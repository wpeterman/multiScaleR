#include <Rcpp.h>

#include <cmath>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

bool is_valid_number(double x) {
  return !NumericVector::is_na(x) && std::isfinite(x);
}

double na_real() {
  return NumericVector::get_na();
}

double min_perimeter(double n_cells) {
  if (n_cells <= 0) {
    return na_real();
  }

  double total_n = std::trunc(std::sqrt(n_cells));
  double total_m = n_cells - total_n * total_n;

  if (total_m == 0) {
    return total_n * 4;
  }
  if (n_cells <= total_n * (1 + total_n)) {
    return 4 * total_n + 2;
  }
  return 4 * total_n + 4;
}

struct CellValue {
  int cell;
  double value;
};

std::vector<CellValue> collect_valid_cells(const NumericVector& d,
                                           const NumericVector& cells,
                                           const NumericMatrix& values,
                                           int column,
                                           double radius) {
  std::vector<CellValue> out;
  out.reserve(d.size());

  for (R_xlen_t i = 0; i < d.size(); ++i) {
    if (!is_valid_number(d[i]) || d[i] > radius ||
        !is_valid_number(cells[i])) {
      continue;
    }

    double value = values(i, column);
    if (!is_valid_number(value)) {
      continue;
    }

    out.push_back(CellValue{static_cast<int>(cells[i]), value});
  }

  return out;
}

std::unordered_map<int, int> index_cells(const std::vector<CellValue>& cells) {
  std::unordered_map<int, int> out;
  out.reserve(cells.size() * 2 + 1);

  for (std::size_t i = 0; i < cells.size(); ++i) {
    out[cells[i].cell] = static_cast<int>(i);
  }

  return out;
}

bool has_neighbor(const std::unordered_map<int, int>& cell_index,
                  int neighbor_cell,
                  int* neighbor_index = nullptr) {
  auto it = cell_index.find(neighbor_cell);
  if (it == cell_index.end()) {
    return false;
  }
  if (neighbor_index != nullptr) {
    *neighbor_index = it->second;
  }
  return true;
}

} // namespace


// [[Rcpp::export]]
NumericVector landscape_composition_metric_cpp(NumericVector d,
                                               NumericMatrix values,
                                               double radius,
                                               std::string metric,
                                               Nullable<NumericVector> weights_,
                                               double base,
                                               Nullable<double> resolution_,
                                               Nullable<double> classes_max_) {
  if (values.nrow() != d.size()) {
    stop("`values` must have one row for each distance in `d`.");
  }
  if (!is_valid_number(radius) || radius <= 0) {
    stop("`radius` must be > 0.");
  }
  if (!is_valid_number(base) || base <= 0 || base == 1) {
    stop("`base` must be > 0 and must not equal 1.");
  }

  bool has_weights = weights_.isNotNull();
  NumericVector weights;
  if (has_weights) {
    weights = NumericVector(weights_);
    if (weights.size() != d.size()) {
      stop("`weights` must have the same length as `d`.");
    }
  }

  bool needs_resolution = metric == "prd" || metric == "ta";
  double resolution = NA_REAL;
  if (needs_resolution) {
    if (resolution_.isNull()) {
      stop("`resolution` is required for PRD and TA.");
    }
    resolution = as<double>(resolution_);
    if (!is_valid_number(resolution) || resolution <= 0) {
      stop("`resolution` must be > 0.");
    }
  }

  double classes_max = NA_REAL;
  bool has_classes_max = classes_max_.isNotNull();
  if (has_classes_max) {
    classes_max = as<double>(classes_max_);
    if (!is_valid_number(classes_max) || classes_max <= 0) {
      stop("`classes_max` must be > 0.");
    }
  }

  NumericVector out(values.ncol(), NA_REAL);
  double log_base = std::log(base);

  for (int col = 0; col < values.ncol(); ++col) {
    std::map<double, double> totals;
    double total_weight = 0;

    for (R_xlen_t i = 0; i < d.size(); ++i) {
      if (!is_valid_number(d[i]) || d[i] > radius) {
        continue;
      }

      double value = values(i, col);
      if (!is_valid_number(value)) {
        continue;
      }

      double weight = has_weights ? weights[i] : 1.0;
      if (!is_valid_number(weight) || weight <= 0) {
        continue;
      }

      totals[value] += weight;
      total_weight += weight;
    }

    if (totals.empty() || total_weight <= 0) {
      continue;
    }

    double shdi = 0;
    double sum_sq = 0;
    for (const auto& kv : totals) {
      double p = kv.second / total_weight;
      if (p > 0) {
        shdi -= p * std::log(p) / log_base;
        sum_sq += p * p;
      }
    }

    double richness = static_cast<double>(totals.size());
    double sidi = 1 - sum_sq;
    double msidi = -std::log(sum_sq);
    double area_ha = needs_resolution ?
      total_weight * resolution * resolution / 10000.0 :
      NA_REAL;

    if (metric == "shdi") {
      out[col] = shdi;
    } else if (metric == "shei") {
      out[col] = richness == 1 ? 0 : shdi / (std::log(richness) / log_base);
    } else if (metric == "sidi") {
      out[col] = sidi;
    } else if (metric == "siei") {
      out[col] = richness <= 1 ? NA_REAL : sidi / (1 - (1 / richness));
    } else if (metric == "msidi") {
      out[col] = msidi;
    } else if (metric == "msiei") {
      out[col] = msidi / std::log(richness);
    } else if (metric == "pr") {
      out[col] = richness;
    } else if (metric == "prd") {
      out[col] = richness / area_ha * 100;
    } else if (metric == "rpr") {
      out[col] = has_classes_max ? richness / classes_max * 100 : NA_REAL;
    } else if (metric == "ta") {
      out[col] = area_ha;
    } else {
      stop("Unsupported composition metric.");
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericVector landscape_edge_metric_cpp(NumericVector d,
                                        NumericMatrix values,
                                        NumericVector cells,
                                        double radius,
                                        double resolution,
                                        int n_cols,
                                        std::string metric) {
  if (values.nrow() != d.size()) {
    stop("`values` must have one row for each distance in `d`.");
  }
  if (cells.size() != d.size()) {
    stop("`cells` must have the same length as `d`.");
  }
  if (!is_valid_number(radius) || radius <= 0) {
    stop("`radius` must be > 0.");
  }
  if (!is_valid_number(resolution) || resolution <= 0) {
    stop("`resolution` must be > 0.");
  }
  if (n_cols <= 0) {
    stop("`n_cols` must be > 0.");
  }

  NumericVector out(values.ncol(), NA_REAL);

  for (int col = 0; col < values.ncol(); ++col) {
    std::vector<CellValue> selected = collect_valid_cells(
      d, cells, values, col, radius
    );
    if (selected.empty()) {
      continue;
    }

    std::unordered_map<int, int> cell_index = index_cells(selected);
    double internal_edges = 0;
    double boundary_edges = 0;

    for (std::size_t i = 0; i < selected.size(); ++i) {
      int cell = selected[i].cell;
      int raster_col = ((cell - 1) % n_cols) + 1;
      int neighbor_index = -1;

      if (raster_col < n_cols &&
          has_neighbor(cell_index, cell + 1, &neighbor_index)) {
        if (selected[i].value != selected[neighbor_index].value) {
          internal_edges += 1;
        }
      }
      if (has_neighbor(cell_index, cell + n_cols, &neighbor_index)) {
        if (selected[i].value != selected[neighbor_index].value) {
          internal_edges += 1;
        }
      }

      int valid_neighbors = 0;
      if (raster_col > 1 && has_neighbor(cell_index, cell - 1)) {
        valid_neighbors += 1;
      }
      if (raster_col < n_cols && has_neighbor(cell_index, cell + 1)) {
        valid_neighbors += 1;
      }
      if (cell > n_cols && has_neighbor(cell_index, cell - n_cols)) {
        valid_neighbors += 1;
      }
      if (has_neighbor(cell_index, cell + n_cols)) {
        valid_neighbors += 1;
      }

      boundary_edges += 4 - valid_neighbors;
    }

    double n_valid = static_cast<double>(selected.size());
    double area_ha = n_valid * resolution * resolution / 10000.0;

    if (metric == "ed") {
      out[col] = internal_edges * resolution / area_ha;
    } else if (metric == "te") {
      out[col] = internal_edges * resolution;
    } else if (metric == "lsi") {
      out[col] = (internal_edges + boundary_edges) / min_perimeter(n_valid);
    } else {
      stop("Unsupported edge metric.");
    }
  }

  return out;
}


// [[Rcpp::export]]
NumericVector landscape_adjacency_metric_cpp(NumericVector d,
                                             NumericMatrix values,
                                             NumericVector cells,
                                             double radius,
                                             int n_cols,
                                             std::string metric) {
  if (values.nrow() != d.size()) {
    stop("`values` must have one row for each distance in `d`.");
  }
  if (cells.size() != d.size()) {
    stop("`cells` must have the same length as `d`.");
  }
  if (!is_valid_number(radius) || radius <= 0) {
    stop("`radius` must be > 0.");
  }
  if (n_cols <= 0) {
    stop("`n_cols` must be > 0.");
  }

  NumericVector out(values.ncol(), NA_REAL);

  for (int col = 0; col < values.ncol(); ++col) {
    std::vector<CellValue> selected = collect_valid_cells(
      d, cells, values, col, radius
    );
    if (selected.empty()) {
      continue;
    }

    std::unordered_map<int, int> cell_index = index_cells(selected);
    std::map<double, int> class_index;

    for (const auto& x : selected) {
      if (class_index.find(x.value) == class_index.end()) {
        int next_index = static_cast<int>(class_index.size());
        class_index[x.value] = next_index;
      }
    }

    int n_classes = static_cast<int>(class_index.size());
    std::vector<double> adjacency(n_classes * n_classes, 0.0);

    for (std::size_t i = 0; i < selected.size(); ++i) {
      int cell = selected[i].cell;
      int raster_col = ((cell - 1) % n_cols) + 1;
      int neighbor_index = -1;

      auto add_pair = [&](int j) {
        int a = class_index[selected[i].value];
        int b = class_index[selected[j].value];
        adjacency[a * n_classes + b] += 1;
        adjacency[b * n_classes + a] += 1;
      };

      if (raster_col < n_cols &&
          has_neighbor(cell_index, cell + 1, &neighbor_index)) {
        add_pair(neighbor_index);
      }
      if (has_neighbor(cell_index, cell + n_cols, &neighbor_index)) {
        add_pair(neighbor_index);
      }
    }

    double total = 0;
    double diag_total = 0;
    for (int i = 0; i < n_classes; ++i) {
      for (int j = 0; j < n_classes; ++j) {
        double value = adjacency[i * n_classes + j];
        total += value;
        if (i == j) {
          diag_total += value;
        }
      }
    }

    if (metric == "pladj") {
      out[col] = total == 0 ? 0 : diag_total / total * 100;
    } else if (metric == "contag") {
      if (n_classes < 2 || total == 0) {
        out[col] = NA_REAL;
      } else {
        double entropy_sum = 0;
        for (double count : adjacency) {
          if (count <= 0) {
            continue;
          }
          double p = count / total;
          entropy_sum += p * std::log(p);
        }
        out[col] = (1 + entropy_sum / (2 * std::log(n_classes))) * 100;
      }
    } else {
      stop("Unsupported adjacency metric.");
    }
  }

  return out;
}
