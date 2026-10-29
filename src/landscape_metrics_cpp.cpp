#include <algorithm>
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

double max_like_adjacencies(double n_cells) {
  if (n_cells <= 0) {
    return na_real();
  }

  double n = std::trunc(std::sqrt(n_cells));
  double m = n_cells - n * n;

  if (m == 0) {
    return 2 * n * (n - 1);
  }
  if (m <= n) {
    return 2 * n * (n - 1) + 2 * m - 1;
  }
  return 2 * n * (n - 1) + 2 * m - 2;
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
                                               Nullable<double> classes_max_,
                                               Nullable<double> focal_class_) {
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

  bool needs_resolution = metric == "prd" || metric == "ta" || metric == "ca";
  double resolution = NA_REAL;
  if (needs_resolution) {
    if (resolution_.isNull()) {
      stop("`resolution` is required for PRD, TA, and CA.");
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

  // Class-level composition metrics (PLAND, CA) target one focal class.
  bool needs_focal = metric == "pland" || metric == "ca";
  double focal_class = NA_REAL;
  if (needs_focal) {
    if (focal_class_.isNull()) {
      stop("`class` is required for the class-level metrics PLAND and CA.");
    }
    focal_class = as<double>(focal_class_);
    if (!is_valid_number(focal_class)) {
      stop("`class` must be a finite categorical value.");
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
    } else if (metric == "pland" || metric == "ca") {
      auto it = totals.find(focal_class);
      double focal_weight = it == totals.end() ? 0.0 : it->second;
      if (metric == "pland") {
        out[col] = focal_weight / total_weight * 100;
      } else {
        out[col] = focal_weight * resolution * resolution / 10000.0;
      }
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


// Shared formulas keep prepared and reference paths statistically identical.
static double adjacency_from_counts(const std::vector<double>& adjacency,
    const std::vector<double>& class_counts, const std::map<double,int>& class_index,
    double total_cells, const std::string& metric, double base, double focal_class) {
  int n_classes = class_counts.size();
  double log_base = std::log(base);
  double result = NA_REAL;
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

    if (metric == "ai") {
      double weighted_ai = 0;

      for (int i = 0; i < n_classes; ++i) {
        double max_adj = max_like_adjacencies(class_counts[i]);
        if (!is_valid_number(max_adj) || max_adj <= 0) {
          continue;
        }

        double like_adjacencies = adjacency[i * n_classes + i] / 2.0;
        double class_ai = like_adjacencies / max_adj * 100.0;
        if (is_valid_number(class_ai)) {
          weighted_ai += class_ai * class_counts[i] / total_cells;
        }
      }

      result = weighted_ai;
    } else if (metric == "pladj") {
      result = total == 0 ? 0 : diag_total / total * 100;
    } else if (metric == "contag") {
      if (n_classes < 2 || total == 0) {
        result = NA_REAL;
      } else {
        double entropy_sum = 0;
        for (double count : adjacency) {
          if (count <= 0) {
            continue;
          }
          double p = count / total;
          entropy_sum += p * std::log(p);
        }
        result = (1 + entropy_sum / (2 * std::log(n_classes))) * 100;
      }
    } else if (metric == "iji") {
      // Interspersion & juxtaposition: entropy of the between-class adjacency
      // distribution, normalised by its maximum. Undefined for < 3 classes.
      if (n_classes < 3) {
        result = NA_REAL;
      } else {
        double between_total = 0;
        for (int a = 0; a < n_classes; ++a) {
          for (int b = a + 1; b < n_classes; ++b) {
            between_total += adjacency[a * n_classes + b];
          }
        }
        if (between_total <= 0) {
          result = NA_REAL;
        } else {
          double e_sum = 0;
          for (int a = 0; a < n_classes; ++a) {
            for (int b = a + 1; b < n_classes; ++b) {
              double e = adjacency[a * n_classes + b];
              if (e > 0) {
                double p = e / between_total;
                e_sum -= p * std::log(p);
              }
            }
          }
          double denom = std::log(0.5 * n_classes * (n_classes - 1));
          result = e_sum / denom * 100;
        }
      }
    } else if (metric == "ent" || metric == "condent" || metric == "joinent" ||
               metric == "mutinf" || metric == "relmutinf") {
      // Information theory on the normalised co-occurrence matrix
      // (Nowosad & Stepinski 2019): marginal entropy `ent`, joint entropy
      // `joinent`, conditional entropy `condent`, mutual information `mutinf`,
      // and relative mutual information `relmutinf`.
      if (total <= 0) {
        result = NA_REAL;
      } else {
        double joinent = 0;
        std::vector<double> row_total(n_classes, 0.0);
        for (int a = 0; a < n_classes; ++a) {
          for (int b = 0; b < n_classes; ++b) {
            double g = adjacency[a * n_classes + b];
            row_total[a] += g;
            if (g > 0) {
              double p = g / total;
              joinent -= p * std::log(p) / log_base;
            }
          }
        }
        double ent = 0;
        for (int a = 0; a < n_classes; ++a) {
          double px = row_total[a] / total;
          if (px > 0) {
            ent -= px * std::log(px) / log_base;
          }
        }
        double condent = joinent - ent;
        double mutinf = ent - condent;
        double relmutinf = mutinf == 0 ? 1.0 : mutinf / ent;
        if (metric == "ent") {
          result = ent;
        } else if (metric == "condent") {
          result = condent;
        } else if (metric == "joinent") {
          result = joinent;
        } else if (metric == "mutinf") {
          result = mutinf;
        } else {
          result = relmutinf;
        }
      }
    } else if (metric == "clumpy") {
      // Clumpiness index for the focal class: like-adjacency proportion given
      // the minimum perimeter, rescaled against the proportion expected under a
      // spatially random distribution. Bounded on [-1, 1].
      auto focal_it = class_index.find(focal_class);
      if (focal_it == class_index.end()) {
        result = NA_REAL;
      } else {
        int fi = focal_it->second;
        double n_i = class_counts[fi];
        double p_i = n_i / total_cells;
        double min_e = min_perimeter(n_i);
        if (p_i >= 1 || !is_valid_number(min_e)) {
          result = NA_REAL;
        } else {
          double like_i = adjacency[fi * n_classes + fi];
          double g_i = like_i / (4.0 * n_i - min_e);
          double clumpy;
          if (g_i >= p_i || p_i >= 0.5) {
            clumpy = (g_i - p_i) / (1 - p_i);
          } else {
            clumpy = (g_i - p_i) / (-p_i);
          }
          result = clumpy;
        }
      }
    } else {
      stop("Unsupported adjacency metric.");
    }
  return result;
}

// [[Rcpp::export]]
NumericVector landscape_adjacency_metric_cpp(NumericVector d,
                                             NumericMatrix values,
                                             NumericVector cells,
                                             double radius,
                                             int n_cols,
                                             std::string metric,
                                             double base,
                                             Nullable<double> focal_class_) {
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
  if (!is_valid_number(base) || base <= 0 || base == 1) {
    stop("`base` must be > 0 and must not equal 1.");
  }

  // CLUMPY is class-level and targets one focal class.
  bool needs_focal = metric == "clumpy";
  double focal_class = NA_REAL;
  if (needs_focal) {
    if (focal_class_.isNull()) {
      stop("`class` is required for the class-level metric CLUMPY.");
    }
    focal_class = as<double>(focal_class_);
    if (!is_valid_number(focal_class)) {
      stop("`class` must be a finite categorical value.");
    }
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
    std::vector<double> class_counts(n_classes, 0.0);
    for (const auto& x : selected) {
      class_counts[class_index[x.value]] += 1;
    }

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

    out[col] = adjacency_from_counts(adjacency, class_counts, class_index,
      static_cast<double>(selected.size()), metric, base, focal_class);
  }

  return out;
}

// Prepare invariant geometry once. Returned ordinary R arrays serialize safely
// to workers; no external pointer or process-local mutable cache is retained.
// [[Rcpp::export]]
List msr_adjacency_prepare_cpp(NumericVector d, NumericVector values,
                              NumericVector cells, int n_cols) {
  if (d.size()!=values.size() || d.size()!=cells.size() || n_cols<1)
    stop("Invalid prepared adjacency dimensions.");
  std::vector<double> distances, labels;
  std::vector<int> ids, codes;
  std::map<double,int> classes;
  for (R_xlen_t i=0;i<d.size();++i) {
    if (!is_valid_number(d[i]) || !is_valid_number(values[i])) continue;
    if (!is_valid_number(cells[i]) || cells[i]<1 || cells[i]>2147483647.0-n_cols ||
        cells[i]!=std::trunc(cells[i])) stop("Invalid prepared cell ID.");
    if (!classes.count(values[i])) {
      if (classes.size() >= 200) stop("Prepared adjacency supports at most 200 classes.");
      int next = classes.size(); classes[values[i]]=next; labels.push_back(values[i]);
    }
    distances.push_back(d[i]); ids.push_back(static_cast<int>(cells[i]));
    codes.push_back(classes[values[i]]);
  }
  std::unordered_map<int,int> index;
  for (size_t i=0;i<ids.size();++i) {
    if (index.count(ids[i])) stop("Duplicate prepared cell ID.");
    index[ids[i]]=i;
  }
  IntegerVector right(ids.size(),-1), down(ids.size(),-1);
  for (size_t i=0;i<ids.size();++i) {
    auto r=index.find(ids[i]+1), b=index.find(ids[i]+n_cols);
    if (ids[i]%n_cols && r!=index.end()) right[i]=r->second;
    if (b!=index.end()) down[i]=b->second;
  }
  return List::create(_["d"]=distances,_["code"]=codes,_["labels"]=labels,
                      _["right"]=right,_["down"]=down);
}

// Only O(classes^2) working storage, independent of the number of selected
// cells and radius proposals. Exact matching radii share one count scan.
// [[Rcpp::export]]
NumericVector msr_adjacency_eval_cpp(List prepared, NumericVector radii,
    CharacterVector metrics, NumericVector bases, NumericVector focals) {
  NumericVector d=prepared["d"], labels=prepared["labels"];
  IntegerVector code=prepared["code"], right=prepared["right"], down=prepared["down"];
  int n=d.size(), k=labels.size(), q=radii.size();
  if (k>200) stop("Prepared adjacency supports at most 200 classes.");
  if(code.size()!=n || right.size()!=n || down.size()!=n || metrics.size()!=q ||
     bases.size()!=q || focals.size()!=q) stop("Invalid prepared adjacency dimensions.");
  for(int i=0;i<n;++i) if(code[i]<0 || code[i]>=k || right[i]<-1 || right[i]>=n || down[i]<-1 || down[i]>=n)
    stop("Invalid prepared adjacency index.");
  NumericVector out(q,NA_REAL);
  std::vector<bool> done(q,false);
  for(int j=0;j<q;++j) {
    if(done[j]) continue;
    double radius=radii[j];
    if(!is_valid_number(radius)||radius<=0) stop("`radius` must be > 0.");
    std::vector<int> local(k,-1);
    std::vector<double> counts;
    std::map<double,int> classes;
    double total_cells=0;
    // Match reference first-occurrence class order at each radius.
    for(int i=0;i<n;++i) if(d[i]<=radius) {
      int c=code[i];
      if(local[c]<0) {local[c]=counts.size();counts.push_back(0);classes[labels[c]]=local[c];}
      counts[local[c]]++;total_cells++;
    }
    int m=counts.size();std::vector<double> adjacency(m*m,0);
    for(int i=0;i<n;++i) if(d[i]<=radius) {
      int a=local[code[i]];
      int neighbors[2]={right[i],down[i]};
      for(int b:neighbors) if(b>=0 && d[b]<=radius) {
        int c=local[code[b]];adjacency[a*m+c]++;adjacency[c*m+a]++;
      }
    }
    for(int z=j;z<q;++z) if(radii[z]==radius) {
      double base=bases[z];std::string metric=as<std::string>(metrics[z]);
      if(!is_valid_number(base)||base<=0||base==1) stop("`base` must be > 0 and must not equal 1.");
      if(metric=="clumpy"&&!is_valid_number(focals[z])) stop("`class` must be a finite categorical value.");
      if(total_cells)out[z]=adjacency_from_counts(adjacency,counts,classes,total_cells,metric,base,focals[z]);
      done[z]=true;
    }
  }
  return out;
}

// Radius-invariant entry thresholds: a pair enters when both cells enter.
// Store only observed class pairs; no cell-by-radius cumulative matrix.
// [[Rcpp::export]]
List msr_adjacency_events_prepare_cpp(List p) {
  NumericVector d=p["d"], labels=p["labels"];
  IntegerVector code=p["code"], right=p["right"], down=p["down"];
  int n=d.size(), k=labels.size();
  if(k>200 || code.size()!=n || right.size()!=n || down.size()!=n)
    stop("Invalid prepared adjacency dimensions.");
  std::vector<std::vector<std::pair<double,int>>> cells(k);
  std::map<int,std::vector<double>> pairs;
  for(int i=0;i<n;++i) {
    if(code[i]<0 || code[i]>=k || right[i]<-1 || right[i]>=n || down[i]<-1 || down[i]>=n)
      stop("Invalid prepared adjacency index.");
  }
  for(int i=0;i<n;++i) {
    cells[code[i]].push_back(std::make_pair(d[i],i));
    int neighbors[2]={right[i],down[i]};
    for(int j:neighbors) if(j>=0) {
      int a=std::min(code[i],code[j]),b=std::max(code[i],code[j]);
      pairs[a*k+b].push_back(std::max(d[i],d[j]));
    }
  }
  std::vector<double> cell_d,pair_d,first_d;
  std::vector<int> first,first_offset(1,0),offset(1,0),pair_key,pair_offset(1,0);
  for(auto &group:cells) {
    std::sort(group.begin(),group.end());
    int earliest=n;
    for(auto &entry:group) {
      cell_d.push_back(entry.first);
      // Only retain changes in the earliest original index, not one index
      // for every cell. This preserves arithmetic order with a small frontier.
      if(entry.second<earliest) {
        earliest=entry.second;first.push_back(earliest);first_d.push_back(entry.first);
      }
    }
    offset.push_back(cell_d.size());
    first_offset.push_back(first.size());
  }
  for(auto &group:pairs) {
    std::sort(group.second.begin(),group.second.end());
    pair_key.push_back(group.first);
    pair_d.insert(pair_d.end(),group.second.begin(),group.second.end());
    pair_offset.push_back(pair_d.size());
  }
  return List::create(_["event_version"]=1,_["labels"]=labels,
    _["cell_d"]=cell_d,_["first"]=first,_["offset"]=offset,
    _["first_d"]=first_d,_["first_offset"]=first_offset,
    _["pair_d"]=pair_d,_["pair_key"]=pair_key,_["pair_offset"]=pair_offset);
}

// [[Rcpp::export]]
NumericVector msr_adjacency_events_eval_cpp(List p, NumericVector radii,
    CharacterVector metrics, NumericVector bases, NumericVector focals) {
  NumericVector labels=p["labels"], cell_d=p["cell_d"], pair_d=p["pair_d"], first_d=p["first_d"];
  IntegerVector first=p["first"], fo=p["first_offset"], offset=p["offset"], keys=p["pair_key"], po=p["pair_offset"];
  int k=labels.size(),q=radii.size(),ng=keys.size();
  if(k>200 || first.size()!=first_d.size() || fo.size()!=k+1 || offset.size()!=k+1 || po.size()!=ng+1 ||
     metrics.size()!=q || bases.size()!=q || focals.size()!=q ||
     offset[0]!=0 || offset[k]!=cell_d.size() || fo[0]!=0 || fo[k]!=first.size() || po[0]!=0 || po[ng]!=pair_d.size())
    stop("Invalid prepared adjacency event dimensions.");
  for(int c=0;c<k;++c) if(offset[c]<0 || offset[c+1]<offset[c] || fo[c]<0 || fo[c+1]<fo[c]) stop("Invalid cell event offset.");
  for(int g=0;g<ng;++g) if(po[g]<0 || po[g+1]<po[g] || keys[g]<0 || keys[g]>=k*k)
    stop("Invalid pair event offset or class.");
  NumericVector out(q,NA_REAL);std::vector<bool> done(q,false);
  for(int j=0;j<q;++j) {
    if(done[j])continue;
    double radius=radii[j];if(!is_valid_number(radius)||radius<=0)stop("`radius` must be > 0.");
    std::vector<double> global_counts(k,0),counts;
    std::vector<std::pair<int,int>> order;
    std::vector<int> local(k,-1);
    std::map<double,int> classes;
    double total_cells=0;
    for(int c=0;c<k;++c) {
      int end=std::upper_bound(cell_d.begin()+offset[c],cell_d.begin()+offset[c+1],radius)-cell_d.begin();
      int count=end-offset[c];global_counts[c]=count;total_cells+=count;
      if(count) {
        int f=std::upper_bound(first_d.begin()+fo[c],first_d.begin()+fo[c+1],radius)-first_d.begin();
        if(f==fo[c])stop("Missing first-occurrence frontier.");
        order.push_back(std::make_pair(first[f-1],c));
      }
    }
    // Preserve the original scan's first-occurrence order, including at ties.
    std::sort(order.begin(),order.end());
    for(auto &entry:order) {
      int c=entry.second;local[c]=counts.size();classes[labels[c]]=local[c];counts.push_back(global_counts[c]);
    }
    int m=counts.size();std::vector<double> adjacency(m*m,0);
    for(int g=0;g<ng;++g) {
      int count=std::upper_bound(pair_d.begin()+po[g],pair_d.begin()+po[g+1],radius)-(pair_d.begin()+po[g]);
      if(count) {
        int a=local[keys[g]/k],b=local[keys[g]%k];
        if(a<0||b<0)stop("Pair events precede their cells.");
        adjacency[a*m+b]+=count;adjacency[b*m+a]+=count;
      }
    }
    for(int z=j;z<q;++z)if(radii[z]==radius) {
      double base=bases[z];std::string metric=as<std::string>(metrics[z]);
      if(!is_valid_number(base)||base<=0||base==1)stop("`base` must be > 0 and must not equal 1.");
      if(metric=="clumpy"&&!is_valid_number(focals[z]))stop("`class` must be a finite categorical value.");
      if(total_cells)out[z]=adjacency_from_counts(adjacency,counts,classes,total_cells,metric,base,focals[z]);
      done[z]=true;
    }
  }
  return out;
}
