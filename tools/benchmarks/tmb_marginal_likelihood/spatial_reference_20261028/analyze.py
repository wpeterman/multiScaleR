"""Consume compact spatial summaries; never fit, pool chains, or alter gates."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROUNDING_HALF_UNIT = 0.5e-12


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def run_qualified(record):
    if not record:
        return False
    return all(record.get(key) is True for key in ("completed", "sampling_completed",
        "derived_pass", "fixed_functional_pass", "own_functional_pass", "source_unchanged",
        "producer_unchanged", "reference_qualified")) and record.get("raw_diagnostics", {}).get("passed") is True


def event_key(row):
    return (row["variable"], row["quantile"], row["threshold"])


def event_agreement(a, b):
    usable = all(row.get("precision_pass") is True and not row.get("constant", True)
        and finite(row.get("estimate")) and 0 <= row["estimate"] <= 1
        and finite(row.get("mcse")) and row["mcse"] > 0 for row in (a, b))
    if not usable:
        return {"state": "unavailable", "agreement": False, "reason": "Nonestimable or failed event precision."}
    difference = abs(a["estimate"] - b["estimate"])
    combined = math.hypot(a["mcse"], b["mcse"])
    # Summaries use jsonlite digits=12. Resolve the original <=3*combined-MCSE
    # screen only when it is unchanged throughout that serialization interval.
    lower_difference = max(0.0, difference - 2 * ROUNDING_HALF_UNIT)
    upper_difference = difference + 2 * ROUNDING_HALF_UNIT
    lower_cutoff = 3 * math.hypot(max(0.0, a["mcse"] - ROUNDING_HALF_UNIT),
                                max(0.0, b["mcse"] - ROUNDING_HALF_UNIT))
    upper_cutoff = 3 * math.hypot(a["mcse"] + ROUNDING_HALF_UNIT, b["mcse"] + ROUNDING_HALF_UNIT)
    state = "agreement" if upper_difference <= lower_cutoff else \
        "disagreement" if lower_difference > upper_cutoff else "serialization_boundary"
    return {"state": state, "agreement": state == "agreement", "absolute_difference": difference,
            "combined_mcse": combined, "three_combined_mcse": 3 * combined}


def compare_pair(a, b, expected_events=228):
    if a is None or b is None:
        return {"qualified": False, "state": "pending", "events": [], "reason": "Both independent runs are required."}
    if a.get("case") != b.get("case") or a.get("input_hash") != b.get("input_hash") or \
            a.get("settings") != b.get("settings") or a.get("sampler_seed") == b.get("sampler_seed"):
        return {"qualified": False, "state": "identity_failure", "events": [],
                "reason": "Matched target/input/settings and distinct planned seeds are required."}
    rows_a = a.get("fixed_functional", [])
    rows_b = b.get("fixed_functional", [])
    mapping_a = {event_key(row): row for row in rows_a}
    mapping_b = {event_key(row): row for row in rows_b}
    if len(rows_a) != expected_events or len(rows_b) != expected_events or \
            len(mapping_a) != expected_events or len(mapping_b) != expected_events or mapping_a.keys() != mapping_b.keys():
        return {"qualified": False, "state": "event_identity_failure", "events": [],
                "reason": "The complete unique fixed-threshold event sets must match exactly."}
    events = []
    for key, row in mapping_a.items():
        result = event_agreement(row, mapping_b[key])
        events.append({"variable": key[0], "quantile": key[1], "threshold": key[2], **result})
    gates = run_qualified(a) and run_qualified(b)
    agreement = all(event["agreement"] for event in events)
    return {"qualified": gates and agreement, "state": "qualified" if gates and agreement else "not_qualified",
            "both_run_gates_passed": gates, "all_fixed_events_agree": agreement,
            "agreement_count": sum(e["state"] == "agreement" for e in events),
            "disagreement_count": sum(e["state"] == "disagreement" for e in events),
            "unavailable_count": sum(e["state"] == "unavailable" for e in events),
            "serialization_boundary_count": sum(e["state"] == "serialization_boundary" for e in events),
            "events": events}


def summarize_record(record):
    if record is None:
        return {"state": "pending", "reference_qualified": False}
    raw = record.get("raw_diagnostics", {})
    fixed = record.get("fixed_functional", [])
    own = record.get("own_functional", [])
    def mean_precision(row):
        return row.get("constant") is False and finite(row.get("mcse")) and row["mcse"] > 0 \
            and finite(row.get("ess_mean")) and row["ess_mean"] > 0
    return {"state": "completed" if record.get("completed") else "failed", "error": record.get("error"),
        "reference_qualified": run_qualified(record), "sampling_completed": record.get("sampling_completed", False),
        "raw_diagnostics": raw, "derived_pass": record.get("derived_pass"),
        "fixed_events_evaluated": len(fixed), "own_events_evaluated": len(own),
        "fixed_precision_failures": sum(row.get("precision_pass") is not True for row in fixed) if fixed else None,
        "own_precision_failures": sum(row.get("precision_pass") is not True for row in own) if own else None,
        "own_mean_precision_nonestimable": sum(not mean_precision(row) for row in own) if own else None,
        "own_mean_mcse_above_target": sum(mean_precision(row) and row["mcse"] > .005 for row in own) if own else None,
        "own_finite_mean_precision_unavailable_rhat": sum(mean_precision(row) and not finite(row.get("rhat")) for row in own) if own else None,
        "historical_laplace_flags": sum(row.get("historical_difference_flag") is True for row in fixed),
        "historical_laplace_unavailable": sum(row.get("historical_difference_flag") is None for row in fixed),
        "geometry": record.get("geometry"), "elapsed_seconds": record.get("elapsed_seconds"),
        "warnings": record.get("warnings", []), "warning_count": record.get("warning_count", 0)}


def analyze():
    plan = read(HERE / "RUNLIST.json")
    records = {}
    for run in plan["runs"]:
        path = HERE / (run["id"] + "-v2.json")
        record = read(path) if path.is_file() else None
        if record is not None and (record.get("run") != run["id"] or record.get("case") != run["case"] or
                record.get("sampler_seed") != run["seed"] or
                (record.get("completed") and record.get("input_hash") != run["input_hash"]) or
                record.get("settings") != plan["settings"]):
            raise ValueError("Run record differs from the frozen scientific plan: " + run["id"])
        records[run["id"]] = record
    pairs = {case: compare_pair(records[case + "-a"], records[case + "-b"]) for case in ("null", "smooth")}
    check = subprocess.run([sys.executable, str(HERE / "verify-repair.py")], capture_output=True, text=True)
    identity_pass = check.returncode == 0
    if not identity_pass:
        for pair in pairs.values():
            pair.update(qualified=False, identity_gate_failed=True)
    return {"planned_runs": 4, "completed_runs": sum(bool(r and r.get("completed")) for r in records.values()),
        "failed_runs": sum(bool(r and not r.get("completed")) for r in records.values()),
        "pending_runs": sum(r is None for r in records.values()), "identity_gate_passed": identity_pass,
        "runs": {key: summarize_record(record) for key, record in records.items()}, "pairs": pairs,
        "prior_pre_sampling_execution_failure_retained": True,
        "scope": "Two selected synthetic datasets with two independent seeds each; no calibration or public spatial readiness.",
        "serialization_note": "New pairs use identical stored thresholds. Near-boundary 12-decimal summary rounding remains unavailable, never rounded into agreement.",
        "historical_comparison_note": "Retained Laplace CDFs were not recomputed at rounded stored thresholds. Flags are descriptive and require qualified reference pairs."}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Save one final compact analysis after all four executions finish")
    args = parser.parse_args()
    result = analyze()
    if args.write:
        panel_complete = any((HERE / name).is_file() for name in ("panel-v2-complete.json", "panel-v3-complete.json"))
        if not panel_complete or result["pending_runs"] or not result["identity_gate_passed"]:
            raise RuntimeError("Final analysis requires complete planned outcomes and unchanged identities")
        path = HERE / "paired-analysis.json"
        if path.exists():
            raise RuntimeError("Retain existing completed analysis")
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
        temporary.replace(path)
    compact = {key: value for key, value in result.items() if key not in ("runs", "pairs")}
    compact["runs"] = {key: {k: v for k, v in row.items() if k not in ("raw_diagnostics", "geometry", "warnings")}
                       for key, row in result["runs"].items()}
    compact["pairs"] = {key: {k: v for k, v in row.items() if k != "events"} for key, row in result["pairs"].items()}
    print(json.dumps(compact, indent=2, allow_nan=False))
