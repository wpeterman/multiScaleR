"""Interpret retained spatial summaries with validated mean-precision semantics.

Never alters an original result, reruns sampling, substitutes Rhat, or pools chains.
"""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ORIGINAL = HERE.parent / "spatial_reference_20261028"
sys.dont_write_bytecode = True
spec = importlib.util.spec_from_file_location("frozen_spatial_consumer", ORIGINAL / "analyze.py")
comparison = importlib.util.module_from_spec(spec)
spec.loader.exec_module(comparison)


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def mean_event(row):
    finite = comparison.finite
    available = row.get("constant") is False and finite(row.get("estimate")) and 0 <= row["estimate"] <= 1 \
        and finite(row.get("threshold")) and finite(row.get("mcse")) and row["mcse"] > 0 \
        and finite(row.get("ess_mean")) and row["ess_mean"] > 0
    # Respect the unchanged <=0.005 criterion and stay conservative at the
    # precision boundary of stored probability-axis summaries.
    boundary = available and abs(row["mcse"] - .005) <= comparison.ROUNDING_HALF_UNIT
    passed = available and not boundary and row["mcse"] <= .005
    return {"original_precision_pass": row.get("precision_pass"), "mean_precision_available": available,
        "declared_mcse_target_met": passed, "serialization_boundary": boundary,
        "indicator_rhat_available": finite(row.get("rhat"))}


def slot(record, expected_events=228):
    fixed = record.get("fixed_functional", [])
    own = record.get("own_functional", [])
    coverage = all(len(rows) == expected_events and len({comparison.event_key(row) for row in rows}) == expected_events
                   for rows in (fixed, own))
    fixed_targets = {(r["variable"], r["quantile"]) for r in fixed}
    own_targets = {(r["variable"], r["quantile"]) for r in own}
    coverage = coverage and len(fixed_targets) == len(own_targets) == expected_events and fixed_targets == own_targets
    fixed_pass = coverage and all(mean_event(row)["declared_mcse_target_met"] for row in fixed)
    own_pass = coverage and all(mean_event(row)["declared_mcse_target_met"] for row in own)
    other_gates = all(record.get(key) is True for key in ("completed", "sampling_completed", "derived_pass",
                                                       "source_unchanged", "producer_unchanged"))
    raw = record.get("raw_diagnostics", {}).get("passed") is True
    differences = []
    for label, rows in (("fixed", fixed), ("own", own)):
        for row in rows:
            interpretation = mean_event(row)
            if interpretation["original_precision_pass"] != interpretation["declared_mcse_target_met"]:
                fields = {key: row.get(key) for key in
                    ("variable", "quantile", "threshold", "estimate", "mcse", "ess_mean", "rhat", "constant")}
                differences.append({"table": label, **fields, **interpretation})
    return {"original_reference_qualified": record.get("reference_qualified"), "raw_gate_passed": raw,
        "other_original_gates_passed": other_gates, "complete_event_coverage": coverage,
        "fixed_declared_mean_precision_passed": fixed_pass, "own_declared_mean_precision_passed": own_pass,
        "declared_selected_reference_gate_passed": raw and other_gates and fixed_pass and own_pass,
        "divergences": record.get("raw_diagnostics", {}).get("divergences"),
        "availability_differences": differences}


def pair(a, b, expected_events=228):
    aligned = a.get("case") == b.get("case") and a.get("input_hash") == b.get("input_hash") and \
        a.get("settings") == b.get("settings") and a.get("sampler_seed") != b.get("sampler_seed")
    left = {comparison.event_key(row): row for row in a.get("fixed_functional", [])}
    right = {comparison.event_key(row): row for row in b.get("fixed_functional", [])}
    aligned = aligned and len(left) == len(right) == expected_events and left.keys() == right.keys()
    events = []
    if aligned:
        for key, first in left.items():
            second = right[key]
            # This temporary view is used only by the frozen mathematical
            # agreement screen. The original row/flag is never changed.
            first_view = dict(first, precision_pass=mean_event(first)["declared_mcse_target_met"])
            second_view = dict(second, precision_pass=mean_event(second)["declared_mcse_target_met"])
            events.append({"variable": key[0], "quantile": key[1], "threshold": key[2],
                **comparison.event_agreement(first_view, second_view)})
    gates = slot(a, expected_events)["declared_selected_reference_gate_passed"] and \
        slot(b, expected_events)["declared_selected_reference_gate_passed"]
    agrees = aligned and bool(events) and all(row["agreement"] for row in events)
    return {"matched_identity_and_event_sets": aligned, "both_declared_run_gates_passed": gates,
        "all_fixed_events_agree": agrees, "agreement_count": sum(row["agreement"] for row in events),
        "declared_selected_reference_gate_passed": gates and agrees,
        "unavailable_or_disagreeing": [row for row in events if not row["agreement"]]}


def build():
    plan = read(ORIGINAL / "RUNLIST.json")
    manifest = read(ORIGINAL / "ARTIFACT_MANIFEST.json")
    if any(digest(ORIGINAL / row["file"]) != row["sha256"] for row in manifest):
        raise ValueError("Frozen spatial evidence changed")
    records = {}
    for run in plan["runs"]:
        record = read(ORIGINAL / (run["id"] + "-v2.json"))
        if record.get("input_hash") != run["input_hash"] or record.get("sampler_seed") != run["seed"] or \
                record.get("settings") != plan["settings"] or record.get("case") != run["case"]:
            raise ValueError("Frozen scientific input/seed/settings identity differs")
        records[run["id"]] = record
    return {"scope": "Separately labeled mean-precision interpretation of all four original fits; no new sampling or overwritten decision.",
        "planned_scientific_slots": 4, "original_reference_flags": {name: r["reference_qualified"] for name, r in records.items()},
        "source_sha256": {name: digest(ORIGINAL / (name + "-v2.json")) for name in records},
        "slots": {name: slot(record) for name, record in records.items()},
        "pairs": {case: pair(records[case + "-a"], records[case + "-b"]) for case in ("null", "smooth")},
        "unchanged_criteria": "Original raw/derived gates; finite nonconstant mean MCSE<=0.005; both independent runs and every fixed-event three-combined-MCSE agreement.",
        "limitations": ["Indicator Rhat remains unavailable where originally undefined; it is never substituted.",
            "Selected-case numerical qualification does not establish calibration, nonstationarity detection or public spatial readiness.",
            "Historical Laplace CDFs retain their original serialization/evaluation qualification.",
            "The null pair retains five/three divergences and remains unqualified."]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    if args.write:
        unit = read(HERE / "sidecar-unit02.json")
        if not unit["passed"] or digest(HERE / "spatial_sidecar.py").upper() != unit["sidecar_sha256"]:
            raise ValueError("Sidecar changed since focused tests")
        subprocess.run([sys.executable, str(HERE / "verify-candidate.py")], check=True)
    result = build()
    if args.write:
        destination = HERE / "spatial_probability_sidecar.json"
        if destination.exists():
            raise ValueError("Retain existing sidecar")
        temporary = destination.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")
        temporary.replace(destination)
    print(json.dumps({"original_reference_flags": result["original_reference_flags"], "pairs": result["pairs"],
        "availability_differences": {key: len(value["availability_differences"]) for key, value in result["slots"].items()}}, indent=2))
