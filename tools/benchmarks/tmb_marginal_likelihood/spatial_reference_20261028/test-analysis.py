"""Focused deterministic consumer tests; no sampling or package mutation."""
import copy
import importlib.util
import math
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location("spatial_analysis", Path(__file__).with_name("analyze.py"))
analysis = importlib.util.module_from_spec(spec)
spec.loader.exec_module(analysis)


def fixture(seed):
    output = {key: True for key in ("completed", "sampling_completed", "derived_pass", "fixed_functional_pass",
        "own_functional_pass", "source_unchanged", "producer_unchanged", "reference_qualified")}
    output.update(case="null", sampler_seed=seed, input_hash="fixed", settings={"iter": 16000},
        raw_diagnostics={"passed": True}, fixed_functional=[{"variable": "field_sd", "quantile": .5,
            "threshold": .1, "estimate": .5, "mcse": .001, "constant": False, "precision_pass": True}])
    return output


class ConsumerContracts(unittest.TestCase):
    def test_matching_independent_qualified_pair(self):
        a, b = fixture(1), fixture(2)
        original = copy.deepcopy((a, b))
        result = analysis.compare_pair(a, b, expected_events=1)
        self.assertTrue(result["qualified"])
        self.assertEqual(result["agreement_count"], 1)
        self.assertEqual((a, b), original)

    def test_missing_run_stays_pending(self):
        self.assertEqual(analysis.compare_pair(fixture(1), None)["state"], "pending")
        self.assertFalse(analysis.run_qualified(None))

    def test_each_original_gate_blocks_qualification(self):
        for key in ("completed", "sampling_completed", "derived_pass", "fixed_functional_pass",
                    "own_functional_pass", "source_unchanged", "producer_unchanged", "reference_qualified"):
            a = fixture(1); a[key] = False
            self.assertFalse(analysis.compare_pair(a, fixture(2), expected_events=1)["qualified"], key)
        a = fixture(1); a["raw_diagnostics"]["passed"] = False
        self.assertFalse(analysis.compare_pair(a, fixture(2), expected_events=1)["qualified"])

    def test_run_identity_and_seed_mismatch(self):
        for key, value in (("case", "smooth"), ("sampler_seed", 2), ("input_hash", "other"), ("settings", {})):
            a = fixture(1); a[key] = value
            self.assertEqual(analysis.compare_pair(a, fixture(2), expected_events=1)["state"], "identity_failure")

    def test_threshold_sets_are_exact_complete_unique(self):
        for mutation in (lambda rows: rows.clear(), lambda rows: rows.append(dict(rows[0])),
                         lambda rows: rows[0].update(threshold=.100001), lambda rows: rows[0].update(quantile=.975)):
            a = fixture(1); mutation(a["fixed_functional"])
            self.assertEqual(analysis.compare_pair(a, fixture(2), expected_events=1)["state"], "event_identity_failure")

    def test_probability_disagreement_is_retained(self):
        a = fixture(1); a["fixed_functional"][0]["estimate"] = .51
        result = analysis.compare_pair(a, fixture(2), expected_events=1)
        self.assertFalse(result["qualified"])
        self.assertEqual(result["disagreement_count"], 1)

    def test_constant_nonfinite_and_failed_precision_unavailable(self):
        for key, value in (("constant", True), ("precision_pass", False), ("mcse", None), ("mcse", 0),
                           ("mcse", float("nan")), ("estimate", float("inf")), ("estimate", -1), ("mcse", True)):
            a = fixture(1); a["fixed_functional"][0][key] = value
            result = analysis.compare_pair(a, fixture(2), expected_events=1)
            self.assertFalse(result["qualified"])
            self.assertEqual(result["unavailable_count"], 1)

    def test_serialization_boundary_never_becomes_agreement(self):
        a = fixture(1)
        a["fixed_functional"][0]["estimate"] += 3 * math.hypot(.001, .001)
        result = analysis.compare_pair(a, fixture(2), expected_events=1)
        self.assertFalse(result["qualified"])
        self.assertEqual(result["serialization_boundary_count"], 1)

    def test_failed_execution_keeps_unmeasured_precision_unavailable(self):
        record = fixture(1)
        record.update(completed=False, error="execution failed", fixed_functional=[], own_functional=[])
        result = analysis.summarize_record(record)
        self.assertEqual(result["state"], "failed")
        self.assertEqual(result["error"], "execution failed")
        self.assertEqual(result["fixed_events_evaluated"], 0)
        self.assertIsNone(result["fixed_precision_failures"])
        self.assertIsNone(result["own_precision_failures"])
        self.assertFalse(result["reference_qualified"])

    def test_finite_mean_precision_and_undefined_rhat_remain_distinct(self):
        record = fixture(1)
        record["own_functional"] = [{"constant": False, "mcse": .002, "ess_mean": 20000,
            "rhat": None, "precision_pass": False}]
        record["reference_qualified"] = False
        result = analysis.summarize_record(record)
        self.assertEqual(result["own_precision_failures"], 1)
        self.assertEqual(result["own_mean_precision_nonestimable"], 0)
        self.assertEqual(result["own_mean_mcse_above_target"], 0)
        self.assertEqual(result["own_finite_mean_precision_unavailable_rhat"], 1)
        self.assertFalse(result["reference_qualified"])


if __name__ == "__main__":
    unittest.main()
