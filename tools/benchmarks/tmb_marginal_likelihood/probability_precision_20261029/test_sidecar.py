import copy
import unittest
import spatial_sidecar as sidecar


def fixture(seed=1):
    event = {"variable": "field_sd", "quantile": .5, "threshold": .1, "estimate": .5,
             "mcse": .001, "ess_mean": 20000, "rhat": None, "constant": False, "precision_pass": False}
    record = {key: True for key in ("completed", "sampling_completed", "derived_pass", "source_unchanged", "producer_unchanged")}
    record.update(reference_qualified=False, raw_diagnostics={"passed": True, "divergences": 0},
        case="smooth", input_hash="fixed", sampler_seed=seed, settings={"chains": 4},
        fixed_functional=[dict(event)], own_functional=[dict(event)])
    return record


class SidecarContracts(unittest.TestCase):
    def test_availability_is_separate_without_mutation(self):
        row = fixture()["own_functional"][0]
        before = copy.deepcopy(row)
        result = sidecar.mean_event(row)
        self.assertTrue(result["mean_precision_available"])
        self.assertTrue(result["declared_mcse_target_met"])
        self.assertFalse(result["indicator_rhat_available"])
        self.assertFalse(result["original_precision_pass"])
        self.assertEqual(row, before)

    def test_original_flag_retained_and_selected_gate_separate(self):
        record = fixture(); before = copy.deepcopy(record)
        result = sidecar.slot(record, 1)
        self.assertFalse(result["original_reference_qualified"])
        self.assertTrue(result["declared_selected_reference_gate_passed"])
        self.assertEqual(len(result["availability_differences"]), 2)
        self.assertEqual(record, before)

    def test_each_original_run_gate_still_required(self):
        for key in ("completed", "sampling_completed", "derived_pass", "source_unchanged", "producer_unchanged"):
            record = fixture(); record[key] = False
            self.assertFalse(sidecar.slot(record, 1)["declared_selected_reference_gate_passed"], key)
        record = fixture(); record["raw_diagnostics"].update(passed=False, divergences=5)
        self.assertFalse(sidecar.slot(record, 1)["declared_selected_reference_gate_passed"])
        self.assertEqual(sidecar.slot(record, 1)["divergences"], 5)

    def test_constant_nonfinite_or_invalid_precision_stays_unavailable(self):
        for key, value in (("constant", True), ("mcse", None), ("mcse", 0), ("mcse", -1),
                           ("mcse", float("inf")), ("ess_mean", None), ("ess_mean", 0),
                           ("estimate", None), ("estimate", -1), ("estimate", 2), ("threshold", None)):
            event = fixture()["own_functional"][0]; event[key] = value
            self.assertFalse(sidecar.mean_event(event)["declared_mcse_target_met"])

    def test_unchanged_mcse_cutoff_and_stored_boundary(self):
        event = fixture()["own_functional"][0]
        event["mcse"] = .006
        self.assertFalse(sidecar.mean_event(event)["declared_mcse_target_met"])
        event["mcse"] = .005
        self.assertTrue(sidecar.mean_event(event)["serialization_boundary"])
        self.assertFalse(sidecar.mean_event(event)["declared_mcse_target_met"])

    def test_complete_unique_event_coverage_required(self):
        record = fixture(); record["own_functional"] = []
        self.assertFalse(sidecar.slot(record, 1)["declared_selected_reference_gate_passed"])
        record = fixture(); record["own_functional"] *= 2
        self.assertFalse(sidecar.slot(record, 1)["declared_selected_reference_gate_passed"])
        record = fixture()
        for table in ("fixed_functional", "own_functional"):
            record[table].append(dict(record[table][0], threshold=.2))
        self.assertFalse(sidecar.slot(record, 2)["declared_selected_reference_gate_passed"])

    def test_independent_matched_pair_passes_declared_screen(self):
        self.assertTrue(sidecar.pair(fixture(1), fixture(2), 1)["declared_selected_reference_gate_passed"])

    def test_pair_cannot_override_failed_raw_gate(self):
        first = fixture(1); first["raw_diagnostics"].update(passed=False, divergences=3)
        result = sidecar.pair(first, fixture(2), 1)
        self.assertTrue(result["all_fixed_events_agree"])
        self.assertFalse(result["declared_selected_reference_gate_passed"])

    def test_pair_identity_and_probability_disagreement_rejected(self):
        for key, value in (("sampler_seed", 1), ("input_hash", "other"), ("case", "null")):
            second = fixture(2); second[key] = value
            self.assertFalse(sidecar.pair(fixture(1), second, 1)["declared_selected_reference_gate_passed"])
        second = fixture(2); second["fixed_functional"][0]["estimate"] = .53
        self.assertFalse(sidecar.pair(fixture(1), second, 1)["declared_selected_reference_gate_passed"])


if __name__ == "__main__":
    unittest.main()
