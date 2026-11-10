"""Test inventory and conservative feedback planning without running any R fit."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("planner", Path(__file__).with_name("plan_tests.py"))
planner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(planner)


class FeedbackPlanTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="msr-feedback-unit-")
        self.root = Path(self.temp.name).resolve()
        self.assertTrue(self.root.is_relative_to(Path(tempfile.gettempdir()).resolve()))
        self.addCleanup(self.temp.cleanup)
        (self.root / "tests/testthat").mkdir(parents=True)
        (self.root / "NAMESPACE").write_text("export(one)\nexport(two)\n")
        for name in ("core", "gdr-refitter", "gdr-tmb-refitter", "cv", "new"):
            (self.root / "tests/testthat" / ("test-" + name + ".R")).write_text("one()\n")
        self.config = dict(focused=dict(core="^core$", gdr="^gdr-", cv="^cv$"),
            smoke_files=["test-core.R"], source_profiles=[
                dict(pattern="^R/gdr", profiles=["gdr", "cv"]),
                dict(pattern="^src/|^tests/testthat/helper", profiles=[], full=True)])

    def test_smoke_is_explicit_partial_feedback(self):
        out = planner.plan(self.root, self.config, ["smoke"])
        self.assertEqual(out["selected_files"], ["test-core.R"])
        self.assertEqual(len(out["excluded_files"]), 4)
        self.assertFalse(out["release_eligible"])
        self.assertTrue(out["full_release_acceptance_required"])

    def test_full_cannot_drop_an_unclassified_test(self):
        out = planner.plan(self.root, self.config, full=True)
        self.assertEqual(out["catalog"]["unclassified"], ["test-new.R"])
        self.assertEqual(out["excluded_files"], [])
        self.assertTrue(out["all_files_selected"])

    def test_source_edges_union_and_exact_regex(self):
        out = planner.plan(self.root, self.config, changed=["R/gdr_tmb_refitter.R"])
        self.assertEqual(out["selected_files"], ["test-cv.R", "test-gdr-refitter.R", "test-gdr-tmb-refitter.R"])
        self.assertEqual(out["inferred_profiles"], ["cv", "gdr"])

    def test_unknown_core_and_changed_helpers_require_full_selection(self):
        for name in ("R/unknown.R", "src/file.cpp", "tests/testthat/helper-new.R"):
            out = planner.plan(self.root, self.config, changed=[name])
            self.assertTrue(out["all_files_selected"], name)

    def test_changed_test_keeps_related_integration_files(self):
        out = planner.plan(self.root, self.config, changed=["tests/testthat/test-gdr-refitter.R"])
        self.assertEqual(out["selected_files"], ["test-gdr-refitter.R", "test-gdr-tmb-refitter.R"])

    def test_field_cv_is_in_explicit_and_source_selected_cv_feedback(self):
        # Exercise the actual profile filters so new field tests cannot be
        # omitted merely because their name lacks the old cv- prefix.
        config = planner.read(Path(__file__).with_name("profiles.json"))
        config = dict(config, smoke_files=[])
        for name in ("field-cv", "cross-validation", "cv-comparison"):
            (self.root / "tests/testthat" / ("test-" + name + ".R")).write_text("one()\n")
        for profile in ("cv", "selection"):
            out = planner.plan(self.root, config, [profile])
            self.assertIn("test-field-cv.R", out["selected_files"])
            self.assertIn("test-cross-validation.R", out["selected_files"])
            self.assertFalse(out["all_files_selected"])
            self.assertFalse(out["release_eligible"])
        for source in ("R/cross_validate_field.R", "R/field_cv_gaussian.R", "R/field_cv_comparison.R"):
            out = planner.plan(self.root, config, changed=[source])
            self.assertIn("test-field-cv.R", out["selected_files"], source)
            self.assertIn("test-cross-validation.R", out["selected_files"], source)
            self.assertFalse(out["all_files_selected"])
            self.assertFalse(out["release_eligible"])

    def test_source_references_are_not_claimed_as_runtime_coverage(self):
        out = planner.inventory(self.root, self.config)
        self.assertEqual(out["exports_without_direct_source_reference"], ["two"])
        self.assertIn("not branch, runtime or assertion coverage", out["coverage_scope"])

    def test_function_valued_do_call_is_an_explicit_source_reference(self):
        (self.root / "tests/testthat/test-new.R").write_text("do.call(two, args)\n")
        out = planner.inventory(self.root, self.config)
        self.assertEqual(out["exports_without_direct_source_reference"], [])

    def test_registered_methods_do_not_claim_dispatch_coverage(self):
        (self.root / "NAMESPACE").write_text("export(one)\nS3method(stats::predict,multiScaleR)\n")
        (self.root / "tests/testthat/test-new.R").write_text("predict(object)\n")
        out = planner.inventory(self.root, self.config)["registered_S3_methods"][0]
        self.assertEqual(out["method"], "predict.multiScaleR")
        self.assertEqual(out["generic_source_references"], ["test-new.R"])
        self.assertFalse(out["runtime_dispatch_coverage_established"])

    def test_multiple_plans_reuse_one_frozen_inventory_without_rereading_sources(self):
        catalog = planner.inventory(self.root, self.config)
        with patch.object(planner, "inventory", side_effect=AssertionError("No repeated inventory")):
            out = planner.plan(self.root, self.config, ["gdr"], catalog=catalog)
        self.assertEqual(len(out["selected_files"]), 2)
        with self.assertRaises(ValueError):
            planner.plan(self.root, dict(self.config, focused={}), ["gdr"], catalog=catalog)

    def test_invalid_profiles_paths_and_stale_smoke_fail(self):
        with self.assertRaises(ValueError): planner.plan(self.root, self.config, ["missing"])
        for path in ("", "../R/foo.R", "C:/R/foo.R", "/R/foo.R", "R/../foo.R"):
            with self.assertRaises(ValueError): planner.plan(self.root, self.config, changed=[path])
        config = dict(self.config, smoke_files=["test-deleted.R"])
        with self.assertRaises(ValueError): planner.inventory(self.root, config)

    def test_records_are_exclusive_and_keep_existing_attempts(self):
        path = self.root / "plan.json"
        planner.write_new(path, dict(accepted=False))
        before = path.read_bytes()
        with self.assertRaises(ValueError): planner.write_new(path, dict(accepted=True))
        self.assertEqual(path.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
