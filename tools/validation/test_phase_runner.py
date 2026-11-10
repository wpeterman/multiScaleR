"""Non-numerical contract tests for the reusable phase runner."""
import copy
import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

spec = importlib.util.spec_from_file_location("phase_runner", Path(__file__).with_name("phase_runner.py"))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)


class PhasePlanTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="msr-runner-unit-")
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)
        self.output = self.root / "tools/benchmarks/example"
        self.output.mkdir(parents=True)
        self.control = self.root / "control.json"
        self.control.write_text(json.dumps(dict(long_launch_until="2026-11-03T09:00:00Z",
            launch_until="2026-11-03T09:30:00Z", worker_stop_at="2026-11-03T09:54:30Z",
            window_end="2026-11-03T10:00:00Z", threads=1, maximum_owned_interpreters=3,
            maximum_owned_fit_workers=2, RAM_reserve_bytes=8 * 2**30,
            disk_runtime_reserve_bytes=5 * 2**30, disk_launch_reserve_bytes=20 * 2**30,
            foreign_jobs_untouched=True, timing_comparison_eligible=False)))
        self.exe = self.root / "Rscript.exe"
        self.exe.write_bytes(b"test fixture, never executed")
        self.input = self.root / "input.R"
        self.input.write_text("stop('unit fixture never launches R')")
        paths = [self.exe, self.input, self.control, Path(runner.__file__),
                 Path(runner.__file__).with_name("run_phase.ps1"),
                 Path(runner.__file__).with_name("resource_control.ps1")]
        identities = [dict(root=str(q.parent), files={q.name: runner.digest(q)}) for q in paths]
        self.plan = dict(schema_version=1, repo=str(self.root), output=str(self.output),
            candidate=str(self.root / "candidate"), branch="inference-dev", head="a" * 40,
            control_file=str(self.control), identities=identities, generated_files=[],
            phases=[dict(name="focus01", executable=str(self.exe), arguments=[str(self.input)],
                         argument_line=subprocess.list2cmdline([str(self.input)]), working_directory=str(self.root))])

    def test_valid_plan_and_input_identities(self):
        runner.validate_plan(self.plan)
        runner.verify_inputs(self.plan)

    def test_unknown_schema_and_branch_are_rejected(self):
        for name, value in (("schema_version", 2), ("branch", "master"), ("head", "short")):
            plan = copy.deepcopy(self.plan); plan[name] = value
            with self.assertRaises(ValueError):
                runner.validate_plan(plan)

    def test_output_escape_and_unsafe_phase_are_rejected(self):
        plan = copy.deepcopy(self.plan); plan["output"] = str(self.root)
        with self.assertRaises(ValueError): runner.validate_plan(plan)
        for name in ("../escape", "focus;echo", "FOCUS", "a/b"):
            plan = copy.deepcopy(self.plan); plan["phases"][0]["name"] = name
            with self.assertRaises(ValueError): runner.validate_plan(plan)

    def test_resource_caps_and_deadline_order_cannot_be_relaxed(self):
        original = json.loads(self.control.read_text())
        for name, value in (("threads", 2), ("maximum_owned_interpreters", 4),
                            ("disk_runtime_reserve_bytes", 1), ("timing_comparison_eligible", True),
                            ("launch_until", "2026-11-03T08:00:00Z")):
            changed = dict(original); changed[name] = value
            self.control.write_text(json.dumps(changed))
            with self.assertRaises(ValueError): runner.validate_plan(self.plan)

    def test_argument_quoting_and_executable_admission(self):
        plan = copy.deepcopy(self.plan); plan["phases"][0]["argument_line"] = "different"
        with self.assertRaises(ValueError): runner.validate_plan(plan)
        plan = copy.deepcopy(self.plan); plan["phases"][0]["executable"] = "cmd.exe"
        with self.assertRaises(ValueError): runner.validate_plan(plan)

    def test_source_and_producer_drift_are_rejected(self):
        self.input.write_text("changed")
        with self.assertRaises(ValueError): runner.verify_inputs(self.plan)
        plan = copy.deepcopy(self.plan); plan["identities"] = plan["identities"][:-1]
        with self.assertRaises(ValueError): runner.validate_plan(plan)

    def test_identity_path_escape_is_rejected(self):
        group = dict(root=str(self.root), files={"../escape": "a" * 64})
        with self.assertRaises(ValueError): list(runner.identity_map(group))

    def test_only_declared_candidate_docs_may_change_after_generation(self):
        doc = self.root / "candidate/man/topic.Rd"
        doc.parent.mkdir(parents=True); doc.write_text("old")
        self.plan["identities"].append(dict(root=str(doc.parent), files={doc.name: runner.digest(doc)}))
        self.plan["generated_files"] = [str(doc)]
        runner.validate_plan(self.plan)
        doc.write_text("generated")
        runner.verify_inputs(self.plan, after=True)
        with self.assertRaises(ValueError): runner.verify_inputs(self.plan)
        self.plan["generated_files"] = [str(self.input)]
        with self.assertRaises(ValueError): runner.validate_plan(self.plan)

    def test_atomic_records_preserve_previous_attempt(self):
        record = self.output / "record.json"
        runner.atomic_json(record, dict(value=1))
        with self.assertRaises(ValueError): runner.atomic_json(record, dict(value=2))
        self.assertEqual(runner.read_json(record), dict(value=1))

    def namespace_plan(self):
        temporary = tempfile.TemporaryDirectory(prefix="msr-doc-candidate-")
        self.addCleanup(temporary.cleanup)
        candidate = Path(temporary.name).resolve()
        namespace = candidate / "NAMESPACE"; namespace.write_text("old")
        input_path = self.root / "document-input.json"
        input_path.write_text(json.dumps(dict(package=str(candidate))))
        script = Path(runner.__file__).with_name("document.R").resolve()
        plan = copy.deepcopy(self.plan); plan["candidate"] = str(candidate)
        plan["generated_files"] = [str(namespace)]
        arguments = ["--vanilla", str(script), str(input_path)]
        plan["phases"][0].update(name="doc01", arguments=arguments,
            argument_line=subprocess.list2cmdline(arguments))
        plan["identities"] += [dict(root=str(q.parent), files={q.name: runner.digest(q)})
                               for q in (namespace, input_path, script)]
        return plan, namespace, input_path

    def test_isolated_roxygen_namespace_has_bound_input_and_starting_identity(self):
        plan, namespace, _ = self.namespace_plan()
        runner.validate_plan(plan); runner.verify_inputs(plan)
        namespace.write_text("generated")
        runner.verify_inputs(plan, after=True)
        with self.assertRaises(ValueError): runner.verify_inputs(plan)
        self.input.write_text("source changed")
        with self.assertRaises(ValueError): runner.verify_inputs(plan, after=True)

    def test_namespace_cannot_relax_live_source_or_phase_guards(self):
        plan, namespace, input_path = self.namespace_plan()
        for field, value in (("name", "focus01"), ("executable", "Rcmd.exe")):
            changed = copy.deepcopy(plan); changed["phases"][0][field] = value
            with self.assertRaises(ValueError): runner.validate_plan(changed)
        changed = copy.deepcopy(plan); changed["candidate"] = str(self.root)
        with self.assertRaises(ValueError): runner.validate_plan(changed)
        changed = copy.deepcopy(plan); changed["phases"] *= 2
        with self.assertRaises(ValueError): runner.validate_plan(changed)
        input_path.write_text(json.dumps(dict(package=str(self.root))))
        with self.assertRaises(ValueError): runner.validate_plan(plan)

    def test_namespace_requires_the_frozen_roxygen_script(self):
        plan, _, _ = self.namespace_plan()
        changed = copy.deepcopy(plan); changed["phases"][0]["arguments"][1] = str(self.input)
        changed["phases"][0]["argument_line"] = subprocess.list2cmdline(changed["phases"][0]["arguments"])
        with self.assertRaises(ValueError): runner.validate_plan(changed)
        changed = copy.deepcopy(plan); changed["identities"] = changed["identities"][:-3] + changed["identities"][-2:]
        with self.assertRaises(ValueError): runner.validate_plan(changed)

    def test_replay_never_reaches_subprocess(self):
        plan_path = self.root / "plan.json"
        plan_path.write_text(json.dumps(self.plan))
        (self.output / "focus01-request.json").write_text("{}")
        with patch.object(runner.subprocess, "run") as execute:
            with self.assertRaises(ValueError): runner.run(plan_path, "focus01")
            execute.assert_not_called()

    def test_input_names_are_not_prior_attempts(self):
        plan_path = self.root / "plan.json"
        plan_path.write_text(json.dumps(self.plan))
        (self.output / "focus01-input.json").write_text("{}")
        with patch.object(runner, "verify_git"), patch.object(runner.subprocess, "run") as execute:
            execute.return_value.returncode = 0
            runner.run(plan_path, "focus01")
            execute.assert_called_once()
        self.assertTrue((self.output / "focus01-verified.json").exists())

    def test_partial_reservation_blocks_replay(self):
        plan_path = self.root / "plan.json"
        plan_path.write_text(json.dumps(self.plan))
        (self.output / "focus01-request.json.tmp").write_text("unfinished")
        with patch.object(runner.subprocess, "run") as execute:
            with self.assertRaises(ValueError): runner.run(plan_path, "focus01")
            execute.assert_not_called()

    def test_new_plans_freeze_runtime_bytes_and_reject_overwriting(self):
        import shutil
        maker_spec = importlib.util.spec_from_file_location("test_plan_maker", Path(__file__).with_name("make_plan.py"))
        maker = importlib.util.module_from_spec(maker_spec); maker_spec.loader.exec_module(maker)
        code = self.root / "tool-code"; code.mkdir()
        for name in maker.RUNTIME_FILES:
            shutil.copyfile(Path(__file__).with_name(name), code / name)
        def git(args, **kwargs):
            if args[1:3] == ["branch", "--show-current"]: return b"inference-dev\n"
            if args[1:3] == ["rev-parse", "HEAD"]: return ("a" * 40).encode()
            if args[1:3] == ["diff", "--cached"]: return b""
            raise AssertionError(args)
        arguments = dict(repo=self.root, output=self.output, control=self.control,
            identities=self.plan["identities"], phases=self.plan["phases"], library_paths=[],
            candidate=self.plan["candidate"])
        plan_path = self.output / "immutable-plan.json"
        with patch.object(maker, "__file__", str(code / "make_plan.py")), patch.object(subprocess, "check_output", side_effect=git):
            frozen = maker.make_plan(plan_path, **arguments)
            original = frozen.read_bytes()
            (code / "phase_runner.py").write_bytes(original + b"\n# later tooling edit\n")
            self.assertEqual(frozen.read_bytes(), original)
            second = maker.make_plan(self.output / "second-plan.json", **arguments)
            self.assertNotEqual(second.parent, frozen.parent)
            with self.assertRaises(ValueError): maker.make_plan(plan_path, **arguments)


if __name__ == "__main__":
    unittest.main(verbosity=2)
