import os
import unittest
from pathlib import Path
import numpy as np
import subprocess

from utils import Dense, Sequential
from utils import collect_input_features, compute_output_vars, write_config, extract_newtit_from_file, extract_unrst_variables
from utils import ABSOLUTE_CASES, RELATIVE_CASES, FEATURE_ENGINEERING_CASES, SCALING_CASES, MULTI_MODEL_CASES, ZERO_NEWTON_CASES, ALL_CASES

from opm.ml.ml_tools.kerasify import export_model

class TestHybridNewton(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        """Run baseline simulation via subprocess and collect timesteps & fluid state."""
        test_dir = Path(__file__).parent
        cls.data_dir = test_dir.parent.parent.parent / "test_data/SPE1CASE1b"
        cls.deck_file = "SPE1CASE1.DATA"

        cls.unrst = {}
        cls.rd_init = {"PERMX": []}
        cls.times = []

        baseline_cmd = [
            os.environ.get("FLOW_BINARY", default="flow"),
            str(cls.data_dir / cls.deck_file),
            "--output-extra-convergence-info=steps,iterations",
            "--newton-min-iterations=0",
            "--output-dir=.",
            "--full-time-step-initially=true",
        ]

        print("Running baseline simulation... ")
        subprocess.run(baseline_cmd, stdout=subprocess.DEVNULL, check=True)

        # Extract Newton iterations
        cls.baseline_iters = extract_newtit_from_file(Path('.') / Path(cls.deck_file).stem)
        print(f"Baseline Newton iterations: {cls.baseline_iters}")

        # Extract fluid state from UNRST for all timesteps
        cls.unrst, cls.times = extract_unrst_variables(".", cls.deck_file)
        cls.rd_init["PERMX"] = [10**2] * len(cls.unrst["PRESSURE"][0])
        cls.n_cells = len(cls.unrst["PRESSURE"][0])

    def _run_cases(self, cases):
        models_dir = Path('.') / "models"
        os.makedirs(models_dir, exist_ok=True)

        for case in cases:
            with self.subTest(case=case["label"]):
                json_path = self._build_case_models(case, models_dir)
                self.assertTrue(json_path and os.path.isfile(json_path),
                                f"Config JSON was not created for case {case['label']}")

                # Run hybrid Newton test
                hybrid_iters = self._run_hybrid_newton_test(json_path)
        return hybrid_iters

    def _build_case_models(self, case, models_dir):
        models_to_write = []

        if case.get("mode", "single") == "single":
            model_entries = [case]
            multi_mode = False
        else:
            model_entries = case.get("models", [])
            multi_mode = True

        for model_case in model_entries:
            start_idx = model_case.get("start_idx", case.get("start_idx"))
            end_idx = model_case.get("end_idx", case.get("end_idx"))
            input_features = model_case["input_features"]
            output_features = model_case["output_features"]
            fe_input = model_case.get("feature_engineering_input")
            fe_output = model_case.get("feature_engineering_output")
            scaling_input = model_case.get("scaling_input")
            scaling_output = model_case.get("scaling_output")

            input_flat, input_scaling_params = collect_input_features(
                input_features, start_idx, self.unrst, self.rd_init, self.times,
                feat_eng=fe_input, scaling=scaling_input
            )
            output_flat, output_scaling_params = compute_output_vars(
                output_features, start_idx, end_idx, self.unrst, self.rd_init,
                feature_engineering=fe_output, scaling=scaling_output,
            )

            input_flat = input_flat.reshape(-1)
            output_flat = output_flat.reshape(-1)

            model = Sequential([
                Dense(input_dim=input_flat.shape[0], output_dim=output_flat.shape[0])
            ])
            model.layers[0].set_weights([np.zeros((input_flat.shape[0], output_flat.shape[0])), output_flat])
            apply_time = self.times[start_idx]

            active_cells = list(range(self.n_cells))

            models_dir.mkdir(exist_ok=True)
            model_name = "_".join(input_features) + "__TO__" + "_".join(output_features) + f"_{start_idx}_{end_idx}.model"
            model_path = models_dir / model_name
            export_model(model, model_path)

            models_to_write.append({
                "input_features": input_features,
                "output_features": output_features,
                "model_path": model_path,
                "active_cells": active_cells,
                "apply_time": apply_time,
                "feature_engineering_input": fe_input,
                "feature_engineering_output": fe_output,
                "scaling_input": scaling_input,
                "scaling_output": scaling_output,
                "input_scaling_params": input_scaling_params,
                "output_scaling_params": output_scaling_params
            })

        if not models_to_write:
            print(f"No valid models for {case['label']}")
            return None

        return write_config(models_to_write, models_dir, multi_mode=multi_mode, case_label=case['label'])

    def _run_hybrid_newton_test(self, json_path):
        """Run hybrid Newton via subprocess and extract Newton iterations."""
        cmd = [
            os.environ.get("FLOW_BINARY", default="flow"),
            str(self.data_dir / self.deck_file),
            f"--hy-ne-config-file={json_path}",
            "--use-hy-ne=true",
            "--output-extra-convergence-info=steps,iterations",
            "--newton-min-iterations=0",
            "--output-dir=.",
            "--full-time-step-initially=true"
        ]

        print(f"Running hybrid Newton simulation with config {json_path}...")
        subprocess.run(cmd,
                       stdout=subprocess.DEVNULL,
                       check=True)

        hybrid_iters = extract_newtit_from_file(Path('.') / Path(self.deck_file).stem)
        return hybrid_iters

    def test_absolute_cases(self):
        self._run_cases(ABSOLUTE_CASES)

    def test_relative_cases(self):
        self._run_cases(RELATIVE_CASES)

    def test_feature_engineering_cases(self):
        self._run_cases(FEATURE_ENGINEERING_CASES)

    def test_scaling_cases(self):
        self._run_cases(SCALING_CASES)

    def test_multi_model_cases(self):
        self._run_cases(MULTI_MODEL_CASES)

    def test_zero_newton_cases(self):
        hybrid_iters = self._run_cases(ZERO_NEWTON_CASES)
        print(f"Hybrid Newton iterations all_features: {hybrid_iters}")

    def test_all_cases(self):
        hybrid_iters = self._run_cases(ALL_CASES)
        print(f"Hybrid Newton iterations all_cases: {hybrid_iters}")


if __name__ == "__main__":
    unittest.main()
