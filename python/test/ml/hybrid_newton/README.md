# Hybrid Newton Test Suite

This test suite validates the integration of **Hybrid Newton (HyNE)** with machine learning models inside the `flow` reservoir simulator.
It ensures that adding ML-based approximations does not degrade solver performance and, in some cases, improves convergence.

---

## Dependencies

In a similar way as for opm-common ML module, these tests require **TensorFlow** to be installed. 
The CPU version is enough for this test and allows faster installation.
```bash
pip install tensorflow-cpu
```

---

## Test Framework

* Implemented with Python’s built-in [`unittest`](https://docs.python.org/3/library/unittest.html).
* Each test case defines a specific **Hybrid Newton configuration** via JSON.
* Simulations are run using `flow`(therefore assuming it is installed):

  * **Baseline** → Standard Newton solver (no ML).
  * **Hybrid Newton** → Newton solver with ML-based updates injected.

---

## Baseline Simulation

Before running hybrid cases:

* Runs a standard `flow` simulation with the given deck.
* Extracts:

  * **Baseline Newton iterations** (`baseline_iters`) from simulation logs.
  * **Reservoir state variables** (`PRESSURE`, `SGAS`, `SWAT`, etc.) from UNRST files.
* Stores `n_cells` and `times` for consistency across tests.

---

## Hybrid Newton Simulation

For each test case:

1. A **minimal Keras model** is built that maps input → output features.

   * Weights are initialized to copy output variables directly.
   * Models are exported using `kerasify.export_model`.
2. A **JSON config** is generated with:

   * Input/output features.
   * Scaling and feature engineering options.
   * Apply time and active cells.
3. `flow` is run with:

   ```bash
   --use-hy-ne=true --hy-ne-config-file=config.json
   ```
4. Hybrid Newton iteration counts are extracted and compared to baseline.

---

## Assertions

Each test validates:

* **Convergence parity**
  Hybrid and baseline must produce the same number of timesteps.
* **Iteration bounds**
  `hybrid_iters[t] <= baseline_iters[t]` for all timesteps.
* **Zero Newton cases**
  Some cases require the first two timesteps to have `0` Newton iterations.
* **Config integrity**
  Config JSONs must exist and be valid.

---

## Test Categories

| Category                    | Description                                           |
| --------------------------- | ----------------------------------------------------- |
| `ABSOLUTE_CASES`            | Tests absolute features.                              |
| `RELATIVE_CASES`            | Tests relative features.                              |
| `FEATURE_ENGINEERING_CASES` | Tests features engineered inputs/outputs (e.g., log). |
| `SCALING_CASES`             | Applies scaling strategies (e.g., min-max, standard). |
| `MULTI_MODEL_CASES`         | Combines multiple models into one simulation.         |
| `ZERO_NEWTON_CASES`         | Use all features, resulting in 0 Newton iterations.   |
| `ALL_CASES`                 | Runs all cases together as a stress test.             |

---

## Utilities

The suite uses helper functions from `utils.py`:

* `collect_input_features` → Prepares input arrays from UNRST data.
* `compute_output_vars` → Builds expected outputs from reservoir state.
* `write_config` → Writes the HyNE JSON configuration.
* `export_model` → Saves Keras models in a format usable by `flow`.
* `extract_newtit_from_file` → Extracts iteration counts from solver logs.
* `extract_unrst_variables` → Reads simulation outputs from UNRST files.

---
