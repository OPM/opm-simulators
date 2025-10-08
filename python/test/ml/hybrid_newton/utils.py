import json
import os
import numpy as np
from pathlib import Path

from opm.io.ecl import ERst, ESmry

# Avoid tensorflow dependency but maintain compability with kerasify export_model method
class Dense:
    def __init__(self, input_dim, output_dim, activation='linear'):
        self.weights = np.zeros((input_dim, output_dim), dtype=np.float32)
        self.biases = np.zeros((output_dim,), dtype=np.float32)
        self.activation = activation

    def set_weights(self, w_b):
        self.weights = w_b[0].astype(np.float32)
        self.biases = w_b[1].astype(np.float32)

    def get_weights(self):
        return [self.weights, self.biases]

    def get_config(self):
        return {'activation': self.activation}

class Sequential:
    def __init__(self, layers):
        self.layers = layers


def extract_unrst_variables(data_dir, deck_file):
    """
    Extracts fluid state variables and timesteps from the .UNRST file
    corresponding to a deck.

    Parameters
    ----------
    data_dir : Path or str
        Directory containing the UNRST file.
    deck_file : str
        Name of the deck file (e.g., SPE1CASE1.DATA).

    Returns
    -------
    unrst_vars : dict
        Dictionary with keys "PRESSURE", "SWAT", "SOIL", "SGAS", "RS"
        Each value is a list of np.ndarray for each report step.
    dts : list
        List of timestep durations (approximated from report steps if needed).
    """
    deck_stem = Path(deck_file).stem
    unrst_file = str(Path(data_dir) / f"{deck_stem}.UNRST")  # convert to string
    unsmry_file = str(Path(data_dir) / f"{deck_stem}.SMSPEC")
    rst = ERst(unrst_file)
    unsmry = ESmry(unsmry_file)
    # Load all report steps
    report_steps = sorted(rst.report_steps)
    for step in report_steps:
        rst.load_report_step(step)

    variables_of_interest = ["PRESSURE", "SWAT", "SOIL", "SGAS", "RS"]
    unrst_vars = {var: [] for var in variables_of_interest}

    for step in report_steps:
        for var in variables_of_interest:
            if (var, step) in rst:
                unrst_vars[var].append(np.array(rst[var, step], copy=True))
            else:
                # Fill with zeros if variable not in UNRST for this step
                n_cells = len(rst["PRESSURE", step])  # assume same size as PRESSURE
                unrst_vars[var].append(np.zeros(n_cells, dtype=np.float32))

    # Compute approximate timesteps (differences between report steps)
    dates = unsmry.dates()
    times = [0.] + [(d - unsmry.start_date).total_seconds() for d in dates]

    return unrst_vars, times


def write_config(models, model_dir, multi_mode=False, case_label=None):
    """
    Write one or multiple model configs to JSON in OPM-supported format.
    """

    model_dir = Path(model_dir)
    if multi_mode and not case_label:
        raise ValueError("case_label must be provided for multi_mode")

    os.makedirs(model_dir, exist_ok=True)
    all_configs = {}

    for idx, m in enumerate(models):
        # Ensure input/output features are lists
        input_features = m["input_features"] if isinstance(m["input_features"], list) else [m["input_features"]]
        output_features = m["output_features"] if isinstance(m["output_features"], list) else [m["output_features"]]

        fe_input = m.get("feature_engineering_input") or [None]*len(input_features)
        fe_output = m.get("feature_engineering_output") or [None]*len(output_features)
        scaling_params_input = m.get("input_scaling_params") or [None]*len(input_features)
        scaling_params_output = m.get("output_scaling_params") or [None]*len(output_features)

        # Build input/output blocks
        input_block = {}
        for i, f in enumerate(input_features):
            d = {}
            if fe_input[i] and fe_input[i].lower() != "none":
                d["feature_engineering"] = fe_input[i].lower()
            if scaling_params_input[i] is not None:
                # convert numpy float32 to native float if needed
                d["scaling_params"] = {k: float(v) for k, v in scaling_params_input[i].items()}
            input_block[f] = d

        output_block = {}
        for i, f in enumerate(output_features):
            d = {}
            if fe_output[i] and fe_output[i].lower() != "none":
                d["feature_engineering"] = fe_output[i].lower()
            if scaling_params_output[i] is not None:
                d["scaling_params"] = {k: float(v) for k, v in scaling_params_output[i].items()}
            output_block[f] = d

        # Active cells
        active_cells = m.get("active_cells")
        if active_cells is None:
            n_cells = getattr(m, "n_cells", None)
            if n_cells is None:
                raise ValueError("Cannot determine active_cells; provide n_cells or active_cells")
            active_cells = list(range(n_cells))

        cells_file = model_dir / f"{os.path.splitext(os.path.basename(m['model_path']))[0]}_active_cells.txt"
        with open(cells_file, "w") as f:
            for cell in active_cells:
                f.write(f"{cell}\n")

        # Apply time(s)
        apply_time = m["apply_time"]
        if isinstance(apply_time, float):
            apply_times_dict = {"0": float(apply_time)}
        elif isinstance(apply_time, list):
            apply_times_dict = {str(i): float(t) for i, t in enumerate(apply_time)}
        else:
            raise ValueError("apply_time must be float or list of floats")

        model_name_rel = f"models/{os.path.basename(m['model_path'])}"
        cells_file_rel = f"models/{cells_file.name}"

        all_configs[str(idx)] = {
            "model_path": model_name_rel,
            "cell_indices_file": cells_file_rel,
            "apply_times": apply_times_dict,
            "features": {"inputs": input_block, "outputs": output_block}
        }

    # JSON filename
    if multi_mode:
        json_name = case_label + ".json"
    else:
        json_name = f"{os.path.splitext(os.path.basename(models[0]['model_path']))[0]}.json"

    json_path = model_dir / json_name
    with open(json_path, "w") as f:
        json.dump(all_configs, f, indent=2)

    return json_path


def seconds_to_days(sec):
    return sec / 86400.0


def apply_feature_engineering(data, method=None):
    if not method or method.lower() == "none":
        return data
    if method.lower() == "log":
        return np.log(data)
    if method.lower() == "log10":
        return np.log10(data)
    if method.lower() == "log1p":
        return np.log1p(data)
    raise ValueError(f"Unsupported feature engineering method: {method}")


def apply_scaling_with_params(data, method):
    data = np.array(data)
    if method == "standard":
        mean, std = np.mean(data), np.std(data)
        return ((data - mean) / std if std != 0 else data - mean), {"mean": mean, "std": std}
    if method == "minmax":
        mn, mx = np.min(data), np.max(data)
        return ((data - mn) / (mx - mn) if mx != mn else data - mn), {"min": mn, "max": mx}
    return data, None


def get_feature_data(name, timestep, unrst, rd_init, times=None, feat_eng=None):
    if name == "PERMX":
        val = np.array(rd_init["PERMX"])
    elif name == "SOIL":
        val = 1.0 - np.array(unrst["SGAS"][timestep]) - np.array(unrst["SWAT"][timestep])
    elif name == "TIMESTEP":
        if times is None:
            raise ValueError("times must be provided for TIMESTEP")
        val = times[timestep + 1] - times[timestep]
    else:
        val = np.array(unrst[name][timestep])
    return apply_feature_engineering(val, feat_eng)


def collect_input_features(input_features, start_idx, unrst, rd_init, times=None, feat_eng=None, scaling=None):
    input_data = []
    scaling_params = []

    for i, f in enumerate(input_features):
        fe = feat_eng[i] if feat_eng else None
        sc = scaling[i] if scaling else None
        val = get_feature_data(f, start_idx, unrst, rd_init, times=times, feat_eng=fe)
        val, params = apply_scaling_with_params(val, sc) if sc else (val, None)
        input_data.append(np.array(val).flatten())
        scaling_params.append(params)

    return np.concatenate(input_data), scaling_params


def compute_output_vars(output_features, start_idx, end_idx, unrst, rd_init, feature_engineering=None, scaling=None):
    output_data = []
    scaling_params = []

    for i, f in enumerate(output_features):
        fe = feature_engineering[i] if feature_engineering else None
        sc = scaling[i] if scaling else None

        if f.startswith("DELTA_"):
            base = f[6:]
            val = get_feature_data(base, end_idx, unrst, rd_init) - get_feature_data(base, start_idx, unrst, rd_init)
        else:
            val = get_feature_data(f, end_idx, unrst, rd_init)

        val = apply_feature_engineering(val, fe)
        val, params = apply_scaling_with_params(val, sc) if sc else (val, None)
        output_data.append(np.array(val).flatten())
        scaling_params.append(params)

    return np.stack(output_data), scaling_params


def extract_newtit_from_file(filename):
    """
    Parse an INFOSTEP file and return the values of the 'NewtIt' column.

    Parameters
    ----------
    filename : str or Path
        Path to the file (with or without .INFOSTEP extension).

    Returns
    -------
    list[int]
        List of NewtIt values for each data row.
    """
    path = Path(filename)
    if path.suffix.lower() != ".infostep":
        path = path.with_suffix(".INFOSTEP")
    if not path.exists():
        raise FileNotFoundError(f"INFOSTEP file not found: {path}")

    with path.open("r") as f:
        lines = [line.strip() for line in f if line.strip()]

    # Locate header with "NewtIt"
    header_idx = next((i for i, line in enumerate(lines) if "NewtIt" in line), None)
    if header_idx is None:
        raise ValueError(f"'NewtIt' column not found in {path}")

    header = lines[header_idx].split()
    idx = header.index("NewtIt")

    newt_values = []
    for row in lines[header_idx + 1 :]:
        cols = row.split()
        if len(cols) > idx:
            try:
                newt_values.append(int(cols[idx]))
            except ValueError:
                continue  # skip malformed rows

    return newt_values

# -------------------------
# Absolute (single timestep predictions)
# -------------------------
ABSOLUTE_CASES = [
    {
        "label": "single_RS",
        "mode": "single",
        "start_idx": 0, "end_idx": 1,
        "input_features": ["TIMESTEP"],
        "output_features": ["RS"],
    },
    {
        "label": "single_PRESSURE",
        "mode": "single",
        "start_idx": 0, "end_idx": 1,
        "input_features": ["PRESSURE"],
        "output_features": ["PRESSURE"],
    },
    {
        "label": "single_PERMX_to_SGAS",
        "mode": "single",
        "start_idx": 0, "end_idx": 1,
        "input_features": ["PERMX"],
        "output_features": ["SGAS"],
    },
    {
        "label": "single_SGAS_to_SWAT",
        "mode": "single",
        "start_idx": 0, "end_idx": 1,
        "input_features": ["SGAS"],
        "output_features": ["SWAT"],
    },
]

# -------------------------
# Relative / delta predictions
# -------------------------
RELATIVE_CASES = [
    {"label": "delta_RS", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["RS"], "output_features": ["DELTA_RS"]},
    {"label": "delta_SGAS", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["SGAS"], "output_features": ["DELTA_SGAS"]},
    {"label": "delta_SWAT", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["SWAT"], "output_features": ["DELTA_SWAT"]},
    {"label": "delta_PRESSURE", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PRESSURE"], "output_features": ["DELTA_PRESSURE"]},
    {"label": "delta_PERMX_to_RS", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PERMX"], "output_features": ["DELTA_RS"]},
    {"label": "delta_PERMX_to_PRESSURE", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PERMX"], "output_features": ["DELTA_PRESSURE"]},
]

# -------------------------
# Feature engineering cases
# -------------------------
FEATURE_ENGINEERING_CASES = [
    {"label": "feat_eng_RS_to_PRESSURE", "mode": "single", "start_idx": 0, "end_idx": 1,
     "input_features": ["RS"], "feature_engineering_input": ["log"],
     "output_features": ["PRESSURE"], "feature_engineering_output": ["log1p"]},
    {"label": "feat_eng_PRESSURE", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["SGAS"], "feature_engineering_input": ["log1p"],
     "output_features": ["PRESSURE"], "feature_engineering_output": ["log10"]},
    {"label": "feat_eng_PERMX_to_PRESSURE", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PERMX"], "feature_engineering_input": ["log"],
     "output_features": ["PRESSURE"], "feature_engineering_output": ["log"]},
]

# -------------------------
# Scaling cases
# -------------------------
SCALING_CASES = [
    {"label": "scaling_RS_minmax", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["RS"], "scaling_input": ["minmax"],
     "output_features": ["RS"], "scaling_output": ["minmax"]},
    {"label": "scaling_PRESSURE_standard", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PRESSURE"], "scaling_input": ["standard"],
     "output_features": ["PRESSURE"], "scaling_output": ["standard"]},
    {"label": "scaling_SGAS_SWAT_standard", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["SGAS", "SWAT"], "scaling_input": ["standard", "none"], 
     "output_features": ["SGAS", "SWAT"], "scaling_output": ["none", "standard"]},
    {"label": "scaling_PRESSURE_SWAT_minmax", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PRESSURE", "SWAT"], "scaling_input": ["minmax", "none"],
     "output_features": ["SWAT", "PRESSURE"], "scaling_output": ["none", "minmax"]},
    {"label": "scaling_PERMX_TIMESTEP_standard", "mode": "single", "start_idx": 1, "end_idx": 2,
     "input_features": ["PRESSURE", "TIMESTEP"], "scaling_input": ["minmax", "none"],
     "output_features": ["SGAS", "SWAT"]},
]

# -------------------------
# Multi-model cases
# -------------------------
MULTI_MODEL_CASES = [
    {"label": "multi_same_index", "mode": "multi",
     "models": [
         {"input_features": ["PERMX"], "output_features": ["RS"], "start_idx": 1, "end_idx": 2},
         {"input_features": ["SGAS"], "output_features": ["SGAS"], "start_idx": 1, "end_idx": 2},
     ]},
    {"label": "multi_diff_index", "mode": "multi",
     "models": [
         {"input_features": ["PERMX"], "output_features": ["RS"], "start_idx": 1, "end_idx": 2},
         {"input_features": ["PRESSURE", "TIMESTEP", "PERMX"], "output_features": ["PRESSURE"], "start_idx": 8, "end_idx": 9},
         {"input_features": ["SGAS"], "output_features": ["SGAS"], "start_idx": 10, "end_idx": 11},
     ]},
]

ZERO_NEWTON_CASES = [
    {
        "label": "all_features",
        "mode": "multi",
        "models": [
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 0, "end_idx": 1},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 1, "end_idx": 2},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 2, "end_idx": 3},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 3, "end_idx": 4},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 4, "end_idx": 5},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 5, "end_idx": 6},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 6, "end_idx": 7},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 7, "end_idx": 8},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 8, "end_idx": 9},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 9, "end_idx": 10},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 10, "end_idx": 11},
            {"input_features": ["PERMX"], "output_features": ["RS", "PRESSURE", "SGAS", "SOIL"], "start_idx": 11, "end_idx": 12},
        ],

    },
]

ALL_CASES = [
    {
        "label": "all_cases",
        "mode": "multi",
        "models": [
            {"input_features": ["PERMX", "TIMESTEP", "RS", "SOIL", "SWAT", "PRESSURE"],
             "scaling_input": ["none", "none", "standard", "standard", "standard","minmax"],
             "feature_engineering_input": ["log", "log10", "log1p", "log1p", "log1p", "log10"],
             "output_features": ["RS", "PRESSURE", "DELTA_SGAS", "DELTA_SOIL"],
             "scaling_output": ["standard", "minmax", "standard", "standard"],
             "feature_engineering_output": ["log", "log10", "log1p", "log1p"],
             "start_idx": 0, "end_idx": 1},
             {"input_features": ["PERMX", "TIMESTEP", "RS", "SOIL", "SWAT", "PRESSURE"],
             "scaling_input": ["none", "none", "standard", "standard", "standard","minmax"],
             "feature_engineering_input": ["log", "log10", "log1p", "log1p", "log1p", "log10"],
             "output_features": ["DELTA_RS", "DELTA_PRESSURE", "DELTA_SGAS", "DELTA_SOIL"],
             "scaling_output": ["standard", "minmax", "standard", "standard"],
             "feature_engineering_output": ["none", "none", "log1p", "log1p"],
             "start_idx": 1, "end_idx": 2},
        ],

    },
]
