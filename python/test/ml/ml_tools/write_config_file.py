import os, json

def write_config(models, model_dir, multi_mode=False, case_label=None):
    """
    Write one or multiple model configs to JSON.

    models: list of dicts, each dict must contain:
        - input_features
        - output_features
        - model_path
        - active_cells
        - apply_time
        - feature_engineering_input (optional)
        - scaling_input (optional)
        - feature_engineering_output (optional)
        - scaling_output (optional)
        - input_scaling_params (optional)
        - output_scaling_params (optional)

    multi_mode: if True, write all configs into a single JSON.
    case_label: used as filename for multi-mode JSON (required if multi_mode=True).
    """

    if multi_mode and not case_label:
        raise ValueError("case_label must be provided for multi_mode=True")

    os.makedirs(model_dir, exist_ok=True)
    all_configs = []

    for m in models:
        # Normalize input/output features
        input_features = m["input_features"]
        output_features = m["output_features"]
        if isinstance(input_features, str):
            input_features = [input_features]
        if isinstance(output_features, str):
            output_features = [output_features]

        n_input = len(input_features)
        n_output = len(output_features)

        # Optional fields (safe defaults)
        fe_input = m.get("feature_engineering_input") or [None]*n_input
        fe_output = m.get("feature_engineering_output") or [None]*n_output
        scaling_input = m.get("scaling_input") or [None]*n_input
        scaling_output = m.get("scaling_output") or [None]*n_output
        input_scaling_params = m.get("input_scaling_params") or [None]*n_input
        output_scaling_params = m.get("output_scaling_params") or [None]*n_output

        # Build input block
        input_block = {}
        for i, fname in enumerate(input_features):
            feat_eng = fe_input[i] if i < len(fe_input) else None
            scale = scaling_input[i] if i < len(scaling_input) else None
            feature_dict = {
                "feature_engineering": feat_eng.lower() if feat_eng and feat_eng.lower() != "none" else None,
                "scaling": scale.lower() if scale and scale.lower() != "none" else None,
            }
            feature_dict = {k: v for k, v in feature_dict.items() if v is not None}
            if i < len(input_scaling_params) and input_scaling_params[i] is not None:
                feature_dict["scaling_params"] = input_scaling_params[i]
            input_block[fname] = feature_dict

        # Build output block
        output_block = {}
        for i, fname in enumerate(output_features):
            feat_eng = fe_output[i] if i < len(fe_output) else None
            scale = scaling_output[i] if i < len(scaling_output) else None
            feature_dict = {
                "feature_engineering": feat_eng.lower() if feat_eng and feat_eng.lower() != "none" else None,
                "scaling": scale.lower() if scale and scale.lower() != "none" else None,
            }
            feature_dict = {k: v for k, v in feature_dict.items() if v is not None}
            if i < len(output_scaling_params) and output_scaling_params[i] is not None:
                feature_dict["scaling_params"] = output_scaling_params[i]
            output_block[fname] = feature_dict

        # Save active cells
        model_base = os.path.splitext(os.path.basename(m["model_path"]))[0]
        cells_file = os.path.join(model_dir, model_base + "_active_cells.txt")
        with open(cells_file, "w") as f:
            for cell in m["active_cells"]:
                f.write(f"{cell}\n")

        # Build config dict
        cfg = {
            "model_path": m["model_path"],
            "cell_indices_file": cells_file,
            "apply_times": [m["apply_time"]],
            "features": {
                "inputs": input_block,
                "outputs": output_block
            }
        }
        all_configs.append(cfg)

    # Determine JSON filename
    if multi_mode:
        if not case_label:
            raise ValueError("case_label is required for multi_mode")
        json_name = case_label + ".json"
    else:
        json_name = os.path.splitext(os.path.basename(models[0]["model_path"]))[0] + ".json"

    json_path = os.path.join(model_dir, json_name)
    
    # Always write a list of configs
    with open(json_path, "w") as f:
        json.dump(all_configs, f, indent=2)

    return json_path