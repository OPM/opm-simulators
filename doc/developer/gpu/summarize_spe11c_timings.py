#!/usr/bin/env python3
"""Summarize three paired full-case timing runs and their transfer counters.

Usage: python3 summarize_spe11c_timings.py BIN_DIRECTORY OUTPUT_PREFIX
"""

import argparse
import json
from pathlib import Path
import re
import statistics

from compare_spe11c import run_statistics, unique_file


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("prefix")
    args = parser.parse_args()
    result = {}
    for mode in ("disabled", "enabled"):
        runs = []
        for repetition in (1, 2, 3):
            directory = args.directory / f"{args.prefix}_timing_{mode}_{repetition}"
            metrics, _ = run_statistics(directory)
            log = unique_file(directory, "PRT").read_text()
            counters = re.findall(r"\[GPU Newton transfers\] (.*)", log)
            if len(counters) != 1:
                raise ValueError(f"Expected one transfer report in {directory}")
            metrics["transfers"] = {name: int(value) for name, value
                                   in re.findall(r"(\w+)=(\d+)", counters[0])}
            metrics["output_directory"] = str(directory)
            runs.append(metrics)
        result[mode] = {
            "runs": runs,
            "median_simulation_seconds": statistics.median(
                run["Simulation time"] for run in runs),
            "median_properties_update_seconds": statistics.median(
                run["Props/update time"] for run in runs),
        }
    result["simulation_speedup"] = (
        result["disabled"]["median_simulation_seconds"] /
        result["enabled"]["median_simulation_seconds"])
    result["properties_update_speedup"] = (
        result["disabled"]["median_properties_update_seconds"] /
        result["enabled"]["median_properties_update_seconds"])
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
