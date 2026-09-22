#!/usr/bin/env python3
"""Correlate rocprofv3 CSV copies with API calls, streams and property kernels.

Usage: python3 summarize_spe11c_trace.py DISABLED_TRACE_DIR ENABLED_TRACE_DIR
This ROCm CSV schema omits copy sizes. Byte counts come from bridge counters;
this script corroborates calls using their API, direction, stream and ordering.
"""

import argparse
import bisect
import collections
import csv
import json
from pathlib import Path


def read_trace(directory, suffix):
    paths = list(directory.rglob(f"*_{suffix}.csv"))
    if len(paths) != 1:
        raise ValueError(f"Expected one {suffix} CSV in {directory}")
    with paths[0].open() as stream:
        return list(csv.DictReader(stream))


def summarize(directory):
    kernels = sorted(read_trace(directory, "kernel_trace"),
                     key=lambda row: int(row["Start_Timestamp"]))
    copies = read_trace(directory, "memory_copy_trace")
    apis = {row["Correlation_Id"]: row for row in read_trace(directory, "hip_api_trace")}
    properties = [row for row in kernels if "dispatcherUpdateAllCellsKernel" in row["Kernel_Name"]]
    newton = [row for row in kernels if "blackoilNewtonUpdateKernel" in row["Kernel_Name"]]
    property_streams = {row["Stream_Id"] for row in properties}
    grouped = collections.Counter((row["Stream_Id"], row["Direction"],
                                   apis.get(row["Correlation_Id"], {}).get("Function", "unmatched"))
                                  for row in copies)
    pv_uploads = [row for row in copies
                  if row["Stream_Id"] in property_streams
                  and row["Direction"] == "MEMORY_COPY_HOST_TO_DEVICE"
                  and apis.get(row["Correlation_Id"], {}).get("Function") == "hipMemcpyAsync"]
    correction_downloads = [row for row in copies
                            if row["Stream_Id"] == "0"
                            and row["Direction"] == "MEMORY_COPY_DEVICE_TO_HOST"
                            and apis.get(row["Correlation_Id"], {}).get("Function") == "hipMemcpyAsync"]
    starts = [int(row["Start_Timestamp"]) for row in kernels]

    def next_property_count(selected):
        count = 0
        for row in selected:
            index = bisect.bisect_right(starts, int(row["End_Timestamp"]))
            if index < len(kernels) and "dispatcherUpdateAllCellsKernel" in kernels[index]["Kernel_Name"]:
                count += 1
        return count

    return {
        "directory": str(directory),
        "requested_collection_period": "5:2:1 seconds (one two-second window after a five-second delay)",
        "copy_sizes_present_in_csv": "Bytes" in copies[0] if copies else False,
        "property_kernels": len(properties),
        "resident_newton_kernels": len(newton),
        "assembly_kernels": sum("void Opm::kernel_linearize<" in row["Kernel_Name"] for row in kernels),
        "property_streams": sorted(property_streams),
        "pv_upload_pattern_calls": len(pv_uploads),
        "pv_uploads_followed_by_property_kernel": next_property_count(pv_uploads),
        "correction_download_pattern_calls": len(correction_downloads),
        "correction_downloads_followed_by_property_kernel": next_property_count(correction_downloads),
        "copies_by_stream_direction_api": [
            {"stream": stream, "direction": direction, "api": api, "calls": count}
            for (stream, direction, api), count in sorted(grouped.items())],
        "api_counts": dict(sorted(collections.Counter(row["Function"] for row in apis.values()).items())),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("disabled", type=Path)
    parser.add_argument("enabled", type=Path)
    args = parser.parse_args()
    print(json.dumps({"disabled": summarize(args.disabled), "enabled": summarize(args.enabled)}, indent=2))


if __name__ == "__main__":
    main()
