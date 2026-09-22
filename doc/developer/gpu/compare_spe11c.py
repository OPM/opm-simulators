#!/usr/bin/env python3
"""Compare SPE11C restart/summary output without third-party Python modules.

Usage: python3 compare_spe11c.py BASELINE_DIRECTORY CANDIDATE_DIRECTORY
Pressure tolerances are applied in Pa; temperature tolerances in K. The
relative tolerance is 1e-6. Missing composition output is reported explicitly.
This checks saved report states; it does not replace per-update shadow checks.
"""

import argparse
import itertools
import json
import math
from pathlib import Path
import re
import struct


def records(path):
    """Read big-endian Eclipse unformatted records, checking block framing."""
    with path.open("rb") as stream:
        def block():
            size = stream.read(4)
            if not size:
                return None
            if len(size) != 4:
                raise ValueError(f"Truncated record in {path}")
            length = struct.unpack(">i", size)[0]
            if length < 0:
                raise ValueError(f"Invalid block length in {path}")
            data = stream.read(length)
            if len(data) != length or stream.read(4) != size:
                raise ValueError(f"Invalid block framing in {path}")
            return data

        while (header := block()) is not None:
            if len(header) != 16:
                raise ValueError(f"Invalid keyword header in {path}")
            name = header[:8].decode("ascii").strip()
            count = struct.unpack(">i", header[8:12])[0]
            kind = header[12:16].decode("ascii")
            widths = {"INTE": 4, "REAL": 4, "DOUB": 8, "LOGI": 4,
                      "CHAR": 8, "MESS": 0}
            width = int(kind[1:]) if kind.startswith("C0") else widths[kind]
            data = bytearray()
            while len(data) < count * width:
                chunk = block()
                if chunk is None or not chunk:
                    raise ValueError(f"Truncated {name} in {path}")
                data.extend(chunk)
            if len(data) != count * width:
                raise ValueError(f"Incorrect size for {name} in {path}")
            if kind in ("INTE", "LOGI", "REAL", "DOUB"):
                code = {"INTE": "i", "LOGI": "i", "REAL": "f", "DOUB": "d"}[kind]
                values = struct.unpack(f">{count}{code}", data)
            elif width:
                values = tuple(data[i:i + width].decode("ascii").strip()
                               for i in range(0, len(data), width))
            else:
                values = ()
            yield name, kind, values


def unique_file(directory, suffix):
    matches = list(directory.glob(f"*.{suffix}"))
    if len(matches) != 1:
        raise ValueError(f"Expected one {suffix} file in {directory}")
    return matches[0]


def numeric_compare(reference, candidate, field, location, scale=1, offset=0, atol=1e-8):
    worst = 0.0
    for index, (a, b) in enumerate(zip(reference, candidate, strict=True)):
        a, b = a * scale + offset, b * scale + offset
        if not math.isfinite(a) or not math.isfinite(b):
            raise ValueError(f"Nonfinite {field} at {location}, index {index}")
        error = abs(a - b)
        tolerance = atol + 1e-6 * abs(a)
        worst = max(worst, error / tolerance)
        if error > tolerance:
            raise ValueError(f"{field} mismatch at {location}, index {index}: "
                             f"reference={a:.17g}, candidate={b:.17g}, "
                             f"error={error:.17g}, tolerance={tolerance:.17g}")
    return worst


def compare_restart(reference, candidate):
    fields = {"PRESSURE", "PCGW", "TEMP", "SWAT", "SGAS", "SOIL",
              "RS", "RV", "RSW", "RVW"}
    observed, maxima, count, step = set(), {}, 0, None
    pressure_scale, temperature_offset = None, None
    for a, b in itertools.zip_longest(records(reference), records(candidate)):
        if a is None or b is None or a[:2] != b[:2] or len(a[2]) != len(b[2]):
            raise ValueError(f"Restart record structure differs at report {step}")
        name, _, values = a
        if name == "SEQNUM":
            if values != b[2]:
                raise ValueError("Restart report sequence differs")
            step = values[0]
            count += 1
        elif name == "INTEHEAD":
            if values[2] != b[2][2]:
                raise ValueError("Restart unit conventions differ")
            if values[2] != 1:
                raise ValueError("Comparator currently supports METRIC restart output only")
            pressure_scale, temperature_offset = 1e5, 273.15
        elif name in fields:
            if pressure_scale is None:
                raise ValueError("Missing restart unit convention")
            pressure = name in {"PRESSURE", "PCGW"}
            temperature = name == "TEMP"
            ratio = numeric_compare(values, b[2], name, f"report {step}",
                                    scale=pressure_scale if pressure else 1,
                                    offset=temperature_offset if temperature else 0,
                                    atol=1e-2 if pressure else 1e-6 if temperature else 1e-8)
            maxima[name] = max(maxima.get(name, 0), ratio)
            observed.add(name)
    required = {"PRESSURE", "TEMP", "SWAT", "SGAS"}
    if not required <= observed:
        raise ValueError(f"Missing required restart fields: {required - observed}")
    return {"reports": count, "fields": sorted(observed),
            "composition_present": bool(observed & {"RS", "RV", "RSW", "RVW"}),
            "maximum_error_fraction_of_tolerance": maxima}


def compare_summary(reference, candidate):
    spec_a = {name: values for name, _, values in records(unique_file(reference, "SMSPEC"))}
    spec_b = {name: values for name, _, values in records(unique_file(candidate, "SMSPEC"))}
    for name in ("KEYWORDS", "WGNAMES", "NUMS", "UNITS"):
        if spec_a.get(name) != spec_b.get(name):
            raise ValueError(f"Summary specification {name} differs")
    units = spec_a["UNITS"]
    maxima, samples = {}, 0
    # These describe execution speed, which is expected to change. They are
    # not simulation state; physical TIME/TIMESTEP and iteration counts remain.
    timing_keywords = {"ELAPSED", "TCPU", "TCPUDAY", "TCPUTS", "TELAPLIN"}
    convergence_keywords = {"MLINEARS", "MSUMLINS", "MSUMNEWT", "NEWTON",
                            "NLINEARS", "NLINSMAX", "NLINSMIN"}
    convergence_differences = {}
    ra = records(unique_file(reference, "UNSMRY"))
    rb = records(unique_file(candidate, "UNSMRY"))
    for a, b in itertools.zip_longest(ra, rb):
        if a is None or b is None or a[:2] != b[:2] or len(a[2]) != len(b[2]):
            raise ValueError(f"Summary structure differs at sample {samples}")
        if a[0] != "PARAMS":
            if a[2] != b[2]:
                raise ValueError(f"Summary metadata {a[0]} differs at sample {samples}")
            continue
        samples += 1
        for i, (x, y, unit) in enumerate(zip(a[2], b[2], units, strict=True)):
            if spec_a["KEYWORDS"][i] in timing_keywords:
                continue
            if spec_a["KEYWORDS"][i] in convergence_keywords:
                keyword = spec_a["KEYWORDS"][i]
                if x != y:
                    entry = convergence_differences.setdefault(keyword, {
                        "first_sample": samples, "first_reference": x,
                        "first_candidate": y, "different_samples": 0,
                        "maximum_absolute_difference": 0})
                    entry["different_samples"] += 1
                    entry["maximum_absolute_difference"] = max(
                        entry["maximum_absolute_difference"], abs(x - y))
                continue
            pressure = unit in {"BAR", "BARA", "BARSA"}
            temperature = unit in {"C", "DEGC", "K"}
            label = f"{i}:{spec_a['KEYWORDS'][i]}:{unit}"
            ratio = numeric_compare([x], [y], label, f"summary sample {samples}",
                                    scale=1e5 if pressure else 1,
                                    offset=273.15 if unit in {"C", "DEGC"} else 0,
                                    atol=1e-2 if pressure else 1e-6 if temperature else 1e-8)
            maxima[label] = max(maxima.get(label, 0), ratio)
    return {"samples": samples, "excluded_execution_timing_keywords": sorted(timing_keywords),
            "convergence_summary_differences": convergence_differences,
            "maximum_error_fraction_of_tolerance": maxima}


def run_statistics(directory):
    log = unique_file(directory, "PRT").read_text()
    if "End of simulation" not in log:
        raise ValueError(f"Incomplete simulation in {directory}")
    result = {}
    for label in ("Number of timesteps", "Simulation time", "Props/update time",
                  "Overall Linearizations", "Overall Newton Iterations", "Overall Linear Iterations"):
        match = re.search(re.escape(label) + r":\s*([0-9.]+)", log)
        if not match:
            raise ValueError(f"Missing {label} in {directory}")
        result[label] = float(match[1])
    history = []
    for line in log.splitlines():
        if line.startswith("Starting time step") or " Newton its=" in line:
            history.append(re.sub(r"\([0-9.]+sec\)", "(timing)", line))
    return result, history


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    stats_a, hist_a = run_statistics(args.reference)
    stats_b, hist_b = run_statistics(args.candidate)
    result = {"reference": str(args.reference), "candidate": str(args.candidate),
              "reference_statistics": stats_a, "candidate_statistics": stats_b,
              "identical_convergence_history": hist_a == hist_b,
              "restart": compare_restart(unique_file(args.reference, "UNRST"),
                                         unique_file(args.candidate, "UNRST")),
              "summary": compare_summary(args.reference, args.candidate)}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
