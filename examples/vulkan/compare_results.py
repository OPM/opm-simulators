#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Compare requested SUMMARY vectors and restart fields at matching report dates.

Requires numpy and resdata. Thresholds are fixed; solver counters are reported
separately. A successful process exit alone does not establish correctness.
"""

import argparse
import json
from contextlib import closing
from pathlib import Path

import numpy as np

ABS_TOL = 2e-2
REL_TOL = 1e-5


def find(directory, suffix):
    matches = list(Path(directory).glob("*" + suffix))
    if len(matches) != 1:
        raise ValueError(f"Expected exactly one {suffix} file in {directory}")
    return matches[0]


def compare(reference: Path, candidate: Path, *, summary_diagnostics: tuple[str, ...] = ()):
    from resdata.resfile import ResdataFile
    from resdata.summary import Summary

    failures, checked = [], 0
    diagnostics = {}

    def check(label, a, b):
        nonlocal checked
        checked += 1
        a, b = np.asarray(a), np.asarray(b)
        # Conservative: absolute OR relative test, rather than adding tolerances.
        if a.shape != b.shape or not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
            failures.append(label + ": shape mismatch or nonfinite values")
            return
        delta = np.abs(a - b)
        valid = (delta <= ABS_TOL) | (delta <= REL_TOL * np.maximum(np.abs(a), np.abs(b)))
        if not np.all(valid):
            failures.append(f"{label}: maximum absolute difference {float(delta.max()):.8g}")

    a, b = Summary(str(find(reference, ".SMSPEC"))), Summary(str(find(candidate, ".SMSPEC")))
    ak, bk = set(a.keys()), set(b.keys())
    if ak != bk:
        failures.append("SUMMARY vector sets differ")
    # Compare every requested vector at identical report dates, without interpolation.
    ad = np.asarray(a.report_dates, dtype="datetime64[ms]")
    bd = np.asarray(b.report_dates, dtype="datetime64[ms]")
    if not np.array_equal(ad, bd):
        failures.append("SUMMARY report dates differ")
    ai, bi = np.searchsorted(a.numpy_dates, ad), np.searchsorted(b.numpy_dates, bd)
    if not np.array_equal(a.numpy_dates[ai], ad) or not np.array_equal(b.numpy_dates[bi], bd):
        failures.append("Exact samples not found at report dates")
    for key in sorted(ak & bk):
        if a.unit(key) != b.unit(key):
            failures.append(f"{key}: units differ")
        if key in summary_diagnostics:
            aa, bb = a.numpy_vector(key)[ai], b.numpy_vector(key)[bi]
            diagnostics[key] = {
                "reference": aa.tolist(),
                "candidate": bb.tolist(),
                "unit": a.unit(key),
                "reason": "Solver performance diagnostic; equality is not a physical criterion",
            }
            if aa.shape != bb.shape or not np.all(np.isfinite(aa)) or not np.all(np.isfinite(bb)):
                failures.append(key + ": diagnóstico con shape mismatch or nonfinite values")
            continue
        check(key, a.numpy_vector(key)[ai], b.numpy_vector(key)[bi])
    with (
        closing(ResdataFile(str(find(reference, ".UNRST")))) as ar,
        closing(ResdataFile(str(find(candidate, ".UNRST")))) as br,
    ):
        if ar.report_list != br.report_list or ar.report_dates != br.report_dates:
            failures.append("Restart reports differ")
        common = set(ar.report_list) & set(br.report_list)
        for step in sorted(common):
            av, bv = ar.restart_view(report_step=step), br.restart_view(report_step=step)
            for key in ("PRESSURE", "SWAT", "SGAS", "RS", "RV"):
                if (key in av) != (key in bv):
                    failures.append(f"{step}/{key}: missing property")
                elif key in av:
                    check(f"{step}/{key}", av[key][0].numpy_view(), bv[key][0].numpy_view())
    return {
        "passed": not failures and checked > 0,
        "checked": checked,
        "failures": failures,
        "absolute_tolerance": ABS_TOL,
        "relative_tolerance": REL_TOL,
        "summary_diagnostics": diagnostics,
        "scope": "Requested SUMMARY and restart properties; not general equivalence with ECLIPSE",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    args = parser.parse_args()
    result = compare(
        args.reference,
        args.candidate,
        summary_diagnostics=("MLINEARS", "MSUMLINS", "NLINEARS", "NLINSMAX", "NLINSMIN", "TCPU"),
    )
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["passed"] else 1)
