#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Create a 120,000-cell SPE1-derived scaling case, not a new field model."""

import argparse
import hashlib
import json
import re
from pathlib import Path


def prepare(source: Path, target: Path):
    raw = source.read_bytes()
    text = raw.decode()

    def replace(keyword, body):
        nonlocal text
        pattern = rf"(?m)^{keyword}[ \t]*\n.*?/"
        text, count = re.subn(
            pattern, lambda _: keyword + "\n" + body + " /", text, count=1, flags=re.S
        )
        if count != 1:
            raise ValueError(f"Expected keyword {keyword}")

    # This transformation is deliberately specific to the public SPE1 case 2.
    stripped = re.sub(r"--[^\n]*", "", text)
    if not re.search(r"\bDIMENS\s+10\s+10\s+3\s*/", stripped):
        raise ValueError("Expected the 10 x 10 x 3 SPE1 case 2 grid")
    text, titles = re.subn(
        r"(?m)^TITLE[ \t]*\n[^\n]*",
        "TITLE\nSPE1-derived 120000-cell refinement",
        text,
        count=1,
    )
    if titles != 1:
        raise ValueError("Expected TITLE")
    replace("DIMENS", "100 100 12")
    replace("WELLDIMS", "2 4 1 2")
    replace("DX", "120000*100")
    replace("DY", "120000*100")
    replace("DZ", "40000*5 40000*7.5 40000*12.5")
    replace("TOPS", "10000*8325")
    replace("PORO", "120000*0.3")
    for keyword in ("PERMX", "PERMY", "PERMZ"):
        replace(keyword, "40000*500 40000*50 40000*200")
    # Replace multi-record keywords separately, including their terminators.
    for keyword, body in (
        ("WELSPECS", "'PROD' 'G1' 100 100 8400 'OIL' /\n'INJ' 'G1' 1 1 8335 'GAS' /\n/"),
        ("COMPDAT", "'PROD' 100 100 9 12 'OPEN' 1* 1* 0.5 /\n'INJ' 1 1 1 4 'OPEN' 1* 1* 0.5 /\n/"),
    ):
        text, count = re.subn(
            rf"(?ms)^{keyword}[ \t]*\n.*?^/[ \t]*$", lambda _: keyword + "\n" + body, text, count=1
        )
        if count != 1:
            raise ValueError(f"Expected multi-record keyword {keyword}")
    replace("TSTEP", "1 2 3 4")
    text = re.sub(r"(?m)^\s*(ECHO|NOECHO)\s*$", r"-- Printing control removed: \1", text)
    text = (
        "-- Derived scaling experiment; original comments below describe the source case.\n" + text
    )
    manifest = {
        "source": str(source.resolve()),
        "source_sha256": hashlib.sha256(raw).hexdigest(),
        "deck_sha256": hashlib.sha256(text.encode()).hexdigest(),
        "grid": [100, 100, 12],
        "active_cells_expected": 120000,
        "days": 10,
        "scope": "Synthetic SPE1-derived refinement; not a validated industrial field model",
        "changes": [
            "10x refinement in x/y and 4x within each original layer; same reservoir extents and pore volume",
            "Preserve PVT, relative permeability, layer permeabilities, porosity and well rate/BHP controls",
            "Move producer to (100,100), complete producer in layers 9-12 and injector in 1-4",
            "Well indices are recomputed on the refined cells; this changes well discretization",
            "Use report steps of 1, 2, 3, 4 days; retain requested SUMMARY and restart output",
            "Remove printing controls only; preserve original dataset license notice",
        ],
    }
    target.mkdir(parents=True, exist_ok=False)
    (target / "SPE1_REFINED.DATA").write_text(text)
    (target / "preparation.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    args = parser.parse_args()
    print(json.dumps(prepare(args.source, args.target), indent=2))
