#!/usr/bin/env python3
"""
Extracts accepted substep end times from an INFOSTEP file and baseline log.
Used by run-timestep-replay-regressionTest.sh for hardcoded timestep replay.
"""

from pathlib import Path
import re
import sys


def extract_timesteps(infostep_path: Path, log_path: Path, output_path: Path):
    if not infostep_path.exists():
        raise FileNotFoundError(f"INFOSTEP file not found: {infostep_path}")
    if not log_path.exists():
        raise FileNotFoundError(f"Simulation log not found: {log_path}")

    lines = [line.strip() for line in infostep_path.read_text().splitlines() if line.strip()]
    header_idx = next((i for i, line in enumerate(lines) if "Time(day)" in line and "Conv" in line), None)
    if header_idx is None:
        raise ValueError(f"Unable to locate INFOSTEP header in {infostep_path}")

    header = lines[header_idx].split()
    time_idx = header.index("Time(day)")
    conv_idx = header.index("Conv")

    accepted_times = []
    for row in lines[header_idx + 1:]:
        cols = row.split()
        if len(cols) <= max(time_idx, conv_idx):
            continue
        if cols[conv_idx] not in {"1", "1.0", "true", "True"}:
            continue
        accepted_times.append(cols[time_idx])

    if not accepted_times:
        raise ValueError(f"No accepted timesteps found in {infostep_path}")

    # INFOSTEP currently omits the last accepted endpoint of the full simulation,
    # so append the total simulation time from the baseline log if needed.
    final_time = None
    for line in reversed(log_path.read_text().splitlines()):
        match = re.search(r"day\s+[^/]+/(\S+)", line)
        if match:
            final_time = match.group(1).rstrip(",;")
            break

    if final_time is None:
        raise ValueError(f"Unable to determine final simulation time from {log_path}")

    if accepted_times[-1] != final_time:
        accepted_times.append(final_time)

    output_path.write_text("\n".join(accepted_times) + "\n")


def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} <infostep_path> <log_path> <output_file>", file=sys.stderr)
        sys.exit(1)

    infostep_path = Path(sys.argv[1])
    log_path = Path(sys.argv[2])
    output_path = Path(sys.argv[3])

    extract_timesteps(infostep_path, log_path, output_path)


if __name__ == "__main__":
    main()
