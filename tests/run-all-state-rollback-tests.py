#!/usr/bin/env python3
"""
Unified State Rollback & Timestep Replay Test Suite for OPM Flow.

Executes replay regression tests across all state-rollback features
(ROCKCOMP hysteresis, WCYCLE timers, MSW, ACTIONX, Group controls, etc.),
extracts simulation convergence metrics, inspects ECLIPSE summary vectors,
and displays formatted comparison tables.

Note on Linear Solver Replay Determinism (--cpr-reuse-setup=0):
In baseline runs with chopped/retried timesteps, Newton iterations attempted
during failed steps advance the linear solver setup counter (solveCount_),
triggering AMG preconditioner recreation at different physical simulation times
than in hardcoded replay runs (which never execute failed iterations).
Passing `--cpr-reuse-setup=0` forces per-solve preconditioner generation,
ensuring exact bitwise trajectory identity and isolating physical state rollback
verification from heuristic linear solver caching cadences.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path


# Color formatting for terminal
class Colors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    CYAN = '\033[96m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

    @classmethod
    def disable(cls):
        cls.HEADER = ''
        cls.BLUE = ''
        cls.CYAN = ''
        cls.GREEN = ''
        cls.YELLOW = ''
        cls.RED = ''
        cls.BOLD = ''
        cls.UNDERLINE = ''
        cls.END = ''


TEST_CASES = [
    {
        "name": "SPE1CASE2_ROCK2DTR",
        "category": "ROCKCOMP / Hysteresis",
        "description": "Rock compaction hysteresis & running max/min tracking (16 failed steps)",
        "deck_dir": "spe1",
        "deck_file": "SPE1CASE2_ROCK2DTR.DATA",
        "args": ["--truncate-time-step-to-float=true", "--full-time-step-initially=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WOPT:PROD", "WGPT:PROD", "WWPT:PROD", "WBHP:PROD"],
        "vector_labels": ["Oil Prod (WOPT)", "Gas Prod (WGPT)", "Water Prod (WWPT)", "BHP (WBHP)"],
    },
    {
        "name": "0_BASE_MODEL2",
        "category": "Well Model / Step Cuts",
        "description": "Multi-well blackoil with convergence failures & step chops",
        "deck_dir": "model2",
        "deck_file": "0_BASE_MODEL2.DATA",
        "args": [
            "--truncate-time-step-to-float=true",
            "--full-time-step-initially=true",
            "--cpr-reuse-setup=0",
            "--newton-max-iterations=8",
            "--tolerance-cnv-relaxed=1e-3",
        ],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "WCYCLE-1",
        "category": "WCYCLE (Well Cycling)",
        "description": "Periodic well cycling open/close timers (well_open_close_time)",
        "deck_dir": "wcycle",
        "deck_file": "WCYCLE-1.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "2_MSWWELL_MODEL3",
        "category": "Multisegment Wells (MSW)",
        "description": "Segmented wells, crossflow, and segment pressure/friction",
        "deck_dir": "model3",
        "deck_file": "2_MSWWELL_MODEL3.DATA",
        "args": ["--truncate-time-step-to-float=true", "--cpr-reuse-setup=0"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "ACTIONX_GCONPROD",
        "category": "ACTIONX Events",
        "description": "Dynamic group production control & well rate adjustments",
        "deck_dir": "actionx",
        "deck_file": "ACTIONX_GCONPROD.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "ACTIONX_UDQ",
        "category": "ACTIONX + UDQ",
        "description": "Dynamic action triggers combined with User-Defined Quantities",
        "deck_dir": "actionx",
        "deck_file": "ACTIONX_UDQ.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "MOD4_GRP",
        "category": "Group Control / NUPCOL",
        "description": "Group hierarchy state machine & dynamic rate allocation",
        "deck_dir": "model4",
        "deck_file": "MOD4_GRP.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["GOPR:MOD4", "GGPR:MOD4", "GWPR:MOD4", "WBHP:OPL1"],
        "vector_labels": ["Group Oil (GOPR)", "Group Gas (GGPR)", "Group Water (GWPR)", "Well BHP (WBHP)"],
    },
    {
        "name": "SPE1CASE1",
        "category": "DRSDT / Base Blackoil",
        "description": "Standard 3-phase blackoil with dissolved gas rate control",
        "deck_dir": "spe1",
        "deck_file": "SPE1CASE1.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WOPT:PROD", "WGPT:PROD", "WWPT:PROD", "WBHP:PROD"],
        "vector_labels": ["Oil Prod (WOPT)", "Gas Prod (WGPT)", "Water Prod (WWPT)", "BHP (WBHP)"],
    },
    {
        "name": "2D_OW_CTAQUIFER",
        "category": "Carter-Tracy Aquifer",
        "description": "Analytical Carter-Tracy aquifer influence functions & water influx",
        "deck_dir": "aquifer-oilwater",
        "deck_file": "2D_OW_CTAQUIFER.DATA",
        "simulator": "flow_oilwater",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "2D_FETKOVICHAQUIFER",
        "category": "Fetkovich Aquifer",
        "description": "Analytical Fetkovich aquifer pressure/flux coupling",
        "deck_dir": "aquifer-fetkovich",
        "deck_file": "2D_FETKOVICHAQUIFER.DATA",
        "simulator": "flow_oilwater",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "2D_THREEPHASE_POLY_HETER",
        "category": "Polymer / Adsorption",
        "description": "Polymer injection, viscosity enhancement & irreversible adsorption",
        "deck_dir": "polymer_simple2D",
        "deck_file": "2D_THREEPHASE_POLY_HETER.DATA",
        "simulator": "flow_polymer",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "SPE1FOAM",
        "category": "Foam / Mobility Reduction",
        "description": "Gas foam injection & mobility reduction factor modeling",
        "deck_dir": "spe1_foam",
        "deck_file": "SPE1FOAM.DATA",
        "simulator": "flow_foam",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WOPT:PROD", "WGPT:PROD", "WWPT:PROD", "WBHP:PROD"],
        "vector_labels": ["Oil Prod (WOPT)", "Gas Prod (WGPT)", "Water Prod (WWPT)", "BHP (WBHP)"],
    },
    {
        "name": "7_HYSTERESIS_MODEL2",
        "category": "Relperm Hysteresis",
        "description": "EHYSTR Killough / Carlson relative permeability scanning curves",
        "deck_dir": "model2",
        "deck_file": "7_HYSTERESIS_MODEL2.DATA",
        "args": ["--truncate-time-step-to-float=true", "--cpr-reuse-setup=0"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "BASE_WT_TRACER",
        "category": "Water Tracers",
        "description": "Sequential water tracer transport & concentration tracking",
        "deck_dir": "tracer",
        "deck_file": "BASE_WT_TRACER.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "SPE1CASE1_BRINE",
        "category": "Brine / Salinity",
        "description": "Brine / salt concentration tracking and water density/viscosity effects",
        "deck_dir": "spe1_brine",
        "deck_file": "SPE1CASE1_BRINE.DATA",
        "simulator": "flow_brine",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WOPT:PROD", "WGPT:PROD", "WWPT:PROD", "WBHP:PROD"],
        "vector_labels": ["Oil Prod (WOPT)", "Gas Prod (WGPT)", "Water Prod (WWPT)", "BHP (WBHP)"],
    },
    {
        "name": "VAPPARS-01",
        "category": "VAPPARS / Rv Control",
        "description": "Vaporized oil in gas phase with DRVDT / VAPPARS rate limits",
        "deck_dir": "vappars",
        "deck_file": "VAPPARS-01.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WOPR:P-1H", "WGPR:P-1H", "WBHP:P-1H", "FPR"],
        "vector_labels": ["Oil Rate (WOPR)", "Gas Rate (WGPR)", "Well BHP (WBHP)", "Avg Press (FPR)"],
    },
    {
        "name": "WAGHYSTR-01",
        "category": "WAG Relperm Hysteresis",
        "description": "Water-Alternating-Gas relative permeability hysteresis cycles",
        "deck_dir": "waghystr",
        "deck_file": "WAGHYSTR-01.DATA",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["WBHP:PROD1", "WBHP:INJ_G1", "WBHP:INJ_W1", "FPR"],
        "vector_labels": ["Prod BHP (PROD1)", "Inj Gas BHP (INJ_G1)", "Inj Water BHP (INJ_W1)", "Avg Press (FPR)"],
    },
    {
        "name": "EXTREME_ROLLBACK_STRESS",
        "category": "Extreme Stress Multi-Physics",
        "description": "Live fluid + VAPOIL + ROCKCOMP + WCYCLE + ACTIONX under forced step cuts",
        "deck_dir": "model_synthetic",
        "deck_file": "EXTREME_ROLLBACK_STRESS.DATA",
        "args": ["--truncate-time-step-to-float=true", "--full-time-step-initially=true", "--newton-max-iterations=6", "--cpr-reuse-setup=0"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "FORCED_ROLLBACK_MULTICUT",
        "category": "Forced Cascading Multi-Cut",
        "description": "Multiple forced Newton convergence failures and cascading cuts on same timestep",
        "deck_dir": "model_synthetic",
        "deck_file": "FORCED_ROLLBACK_MULTICUT.DATA",
        "args": ["--truncate-time-step-to-float=true", "--full-time-step-initially=true", "--newton-max-iterations=4", "--cpr-reuse-setup=0"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
    {
        "name": "1D_COMP",
        "category": "Compositional Flow & Wells",
        "description": "3-component Peng-Robinson EOS compositional flow with CompWellModel",
        "deck_dir": "compositional",
        "deck_file": "1D_COMP.DATA",
        "simulator": "flow_comp",
        "args": ["--truncate-time-step-to-float=true"],
        "replay_args": ["--initial-time-step-in-days=11111111"],
        "key_vectors": ["FOPT", "FGPT", "FWPT", "FPR"],
        "vector_labels": ["Oil Prod (FOPT)", "Gas Prod (FGPT)", "Water Prod (FWPT)", "Avg Press (FPR)"],
    },
]


def find_default_paths():
    """Locates repo directories, build directories, and binaries."""
    script_dir = Path(__file__).resolve().parent
    # Check potential repo root paths
    candidates = [
        script_dir.parent.parent,
        Path("/home/ahmed/Documents/Projects/opm/src"),
        Path.cwd(),
    ]
    repo_root = None
    for cand in candidates:
        if (cand / "opm-simulators").exists() and (cand / "opm-tests").exists():
            repo_root = cand.resolve()
            break
    if not repo_root:
        repo_root = script_dir.parent.parent.resolve()

    tests_dir = repo_root / "opm-tests"
    if not tests_dir.exists():
        tests_dir = repo_root.parent / "src" / "opm-tests"

    build_candidates = [
        repo_root.parent / "build-mpi",
        repo_root / "build-mpi",
        repo_root / "cmake-build",
        repo_root / "build",
    ]
    build_dir = None
    for bcand in build_candidates:
        if (bcand / "opm-simulators" / "bin" / "flow_blackoil").exists():
            build_dir = bcand.resolve()
            break

    if not build_dir:
        build_dir = repo_root.parent / "build-mpi"

    return repo_root, tests_dir, build_dir


def parse_log_stats(log_path):
    """Extracts execution metrics and wasted Newton iterations from Flow log."""
    stats = {
        "timesteps": 0,
        "sim_time": 0.0,
        "newton_total": 0,
        "newton_wasted": 0,
        "linear_total": 0,
        "linear_wasted": 0,
        "failed_steps": 0,
    }
    if not log_path.exists():
        return stats

    content = log_path.read_text(errors="ignore")

    # Count convergence failures / chops
    stats["failed_steps"] = len(re.findall(r"Cutting time step|Oscillating behavior detected|NonlinearSolver\.hpp", content))

    for line in content.splitlines():
        if "Number of timesteps:" in line:
            m = re.search(r"Number of timesteps:\s+(\d+)", line)
            if m:
                stats["timesteps"] = int(m.group(1))
        elif "Simulation time:" in line:
            m = re.search(r"Simulation time:\s+([\d\.]+)\s*s", line)
            if m:
                stats["sim_time"] = float(m.group(1))
        elif "Overall Newton Iterations:" in line:
            m = re.search(r"Overall Newton Iterations:\s+(\d+)(?:\s+\(Wasted:\s+(\d+))?", line)
            if m:
                stats["newton_total"] = int(m.group(1))
                if m.group(2):
                    stats["newton_wasted"] = int(m.group(2))
        elif "Overall Linear Iterations:" in line:
            m = re.search(r"Overall Linear Iterations:\s+(\d+)(?:\s+\(Wasted:\s+(\d+))?", line)
            if m:
                stats["linear_total"] = int(m.group(1))
                if m.group(2):
                    stats["linear_wasted"] = int(m.group(2))

    return stats


def extract_summary_values(summary_bin, smspec_file, vectors):
    """Uses opm-common `summary` binary to read final report-step vector values."""
    results = {}
    if not summary_bin.exists() or not smspec_file.exists():
        return results

    cmd = [str(summary_bin), "-n", "-r", str(smspec_file)] + vectors
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=10)
        if proc.returncode == 0 and proc.stdout.strip():
            last_line = proc.stdout.strip().splitlines()[-1]
            tokens = last_line.split()
            if len(tokens) == len(vectors):
                for vec, val in zip(vectors, tokens):
                    try:
                        results[vec] = float(val)
                    except ValueError:
                        results[vec] = val
    except Exception:
        pass

    return results


def format_table(headers, rows, alignments=None):
    """Formats a neat, aligned ASCII/Unicode table."""
    if not rows:
        return ""
    col_count = len(headers)
    if not alignments:
        alignments = ['left'] * col_count

    col_widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            clean_cell = re.sub(r'\033\[[0-9;]*m', '', str(cell))
            if len(clean_cell) > col_widths[i]:
                col_widths[i] = len(clean_cell)

    # Padding
    pad = 2
    col_widths = [w + pad for w in col_widths]

    def render_row(cells, aligns):
        out = ""
        for i, cell in enumerate(cells):
            clean_len = len(re.sub(r'\033\[[0-9;]*m', '', str(cell)))
            spaces = col_widths[i] - clean_len
            if aligns[i] == 'right':
                out += " " * spaces + str(cell)
            elif aligns[i] == 'center':
                left_spaces = spaces // 2
                right_spaces = spaces - left_spaces
                out += " " * left_spaces + str(cell) + " " * right_spaces
            else:
                out += str(cell) + " " * spaces
        return out

    sep = "-" * sum(col_widths)
    header_str = render_row(headers, alignments)

    res = [sep, header_str, sep]
    for r in rows:
        res.append(render_row(r, alignments))
    res.append(sep)
    return "\n".join(res)


def run_test_case(case, repo_root, tests_dir, build_dir, result_root, use_mpi=True):
    """Runs a single test case through the replay regression driver."""
    replay_driver = repo_root / "opm-simulators" / "tests" / "run-timestep-replay-regressionTest.sh"
    sim_name = case.get("simulator", "flow_blackoil")
    sim_bin = build_dir / "opm-simulators" / "bin" / sim_name
    compare_bin = build_dir / "opm-common" / "bin" / "compareECL"
    summary_bin = build_dir / "opm-common" / "bin" / "summary"
    rst_bin = build_dir / "opm-common" / "bin" / "rst_deck"

    deck_path = tests_dir / case["deck_dir"]
    case_result_dir = result_root / case["name"]
    baseline_log = case_result_dir / f"{case['name']}.baseline.log"
    smspec_file = case_result_dir / "baseline" / f"{case['name']}.SMSPEC"

    cmd = []
    if use_mpi and shutil.which("mpirun"):
        cmd.extend(["mpirun", "-np", "1"])

    cmd.extend([
        str(replay_driver),
        "-i", str(deck_path),
        "-r", str(case_result_dir),
        "-b", str(sim_bin.parent),
        "-f", case["name"],
        "-a", "2e-14",
        "-t", "2e-14",
        "-c", str(compare_bin),
        "-d", str(rst_bin),
        "-u", sim_name,
        "-e", str(sim_bin),
    ])

    for rarg in case.get("replay_args", []):
        cmd.extend(["-y", rarg])

    cmd.append("--")
    cmd.extend(case.get("args", []))

    start_time = time.time()
    try:
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, timeout=120)
        elapsed = time.time() - start_time
        success = (proc.returncode == 0)
        output = proc.stdout + proc.stderr
    except subprocess.TimeoutExpired:
        elapsed = time.time() - start_time
        success = False
        output = "TIMEOUT EXPIRED (>120s)"
    except Exception as e:
        elapsed = time.time() - start_time
        success = False
        output = f"Execution error: {e}"

    stats = parse_log_stats(baseline_log)
    summary_vals = extract_summary_values(summary_bin, smspec_file, case.get("key_vectors", []))

    return {
        "name": case["name"],
        "category": case["category"],
        "description": case["description"],
        "success": success,
        "elapsed": elapsed,
        "stats": stats,
        "summary": summary_vals,
        "vector_labels": case.get("vector_labels", []),
        "key_vectors": case.get("key_vectors", []),
        "output": output,
    }


def main():
    parser = argparse.ArgumentParser(description="OPM Flow State Rollback & Timestep Replay Test Suite")
    parser.add_argument("--filter", "-k", type=str, help="Run only tests matching pattern (substring)")
    parser.add_argument("--build-dir", "-b", type=str, help="Path to build directory containing binaries")
    parser.add_argument("--test-dir", "-t", type=str, help="Path to opm-tests directory")
    parser.add_argument("--no-color", action="store_true", help="Disable colored output")
    parser.add_argument("--results-dir", "-r", type=str, default="/tmp/opm_state_rollback_test_results", help="Directory to store run outputs")
    parser.add_argument("--no-mpi", action="store_true", help="Do not wrap simulator execution with mpirun")
    args = parser.parse_args()

    if args.no_color:
        Colors.disable()

    repo_root, default_tests_dir, default_build_dir = find_default_paths()
    tests_dir = Path(args.test_dir) if args.test_dir else default_tests_dir
    build_dir = Path(args.build_dir) if args.build_dir else default_build_dir
    results_dir = Path(args.results_dir)

    print(f"{Colors.BOLD}{Colors.HEADER}========================================================================={Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}    OPM Flow Simulator — State Rollback & Replay Verification Suite      {Colors.END}")
    print(f"{Colors.BOLD}{Colors.HEADER}========================================================================={Colors.END}")
    print(f"Repo Root  : {repo_root}")
    print(f"Tests Dir  : {tests_dir}")
    print(f"Build Dir  : {build_dir}")
    print(f"Results Dir: {results_dir}")
    print(f"Comparison : Deterministic Bitwise Identity @ Machine Precision (2e-14)")
    print()

    active_cases = TEST_CASES
    if args.filter:
        active_cases = [c for c in TEST_CASES if re.search(args.filter, c["name"], re.IGNORECASE) or re.search(args.filter, c["category"], re.IGNORECASE)]

    if not active_cases:
        print(f"{Colors.RED}No tests matching filter '{args.filter}' found.{Colors.END}")
        sys.exit(1)

    print(f"Running {len(active_cases)} targeted test case(s)...\n")

    results = []
    total_start = time.time()

    for idx, case in enumerate(active_cases, 1):
        print(f"[{idx}/{len(active_cases)}] Running {Colors.BOLD}{case['name']}{Colors.END} ({case['category']})... ", end="", flush=True)
        res = run_test_case(case, repo_root, tests_dir, build_dir, results_dir, use_mpi=(not args.no_mpi))
        results.append(res)
        if res["success"]:
            print(f"{Colors.GREEN}[PASSED]{Colors.END} in {res['elapsed']:.2f}s")
        else:
            print(f"{Colors.RED}[FAILED]{Colors.END} in {res['elapsed']:.2f}s")

    total_elapsed = time.time() - total_start
    passed_count = sum(1 for r in results if r["success"])
    failed_count = len(results) - passed_count

    print(f"\n{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}                       TABLE 1: EXECUTION & REPLAY STATUS               {Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")

    t1_headers = ["Test Case", "Feature Tested", "Timesteps", "Failed / Cuts", "Newton Its (Wasted)", "Time (s)", "Replay Status"]
    t1_align = ["left", "left", "right", "right", "right", "right", "center"]
    t1_rows = []
    for r in results:
        status_str = f"{Colors.GREEN}MATCH (2e-14){Colors.END}" if r["success"] else f"{Colors.RED}MISMATCH / FAIL{Colors.END}"
        st = r["stats"]
        newton_str = f"{st['newton_total']} ({st['newton_wasted']})" if st['newton_total'] > 0 else "-"
        t1_rows.append([
            r["name"],
            r["category"],
            str(st["timesteps"]) if st["timesteps"] > 0 else "-",
            str(st["failed_steps"]),
            newton_str,
            f"{r['elapsed']:.2f}",
            status_str,
        ])
    print(format_table(t1_headers, t1_rows, t1_align))

    print(f"\n{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}                   TABLE 2: ECLIPSE SUMMARY TOTALS INSPECTION           {Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")

    t2_headers = ["Test Case", "Cumulative Oil", "Cumulative Gas", "Cumulative Water", "Pressure / BHP"]
    t2_align = ["left", "right", "right", "right", "right"]
    t2_rows = []
    for r in results:
        vals = []
        for vec in r["key_vectors"]:
            v = r["summary"].get(vec, "-")
            if isinstance(v, float):
                vals.append(f"{v:.4e}")
            else:
                vals.append(str(v))
        # Pad to 4 columns if needed
        while len(vals) < 4:
            vals.append("-")
        t2_rows.append([r["name"]] + vals[:4])
    print(format_table(t2_headers, t2_rows, t2_align))

    print(f"\n{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}                TABLE 3: STATE ROLLBACK MECHANISM MATRIX                 {Colors.END}")
    print(f"{Colors.BOLD}{Colors.CYAN}========================================================================={Colors.END}")

    t3_headers = ["Feature", "Snapshotted State Members", "Rollback Trigger", "Impact of Rollback"]
    t3_align = ["left", "left", "left", "left"]
    t3_rows = [
        ["ROCKCOMP (Hysteresis)", "max_p_, minRefPressure_, rockCompTransMultVal_", "updateFailed()", "Prevents unphysical pore compaction & scanning curve drift"],
        ["VAPPARS / DRSDT", "max_oil_saturation, max_water_saturation, mixControls", "updateFailed()", "Reverts trial saturation extrema to previous converged state"],
        ["WCYCLE (Well Cycling)", "well_open_close_time_", "updateFailed()", "Reverts periodic open/close timers when well steps chop"],
        ["MSW (Multisegment Wells)", "active_wgstate_, closed_this_timestep_", "updateFailed()", "Restores segment pressure/rate states & well control modes"],
        ["ACTIONX / UDQ Triggers", "Deferred to converged endTimeStep()", "convergenceFailed", "Actions evaluate strictly on converged states; no corruption"],
        ["Group Control / NUPCOL", "nupcol_wgstate_, group_state_, network_state_", "updateFailed()", "Prevents guide rate / group target contamination across retries"],
        ["Carter-Tracy Aquifer", "W_flux_ deferred; influence tables dynamic", "updateFailed()", "Prevents non-physical aquifer influx accumulation on retry"],
        ["Fetkovich Aquifer", "pressure_previous_, aquifer flux coupling", "updateFailed()", "Ensures aquifer pressure decay remains synchronized with step size"],
        ["Polymer / Adsorption", "max_poly_ads_ (irreversible adsorption)", "updateFailed()", "Prevents artificial polymer stripping across divergent trial steps"],
        ["Foam / Mobility Reduction", "Dynamic mobility reduction factors (Fw, Fs, Fg)", "updateFailed()", "Re-evaluated cleanly from primary variables on retried step"],
        ["Relperm Hysteresis", "Historical turning points & scanning curve flags", "updateFailed()", "Prevents locking into incorrect drainage/imbibition scanning curves"],
        ["Water Tracers", "Tracer transport deferred to converged step", "endTimeStep()", "Sequential solver executes only on accepted flow fields"],
        ["Brine / Salinity", "Instantaneous salinity state functions", "updateFailed()", "Salt concentration reverts with primary solution vector u(t)"],
    ]
    print(format_table(t3_headers, t3_rows, t3_align))

    print(f"\n{Colors.BOLD}Suite Summary:{Colors.END} {passed_count}/{len(results)} Passed ({failed_count} Failed) in {total_elapsed:.2f}s total wall time.")

    if failed_count > 0:
        print(f"\n{Colors.RED}Failure details:{Colors.END}")
        for r in results:
            if not r["success"]:
                print(f"\n--- {r['name']} Failure Output ---\n{r['output']}\n")
        sys.exit(1)
    else:
        print(f"{Colors.GREEN}All tested rollback features demonstrated 100% deterministic replay at machine precision!{Colors.END}\n")
        sys.exit(0)


if __name__ == "__main__":
    main()
