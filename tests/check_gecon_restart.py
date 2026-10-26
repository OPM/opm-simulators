#!/usr/bin/env python3
"""Check GECON-02 MPI outputs or UNRST reconstruction.

Use tight nonlinear convergence and Flow's extended double-precision UNRST
output to separate restart/parallel consistency from solve and file-quantization
errors.
Comparison tolerances remain the normal .02 absolute / 1e-5 relative values.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

from check_gecon_economic_events import read_summary, require, validate

OPTIONS = ['--threads-per-process=1', '--group-economic-wgr-step-tolerance=0.001',
           '--initial-time-step-in-days=1', '--tolerance-cnv=1e-7',
           '--tolerance-cnv-relaxed=1e-7', '--tolerance-mb=1e-10',
           '--tolerance-mb-relaxed=1e-10', '--tolerance-wells=1e-8',
           '--enable-opm-rst-file=true', '--ecl-output-double-precision=true']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('flow', 'summary', 'compare', 'deck', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--mode', choices=('mpi', 'unrst'), required=True)
    parser.add_argument('--rst-deck', type=Path)
    parser.add_argument('--mpi-exec', type=Path)
    parser.add_argument('--mpi-numproc-flag', default='-n')
    args = parser.parse_args()
    if args.mode == 'mpi' and not args.mpi_exec:
        parser.error('--mode mpi requires --mpi-exec')
    if args.mode == 'unrst' and not args.rst_deck:
        parser.error('--mode unrst requires --rst-deck')
    args.output.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix=args.mode + '-', dir=str(args.output.resolve())))
    print('GECON consistency artifacts: ' + str(root), flush=True)
    commands = []

    def command(argv, log, env=None):
        argv = list(map(str, argv))
        commands.append(argv)
        (root / 'commands.json').write_text(json.dumps(commands, indent=2))
        with log.open('w') as stream:
            subprocess.run(argv, stdout=stream, stderr=subprocess.STDOUT,
                           env=env, check=True, timeout=600)

    source = root / 'input'
    source.mkdir()
    text = args.deck.read_text()
    text = re.sub(r"(INCLUDE\s*\n\s*)'([^']+)'",
                  lambda m: m[1] + "'" + str((args.deck.parent / m[2]).resolve()) + "'",
                  text, flags=re.IGNORECASE)
    # Add report boundaries immediately before the two economic closures.
    for marker, extra in [("  20 'JAN' 2020 /", " 24 'JAN' 2020 /"),
                          ("  01 'MAR' 2020 /", " 8 'MAR' 2020 /")]:
        require(text.count(marker) == 1, 'Missing GECON report boundary: ' + marker)
        if args.mode == 'unrst':
            extra = extra[:-1] + ("'06:00:00' /" if 'JAN' in extra else "'20:24:00' /")
        text = text.replace(marker, marker + '\n' + extra)
    deck = source / 'GECON-02.DATA'
    deck.write_text(text)
    env = os.environ.copy()
    env['HWLOC_COMPONENTS'] = '-opencl'
    env['OMP_NUM_THREADS'] = '1'

    def run(name, input_deck, mpi=False, extra=()):
        output = root / name
        output.mkdir()
        run_env = env.copy()
        if mpi:
            run_env.pop('OMPI_MCA_btl', None)
            prefix = [args.mpi_exec, args.mpi_numproc_flag, '2']
        else:
            run_env['OMPI_MCA_btl'] = 'self'
            prefix = []
        argv = prefix + [args.flow, input_deck, '--output-dir=' + str(output)] + OPTIONS + list(extra)
        command(argv, output / 'console.log', run_env)
        return output

    flow_hash = hashlib.sha256(args.flow.read_bytes()).hexdigest()
    baseline = run('continuous', deck)
    rows = read_summary(args.summary, baseline)
    validate(rows, (baseline / 'console.log').read_text(), False)
    (root / 'provenance.json').write_text(json.dumps({
        'flow': str(args.flow), 'sha256': flow_hash,
        'options': OPTIONS, 'absolute': .02, 'relative': 1e-5}, indent=2))

    def compare(output):
        command([args.compare, '-d', '-t', 'SMRY', baseline / 'GECON-02', output / 'GECON-02',
                 '.02', '1e-5'], output / 'compare-summary.log')
        command([args.compare, '-t', 'UNRST', baseline / 'GECON-02',
                 output / 'GECON-02', '.02', '1e-5'], output / 'compare-restart.log')

    if args.mode == 'mpi':
        parallel = run('mpi2', deck, mpi=True)
        validate(read_summary(args.summary, parallel),
                 (parallel / 'console.log').read_text(), False)
        compare(parallel)
    elif args.mode == 'unrst':
        # Reproduce the lost tiny-phase / unconfigured economic closure defect.
        # Compare the complete first post-restart report. Later PID trajectories
        # differ because UNRST files omit controller history.
        input_dir = root / 'input-resume'
        input_dir.mkdir()
        restarted_deck = input_dir / 'GECON-02.DATA'
        command([args.rst_deck, '-m', 'inline', '-s', deck,
                 str(baseline / 'GECON-02.UNRST') + ':4', restarted_deck],
                input_dir / 'rst-deck.log')
        output = run('unrst-resume', restarted_deck, extra=['--sched-restart=true'])
        resumed = read_summary(args.summary, output)
        require(bool(resumed) and resumed[-1]['TIME'] == 213, 'UNRST restart did not finish')
        events = re.findall(r'Well (\S+) will be shut because:',
                            (output / 'console.log').read_text())
        require(events == ['B-1H'], 'UNRST restart lost GECON closure state: ' + str(events))
        closure = re.search(r'At time = ([\d.]+) DAYS[^\n]*Well B-1H will be shut',
                            (output / 'console.log').read_text())
        require(closure is not None and abs(float(closure[1]) - 67.86) < .125,
                'UNRST restart changed second economic closure timing')
        require(all(row['WGPR:B-3H'] == 0 and row['WWPR:B-3H'] == 0
                    for row in resumed), 'UNRST restart reopened B-3H')
        require(abs(resumed[-1]['FWPT'] - rows[-1]['FWPT']) < 1.5,
                'UNRST restart changed water production outside the economic envelope')
        command([args.compare, '-r', '5', '-t', 'UNRST', baseline / 'GECON-02',
                 output / 'GECON-02', '.02', '1e-5'], output / 'compare-first-report.log')
    require(hashlib.sha256(args.flow.read_bytes()).hexdigest() == flow_hash,
            'Flow executable changed during validation')
    (root / 'result.json').write_text(json.dumps({'mode': args.mode, 'result': 'PASS',
        'flow_sha256': flow_hash, 'absolute': .02, 'relative': 1e-5,
        'scope': {'mpi': 'all report summaries and all restart fields',
                  'unrst': 'first post-restart report and economic closure semantics'}[args.mode]}, indent=2))
    print('GECON {} consistency passed'.format(args.mode), flush=True)


if __name__ == '__main__':
    main()
