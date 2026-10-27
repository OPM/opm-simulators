#!/usr/bin/env python3
"""Check explicit well-level reopening at the first report step after UNRST restart.

A SPE1 producer is economically closed on day 1, remains closed through day 2,
then receives WELOPEN OPEN with its economic limit disabled. Restart at day 2
must honor that event just as the continuous run does. A long WTEST interval
preserves the closure on common master without triggering automatic reopening.
Uses only Python's standard library and the common summary/rst_deck utilities.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('flow', 'summary', 'rst-deck', 'deck', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve() / 'reopening'
    if root.exists():
        shutil.rmtree(root)
    root.mkdir(parents=True)
    print('Well reopening restart artifacts: ' + str(root), flush=True)
    commands = []
    env = dict(os.environ, OMP_NUM_THREADS='1', HWLOC_COMPONENTS='-opencl',
               OMPI_MCA_btl='self')

    def command(argv, log):
        commands.append(list(map(str, argv)))
        (root / 'commands.json').write_text(json.dumps(commands, indent=2))
        with log.open('w') as stream:
            subprocess.run(commands[-1], stdout=stream, stderr=subprocess.STDOUT,
                           env=env, check=True, timeout=120)

    text = args.deck.read_text()
    matches = list(re.finditer(r'^TSTEP\s*$', text, re.MULTILINE))
    require(len(matches) == 1, 'Expected the SPE1 fixture to have one TSTEP section')
    text = text[:matches[0].start()] + """WTEST
 'PROD' 1000 'E' 1 /
/
WECON
 'PROD' 1e10 4* 'WELL' /
/
TSTEP
 1 1 /
WECON
 'PROD' 0 4* 'NONE' /
/
WELOPEN
 'PROD' 'OPEN' /
/
TSTEP
 1 /
END
"""
    deck = root / 'SPE1CASE1.DATA'
    deck.write_text(text)
    digest = hashlib.sha256(args.flow.read_bytes()).hexdigest()

    def run(name, deck, extra=()):
        output = root / name
        output.mkdir()
        command([args.flow, deck, '--output-dir=' + str(output),
                 '--threads-per-process=1', '--enable-opm-rst-file=true',
                 '--ecl-output-double-precision=true', *extra], output / 'console.log')
        data = subprocess.check_output(
            [str(args.summary), str(output / (deck.stem + '.SMSPEC')),
             'TIME', 'WOPR:PROD', 'WOPT:PROD'], text=True, env=env, timeout=30)
        (output / 'summary.txt').write_text(data)
        rows = [list(map(float, line.split())) for line in data.splitlines()
                if line.split() and line.split()[0] != 'TIME']
        require(bool(rows) and rows[-1][0] == 3, name + ' did not reach day 3')
        return output, rows

    baseline, continuous = run('continuous', deck)
    day2 = next(row for row in continuous if row[0] == 2)
    require(day2[1] == 0 and day2[2] > 0, 'Producer must be closed at the checkpoint')
    require(continuous[-1][1] > 0 and continuous[-1][2] > day2[2],
            'Continuous run did not reopen and produce after day 2')
    restarted_deck = root / 'RESUME.DATA'
    command([args.rst_deck, '-m', 'inline', '-s', deck,
             str(baseline / 'SPE1CASE1.UNRST') + ':2', restarted_deck],
            root / 'rst-deck.log')
    resumed_output, resumed = run('resumed', restarted_deck, ['--sched-restart=true'])
    require(all(row[0] > 2 for row in resumed), 'Unexpected pre-restart output')
    require(resumed[-1][1] > 0 and resumed[-1][2] > day2[2],
            'Restored closure blocked the explicit OPEN request at day 2')
    for index, name in ((1, 'oil rate'), (2, 'cumulative oil')):
        before, after = continuous[-1][index], resumed[-1][index]
        require(abs(before - after) <= .02 + 1e-5 * max(abs(before), abs(after)),
                'Continuous/restarted ' + name + ' differs')
    # A restored closure must remain closed when no explicit OPEN is requested.
    # This prevents fixing reopening by simply dropping all restored closures.
    no_open_source = root / 'NOOPEN.DATA'
    no_open_source.write_text(text.replace("WELOPEN\n 'PROD' 'OPEN' /\n/\n", ""))
    no_open_deck = root / 'REMAIN_CLOSED.DATA'
    command([args.rst_deck, '-m', 'inline', '-s', no_open_source,
             str(baseline / 'SPE1CASE1.UNRST') + ':2', no_open_deck],
            root / 'rst-deck-no-open.log')
    _, still_closed = run('no-open', no_open_deck, ['--sched-restart=true'])
    require(all(row[1] == 0 for row in still_closed),
            'Restored well reopened without an explicit OPEN request')
    require(abs(still_closed[-1][2] - day2[2]) <= .02,
            'Closed well accumulated production without OPEN')
    require(hashlib.sha256(args.flow.read_bytes()).hexdigest() == digest,
            'Flow executable changed during the test')
    (root / 'result.json').write_text(json.dumps({
        'result': 'PASS', 'flow_sha256': digest, 'restart_day': 2,
        'closed_checkpoint': day2, 'continuous_day3': continuous[-1],
        'resumed_day3': resumed[-1], 'no_open_day3': still_closed[-1],
        'scope': 'Explicit WELOPEN OPEN at first restarted report step; production resumes'
    }, indent=2) + '\n')
    print('Restart-boundary explicit well reopening passed', flush=True)


if __name__ == '__main__':
    main()
