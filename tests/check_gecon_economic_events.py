#!/usr/bin/env python3
"""Run GECON-02 with resolved timesteps and check economic-event semantics.

Uses only Python's standard library and opm-common's summary executable.
Checks closure times within 0.125 days of EXPECTED_CLOSURES and day-213
water production within 316.0 +/- 1.5 m3, plus closure semantics and, for
adaptive runs, accepted normalized WGR changes and initial-step agreement.
"""
import argparse
import math
from pathlib import Path
import re
import subprocess
import sys
import tempfile

WELLS = ('B-1H', 'B-2H', 'B-3H', 'C-1H', 'C-2H')
VECTORS = ('TIME', 'FWPR', 'FGPR', 'FWPT') + tuple(
    '{}:{}'.format(vector, well) for well in WELLS for vector in ('WWPR', 'WGPR'))
STEP_DAYS = 0.0625
ENDPOINTS = tuple(i * STEP_DAYS for i in range(1, 1457)) + (121, 152, 182, 213)
EXPECTED_CLOSURES = (('B-3H', 23.25), ('B-1H', 67.875))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_summary(executable, output):
    result = subprocess.run(
        [str(executable), str(output / 'GECON-02.SMSPEC'), *VECTORS],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
        check=True, timeout=60)
    rows = []
    for line in result.stdout.splitlines():
        fields = line.split()
        if not fields or fields[0] == 'TIME':
            continue
        require(len(fields) == len(VECTORS), 'Malformed summary row: ' + line)
        values = tuple(map(float, fields))
        require(all(math.isfinite(value) for value in values),
                'Non-finite summary value: ' + line)
        rows.append(dict(zip(VECTORS, values)))
    return rows


def validate(rows, log, prescribed_schedule=True, wgr_tolerance=0.001):
    require(bool(rows), 'Missing summary timesteps')
    require(rows[-1]['TIME'] == 213, 'Simulation did not reach day 213')
    require(all(a['TIME'] < b['TIME'] for a, b in zip(rows, rows[1:])),
            'Summary times do not increase')
    if prescribed_schedule:
        require(len(rows) == len(ENDPOINTS), 'Missing or unexpected summary timesteps')
        require(all(abs(row['TIME'] - endpoint) < 1e-8
                    for row, endpoint in zip(rows, ENDPOINTS)),
                'Summary endpoints differ from the prescribed schedule')
        require('Timestep chopped to' not in log, 'Unexpected timestep chop')
    else:
        require('Rejecting timestep for group economic WGR change' in log,
                'Economic accuracy rejection path was not exercised')
    events = re.findall(
        r'At time = ([\d.]+) DAYS[^\n]*Well (\S+) will be shut because:\n([^\n]+)', log)
    require(len(events) == len(EXPECTED_CLOSURES),
            'Expected exactly two economic shut-ins, got {}'.format(events))
    closures = []
    for (printed_time, well, reason), (expected_well, expected_time) in zip(
            events, EXPECTED_CLOSURES):
        require(well == expected_well, 'Wrong well/order: ' + well)
        require('Water-gas ratio for group FIELD' in reason,
                'Unexpected economic-limit reason: ' + reason)
        index = min(range(len(rows)), key=lambda i: abs(rows[i]['TIME'] - float(printed_time)))
        at = rows[index]
        require(abs(at['TIME'] - float(printed_time)) <= 0.0051,
                'Logged shut-in does not match a summary endpoint')
        require(abs(at['TIME'] - expected_time) <= 0.125,
                '{} shut at day {}, expected {} +/- 0.125'.format(
                    well, at['TIME'], expected_time))
        require(0 < index < len(rows) - 1, 'Missing event neighbours')
        before, after = rows[index - 1], rows[index + 1]
        require(before['FGPR'] > 0 and at['FGPR'] > 0, 'Non-positive group gas rate')
        require(before['FWPR'] / before['FGPR'] <= 1e-6 < at['FWPR'] / at['FGPR'],
                'FIELD WGR does not cross its 1e-6 limit at shut-in')
        ratios = {w: at['WWPR:' + w] / at['WGPR:' + w]
                  for w in WELLS if at['WGPR:' + w] > 0}
        require(bool(ratios) and max(ratios, key=ratios.get) == well,
                'Closed well is not the worst producing-well WGR')
        require(at['WWPR:' + well] > 0 and at['WGPR:' + well] > 0,
                'Selected well was not producing at the completed step')
        # A group-limit closure affects subsequent steps, not the completed one.
        require(all(row['WWPR:' + well] == 0 and row['WGPR:' + well] == 0
                    for row in rows[index + 1:]),
                'Closed well continues producing or reopens unexpectedly')
        closures.append(at['TIME'])
        print('{} shut at day {}; FIELD WGR={:.9g}'.format(
            well, at['TIME'], at['FWPR'] / at['FGPR']))
    require(abs(rows[-1]['FWPT'] - 316.0) <= 1.5,
            'Day-213 water production {} outside 316.0 +/- 1.5 m3'.format(rows[-1]['FWPT']))
    print('Day-213 water production: {} m3'.format(rows[-1]['FWPT']))
    if not prescribed_schedule:
        # Verify the actual accepted changes, not just the presence of a retry
        # marker. A closure belongs to the preceding completed step; remove its
        # rates from that step before comparing against the next one.
        closed = set()
        maximum_change = 0.0
        for before, after in zip(rows, rows[1:]):
            closed.update(w for w in WELLS
                          if before['WGPR:' + w] > 0 and after['WGPR:' + w] == 0)
            water = before['FWPR'] - sum(before['WWPR:' + w] for w in closed)
            gas = before['FGPR'] - sum(before['WGPR:' + w] for w in closed)
            require(gas > 0 and after['FGPR'] > 0, 'Missing gas in adaptive WGR check')
            change = abs(after['FWPR'] / after['FGPR'] - water / gas) / 1e-6
            maximum_change = max(maximum_change, change)
            # summary writes floating-point rates with limited precision.
            require(change <= wgr_tolerance + 2e-6,
                    'Accepted normalized WGR change {} exceeds {} at day {}'.format(
                        change, wgr_tolerance, after['TIME']))
        print('Maximum accepted normalized WGR change: {}'.format(maximum_change))
    return closures, rows[-1]['FWPT']


def check_reopening(args):
    output = Path(tempfile.mkdtemp(prefix='reopen-', dir=str(args.output.resolve())))
    print('GECON reopening artifacts: {}'.format(output), flush=True)
    text = args.deck.read_text()
    text = re.sub(r"(INCLUDE\s*\n\s*)'([^']+)'",
                  lambda match: match[1] + "'" + str((args.deck.parent / match[2]).resolve()) + "'",
                  text, flags=re.IGNORECASE)
    marker = "  01 'FEB' 2020 /"
    require(text.count(marker) == 1, 'Cannot locate day-31 GECON report boundary')
    text = text.replace(marker, marker + "\n/\nWELOPEN\n 'B-3H' 'OPEN' /\n/\nDATES")
    deck = output / 'GECON-02.DATA'
    deck.write_text(text)
    with (output / 'console.log').open('w') as log:
        subprocess.run(
            [str(args.flow), str(deck), '--output-dir=' + str(output),
             '--threads-per-process=1', '--initial-time-step-in-days=1',
             '--group-economic-wgr-step-tolerance=0.001'],
            stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
    rows = read_summary(args.summary, output)
    log = (output / 'console.log').read_text()
    require(bool(rows) and rows[-1]['TIME'] == 213, 'Reopening run did not finish')
    events = re.findall(r'At time = ([\d.]+) DAYS[^\n]*Well (\S+) will be shut', log)
    require([well for _, well in events] == ['B-3H', 'B-3H', 'B-1H'],
            'Unexpected reopening closure sequence: ' + str(events))
    require('Skipping group economic WGR change check for changed well-test status in group FIELD' in log,
            'Reopening did not exercise discrete-status handling')
    index = min(range(len(rows)), key=lambda i: abs(rows[i]['TIME'] - float(events[1][0])))
    at = rows[index]
    require(31 < at['TIME'] < 60 and at['WGPR:B-3H'] > 0 and at['WWPR:B-3H'] > 0,
            'Scheduled reopening did not restore production')
    require(at['FWPR'] / at['FGPR'] > 1e-6, 'Reclosed group does not violate WGR')
    require(all(row['WGPR:B-3H'] == 0 and row['WWPR:B-3H'] == 0 for row in rows[index + 1:]),
            'Reclosed well continues producing')
    print('Scheduled reopening and economic reclosure passed')


def check_control_change(args):
    """A later WCONPROD edit must not pin retries to the old solved rates."""
    output = Path(tempfile.mkdtemp(prefix='control-change-', dir=str(args.output.resolve())))
    print('GECON control-change artifacts: {}'.format(output), flush=True)
    text = args.deck.read_text()
    text = re.sub(r"(INCLUDE\s*\n\s*)'([^']+)'",
                  lambda match: match[1] + "'" + str((args.deck.parent / match[2]).resolve()) + "'",
                  text, flags=re.IGNORECASE)
    marker = "  20 'JAN' 2020 /"
    require(text.count(marker) == 1, 'Cannot locate day-19 GECON report boundary')
    # Reduce the wet producer's target by a factor of ten. Keep this run short:
    # only the first accepted step after the edit is relevant to this regression.
    text = text[:text.index(marker)] + marker + (
        "\n/\nWCONPROD\n 'B-3H' 'OPEN' 'GRAT' 2* 3.0E5 2* 100 /\n/\n"
        "TSTEP\n 1 /\nEND\n")
    deck = output / 'GECON-02.DATA'
    deck.write_text(text)
    with (output / 'console.log').open('w') as log:
        subprocess.run(
            [str(args.flow), str(deck), '--output-dir=' + str(output),
             '--threads-per-process=1', '--initial-time-step-in-days=1',
             '--group-economic-wgr-step-tolerance=0.001'],
            stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
    rows = read_summary(args.summary, output)
    require(bool(rows) and abs(rows[-1]['TIME'] - 20) < 1e-8,
            'Control-change run did not finish')
    before = [row for row in rows if abs(row['TIME'] - 19) < 1e-8]
    after = [row for row in rows if row['TIME'] > 19]
    require(len(before) == 1 and bool(after), 'Missing control-change neighbours')
    require(before[0]['WGPR:B-3H'] > 1e6,
            'Producer was not on its original gas target before the edit')
    require(0 < after[0]['WGPR:B-3H'] <= 3.0e5 + 1,
            'Scheduled WCONPROD target edit did not take effect')
    require(after[0]['WWPR:B-3H'] > 0,
            'Control-change regression did not exercise a wet producer')
    require(before[0]['FGPR'] > 0 and after[0]['FGPR'] > 0,
            'Missing group gas rate around the target edit')
    jump = abs(after[0]['FWPR'] / after[0]['FGPR']
               - before[0]['FWPR'] / before[0]['FGPR']) / 1e-6
    require(jump > 0.001,
            'Target edit did not exceed the stale-baseline WGR tolerance')
    print('Scheduled WCONPROD change and first-step acceptance passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--flow', type=Path)
    parser.add_argument('--summary', type=Path, required=True)
    parser.add_argument('--deck', type=Path)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--adaptive-control', action='store_true',
                        help='Exercise WGR accuracy control with adaptive initial steps of 1 and 8 days')
    parser.add_argument('--check-only', action='store_true',
                        help='Validate an existing run without starting Flow')
    args = parser.parse_args()
    if args.adaptive_control:
        if args.check_only or not (args.flow and args.deck):
            parser.error('--adaptive-control requires --flow and --deck, without --check-only')
        args.output.mkdir(parents=True, exist_ok=True)
        results = []
        for initial in (1, 8):
            output = Path(tempfile.mkdtemp(prefix='adaptive-{}-'.format(initial),
                                          dir=str(args.output.resolve())))
            print('GECON artifacts: {}'.format(output), flush=True)
            with (output / 'console.log').open('w') as log:
                subprocess.run(
                    [str(args.flow), str(args.deck.resolve()), '--output-dir=' + str(output),
                     '--threads-per-process=1', '--initial-time-step-in-days=' + str(initial),
                     '--group-economic-wgr-step-tolerance=0.001'],
                    stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
            results.append(validate(read_summary(args.summary, output),
                                    (output / 'console.log').read_text(), False))
        require(all(abs(a - b) <= 0.0625 for a, b in zip(results[0][0], results[1][0])),
                'Closure timing depends on the initial timestep by more than 0.0625 days')
        require(abs(results[0][1] - results[1][1]) <= 0.5,
                'Water production depends on the initial timestep by more than 0.5 m3')
        check_reopening(args)
        check_control_change(args)
        print('Adaptive GECON WGR accuracy regression passed')
        return
    if not args.check_only:
        if not (args.flow and args.deck):
            parser.error('--flow and --deck are required')
        args.output.mkdir(parents=True, exist_ok=True)
        args.output = Path(tempfile.mkdtemp(prefix='run-', dir=str(args.output.resolve())))
        print('GECON artifacts: {}'.format(args.output), flush=True)
        steps = args.output / 'timesteps'
        steps.write_text('\n'.join(map(str, ENDPOINTS)) + '\n')
        with (args.output / 'console.log').open('w') as log:
            subprocess.run(
                [str(args.flow), str(args.deck.resolve()),
                 '--output-dir=' + str(args.output.resolve()), '--threads-per-process=1',
                 '--time-step-control=hardcoded', '--time-step-control-file-name=' + str(steps)],
                stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
    validate(read_summary(args.summary, args.output), (args.output / 'console.log').read_text())
    print('GECON economic-event regression passed')


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        print('GECON economic-event regression FAILED: {}'.format(error), file=sys.stderr)
        sys.exit(1)
