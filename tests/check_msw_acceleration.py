#!/usr/bin/env python3
"""Check MSW segment pressures against a momentum balance, for both flow directions.

An injector and a producer have identical horizontal pipes, so the only pressure
terms are acceleration and friction. Fluid leaving through the wall keeps its axial
momentum, so p + rho v^2 / 2 is conserved in the injector. Fluid entering has none,
so p + rho v^2 is conserved in the producer. INJB and PRODB split into two branches
at segment 3. Where flow divides, the branches enter as their flow-weighted mean.
PRODX is a producer whose segment 3 loses fluid through the wall.
"""
import argparse
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys

DENSITY = 1000.0  # kg/m3, tests/msw_acceleration.DATA
DIAMETER = 0.1  # m
TOLERANCE = 1e-3  # bar
TARGET_RATE = 10000.0  # sm3/day, tests/msw_acceleration.DATA
RATE_TOLERANCE = 10.0  # sm3/day

# Segment: (outlet, inlets). Segment 2 is blank.
PIPE = {2: (1, [3]), 3: (2, [4]), 4: (3, [])}
BRANCHED = {2: (1, [3]), 3: (2, [4, 5]), 4: (3, []), 5: (3, [])}
# Producer segments with net outflow through the wall.
OUTFLOW_SEGMENTS = {'PRODX': (3,)}


def momentum_flux(rate):
    """rho v^2 in bar for a rate in sm3/day; Bw = 1."""
    area = math.pi * DIAMETER * DIAMETER / 4
    return DENSITY * (rate / 86400.0 / area) ** 2 / 1e5


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('flow', 'summary', 'deck', 'output'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()

    if args.output.exists():
        shutil.rmtree(args.output)
    args.output.mkdir(parents=True)
    env = dict(os.environ, OMP_NUM_THREADS='1')
    subprocess.run([args.flow, f'--output-dir={args.output}', args.deck],
                   env=env, check=True, stdout=subprocess.DEVNULL)

    smspec = args.output / (args.deck.stem.upper() + '.SMSPEC')
    failed = False
    for well, coefficient, segments in (('INJ', 0.5, PIPE), ('PROD', 1.0, PIPE),
                                        ('INJB', 0.5, BRANCHED), ('PRODB', 1.0, BRANCHED),
                                        ('PRODX', 1.0, PIPE)):
        keys = [f'{v}:{well}:{s}' for v in ('SPR', 'SPRDF', 'SWFR') for s in (1, *segments)]
        out = subprocess.run([args.summary, smspec, *keys],
                             env=env, check=True, capture_output=True, text=True).stdout
        last = [float(x) for x in out.strip().splitlines()[-1].split()]
        val = dict(zip(keys, last[-len(keys):]))

        # A shut well has zero flux and friction, so it satisfies every balance below.
        rate = abs(val[f'SWFR:{well}:2'])
        print(f'{well} rate: {rate:.1f} target {TARGET_RATE:.1f}')
        if abs(rate - TARGET_RATE) > RATE_TOLERANCE:
            print(f'{well}: does not deliver its target rate')
            failed = True

        for seg, (outlet, inlets) in segments.items():
            rates = [abs(val[f'SWFR:{well}:{inlet}']) for inlet in inlets]
            inlet_flux = sum(momentum_flux(r) for r in rates)
            if coefficient < 1.0 and len(inlets) > 1:
                inlet_flux = sum(r * momentum_flux(r) for r in rates) / sum(rates)
            flux = momentum_flux(abs(val[f'SWFR:{well}:{seg}'])) - inlet_flux
            friction = val[f'SPRDF:{well}:{seg}']
            segment_coefficient = 0.5 if seg in OUTFLOW_SEGMENTS.get(well, ()) else coefficient
            expected = segment_coefficient * flux + friction
            measured = val[f'SPR:{well}:{seg}'] - val[f'SPR:{well}:{outlet}']
            print(f'{well} {seg}-{outlet}: measured {measured:.4f} expected {expected:.4f}'
                  f' flux {flux:.4f} friction {friction:.4f}')
            if abs(measured - expected) > TOLERANCE:
                failed = True
            # The balance only tests the acceleration term if it outweighs friction.
            if seg != 2 and abs(segment_coefficient * flux) < 10 * abs(friction):
                print(f'{well} {seg}: acceleration does not dominate')
                failed = True
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
