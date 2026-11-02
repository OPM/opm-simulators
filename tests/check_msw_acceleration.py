#!/usr/bin/env python3
"""Check MSW segment pressures against a momentum balance, for both flow directions.

An injector and a producer have identical horizontal pipes, so the only pressure
terms are acceleration and friction.
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
    for well in ('INJ', 'PROD'):
        keys = [f'{v}:{well}:{s}' for v in ('SPR', 'SPRDF', 'SWFR') for s in (1, 2, 3, 4)]
        out = subprocess.run([args.summary, smspec, *keys],
                             env=env, check=True, capture_output=True, text=True).stdout
        last = [float(x) for x in out.strip().splitlines()[-1].split()]
        val = dict(zip(keys, last[-len(keys):]))

        # Segment 4 is the closed end, segment 2 is blank.
        for seg, inlet in ((4, None), (3, 4), (2, 3)):
            flux = momentum_flux(abs(val[f'SWFR:{well}:{seg}']))
            if inlet is not None:
                flux -= momentum_flux(abs(val[f'SWFR:{well}:{inlet}']))
            friction = val[f'SPRDF:{well}:{seg}']
            expected = flux + friction
            measured = val[f'SPR:{well}:{seg}'] - val[f'SPR:{well}:{seg - 1}']
            print(f'{well} {seg}-{seg - 1}: measured {measured:.4f} expected {expected:.4f}'
                  f' flux {flux:.4f} friction {friction:.4f}')
            if abs(measured - expected) > TOLERANCE:
                failed = True
            # Keeps the check sensitive to the acceleration term.
            if seg != 2 and abs(flux) < 10 * abs(friction):
                print(f'{well} {seg}: acceleration does not dominate')
                failed = True
    return 1 if failed else 0


if __name__ == '__main__':
    sys.exit(main())
