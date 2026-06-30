#!/usr/bin/env python3

# Generates a PDF with plots of all summary curves from a reference
# case and a 'new' simulation.

import argparse
from filelock import FileLock
from functools import wraps
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from opm.io.ecl import ESmry
import os
import pickle
from scipy import integrate, stats
import sys


def with_file_lock(method):
    @wraps(method)
    def wrapper(self, *args, **kwargs):
        with FileLock(self.lock_file):
            return method(self, *args, **kwargs)

    return wrapper


class WellComparisonManager:
    def __init__(self, max_devs_file='max_devs.pkl', lock_file='max_devs.lck'):
        self.max_devs_file = max_devs_file
        self.lock_file = lock_file
        plt.rcParams['font.size'] = 8

    @staticmethod
    def _is_curve_key(key):
        return key not in ('TIME', 'YEARS')

    @staticmethod
    def _has_nonzero_data(values):
        return len(values) > 0 and np.any(values)

    def _write_pdf(self, test_name, plot_entries):
        with PdfPages(f'{test_name}.pdf') as pdf:
            for entry in plot_entries:
                fig, ax = plt.subplots()
                for series in entry['series']:
                    ax.plot(series['time'], series['values'], **series['style'])

                ax.legend(entry['legend'])
                plt.title(entry['name'])

                if entry['unit']:
                    plt.ylabel(entry['unit'])

                myFmt = DateFormatter("%Y-%b")
                ax.xaxis.set_major_formatter(myFmt)
                ax.xaxis.set_major_locator(plt.MaxNLocator(20))
                plt.grid()
                fig.autofmt_xdate()
                fig.savefig(pdf, format='pdf')
                plt.close(fig)

    @with_file_lock
    def _update_max_deviation(self, test_name, deviation_value):
        if os.path.exists(self.max_devs_file):
            try:
                with open(self.max_devs_file, 'rb') as f:
                    max_deviations = pickle.load(f)
            except (EOFError, pickle.UnpicklingError, OSError) as err:
                print(
                    f"Cannot read ranking data from {self.max_devs_file}: {err}",
                    file=sys.stderr,
                )
                return
        else:
            max_deviations = {}

        max_deviations[test_name] = deviation_value

        with open(self.max_devs_file, 'wb') as f:
            pickle.dump(max_deviations, f)

    # Run analysis of a test.
    # Calculate the deviation for each curve
    # and generate a pdf with plots ordered according to deviation
    def run_analysis(self, ref_smspec, sim_smspec, test_name, ref_name, sim_name):
        ref_file = ESmry(ref_smspec)
        sim_file = ESmry(sim_smspec)

        ref_time = ref_file.dates()
        sim_time = sim_file.dates()

        ref_time_in_secs = [(v - ref_time[0]).total_seconds() for v in ref_time]
        sim_time_in_secs = [(v - sim_time[0]).total_seconds() for v in sim_time]

        # Attempt at sorting the graphs in descending eyeball norm order.
        # - Normalize by inf-norm to get the same range in each graph, ie (0,1).
        # - Convert graphs to probability distributions (ie integral under curve should be 1).
        # - Use the wasserstein distance scaled by area under reference curve.
        deviation = {}
        for curve_name in ref_file.keys():
            if not self._is_curve_key(curve_name):
                continue

            try:
                ref = ref_file[curve_name]
                sim = sim_file[curve_name]
            except Exception:
                continue

            if len(ref) == 0 and len(sim) == 0:
                continue
            if len(ref) == 0 or len(sim) == 0:
                deviation[curve_name] = float('inf')
                continue

            if not (np.any(ref) or np.any(sim)):
                continue

            ref_norm = np.linalg.norm(ref, np.inf)
            sim_norm = np.linalg.norm(sim, np.inf)
            if ref_norm == 0.0 or sim_norm == 0.0:
                deviation[curve_name] = float('inf')
                continue

            ref_scaled = ref / ref_norm
            sim_scaled = sim / sim_norm

            A_ref = integrate.trapezoid(ref_scaled, ref_time_in_secs, 0.0)
            A_sim = integrate.trapezoid(sim_scaled, sim_time_in_secs, 0.0)

            if A_ref == 0.0 and A_sim == 0.0:
                deviation[curve_name] = 0.0
                continue

            if A_ref == 0.0 or A_sim == 0.0:
                deviation[curve_name] = float('inf')
                continue

            deviation[curve_name] = stats.wasserstein_distance(ref_scaled / A_ref, sim_scaled / A_sim) * A_ref

        plot_entries = []
        for curve_name in sorted(deviation, key=lambda x: deviation[x], reverse=True):
            try:
                ref = ref_file[curve_name]
                sim = sim_file[curve_name]
            except Exception:
                continue

            if len(ref) == 0 or len(sim) == 0:
                continue

            plot_entries.append(
                {
                    'name': curve_name,
                    'unit': ref_file.units(curve_name),
                    'legend': [ref_name, sim_name],
                    'series': [
                        {
                            'time': ref_time,
                            'values': ref,
                            'style': {
                                'linestyle': 'dashed',
                                'linewidth': 0.5,
                                'marker': 'o',
                                'markersize': 1.0,
                            },
                        },
                        {
                            'time': sim_time,
                            'values': sim,
                            'style': {
                                'linewidth': 0.5,
                                'marker': 'x',
                                'markersize': 1.0,
                            },
                        },
                    ],
                }
            )

        self._write_pdf(test_name, plot_entries)

        # Deviation is infinite in the following cases:
        # - One of the curves is empty and the other is not.
        # - One of the curves is zero and the other is not.
        #
        # This is intentional and ensures that new cases are
        # ranked above (reverse=True in sorted() call) any
        # case with a reference solution.
        max_deviation = max(deviation.values(), default=0.0)
        self._update_max_deviation(test_name, max_deviation)

    # Analyze a case where no reference solution exists.
    # Generate plots for all non-empty curves from the simulation only.
    def analyze_new_case(self, sim_smspec, test_name, sim_name):
        sim_file = ESmry(sim_smspec)
        sim_time = sim_file.dates()

        plot_entries = []
        for curve_name in sim_file.keys():
            if not self._is_curve_key(curve_name):
                continue

            try:
                sim = sim_file[curve_name]
            except Exception:
                continue

            if not self._has_nonzero_data(sim):
                continue

            plot_entries.append(
                {
                    'name': curve_name,
                    'unit': sim_file.units(curve_name),
                    'legend': [sim_name],
                    'series': [
                        {
                            'time': sim_time,
                            'values': sim,
                            'style': {
                                'linewidth': 0.5,
                                'marker': 'x',
                                'markersize': 1.0,
                            },
                        }
                    ],
                }
            )

        self._write_pdf(test_name, plot_entries)

        # No reference exists, so register zero deviation for ranking.
        self._update_max_deviation(test_name, 0.0)

    # Rename files to rank them according to maximum deviations
    @with_file_lock
    def reorder_files(self):
        if not os.path.exists(self.max_devs_file):
            print(
                f"Case ranking data file {self.max_devs_file} does not exist",
                file=sys.stderr,
            )
            return False

        try:
            with open(self.max_devs_file, 'rb') as f:
                max_deviations = pickle.load(f)
        except (EOFError, pickle.UnpicklingError, OSError) as err:
            print(
                f"Cannot read ranking data from {self.max_devs_file}: {err}",
                file=sys.stderr,
            )
            return False

        c = 1
        for test in sorted(max_deviations, key=lambda x: max_deviations[x], reverse=True):
            src = f'{test}.pdf'
            if not os.path.exists(src):
                continue

            os.replace(src, f'{c:02d}_{test}.pdf')
            c += 1

        return True

# Main code
parser = argparse.ArgumentParser('plot_well_comparison.py')

parser.add_argument('-c', help='Name of test to process', dest='test_name')
parser.add_argument('-r', help='Reference file', dest='ref_file')
parser.add_argument('-s', help='Simulation file', dest='sim_file')
parser.add_argument('-t', help='Name for first simulation', dest='ref_name', default='Reference')
parser.add_argument('-u', help='Name for second simulation', dest='sim_name', default='New simulation')
parser.add_argument('-o', choices=['plot', 'rename'], help='Operation to do', required=True, dest='operation')
args = parser.parse_args()
manager = WellComparisonManager()

if args.operation == 'plot':
    print(f"Processing {args.test_name}")
    sim_smspec = args.sim_file + '.SMSPEC'
    ref_smspec = (args.ref_file + '.SMSPEC') if args.ref_file else None

    if not os.path.exists(sim_smspec):
        print(f"Cannot process {args.test_name}: missing simulation summary {sim_smspec}", file=sys.stderr)
        sys.exit(1)

    if ref_smspec and os.path.exists(ref_smspec):
        manager.run_analysis(ref_smspec, sim_smspec, args.test_name, args.ref_name, args.sim_name)
    else:
        manager.analyze_new_case(sim_smspec, args.test_name, args.sim_name)
else:
    if not manager.reorder_files():
        sys.exit(1)
