#!/usr/bin/env python3

# Generates a PDF with plots of all summary curves from a reference
# case and a 'new' simulation.

import argparse
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from opm.io.ecl import ESmry
import os
import pickle
from scipy import integrate, stats

# Run analysis of a test.
# Calculate the deviation for each curve
# and generate a pdf with plots ordered according to deviation
def run_analysis(ref_file_name, sys_file_name, test_name, ref_name, sim_name):
    ref_file = ESmry(ref_file_name + '.SMSPEC')
    sim_file = ESmry(sys_file_name + '.SMSPEC')

    ref_time = ref_file.dates()
    sim_time = sim_file.dates()

    ref_time_in_secs = [(v - ref_time[0]).total_seconds() for v in ref_time]
    sim_time_in_secs = [(v - sim_time[0]).total_seconds() for v in sim_time]

    plt.rcParams['font.size'] = 8

    # Attempt at sorting the graphs in descending eyeball norm order.
    # - Normalize by inf-norm to get the same range in each graph, ie (0,1).
    # - Convert graphs to probability distributions (ie integral under curve should be 1).
    # - Use the wasserstein distance scaled by area under reference curve.
    deviation={}
    for r in ref_file.keys():
        if r == 'TIME' or r == 'YEARS':
            continue
        try:
            ref = ref_file[r]
            sim = sim_file[r]
        except:
            continue

        if len(ref) == 0 and len(sim) == 0:
            continue

        if not (any(ref) or any(sim)):
            continue

        ref /= np.linalg.norm(ref, np.inf)
        sim /= np.linalg.norm(sim, np.inf)

        A_ref = integrate.trapezoid(ref, ref_time_in_secs, 0.0)
        A_sim = integrate.trapezoid(sim, sim_time_in_secs, 0.0)

        deviation[r] = stats.wasserstein_distance(ref / A_ref, sim / A_sim) * A_ref

    p = PdfPages(f'{test_name}.pdf')
    for r in sorted(deviation, key = lambda x: deviation[x], reverse=True):
        try:
            ref = ref_file[r]
            sim = sim_file[r]
        except:
            continue

        fig, ax = plt.subplots()
        ax.plot(ref_time, ref, linestyle='dashed', linewidth=0.5, marker='o', markersize=1.0)
        ax.plot(sim_time, sim, linewidth=0.5, marker='x', markersize=1.0)
        ax.legend([ref_name, sim_name])
        plt.title(r)
        u = ref_file.units(r)
        if u:
            plt.ylabel(u)
        myFmt = DateFormatter("%Y-%b")
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_major_locator(plt.MaxNLocator(20))
        plt.grid()
        fig.autofmt_xdate()
        fig.savefig(p, format='pdf')
        plt.close(fig)

    p.close()

    if os.path.exists('max_devs.pkl'):
        with open('max_devs.pkl', 'rb') as f:
            max_deviations = pickle.load(f)
    else:
        max_deviations = {}

    max_dev = max(deviation, key = lambda x: deviation[x])
    max_deviations[test_name] = deviation[max_dev]

    with open('max_devs.pkl', 'wb') as f:
        pickle.dump(max_deviations, f)

# Rename files to rank them according to maximum deviations
def reorder_files():
    with open('max_devs.pkl', 'rb') as f:
        max_deviations = pickle.load(f)

    c = 1
    for file in sorted(max_deviations, key = lambda x: max_deviations[x], reverse=True):
        os.rename(f'{file}.pdf', f'{c:02d}_{file}.pdf')
        c += 1

# Main code
parser = argparse.ArgumentParser('plot_well_comparison.py')

parser.add_argument('-c', help='Name of test to process', dest='test_name')
parser.add_argument('-r', help='Reference file', dest='ref_file')
parser.add_argument('-s', help='Simulation file', dest='sim_file')
parser.add_argument('-t', help='Name for first simulation', dest='ref_name', default='Reference')
parser.add_argument('-u', help='Name for second simulation', dest='sim_name', default='New simulation')
parser.add_argument('-o', choices=['plot', 'rename'], help='Operation to do', required=True, dest='operation')
args = parser.parse_args()

if args.operation == 'plot':
    run_analysis(args.ref_file, args.sim_file, args.test_name, args.ref_name, args.sim_name)
else:
    reorder_files()
