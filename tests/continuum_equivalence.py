#!/usr/bin/env python3
"""Compare one continuum of a dual-continuum run with a single porosity run.

With no flow between the continua, the cells of one continuum behave as a single porosity
grid made of those cells. This compares the INIT arrays, the restart arrays and the summary
vectors of the wells of such a run (the dual-continuum case) with the single porosity
reference. Only Python's standard library is used.

usage: continuum_equivalence.py <single porosity dir> <case> <dual dir> <case> fracture|matrix
"""
import struct
import sys

ABS_TOL = 1.0e-4
REL_TOL = 1.0e-5
INIT_ABS_TOL = 1.0e-5
INIT_REL_TOL = 1.0e-6

SUMMARY_VECTORS = {'FOPR', 'FWPR', 'FGPR', 'FOPT', 'FWPT', 'FGPT', 'FWIR', 'FGIR', 'FWIT', 'FGIT',
                   'WOPR', 'WWPR', 'WGPR', 'WWIR', 'WGIR', 'WBHP', 'WOPT', 'WWPT', 'WGPT'}


def read_blocks(path):
    """The blocks (name, values) of an unformatted Eclipse file."""
    out = []
    with open(path, 'rb') as f:
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            name = f.read(8).decode().strip()
            (n,) = struct.unpack('>i', f.read(4))
            typ = f.read(4).decode()
            f.read(4)
            elem = {'INTE': ('i', 4), 'REAL': ('f', 4), 'DOUB': ('d', 8), 'LOGI': ('i', 4),
                    'CHAR': ('s', 8), 'MESS': ('s', 0)}[typ]
            data = []
            remaining = n
            while remaining > 0:
                (bsize,) = struct.unpack('>i', f.read(4))
                count = bsize // elem[1] if elem[1] else 0
                raw = f.read(bsize)
                f.read(4)
                if elem[0] == 's':
                    data += [raw[i * 8:(i + 1) * 8].decode().strip() for i in range(count)]
                else:
                    data += list(struct.unpack('>' + str(count) + elem[0], raw))
                remaining -= count
            out.append((name, data))
    return out


def differ(reference, value, abs_tol, rel_tol):
    """As compareECL: a value fails only if both tolerances are exceeded."""
    d = abs(reference - value)
    return d > abs_tol and d > rel_tol * abs(reference)


def main():
    sp_dir, sp_case, du_dir, du_case, which = sys.argv[1:6]
    failures = []

    # ------------------------------------------------------------------ INIT
    sp_init = read_blocks('%s/%s.INIT' % (sp_dir, sp_case))
    du_init = dict(read_blocks('%s/%s.INIT' % (du_dir, du_case)))
    n_sp = dict(sp_init)['INTEHEAD'][11]       # number of active cells
    n_du = du_init['INTEHEAD'][11]

    def part(values):
        return values[:n_sp] if which == 'matrix' else values[len(values) - n_sp:]

    compared = 0
    for name, data in sp_init:
        if len(data) != n_sp or name not in du_init or len(du_init[name]) != n_du or name == 'ACTNUM':
            continue
        sliced = part(du_init[name])
        compared += 1
        bad = [i for i, (a, b) in enumerate(zip(data, sliced)) if differ(a, b, INIT_ABS_TOL, INIT_REL_TOL)]
        if bad:
            failures.append('INIT %s differs in %d cells, first %d: %g vs %g'
                            % (name, len(bad), bad[0], data[bad[0]], sliced[bad[0]]))
    print('INIT: %d arrays of %d active cells compared with the %s continuum of %d cells'
          % (compared, n_sp, which, n_du))
    if compared < 10:
        failures.append('too few INIT arrays compared: %d' % compared)

    # ----------------------------------------------------------------- UNRST
    def steps(blocks):
        out, cur = [], {}
        for name, data in blocks:
            if name == 'SEQNUM':
                if cur:
                    out.append(cur)
                cur = {}
            elif name in ('PRESSURE', 'SWAT', 'SGAS', 'RS', 'RV'):
                cur[name] = data
        if cur:
            out.append(cur)
        return out

    s_sp = steps(read_blocks('%s/%s.UNRST' % (sp_dir, sp_case)))
    s_du = steps(read_blocks('%s/%s.UNRST' % (du_dir, du_case)))
    if len(s_sp) != len(s_du):
        failures.append('report steps: %d against %d' % (len(s_sp), len(s_du)))
    for key in ('PRESSURE', 'SWAT', 'SGAS', 'RS', 'RV'):
        worst = 0.0
        for a, b in zip(s_sp, s_du):
            if key not in a or key not in b:
                continue
            for x, y in zip(a[key], part(b[key])):
                if differ(x, y, ABS_TOL, REL_TOL):
                    failures.append('UNRST %s: %g vs %g' % (key, x, y))
                    break
            worst = max([worst] + [abs(x - y) for x, y in zip(a[key], part(b[key]))])
        print('UNRST %-8s largest difference %.2e' % (key, worst))

    # --------------------------------------------------------------- SUMMARY
    def summary(case_dir, case):
        spec = dict(read_blocks('%s/%s.SMSPEC' % (case_dir, case)))
        keys = list(zip(spec['KEYWORDS'], spec['WGNAMES']))
        rows = [d for n, d in read_blocks('%s/%s.UNSMRY' % (case_dir, case)) if n == 'PARAMS']
        return {k: i for i, k in enumerate(keys)}, rows

    i_sp, r_sp = summary(sp_dir, sp_case)
    i_du, r_du = summary(du_dir, du_case)
    time_key = ('TIME', ':+:+:+:+')
    t_sp = {round(row[i_sp[time_key]], 6): row for row in r_sp}
    t_du = {round(row[i_du[time_key]], 6): row for row in r_du}
    common = sorted(set(t_sp) & set(t_du))
    if len(common) != len(t_sp) or len(common) != len(t_du):
        failures.append('summary times: %d and %d, %d in common' % (len(t_sp), len(t_du), len(common)))
    compared = 0
    for key, i in i_sp.items():
        if key[0] not in SUMMARY_VECTORS or key not in i_du:
            continue
        compared += 1
        for t in common:
            if differ(t_sp[t][i], t_du[t][i_du[key]], ABS_TOL, REL_TOL):
                failures.append('SMRY %s:%s at %g: %g vs %g' % (key[0], key[1], t, t_sp[t][i], t_du[t][i_du[key]]))
                break
    print('SUMMARY: %d vectors at %d times' % (compared, len(common)))
    if compared < 10:
        failures.append('too few summary vectors compared: %d' % compared)

    for line in failures:
        print('FAIL', line)
    return 1 if failures else 0


sys.exit(main())
