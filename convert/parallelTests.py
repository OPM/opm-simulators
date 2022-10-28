#!/usr/bin/python3

import json
import os
import re

import shutil

try:
    shutil.rmtree("definitions/parallel")
    os.makedirs('definitions/parallel')
except:
    pass

with open('../parallelTests.cmake', 'r') as f:
    s = f.read()

for match in re.finditer('^add_test_compare_parallel_simulation\(([^\)]*)\)', s, re.MULTILINE):
    params = match.group(1).split()
    i = 0
    mapping = dict()
    mapping['test_cases'] = [dict()]
    while i < len(params):
        if params[i] == 'TEST_ARGS':
            j = 1
            val = ""
            while i + j < len(params) and not params[i+j][0].isupper():
                val = val + ' ' + params[i+j]
                j = j + 1
            mapping['test_cases'][0][params[i].lower()] = val.strip()
            i = i + j
        else:
            if params[i] == 'ABS_TOL' and params[i+1] == '${abs_tol}':
                params[i+1] = 2e-2
            if params[i] == 'REL_TOL' and params[i+1] == '${rel_tol}':
                params[i+1] = 1e-5
            if params[i] == 'REL_TOL' and params[i+1] == '${coarse_rel_tol}':
                params[i+1] = 1e-2
            if params[i] == 'ABS_TOL' and params[i+1] == '${abs_tol_parallel}':
                params[i+1] = 0.02
            if params[i] == 'REL_TOL' and params[i+1] == '${rel_tol_parallel}':
                params[i+1] = 8e-5
            if params[i] == 'REL_TOL' and params[i+1] == '${coarse_rel_tol_parallel}':
                params[i+1] = 1e-2
            if params[i] == 'CASENAME' or params[i] == 'FILENAME' or params[i] == 'RESTART_STEP' or params[i] == 'RESTART_SCHEDULE':
                mapping['test_cases'][0][params[i].lower()] = params[i+1]
            else:
                mapping[params[i].lower()] = params[i+1]
            i = i + 2
    with open('definitions/parallel/{}.json'.format(mapping['test_cases'][0]['casename']),'w') as f:
        f.write(json.dumps(mapping, indent=4))
