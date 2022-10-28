#!/usr/bin/python3

from pathlib import Path
import json
import sys

p = Path(sys.argv[1])
for file in p.iterdir():
    with open(file, 'r') as f:
        data = f.read()
    js = json.loads(data)
    for test in js:
        if 'dir' in test:
            dir_name = test['dir']
        else:
            dir_name = test['casename']
        print(test['casename'])
        print('{} {} {}'.format(test['simulator'],dir_name,test['casename']))
