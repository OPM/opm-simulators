#!/usr/bin/env python
import subprocess
import sys
import os.path

from ert.ecl import EclFile
from ert.test import TestAreaContext


def compare_files( flow_file , ref_file , kw_list):
    abs_epsilon = 0
    rel_epsilon = 1e-2
    flow = EclFile( flow_file )
    ref = EclFile( ref_file )

    for kw in kw_list:
        flow_kw = flow[kw][0]
        ref_kw = ref[kw][0]

        if not flow_kw.equal_numeric( ref_kw , abs_epsilon = abs_epsilon , rel_epsilon = rel_epsilon):
            first_different = ref_kw.firstDifferent( flow_kw , abs_epsilon = abs_epsilon , rel_epsilon = rel_epsilon)
            sys.exit("Keyword:%s was different in flow simulation and reference. Index of first difference: %d" % (kw , first_different))

            

        




flow = sys.argv[1]
data_file = sys.argv[2]
ref_init_file = sys.argv[3]
kw_list = sys.argv[4:]

with TestAreaContext("flow_init") as ta:
    subprocess.check_call( [flow , "nosim=true" , data_file ] )

    compare_files( os.path.splitext( os.path.basename( data_file ) )[0] + ".INIT" , ref_init_file , kw_list)


