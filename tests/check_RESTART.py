#!/usr/bin/env python
import subprocess
import sys
import os.path

from ert.ecl import EclFile
from ert.test import TestAreaContext


def compare_files( flow_file , ref_file):
    flow = EclFile( flow_file )
    ref = EclFile( ref_file )

    for kw in ["PRESSURE" , "SWAT" , "SGAS" , "RS" , "RV"]:
        flow_kw = flow[kw][-1]
        ref_kw = ref[kw][-1]

        if not flow_kw.equal_numeric( ref_kw , abs_epsilon = 0, rel_epsilon = 2.5e-2 ):
            first_different = flow_kw.firstDifferent( ref_kw , abs_epsilon = 0, rel_epsilon = 2.5e-2 )
            sys.exit("Keyword:%s was different in flow simulation and reference. First difference in index:%d" % (kw , first_different))


#-----------------------------------------------------------------
# Small script running flow twice, once with restart and normal run
# from step 0. The final state of the two runs is compared.

flow = sys.argv[1]
full_deck = sys.argv[2]
restart_deck = sys.argv[3]

with TestAreaContext("flow_init") as ta:
    ta.copy_file( full_deck )
    ta.copy_file( restart_deck )

    subprocess.check_call( [flow , os.path.basename( full_deck ) ] )
    subprocess.check_call( [flow , os.path.basename( restart_deck ) ] )

    compare_files( os.path.splitext( os.path.basename( full_deck ) )[0] + ".UNRST" , os.path.splitext( os.path.basename( restart_deck ) )[0] + ".UNRST")


