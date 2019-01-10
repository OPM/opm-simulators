#!/usr/bin/env python3
import sys
import os
import argparse
import numpy
def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--outputdir')
    parser.add_argument('--simulator')
    parser.add_argument('--datadir')
    parser.add_argument('--case');
    parser.add_argument('--max_diff',default=1e-6);
    args=parser.parse_args()
    print args
    
    options = "--enable-adaptive-time-stepping=false \
    --solve-welleq-initially=true --use-adjoint=true \
    --linear-solver-reduction=1e-10 \
    --tolerance-cnv=1e-6 \
    --tolerance-cnv-relaxed=1e-5\
    --tolerance-mb=1e-6 \
    --tolerance-pressure-ms-wells=1000 \
    --tolerance-well-control=1e-7 \
    --tolerance-wells=1e-4"
    
    datadir = "%s/%s" % (args.datadir,args.case)
    deckfile = "%s/inputfiles/%s.DATA" % (datadir,args.case)
    adjoint_test_file = "%s/adjoint_results.txt" % (args.outputdir)
    if os.path.exists(adjoint_test_file):
        os.remove(adjoint_test_file)
    print "************************************************************"
    scommand= "%s --output-dir=%s %s %s  " % (args.simulator, args.outputdir, options, deckfile)
    print scommand 
    print "************************************************************"
    ok = os.system(scommand)
    adjoint_ref_file = "%s/adjoint_results.txt" % (datadir)
    
    adj_ref = numpy.genfromtxt(adjoint_ref_file,comments='%')
    adj_test = numpy.genfromtxt(adjoint_test_file,comments='%')
    diff=numpy.max(numpy.abs(adj_ref-adj_test))
    print "max difference %e" % diff
    ok = (diff < args.max_diff) and ok==0
    if ok:
        return 0
    else:
        return 1
       
if __name__ == "__main__":
    sys.exit(main())
else:
    main()
