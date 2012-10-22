#! /usr/bin/python

import argparse
import re

# fuzzy compare two VTK files
def isFuzzyEqual(vtkFile1, vtkFile2, absTol, relTol):
    inField = False
    curFieldName = ""
    for curLine1 in vtkFile1.xreadlines():
        curLine2 = vtkFile2.readline()

        if curLine1.find("</DataArray>") >= 0:
            inField = False
            continue

        m = re.search(r'Name="([a-zA-Z0-9 _\-]*)"', curLine1)
        if m:
            curFieldName = m.group(1)
            inField = True
            continue

        if not inField:
            continue

        curVals1 = map(float, curLine1.split())
        curVals2 = map(float, curLine2.split())
        
        if len(curVals1) != len(curVals2):
            print "Length of field '%s' is different"%curFieldName
            return False
        
        for i in range(0, len(curVals1)):
            number1 = curVals1[i]
            number2 = curVals2[i]
            if abs(number1 - number2) > absTol and number2 != 0 and abs(number1/number2 - 1) > relTol:
                print 'Difference between %f and %f too large (%f%%) in data field "%s"'%(number1,number2,abs(number1/number2 - 1)*100, curFieldName)
                return False
    return True

# main programm
# handle arguments and print help message
parser = argparse.ArgumentParser(description='Fuzzy compare two VTK\
    (Visualization Toolkit) files. The files are accepted if for every\
    value the difference is below the absTol error or below the\
    relTol error or below both.')
parser.add_argument('vtu_file_1', type=open,
    help='first file to compare')
parser.add_argument('vtu_file_2', type=open,
    help='second file to compare')
parser.add_argument('-r', '--relTol', type=float, default=1e-2,
    help='maximum tolerated absolute error (default=1e-2)')
parser.add_argument('-a', '--absTol', type=float, default=1e-6,
    help='maximum tolerated relative error (default=1e-6)')
args = parser.parse_args()

# fuzzy compare
if isFuzzyEqual(args.vtu_file_1, args.vtu_file_2, args.absTol, args.relTol):
    exit(0)
else:
    exit(1)
