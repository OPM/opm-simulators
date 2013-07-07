#! /usr/bin/python
import re
import sys

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
            if curFieldName.find("saturation") >= 0:
                if abs(number1 - number2) > 0.1:
                    print 'Difference between %f and %f too large in data field "%s": %s'%(number1,number2,curFieldName,abs(number1 - number2))
                    return False
            if curFieldName.find("mole") >= 0 or curFieldName.find("mass") >= 0:
                if abs(number1 - number2) > 0.1:
                    print 'Difference between %f and %f too large in data field "%s": %s'%(number1,number2,curFieldName,abs(number1 - number2))
                    return False
            elif curFieldName.find("velocity") >= 0:
                if abs(number1 - number2) > 0.02:
                    print 'Difference between %f and %f too large in data field "%s": %s'%(number1,number2,curFieldName,abs(number1 - number2))
                    return False
            elif curFieldName.find("pressure") >= 0:
                if abs(number1 - number2) > 0.1*abs(number1 + number2):
                    print 'Difference between %f and %f too large in data field "%s": %s'%(number1,number2,curFieldName,abs(number1 - number2))
                    return False
            else:
                # don't compare any other fields
                pass
    return True

if len(sys.argv) != 3:
    print sys.argv[0], "CURRENT_RESULT REFERENCE_RESULT"
    exit(1);

if isFuzzyEqual(open(sys.argv[1]), open(sys.argv[2]), absTol=1e-6, relTol=1e-2):
    exit(0)
else:
    exit(1)
