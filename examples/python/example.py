#!/usr/bin/env python
import sys
from trans_graph import TransGraph

# This file just contains some example functions for how one can poke
# around in the TransGraph datastructure.


def direction_count(tg):
    dir_count = {"X" : 0 , 
                 "Y" : 0 ,
                 "Z" : 0 ,
                 "NNC" : 0 }
    
    for cell in tg:
        if cell:
            for conn in cell:
                dir_count[ conn.dir ] += 1
    
    dir_count["Total"] = dir_count["X"] + dir_count["Y"] + dir_count["Z"] + dir_count["NNC"]
    return dir_count
                



def print_cell( prefix , cell ):
    print "%s: Cell: (%d,%d,%d) " % (prefix , cell.i , cell.j , cell.k)
    for conn in cell:
        print "    Connection => (%3d,%3d,%3d)   Transmissibility: %g     Direction: %s" % (conn.i , conn.j , conn.k , conn.T , conn.dir)
    
    
    print "    cell[\"X\"] => %s" % cell["X"]
    print "    cell[\"Y\"] => %s" % cell["Y"]
    print "    cell[\"Z\"] => %s" % cell["Z"]
    

def connection_to_count(tg , i,j,k):
    count = 0
    for cell in tg:
        if cell.connectsWith(i,j,k):
            count += 1
    
    return count
        











#-----------------------------------------------------------------

def direction_example( opm_tg , ecl_tg ):
    opm_count = direction_count( opm_tg )
    ecl_count = direction_count( ecl_tg )

    print "OPM: %s" % opm_count
    print "ECL: %s" % ecl_count
    print



def cell_example(opm_tg , ecl_tg ):
    opm_cell = opm_tg[21,27,10]
    ecl_cell = ecl_tg[21,27,10]

    print_cell( "OPM: " , opm_cell )
    print_cell( "ECL: " , ecl_cell )

        


def count_example(opm_tg , ecl_tg ):
    i = 5
    j = 10
    k = 7
    print "Opm connections to (%d,%d,%d): %d" % (i,j,k , connection_to_count(opm_tg , i,j,k))
    print "Ecl connections to (%d,%d,%d): %d" % (i,j,k , connection_to_count(ecl_tg , i,j,k))
    print 


def xtrace_example( opm_tg , ecl_tg , j , k):
    opm_trace = opm_tg.getXTrace( j , k)
    ecl_trace = ecl_tg.getXTrace( j , k)

    print "OPM: %s" % opm_trace
    print "ECL: %s" % ecl_trace
    print


def ytrace_example( opm_tg , ecl_tg , i , k):
    opm_trace = opm_tg.getYTrace(  i , k )
    ecl_trace = ecl_tg.getYTrace(  i , k )

    print "OPM: %s" % opm_trace
    print "ECL: %s" % ecl_trace
    print


def ztrace_example( opm_tg , ecl_tg , i , j):
    opm_trace = opm_tg.getZTrace( i , j )
    ecl_trace = ecl_tg.getZTrace( i , j )

    print "OPM: %s" % opm_trace
    print "ECL: %s" % ecl_trace
    print
    

#-----------------------------------------------------------------

if len(sys.argv) < 3:
    sys.exit("example.py  opm_trans.json  eclipse_trans.json")

opm_tg = TransGraph.load(sys.argv[1])
ecl_tg = TransGraph.load(sys.argv[2])


direction_example( opm_tg , ecl_tg )
cell_example( opm_tg , ecl_tg )
count_example( opm_tg , ecl_tg )


xtrace_example( opm_tg , ecl_tg , 20 ,20)
ytrace_example( opm_tg , ecl_tg , 10 ,10)

ztrace_example( opm_tg , ecl_tg , 10 ,10)
ztrace_example( opm_tg , ecl_tg , 10 ,20)
ztrace_example( opm_tg , ecl_tg , 20 ,10)
ztrace_example( opm_tg , ecl_tg , 20 ,20)
ztrace_example( opm_tg , ecl_tg , 30 ,70)
ztrace_example( opm_tg , ecl_tg , 30 ,60)
