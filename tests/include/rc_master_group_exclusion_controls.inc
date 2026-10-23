-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- The slaves and the group controls shared by the decks
-- RC_MASTER_GROUP_EXCLUSION*.DATA for the unit tests in
-- tests/rescoup/test_rescoup_constraints.cpp.  The decks differ only in the
-- group hierarchy (GRUPTREE), which must come before this file.

SLAVES
-- name  DATA file  host  directory  processes
 'S1'   'S1'       1*    'slave1'   1 /
 'S2'   'S2'       1*    'slave2'   1 /
/

GRUPMAST
-- master group  slave  slave group
 'M1'  'S1'  'G1' /
 'M2'  'S1'  'G2' /
 'M3'  'S2'  'G3' /
 'M4'  'S2'  'G4' /
/

-- M1, M2 and M3 are on their own (non-binding) ORAT limits, but available for
-- higher-level control, so they take part in FIELD's guide-rate distribution.
-- FIELD also has an LRAT limit of 2000, so that a test can check the limits
-- the master sends besides the ORAT target.
GCONPROD
-- group  cmode   oil     water  gas  liquid  action  item8  guide rate  phase
 'FIELD' 'ORAT'  1000.0  1*     1*   2000.0  'RATE' /
 'M1'    'ORAT'  5000.0  1*     1*   1*      'RATE'  'YES'  100.0       'OIL' /
 'M2'    'ORAT'  5000.0  1*     1*   1*      'RATE'  'YES'  200.0       'OIL' /
 'M3'    'ORAT'  5000.0  1*     1*   1*      'RATE'  'YES'  300.0       'OIL' /
 'M4'    'ORAT'  50.0    1*     1*   1*      'RATE'  'NO'   1*          1* /
/

GCONINJE
-- group  phase    cmode   surface  resv  reinj  voidage  item8  guide rate  def
 'FIELD' 'WATER'  'RATE'  1000.0 /
 'FIELD' 'GAS'    'RATE'  1000.0 /
 'M1'    'WATER'  'RATE'  5000.0   1*    1*     1*       'YES'  100.0       'RATE' /
 'M2'    'WATER'  'RATE'  5000.0   1*    1*     1*       'YES'  200.0       'RATE' /
 'M3'    'WATER'  'RATE'  5000.0   1*    1*     1*       'YES'  300.0       'RATE' /
 'M4'    'WATER'  'RATE'  50.0     1*    1*     1*       'NO' /
 'M1'    'GAS'    'RATE'  5000.0   1*    1*     1*       'YES'  100.0       'RATE' /
 'M2'    'GAS'    'RATE'  5000.0   1*    1*     1*       'YES'  200.0       'RATE' /
 'M3'    'GAS'    'RATE'  5000.0   1*    1*     1*       'YES'  300.0       'RATE' /
 'M4'    'GAS'    'RATE'  50.0     1*    1*     1*       'NO' /
/
