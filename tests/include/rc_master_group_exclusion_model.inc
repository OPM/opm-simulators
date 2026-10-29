-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- The one-cell model shared by the decks RC_MASTER_GROUP_EXCLUSION*.DATA for
-- the unit tests in tests/rescoup/test_rescoup_constraints.cpp: the dimensions,
-- the grid, the fluid properties and the initial state.  There are no wells: the
-- master groups get their rates from the data the test sets for the slaves.
-- WELLDIMS allows for FIELD, the four master groups and the plain group P of
-- the nested decks.

DIMENS
   1 1 1 /

OIL
WATER
GAS
DISGAS

METRIC

START
 1 'JAN' 2020 /

WELLDIMS
-- max wells, max cons/well, max groups, max wells/group
   1 1 7 5 /

TABDIMS
/

EQLDIMS
/

UNIFIN
UNIFOUT

-- =============================================================================
GRID
-- =============================================================================

DXV
100.0 /

DYV
100.0 /

DZV
10.0 /

TOPS
2000 /

PORO
0.25 /

PERMX
500 /

COPY
  PERMX PERMY /
  PERMX PERMZ /
/

-- =============================================================================
PROPS
-- =============================================================================

PVTW
-- ref.pres  Bw      Cw          Vw    Vw.visc
   250.0    1.02   4.0E-05      0.5    0.0 /

ROCK
-- ref.pres  compressibility
   250.0    4.0E-05 /

DENSITY
-- oil    water   gas
   850    1014    1.0 /

PVDO
-- P       Bo      Vo
  100.0   1.05    2.0
  250.0   1.02    2.5
  400.0   1.00    3.0 /

PVDG
-- P       Bg       Vg
  100.0   0.010    0.015
  250.0   0.005    0.020
  400.0   0.003    0.025 /

SWOF
-- Sw      Krw      Kro      Pcow
  0.2     0.0      1.0      0.0
  0.3     0.05     0.6      0.0
  0.5     0.15     0.3      0.0
  0.7     0.35     0.1      0.0
  0.8     0.5      0.0      0.0 /

SGOF
-- Sg      Krg      Kro      Pcog
  0.0     0.0      1.0      0.0
  0.1     0.02     0.6      0.0
  0.3     0.10     0.3      0.0
  0.5     0.25     0.1      0.0
  0.7     0.5      0.0      0.0 /

-- =============================================================================
SOLUTION
-- =============================================================================

EQUIL
-- datum  pressure  OWC    Pc@OWC  GOC    Pc@GOC  RSVD  RVVD  SOLN
   2000   250.0     2100   0.0     1900   0.0     1     0     0 /

RSVD
   1900   100.0
   2100   100.0 /
