#!/bin/bash
#
# run-parallel-init-regressionTest.sh
# ---------------------------------------------------------------------------
# Synthetic serial-vs-parallel regression for the parallel LGR .INIT / .EGRID output.
#
# The LGR deck (opm-common's LGR-WELL-3X3-1LGR: 3x3x1, one CARFIN box) is embedded below as a
# heredoc -- a self-contained synthetic case, no external .DATA dependency. It is written to a
# temp file and run SERIAL and PARALLEL with --enable-dry-run=true (NOSIM: the INIT+EGRID are
# written at init and the simulation is skipped, so this does not depend on the per-step gather).
# The two runs' .EGRID and .INIT are compared with opm-common's compareECL; -x on the INIT
# ignores the parallel-only MPI_RANK keyword.
#
# Pre-fix (the coarse per-rank transmissibility is handed to the LGR INIT writer) the parallel
# run deadlocks in finishInit -> this test times out. With the refined-transmissibility fix the
# parallel run completes and its INIT/EGRID match the serial reference.
#
if test $# -eq 0
then
  echo -e "Usage:\t$0 -r <result> -a <abstol> -t <reltol> -c <compareECL> -e <exe> [-n <procs>]"
  exit 1
fi

MPI_PROCS=4
OPTIND=1
while getopts "i:r:f:a:t:c:e:n:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    n) MPI_PROCS=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

CASE=LGR_WELL_3X3_1LGR
rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}/mpi
DECK_FILE=${RESULT_PATH}/${CASE}.DATA

# ---- synthetic LGR deck (embedded; written to a temp file) ----
cat > ${DECK_FILE} <<'OPM_SYNTHETIC_DECK_EOF'
-- This reservoir simulation deck is made available under the Open Database
-- License: http://opendatacommons.org/licenses/odbl/1.0/. Any rights in
-- individual contents of the database are licensed under the Database Contents
-- License: http://opendatacommons.org/licenses/dbcl/1.0/

-- Copyright (C) 2023 Equinor

-- Minimal dummy case based on SPE1 (Odeh 1981) physics.
-- Purpose: simple single-LGR development case for PR00 (Fixes 1, 2, 3).
--
-- Global grid: 3 x 3 x 1
--   +-------+-------+-------+
--   | (1,3) | (2,3) | (3,3) | <- PROD here (global grid, no LGR)
--   +-------+-------+-------+
--   | (1,2) | (2,2) | (3,2) |
--   +-------+-------+-------+
--   |(1,1)  | (2,1) | (3,1) | <- LGR1 host (INJ)
--   +-------+-------+-------+
--
-- LGR1: global cell (1,1,1) refined into 3x3x1 fine cells
--   Fine cell layout:
--     (1,3)(2,3)(3,3)
--     (1,2)(2,2)(3,2)  <- INJ at centre (2,2,1)
--     (1,1)(2,1)(3,1)
--
-- Wells:
--   INJ  - gas injector - LGR1 fine cell (2,2,1)  — uses WELSPECL/COMPDATL
--   PROD - oil producer - global cell   (3,3,1)   — uses WELSPECS/COMPDAT
--
-- PR00 fix coverage:
--   Fix 1: all three second-letter branches (LW* -> Well, LC* -> Connection,
--           LB* -> Block) are exercised via INJ's LGR summary vectors.
--   Fix 2: lgr_ field is populated and converted for INJ; PROD has lgr_={}
--           confirming non-LGR wells are unaffected (regression).
--   Fix 3: serializeOp / operator== on lgr_info — covered by unit tests.
--   Fix 4: NOT the focus here (only one LGR, no dedup scenario possible).
--           Use LGR-WELL-3X3-2LGR.DATA to stress-test Fix 4.

---------------------------------------------------------------------------
------- 3x3x1 1LGR (INJ in LGR1, PROD in global) — development case ------
---------------------------------------------------------------------------

RUNSPEC
-- -----------------------------------------------------------------------

TITLE
   3x3x1 1LGR - LGR1(1,1,1) INJ, PROD in global (3,3,1)

DIMENS
   3 3 1 /

EQLDIMS
/

TABDIMS
/

OIL
GAS
WATER
DISGAS

FIELD

START
   1 'JAN' 2015 /

WELLDIMS
-- MaxWells  MaxCompl  MaxGroups  MaxWellsInGroup
        2         1         1          2 /

UNIFOUT

GRID
-- -----------------------------------------------------------------------

-- LGR1: refine global cell (1,1,1) into 3x3x1 fine cells
CARFIN
-- NAME  I1 I2  J1 J2  K1 K2  NX  NY  NZ
'LGR1'   1  1   1  1   1  1   3   3   1 /
ENDFIN

INIT

DX
   9*1000 /

DY
   9*1000 /

DZ
   9*50 /

TOPS
   9*8325 /

PORO
   9*0.3 /

PERMX
   9*500 /

PERMY
   9*200 /

PERMZ
   9*200 /

ECHO

PROPS
-- -----------------------------------------------------------------------

PVTW
-- PRef(psia)  Bw(rb/stb)  Cw(psi-1)   Visc(cP)  Viscosibility
    4017.55     1.038       3.22E-6      0.318      0.0 /

ROCK
-- PRef(psia)  Cr(psi-1)
    14.7        3E-6 /

SWOF
0.12   0.000000000   1.000   0
0.18   4.64876E-008  1.000   0
0.24   1.86000E-007  0.997   0
0.30   4.18388E-007  0.980   0
0.36   7.43802E-007  0.700   0
0.42   1.16219E-006  0.350   0
0.48   1.67355E-006  0.200   0
0.54   2.27789E-006  0.090   0
0.60   2.97521E-006  0.021   0
0.66   3.76550E-006  0.010   0
0.72   4.64876E-006  0.001   0
0.78   5.62500E-006  0.0001  0
0.84   6.69421E-006  0.000   0
0.91   8.05914E-006  0.000   0
1.00   0.000010000   0.000   0 /

SGOF
0.000  0.000  1.000  0
0.001  0.000  1.000  0
0.020  0.000  0.997  0
0.050  0.005  0.980  0
0.120  0.025  0.700  0
0.200  0.075  0.350  0
0.250  0.125  0.200  0
0.300  0.190  0.090  0
0.400  0.410  0.021  0
0.450  0.600  0.010  0
0.500  0.720  0.001  0
0.600  0.870  0.000  0
0.700  0.940  0.000  0
0.850  0.980  0.000  0
0.880  0.984  0.000  0 /

DENSITY
-- Oil(lb/ft3)  Water(lb/ft3)  Gas(lb/ft3)
   53.66        64.49          0.0533 /

PVDG
--  P(psia)   Bg(rb/Mscf)   Visc(cP)
     14.700   166.666        0.00800
    264.700    12.093        0.00960
    514.700     6.274        0.01120
   1014.700     3.197        0.01400
   2014.700     1.614        0.01890
   2514.700     1.294        0.02080
   3014.700     1.080        0.02280
   4014.700     0.811        0.02680
   5014.700     0.649        0.03090
   9014.700     0.386        0.04700 /

PVTO
--  Rs(Mscf/stb)  Pbub(psia)  Bo(rb/stb)  Visc(cP)
    0.0010         14.7        1.0620       1.0400 /
    0.0905        264.7        1.1500       0.9750 /
    0.1800        514.7        1.2070       0.9100 /
    0.3710       1014.7        1.2950       0.8300 /
    0.6360       2014.7        1.4350       0.6950 /
    0.7750       2514.7        1.5000       0.6410 /
    0.9300       3014.7        1.5650       0.5940 /
    1.2700       4014.7        1.6950       0.5100
               9014.7        1.5790       0.7400 /
    1.6180       5014.7        1.8270       0.4490
               9014.7        1.7370       0.6310 /
/

SOLUTION
-- -----------------------------------------------------------------------

RPTSOL
RESTART=1 SOIL SWAT SGAS /

EQUIL
-- Datum(ft)  P@datum(psia)  WOC(ft)  Pcow@WOC  GOC(ft)  Pcog@GOC  RSVD  RVVD  NRPVT
      8400       4800          8450      0         8300      0        1     0     0 /

RSVD
-- Depth(ft)  Rs(Mscf/stb)
   8300        1.270
   8450        1.270 /

SUMMARY
-- -----------------------------------------------------------------------

-- =========================================================
-- FIELD LEVEL
-- =========================================================
FOPR
FGOR
FWCT
FWIR

-- =========================================================
-- GLOBAL WELL LEVEL  (W*)
-- Both wells at global level. Cross-check: W* must equal LW* for INJ.
-- PROD has no LW* — confirms non-LGR wells are unaffected (Fix 2 regression).
-- =========================================================
WBHP
  'INJ'
  'PROD'
/
WOPR
  'INJ'
  'PROD'
/
WGIR
  'INJ'
  'PROD'
/
WGOR
  'PROD'
/
WWPR
  'INJ'
  'PROD'
/

-- =========================================================
-- GLOBAL BLOCK LEVEL  (B*)
-- Cell (1,1,1) is the LGR1 host — coarse aggregate.
-- Cell (3,3,1) holds PROD. Cell (2,2,1) is the undisturbed reference.
-- Cross-check: BPR(1,1,1) ≈ LBPR(LGR1,2,2,1).
-- =========================================================
BPR
  1  1  1 /
  2  2  1 /
  3  3  1 /
/

BGSAT
  1  1  1 /
  2  2  1 /
  3  3  1 /
/

-- =========================================================
-- LGR WELL LEVEL  (LW*)
-- Syntax: LWkeyword / 'LGRName'  'WellName' / /
-- Fix 1: exercises the 'W' branch of case 'L' in category_from_keyword().
-- Fix 2: lgr_ = {LGR1, {0,0,0}} populated and converted to EclIO::SummaryNode.
-- =========================================================

LWBHP
  'LGR1'  'INJ' /
/

LWGIR
  'LGR1'  'INJ' /
/

LWGIT
  'LGR1'  'INJ' /
/

LWVIR
  'LGR1'  'INJ' /
/

LWVIT
  'LGR1'  'INJ' /
/

-- Regression sentinel: must be zero for a pure gas injector
LWWIR
  'LGR1'  'INJ' /
/

LWSTAT
  'LGR1'  'INJ' /
/

-- =========================================================
-- LGR CONNECTION LEVEL  (LC*)
-- Syntax: LCkeyword / 'LGRName'  'WellName'  I  J  K / /
-- INJ perforated at LGR1 fine cell (2,2,1).
-- Fix 1: exercises the 'C' branch of case 'L' in category_from_keyword().
-- =========================================================

LCGFR
  'LGR1'  'INJ'  2  2  1 /
/

-- Regression sentinel: must be zero for a pure gas injector
LCOFR
  'LGR1'  'INJ'  2  2  1 /
/

LCWFR
  'LGR1'  'INJ'  2  2  1 /
/

LCPR
  'LGR1'  'INJ'  2  2  1 /
/

-- =========================================================
-- LGR BLOCK LEVEL  (LB*)
-- Syntax: LBkeyword / 'LGRName'  I  J  K / /
-- Three cells: corners (1,1,1), (3,3,1) and centre (2,2,1).
-- Fix 1: exercises the 'B' branch of case 'L' in category_from_keyword().
-- =========================================================

-- Block pressure gradient across LGR1
LBPR
  'LGR1'  1  1  1 /
  'LGR1'  2  2  1 /
  'LGR1'  3  3  1 /
/

-- Gas saturation grows from injector cell (2,2,1) outward
LBGSAT
  'LGR1'  1  1  1 /
  'LGR1'  2  2  1 /
  'LGR1'  3  3  1 /
/

-- Saturation accounting at injector cell (Sw + So + Sg = 1)
LBSOIL
  'LGR1'  2  2  1 /
/

LBWSAT
  'LGR1'  2  2  1 /
/

SCHEDULE
-- -----------------------------------------------------------------------

RPTSCHED
   'PRES' 'SGAS' 'RS' 'WELLS' /

RPTRST
   'BASIC=1' /

DRSDT
 0 /

-- INJ in LGR1 -> WELSPECL + COMPDATL
-- PROD in global grid -> WELSPECS + COMPDAT
WELSPECL
-- WellName  Group  LGRName   I  J  RefDepth(ft)  Phase
   'INJ'    'G1'   'LGR1'    2  2   8335         'GAS' /
/

COMPDATL
-- WellName  LGRName   I  J  K1 K2  Flag    Tmult  PI     WellRad(ft)
   'INJ'    'LGR1'     2  2  1  1  'OPEN'   1*     1*     0.5 /
/

WELSPECS
-- WellName  Group  I  J  RefDepth(ft)  Phase
   'PROD'   'G1'   3  3   8400        'OIL' /
/

COMPDAT
-- WellName  I  J  K1 K2  Flag    Tmult  PI     WellRad(ft)
   'PROD'   3  3  1  1  'OPEN'   1*     1*     0.5 /
/

WCONPROD
-- WellName  Flag   Control  OilRate(stb/d)  ...  BHPmin(psia)
   'PROD'   'OPEN' 'ORAT'    20000           4*    1000 /
/

WCONINJE
-- WellName  Phase  Flag    Control  Rate(Mscf/d)  1*  BHPmax(psia)
   'INJ'    'GAS' 'OPEN'  'RATE'    100000         1*   9014 /
/

TSTEP
-- 12 monthly steps (Jan-Dec 2015)
   31 28 31 30 31 30 31 31 30 31 30 31 /

END
OPM_SYNTHETIC_DECK_EOF
# ---- end synthetic deck ----

COMMON_ARGS="${TEST_ARGS} --parsing-strictness=low --enable-dry-run=true --enable-ecl-output=true"

echo "=== serial run (NOSIM: writes INIT + EGRID) ==="
"${EXE_NAME}" ${COMMON_ARGS} --output-dir=${RESULT_PATH} ${DECK_FILE}
test $? -eq 0 || { echo "FAIL: serial run errored"; exit 1; }

echo "=== parallel run (NOSIM, ${MPI_PROCS} ranks) ==="
mpirun -np ${MPI_PROCS} "${EXE_NAME}" ${COMMON_ARGS} --output-dir=${RESULT_PATH}/mpi ${DECK_FILE}
test $? -eq 0 || { echo "FAIL: parallel run errored (deadlock/crash?)"; exit 1; }

S=${RESULT_PATH}/${CASE}
P=${RESULT_PATH}/mpi/${CASE}
ecode=0

echo "=== compareECL EGRID (serial vs parallel) ==="
${COMPARE_ECL_COMMAND} -t EGRID ${S} ${P} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]; then ecode=1; ${COMPARE_ECL_COMMAND} -a -t EGRID ${S} ${P} ${ABS_TOL} ${REL_TOL}; fi

echo "=== compareECL INIT (serial vs parallel; -x ignores parallel-only MPI_RANK) ==="
${COMPARE_ECL_COMMAND} -t INIT -x ${S} ${P} ${ABS_TOL} ${REL_TOL}
if [ $? -ne 0 ]; then ecode=1; ${COMPARE_ECL_COMMAND} -a -t INIT -x ${S} ${P} ${ABS_TOL} ${REL_TOL}; fi

exit $ecode
