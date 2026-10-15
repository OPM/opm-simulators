# Copyright 2026 Equinor
#
# This file is part of the Open Porous Media project (OPM).
#
# OPM is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# OPM is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with OPM.  If not, see <http://www.gnu.org/licenses/>.

# Reservoir coupling regression tests: run a coupled case, where the master
# spawns its slaves, and compare the master's and the slaves' results with
# reference results in opm-tests.
# Included from rescoupTests.cmake when opm-tests is found, which provides
# the MPI launcher in _rescoup_mpiexec.

set(abs_tol 2e-2)
set(rel_tol 1e-5)

# SLAVES and GRUPMAST are flagged as unsupported keywords, so the coupled
# decks need --parsing-strictness=low. The slaves receive the master's
# options, including these.
add_test_compare_rescoup_simulation(
  CASENAME
    rc_m1
  DIR
    rescoup/rc_m1
  MASTER_DIR
    master
  FILENAME
    RC-01_MAST_PRED
  SLAVE_FILENAMES
    RC-01_MOD1_PRED
    RC-01_MOD2_PRED
  # SLAVES item 5 of each slave above
  SLAVE_PROCS
    4
    1
  MASTER_PROCS
    1
  MPIEXEC
    ${_rescoup_mpiexec}
  SIMULATOR
    flow
  DEV_SIMULATOR
    flow_blackoil
  ABS_TOL
    ${abs_tol}
  REL_TOL
    ${rel_tol}
  TEST_ARGS
    --parsing-strictness=low
    --enable-tuning=true
)
