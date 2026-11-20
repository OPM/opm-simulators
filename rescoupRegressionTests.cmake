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
  # Data files that rc_m1 shares with other cases
  SHARED_DIRS
    rescoup/common
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

# The rc_m2 cases: variants of one coupled model. Each master spawns RES-1
# (4 ranks) and RES-2 (1 rank). The cases take RES-2, and some take RES-1,
# from rc_m2/t2a, so that directory is copied with each case.
#   case_dir: the case's directory under rescoup/rc_m2
#   master:   the master deck (no extension)
#   res1:     the deck of RES-1 (RES-2 is always T2A_MOD2_FLOW)
# Further arguments are passed on, e.g. REFERENCE_SUBDIR.
function(add_rc_m2_test case_dir master res1)
  string(TOLOWER "${master}" master_lower)
  add_test_compare_rescoup_simulation(
    CASENAME
      rc_m2_${case_dir}_${master_lower}
    DIR
      rescoup/rc_m2/${case_dir}
    MASTER_DIR
      master
    FILENAME
      ${master}
    SLAVE_FILENAMES
      ${res1}
      T2A_MOD2_FLOW
    SLAVE_PROCS
      4
      1
    MASTER_PROCS
      1
    SHARED_DIRS
      rescoup/common
      rescoup/rc_m2/t2a
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
    ${ARGN}
  )
endfunction()

add_rc_m2_test(t2a TEST2A_FLOW T2A_MOD1_FLOW)
add_rc_m2_test(injcap TEST2A_IC_FLOW T2A_MOD1_FLOW)
add_rc_m2_test(netroot TEST2A_FLOW_NETROOT T2A_MOD1_FLOW)
add_rc_m2_test(slavenet TEST2A_SN_FLOW T2A_SN_MOD1_FLOW)
add_rc_m2_test(slavenet2 TEST2A_SN2_FLOW T2A_SN2_MOD1_FLOW)
add_rc_m2_test(slavenet3 TEST2A_SN3_FLOW T2A_SN3_MOD1_FLOW)
# A slave group with a gas limit of its own (GRUPSLAV flag BOTH), which it breaks while
# the master's target is an oil rate, and the oil flag BOTH without an oil limit of its own.
add_rc_m2_test(slavenet TEST2A_SN_OWNLIM_FLOW T2A_SN_OWNLIM_MOD1_FLOW
               REFERENCE_SUBDIR TEST2A_SN_OWNLIM_FLOW)

# Gas lift: the four masters share one case directory and spawn the same
# slave decks, so each keeps its reference files in its own sub-directory.
foreach(master TEST2A_GL_FLOW TEST2A_GL_MGL_LEAF_FLOW TEST2A_GL_MGL_MOD1_FLOW)
  add_rc_m2_test(gaslift ${master} T2A_GL_MOD1_FLOW REFERENCE_SUBDIR ${master})
endforeach()
add_rc_m2_test(gaslift TEST2A_GL_NOGRP_FLOW T2A_GL_NOGRP_MOD1_FLOW
               REFERENCE_SUBDIR TEST2A_GL_NOGRP_FLOW)
# The GRUPSLAV oil flag SLAV, while the master's targets change between oil and gas.
add_rc_m2_test(gaslift TEST2A_GL_SLAV_FLOW T2A_GL_SLAV_MOD1_FLOW
               REFERENCE_SUBDIR TEST2A_GL_SLAV_FLOW)
# The master withdraws the production target of a slave group that has a gas limit
# of its own but no oil limit.
add_rc_m2_test(gaslift TEST2A_GL_WD_OWN_FLOW T2A_GL_WD_OWN_MOD1_FLOW
               REFERENCE_SUBDIR TEST2A_GL_WD_OWN_FLOW)
# The master withdraws the water injection targets of all its groups, and a slave
# group has the water injection flag BOTH without a GCONINJE of its own.
add_rc_m2_test(gaslift TEST2A_GL_WD_INJ_FLOW T2A_GL_WD_INJ_MOD1_FLOW
               REFERENCE_SUBDIR TEST2A_GL_WD_INJ_FLOW)
