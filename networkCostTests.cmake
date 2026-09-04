# What the network solve costs, as opposed to what it computes.
#
# The regression tests compare the answer and would pass a run that reaches it
# by falling back to the relaxed update on every step, or by re-solving one
# well thousands of times because the network and the well model disagree
# about whether it flows. Both have happened on these two decks. The bounds
# are roughly 25 % above what the simultaneous solve does today; they are
# there to catch a regression, not to pin an exact count.

function(add_network_cost_test)
  set(oneValueArgs CASENAME FILENAME DIR MAX_FALLBACKS MAX_NEWTON MAX_WELL_SOLVES MAX_WELL_SHARE)
  set(multiValueArgs TEST_ARGS)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  add_test(NAME network_cost_${PARAM_CASENAME}
           COMMAND ${PROJECT_SOURCE_DIR}/tests/run-network-cost-test.sh
                   -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                   -f ${PARAM_FILENAME}
                   -r ${PROJECT_BINARY_DIR}/tests/results/network_cost_${PARAM_CASENAME}
                   -e ${PROJECT_BINARY_DIR}/bin/flow_blackoil
                   -F ${PARAM_MAX_FALLBACKS}
                   -N ${PARAM_MAX_NEWTON}
                   -W ${PARAM_MAX_WELL_SOLVES}
                   -S ${PARAM_MAX_WELL_SHARE}
                   -- ${PARAM_TEST_ARGS})
  set_tests_properties(network_cost_${PARAM_CASENAME} PROPERTIES LABELS "network_cost")
endfunction()

# An autochoke under the simultaneous solve. Measured 2026-08-24: 0 fallbacks,
# 199 Newton, 4603 well solves, B-1H 45 % of them. Legacy on the same deck
# needs 746 Newton and 70040 well solves for the same answer.
add_network_cost_test(
  CASENAME autochoke_complementarity
  FILENAME NETWORK_MODEL5_STDW_AUTOCHK
  DIR network
  MAX_FALLBACKS 0
  MAX_NEWTON 260
  MAX_WELL_SOLVES 6000
  MAX_WELL_SHARE 60
  TEST_ARGS --network-solver=newton --network-analytic-jacobian=true
            --network-group-control=true --network-autochoke=true
            --network-complementarity=true
)

# Gas lift answered by the network. Measured 2026-08-24: 0 fallbacks, 521
# Newton, 9051 well solves, B-3H 32 % of them. Legacy: 2607 and 49320.
add_network_cost_test(
  CASENAME gaslift_complementarity
  FILENAME GASLIFT-13
  DIR gaslift
  MAX_FALLBACKS 0
  MAX_NEWTON 680
  MAX_WELL_SOLVES 11500
  MAX_WELL_SHARE 45
  TEST_ARGS --network-solver=newton --network-analytic-jacobian=true
            --network-group-control=true --gas-lift-network-response=true
            --network-complementarity=true
)

# The simultaneous solve reaches the same balance as the relaxed update by a
# different route, so on a deck where the relaxed update converges the two must
# agree. Checkable without any external reference, which matters because these
# decks have none.
#
# The tolerance is looser than the regression suite's (which compares a run
# against its own stored reference): two solution paths through the same
# schedule differ at the sub-percent level. It still discriminates -- where the
# formulations genuinely disagree the gap is 4-6 % (the active set against
# legacy on autochoke, and the gas-lift decks, which is why those are covered
# by cost tests instead).
set(network_agreement_abs_tol 2e-2)
set(network_agreement_rel_tol 2e-2)

function(add_network_agreement_test)
  set(oneValueArgs CASENAME FILENAME DIR)
  set(multiValueArgs BOTH_ARGS NEWTON_ARGS)
  cmake_parse_arguments(PARAM "" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
  string(REPLACE ";" " " _both "${PARAM_BOTH_ARGS}")
  add_test(NAME network_agreement_${PARAM_CASENAME}
           COMMAND ${PROJECT_SOURCE_DIR}/tests/run-network-agreement-test.sh
                   -i ${OPM_TESTS_ROOT}/${PARAM_DIR}
                   -f ${PARAM_FILENAME}
                   -r ${PROJECT_BINARY_DIR}/tests/results/network_agreement_${PARAM_CASENAME}
                   -e ${PROJECT_BINARY_DIR}/bin/flow_blackoil
                   -c $<TARGET_FILE:compareECL>
                   -a ${network_agreement_abs_tol}
                   -t ${network_agreement_rel_tol}
                   -b "${_both}"
                   -- ${PARAM_NEWTON_ARGS})
  set_tests_properties(network_agreement_${PARAM_CASENAME}
                       PROPERTIES LABELS "network_agreement")
endfunction()

# Every deck where the production solve engages and the relaxed update also
# converges. Measured 2026-08-24 (legacy FOPT vs simultaneous, well solves):
#   NETWORK-01              538142 / 449   vs 538842 / 379
#   NETWORK-01_STANDARD     538142 / 449   vs 538842 / 379
#   NETWORK-01-REROUTE      543737 / 756   vs 544171 / 493
#   NETWORK-01-WTEST        674921 / 1016  vs 675875 / 549
#   NETWORK-01-WEFAC-...    128788 / 430   vs 128772 / 365
#   NETWORK-01-MULTIROOT    778152 / 407   vs 778188 / 380
#   6_UDA_MODEL5_STDW       317225 / 1781  vs 317225 / 1788  (exact)
set(_agree_newton --network-solver=newton --network-analytic-jacobian=true)

foreach(case NETWORK-01 NETWORK-01_STANDARD NETWORK-01-WTEST
             NETWORK-01-WEFAC-GEFAC-ITEM3NO NETWORK-01-MULTIROOT)
  string(TOLOWER ${case} _lc)
  string(REPLACE "-" "_" _lc ${_lc})
  add_network_agreement_test(
    CASENAME ${_lc}
    FILENAME ${case}
    DIR network
    BOTH_ARGS --enable-tuning=true
    NEWTON_ARGS ${_agree_newton}
  )
endforeach()

add_network_agreement_test(
  CASENAME network_01_reroute
  FILENAME NETWORK-01-REROUTE
  DIR network
  BOTH_ARGS --enable-tuning=true --local-well-solve-control-switching=true
  NEWTON_ARGS ${_agree_newton}
)

add_network_agreement_test(
  CASENAME uda_model5_stdw
  FILENAME 6_UDA_MODEL5_STDW
  DIR model5
  BOTH_ARGS --enable-tuning=true
  NEWTON_ARGS ${_agree_newton} --network-group-control=true
)
