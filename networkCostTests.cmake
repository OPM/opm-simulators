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
