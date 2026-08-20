#!/bin/bash

# This performs a serial and a parallel for a simulator,
# then compares the summary and restart files from the two runs.
# Meant to track regression in parallel simulators.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\tOptional options:"
  echo -e "\t\t -n <procs>    Number of MPI processes to use"
  echo -e "\t\t -m <mode>     Comparison mode: summary (default, compares SMRY+UNRST"
  echo -e "\t\t               against a normal run) or init (dry-run, compares"
  echo -e "\t\t               EGRID+INIT only -- for tracking parallel INIT-file"
  echo -e "\t\t               regressions, e.g. parallel LGR transmissibility)"
  exit 1
fi

MPI_PROCS=4
MODE=summary
OPTIND=1

while getopts "i:r:f:a:t:c:e:n:m:" OPT
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
    m) MODE=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

case "${MODE:-summary}" in
  summary|init) ;;
  *) echo "Unknown mode '${MODE}' (expected 'summary' or 'init')"; exit 1 ;;
esac

if [ "${MODE}" = "init" ]
then
  RUN_FLAG="--enable-dry-run=true --enable-ecl-output=true"
else
  RUN_FLAG="--enable-opm-rst-file=true"
fi

rm -Rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
"${EXE_NAME}" ${TEST_ARGS} ${RUN_FLAG} --output-dir=${RESULT_PATH}

test $? -eq 0 || exit 1
mkdir mpi
cd mpi
mpirun -np ${MPI_PROCS} "${EXE_NAME}" ${TEST_ARGS} ${RUN_FLAG} --output-dir=${RESULT_PATH}/mpi
test $? -eq 0 || exit 1
cd ..

ecode=0
if [ "${MODE}" = "init" ]
then
  echo "=== Executing comparison for EGRID file ==="
  ${COMPARE_ECL_COMMAND} -t EGRID ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -t EGRID ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  fi

  # -x ignores the parallel-only MPI_RANK keyword, which has no serial counterpart to compare against.
  echo "=== Executing comparison for INIT file (ignoring parallel-only MPI_RANK) ==="
  ${COMPARE_ECL_COMMAND} -t INIT -x ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -t INIT -x ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  fi
else
  echo "=== Executing comparison for summary file ==="
  ${COMPARE_ECL_COMMAND} -t SMRY -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} -t SMRY -a -R ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  fi

  echo "=== Executing comparison for restart file ==="
  ${COMPARE_ECL_COMMAND} -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} -a -l -t UNRST ${RESULT_PATH}/${FILENAME} ${RESULT_PATH}/mpi/${FILENAME} ${ABS_TOL} ${REL_TOL}
  fi
fi

exit $ecode
