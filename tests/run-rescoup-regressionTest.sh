#!/bin/bash

# This runs a reservoir coupling simulation, where the master spawns its
# slaves, then compares the summary, restart and init files of the master
# and of every slave against a reference.
#
# The case directory is copied into the result directory first, and the
# master is run from its copy. The slaves' deck directories are given
# relative to the master's working directory (SLAVES item 4), and the
# slaves' log files are written there, so running from the copy keeps the
# test from writing into the test data.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to the case directory (holding the master's and the slaves' deck directories)"
  echo -e "\t\t -m <dir>      The master deck's sub-directory of the case directory"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -f <filename> Master deck file name"
  echo -e "\t\t -s <names>    Comma-separated base names of all result files to compare (master and slaves)"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -x <path>     MPI launcher"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\tOptional options:"
  echo -e "\t\t -u <filename> Simulator binary for reference data"
  exit 1
fi

OPTIND=1
while getopts "i:m:r:f:s:a:t:c:x:e:u:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    m) MASTER_DIR=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    s) RESULT_NAMES=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    x) MPIEXEC=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    u) REF_EXE_NAME=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

BASE_EXE_NAME=$(basename -- "${EXE_NAME}")
REF_EXE_NAME=${REF_EXE_NAME:-${BASE_EXE_NAME}}

if [ ! -x "${EXE_NAME}" ]; then
    echo "ERROR: Simulator binary not found: ${EXE_NAME}"
    echo ""
    echo "To build this binary, run one of:"
    echo "  ninja ${BASE_EXE_NAME}"
    echo "  make ${BASE_EXE_NAME}"
    echo "  cmake --build . --target ${BASE_EXE_NAME}"
    exit 1
fi

rm -rf ${RESULT_PATH}
mkdir -p ${RESULT_PATH}/case
tar -C ${INPUT_DATA_PATH} --exclude=opm-simulation-reference -cf - . | \
  tar -C ${RESULT_PATH}/case -xf -
test $? -eq 0 || exit 1

cd ${RESULT_PATH}/case/${MASTER_DIR} || exit 1
"${MPIEXEC}" -np 1 "${EXE_NAME}" ${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}
test $? -eq 0 || exit 1
cd ${RESULT_PATH}

ignore_extra_kw=""
if grep -q "ignore_extra" <<< $ghprbCommentBody
then
    ignore_extra_kw="-x"
fi

type=""
if grep -q "only_summary" <<< $ghprbCommentBody
then
  type="-t SMRY"
  echo "=== Executing comparison for UNSMRY files if these exists in reference folder ==="
else
  echo "=== Executing comparison for EGRID, INIT, UNRST, UNSMRY and RFT files if these exists in reference folder ==="
fi

ecode=0
for name in ${RESULT_NAMES//,/ }
do
  echo "=== Comparing ${name} ==="
  ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${type} ${INPUT_DATA_PATH}/opm-simulation-reference/${REF_EXE_NAME}/${name} ${RESULT_PATH}/${name} ${ABS_TOL} ${REL_TOL}
  if [ $? -ne 0 ]
  then
    ecode=1
    ${COMPARE_ECL_COMMAND} ${ignore_extra_kw} ${type} -a ${INPUT_DATA_PATH}/opm-simulation-reference/${REF_EXE_NAME}/${name} ${RESULT_PATH}/${name} ${ABS_TOL} ${REL_TOL}
  fi
done

exit $ecode
