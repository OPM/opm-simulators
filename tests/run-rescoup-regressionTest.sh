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
# Directories given with -d (data files that several cases share, or the deck
# directories of slaves that several cases use) are copied together with the
# case directory.  The copies keep their paths relative to the root of the test
# data (-o), so that the decks can include their files with the same relative
# paths as in the test data, whatever the depth of the case directory.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to the case directory (holding the master's and the slaves' deck directories)"
  echo -e "\t\t -o <path>     Root of the test data (opm-tests), which the case directory and the -d directories are under"
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
  echo -e "\t\t -n <procs>    Number of MPI processes of the master (default 1); the"
  echo -e "\t\t               master spawns the slaves with the ranks given in SLAVES"
  echo -e "\t\t -u <filename> Simulator binary for reference data"
  echo -e "\t\t -d <path>     A directory under the root of the test data to copy together with the case directory, for"
  echo -e "\t\t               data files that several cases share. Its path relative to the root is kept in the copy."
  echo -e "\t\t               Can be given several times. A missing directory is skipped"
  exit 1
fi

MASTER_PROCS=1
OPTIND=1
SHARED_DIRS=()
while getopts "i:o:m:r:f:s:a:t:c:x:e:n:u:d:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    o) TESTS_ROOT=${OPTARG%/} ;;
    m) MASTER_DIR=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    s) RESULT_NAMES=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) COMPARE_ECL_COMMAND=${OPTARG} ;;
    x) MPIEXEC=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    n) MASTER_PROCS=${OPTARG} ;;
    u) REF_EXE_NAME=${OPTARG} ;;
    d) SHARED_DIRS+=("${OPTARG}") ;;
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

# The case directory and the shared directories are copied to the same paths
# relative to the root of the test data under ${RESULT_PATH}/tree, so that the
# decks find each other's files with the same relative paths as in the test data.
# The reference results stay where they are.
if [ -z "${TESTS_ROOT}" ]; then
    echo "ERROR: The root of the test data (-o) is needed"
    exit 1
fi
under_root() {
  case "$1" in
    "${TESTS_ROOT}"/*) return 0 ;;
    *) return 1 ;;
  esac
}
if ! under_root "${INPUT_DATA_PATH}"; then
    echo "ERROR: The case directory ${INPUT_DATA_PATH} is not under the root of the test data ${TESTS_ROOT}"
    exit 1
fi

rm -rf ${RESULT_PATH}
TREE_PATH=${RESULT_PATH}/tree
CASE_PATH=${TREE_PATH}/${INPUT_DATA_PATH#${TESTS_ROOT}/}
mkdir -p ${CASE_PATH}
tar -C ${INPUT_DATA_PATH} --exclude=opm-simulation-reference -cf - . | \
  tar -C ${CASE_PATH} -xf -
test $? -eq 0 || exit 1

for SHARED_DIR in "${SHARED_DIRS[@]}"
do
  if ! under_root "${SHARED_DIR}"; then
    echo "ERROR: The shared directory ${SHARED_DIR} is not under the root of the test data ${TESTS_ROOT}"
    exit 1
  fi
  if [ -d "${SHARED_DIR}" ]; then
    SHARED_COPY=${TREE_PATH}/${SHARED_DIR#${TESTS_ROOT}/}
    mkdir -p ${SHARED_COPY}
    tar -C ${SHARED_DIR} --exclude=opm-simulation-reference -cf - . | tar -C ${SHARED_COPY} -xf -
    test $? -eq 0 || exit 1
  else
    echo "Skipping the shared directory ${SHARED_DIR}: it does not exist"
  fi
done

cd ${CASE_PATH}/${MASTER_DIR} || exit 1
"${MPIEXEC}" -np ${MASTER_PROCS} "${EXE_NAME}" ${FILENAME}.DATA ${TEST_ARGS} --output-dir=${RESULT_PATH}
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
