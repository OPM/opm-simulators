#!/bin/bash

# This runs a simulator, then compares the damaris output files
# against a reference.

if test $# -eq 0
then
  echo -e "Usage:\t$0 <options> -- [additional simulator options]"
  echo -e "\tMandatory options:"
  echo -e "\t\t -i <path>     Path to read deck from"
  echo -e "\t\t -r <path>     Path to store results in"
  echo -e "\t\t -b <path>     Path to simulator binary"
  echo -e "\t\t -f <filename> Deck file name"
  echo -e "\t\t -a <tol>      Absolute tolerance in comparison"
  echo -e "\t\t -t <tol>      Relative tolerance in comparison"
  echo -e "\t\t -c <path>     Path to comparison tool"
  echo -e "\t\t -e <filename> Simulator binary to use"
  echo -e "\tOptional options:"
  echo -e "\t\t -n <procs>    Number of MPI processes to use"
  exit 1
fi

MPI_PROCS=4
OPTIND=1
while getopts "i:r:b:f:a:t:c:e:n:" OPT
do
  case "${OPT}" in
    i) INPUT_DATA_PATH=${OPTARG} ;;
    r) RESULT_PATH=${OPTARG} ;;
    b) BINPATH=${OPTARG} ;;
    f) FILENAME=${OPTARG} ;;
    a) ABS_TOL=${OPTARG} ;;
    t) REL_TOL=${OPTARG} ;;
    c) H5DIFF_COMMAND=${OPTARG} ;;
    e) EXE_NAME=${OPTARG} ;;
    n) MPI_PROCS=${OPTARG} ;;
  esac
done
shift $(($OPTIND-1))
TEST_ARGS="$@"

mkdir -p ${RESULT_PATH}
cd ${RESULT_PATH}
mpirun -np ${MPI_PROCS} ${BINPATH}/${EXE_NAME} ${INPUT_DATA_PATH}/${FILENAME} ${TEST_ARGS} --output-dir=${RESULT_PATH} --damaris-dedicated-cores=1 --damaris-save-mesh-to-hdf=true --enable-damaris-output=true --enable-ecl-output=false --damaris-output-hdf-collective=1 --damaris-save-to-hdf=1 --damaris-sim-name="${FILENAME}"
test $? -eq 0 || exit 1
cd ..

echo "=== Executing comparison for files if these exists in reference folder ==="
for fname in `ls ${RESULT_PATH}/${FILENAME}*.h5`
do
  fname=`basename $fname`
  if test -f ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${fname}
  then
    echo -e "\t - ${fname}"
    ${H5DIFF_COMMAND} -d ${ABS_TOL} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${fname} ${RESULT_PATH}/${fname}
    test $? -eq 0 || exit 1
    ${H5DIFF_COMMAND} -p ${REL_TOL} ${INPUT_DATA_PATH}/opm-simulation-reference/${EXE_NAME}/${fname} ${RESULT_PATH}/${fname}
    test $? -eq 0 || exit 1
  fi
done
