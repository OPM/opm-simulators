#!/bin/bash

# Generates a summary plot comparison for each failed test case

OPM_TESTS_ROOT=$1
BUILD_DIR=$2
RESULT_DIR=$3
SOURCE_DIR=`dirname "$0"`
JOBS=${BUILDTHREADS:-16}

# ctest only writes this log when tests fail
FAILED_TESTS=$(cat $BUILD_DIR/Testing/Temporary/LastTestsFailed*.log 2>/dev/null)

mkdir -p $BUILD_DIR/failure_report
cd $BUILD_DIR/failure_report
rm -f *

JOBLIST=""
for failed_test in $FAILED_TESTS
do
  if grep -q -E "compareECLFiles" <<< $failed_test
  then
    failed_test=$(echo "${failed_test}" | sed -e 's/.*://')
    # Extract test properties
    binary=$(awk -v test="${failed_test}" -v prop="SIMULATOR" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    dir_name=$(awk -v test="${failed_test}" -v prop="DIRNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    file_name=$(awk -v test="${failed_test}" -v prop="FILENAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    test_name=$(awk -v test="${failed_test}" -v prop="TESTNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    ref_sim=$(awk -v test="${failed_test}" -v prop="REFERENCE_SIMULATOR" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    JOBLIST+="-r $OPM_TESTS_ROOT/${dir_name}/opm-simulation-reference/${ref_sim:-$binary}/${file_name} -s $RESULT_DIR/tests/results/$binary+$test_name/$file_name -c $test_name -o plot\\n"
  elif grep -q -E "compareSeparateECLFiles" <<< $failed_test
  then
    failed_test=$(echo "${failed_test}" | sed -e 's/.*://')
    # Extract test properties
    binary=$(awk -v test="${failed_test}" -v prop="SIMULATOR" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    dir_name=$(awk -v test="${failed_test}" -v prop="DIRNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    file_name1=$(awk -v test="${failed_test}" -v prop="FILENAME1" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    file_name2=$(awk -v test="${failed_test}" -v prop="FILENAME2" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    test_name=$(awk -v test="${failed_test}" -v prop="TESTNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
    JOBLIST+="-r $RESULT_DIR/tests/results/$binary+$test_name/$file_name1 -s $RESULT_DIR/tests/results/$binary+$test_name/$file_name2 -c $test_name -o plot -t $file_name1 -u $file_name2\\n"
  fi
done

if test -n "$JOBLIST"
then
  echo -e $JOBLIST | xargs -L1 -P${JOBS} $SOURCE_DIR/plot_well_comparison.py
  $SOURCE_DIR/plot_well_comparison.py  -o rename
else
  echo "No failed regression tests to report"
fi
