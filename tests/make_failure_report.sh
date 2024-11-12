#!/bin/bash

# Generates a summary plot comparison for each failed test case

OPM_TESTS_ROOT=$1
BUILD_DIR=$2
RESULT_DIR=$3
SOURCE_DIR=`dirname "$0"`

FAILED_TESTS=`cat $BUILD_DIR/Testing/Temporary/LastTestsFailed*.log`

mkdir -p $BUILD_DIR/failure_report
cd $BUILD_DIR/failure_report
rm -f *

for failed_test in $FAILED_TESTS
do
  grep -q -E "compareECLFiles" <<< $failed_test
  test $? -ne 0 && continue

  failed_test=`echo $failed_test | sed -e 's/.*://g' -e 's/\+/./g'`
  # Extract test properties
  binary=$(awk -v search="set_tests_properties\\\($failed_test\$" -v prop="SIMULATOR" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
  dir_name=$(awk -v search="set_tests_properties\\\($failed_test\$" -v prop="DIRNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
  file_name=$(awk -v search="set_tests_properties\\\($failed_test\$" -v prop="FILENAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
  test_name=$(awk -v search="set_tests_properties\\\($failed_test\$" -v prop="TESTNAME" -f ${SOURCE_DIR}/getprop.awk $RESULT_DIR/CTestTestfile.cmake)
  echo "Processing ${test_name}"
  $SOURCE_DIR/plot_well_comparison.py -r $OPM_TESTS_ROOT/$dir_name/opm-simulation-reference/$binary/$file_name -s $RESULT_DIR/tests/results/$binary+$test_name/$file_name -c $test_name -o plot
done

if test -n "$FAILED_TESTS"
then
  $SOURCE_DIR/plot_well_comparison.py  -o rename
fi
