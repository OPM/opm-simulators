#! /bin/bash
#
# Runs a test from the test directory and compare the resulting VTU files.
#
# Usage:
#
# runTest.sh REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY TEST_ARGS
#

function usage() {
    echo "Usage:"
    echo
    echo "runTest.sh REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY [TEST_ARGS]"
};

REFERENCE_RESULT="$1"
TEST_RESULT="$2"
TEST_BINARY="$3"
TEST_ARGS="${@:4:100}"

# make sure we have at least 3 parameters
if test "$#" -lt 3; then
    echo "Wrong number of parameters"
    echo
    usage
    exit 1
fi

# make sure the reference result exists
if ! test -r "$REFERENCE_RESULT"; then
    echo "File $REFERENCE_RESULT does not exist or is not readable"
    echo
    usage
    exit 1
fi
    
# make sure the binary is of the test is present
if ! test -x "$TEST_BINARY"; then
    echo "$TEST_BINARY does not exist or is not executable"
    echo
    usage
    exit 1
fi

#run the test
echo "######################"
echo "# Running test"
echo "######################"
if ! "$TEST_BINARY" $TEST_ARGS; then
    echo "Executing the binary failed!"
    exit 1
fi

# compare the results
echo "######################"
echo "# Comparing results"
echo "######################"
if ! test -r "$TEST_RESULT"; then
    echo "File $TEST_RESULT does not exist or is not readable"
    exit 1
fi

if ! python bin/fuzzycomparevtu.py "$REFERENCE_RESULT" "$TEST_RESULT"; then
    echo "The files \"$TEST_RESULT\" and \"$REFERENCE_RESULT\" are different."
    echo "Make sure the contents of \"$TEST_RESULT\" are still valid and "
    echo "make it the reference result if necessary."
    exit 1
fi

# SUCCESS!!!!!!
echo "Result and reference result are identical" 
exit 0
