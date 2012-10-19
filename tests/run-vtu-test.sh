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
    echo "runTest.sh TEST_TYPE TEST_BINARY [TEST_ARGS]"
    echo "where TEST_TYPE can either be --plain or --simulation"
};

TEST_TYPE="$1"
TEST_NAME="$2"
TEST_ARGS="${@:3:100}"

# make sure we have at least 2 parameters
if test "$#" -lt 2; then
    echo "Wrong number of parameters"
    echo
    usage
    exit 1
fi

# find the binary in the its folder
TEST_BINARY=$(find -type f -executable -name "$TEST_NAME")
    
# make sure the binary is of the test is present
if ! test -x "$TEST_BINARY"; then
    echo "$TEST_NAME does not exist or is not executable"
    echo
    usage
    exit 1
fi

#run the test
echo "######################"
echo "# Running test '$TEST_NAME'"
echo "######################"

if ! "$TEST_BINARY" $TEST_ARGS; then
    echo "Executing the binary failed!"
    exit 1
fi

case "$TEST_TYPE" in
    "--simulation")

        # compare the results
        echo "######################"
        echo "# Comparing results"
        echo "######################"
        TEST_RESULT=$(ls --sort=time *.vtu *.vtp | head -n 1 2> /dev/null)
        if ! test -r "$TEST_RESULT"; then
            echo "File $TEST_RESULT does not exist or is not readable"
            exit 1
        fi

        REFERENCE_RESULT="referencesolutions/$TEST_RESULT"
        if ! test -r "$REFERENCE_RESULT"; then
            echo "File $REFERENCE_RESULT does not exist or is not readable"
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

        ;;
    *)
        exit 0
        ;;
esac
