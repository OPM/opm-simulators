#! /bin/bash
#
# Runs a test from the test directory and compare the resulting VTU files.
#
# Usage:
#
# runTest.sh REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY TEST_ARGS
#
MY_DIR="$(dirname "$0")"

usage() {
    echo "Usage:"
    echo
    echo "runTest.sh TEST_TYPE TEST_BINARY [TEST_ARGS]"
    echo "where TEST_TYPE can either be --plain or --simulation=\$NUM_CORES (is '$TEST_TYPE')."
};

validateResults() {
    local OUTPUT_FILE="$1"
    local SIM_NAME="$2"

    for REFERENCE_RESULT in ${MY_DIR}/../tests/referencesolutions/$SIM_NAME*; do
        echo "Comparing with \"$REFERENCE_RESULT\"... "
        if python "${MY_DIR}/fuzzycomparevtu.py" "$REFERENCE_RESULT" "$OUTPUT_FILE"; then
            # SUCCESS!!!!!!
            echo "Result file '$OUTPUT_FILE' and reference '$REFERENCE_RESULT' are identical" 
            return 0
        fi
    done
    
    echo "There are no reference results which are are identical to \"$TEST_RESULT\"."
    echo "Make sure the contents of \"$TEST_RESULT\" are still valid and "
    echo "if necessary, add a reference result."
    exit 1
}


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
NUM_BINARIES=$(echo "$TEST_BINARY" | wc -w)


if test "$NUM_BINARIES" != "1"; then
    echo "No binary file found or binary file is non-unique (is: $TEST_BINARY)"
    echo
    usage
    exit 1
fi

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


RND="$(dd if=/dev/urandom bs=20 count=1 2> /dev/null | md5sum | cut -d" " -f1)"
case "$TEST_TYPE" in
    "--simulation")
        "$TEST_BINARY" $TEST_ARGS | tee "test-$RND.log"
        RET="${PIPESTATUS[0]}"
        if test "$RET" != "0"; then
            echo "Executing the binary failed!"
            rm "test-$RND.log"
            exit 1
        fi

        # compare the results
        echo "######################"
        echo "# Comparing results"
        echo "######################"
        echo "RND: '$RND'"

        SIM_NAME=$(grep "Initializing problem" "test-$RND.log" | sed "s/.*\"\(.*\)\".*/\1/" | head -n1)
        NUM_TIMESTEPS=$(grep "Writing result" "test-$RND.log" | wc -l)
        TEST_RESULT=$(printf "%s-%05i" "$SIM_NAME" "$NUM_TIMESTEPS")
        TEST_RESULT=$(ls "$TEST_RESULT".*)
        rm "test-$RND.log"
        if ! test -r "$TEST_RESULT"; then
            echo "File $TEST_RESULT does not exist or is not readable"
            exit 1
        fi

        echo "Simulation name: '$SIM_NAME'"
        echo "Number of timesteps: '$NUM_TIMESTEPS'"
        echo "Test result file: '$TEST_RESULT'"

        validateResults "$TEST_RESULT" "$SIM_NAME"
        exit 0
        ;;

    "--parallel-simulation="*)
        NUM_PROCS="${TEST_TYPE/--parallel-simulation=/}"

        mpirun -np "$NUM_PROCS" "$TEST_BINARY" $TEST_ARGS | tee "test-$RND.log"
        RET="${PIPESTATUS[0]}"
        if test "$RET" != "0"; then
            echo "Executing the binary failed!"
            rm "test-$RND.log"
            exit 1
        fi

        grep "Initializing problem" "test-$RND.log"
        SIM_NAME=$(grep "Initializing problem" "test-$RND.log" | sed "s/.*\"\(.*\)\".*/\1/" | head -n1)
        NUM_TIMESTEPS=$(grep "Writing result" "test-$RND.log" | wc -l)
        rm "test-$RND.log"

        echo "Simulation name: '$SIM_NAME'"
        echo "Number of timesteps: '$NUM_TIMESTEPS'"
        for PROC_NUM in 0 1 2 3; do
            REF_FILE=$(printf "s%04d-p%04d-%s" "$NUM_PROCS" "$PROC_NUM" "$SIM_NAME")
            TEST_RESULT=$(printf "s%04d-p%04d-%s-%05i" "$NUM_PROCS" "$PROC_NUM" "$SIM_NAME" "$NUM_TIMESTEPS")
            TEST_RESULT=$(ls "$TEST_RESULT".*)
            if ! test -r "$TEST_RESULT"; then
                echo "File $TEST_RESULT does not exist or is not readable"
                exit 1
            fi

            echo "Validate result for process $PROC_NUM (file: $TEST_RESULT)"

            validateResults "$TEST_RESULT" "$REF_FILE"
        done
        exit 0
        ;;

    "--restart")
        "$TEST_BINARY" $TEST_ARGS | tee "test-$RND.log"
        RET="${PIPESTATUS[0]}"
        if test "$RET" != "0"; then
            echo "Executing the binary failed!"
            rm "test-$RND.log"
            exit 1
        fi
        RESTART_TIME=$(grep "Serialize" "test-$RND.log" | tail -n 1 | sed "s/.*time=\([0-9.e+\-]*\).*/\1/")
        rm "test-$RND.log"
        
        if ! "$TEST_BINARY" $TEST_ARGS --restart-time="$RESTART_TIME" --newton-write-convergence=true; then
            echo "Restarting $TEST_BINARY failed"
            exit 1;
        fi
        exit 0
        ;;        

    "--parameters")
        HELP_MSG="$($TEST_BINARY --help)"
        if test "$(echo "$HELP_MSG" | grep -i usage)" == ''; then
            echo "$TEST_BINARY did not accept '--help' parameter"
            exit 1
        fi
        HELP_MSG2="$($TEST_BINARY -h)"
        if test "$HELP_MSG" != "$HELP_MSG2"; then
            echo "Output of $TEST_BINARY different when passing '--help' and '-h'"
            exit 1
        fi

        cat > "paramfile-$RND.ini" <<EOF
EndTime=100
InitialTimeStepSize=100
UndefinedParam="blubb"
EOF
        if ! $TEST_BINARY --parameter-file="paramfile-$RND.ini" 2>&1 > /dev/null; then
            echo "$TEST_BINARY does not correctly read a parameter file"
            exit 1
        elif $TEST_BINARY --parameter-file="foobar.ini" 2>&1 > /dev/null; then
            echo "$TEST_BINARY does not abort even though the specified parameter file does not exist"
            exit 1
        elif ! $TEST_BINARY --foo --end-time=1 2>&1 > /dev/null; then
            echo "$TEST_BINARY does not accept a flag parameters"
            exit 1
        fi

        # test some invalid parameter names
        for PARAM in foo -- -0foo --0foo --foo--bar --foo- -foo --foo-bar§=abc ; do
            if $TEST_BINARY "$PARAM" --end-time=100 2>&1 > /dev/null; then
                echo "$TEST_BINARY accepted invalid command line option '$PARAM'"
                exit 1
            fi
        done
        echo "Test successful"
        exit 0

        ;;

    "--plain")
        if ! "$TEST_BINARY" $TEST_ARGS; then
            exit 1
        fi

        exit 0
        ;;
esac
