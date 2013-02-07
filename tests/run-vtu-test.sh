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
    echo "where TEST_TYPE can either be --plain or --simulation (is '$TEST_TYPE')."
};

function validateResults() {
    OUTPUT_FILE="$1"
    SIM_NAME="$2"

    for REFERENCE_RESULT in referencesolutions/$SIM_NAME*; do
        echo "Comparing with \"$REFERENCE_RESULT\"... "
        if python bin/fuzzycomparevtu.py "$REFERENCE_RESULT" "$OUTPUT_FILE"; then
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
        TEST_RESULT=$(ls $TEST_RESULT.*)
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

    "--simulation-diffusion")
        if ! "$TEST_BINARY" $TEST_ARGS ; then
            echo "Executing the binary failed!"
            exit 1
        fi
        validateResults mimeticdiffusion-00001.vtu mimeticdiffusion
        exit 0
        ;;
    
    "--parallel-simulation")
        mpirun -np 4 "$TEST_BINARY" $TEST_ARGS | tee "test-$RND.log"
        RET="${PIPESTATUS[0]}"
        if test "$RET" != "0"; then
            echo "Executing the binary failed!"
            rm "test-$RND.log"
            exit 1
        fi
        rm "test-$RND.log"

        # # compare the results
        # echo "######################"
        # echo "# Comparing results"
        # echo "######################"
        # SIM_NAME=$(grep "Writing result file for" test-$RND.log | sed "s/.*\"\(.*\)\".*/\1/" | head -n1)
        # TEST_RESULT=$(ls $SIM_NAME*.vtu $SIM_NAME*.vtp 2> /dev/null | sort | tail -n 1)
        # rm "test-$RND.log"
        # if ! test -r "$TEST_RESULT"; then
        #     echo "File $TEST_RESULT does not exist or is not readable"
        #     exit 1
        # fi

        # REFERENCE_RESULT="referencesolutions/$TEST_RESULT"
        # if ! test -r "$REFERENCE_RESULT"; then
        #     echo "File $REFERENCE_RESULT does not exist or is not readable"
        #     exit 1
        # fi
        

        # if ! python bin/fuzzycomparevtu.py "$REFERENCE_RESULT" "$TEST_RESULT"; then
        #     echo "The files \"$TEST_RESULT\" and \"$REFERENCE_RESULT\" are different."
        #     echo "Make sure the contents of \"$TEST_RESULT\" are still valid and "
        #     echo "make it the reference result if necessary."
        #     exit 1
        # fi
        
        # # SUCCESS!!!!!!
        # echo "Result and reference result are identical" 
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
        if test "$(echo $HELP_MSG | grep -i usage)" == ''; then
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
        if ! $TEST_BINARY --parameter-file="paramfile-$RND.ini" > /dev/null; then
            echo "$TEST_BINARY does not correctly read a parameter file"
            exit 1
        elif $TEST_BINARY --parameter-file="foobar.ini" > /dev/null; then
            echo "$TEST_BINARY does not abort even though the specified parameter file does not exist"
            exit 1
        elif ! $TEST_BINARY --foo --end-time=1 > /dev/null; then
            echo "$TEST_BINARY des not accept a flag parameters"
            exit 1
        fi

        # test some invalid parameter names
        for PARAM in foo -- -0foo --0foo --foo--bar --foo- -foo --foo-bar§=abc ; do
            if $TEST_BINARY "$PARAM" --end-time=100 > /dev/null; then
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
