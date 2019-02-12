#! /bin/bash
#
# Runs a test from the test directory and compare the resulting VTU files.
#
# Usage:
#
# runTest.sh TEST_TYPE [TEST_ARGS]
#
MY_DIR="$(dirname "$0")"

usage() {
    echo "Usage:"
    echo
    echo "runTest.sh TEST_TYPE [TEST_ARGS]"
    echo "where TEST_TYPE can either be --plain, --simulation, --spe1 or --parallel-simulation=\$NUM_CORES (is '$TEST_TYPE')."
};

# this function clips the help message printed by an ewoms simulation
# to what is actually printed, throwing away all garbage which is
# printed before or after the "meat"
clipToHelpMessage()
{
    STATUS="not started"
    while read CUR_LINE; do
        if echo $CUR_LINE | grep -q "Usage: "; then
            STATUS="started"
        elif test "$STATUS" = "started" && echo $CUR_LINE | grep -q "^--"; then
            STATUS="params"
        elif test "$STATUS" = "params" && echo $CUR_LINE | grep -q "^[^-]"; then
            STATUS="stopped"
        fi

        if test "$STATUS" != "not started" && test "$STATUS" != "stopped"; then
            echo "$CUR_LINE"
        fi
    done
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

if test "$TEST_TYPE" != "--spe1"; then
    # find the binary in the its folder
    TEST_BINARY=$(find . -type f -perm -0111 -name "$TEST_NAME")
    NUM_BINARIES=$(echo "$TEST_BINARY" | wc -w | tr -d '[:space:]')

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
fi

#run the test
echo "######################"
echo "# Running test '$TEST_NAME'"
echo "######################"


RND="$(dd if=/dev/urandom bs=20 count=1 2> /dev/null | md5sum | cut -d" " -f1)"
case "$TEST_TYPE" in
    "--simulation")
        echo "executing \"$TEST_BINARY $TEST_ARGS\""
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

        SIM_NAME=$(grep "Applying the initial solution of the" "test-$RND.log" | sed "s/.*\"\(.*\)\".*/\1/" | head -n1)
        NUM_TIMESTEPS=$(( $(grep "Time step [0-9]* done" "test-$RND.log" | wc -l)))
        TEST_RESULT=$(printf "%s-%05i" "$SIM_NAME" "$NUM_TIMESTEPS")
        TEST_RESULT=$(ls -- "$TEST_RESULT".*)
        rm "test-$RND.log"
        if ! test -r "$TEST_RESULT"; then
            echo "File $TEST_RESULT does not exist or is not readable"
            exit 1
        fi

        echo "Simulation name: '$SIM_NAME'"
        echo "Number of timesteps: '$NUM_TIMESTEPS'"

        exit 0
        ;;

    "--parallel-simulation="*)
        NUM_PROCS="${TEST_TYPE/--parallel-simulation=/}"

        echo "executing \"mpirun -np \"$NUM_PROCS\" $TEST_BINARY $TEST_ARGS\""
        mpirun -np "$NUM_PROCS" "$TEST_BINARY" $TEST_ARGS | tee "test-$RND.log"
        RET="${PIPESTATUS[0]}"
        if test "$RET" != "0"; then
            echo "Executing the binary failed!"
            rm "test-$RND.log"
            exit 1
        fi

        SIM_NAME=$(grep "Applying the initial solution of the" "test-$RND.log" | sed "s/.*\"\(.*\)\".*/\1/" | head -n1)
        NUM_TIMESTEPS=$(( $(grep "Time step [0-9]* done" "test-$RND.log" | wc -l)))
        rm "test-$RND.log"

        echo "Simulation name: '$SIM_NAME'"
        echo "Number of timesteps: '$NUM_TIMESTEPS'"
        for PROC_NUM in 0 1 2 3; do
            REF_FILE=$(printf "s%04d-p%04d-%s" "$NUM_PROCS" "$PROC_NUM" "$SIM_NAME")
            TEST_RESULT=$(printf "s%04d-p%04d-%s-%05i" "$NUM_PROCS" "$PROC_NUM" "$SIM_NAME" "$NUM_TIMESTEPS")
            TEST_RESULT=$(ls -- "$TEST_RESULT".*)
            if ! test -r "$TEST_RESULT"; then
                echo "File $TEST_RESULT does not exist or is not readable"
                exit 1
            fi
        done
        exit 0
        ;;

    "--spe1")
        echo "Running the ebos simulator for SPE1CASE1"

        EBOS_COMMAND=$(find . -type f -perm -0111 -name "ebos")

        NUM_BINARIES=$(echo "$EBOS_COMMAND" | wc -w | tr -d '[:space:]')
        if test "$NUM_BINARIES" != "1"; then
            echo "No ebos executable found (is: $EBOS_COMMAND)"
            echo
            usage
            exit 1
        fi

        #########
        # Run the simulator
        if ! "$EBOS_COMMAND" "data/SPE1CASE1" ; then
            exit 1
        fi
        #########

        # make sure that the simulator produced non-empty output files
        for EXT in "INIT" "EGRID" "UNRST" "UNSMRY" "SMSPEC"; do
            if ! test -s "SPE1CASE1.$EXT"; then
                echo "The simulator did not produce a non-empty SPE1CASE1.$EXT file"
                exit 1
            fi
        done

        exit 0
        ;;

    "--restart")
        echo "executing \"$TEST_BINARY $TEST_ARGS\""
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
        HELP_MSG="$($TEST_BINARY --help | clipToHelpMessage)"
        if test "$(echo "$HELP_MSG" | grep -i usage)" == ''; then
            echo "$TEST_BINARY did not accept '--help' parameter"
            exit 1
        fi
        HELP_MSG2="$($TEST_BINARY -h | clipToHelpMessage)"
        if test "$HELP_MSG" != "$HELP_MSG2"; then
            echo "Output of $TEST_BINARY different when passing '--help' and '-h'"
            exit 1
        fi

        cat > "paramfile-$RND.ini" <<EOF
EndTime=100   
  
InitialTimeStepSize=100   # first in-line comment

UndefinedParam  =  "blubb"; second in-line comment
  # full line comment
; also a full line comment
EOF
        if ! $TEST_BINARY --parameter-file="paramfile-$RND.ini" 2>&1 > /dev/null; then
            echo "$TEST_BINARY does not correctly read a parameter file"
            exit 1
        elif $TEST_BINARY --parameter-file="foobar.ini" 2>&1 > /dev/null; then
            echo "$TEST_BINARY does not abort even though the specified parameter file does not exist"
            exit 1
        elif $TEST_BINARY --foo 2>&1 > /dev/null; then
            echo "$TEST_BINARY accepts flag parameters besides --help"
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
        echo "executing \"$TEST_BINARY $TEST_ARGS\""
        if ! "$TEST_BINARY" $TEST_ARGS; then
            exit 1
        fi

        exit 0
        ;;
esac
