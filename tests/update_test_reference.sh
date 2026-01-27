#!/bin/bash -norc

# Location of this script.
dir=$(dirname $0)

# Tests that are ultimately updated in this script run. Aggregated
# mostly in order to create the final commit message.
changed_tests=""

# Collection of properties for the test currently being processed.
# Populated by function extractTestProperties().
declare -A testProperty

# ===========================================================================
# Helper functions.

# Generate textual diffs of formatted file contents.
#
# Appends to ${WORKSPACE}/data_diff.
#
# $1 is location of reference results
# $2 is location of new simulation results
# $3 is name of result file.
compareResultFileContents () {
    local reference_dir="${1}"
    local simoutput_dir="${2}"
    local result_file="${3}"

    # Recall: convertECL prints a line of the form
    #
    #    converting path/to/CASE.UNRST -> path/to/CASE.FUNRST
    #
    # so we use this output to infer the result file name.

    local tmp2=$(mktemp -d)
    mkdir -p "${tmp2}"/{orig,new}

    if [ -s "${reference_dir}/${result_file}" ]
    then
        cp -a "${reference_dir}/${result_file}" "${tmp2}/orig"
        ref=$(${CONVERT_ECL} "${tmp2}/orig/${result_file}" | \
                  awk '/converting/{print $NF}')
    else
        # Reference file does not exist. New file type or new case.
        ref=/dev/null
    fi

    if [ -s "${simoutput_dir}/${result_file}" ]
    then
        cp -a "${simoutput_dir}/${result_file}" "${tmp2}/new"
        new=$(${CONVERT_ECL} "${tmp2}/new/${result_file}" | \
                  awk '/converting/{print $NF}')
    else
        # New simulation file does not exist. Typically when the file
        # type is removed. Unexpected.
        new=/dev/null
    fi


    {
        flock -x 200
        if [ "${new}" != "${ref}" ]
        then
            diff -u "${ref}" "${new}" >> "${WORKSPACE}/data_diff"
        else
            # Neither reference nor updated file exists. This really
            # shouldn't happen.
            echo "No difference between reference and update for ${result_file}" \
                 >> "${WORKSPACE}/data_diff"
        fi
    } 200>|$TMPDIR/opm_diff_lock

    rm -rf ${tmp2}
}

# Extract individual test properties of a failed test
#
# Input argument $1 is the test description from 'LastTestsFailed.log'.
#
# Populates 'testProperty' in such a way that
#
#   ${testProperty["binary"]}    is the simulation tool used to run the test.
#   ${testProperty["dir_name"]}  is the source directory in opm-tests.
#   ${testProperty["file_name"]} is the base name of the simulation input file.
#   ${testProperty["test_name"]} is the test suite's name of the particular test.
extractTestProperties () {
    local failed_test

    failed_test=$(echo "$1" | sed -e 's/.*://' -e 's/\+/./g')

    testProperty["binary"]=$(awk -v search="set_tests_properties\\\(.*${failed_test}.*\$" \
                             -v prop="SIMULATOR" -f "${dir}/getprop.awk" \
                             "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["dir_name"]=$(awk -v search="set_tests_properties\\\(.*${failed_test}.*\$" \
                               -v prop="DIRNAME" -f "${dir}/getprop.awk" \
                               "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["file_name"]=$(awk -v search="set_tests_properties\\\(.*${failed_test}.*\$" \
                                -v prop="FILENAME" -f "${dir}/getprop.awk" \
                                "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["test_name"]=$(awk -v search="set_tests_properties\\\(.*${failed_test}.*\$" \
                                -v prop="TESTNAME" -f "${dir}/getprop.awk" \
                                "${BUILD_DIR}/CTestTestfile.cmake")
}

# Copy results from a test run to reference dir
#
# $1 = source directory to copy data from
# $2 = destination directory to copy data to
# $3 = base file name for files to copy
# $4...$@ = file types to copy
copyToReferenceDir () {
    SRC_DIR=$1
    DST_DIR=$2
    STEM=$3
    FILETYPES=${@:4}
    mkdir -p "${DST_DIR}"

    DIFF=1
    for filetype in $FILETYPES
    do
        result_file="${STEM}.${filetype}"

        # Don't flag results as changed if neither the reference nor
        # the result directory have a specific file type. This means
        # we transparently handle result files like *.RFT which are
        # created only if the simulation model explicitly requests the
        # file type.
        if [ ! -f "${SRC_DIR}/${result_file}" -a ! -f "${DST_DIR}/${result_file}" ]
        then
            continue
        fi

        if ! cmp -s "${SRC_DIR}/${result_file}" "${DST_DIR}/${result_file}"
        then
            if  [ -n "$CONVERT_ECL" ]
            then
                compareResultFileContents "${DST_DIR}" "${SRC_DIR}" "${result_file}"
            fi

            cp -a "${SRC_DIR}/${result_file}" "${DST_DIR}"
            DIFF=0
        fi
    done

    return $DIFF
}

# Copy damaris results from a test run to reference dir
#
# $1 = source directory to copy data from
# $2 = destination directory to copy data to
# $3 = base file name for files to copy
copyDamarisToReferenceDir () {
  SRC_DIR=$1
  DST_DIR=$2
  STEM=$3
  mkdir -p $DST_DIR

  FIRST_FILE=`ls -v1 $SRC_DIR/$STEM*.h5 | head -n1`
  LAST_FILE=`ls -v1 $SRC_DIR/$STEM*.h5 | tail -n1`

  for file in $FIRST_FILE $LAST_FILE
  do
    {
        flock -x 200
        h5diff -v "$file" "$DST_DIR/`basename $file`" >> $WORKSPACE/data_diff
    } 200>|$TMPDIR/opm_diff_lock
    cp "$file" $DST_DIR
  done
}

# Generate reference solution updates for restarted simulation runs.
#
# $1 is name of simulator binary
# $2 is location of reference solutions in opm-tests
# $3 is base name of simulation input file
# $4 is test suite's name of test base run
updateRestartResults () {
    local binary="${1}"
    local dir_name="${2}"
    local file_name="${3}"
    local test_name="${4}"

    local rst_steps
    local rst_step
    local result
    local case_result

    rst_steps=$(echo ${BUILD_DIR}/tests/results/${binary}+${test_name}/restart/*RESTART_*.UNRST | \
                    sed -E -e 's;[^ ]*RESTART_([0-9]+)\.UNRST;\1;g')

    result=0
    for rst_step in ${rst_steps}
    do
        copyToReferenceDir \
            "${BUILD_DIR}/tests/results/${binary}+${test_name}/restart/" \
            "${OPM_TESTS_ROOT}/${dir_name}/opm-simulation-reference/${binary}/restart" \
            "${file_name}_RESTART_${rst_step}" \
            EGRID INIT RFT SMSPEC UNRST UNSMRY

        case_result=$?
        if [ ${result} -eq 0 ]
        then
            result=${case_result}
        fi
    done

    if [ ${result} -eq 0 ]
    then
        changed_tests+=" ${test_name}(restart)"
    fi
}

# Generate reference solution updates for base and restarted
# simulation runs.
#
# Uses properties of failing test being processed (${testProperty}).
updateFullSimulationResults () {
    local binary=${testProperty["binary"]}
    local dir_name=${testProperty["dir_name"]}
    local file_name=${testProperty["file_name"]}
    local test_name=${testProperty["test_name"]}

    if copyToReferenceDir \
           "${BUILD_DIR}/tests/results/${binary}+${test_name}" \
           "${OPM_TESTS_ROOT}/${dir_name}/opm-simulation-reference/${binary}" \
           "${file_name}" \
           EGRID INIT RFT SMSPEC UNRST UNSMRY
    then
        changed_tests+=" ${test_name}"
    fi

    if [ -d "${BUILD_DIR}/tests/results/${binary}+${test_name}/restart" ]
    then
        updateRestartResults "${binary}" "${dir_name}" "${file_name}" "${test_name}"
    fi
}

# Generate reference solution updates for Damaris result files.
#
# Uses properties of failing test being processed (${testProperty}).
updateDamarisResults () {
    local binary=${testProperty["binary"]}
    local dir_name=${testProperty["dir_name"]}
    local file_name=${testProperty["file_name"]}
    local test_name=${testProperty["test_name"]}

    copyDamarisToReferenceDir \
        "${BUILD_DIR}/tests/results/${binary}+${test_name}" \
        "${OPM_TESTS_ROOT}/${dir_name}/opm-simulation-reference/${binary}" \
        "${file_name}"

    changed_tests+=" ${test_name}"
}

# Generate reference solution updates for cases that don't run a full
# simulation (e.g., those that use the 'NOSIM' keyword).
#
# Updates the .EGRID and .INIT files only.
updateInitFileResults () {
    local binary=${testProperty["binary"]}
    local dir_name=${testProperty["dir_name"]}
    local file_name=${testProperty["file_name"]}
    local test_name=${testProperty["test_name"]}

    if copyToReferenceDir \
           "${BUILD_DIR}/tests/results/init/${binary}+${test_name}" \
           "${OPM_TESTS_ROOT}/${dir_name}/opm-simulation-reference/${binary}" \
           "${file_name}" \
           EGRID INIT
    then
        changed_tests+=" ${test_name}"
    fi
}

changed_tests=""

# If we get here, the ${failed_test} is one of the compare*Files_*
# tests.  Extract specific information about the test itself and
# update its reference solutions accordingly.

extractTestProperties "$@"

case "$@" in
  *compareECLFiles* )
    updateFullSimulationResults ;;

  *compareECLInitFiles* )
    updateInitFileResults ;;

  *compareDamarisFiles* )
    updateDamarisResults ;;
esac

{
    flock -x 200
    echo $changed_tests >> $TMPDIR/changed_tests
} 200>|$TMPDIR/opm_diff_lock
