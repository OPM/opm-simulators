#!/bin/bash -norc

# Generate OPM-Tests repository PR to update one or more reference
# solution files as a result of a known improvement to the simulator
# or a new regression test being added into the test suite.
#
# Processes 'LastTestsFailed.log' in OPM-Simulators' build directory
# to infer applicable failing tests.
#
# Positional arguments:
#   - $1 Location of OPM-Tests root directory.
#   - $2 Location of OPM-Simulators' build directory.
#   - $3 Full path to 'convertECL' utility, including the 'convertECL'
#        utility name. Needed only if caller wants 'diff -u' style
#        content differences between existing reference solutions and
#        new candidate solutions. Pass empty string if not.
#
# Environment:
#   - $WORKSPACE Location of Jenkins' build work space.
#   - $configuration Name of Jenkins build configuration.
#   - $REASON Underlying reason for this data update. Typically one or
#     more PR names. An empty REASON will stop commit process and
#     prompt user for an appropriate commit message.
#   - $BRANCH_BASE Branch relative to which to start this data
#     update. Typically 'master'.
#   - $BRANCH_NAME Name of data update PR branch in OPM-Tests repo.

OPM_TESTS_ROOT=$1
BUILD_DIR=$2
CONVERT_ECL=$3

TMPDIR=$(mktemp -d)
mkdir -p "${TMPDIR}"/{orig,new}

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

    if [ -s "${reference_dir}/${result_file}" ]
    then
        cp -a "${reference_dir}/${result_file}" "${TMPDIR}/orig"
        ref=$(${CONVERT_ECL} "${TMPDIR}/orig/${result_file}" | \
                  awk '/converting/{print $NF}')
    else
        # Reference file does not exist. New file type or new case.
        ref=/dev/null
    fi

    if [ -s "${simoutput_dir}/${result_file}" ]
    then
        cp -a "${simoutput_dir}/${result_file}" "${TMPDIR}/new"
        new=$(${CONVERT_ECL} "${TMPDIR}/new/${result_file}" | \
                  awk '/converting/{print $NF}')
    else
        # New simulation file does not exist. Typically when the file
        # type is removed. Unexpected.
        new=/dev/null
    fi

    if [ "${new}" != "${ref}" ]
    then
        diff -u "${ref}" "${new}" >> "${WORKSPACE}/data_diff"
    else
        # Neither reference nor updated file exists. This really
        # shouldn't happen.
        echo "No difference between reference and update for ${result_file}" \
             >> "${WORKSPACE}/data_diff"
    fi
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

    testProperty["binary"]=$(awk -v search="set_tests_properties\\\(${failed_test}\$" \
                             -v prop="SIMULATOR" -f "${dir}/getprop.awk" \
                             "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["dir_name"]=$(awk -v search="set_tests_properties\\\(${failed_test}\$" \
                               -v prop="DIRNAME" -f "${dir}/getprop.awk" \
                               "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["file_name"]=$(awk -v search="set_tests_properties\\\(${failed_test}\$" \
                                -v prop="FILENAME" -f "${dir}/getprop.awk" \
                                "${BUILD_DIR}/CTestTestfile.cmake")

    testProperty["test_name"]=$(awk -v search="set_tests_properties\\\(${failed_test}\$" \
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
    h5diff -v "$file" "$DST_DIR/`basename $file`" >> $WORKSPACE/data_diff
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

    if [ -d "${configuration}/build-opm-simulators/tests/results/${binary}+${test_name}/restart" ]
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

# ===========================================================================
# Main data update loop.
#
# Processes one line of 'LastTestsFailed.log' at a time and updates
# applicable reference solutions.

for logfile in $(echo "${BUILD_DIR}/Testing/Temporary/LastTestsFailed*.log")
do
    if [ ! -s "${logfile}" ]
    then
        continue
    fi

    while read failed_test
    do
        if ! grep -q 'compare.*Files_.*' <<< ${failed_test}
        then
            # The ${failed_test} is not among the compare*Files_* tests to
            # which this update script applies.  Nothing to do.
            continue
        fi

        # If we get here, the ${failed_test} is one of the compare*Files_*
        # tests.  Extract specific information about the test itself and
        # update its reference solutions accordingly.

        extractTestProperties "${failed_test}"

        case "${failed_test}" in
            *compareECLFiles* )
                updateFullSimulationResults ;;

            *compareECLInitFiles* )
                updateInitFileResults ;;

            *compareDamarisFiles* )
                updateDamarisResults ;;
        esac
    done < "${logfile}"
done

if [ -z "${changed_tests}" ]
then
    exit 5
fi

# ===========================================================================
# Create data update PR.

# 1) Create commit message (or commit message template, depending on
# the contents of ${REASON}).  Empty ${REASON} ultimately starts an
# interactive commit.
changed_tests=`echo $changed_tests | xargs`
echo -e "Automatic Reference Data Update for ${REASON:-(Unknown)}\n" > /tmp/cmsg
if [ -z "$REASON" ]
then
  echo -e "Reason: fill in this\n" >> /tmp/cmsg
else
  echo -e "Reason: $REASON\n" >> /tmp/cmsg
fi
if [ -n "$CONVERT_ECL" ]
then
  for dep in opm-common opm-grid
  do
    pushd $WORKSPACE/deps/$dep > /dev/null
    name=`printf "%-14s" $dep`
    rev=`git rev-parse HEAD`
    echo -e "$name = $rev" >> /tmp/cmsg
    popd > /dev/null
  done
  echo -e "opm-simulators = `git rev-parse HEAD`" >> /tmp/cmsg
fi

echo -e "\n### Changed Tests ###\n" >> /tmp/cmsg
printf "  * %s\n" ${changed_tests} >> /tmp/cmsg

# ---------------------------------------------------------------------------

# 2) Create branch for new commit.
cd $OPM_TESTS_ROOT
if [ -n "$BRANCH_NAME" ]
then
  git checkout -b $BRANCH_NAME $BRANCH_BASE
fi

# ---------------------------------------------------------------------------

# 3) Add new files as needed.
untracked=`git status --porcelain | awk '$1~/\?/{print $2}'`
if [ -n "$untracked" ]
then
  git add $untracked
fi

# ---------------------------------------------------------------------------

# 4) Commit reference solution update.
if [ -z "$REASON" ]
then
  git commit -a -t /tmp/cmsg
else
  git commit -a -F /tmp/cmsg
fi

# ===========================================================================
# Clean up intermediate files.

rm -rf $TMPDIR
