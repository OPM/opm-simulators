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

    tests=""
    while read failed_test
    do
        if ! grep -q 'compare.*Files_.*' <<< ${failed_test}
        then
            # The ${failed_test} is not among the compare*Files_* tests to
            # which this update script applies.  Nothing to do.
            continue
        fi

        tests+="${failed_test}\\n"
    done < "${logfile}"

    JOBS=${TESTTHREADS:-16}
    echo -e $tests | \
      OPM_TESTS_ROOT=$OPM_TESTS_ROOT \
      BUILD_DIR=$BUILD_DIR \
      CONVERT_ECL=$CONVERT_ECL \
      TMPDIR=$TMPDIR \
        xargs -P${JOBS} -L1 $(dirname $0)/update_test_reference.sh
    changed_tests=$(cat $TMPDIR/changed_tests)
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
MSGFILE=$WORKSPACE/deps/cmsg
echo -e "Automatic Reference Data Update for ${REASON:-(Unknown)}\n" > $MSGFILE
if [ -z "$REASON" ]
then
  echo -e "Reason: fill in this\n" >> $MSGFILE
else
  echo -e "Reason: $REASON\n" >> $MSGFILE
fi
if [ -n "$CONVERT_ECL" ]
then
  for dep in opm-common opm-grid opm-simulators
  do
    pushd $WORKSPACE/deps/$dep > /dev/null
    name=`printf "%-14s" $dep`
    rev=`git rev-parse HEAD`
    echo -e "$name = $rev" >> $MSGFILE
    popd > /dev/null
  done
fi

echo -e "\n### Changed Tests ###\n" >> $MSGFILE
printf "  * %s\n" ${changed_tests} >> $MSGFILE

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
  git commit -a -t $MSGFILE
else
  git commit -a -F $MSGFILE
fi

# ===========================================================================
# Clean up intermediate files.

rm -rf $TMPDIR
