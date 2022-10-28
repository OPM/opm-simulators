#!/bin/bash

OPM_TESTS_ROOT=$1
BUILD_DIR=$2
CONVERT_ECL=$3

TMPDIR=`mktemp -d`
mkdir $TMPDIR/orig
mkdir $TMPDIR/new

# Copy results from a test run to refence dir
# $1 = source directory to copy data from
# $2 = destination directory to copy data to
# $3 = base file name for files to copy
# $4...$@ = file types to copy
copyToReferenceDir () {
  SRC_DIR=$1
  DST_DIR=$2
  STEM=$3
  FILETYPES=${@:4}
  mkdir -p $DST_DIR

  DIFF=1
  for filetype in $FILETYPES
  do
    # Don't flag as changed if both reference and result dir lack a file type
    # In particular to handle the optional RFT's
    if [ ! -f $SRC_DIR/$STEM.$filetype ] && [ ! -f $DST_DIR/$STEM.$filetype ]
    then
      continue
    fi
    diff -q "$SRC_DIR/$STEM.$filetype" "$DST_DIR/$STEM.$filetype"
    res=$?
    if test $res -ne 0 && test -n "$CONVERT_ECL"
    then
      cp $SRC_DIR/$STEM.$filetype $TMPDIR/new
      $CONVERT_ECL $TMPDIR/new/$STEM.$filetype
      cp $DST_DIR/$STEM.$filetype $TMPDIR/orig
      $CONVERT_ECL $TMPDIR/orig/$STEM.$filetype
      diff -u $TMPDIR/orig/$STEM.F$filetype $TMPDIR/new/$STEM.F$filetype >> $WORKSPACE/data_diff
    fi
    if test $res -ne 0
    then
      cp "$SRC_DIR/$STEM.$filetype" $DST_DIR
      DIFF=0
    fi
  done

  return $DIFF
}

declare -A tests
# Read test definitions from json files

DIRECTORY=`dirname $0`

input=`python3 $DIRECTORY/gen_test_data.py $DIRECTORY/definitions/regression`

while IFS= read -r line; do
  read -r line2
  tests[$line]="$line2"
done <<< "$input"

changed_tests=""

# Read failed tests
FAILED_TESTS=`cat $BUILD_DIR/Testing/Temporary/LastTestsFailed*.log`

test -z "$FAILED_TESTS" && exit 5

for failed_test in $FAILED_TESTS
do
  failed=`echo $failed_test | sed -e 's/.*:compareECLFiles_//g'`
  for test_name in ${!tests[*]}
  do
    binary=`echo ${tests[$test_name]} | awk -F ' ' '{print $1}'`
    dirname=`echo ${tests[$test_name]} | awk -F ' ' '{print $2}'`
    casename=`echo ${tests[$test_name]} | awk -F ' ' '{print $3}'`
    if grep -q "$failed" <<< "$binary+$casename"
    then
      copyToReferenceDir \
          $BUILD_DIR/tests/results/$binary+$test_name \
          $OPM_TESTS_ROOT/$dirname/opm-simulation-reference/$binary \
          $casename \
          EGRID INIT RFT SMSPEC UNRST UNSMRY
      test $? -eq 0 && changed_tests="$changed_tests $test_name"

      if [ -d $configuration/build-opm-simulators/tests/results/$binary+$test_name/restart ]
      then

        RSTEPS=`ls -1 $BUILD_DIR/tests/results/$binary+$test_name/restart/*.UNRST | sed -e 's/.*RESTART_*//' | sed 's/[.].*//' `
        result=0
        for RSTEP in $RSTEPS
        do
          copyToReferenceDir \
              $BUILD_DIR/tests/results/$binary+$test_name/restart/ \
              $OPM_TESTS_ROOT/$dirname/opm-simulation-reference/$binary/restart \
              ${casename}_RESTART_${RSTEP} \
              EGRID INIT RFT SMSPEC UNRST UNSMRY
          res=$?
          test $result -eq 0 || result=$res
        done
        test $result -eq 0 && changed_tests="$changed_tests $test_name(restart)"
      fi
    fi
  done
done

# special tests
copyToReferenceDir \
      $BUILD_DIR/tests/results/init/flow+norne \
      $OPM_TESTS_ROOT/norne/opm-simulation-reference/flow \
      NORNE_ATW2013 \
      EGRID INIT
test $? -eq 0 && changed_tests="$changed_tests norne_init"

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
  for dep in opm-common opm-grid opm-material opm-models
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
for t in ${changed_tests}
do
  echo "  * ${t}" >> /tmp/cmsg
done

cd $OPM_TESTS_ROOT
if [ -n "$BRANCH_NAME" ]
then
  git checkout -b $BRANCH_NAME $BRANCH_BASE
fi

# Add potential new files
untracked=`git status --porcelain | awk '$1~/\?/{print $2}'`
if [ -n "$untracked" ]
then
  git add $untracked
fi

if [ -z "$REASON" ]
then
  git commit -a -t /tmp/cmsg
else
  git commit -a -F /tmp/cmsg
fi

rm -rf $TMPDIR
