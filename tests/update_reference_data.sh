#!/bin/bash

OPM_TESTS_ROOT=$1

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
    diff -q "$WORKSPACE/$SRC_DIR$STEM.$filetype" "$DST_DIR/$STEM.$filetype"
    if test $? -ne 0
    then
      cp "$WORKSPACE/$SRC_DIR$STEM.$filetype" $DST_DIR
      DIFF=0
    fi
  done

  return $DIFF
}

declare -A tests
# binary dirname casename [testname]
# you only have to specify testname if it differs from dirname
tests[spe1]="flow spe1 SPE1CASE1"
tests[spe12]="flow spe1 SPE1CASE2"
tests[spe12p]="flow spe1 SPE1CASE2_2P spe1_2p"
tests[spe1oilgas]="flow spe1 SPE1CASE2_OILGAS spe1_oilgas"
tests[spe1nowells]="flow spe1 SPE1CASE2_NOWELLS spe1_nowells"
tests[spe1thermal]="flow spe1 SPE1CASE2_THERMAL spe1_thermal"
tests[ctaquifer_2d_oilwater]="flow aquifer-oilwater 2D_OW_CTAQUIFER ctaquifer_2d_oilwater"
tests[fetkovich_2d]="flow aquifer-fetkovich 2D_FETKOVICHAQUIFER fetkovich_2d"
tests[msw_2d_h]="flow msw_2d_h 2D_H__"
tests[msw_3d_hfa]="flow msw_3d_hfa 3D_MSW"
tests[polymer_oilwater]="flow polymer_oilwater 2D_OILWATER_POLYMER"
tests[polymer2d]="flow polymer_simple2D 2D_THREEPHASE_POLY_HETER"
tests[spe3]="flow spe3 SPE3CASE1"
tests[spe5]="flow spe5 SPE5CASE1"
tests[spe9group]="flow spe9group SPE9_CP_GROUP"
tests[spe9]="flow spe9 SPE9_CP_SHORT"
tests[wecon_wtest]="flow wecon_wtest 3D_WECON"
tests[spe1_metric_vfp1]="flow vfpprod_spe1 SPE1CASE1_METRIC_VFP1 spe1_metric_vfp1"
tests[base_model_1]="flow model1 BASE_MODEL_1 base_model_1"
tests[msw_model_1]="flow model1 MSW_MODEL_1 msw_model_1"
tests[faults_model_1]="flow model1 FAULTS_MODEL_1 faults_model_1"
tests[polymer_injectivity]="flow polymer_injectivity 2D_POLYMER_INJECTIVITY"
tests[nnc]="flow editnnc NNC_AND_EDITNNC nnc"

changed_tests=""
for test_name in ${!tests[*]}
do
  binary=`echo ${tests[$test_name]} | awk -F ' ' '{print $1}'`
  dirname=`echo ${tests[$test_name]} | awk -F ' ' '{print $2}'`
  casename=`echo ${tests[$test_name]} | awk -F ' ' '{print $3}'`
  tname=`echo ${tests[$test_name]} | awk -F ' ' '{print $4}'`
  test -z "$tname" && tname=$dirname
  copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/$binary+$tname/ \
      $OPM_TESTS_ROOT/$dirname/opm-simulation-reference/$binary \
      $casename \
      EGRID INIT SMSPEC UNRST UNSMRY
  test $? -eq 0 && changed_tests="$changed_tests $test_name"
done

# special tests
copyToReferenceDir \
      $configuration/build-opm-simulators/tests/results/init/flow+norne/ \
      $OPM_TESTS_ROOT/norne/opm-simulation-reference/flow \
      NORNE_ATW2013 \
      EGRID INIT
test $? -eq 0 && changed_tests="$changed_tests norne_init"

changed_tests=`echo $changed_tests | xargs`
echo -e "update reference data for $changed_tests\n" > /tmp/cmsg
if [ -z "$REASON" ]
then
  echo -e "Reason: fill in this\n" >> /tmp/cmsg
else
  echo -e "Reason: $REASON\n" >> /tmp/cmsg
fi
for dep in libecl opm-common opm-grid opm-material ewoms
do
  pushd $WORKSPACE/deps/$dep > /dev/null
  name=`printf "%-14s" $dep`
  rev=`git rev-parse HEAD`
  echo -e "$name = $rev" >> /tmp/cmsg
  popd > /dev/null
done
echo -e "opm-simulators = `git rev-parse HEAD`" >> /tmp/cmsg

cd $OPM_TESTS_ROOT
if [ -n "$BRANCH_NAME" ]
then
  git checkout -b $BRANCH_NAME origin/master
fi

# Add potential new files
untracked=`git status | sed '1,/Untracked files/d' | tail -n +3 | head -n -2`
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
