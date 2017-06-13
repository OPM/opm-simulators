#!/bin/bash
# This will perform some analysis on a failed restart/init
# comparisons in a regression tests.

# Analyze restart/init test failure.
# - Print failed keywords
# - Print maximum deviations for each keyword
analyzeEclFailure() {
kwds=`cat $KWTMPFILE | tr '\n' ' ' | sed -e 's/.*Common keywords for the two cases:\(.*\)Uncommon.*/\1/'`
TMPFILE=$(mktemp)
allok=1
for kwd in $kwds
do
  ${COMPARE_ECL_COMMAND} -n -l -k ${kwd} ${FILE1} ${FILE2} ${ABS_TOL} ${REL_TOL} &> ${TMPFILE}
  nfailure=`cat ${TMPFILE} | grep "Deviations exceed tolerances" | wc -l`
  if [ ${nfailure} -ne 0 ]
  then
    allok=0
    echo "Failure for keyword ${kwd}"
    echo -e "\t Fails for ${nfailure} entries"
    alines=`cat $TMPFILE | grep "The absolute deviation is"`
    IFS=$'\n'
    abs=0
    for line in $alines
    do
      abs_new=`echo $line| awk -F ' ' '{print $5}'`
      abs_new=`echo ${abs_new} | sed -e 's/\.$//g' | awk '{printf sprintf("%.16f", $1); }'`
      if [ `bc <<< "$abs_new>$abs"` -eq 1 ]
      then
        abs=$abs_new
      fi
    done
    rlines=`cat $TMPFILE | grep "The relative deviation is"`
    rel=0
    for line in $rlines
    do
      rel_new=`echo $line| awk -F ' ' '{print $5}'`
      rel_new=`echo ${rel_new} | sed -e 's/\.$//g' | awk '{printf sprintf("%.16f", $1); }'`
      if [ `bc <<< "$rel_new>$rel"` -eq 1 ]
      then
        rel=$rel_new
      fi
    done
    echo -e "\t Largest absolute deviation: `echo $abs | awk '{printf sprintf("%e", $1); }'`"
    echo -e "\t Largest relative deviation: `echo $rel | awk '{printf sprintf("%e", $1); }'`"
  fi
done
if [ $allok -eq 1 ]
then
  echo "Comparisons pass for all common keywords."
fi
}

COMPARE_ECL_COMMAND=$1
TYPE=$2
FILE1=$3
FILE2=$4
ABS_TOL=$5
REL_TOL=$6

KWTMPFILE=$(mktemp)
${COMPARE_ECL_COMMAND} -t ${TYPE} -l -P ${FILE1} ${FILE2} ${ABS_TOL} ${REL_TOL} &> $KWTMPFILE

cat $KWTMPFILE
echo ""

analyzeEclFailure

rm -f $KWTMPFILE
