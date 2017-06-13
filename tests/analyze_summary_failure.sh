#!/bin/bash
# This will perform some analysis on a failed summary
# comparisons in a regression tests.

# Analyze summary test failure.
# - Print number of failed keywords
# - Print number of levels each keyword failed at
# - Print maximum deviations for each keyword
analyzeSummaryFailure() {
  lines=`cat $TMPFILE | grep "For keyword"`
  IFS=$'\n'
  kwds=""
  for line in $lines
  do
    kwds+="`echo $line | awk -F ' ' '{print $3}'`\n"
  done
  unique_kwds=`echo -e $kwds | uniq`
  kws_failed=`echo -e "$unique_kwds" | wc -l`
  echo "$kws_failed summary keyword(s) exhibit failures"
  numsteps=`cat $TMPFILE | grep "Comparing " | grep steps | awk -F ' ' '{print $2}'`

  for kwd in $unique_kwds
  do
    lines=`cat $TMPFILE | grep -n "For keyword $kwd$"`
    echo -e "\t $kwd"
    echo -e "\t \t Fails for: `cat $TMPFILE | grep -n "$kwd$" | wc -l` / $numsteps steps."
    abs=0
    rel=0
    for line in $lines
    do
      ln=`echo $line | awk -F ':' '{print $1'}`
      abs_line=$(($ln+2))
      rel_line=$(($ln+3))
      abs_new=`sed "${abs_line}q;d" $TMPFILE | awk -F ' ' '{print $5}'`
      rel_new=`sed "${rel_line}q;d" $TMPFILE | awk -F ' ' '{print $5}'`
      abs_new=`echo ${abs_new} | sed -e 's/\.$//g' | awk '{printf sprintf("%.16f", $1); }'`
      rel_new=`echo ${rel_new} | sed -e 's/\.$//g' | awk '{printf sprintf("%.16f", $1); }'`
      if [ `bc <<< "$abs_new>$abs"` -eq 1 ]
      then
        abs=$abs_new
      fi
      if [ `bc <<< "$rel_new>$rel"` -eq 1 ]
      then
        rel=$rel_new
      fi
    done
    echo -e "\t\t Largest absolute error: `echo $abs | awk '{printf sprintf("%e", $1); }'`"
    echo -e "\t\t Largest relative error: `echo $rel | awk '{printf sprintf("%e", $1); }'`"
  done
}

COMPARE_SUMMARY_COMMAND=$1
PARAM=$2
FILE1=$3
FILE2=$4
ABS_TOL=$5
REL_TOL=$6

TMPFILE=$(mktemp)
${COMPARE_SUMMARY_COMMAND} -p -n ${PARAM} ${FILE1} ${FILE2} ${ABS_TOL} ${REL_TOL} &> $TMPFILE

analyzeSummaryFailure
rm -f $TMPFILE
