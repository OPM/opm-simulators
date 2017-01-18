#!/bin/bash

#Stop on first error
set -e

#Echo all commands to terminal
#set -x


# Get the absolute path from a (potentially) relative path
# $1 = relative path
# Prints the absolute path to stdout for capture
getAbsPath() {
  local MY_PATH=$1
  local ORIG_PATH=`pwd`

  MY_PATH="`( cd \"$MY_PATH\" && pwd )`"
  cd "$ORIG_PATH"

  if [ -z "$MY_PATH" ] ; then
    echo "Script path became something strange: '$SCRIPT_PATH'"
    exit 1  # fail
  fi
  echo "$MY_PATH"
}


#Print usage if requested
usage() {
  #Grep through this source file for the options
  echo "$0 usage:" && grep "    .)\ # " $0
  exit 1
}



#Get the script absolute path to reference other scripts
SCRIPT_PATH="`dirname \"$0\"`"
SCRIPT_PATH=`getAbsPath "$SCRIPT_PATH"`



#Get options
OPM_DATA=
WORKSPACE="$SCRIPT_PATH/.."

# Make sure the working tree is clean
pushd . > /dev/null
cd $WORKSPACE
if [[ -n $(git diff-files) ]]
then
  echo "Cannot run with working tree changes. Commit, drop or stash them."
  exit 1
fi
popd > /dev/null


if [ "$OSTYPE" == "linux-gnu" ]
then
  if which nproc > /dev/null
  then
    BUILDTHREADS=`nproc`
  fi
else # OSX
  BUILDTHREADS=`sysctl hw.ncpu | awk -F ' ' '{print $1}'`
fi
test -z "$BUILDTHREADS" && BUILDTHREADS=1

PULL_REQUESTS=
[ $# -eq 0 ] && usage
while getopts "d:w:p:t:o:h" arg; do
  case $arg in
    d) # OPM Data path to use
      echo "OPM_DATA ${OPTARG}"
      OPM_DATA=${OPTARG}
      ;;
    w) # Workspace directory to use for opm-simulators
      echo "WORKSPACE ${OPTARG}"
      WORKSPACE=${OPTARG}
      ;;
    t) # Build threads
      echo "BUILDTHREADS ${OPTARG}"
      BUILDTHREADS=${OPTARG}
      ;;
    p) # Select pull requests to build (e.g., opm-core=1100)
      echo "PULL_REQUESTS ${OPTARG}"
      PULL_REQUESTS=${OPTARG}
      ;;
    o) # Use given toolchain file
      CMAKE_TOOLCHAIN_FILES=`getAbsPath "$OPTARG"`
      echo "Using ${CMAKE_TOOLCHAIN_FILES} toolchain"
      ;;
    h) # Display help.
      usage
      exit 0
      ;;
  esac
done


if [ -z "$OPM_DATA" ]; then
  echo "No OPM_DATA dir, will clone using network (slow)";
fi
WORKSPACE=`getAbsPath "$WORKSPACE"`






#The BUILD_SCRIPT uses environment variables, so let's expose some arguments
export OPM_DATA_ROOT_PREDEFINED="$OPM_DATA"
export OPM_DATA_ROOT="$OPM_DATA"
export WORKSPACE
export BUILDTHREADS
export BTYPES=serial
export CMAKE_TOOLCHAIN_FILES
export ghprbCommentBody="$PULL_REQUESTS"


#Build the requested configuration
$SCRIPT_PATH/../jenkins/build.sh


#Now update reference data-command
echo " "
echo " "
echo " "
echo ">>> To update opm-data, please use the following command <<<"
echo ">>> This creates a commit in $OPM_DATA_ROOT. You can then push the result to a branch in your repository on github, and create a pull request."
echo " "
echo "configuration=serial WORKSPACE=$WORKSPACE $SCRIPT_PATH/update_reference_data.sh $OPM_DATA_ROOT [spe11] [spe12] [spe3] [spe9]"
echo " "

echo ">>> To plot difference between new and old data, usse the following <<<"
echo "SUMMARY_X=<path-to-summary.x> deps/opm-data/norne/plothelper.sh -o compare_runs.pdf -d [deck_name] -r serial/build-opm/simulators/tests/results/<casename> -r $OPM_DATA_ROOT/<casename>/opm-simulation_reference/ -v WBHP -v WOPR -v WGFR -v WWPR -c"
