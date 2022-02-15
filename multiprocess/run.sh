#!/bin/bash -e

source /etc/profile.d/umd_ui.sh #to enable proxy

if [[ $8 != "" ]]
then
   # ~xrdcp "$8" $TMP   #copy to node
   dccp -H "$8" $TMP
fi

submitDir=${10}/multiprocess # framework base dir
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${9} # CMSSW/src dir
eval `scramv1 runtime -sh`

cd $submitDir

cd ../build

export ANALYSIS_YEAR_CONFIG="$4"

./run.x "$1" "$2" "$3" "$5" "$6" "$7"

if [[ $8 != "" ]]
then
   rm -v $TMP/*.root    #remove input from node 
fi
