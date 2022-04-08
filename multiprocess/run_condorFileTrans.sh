#!/bin/bash -e
# this bash script requires the input file to be transfered with condor file transfer (if not use run.sh)

inputFileDir=`pwd`

source /etc/profile.d/umd_ui.sh #to enable proxy

submitDir=${8}/multiprocess # framework base dir
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${7} # CMSSW/src dir
eval `scramv1 runtime -sh`

cd $submitDir

cd ../build

export ANALYSIS_YEAR_CONFIG="$4"

./run.x "$1" "$2" "$3" "$5" "$6" -d $inputFileDir/
