#!/bin/bash -e
# ~inputFileDir=`pwd`
# ~export PATH=$PATH:/net/software_t2k/tools/cccp/
# ~cccp "$7" ./

submitDir=/home/home4/institut_1b/dmeuser/top_analysis/framework/multiprocess
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/home4/institut_1b/dmeuser/CMSSW_10_5_0/CMSSW_10_5_0/src
eval `scramv1 runtime -sh`

cd $submitDir

cd ../build

export ANALYSIS_YEAR_CONFIG="$4"

# ~./run.x "$1" "$2" "$3" "$5" "$6" -d "$inputFileDir/"
./run.x "$1" "$2" "$3" "$5" "$6" 

# ~cd $inputFileDir
# ~rm ./*.root    #Remove copied files
