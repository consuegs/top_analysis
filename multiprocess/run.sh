#!/bin/bash -e
TEMPVAR=`pwd`
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /home/home4/institut_1b/dmeuser/CMSSW_10_5_0/CMSSW_10_5_0/src
eval `scramv1 runtime -sh`

cd $TEMPVAR

cd ../build

# ~./run.x -f0.01 distributions
./run.x "$1" "$2" "$3"
