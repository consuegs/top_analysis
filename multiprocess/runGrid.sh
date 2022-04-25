#!/bin/bash

echo $(date)
sourceDir=`pwd`

echo $(date)
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo $(date)
dccp dcap://grid-dcap.physik.rwth-aachen.de"$7"/gridJobInputs/CMSSW_10_5_0.tgz ./
tar -xf CMSSW_10_5_0.tgz
cd CMSSW_10_5_0/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

cd $sourceDir

echo $(date)
dccp dcap://grid-dcap.physik.rwth-aachen.de"$7"/gridJobInputs/FW.tgz ./
tar -xf FW.tgz
mkdir build
cd build
cmake ..
make

echo $(date)

#random sleep to avoid to many simultaneous copy jobs
duration=$[ ( $RANDOM % 5 )  + 1 ]
echo "Sleeping for "$duration" minutes"
sleep $duration"m"

dccp "$8" ../

echo $(date)

export ANALYSIS_YEAR_CONFIG="$4"
export HOSTNAME=$HOSTNAME

./run.x "$1" "$2" "$3" "$5" "$6" -d ../
