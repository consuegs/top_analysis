#!/bin/bash -e

source /etc/profile.d/umd_ui.sh #to enable proxy

submitDir=/home/home4/institut_1b/dmeuser/top_analysis/framework/multiprocess
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh

#random sleep to avoid to many simultaneous copy jobs
duration=$[ ( $RANDOM % 15 )  + 1 ]
echo "Sleeping for "$duration" minutes"
sleep $duration"m"

# copy inputs for Unfolding from dCache
syst="$4"
syst=${syst:2}
mkdir $TMP/100.0
cd "$TMP/100.0"
eval `scram unsetenv -sh`;
gfal-copy -r srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="${6}"/minTrees/"$3"/"$5"/minTrees/100.0/Nominal ./Nominal
if [[ $syst != "Nominal" ]]
then
   gfal-copy -r srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="${6}"/minTrees/"$3"/"$5"/minTrees/100.0/$syst ./$syst
fi
cd /home/home4/institut_1b/dmeuser/CMSSW_10_5_0/CMSSW_10_5_0/src
eval `scramv1 runtime -sh`

cd $submitDir

cd ../build

export ANALYSIS_YEAR_CONFIG="$3"
export HOSTNAME=$HOSTNAME

./run.x "$1" "$2" "$4" -m$TMP/

if [[ $8 != "" ]]
then
   rm -vr $TMP/Nominal    #remove input from node 
   rm -vr $TMP/$syst    #remove input from node 
fi
