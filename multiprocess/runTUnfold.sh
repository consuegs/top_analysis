#!/bin/bash -e

source /cvmfs/cms.cern.ch/cmsset_default.sh

submitDir=${9}/multiprocess

#random sleep to avoid to many simultaneous copy jobs
duration=$[ ( $RANDOM % 15 )  + 1 ]
echo "Sleeping for "$duration" minutes"
sleep $duration"m"


# copy inputs for Unfolding from dCache
syst="$4"
syst=${syst:2}
mkdir $TMP/100.0
cd "$TMP/100.0"
eval `scram unsetenv -sh`
gfal-copy -r -t 20000 davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"$7"/minTrees/"$3"/"$5"/minTrees/100.0/Nominal ./Nominal
if [[ $syst != "Nominal" ]]
then
   gfal-copy -r -t 20000 davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"$7"/minTrees/"$3"/"$5"/minTrees/100.0/$syst ./$syst
fi

# ~cmssw-cc7 -- "export SCRAM_ARCH=slc7_amd64_gcc820 && source /cvmfs/cms.cern.ch/cmsset_default.sh && cd ${8} && cmsenv && cd $submitDir && cd ../build && export ANALYSIS_YEAR_CONFIG=$3 && export HOSTNAME=$HOSTNAME && ./run.x $1 $2 $4 -m$TMP/"

if [[ $8 != "" ]]
then
   rm -vr $TMP/100.0/Nominal    #remove input from node
   if [[ $syst != "Nominal" ]]
   then 
      rm -vr $TMP/100.0/$syst    #remove input from node
   fi
fi
