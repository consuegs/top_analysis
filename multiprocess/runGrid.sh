#!/bin/bash -e

echo $(date)
sourceDir=`pwd`

echo $(date)
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh

#get CMSSW from dCache and setup environment
echo $(date)
dccp dcap://grid-dcap.physik.rwth-aachen.de"$7"/gridJobInputs/CMSSW_10_5_0.tgz ./
tar -xf CMSSW_10_5_0.tgz
cd CMSSW_10_5_0/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

cd $sourceDir

#get FW from dCache and compile
echo $(date)
dccp dcap://grid-dcap.physik.rwth-aachen.de"$7"/gridJobInputs/FW.tgz ./
tar -xf FW.tgz
mkdir build
cd build
cmake ..
make

echo $(date)

#random sleep to avoid to many simultaneous copy jobs
duration=$[ ( $RANDOM % 15 )  + 1 ]
echo "Sleeping for "$duration" minutes"
sleep $duration"m"

#get inputs from dCache
if [[ $2 != "TUnfold_binning" ]]
then
   dccp "$8" ../
else
   syst="$5"
   syst=${syst:2}
   mkdir ../100.0
   cd "../100.0"
   eval `scram unsetenv -sh`;
   gfal-copy -r -t 7200 srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="${10}"/minTrees/"$4"/"$9"/minTrees/100.0/Nominal ./Nominal
   if [[ $syst != "Nominal" ]]
   then
      gfal-copy -r -t 7200 srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="${10}"/minTrees/"$4"/"$9"/minTrees/100.0/$syst ./$syst
   fi
   cd ../CMSSW_10_5_0/src/
   eval `scramv1 runtime -sh`
   cd $sourceDir
   cd build
fi

echo $(date)

export ANALYSIS_YEAR_CONFIG="$4"
export HOSTNAME=$HOSTNAME
echo $HOSTNAME

#run framework
if [[ $2 != "TUnfold_binning" ]]
then
   ./run.x "$1" "$2" "$3" "$5" "$6" -d ../
else
   ./run.x "$1" "$2" "$5" -m ../
fi

#check if output directory exist
FILE=../output_framework
if [ -d "$FILE" ]; then
   echo "$FILE exists."
else
   echo "$FILE does not exists, job fails."
   exit 1
fi


#copy mintree to dCache if distribution module is running
if [[ $2 == "distributions" ]]
then
   echo "Copy minTree output to dCache"
   eval `scram unsetenv -sh`; gfal-copy -r -f ../minTrees/ srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="$7"/minTrees/"$4"/"$9"/minTrees/
   rm -r ../minTrees/
fi

exit 0

