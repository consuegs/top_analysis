#!/bin/bash -e

echo $(date)
sourceDir=`pwd`

echo $(date)
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
source /cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/cmake/3.10.2/etc/profile.d/init.sh
export HOSTNAME=$HOSTNAME
echo $HOSTNAME

#Get CMSSW and FW from dCache
echo $(date): "Get CMSSW and FW from dCache"
gfal-copy davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"${10}"/gridJobInputs/CMSSW_10_5_0.tgz ./
gfal-copy davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"${10}"/gridJobInputs/FW.tgz ./

#setup CMSSW environment
echo $(date)
tar -xf CMSSW_10_5_0.tgz
cd CMSSW_10_5_0/src/
scramv1 b ProjectRename
eval `scramv1 runtime -sh`

cd $sourceDir

#compile FW
echo $(date)
tar -xf FW.tgz
mkdir build
cd build
# ~cmake ..
cmake .. -DPYTHON_INCLUDE_DIR=$(python -c "import sysconfig; print(sysconfig.get_path('include'))") -DPYTHON_LIBRARY=$(python -c "import sysconfig; print(sysconfig.get_config_var('LIBDIR'))")
make

echo $(date)

cd $sourceDir

#random sleep to avoid to many simultaneous copy jobs
# ~duration=$[ ( $RANDOM % 15 )  + 1 ]
# ~echo "Sleeping for "$duration" minutes"
# ~sleep $duration"m"

#get inputs from dCache
if [[ $2 != "TUnfold_binning" ]]
then
   eval `scram unsetenv -sh`;
   gfal-copy -r -t 3600 "$8" ./
   cd CMSSW_10_5_0/src/
   eval `scramv1 runtime -sh`
   cd $sourceDir
   cd build
else
   syst="$5"
   syst=${syst:2}
   mkdir 100.0
   cd "100.0"
   eval `scram unsetenv -sh`;
   echo grid-webdav.physik.rwth-aachen.de:2889/store/user/"${10}"/minTrees/"$4"/"$9"/minTrees/100.0/Nominal
   gfal-copy -r -t 3600 davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"${10}"/minTrees/"$4"/"$9"/minTrees/100.0/Nominal ./Nominal
   if [[ $syst != "Nominal" ]]
   then
      gfal-copy -r -t 3600 davs://grid-webdav.physik.rwth-aachen.de:2889/store/user/"${10}"/minTrees/"$4"/"$9"/minTrees/100.0/$syst ./$syst
   fi
   cd ../CMSSW_10_5_0/src/
   eval `scramv1 runtime -sh`
   cd $sourceDir
   cd build
fi

echo $(date)

export ANALYSIS_YEAR_CONFIG="$4"

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

#remove TUnfold_binning inputs
if [[ $2 == "TUnfold_binning" ]]
then
   echo "Removing TUnfold inputs"
   rm -r ../100.0/
fi

exit 0

