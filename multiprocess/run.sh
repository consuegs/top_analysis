#!/bin/bash -e

source /etc/profile.d/umd_ui.sh #to enable proxy

export HOSTNAME=$HOSTNAME

#random sleep to avoid to many simultaneous copy jobs
duration=$[ ( $RANDOM % 10 )  + 1 ]
echo "Sleeping for "$duration" minutes"
sleep $duration"m"

if [[ $8 != "" ]]
then
   # ~xrdcp "$8" $TMP   #copy to node
   dccp "$8" $TMP
fi

submitDir=${10}/multiprocess # framework base dir
export SCRAM_ARCH=slc7_amd64_gcc820
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${9} # CMSSW/src dir
eval `scramv1 runtime -sh`

cd $submitDir

cd ../build

export ANALYSIS_YEAR_CONFIG="$4"

./run.x "$1" "$2" "$3" "$5" "$6" "$7" 2>&1 | tee $TMP/OutFile.txt

#copy minTrees to dCache
if [[ $1 != "distributions" ]]
then
   minTreePath=$(cat $TMP/OutFile.txt | grep minTree | tr "'" "\n"| sed -n '4p')
   echo $minTreePath
   minTreePathHelper="${minTreePath/output_framework/}"
   minTreePathHelper="${minTreePathHelper/../}"
   minTreePathHelper="${minTreePathHelper/top_analysis/+minTrees}"
   minTreePathDcache=$(echo $minTreePathHelper | tr "+" "\n"| sed -n '2p')
   echo $minTreePathDcache
   echo "Copy minTree output to dCache"
   eval `scram unsetenv -sh`; gfal-copy -r -f $minTreePath srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN="${11}"/$minTreePathDcache
fi

if [[ $8 != "" ]]
then
   rm -v $TMP/*.root    #remove input from node 
   rm -v $TMP/*.txt    #remove input from node 
fi
