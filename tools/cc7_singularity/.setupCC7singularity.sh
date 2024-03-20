#!/bin/bash
TEMPVAR=`pwd`
source /cvmfs/cms.cern.ch/cmsset_default.sh

LC_ALL="en_US.UTF-8"
export LC_ALL

export PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external/cmake/3.10.2/bin:$PATH"
export PATH="/cvmfs/cms.cern.ch/slc7_amd64_gcc820/external:$PATH"
export PATH="//usr/lib64/python2.7/site-packages:$PATH"

cd $TEMPVAR
