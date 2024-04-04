## Requirements ##
- ROOT6
- C++17 compiler
- CMake >= 2.8 (lx-cluster: cmake28 executable)
- Boost
- Scientific Linux 7 (use singularity)

## Setup scientific linux 7 singularity ##
Before setting up the singularity, the repository has to be cloned. Add the following lines to your `.bashrc` to start the singularity each time you run on the Alma9 node (standard node on the physik cluster). To work properly the repository has to be clone into your home directory.

    # Run singularity if not on cc7 node
    linuxVersion=$( cat /etc/redhat-release )
    hostname=$( hostname )
    check="portal"
    
    if echo "$hostname" | grep -q "$check"; then
      export RUN_ON_NODE=1;
    else
      export RUN_ON_NODE=0;
    fi
    
    if [[ $linuxVersion == "AlmaLinux release 9.3 (Shamrock Pampas Cat)" ]]
    then
        if [[ $RUN_ON_NODE == 0 ]]
        then
    	echo "Start singularity"
    	export RUN_IN_SINGULARITY=1
    	cd ~
    	./top_analysis/tools/cc7_singularity/cmssw-cc7_modified
        fi
    fi

To exit the singularity, just use `exit`.

## Setup CMSSW Environmet ##
Before setting up the CMSSW environment, run `source /cvmfs/cms.cern.ch/cmsset_default.sh` or, even better, add the command to your `.bashrc`.

    export SCRAM_ARCH=slc7_amd64_gcc820
    cmsrel CMSSW_10_5_0
    cd CMSSW_10_5_0/src/
    cmsenv
    # Following lines enable use of CombineHarvester in framework
    git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
    cd $CMSSW_BASE/src/HiggsAnalysis/CombinedLimit
    git fetch origin
    git checkout v8.2.0
    cd $CMSSW_BASE/src/
    git clone --branch 102x https://github.com/cms-analysis/CombineHarvester.git CombineHarvester
    scram b -j7

The CMSSW version has to be sourced before running the framework. Would recommend to include an alias.

## Build framework ##
    cd ~/top_analysis
    mkdir build; cd build
    cmake ..
    make

## Set environment variable for data taking period ##
The framework is steered based on the entry given in `$ANALYSIS_YEAR_CONFIG`. For convenience just add the following aliases to your `.bashrc`

    alias set2016preVFP="export ANALYSIS_YEAR_CONFIG=2016_preVFP"
    alias set2016postVFP="export ANALYSIS_YEAR_CONFIG=2016_postVFP"
    alias set2017="export ANALYSIS_YEAR_CONFIG=2017"
    alias set2018="export ANALYSIS_YEAR_CONFIG=2018"
    

## Run ##
Execute `run.x` (created in the build directory).
Meaningful modules can be chosen as arguments:

    $ ./run.x --help
        Running on year 2018 (to change set ANALYSIS_YEAR_CONFIG variable)
        Usage: ./run.x module1[ module2[...]] [options]
        Allowed options:
          -h [ --help ]                         produce help message
          -f [ --fraction ] arg (=1)            Fraction of events to process
          --mc_dataset arg                      Single MC dataset to be processed, if 
                                                empty datasets are taken from 
                                                config.ini
          --data_dataset arg                    Single data dataset to be processed, if
                                                empty datasets are taken from 
                                                config.ini
          --signal_dataset arg                  Single signal dataset to be processed, 
                                                if empty datasets are taken from 
                                                config.ini
          --fileNR arg (=0)                     FileNR to be processed, if zero all 
                                                files are processed (filelist defined 
                                                in .ini)
          -s [ --systematic ] arg (=Nominal)    Systematic uncertainty to be applied 
                                                (only for distributions)
          -d [ --dataBasePath ] arg (=/net/data_cms1b/user/dmeuser/top_analysis/2018/v08/nTuple/)
                                                Alternative dataBasePath, which can be 
                                                used for condor submission with 
                                                inputFileTransfer
          -m [ --minTreePath ] arg (=/net/data_cms1b/user/dmeuser/top_analysis/2018/v08/minTrees/)
                                                Alternative minTreePath, which can be 
                                                used for condor submission with 
                                                inputFileTransfer
          --release                             Release mode (don't draw version 
                                                labels)


## Adding new modules ##
To add a new module, create a new `<module>.cpp` file in `src/modules`.
The entry function has to be `extern "C" void run() {...}`.
Add `<module>` to the list in `src/modules/CMakeLists.txt`

## Running on condor: general information ##
The `distributions`, `bTagEff` and `triggerEff` modules can be executed in multiple condor jobs. One job corresponds to one sample file (defined in `config*.ini`). The scripts for submitting the jobs can be found in `multiprocess/`. 
Attention: The output hists are saved in an additional folder `output_framework/multiHists/` (in case of the distributions module) and have to be merged before plotting. In addition, the `TUnfold_binning` module can be submitted for individual systematics. For a new user to use the submission option a corresponding user config in `users/` has to be created and a grid certificate has to be present to setup a voms proxy.

## Running on condor: site dependencies ##
Currently there are three different options to submit condor jobs using `multiprocess/process.py`. Two options can be used from the lxcluster nodes:
* **Submitting to `lxbatch*.physik.rwth-aachen.de` nodes:** Jobs are submitted to local cluster nodes, with direct network connection to the dCache. This option is selcted based on `--copyDCache`
* **Direct Grid submission**: Jobs are directly submitted to `grid-ce-1(2)-rwth.gridka.de`. This submission is currently only possible from the remaining cc7 cluster nodes, which can be accessed via `ssh dmeuser@lxv7.physik.rwth-aachen.de -i <public key>`. This option is expected to be no longer present once the CEs migrate to a newer HTcondor version (around summer 2024).

The third option is based on the [cmsconnect service](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCMSConnect). To use this service, follow the instructions given in the twiki. Based on this service the jobs can be submitted to the CMS global pool using via a specifig login node. Since there is no shared filesystem with the lxcluster nodes, the framework has to be setup on the login node. In addition, the storage of the outputs has to be properly solved (currently using `/scratch(<username>` seems to be an option). Furthermore, a new user config `<username>_cmsconnect.ini` has to be created. Once the framework is setup, the jobs can be submitted using the additional `--cmsconnect` option.
...
