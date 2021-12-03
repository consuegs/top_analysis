#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess
import sys
import configparser

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

def get_fileNR(dataset,year):    #checks config for number of files for given dataset
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   return len(config[dataset]["files"].split(","))
   
def get_fileName(dataset,year,fileNR):      #checks config for name of file
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   return config[dataset]["files"].split(",")[fileNR]

def get_dataBasePath_dCache(year):      #return dataBasePath on dCache for given year
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   #  ~return "root://xrootd-cms.infn.it///store/user/dmeuser/mergedNtuple/{0}/{1}/".format(year,config["input"]["version"])
   return "root://grid-cms-xrootd.physik.rwth-aachen.de///store/user/dmeuser/mergedNtuple/{0}/{1}/".format(year,config["input"]["version"])

if __name__ == "__main__":
   
   #############################################
   # Select datasets to process
   #############################################
   #  ~toProcess_mc=["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
   #  ~toProcess_mc=["TTbar_diLepton"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_ERDON"]
   #  ~toProcess_mc=["TTbar_diLepton_CR1"]
   #  ~toProcess_mc=["TTbar_diLepton_CR2"]
   #  ~toProcess_mc=["TTbar_diLepton_MTOP169p5"]
   #  ~toProcess_mc=["TTbar_diLepton_MTOP175p5"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_UP","TTbar_diLepton_tau_UETUNE_UP","TTbar_singleLepton_UETUNE_UP","TTbar_hadronic_UETUNE_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_DOWN","TTbar_diLepton_tau_UETUNE_DOWN","TTbar_singleLepton_UETUNE_DOWN","TTbar_hadronic_UETUNE_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_UP","TTbar_diLepton_tau_MATCH_UP","TTbar_singleLepton_MATCH_UP","TTbar_hadronic_MATCH_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_DOWN","TTbar_diLepton_tau_MATCH_DOWN","TTbar_singleLepton_MATCH_DOWN","TTbar_hadronic_MATCH_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_ERDON","TTbar_diLepton_tau_ERDON","TTbar_singleLepton_ERDON","TTbar_hadronic_ERDON"]
   #  ~toProcess_mc=["TTbar_diLepton_CR1","TTbar_diLepton_tau_CR1","TTbar_singleLepton_CR1","TTbar_hadronic_CR1"]
   #  ~toProcess_mc=["TTbar_diLepton_CR2","TTbar_diLepton_tau_CR2","TTbar_singleLepton_CR2","TTbar_hadronic_CR2"]
   #  ~toProcess_mc=["TTbar_diLepton_MTOP169p5","TTbar_diLepton_tau_MTOP169p5","TTbar_singleLepton_MTOP169p5","TTbar_hadronic_MTOP169p5"]
   toProcess_mc=["TTbar_diLepton_MTOP175p5","TTbar_diLepton_tau_MTOP175p5","TTbar_singleLepton_MTOP175p5","TTbar_hadronic_MTOP175p5"]
   
   #  ~toProcess_data=["DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"]
   #  ~toProcess_data=["DoubleMuon","MuonEG","SingleMuon","EGamma"]      #2018
   #  ~toProcess_data=["MET"]      
   toProcess_data=[]
         
   #  ~toProcess_signal=["T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"]
   #  ~toProcess_signal=["T1tttt_1200_800","T1tttt_1500_100","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10"]
   #  ~toProcess_signal=["T1tttt_1200_800"]
   toProcess_signal=[]

   # define command line arguments
   parser = argparse.ArgumentParser()
   parser.add_argument('-m', type=str, default="distributions", help="module, default is distributions")
   parser.add_argument('-f', type=float, choices=[Range(0.0, 1.0)], default=1.0, help="process fraction")
   parser.add_argument('-y', type=str, help="year to be set as ANALYSIS_YEAR_CONFIG",required=True)
   parser.add_argument('-s', type=str, default="Nominal", help="systematic shift")
   parser.add_argument('--scratchInput', action='store_true', default=False, help="Use nTuple stored on scratch, otherwise dCache Input is used.")
   parser.add_argument('--SingleSubmit', action='store_true' )

   args = parser.parse_args()

   print "Running "+args.m

   print "Process "+str(args.f*100)+"% of events"

   print toProcess_mc
   print toProcess_data
   print toProcess_signal
   
   # set dataBasePath depending on dCache options
   dataBasePath = "-d"+get_dataBasePath_dCache(args.y) if args.scratchInput==False else ""
   if args.scratchInput:
      print "Use scratch input"
   else:
       print "Use dCache input"
   
   # Ask if selected settings are correct
   correctSamples = input("If you want to continue with the selected setting, enter 1:\n")
   if (correctSamples != 1):
      print "Abort Submission"
      exit(97)

   # For single submit check if only one dataset is selected and then ask for file nr
   if (args.SingleSubmit and ((len(toProcess_mc)+len(toProcess_data)+len(toProcess_signal))==1)):
      singleFileNR = input("Please input the filenumber to submit:\n")
   elif args.SingleSubmit:
      print "SingleSubmit can only be used if one dataset is selected"
      exit(98)

   # not data processing if systematic shift is choosen
   if (toProcess_data and args.s!="Nominal"):
      toProcess_data=[]
      print "!!!!!!!!!!!!!!!!Data is not processed with systematic shift!!!!!!!!!!!!!!!!!!!!"
   
   # bTagEff only submitted with TTbar_dilepton
   if (args.m == "bTagEff"):
      if (any("TTbar_diLepton" in s for s in toProcess_mc)==False or len(toProcess_mc)>1 or toProcess_data or toProcess_signal):
         print "Error: bTagEff should only be submitted with TTbar_dilepton"
         sys.exit(1)
   
   # triggerEff only submitted with TTbar_dilepton and MET
   if (args.m == "triggerEff"):
      if ("TTbar_diLepton" not in toProcess_mc or len(toProcess_mc)>1 or toProcess_signal or "MET" not in toProcess_data):
         print "Error: triggerEff should only be submitted with TTbar_dilepton and MET"
         sys.exit(1)
   
   # create logpath if not existing
   logpath="logs/"+args.y+"/"+args.s+"/"+str(args.f)+"/"+args.m
   if not os.path.exists(logpath):
      try:
         os.makedirs(logpath)
      except OSError as exc: # Guard against race condition
         if exc.errno != errno.EEXIST:
            raise

   requ_mem=1500   #standard value, allocated if not defined
   
   # loop over selected datasets and submit corresponding jobs
   sampleStr = ""
   startFileNR = 1
   for sel in [[toProcess_mc,"--mc_dataset="],[toProcess_data,"--data_dataset="],[toProcess_signal,"--signal_dataset="]]:
       for x in sel[0]:
               sampleStr=sel[1]+x   # create dataset argument for run.x depending on sample type
               totalFiles=get_fileNR(x,args.y)     # get total number of files for dataset
               
               if (args.SingleSubmit):    # set fileNR if SingleSubmit option is selected
                  if (singleFileNR<=totalFiles):
                     totalFiles=singleFileNR
                     startFileNR=singleFileNR
                  else:
                     print "Selected SingleFileNR is not available"
                     exit(98)
               
               for fileNR in range(startFileNR-1,totalFiles):     # loop over files for given dataset
                                    
                  for logFile in glob.glob(logpath+"/"+args.m+"_"+x+"_"+str(fileNR+1)+"*"):  # Remove old logs
                     if os.path.exists(logFile):
                        os.remove(logFile)
                  
                  submitFile = logpath+"/"+args.m+"_"+x+"_"+str(fileNR+1)+".submit"    #define submit file
                                    
                  with open(submitFile,"w") as f:   # write condor submit
                     f.write("""
Universe                = vanilla
Executable              = run.sh
Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {8}
Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
use_x509userproxy       = true
Request_Memory          = {4} Mb
Requirements            = (TARGET.CpuFamily > 6) && (TARGET.Machine != "lxcip16.physik.rwth-aachen.de")  {9}
Queue
""".format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1),dataBasePath,"\nRank = CpuFamily" if(x=="TTbar_diLepton") else ""),)
                  subprocess.call(["condor_submit", submitFile])
