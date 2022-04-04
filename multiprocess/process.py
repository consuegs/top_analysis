#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess
import sys
import configparser
import sys
sys.path.append("../users")
from getPath import getPath

bTagEff_sample_syst_dict = {
      "Nominal" : "TTbar_diLepton",
      "JESTotal_UP" : "TTbar_diLepton",
      "JESTotal_DOWN" : "TTbar_diLepton",
      "JER_UP" : "TTbar_diLepton",
      "JER_DOWN" : "TTbar_diLepton",
      "UETUNE_UP" : "TTbar_diLepton_UETUNE_UP",
      "UETUNE_DOWN" : "TTbar_diLepton_UETUNE_DOWN",
      "MATCH_UP" : "TTbar_diLepton_MATCH_UP",
      "MATCH_DOWN" : "TTbar_diLepton_MATCH_DOWN",
      "ERDON" : "TTbar_diLepton_ERDON",
      "CR1" : "TTbar_diLepton_CR1",
      "CR2" : "TTbar_diLepton_CR2",
      "MTOP169p5" : "TTbar_diLepton_MTOP169p5",
      "MTOP175p5" : "TTbar_diLepton_MTOP175p5",
      "TOP_PT" : "TTbar_diLepton",
      "MESCALE_UP" : "TTbar_diLepton",
      "MESCALE_DOWN" : "TTbar_diLepton",
      "MERENSCALE_UP" : "TTbar_diLepton",
      "MERENSCALE_DOWN" : "TTbar_diLepton",
      "MEFACSCALE_UP" : "TTbar_diLepton",
      "MEFACSCALE_DOWN" : "TTbar_diLepton",
      "PSISRSCALE_UP" : "TTbar_diLepton",
      "PSISRSCALE_DOWN" : "TTbar_diLepton",
      "PSFSRSCALE_UP" : "TTbar_diLepton",
      "PSFSRSCALE_DOWN" : "TTbar_diLepton",
      "BFRAG_UP" : "TTbar_diLepton",
      "BFRAG_DOWN" : "TTbar_diLepton",
      "BSEMILEP_UP" : "TTbar_diLepton",
      "BSEMILEP_DOWN" : "TTbar_diLepton",
      "PDF_ALPHAS_UP" : "TTbar_diLepton",
      "PDF_ALPHAS_DOWN" : "TTbar_diLepton",
      
      "JESAbsoluteMPFBias_UP" : "TTbar_diLepton",
      "JESAbsoluteMPFBias_DOWN" : "TTbar_diLepton",
      "JESAbsoluteScale_UP" : "TTbar_diLepton",
      "JESAbsoluteScale_DOWN" : "TTbar_diLepton",
      "JESAbsoluteStat_UP" : "TTbar_diLepton",
      "JESAbsoluteStat_DOWN" : "TTbar_diLepton",
      "JESFlavorQCD_UP" : "TTbar_diLepton",
      "JESFlavorQCD_DOWN" : "TTbar_diLepton",
      "JESFragmentation_UP" : "TTbar_diLepton",
      "JESFragmentation_DOWN" : "TTbar_diLepton",
      "JESPileUpDataMC_UP" : "TTbar_diLepton",
      "JESPileUpDataMC_DOWN" : "TTbar_diLepton",
      "JESPileUpPtBB_UP" : "TTbar_diLepton",
      "JESPileUpPtBB_DOWN" : "TTbar_diLepton",
      "JESPileUpPtEC1_UP" : "TTbar_diLepton",
      "JESPileUpPtEC1_DOWN" : "TTbar_diLepton",
      "JESPileUpPtRef_UP" : "TTbar_diLepton",
      "JESPileUpPtRef_DOWN" : "TTbar_diLepton",
      "JESRelativeBal_UP" : "TTbar_diLepton",
      "JESRelativeBal_DOWN" : "TTbar_diLepton",
      "JESRelativeFSR_UP" : "TTbar_diLepton",
      "JESRelativeFSR_DOWN" : "TTbar_diLepton",
      "JESRelativeJEREC1_UP" : "TTbar_diLepton",
      "JESRelativeJEREC1_DOWN" : "TTbar_diLepton",
      "JESRelativePtBB_UP" : "TTbar_diLepton",
      "JESRelativePtBB_DOWN" : "TTbar_diLepton",
      "JESRelativePtEC1_UP" : "TTbar_diLepton",
      "JESRelativePtEC1_DOWN" : "TTbar_diLepton",
      "JESRelativeSample_UP" : "TTbar_diLepton",
      "JESRelativeSample_DOWN" : "TTbar_diLepton",
      "JESRelativeStatEC_UP" : "TTbar_diLepton",
      "JESRelativeStatEC_DOWN" : "TTbar_diLepton",
      "JESRelativeStatFSR_UP" : "TTbar_diLepton",
      "JESRelativeStatFSR_DOWN" : "TTbar_diLepton",
      "JESSinglePionECAL_UP" : "TTbar_diLepton",
      "JESSinglePionECAL_DOWN" : "TTbar_diLepton",
      "JESSinglePionHCAL_UP" : "TTbar_diLepton",
      "JESSinglePionHCAL_DOWN" : "TTbar_diLepton",
      "JESTimePtEta_UP" : "TTbar_diLepton",
      "JESTimePtEta_DOWN" : "TTbar_diLepton",
      
      "JESFlavorRealistic_UP" : "TTbar_diLepton",
      "JESFlavorRealistic_DOWN" : "TTbar_diLepton",
      "JESFlavorPureGluon_UP" : "TTbar_diLepton",
      "JESFlavorPureGluon_DOWN" : "TTbar_diLepton",
      "JESFlavorPureQuark_UP" : "TTbar_diLepton",
      "JESFlavorPureQuark_DOWN" : "TTbar_diLepton",
      "JESFlavorPureCharm_UP" : "TTbar_diLepton",
      "JESFlavorPureCharm_DOWN" : "TTbar_diLepton",
      "JESFlavorPureBottom_UP" : "TTbar_diLepton",
      "JESFlavorPureBottom_DOWN" : "TTbar_diLepton",
      
      "JESRelativeBalreg_UP" : "TTbar_diLepton",
      "JESRelativeBalreg_DOWN" : "TTbar_diLepton",
      "JESFlavorQCDreg_UP" : "TTbar_diLepton",
      "JESFlavorQCDreg_DOWN" : "TTbar_diLepton",
      "JESRelativeSampleYear_UP" : "TTbar_diLepton",
      "JESRelativeSampleYear_DOWN" : "TTbar_diLepton",
      "JESAbsoluteYear_UP" : "TTbar_diLepton",
      "JESAbsoluteYear_DOWN" : "TTbar_diLepton",
      "JESAbsolute_UP" : "TTbar_diLepton",
      "JESAbsolute_DOWN" : "TTbar_diLepton",
      "JESBBEC1Year_UP" : "TTbar_diLepton",
      "JESBBEC1Year_DOWN" : "TTbar_diLepton",
      "JESBBEC1_UP" : "TTbar_diLepton",
      "JESBBEC1_DOWN" : "TTbar_diLepton",
   }

#  ~allMC = ["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_ee","DrellYan_mumu","DrellYan_tautau","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
allMC = ["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]

allData2018 = ["DoubleMuon","MuonEG","SingleMuon","EGamma"] 
allData2017 = ["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"] 
allData2016 = ["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"] 

sample_allSyst_dict = {
      "JESTotal_UP" : allMC,
      "JESTotal_DOWN" : allMC,
      "JER_UP" : allMC,
      "JER_DOWN" : allMC,
      "BTAGBC_UP" : allMC,
      "BTAGBC_DOWN" : allMC,
      "BTAGL_UP" : allMC,
      "BTAGL_DOWN" : allMC,
      "ELECTRON_ID_UP" : allMC,
      "ELECTRON_ID_DOWN" : allMC,
      "ELECTRON_RECO_UP" : allMC,
      "ELECTRON_RECO_DOWN" : allMC,
      "ELECTRON_SCALESMEARING_UP" : allMC,
      "ELECTRON_SCALESMEARING_DOWN" : allMC,
      "MUON_ID_UP" : allMC,
      "MUON_ID_DOWN" : allMC,
      "MUON_ISO_UP" : allMC,
      "MUON_ISO_DOWN" : allMC,
      "MUON_SCALE_UP" : allMC,
      "MUON_SCALE_DOWN" : allMC,
      "PU_UP" : allMC,
      "PU_DOWN" : allMC,
      "UNCLUSTERED_UP" : allMC,
      "UNCLUSTERED_DOWN" : allMC,
      "UETUNE_UP" : ["TTbar_diLepton_UETUNE_UP","TTbar_diLepton_tau_UETUNE_UP","TTbar_singleLepton_UETUNE_UP","TTbar_hadronic_UETUNE_UP"],
      "UETUNE_DOWN" : ["TTbar_diLepton_UETUNE_DOWN","TTbar_diLepton_tau_UETUNE_DOWN","TTbar_singleLepton_UETUNE_DOWN","TTbar_hadronic_UETUNE_DOWN"],
      "MATCH_UP" : ["TTbar_diLepton_MATCH_UP","TTbar_diLepton_tau_MATCH_UP","TTbar_singleLepton_MATCH_UP","TTbar_hadronic_MATCH_UP"],
      "MATCH_DOWN" : ["TTbar_diLepton_MATCH_DOWN","TTbar_diLepton_tau_MATCH_DOWN","TTbar_singleLepton_MATCH_DOWN","TTbar_hadronic_MATCH_DOWN"],
      "MTOP169p5" : ["TTbar_diLepton_MTOP169p5","TTbar_diLepton_tau_MTOP169p5","TTbar_singleLepton_MTOP169p5","TTbar_hadronic_MTOP169p5"],
      "MTOP175p5" : ["TTbar_diLepton_MTOP175p5","TTbar_diLepton_tau_MTOP175p5","TTbar_singleLepton_MTOP175p5","TTbar_hadronic_MTOP175p5"],
      "CR1" : ["TTbar_diLepton_CR1","TTbar_diLepton_tau_CR1","TTbar_singleLepton_CR1","TTbar_hadronic_CR1"],
      "CR2" : ["TTbar_diLepton_CR2","TTbar_diLepton_tau_CR2","TTbar_singleLepton_CR2","TTbar_hadronic_CR2"],
      "ERDON" : ["TTbar_diLepton_ERDON","TTbar_diLepton_tau_ERDON","TTbar_singleLepton_ERDON","TTbar_hadronic_ERDON"],
      "TRIG_UP" : allMC,
      "TRIG_DOWN" : allMC,
      "MERENSCALE_UP" : allMC,
      "MERENSCALE_DOWN" : allMC,
      "MEFACSCALE_UP" : allMC,
      "MEFACSCALE_DOWN" : allMC,
      "PSISRSCALE_UP" : allMC,
      "PSISRSCALE_DOWN" : allMC,
      "PSFSRSCALE_UP" : allMC,
      "PSFSRSCALE_DOWN" : allMC,
      "BFRAG_UP" : allMC,
      "BFRAG_DOWN" : allMC,
      "BSEMILEP_UP" : allMC,
      "BSEMILEP_DOWN" : allMC,
      "PDF_ALPHAS_UP" : allMC,
      "PDF_ALPHAS_DOWN" : allMC,
      "TOP_PT" : allMC,
      "L1PREFIRING_UP" : allMC,
      "L1PREFIRING_DOWN" : allMC,
      
      "JESAbsoluteMPFBias_UP" : allMC,
      "JESAbsoluteMPFBias_DOWN" : allMC,
      "JESAbsoluteScale_UP" : allMC,
      "JESAbsoluteScale_DOWN" : allMC,
      "JESAbsoluteStat_UP" : allMC,
      "JESAbsoluteStat_DOWN" : allMC,
      "JESFlavorQCD_UP" : allMC,
      "JESFlavorQCD_DOWN" : allMC,
      "JESFragmentation_UP" : allMC,
      "JESFragmentation_DOWN" : allMC,
      "JESPileUpDataMC_UP" : allMC,
      "JESPileUpDataMC_DOWN" : allMC,
      "JESPileUpPtBB_UP" : allMC,
      "JESPileUpPtBB_DOWN" : allMC,
      "JESPileUpPtEC1_UP" : allMC,
      "JESPileUpPtEC1_DOWN" : allMC,
      "JESPileUpPtRef_UP" : allMC,
      "JESPileUpPtRef_DOWN" : allMC,
      "JESRelativeBal_UP" : allMC,
      "JESRelativeBal_DOWN" : allMC,
      "JESRelativeFSR_UP" : allMC,
      "JESRelativeFSR_DOWN" : allMC,
      "JESRelativeJEREC1_UP" : allMC,
      "JESRelativeJEREC1_DOWN" : allMC,
      "JESRelativePtBB_UP" : allMC,
      "JESRelativePtBB_DOWN" : allMC,
      "JESRelativePtEC1_UP" : allMC,
      "JESRelativePtEC1_DOWN" : allMC,
      "JESRelativeSample_UP" : allMC,
      "JESRelativeSample_DOWN" : allMC,
      "JESRelativeStatEC_UP" : allMC,
      "JESRelativeStatEC_DOWN" : allMC,
      "JESRelativeStatFSR_UP" : allMC,
      "JESRelativeStatFSR_DOWN" : allMC,
      "JESSinglePionECAL_UP" : allMC,
      "JESSinglePionECAL_DOWN" : allMC,
      "JESSinglePionHCAL_UP" : allMC,
      "JESSinglePionHCAL_DOWN" : allMC,
      "JESTimePtEta_UP" : allMC,
      "JESTimePtEta_DOWN" : allMC,
      
      "JESFlavorRealistic_UP" : allMC,
      "JESFlavorRealistic_DOWN" : allMC,
      "JESFlavorPureGluon_UP" : allMC,
      "JESFlavorPureGluon_DOWN" : allMC,
      "JESFlavorPureQuark_UP" : allMC,
      "JESFlavorPureQuark_DOWN" : allMC,
      "JESFlavorPureCharm_UP" : allMC,
      "JESFlavorPureCharm_DOWN" : allMC,
      "JESFlavorPureBottom_UP" : allMC,
      "JESFlavorPureBottom_DOWN" : allMC,
      
      "JESRelativeBalreg_UP" : allMC,
      "JESRelativeBalreg_DOWN" : allMC,
      "JESFlavorQCDreg_UP" : allMC,
      "JESFlavorQCDreg_DOWN" : allMC,
      "JESRelativeSampleYear_UP" : allMC,
      "JESRelativeSampleYear_DOWN" : allMC,
      "JESAbsoluteYear_UP" : allMC,
      "JESAbsoluteYear_DOWN" : allMC,
      "JESAbsolute_UP" : allMC,
      "JESAbsolute_DOWN" : allMC,
      "JESBBEC1Year_UP" : allMC,
      "JESBBEC1Year_DOWN" : allMC,
      "JESBBEC1_UP" : allMC,
      "JESBBEC1_DOWN" : allMC
}

tunfold_syst = [
"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP169p5","MTOP175p5","CR1","CR2","ERDON","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN",

"JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN",

"JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFlavorPureGluon_UP","JESFlavorPureGluon_DOWN","JESFlavorPureQuark_UP","JESFlavorPureQuark_DOWN","JESFlavorPureCharm_UP","JESFlavorPureCharm_DOWN","JESFlavorPureBottom_UP","JESFlavorPureBottom_DOWN",

]


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

def get_dataBasePath_dCache(year,dcap=False):      #return dataBasePath on dCache for given year
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   if dcap:
      return "dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname"),year,config["input"]["version"])
   else:
      return "root://grid-cms-xrootd.physik.rwth-aachen.de///store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname"),year,config["input"]["version"])

def printSubmitInfo(args):
   print "Running "+args.m
   print "Systematic: "+args.s
   print "Process "+str(args.f*100)+"% of events"
   
def createLogPath(args):
   logpath="logs/"+args.y+"/"+args.s+"/"+str(args.f)+"/"+args.m
   if not os.path.exists(logpath):
      try:
         os.makedirs(logpath)
      except OSError as exc: # Guard against race condition
         if exc.errno != errno.EEXIST:
            raise
   return logpath
   
def submit(args,toProcess_mc,toProcess_data,toProcess_signal,disableConfirm=False):
   
   printSubmitInfo(args)
   
   # get correct sample for bTagEff
   if (args.m == "bTagEff"):
      toProcess_mc = [bTagEff_sample_syst_dict[args.s]]
      toProcess_data = []
      toProcess_signal = []
      if (any("TTbar_diLepton" in s for s in toProcess_mc)==False or len(toProcess_mc)>1 or toProcess_data or toProcess_signal):
         print "Error: bTagEff should only be submitted with TTbar_dilepton"
         sys.exit(1)

   print toProcess_mc
   print toProcess_data
   print toProcess_signal
   
   # set dataBasePath depending on dCache options
   if args.scratchInput==False and args.copyDCache==False:
      print "Use dCache input"
      dataBasePath = "-d"+get_dataBasePath_dCache(args.y)
   elif args.copyDCache:
      print "Copy input from dCache to condor node"
      dataBasePath = "-d$TMP/"
   else:
      print "Use scratch input"
      dataBasePath = ""
   
   # Ask if selected settings are correct
   runSubmit = True
   if disableConfirm==False:
      correctSamples = input("If you want to continue with the selected setting, enter 1:\n")
      if (correctSamples != 1):
         print "Abort Submission"
         runSubmit = False
   
   if runSubmit:
      # For single submit check if only one dataset is selected and then ask for file nr
      if (args.SingleSubmit and ((len(toProcess_mc)+len(toProcess_data)+len(toProcess_signal))==1)):
         singleFileNR = input("Please input the filenumber to submit:\n")
      elif args.SingleSubmit:
         print "SingleSubmit can only be used if one dataset is selected"
         exit(98)

      # not data processing if systematic shift is choosen
      if (toProcess_data and args.s!="Nominal" and args.s!="met40Cut"):
         toProcess_data=[]
         print "!!!!!!!!!!!!!!!!Data is not processed with systematic shift!!!!!!!!!!!!!!!!!!!!"
      
      # triggerEff only submitted with TTbar_dilepton and MET
      if (args.m == "triggerEff"):
         if ("TTbar_diLepton" not in toProcess_mc or len(toProcess_mc)>1 or toProcess_signal or "MET" not in toProcess_data):
            print "Error: triggerEff should only be submitted with TTbar_dilepton and MET"
            sys.exit(1)
      
      # create logpath if not existing
      logpath = createLogPath(args)

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
                     
                     if args.copyDCache:     #set path to dCache input in case of copying to condor node
                        inputPath = get_dataBasePath_dCache(args.y,True)+get_fileName(x,args.y,fileNR)
                        inputPath = inputPath.replace(" ", "")
                     else:
                        inputPath = ""
                                       
                     with open(submitFile,"w") as f:   # write condor submit
                        #  ~f.write("""
#  ~Universe                = vanilla
#  ~Executable              = run.sh
#  ~Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {8} {9} {11} {12}
#  ~Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
#  ~Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
#  ~Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
#  ~use_x509userproxy       = true
#  ~Request_Memory          = {4} Mb
#  ~Requirements            = (TARGET.CpuFamily > 6) && (TARGET.Machine != "lxcip16.physik.rwth-aachen.de")  {10}
#  ~Queue
   #  ~""".format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1),dataBasePath,inputPath,"\nRank = CpuFamily" if(x=="TTbar_diLepton") else "", getPath("cmsswBasePath"), getPath("frameworkBasePath")),)
                        f.write("""
Universe                = vanilla
Executable              = run.sh
Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {8} {9} {11} {12}
Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
use_x509userproxy       = true
Request_Memory          = {4} Mb
Requirements            = (TARGET.Machine == "lxblade33.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade34.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade35.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade36.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade37.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade38.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade39.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade40.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade41.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade42.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade43.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade44.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade45.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade46.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade47.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade48.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade49.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade50.physik.rwth-aachen.de") || (substr(TARGET.Machine,0,4) == "lx1b")
Queue
   """.format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1),dataBasePath,inputPath,"\nRank = CpuFamily" if(x=="TTbar_diLepton") else "", getPath("cmsswBasePath"), getPath("frameworkBasePath")),)
                     subprocess.call(["condor_submit", submitFile])


def submitTUnfold(args,systematics):
   for syst in systematics:
      
      # override systematic
      args.s = syst
      
      printSubmitInfo(args)
      
      # create logpath if not existing
      logpath = createLogPath(args)
      
      # loop over selected systematics and submit corresponding jobs
                           
      for logFile in glob.glob(logpath+"/*"):  # Remove old logs
         if os.path.exists(logFile):
            os.remove(logFile)
      
      submitFile = logpath+"/"+args.m+".submit"    #define submit file
                        
      with open(submitFile,"w") as f:   # write condor submit
         f.write("""
Universe                = vanilla
Executable              = runTUnfold.sh
Arguments               = -f{0} {1} {2} -s{3}
Log                     = logs/{2}/{3}/{0}/{1}/{1}.log
Output                  = logs/{2}/{3}/{0}/{1}/{1}.out
Error                   = logs/{2}/{3}/{0}/{1}/{1}.error
use_x509userproxy       = true
Queue
""".format(str(args.f),args.m,args.y,args.s),)
      subprocess.call(["condor_submit", submitFile])
   
   

if __name__ == "__main__":
   
   #############################################
   # Select datasets to process
   #############################################
   #  ~toProcess_mc=["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
   #  ~toProcess_mc=["ZZ"]
   #  ~toProcess_mc=["DrellYan_NLO"]
   #  ~toProcess_mc=["TTbar_diLepton"]
   #  ~toProcess_mc=["TTbar_diLepton_tau_MATCH_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_UP","TTbar_diLepton_tau_UETUNE_UP","TTbar_singleLepton_UETUNE_UP","TTbar_hadronic_UETUNE_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_UETUNE_DOWN","TTbar_diLepton_tau_UETUNE_DOWN","TTbar_singleLepton_UETUNE_DOWN","TTbar_hadronic_UETUNE_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_UP","TTbar_diLepton_tau_MATCH_UP","TTbar_singleLepton_MATCH_UP","TTbar_hadronic_MATCH_UP"]
   #  ~toProcess_mc=["TTbar_diLepton_MATCH_DOWN","TTbar_diLepton_tau_MATCH_DOWN","TTbar_singleLepton_MATCH_DOWN","TTbar_hadronic_MATCH_DOWN"]
   #  ~toProcess_mc=["TTbar_diLepton_ERDON","TTbar_diLepton_tau_ERDON","TTbar_singleLepton_ERDON","TTbar_hadronic_ERDON"]
   #  ~toProcess_mc=["TTbar_hadronic_CR1"]
   #  ~toProcess_mc=["TTbar_diLepton_CR1","TTbar_diLepton_tau_CR1","TTbar_singleLepton_CR1","TTbar_hadronic_CR1"]
   #  ~toProcess_mc=["TTbar_diLepton_CR2","TTbar_diLepton_tau_CR2","TTbar_singleLepton_CR2","TTbar_hadronic_CR2"]
   #  ~toProcess_mc=["TTbar_diLepton_MTOP169p5","TTbar_diLepton_tau_MTOP169p5","TTbar_singleLepton_MTOP169p5","TTbar_hadronic_MTOP169p5"]
   #  ~toProcess_mc=["TTbar_diLepton_MTOP175p5","TTbar_diLepton_tau_MTOP175p5","TTbar_singleLepton_MTOP175p5","TTbar_hadronic_MTOP175p5"]
   #  ~toProcess_mc=[]
   toProcess_mc=allMC
   
   #  ~toProcess_data=["DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"]
   #  ~toProcess_data=["DoubleMuon","MuonEG","SingleMuon","EGamma"]      #2018
   #  ~toProcess_data=["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"]       #2017, 2016
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
   parser.add_argument('--copyDCache', action='store_true', default=True, help="Copy nTuples stored on DCache to node before running the code.")
   parser.add_argument('--SingleSubmit', action='store_true' )
   parser.add_argument('--bTagEff_complete', action='store_true', default=False, help="Submits bTagEff jobs with all relevant systematics (use with care!)")
   parser.add_argument('--distributions_complete', action='store_true', default=False, help="Submits distributions jobs with all relevant systematics (use with care!)")
   parser.add_argument('--tunfold_complete', action='store_true', default=False, help="Submits tunfold jobs with all relevant systematics (use with care!)")
   parser.add_argument('--tunfold_pdf_complete', action='store_true', default=False, help="Submits tunfold jobs with all pdf shifts (use with care!)")
   parser.add_argument('--pdf_complete', action='store_true', default=False, help="Submits distributions jobs with all pdf shifts (use with care!)")
   parser.add_argument('--noConfirmation', action='store_true', default=False, help="Disables keyboard input befor submission")

   args = parser.parse_args()
   
   if (args.bTagEff_complete == False and args.distributions_complete==False and args.pdf_complete==False):
      if (args.m == "TUnfold_binning"):
         if (args.tunfold_complete):
            submitTUnfold(args,tunfold_syst)
         elif (args.tunfold_pdf_complete):
            for varNumber in xrange(1,51):
               for var in ["UP","DOWN"]:
                  args.s = "PDF_{}_{}".format(varNumber,var)
                  submitTUnfold(args,{args.s})
         else:
            submitTUnfold(args,{args.s})
      else:
         submit(args,toProcess_mc,toProcess_data,toProcess_signal)
   elif (args.bTagEff_complete):
      if (args.m == "bTagEff"):
         for syst in bTagEff_sample_syst_dict.keys():
            args.s = syst
            submit(args,toProcess_mc,toProcess_data,toProcess_signal,args.noConfirmation)
      else:
         print "bTagEff_complete can only be used if bTagEff is selected as module!"
         exit(98)
   elif (args.distributions_complete):
      if (args.m == "" or args.m == "distributions"):
         print "Submit nominal"
         if (args.y=="2018"):
            submit(args,allMC,allData2018,[])
         elif (args.y=="2017"):
            submit(args,allMC,allData2017,[])
         elif (args.y=="2016_preVFP" or args.y == "2016_postVFP"):
            submit(args,allMC,allData2016,[])
         else:
            print args.y+" does not match correct year"
         for syst in sample_allSyst_dict.keys():
            args.s = syst
            submit(args,sample_allSyst_dict[syst],[],[],args.noConfirmation)
      else:
         print "distributions_complete can only be used if distributions is selected as module!"
         exit(98)
   elif (args.pdf_complete):
      if (args.m == "" or args.m == "distributions"):
         for varNumber in xrange(1,51):
            for var in ["UP","DOWN"]:
               args.s = "PDF_{}_{}".format(varNumber,var)
               submit(args,["TTbar_diLepton"],[],[],args.noConfirmation)
      else:
         print "pdf_complete can only be used if distributions is selected as module!"
         exit(98)
      
   
