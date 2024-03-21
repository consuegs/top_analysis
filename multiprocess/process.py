#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess
import sys
import configparser
import sys
import re
sys.path.append("../users")
from getPath import getPath

bTagEff_sample_syst_dict = {
      "Nominal" : "TTbar_diLepton",
      "BSEMILEP_UP" : "TTbar_diLepton",
      "BSEMILEP_DOWN" : "TTbar_diLepton",
      "CR1" : "TTbar_diLepton_CR1",
      "CR2" : "TTbar_diLepton_CR2",
      "ERDON" : "TTbar_diLepton_ERDON",
      "JEREta0_UP" : "TTbar_diLepton",
      "JEREta0_DOWN" : "TTbar_diLepton",
      "JEREta1_UP" : "TTbar_diLepton",
      "JEREta1_DOWN" : "TTbar_diLepton",
      "JESAbsolute_UP" : "TTbar_diLepton",
      "JESAbsolute_DOWN" : "TTbar_diLepton",
      "JESAbsoluteYear_UP" : "TTbar_diLepton",
      "JESAbsoluteYear_DOWN" : "TTbar_diLepton",
      "JESBBEC1_UP" : "TTbar_diLepton",
      "JESBBEC1_DOWN" : "TTbar_diLepton",
      "JESBBEC1Year_UP" : "TTbar_diLepton",
      "JESBBEC1Year_DOWN" : "TTbar_diLepton",
      "JESFlavorRealistic_UP" : "TTbar_diLepton",
      "JESFlavorRealistic_DOWN" : "TTbar_diLepton",
      "JESRelativeBalreg_UP" : "TTbar_diLepton",
      "JESRelativeBalreg_DOWN" : "TTbar_diLepton",
      "JESRelativeSampleYear_UP" : "TTbar_diLepton",
      "JESRelativeSampleYear_DOWN" : "TTbar_diLepton",
      "JETPILEUPID_UP" : "TTbar_diLepton",
      "JETPILEUPID_DOWN" : "TTbar_diLepton",
      "MATCH_UP" : "TTbar_diLepton_MATCH_UP",
      "MATCH_DOWN" : "TTbar_diLepton_MATCH_DOWN",
      "MATCH_DCTR_UP" : "TTbar_diLepton",
      "MATCH_DCTR_DOWN" : "TTbar_diLepton",
      "MERENSCALE_UP" : "TTbar_diLepton",
      "MERENSCALE_DOWN" : "TTbar_diLepton",
      "MEFACSCALE_UP" : "TTbar_diLepton",
      "MEFACSCALE_DOWN" : "TTbar_diLepton",
      "MESCALE_UP" : "TTbar_diLepton",
      "MESCALE_DOWN" : "TTbar_diLepton",
      "MTOP169p5" : "TTbar_diLepton_MTOP169p5",
      "MTOP175p5" : "TTbar_diLepton_MTOP175p5",
      "PDF_ALPHAS_UP" : "TTbar_diLepton",
      "PDF_ALPHAS_DOWN" : "TTbar_diLepton",
      "PSISRSCALE_UP" : "TTbar_diLepton",
      "PSISRSCALE_DOWN" : "TTbar_diLepton",
      "PSFSRSCALE_UP" : "TTbar_diLepton",
      "PSFSRSCALE_DOWN" : "TTbar_diLepton",
      "TOP_PT" : "TTbar_diLepton",
      "UETUNE_UP" : "TTbar_diLepton_UETUNE_UP",
      "UETUNE_DOWN" : "TTbar_diLepton_UETUNE_DOWN",

      
      #  ~"JESAbsoluteMPFBias_UP" : "TTbar_diLepton",
      #  ~"JESAbsoluteMPFBias_DOWN" : "TTbar_diLepton",
      #  ~"JESAbsoluteScale_UP" : "TTbar_diLepton",
      #  ~"JESAbsoluteScale_DOWN" : "TTbar_diLepton",
      #  ~"JESAbsoluteStat_UP" : "TTbar_diLepton",
      #  ~"JESAbsoluteStat_DOWN" : "TTbar_diLepton",
      #  ~"JESFlavorQCD_UP" : "TTbar_diLepton",
      #  ~"JESFlavorQCD_DOWN" : "TTbar_diLepton",
      #  ~"JESFragmentation_UP" : "TTbar_diLepton",
      #  ~"JESFragmentation_DOWN" : "TTbar_diLepton",
      #  ~"JESPileUpDataMC_UP" : "TTbar_diLepton",
      #  ~"JESPileUpDataMC_DOWN" : "TTbar_diLepton",
      #  ~"JESPileUpPtBB_UP" : "TTbar_diLepton",
      #  ~"JESPileUpPtBB_DOWN" : "TTbar_diLepton",
      #  ~"JESPileUpPtEC1_UP" : "TTbar_diLepton",
      #  ~"JESPileUpPtEC1_DOWN" : "TTbar_diLepton",
      #  ~"JESPileUpPtRef_UP" : "TTbar_diLepton",
      #  ~"JESPileUpPtRef_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeBal_UP" : "TTbar_diLepton",
      #  ~"JESRelativeBal_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeFSR_UP" : "TTbar_diLepton",
      #  ~"JESRelativeFSR_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeJEREC1_UP" : "TTbar_diLepton",
      #  ~"JESRelativeJEREC1_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativePtBB_UP" : "TTbar_diLepton",
      #  ~"JESRelativePtBB_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativePtEC1_UP" : "TTbar_diLepton",
      #  ~"JESRelativePtEC1_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeSample_UP" : "TTbar_diLepton",
      #  ~"JESRelativeSample_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeStatEC_UP" : "TTbar_diLepton",
      #  ~"JESRelativeStatEC_DOWN" : "TTbar_diLepton",
      #  ~"JESRelativeStatFSR_UP" : "TTbar_diLepton",
      #  ~"JESRelativeStatFSR_DOWN" : "TTbar_diLepton",
      #  ~"JESSinglePionECAL_UP" : "TTbar_diLepton",
      #  ~"JESSinglePionECAL_DOWN" : "TTbar_diLepton",
      #  ~"JESSinglePionHCAL_UP" : "TTbar_diLepton",
      #  ~"JESSinglePionHCAL_DOWN" : "TTbar_diLepton",
      #  ~"JESTimePtEta_UP" : "TTbar_diLepton",
      #  ~"JESTimePtEta_DOWN" : "TTbar_diLepton",
      #  ~"JESFlavorRealistic_UP" : "TTbar_diLepton",
      #  ~"JESFlavorRealistic_DOWN" : "TTbar_diLepton",
      
      # Non-split uncertainties (for combine)
      #  ~"JER_UP" : "TTbar_diLepton",
      #  ~"JER_DOWN" : "TTbar_diLepton",
      
      # Additional studies
      #  ~"applyJerMET" : "TTbar_diLepton",
      #  ~"JERMET_UP" : "TTbar_diLepton",
      #  ~"JERMET_DOWN" : "TTbar_diLepton",
   }

#  ~allMC = ["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_ee","DrellYan_mumu","DrellYan_tautau","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
#  ~allMC = ["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","DrellYan_M10to50_NLO","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
allMC = ["SingleTop_DS"]

allData2018 = ["DoubleMuon","MuonEG","SingleMuon","EGamma"] 
allData2017 = ["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"] 
allData2016 = ["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"]

nominalTypeSysts = ["Nominal","met40Cut","removeMLLcut","jetPileupIDapplied","applyJerMET","applyJetVetoMaps","applyJetVetoMaps_subleading","applyJetVetoMaps_leading","applyJetVetoMaps_loose","applyJetVetoMaps_cleanedJets","applyJetVetoMaps_HEM1516","applyGenLevel_DeltaRcut","useDNNnoMETcut","useDNNmumu","useDNNnoMetCutDY"]

sample_allSyst_dict = {
      "BSEMILEP_UP" : allMC,
      "BSEMILEP_DOWN" : allMC,
      "BTAGBC_CORR_UP" : allMC,
      "BTAGBC_CORR_DOWN" : allMC,
      "BTAGBC_UNCORR_UP" : allMC,
      "BTAGBC_UNCORR_DOWN" : allMC,
      "BTAGL_CORR_UP" : allMC,
      "BTAGL_CORR_DOWN" : allMC,
      "BTAGL_UNCORR_UP" : allMC,
      "BTAGL_UNCORR_DOWN" : allMC,
      "CR1" : ["TTbar_diLepton_CR1","TTbar_diLepton_tau_CR1","TTbar_singleLepton_CR1","TTbar_hadronic_CR1"],
      "CR2" : ["TTbar_diLepton_CR2","TTbar_diLepton_tau_CR2","TTbar_singleLepton_CR2","TTbar_hadronic_CR2"],
      "ERDON" : ["TTbar_diLepton_ERDON","TTbar_diLepton_tau_ERDON","TTbar_singleLepton_ERDON","TTbar_hadronic_ERDON"],
      "ELECTRON_ID_UP" : allMC,
      "ELECTRON_ID_DOWN" : allMC,
      "ELECTRON_RECO_UP" : allMC,
      "ELECTRON_RECO_DOWN" : allMC,
      "ELECTRON_SCALESMEARING_UP" : allMC,
      "ELECTRON_SCALESMEARING_DOWN" : allMC,
      "JEREta0_UP" : allMC,
      "JEREta0_DOWN" : allMC,
      "JEREta1_UP" : allMC,
      "JEREta1_DOWN" : allMC,
      "JESAbsolute_UP" : allMC,
      "JESAbsolute_DOWN" : allMC,
      "JESAbsoluteYear_UP" : allMC,
      "JESAbsoluteYear_DOWN" : allMC,
      "JESBBEC1_UP" : allMC,
      "JESBBEC1_DOWN" : allMC,
      "JESBBEC1Year_UP" : allMC,
      "JESBBEC1Year_DOWN" : allMC,
      "JESFlavorRealistic_UP" : allMC,
      "JESFlavorRealistic_DOWN" : allMC,
      "JESRelativeBalreg_UP" : allMC,
      "JESRelativeBalreg_DOWN" : allMC,
      "JESRelativeSampleYear_UP" : allMC,
      "JESRelativeSampleYear_DOWN" : allMC,
      "JETPILEUPID_UP" : allMC,
      "JETPILEUPID_DOWN" : allMC,
      "L1PREFIRING_UP" : allMC,
      "L1PREFIRING_DOWN" : allMC,
      "MATCH_UP" : ["TTbar_diLepton_MATCH_UP","TTbar_diLepton_tau_MATCH_UP","TTbar_singleLepton_MATCH_UP","TTbar_hadronic_MATCH_UP"],
      "MATCH_DOWN" : ["TTbar_diLepton_MATCH_DOWN","TTbar_diLepton_tau_MATCH_DOWN","TTbar_singleLepton_MATCH_DOWN","TTbar_hadronic_MATCH_DOWN"],
      "MATCH_DCTR_UP" : ["TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"],
      "MATCH_DCTR_DOWN" : ["TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"],
      "MERENSCALE_UP" : allMC,
      "MERENSCALE_DOWN" : allMC,
      "MEFACSCALE_UP" : allMC,
      "MEFACSCALE_DOWN" : allMC,
      "MESCALE_UP" : allMC,
      "MESCALE_DOWN" : allMC,
      "MTOP169p5" : ["TTbar_diLepton_MTOP169p5","TTbar_diLepton_tau_MTOP169p5","TTbar_singleLepton_MTOP169p5","TTbar_hadronic_MTOP169p5"],
      "MTOP175p5" : ["TTbar_diLepton_MTOP175p5","TTbar_diLepton_tau_MTOP175p5","TTbar_singleLepton_MTOP175p5","TTbar_hadronic_MTOP175p5"],
      "MUON_ID_STAT_UP" : allMC,
      "MUON_ID_STAT_DOWN" : allMC,
      "MUON_ID_SYST_UP" : allMC,
      "MUON_ID_SYST_DOWN" : allMC,
      "MUON_ISO_STAT_UP" : allMC,
      "MUON_ISO_STAT_DOWN" : allMC,
      "MUON_ISO_SYST_UP" : allMC,
      "MUON_ISO_SYST_DOWN" : allMC,
      "MUON_SCALE_UP" : allMC,
      "MUON_SCALE_DOWN" : allMC,
      "PDF_ALPHAS_UP" : allMC,
      "PDF_ALPHAS_DOWN" : allMC,
      "PSISRSCALE_UP" : allMC,
      "PSISRSCALE_DOWN" : allMC,
      "PSFSRSCALE_UP" : allMC,
      "PSFSRSCALE_DOWN" : allMC,
      "PU_UP" : allMC,
      "PU_DOWN" : allMC,
      "TOP_PT" : allMC,
      "TRIG_UP" : allMC,
      "TRIG_DOWN" : allMC,
      "TWDS" : ["SingleTop_TWDS"],
      "UETUNE_UP" : ["TTbar_diLepton_UETUNE_UP","TTbar_diLepton_tau_UETUNE_UP","TTbar_singleLepton_UETUNE_UP","TTbar_hadronic_UETUNE_UP"],
      "UETUNE_DOWN" : ["TTbar_diLepton_UETUNE_DOWN","TTbar_diLepton_tau_UETUNE_DOWN","TTbar_singleLepton_UETUNE_DOWN","TTbar_hadronic_UETUNE_DOWN"],
      "UNCLUSTERED_UP" : allMC,
      "UNCLUSTERED_DOWN" : allMC,
            
      #  ~"JESAbsoluteMPFBias_UP" : allMC,
      #  ~"JESAbsoluteMPFBias_DOWN" : allMC,
      #  ~"JESAbsoluteScale_UP" : allMC,
      #  ~"JESAbsoluteScale_DOWN" : allMC,
      #  ~"JESAbsoluteStat_UP" : allMC,
      #  ~"JESAbsoluteStat_DOWN" : allMC,
      #  ~"JESFlavorQCD_UP" : allMC,
      #  ~"JESFlavorQCD_DOWN" : allMC,
      #  ~"JESFragmentation_UP" : allMC,
      #  ~"JESFragmentation_DOWN" : allMC,
      #  ~"JESPileUpDataMC_UP" : allMC,
      #  ~"JESPileUpDataMC_DOWN" : allMC,
      #  ~"JESPileUpPtBB_UP" : allMC,
      #  ~"JESPileUpPtBB_DOWN" : allMC,
      #  ~"JESPileUpPtEC1_UP" : allMC,
      #  ~"JESPileUpPtEC1_DOWN" : allMC,
      #  ~"JESPileUpPtRef_UP" : allMC,
      #  ~"JESPileUpPtRef_DOWN" : allMC,
      #  ~"JESRelativeBal_UP" : allMC,
      #  ~"JESRelativeBal_DOWN" : allMC,
      #  ~"JESRelativeFSR_UP" : allMC,
      #  ~"JESRelativeFSR_DOWN" : allMC,
      #  ~"JESRelativeJEREC1_UP" : allMC,
      #  ~"JESRelativeJEREC1_DOWN" : allMC,
      #  ~"JESRelativePtBB_UP" : allMC,
      #  ~"JESRelativePtBB_DOWN" : allMC,
      #  ~"JESRelativePtEC1_UP" : allMC,
      #  ~"JESRelativePtEC1_DOWN" : allMC,
      #  ~"JESRelativeSample_UP" : allMC,
      #  ~"JESRelativeSample_DOWN" : allMC,
      #  ~"JESRelativeStatEC_UP" : allMC,
      #  ~"JESRelativeStatEC_DOWN" : allMC,
      #  ~"JESRelativeStatFSR_UP" : allMC,
      #  ~"JESRelativeStatFSR_DOWN" : allMC,
      #  ~"JESSinglePionECAL_UP" : allMC,
      #  ~"JESSinglePionECAL_DOWN" : allMC,
      #  ~"JESSinglePionHCAL_UP" : allMC,
      #  ~"JESSinglePionHCAL_DOWN" : allMC,
      #  ~"JESTimePtEta_UP" : allMC,
      #  ~"JESTimePtEta_DOWN" : allMC,
      
      # Non-split uncertainties (for combine)
      #  ~"BTAGBC_UP" : allMC,
      #  ~"BTAGBC_DOWN" : allMC,
      #  ~"BTAGL_UP" : allMC,
      #  ~"BTAGL_DOWN" : allMC,
      #  ~"JER_UP" : allMC,
      #  ~"JER_DOWN" : allMC,
      #  ~"MUON_ID_UP" : allMC,
      #  ~"MUON_ID_DOWN" : allMC,
      #  ~"MUON_ISO_UP" : allMC,
      #  ~"MUON_ISO_DOWN" : allMC
      
      #  ~"JERMET_UP" : allMC,
      #  ~"JERMET_DOWN" : allMC,
}

#  ~tunfold_syst = [
#  ~"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"
#  ~]
tunfold_syst = [
"Nominal"
]
#  ~tunfold_syst = ["CR1","CR2","ERDON"]
#  ~tunfold_syst = ["MATCH_DCTR_UP","MATCH_DCTR_DOWN"]

#  ~tunfold_syst = ["JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN"]


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

def get_version(year):      #checks config for name of file
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   return config["input"]["version"]

def get_dataBasePath_dCache(year,curl=True,dcap=False,extern=False):      #return dataBasePath on dCache for given year
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   if (args.y.find("2016")>=0):
      user="teroerde"
   else:
      user=None
   if dcap:
      if extern:
         return "dcap://grid-dcap-extern.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname",user),year,config["input"]["version"])
      else:
         return "dcap://grid-dcap.physik.rwth-aachen.de/pnfs/physik.rwth-aachen.de/cms/store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname",user),year,config["input"]["version"])
   elif curl:
      return "curlsimple://grid-webdav.physik.rwth-aachen.de:2889/store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname",user),year,config["input"]["version"])
   else:
      return "root://grid-cms-xrootd.physik.rwth-aachen.de///store/user/{0}/mergedNtuple/{1}/{2}/".format(getPath("gridname",user),year,config["input"]["version"])

def printSubmitInfo(args):
   print("Running "+args.m)
   print("Systematic: "+args.s)
   print("Process "+str(args.f*100)+"% of events")
   
def createLogPath(args):
   logpath = "logs/"+args.y+"/"+args.s+"/"+str(args.f)+"/"+args.m
   if not os.path.exists(logpath):
      try:
         os.makedirs(logpath)
      except OSError as exc: # Guard against race condition
         if exc.errno != errno.EEXIST:
            raise
   return logpath

def removeJES_HEM(year):
   if year != "2018":
      key = "JESUserDefinedHEM1516_DOWN"
      if key in tunfold_syst:
         tunfold_syst.remove(key)
      if key in bTagEff_sample_syst_dict:
         del bTagEff_sample_syst_dict[key]
      if key in sample_allSyst_dict:
         del sample_allSyst_dict[key]

def checkCEjobs():      #check how many jobs are present in both CE and choose with CE to submit to
   out = subprocess.check_output(["condor_status", "-grid"])
   jobs_ce1 = 0
   jobs_ce2 = 0
   for line in out.split("\n"):
      if "ce-1" in line:
         line = re.sub(" +"," ",line)
         jobs_ce1 += int(line.split(" ")[3])
      elif "ce-2" in line:
         line = re.sub(" +"," ",line)
         jobs_ce2 += int(line.split(" ")[3])
   
   if jobs_ce1>jobs_ce2:
      return "ce-2"
   else:
      return "ce-1"

def uploadCompressedCMSSW():    #compress and upload CMSSW to dCache to run on grid
   runDir = os.getcwd()
   filePath = runDir+"/inputs/"
   if not os.path.exists(filePath):
      try:
         os.makedirs(filePath)
      except OSError as exc: # Guard against race condition
         if exc.errno != errno.EEXIST:
            raise
   CMSSW_path = os.getenv('CMSSW_BASE')
   CMSSW_version = CMSSW_path.split("/")[-1]
   filePath+=CMSSW_version+".tgz"
   
   print("Compressing "+CMSSW_version)
   try:
      os.chdir(CMSSW_path+"/../")
      subprocess.call(["tar","-zcvf",filePath,CMSSW_version],stdout=open(os.devnull, 'wb'))
      os.chdir(runDir)
   except Exception as error:
      print("Compressing failed")
      errno, errstr = error.args
      print(errstr)
      sys.exit(20)
   
   print("Uploading "+CMSSW_version+" to dCache")
   targetPath = "gridJobInputs/"+CMSSW_version+".tgz"
   command = "eval `scram unsetenv -sh`;","gfal-copy",filePath,"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/".format(getPath("dCacheBasePath"))+targetPath,"-f","-r"
   command = " ".join(command)
   try:
      subprocess.call(command,shell=True)
   except Exception as error:
      print("Uploading failed")
      errno, errstr = error.args
      print(errstr)
      sys.exit(21)

def uploadCompressedFW():    #compress and upload local Framework to dCache to run on grid
   runDir = os.getcwd()
   filePath = runDir+"/inputs/"
   if not os.path.exists(filePath):
      try:
         os.makedirs(filePath)
      except OSError as exc: # Guard against race condition
         if exc.errno != errno.EEXIST:
            raise
   
   FW_files = ["CMakeLists.txt","cmake","config2018.ini","config2017.ini","config2016_preVFP.ini","config2016_postVFP.ini","data","src"]
   filePath+="FW.tgz"
   FW_path = getPath("frameworkBasePath")
   
   print("Compressing framework from "+FW_path)
   try:
      os.chdir(FW_path)
      subprocess.call(["tar","-hzcvf",filePath]+FW_files,stdout=open(os.devnull, 'wb'))
      os.chdir(runDir)
   except Exception as error:
      print("Compressing failed")
      errno, errstr = error.args
      print(errstr)
      sys.exit(30)
   
   print("Uploading framework to dCache")
   targetPath = "gridJobInputs/FW.tgz"
   command = "eval `scram unsetenv -sh`;","gfal-copy",filePath,"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/".format(getPath("dCacheBasePath"))+targetPath,"-f","-r"
   command = " ".join(command)
   try:
      subprocess.call(command,shell=True)
   except Exception as error:
      print("Uploading failed")
      errno, errstr = error.args
      print(errstr)
      sys.exit(31)
   
def submit(args,toProcess_mc,toProcess_data,toProcess_signal,disableConfirm=False):
   
   printSubmitInfo(args)
   
   # get correct sample for bTagEff
   if (args.m == "bTagEff"):
      toProcess_mc = [bTagEff_sample_syst_dict[args.s]]
      toProcess_data = []
      toProcess_signal = []
      if (any("TTbar_diLepton" in s for s in toProcess_mc)==False or len(toProcess_mc)>1 or toProcess_data or toProcess_signal):
         print("Error: bTagEff should only be submitted with TTbar_dilepton")
         sys.exit(1)

   print(toProcess_mc)
   print(toProcess_data)
   print(toProcess_signal)
   
   # set dataBasePath depending on dCache options
   if args.scratchInput==False and args.copyDCache==False and args.condFileTransfer==False:
      print("Use dCache input and run on grid")
      dataBasePath = "-d"+get_dataBasePath_dCache(args.y)
   elif args.copyDCache:
      print("Copy input from dCache to condor lx-node and run on lx-node")
      dataBasePath = "-d$TMP/"
   elif args.condFileTransfer:
      print("Use Condor file transfer to copy from dCache to node")
      dataBasePath = ""
   else:
      print("Use scratch input")
      dataBasePath = ""
   
   # Ask if selected settings are correct
   runSubmit = True
   if disableConfirm==False:
      correctSamples = input("If you want to continue with the selected setting, enter 1:\n")
      if (correctSamples != 1):
         print("Abort Submission")
         runSubmit = False
   
   if runSubmit:
      # For single submit check if only one dataset is selected and then ask for file nr
      if (args.SingleSubmit and ((len(toProcess_mc)+len(toProcess_data)+len(toProcess_signal))==1)):
         singleFileNR = input("Please input the filenumber to submit:\n")
      elif args.SingleSubmit:
         print("SingleSubmit can only be used if one dataset is selected")
         exit(98)

      # not data processing if systematic shift is choosen
      if (toProcess_data and args.s not in nominalTypeSysts):
         toProcess_data=[]
         print("!!!!!!!!!!!!!!!!Data is not processed with systematic shift!!!!!!!!!!!!!!!!!!!!")
      
      # triggerEff only submitted with TTbar_dilepton and MET
      if (args.m == "triggerEff"):
         if ("TTbar_diLepton" not in toProcess_mc or len(toProcess_mc)>1 or toProcess_signal or "MET" not in toProcess_data):
            print("Error: triggerEff should only be submitted with TTbar_dilepton and MET")
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
                  print("Selected SingleFileNR is not available")
                  exit(98)
            
            for fileNR in range(startFileNR-1,totalFiles):     # loop over files for given dataset
                                 
               for logFile in glob.glob(logpath+"/"+args.m+"_"+x+"_"+str(fileNR+1)+"*"):  # Remove old logs
                  if os.path.exists(logFile):
                     os.remove(logFile)
               
               submitFile = logpath+"/"+args.m+"_"+x+"_"+str(fileNR+1)+".submit"    #define submit file
               
               if args.copyDCache:     #set path to dCache input in case of copying to condor node
                  inputPath = get_dataBasePath_dCache(args.y,False,True,True)+get_fileName(x,args.y,fileNR)
                  inputPath = inputPath.replace(" ", "")
               elif args.condFileTransfer:    #set path to dCache input in case of copying to condor node
                  inputPath = get_dataBasePath_dCache(args.y,True,False)+get_fileName(x,args.y,fileNR)
                  inputPath = inputPath.replace(" ", "")
               else:
                  inputPath = get_dataBasePath_dCache(args.y,False,True)+get_fileName(x,args.y,fileNR)
                  inputPath = inputPath.replace(" ", "")
                                 
               with open(submitFile,"w") as f:   # write condor submit
                  if args.copyDCache:
                     f.write("""
   Universe                = vanilla
   Executable              = run.sh
   Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {8} {9} {11} {12} {13}
   Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
   Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
   Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
   use_x509userproxy       = true
   Request_Memory          = {4} Mb
   Requirements            = (TARGET.Machine == "lxblade33.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade34.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade35.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade36.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade37.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade38.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade39.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade40.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade41.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade42.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade43.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade44.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade45.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade46.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade47.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade48.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade49.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade50.physik.rwth-aachen.de")
   Queue
      """.format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1),dataBasePath,inputPath,"\nRank = CpuFamily" if(x=="TTbar_diLepton") else "", getPath("cmsswBasePath"), getPath("frameworkBasePath"),getPath("dCacheBasePath")),)
                  elif args.condFileTransfer:
                     f.write("""
   Universe                = vanilla
   Executable              = run_condorFileTrans.sh
   Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {10} {11}
   Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
   Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
   Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
   stream_output           = true
   use_x509userproxy       = true
   should_transfer_files   = YES
   transfer_plugins        = curlsimple=curlSimple_plugin.py
   transfer_input_files    = {8}
   Request_Memory          = {4} Mb
   Requirements            = (TARGET.CpuFamily > 6) && (TARGET.Machine != "lxcip16.physik.rwth-aachen.de")  {9}
   Queue
      """.format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1),inputPath,"\nRank = CpuFamily" if(x=="TTbar_diLepton") else "", getPath("cmsswBasePath"), getPath("frameworkBasePath")),)
                  else:    # Grid submission
                     f.write("""
   universe                = grid
   grid_resource           = condor grid-{12}-rwth.gridka.de grid-{12}-rwth.gridka.de:9619
   Executable              = runGrid.sh
   Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {10} {11} {9}
   Log                     = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.log
   Output                  = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.out
   Error                   = logs/{5}/{6}/{0}/{1}/{1}_{3}_{7}.error
   use_x509userproxy       = true
   should_transfer_files   = YES
   transfer_output_files   = output_framework
   transfer_output_remaps  = "output_framework = {8}/{5}/{9}/output_framework"
   when_to_transfer_output = ON_SUCCESS
   success_exit_code      = 0
   getenv                  = yes
   Queue
      """.format(str(args.f),args.m,sampleStr,x,str(requ_mem),args.y,args.s,str(fileNR+1), getPath("scratchBasePath"),get_version(args.y),getPath("dCacheBasePath"),inputPath,checkCEjobs()))
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
      
      # get dCache path (currently 2016 is stored in Marius folder)
      #  ~if (args.y.find("2016")>=0):
         #  ~dCachePath = getPath("dCacheBasePath","teroerde")
      #  ~else:
         #  ~dCachePath = getPath("dCacheBasePath")
      
      dCachePath = getPath("dCacheBasePath")
                        
      with open(submitFile,"w") as f:   # write condor submit
         if args.copyDCache:
            f.write("""
      Universe                = vanilla
      Executable              = runTUnfold.sh
      Arguments               = -f{0} {1} {2} -s{3} {4} {5}
      Log                     = logs/{2}/{3}/{0}/{1}/{1}.log
      Output                  = logs/{2}/{3}/{0}/{1}/{1}.out
      Error                   = logs/{2}/{3}/{0}/{1}/{1}.error
      use_x509userproxy       = true
      RequestMemory           = 5000
      Requirements            = (TARGET.Machine == "lxblade33.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade34.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade35.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade36.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade37.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade38.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade39.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade40.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade41.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade42.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade43.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade44.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade45.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade46.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade47.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade48.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade49.physik.rwth-aachen.de") || (TARGET.Machine == "lxblade50.physik.rwth-aachen.de")
      Queue
      """.format(str(args.f),args.m,args.y,args.s,get_version(args.y),dCachePath),)
         else:    #grid submission
            f.write("""
      universe                = grid
      grid_resource           = condor grid-{12}-rwth.gridka.de grid-{12}-rwth.gridka.de:9619
      Executable              = runGrid.sh
      Arguments               = -f{0} {1} {2} {5} -s{6} --fileNR={7} {10} {11} {9} {13}
      Log                     = logs/{5}/{6}/{0}/{1}/{1}.log
      Output                  = logs/{5}/{6}/{0}/{1}/{1}.out
      Error                   = logs/{5}/{6}/{0}/{1}/{1}.error
      use_x509userproxy       = true
      should_transfer_files   = YES
      transfer_output_files   = output_framework
      transfer_output_remaps  = "output_framework = {8}/{5}/{9}/output_framework"
      when_to_transfer_output = ON_SUCCESS
      success_exit_code       = 0
      RequestMemory           = 10000
      getenv                  = yes
      Queue
         """.format(str(args.f),args.m,"placeholder","placeholder","placeholder",args.y,args.s,"placeholder", getPath("scratchBasePath"),get_version(args.y),getPath("dCacheBasePath"),"placeholder",checkCEjobs(),getPath("gridname")))
      subprocess.call(["condor_submit", submitFile])
   
   

if __name__ == "__main__":
   
   #############################################
   # Select datasets to process
   #############################################
   #  ~toProcess_mc=["TTbar_diLepton","TTbar_amcatnlo","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"]
   #  ~toProcess_mc=["ZZ"]
   #  ~toProcess_mc=["DrellYan_NLO"]
   #  ~toProcess_mc=["TTbar_diLepton"]
   #  ~toProcess_mc=["TTbar_diLepton","TTbar_diLepton_tau","TTbar_hadronic","TTbar_singleLepton"]
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
   #  ~toProcess_mc=["TTbar_diLepton_tau_MTOP175p5"]
   #  ~toProcess_mc=["DrellYan_M10to50"]
   toProcess_mc=["DrellYan"]
   #  ~toProcess_mc=["SingleTop","SingleTop_DS"]
   #  ~toProcess_mc=["bb4l"]
   #  ~toProcess_mc=["bb4l_new"]
   #  ~toProcess_mc=[]
   #  ~toProcess_mc=allMC
   
   #  ~toProcess_data=["DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"]
   #  ~toProcess_data=["DoubleMuon","MuonEG","SingleMuon","EGamma"]      #2018
   #  ~toProcess_data=["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"]       #2017, 2016
   #  ~toProcess_data=["MET"]      
   toProcess_data=[]
         
   #  ~toProcess_signal=["T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"]
   #  ~toProcess_signal=["T1tttt_1200_800","T1tttt_1500_100","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10"]
   #  ~toProcess_signal=["T1tttt_1200_800"]
   #  ~toProcess_signal=["T2tt_525_438","T2tt_525_350"]
   toProcess_signal=[]

   # define command line arguments
   parser = argparse.ArgumentParser()
   parser.add_argument('-m', type=str, default="distributions", help="module, default is distributions")
   parser.add_argument('-f', type=float, choices=[Range(0.0, 1.0)], default=1.0, help="process fraction")
   parser.add_argument('-y', type=str, help="year to be set as ANALYSIS_YEAR_CONFIG (e.g. 2016_preVFP)",required=True)
   parser.add_argument('-s', type=str, default="Nominal", help="systematic shift")
   parser.add_argument('--scratchInput', action='store_true', default=False, help="Use nTuple stored on scratch, otherwise dCache Input is used.")
   parser.add_argument('--copyDCache', action='store_true', default=False, help="Copy nTuples stored on dCache to lx-node before running the code on lx-node.")
   parser.add_argument('--condFileTransfer', action='store_true', default=False, help="Use condor file transfer to copy from dCache to node before running the code.")
   parser.add_argument('--SingleSubmit', action='store_true' )
   parser.add_argument('--bTagEff_complete', action='store_true', default=False, help="Submits bTagEff jobs with all relevant systematics (use with care!)")
   parser.add_argument('--distributions_complete', action='store_true', default=False, help="Submits distributions jobs with all relevant systematics (use with care!)")
   parser.add_argument('--tunfold_complete', action='store_true', default=False, help="Submits tunfold jobs with all relevant systematics (use with care!)")
   parser.add_argument('--tunfold_pdf_complete', action='store_true', default=False, help="Submits tunfold jobs with all pdf shifts (use with care!)")
   parser.add_argument('--pdf_complete', action='store_true', default=False, help="Submits distributions jobs with all pdf shifts (use with care!)")
   parser.add_argument('--noConfirmation', action='store_true', default=False, help="Disables keyboard input befor submission")

   args = parser.parse_args()
      
   # Upload CMSSW and FW to dCache if running on grid
   if args.scratchInput==False and args.copyDCache==False and args.condFileTransfer==False:
      uploadCompressedCMSSW()
      uploadCompressedFW()
   
   # Remove HEM syst uncertainty for all years except 2018
   removeJES_HEM(args.y)
   
   # Set correct module for tunfold_complete
   if (args.tunfold_complete or args.tunfold_pdf_complete):
      args.m = "TUnfold_binning"
      
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
         print("bTagEff_complete can only be used if bTagEff is selected as module!")
         exit(98)
   elif (args.distributions_complete):
      if (args.m == "" or args.m == "distributions"):
         print("Submit nominal")
         if (args.y=="2018"):
            submit(args,allMC,allData2018,[])
         elif (args.y=="2017"):
            submit(args,allMC,allData2017,[])
         elif (args.y=="2016_preVFP" or args.y == "2016_postVFP"):
            submit(args,allMC,allData2016,[])
         else:
            print(args.y+" does not match correct year")
         for syst in sample_allSyst_dict.keys():
            args.s = syst
            submit(args,sample_allSyst_dict[syst],[],[],args.noConfirmation)
      else:
         print("distributions_complete can only be used if distributions is selected as module!")
         exit(98)
   elif (args.pdf_complete):
      if (args.m == "" or args.m == "distributions"):
         for varNumber in xrange(1,51):
            for var in ["UP","DOWN"]:
               args.s = "PDF_{}_{}".format(varNumber,var)
               submit(args,["TTbar_diLepton"],[],[],args.noConfirmation)
      else:
         print("pdf_complete can only be used if distributions is selected as module!")
         exit(98)
      
   
