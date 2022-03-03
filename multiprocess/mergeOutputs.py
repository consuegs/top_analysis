#!/usr/bin/env python2
import process
import condor_status
import argparse
import subprocess as sp
import re
import configparser
import glob
import os
import sys
import shutil
import shlex
import ROOT
from ROOT import *
import utilities
sys.path.append("../users")
from getPath import getPath


def getNTupleVersion(year):    # checks config for number of files for given dataset
   config = configparser.ConfigParser()
   config.read("../config"+year+".ini")
   return config["input"]["version"]

def getInfoFromLogPath(logPath):
    info = {}
    info["year"] = logPath.split("/")[1].replace("/","")
    info["syst"] = logPath.split("/")[2].replace("/","")
    info["fraction"] = logPath.split("/")[3].replace("/","")
    info["module"] = logPath.split("/")[4].replace("/","")
    return info

def isComplete(logPath,dataset):    # checks if files in finished.txt match expected files defined in config[].ini
    file = open(logPath+"finished.txt")
    data = file.read()
    pattern = dataset+"_\d+"
    producedFiles = len(re.findall(pattern,data))
    expectedFiles = process.get_fileNR(dataset,getInfoFromLogPath(logPath)["year"])
    
    return producedFiles == expectedFiles

def getSplitSamples(dataset,logPath):   # extracts sample names from finished.txt
    file = open(logPath+"finished.txt")
    data = file.read()
    pattern = dataset+"_\d+"
    return re.findall(pattern,data)

def listAllObjectPaths(file_,output,path=""):
    file_.cd(path)
    for key in gDirectory.GetListOfKeys():
        if key.IsFolder():
            listAllObjectPaths(file_,output,path+"/"+key.GetName())
        else:
            output.append(path+"/"+key.GetName())

def mergeTree(dataset,logPath,treePath):    # merges minTrees of given dataset
    splitSamples = []
    outfile = dataset+"_merged.root"
    for splitSample in getSplitSamples(dataset,logPath):
        splitSamples.append(splitSample+".root")
    runDir = os.getcwd()
    os.chdir(treePath)
    if os.path.exists(outfile):     # remove old output if existing
        os.remove(outfile)
    for splitSample in splitSamples:
        sp.call(["rootrm",splitSample+":ttbar_res"])    # remove buggy parts of output (probably feature of older root version)
    if sp.call(["hadd","-f",outfile]+splitSamples):     # merging command
        sys.exit(1)
    sp.call(["rootrm",outfile+":ttbar_res100.0/ttbar_res"])
    os.chdir(runDir)

def mergeHist(dataset,logPath,histPath):    # merges histograms of given dataset
    nTupleVersion = getNTupleVersion(getInfoFromLogPath(logPath)["year"])
    splitSamples = []
    fileStart = "histograms_" if getInfoFromLogPath(logPath)["module"] == "distributions" else ""
    outfile = fileStart+dataset+"_merged_"+nTupleVersion+".root"
    for splitSample in getSplitSamples(dataset,logPath):
        splitSamples.append(fileStart+splitSample+"_"+nTupleVersion+".root")
    runDir = os.getcwd()
    os.chdir(histPath)
    if os.path.exists(outfile):
        os.remove(outfile)
    if sp.call(["hadd","-f",outfile]+splitSamples,stdout=open(os.devnull, 'wb')):
        sys.exit(1)
    print "Created",outfile
    os.chdir(runDir)
    
def mergeHists_forSamples(logPath,histPath,datasetList,samplesToMerge,outputName):        # merges histograms for different samples (mainly needed for datacard production)
    
    for sample in samplesToMerge:
        if sample not in datasetList:
            print "No merging for:",samplesToMerge
            return
    
    nTupleVersion = getNTupleVersion(getInfoFromLogPath(logPath)["year"])
    fileStart = "histograms_" if getInfoFromLogPath(logPath)["module"] == "distributions" else ""
    outfile = fileStart+outputName+"_merged_"+nTupleVersion+".root"
    inputSamples = [fileStart+x+"_merged_"+nTupleVersion+".root" for x in samplesToMerge]
    tempSamples = []
    
    if os.path.exists(histPath+"/temp/") == False:      #create temp folder for files with renamed hists
        os.mkdir(histPath+"/temp/")
    runDir = os.getcwd()
    os.chdir(histPath)
    for file_ in inputSamples:  #loop over inputSamples
        f = TFile(file_,"READ")
        hists = {}
        histNames = []
        listAllObjectPaths(f,histNames)     # get all paths in root file
        for histName in histNames:
            hist = f.Get(histName)
            hists[histName] = hist
        f_out = TFile(histPath+"/temp/"+file_,"RECREATE")
        tempSamples.append(histPath+"/temp/"+file_)
        f_out.cd()
        for histName in hists.keys():
            folderName = histName.rsplit("/",1)[0]
            folderName = folderName[1:] 
            with utilities.Quiet(kError + 1):   # avoid spaming from root
                if f_out.cd(folderName) == False:
                    f_out.mkdir(folderName)
                    f_out.cd(folderName)
            hists[histName].Write(outputName)
        f_out.Close()
        f.Close()
    
    if os.path.exists(outfile):
        os.remove(outfile)
    if sp.call(["hadd","-f",outfile]+tempSamples,stdout=open(os.devnull, 'wb')):    # combine temp hists with hadd
        sys.exit(1)
    if os.path.exists(histPath+"/temp/"):   # cleaning
        shutil.rmtree(histPath+"/temp/")
    
    datasetList.append(outputName)      # add merged samples to datasetList
    print "Created",outfile
    os.chdir(runDir)

def mergeAllHists(datasets,logPath,histPath):       # merge all histograms into one file used for plotting scripts
    nTupleVersion = getNTupleVersion(getInfoFromLogPath(logPath)["year"])
    datasetFiles = []
    outfile = "histograms_merged_"+nTupleVersion+".root"
    fileStart = "histograms_" if getInfoFromLogPath(logPath)["module"] == "distributions" else ""
    for dataset in datasets:
        datasetFiles.append(fileStart+dataset+"_merged_"+nTupleVersion+".root")
    runDir = os.getcwd()
    os.chdir(histPath)
    for datasetFile in datasetFiles:
        if os.path.exists(datasetFile)==False:
            print "Error:", datasetFile, "is missing"
            sys.exit(98)
    if os.path.exists(outfile):
        os.remove(outfile)
    if sp.call(["hadd","-f",outfile]+datasetFiles,stdout=open(os.devnull, 'wb')):
        sys.exit(1)
    print "Created",outfile
    os.chdir(runDir)

def getDatasetList(logPath):    # get list of dataset, where jobs have been submitted
    datasets = []
    for out in glob.glob(logPath+"*.out"):
        fileName = condor_status.getNameFromOutFile(out)
        datasetName = fileName.split("_")
        datasetName.pop()
        datasetName = "_".join(datasetName)
        if datasetName not in datasets:
            datasets.append(datasetName)
    return datasets
    
def getTreePath(logPath):       # get correct treePath from logPath
    info = getInfoFromLogPath(logPath)
    treePath = "{}/{}/{}/minTrees/{}/{}/".format(getPath("scratchBasePath"),info["year"],getNTupleVersion(info["year"]),str(float(info["fraction"])*100),info["syst"])
    return treePath

def getHistPath(logPath):       # get correct histPath from logPath
    info = getInfoFromLogPath(logPath)
    treePath = "{}/{}/{}/output_framework/multiHists/{}/".format(getPath("scratchBasePath"),info["year"],getNTupleVersion(info["year"]),info["syst"])
    return treePath
    
def getHistPath_module(logPath):       # get correct histPath from logPath for modules other than distributions
    info = getInfoFromLogPath(logPath)
    treePath = "{}/{}/{}/output_framework/{}/{}/".format(getPath("scratchBasePath"),info["year"],getNTupleVersion(info["year"]),info["module"],info["syst"])
    return treePath
    
def mergeRequired(logPath,isDistributions=True,useTrees=True):       # check if there are more recent individual outputs than the merged one
    if useTrees and isDistributions:
        fileList = glob.glob(getTreePath(logPath)+"/*")
    elif isDistributions:
        fileList = glob.glob(getHistPath(logPath)+"/*")
    else:
        fileList = glob.glob(getHistPath_module(logPath)+"/*")
    latestPath = max(fileList, key=os.path.getctime)
    if latestPath.find("merged")>0:
        return False
    else:
        return True

def getSytNameCombine(systName):    # correct syst name for combine (needs "Up" and "Down")
    
    # fix naming for syst with own samples
    if systName=="Nominal" :
        return ""
    elif systName.find("UP")>0 :
        return "_"+systName.replace("_UP","Up")
    elif systName.find("DOWN")>0 :
        return "_"+systName.replace("_DOWN","Down")
    elif systName in ["CR1","CR2","ERDON","MTOP169p5","MTOP175p5"]:
        return ""
    elif systName=="TOP_PT":
        return "_"+systName+"Up"
    else:
        return "_"+systName
        

def getSampleNameCombine(sampleName):       # needed to correctly handle SampleNames of alternative Samples for combine
    if sampleName.find("UETUNE")>0 :
        if sampleName.find("UP")>0 :
            return sampleName.replace("_UETUNE_UP","")
        else:
            return sampleName.replace("_UETUNE_DOWN","")
    elif sampleName.find("MATCH")>0 :
        if sampleName.find("UP")>0 :
            return sampleName.replace("_MATCH_UP","")
        else:
            return sampleName.replace("_MATCH_DOWN","")
    elif sampleName.find("MTOP169p5")>0 :
        return sampleName.replace("MTOP169p5","")+"MTOPDown"
    elif sampleName.find("MTOP175p5")>0 :
        return sampleName.replace("MTOP175p5","")+"MTOPUp"
    elif any(x in sampleName for x in["CR1","CR2","ERDON"]):    # one sided uncertainties have Up=Shift and Down=Nominal
        return sampleName+"Up"
    else:
        return sampleName

def mergeForCombine(logPath,histPath,sampleList):
    info = getInfoFromLogPath(logPath)
    nTupleVersion = getNTupleVersion(info["year"])
    #  ~outputDir = histPath+"../combine/"
    outputDir = histPath+"../combine_new/"
    outfile = "{}combineInput_{}.root".format(outputDir,nTupleVersion)
    
    filesToMerge = []    # keep track of files to merge
    
    if os.path.exists(outputDir) == False:      #create combine folder for output files
        os.mkdir(outputDir)
    
    for systPath in glob.glob("logs/"+info["year"]+"/*"):   # loop for renaming and merging histograms
        distrLogPath = systPath+"/1.0/distributions/"
        if(os.path.exists(distrLogPath)):
            info_syst = getInfoFromLogPath(distrLogPath)
            print "-------Renaming {}---------------".format(info_syst["syst"])
            mergedHistPath = getHistPath(distrLogPath)+"histograms_merged_{}.root".format(nTupleVersion)
            f = TFile(mergedHistPath,"READ")
            hists = {}
            histNames = []
            listAllObjectPaths(f,histNames)     # get all paths in root file
            for histName in histNames:
                histNameSplitted = histName.split("/")
                if histNameSplitted[1] != "distributions100.0" or histNameSplitted[2] != "baseline":     # only store relevant histograms with full stat.
                    continue
                sampleName = histName.split("/")[-1]
                strippedSampleName = sampleName.replace("_"+info_syst["syst"],"")   # only store processes needed for combine
                if strippedSampleName not in sampleList:
                    continue
                else:
                    hist = f.Get(histName)
                    hists[histName] = hist
            
            tempFileNames = {}
            tempFileNames[info_syst["syst"]] = outputDir+info_syst["syst"]+".root"
            if (info_syst["syst"]=="Nominal"):      # one sided syst need Down=Nominal
                tempFileNames[info_syst["syst"]] = outputDir+info_syst["syst"]+".root"
                for syst in ["CR1","CR2","ERDON","TOP_PT"]:
                    tempFileNames[syst] = outputDir+syst+"_DOWN.root"
            
            for syst,tempFileName in tempFileNames.items():
                if (info_syst["syst"]!=syst):   # for one sided down shift
                    systNameCombine = "_"+syst+"Down"
                else:
                    systNameCombine = getSytNameCombine(info_syst["syst"])  # get correct name of systematic for combine
                f_out = TFile(tempFileName,"RECREATE")
                filesToMerge.append(tempFileName)
                f_out.cd()
                for histName in hists.keys():
                    folderName = histName.rsplit("/",1)[0]
                    folderName = folderName[1:]
                    sampleName = histName.split("/")[-1]
                    outputHistName = "{}{}".format(getSampleNameCombine(sampleName),systNameCombine)     # add syst name to histName
                    
                    with utilities.Quiet(kError + 1):   # avoid spaming from root
                        if f_out.cd(folderName) == False:
                            f_out.mkdir(folderName)
                            f_out.cd(folderName)
                    hists[histName].Write(outputHistName)
                f_out.Close()
            
            f.Close()
    
    print "-------Merging renamed histograms---------------"
    if os.path.exists(outfile):
        os.remove(outfile)
    if sp.call(["hadd","-f",outfile]+filesToMerge,stdout=open(os.devnull, 'wb')):    # combine temp hists with hadd
        sys.exit(1)
    if sp.call(["chmod","a+rx",outfile],stdout=open(os.devnull, 'wb')):     # set rights for other to excess output file
        sys.exit(1)
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("logPath", type=str, nargs='?', default="logs/2018/Nominal/1.0/distributions/", help='<Required> Path to log directory of jobs to merge')
    parser.add_argument("--onlyTrees", action='store_true', help="Only minTrees are merged")
    parser.add_argument("--onlyHists", action='store_true', help="Only histograms are merged")
    parser.add_argument('--singleDataset', type=str, default="", help="Defines single dataset to merge (otherwise all samples are merged)")
    parser.add_argument("--mergeAllHists", action='store_true', default=True, help="Merge all hist of different datasets to one file")
    parser.add_argument("--mergeForCombine", action='store_true', default=False, help="Merge all hist of different datasets and systematics for combine input")
    parser.add_argument('-y', type=str, help="Year to be merged, currently only used in combination with --mergeForCombine",default="2018")
    args = parser.parse_args()
    
    isDistributions = (args.logPath.find("distributions") > 0)
    isTriggerEff = (args.logPath.find("triggerEff") > 0)
        
    if (isDistributions == False and (args.onlyHists == False)):
        print "Error: For other modules than distributions Trees can not be merged!"
        print "Continuing with --onlyHists option"
        args.onlyHists = True
    
    treePath = getTreePath(args.logPath)
    histPath = getHistPath(args.logPath)
        
    if(isDistributions == False):
        histPath = getHistPath_module(args.logPath)
    
    if args.singleDataset == "":
        datasetList = getDatasetList(args.logPath)
    else:
        if args.mergeAllHists or args.onlyTrees:
            print "Error: mergeAllHists option not compatible with singeSample or onlyTrees"
        if args.singleDataset in getDatasetList(args.logPath):
            datasetList = [args.singleDataset]
        else:
            print "Error: Selected sample not in logPath"
            sys.exit(1)
        
    for dataset in datasetList:     # check if all datasets, which are intended for merging, are finished
            if isComplete(args.logPath,dataset) == False:
                print "Error:",dataset,"is not complete"
                sys.exit(1)
    
    if (args.mergeForCombine==False):
        for dataset in datasetList:
                if isComplete(args.logPath,dataset):
                    print "-------------Start merging", dataset, "----------------------"
                    if args.onlyHists == False:
                        mergeTree(dataset,args.logPath,treePath)
                    if args.onlyTrees == False:
                        mergeHist(dataset,args.logPath,histPath)
        
        
        mergeHists_forSamples(args.logPath,histPath,datasetList,["DrellYan_M10to50","DrellYan_NLO"],"DrellYan_comb")    # merge samples for datacards
        mergeHists_forSamples(args.logPath,histPath,datasetList,["TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"],"TTbar_other")
        mergeHists_forSamples(args.logPath,histPath,datasetList,["WJetsToLNu","WW","WZ","ZZ","ttZ_QQ","ttZ_2L","ttW"],"otherBKG")
        mergeHists_forSamples(args.logPath,histPath,datasetList,["DoubleMuon","EGamma","MuonEG","SingleMuon"],"data_obs") # 2018
        mergeHists_forSamples(args.logPath,histPath,datasetList,["DoubleMuon","MuonEG","SingleMuon","DoubleEG","SingleElectron"],"data_obs") #2016/2017
        
        altSampleNames = ["CR1","CR2","ERDON","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP169p5","MTOP175p5"]    # handle alt. sample for datacards
        if any(x in args.logPath for x in altSampleNames):
            info = getInfoFromLogPath(args.logPath)
            mergeHists_forSamples(args.logPath,histPath,datasetList,["TTbar_diLepton_tau_"+info["syst"],"TTbar_singleLepton_"+info["syst"],"TTbar_hadronic_"+info["syst"]],"TTbar_other_"+info["syst"])
        
        if args.mergeAllHists and (isDistributions or isTriggerEff):
            print "-------------Start merging all datasets to one histFile----------------------"
            mergeAllHists(datasetList,args.logPath,histPath)
    elif(args.mergeForCombine):
        args.logPath = "logs/{0}/Nominal/1.0/distributions/".format(args.y)
        histPath = getHistPath(args.logPath)
        mergeForCombine(args.logPath,histPath,["TTbar_diLepton","DrellYan_comb","TTbar_other","otherBKG","SingleTop","data_obs"])
