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
import shlex

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
    info = getInfoFromLogPath(args.logPath)
    treePath = "/net/data_cms1b/user/dmeuser/top_analysis/{}/{}/minTrees/{}/{}/".format(info["year"],getNTupleVersion(info["year"]),str(float(info["fraction"])*100),info["syst"])
    return treePath

def getHistPath(logPath):       # get correct histPath from logPath
    info = getInfoFromLogPath(args.logPath)
    treePath = "/net/data_cms1b/user/dmeuser/top_analysis/{}/{}/output_framework/multiHists/{}/".format(info["year"],getNTupleVersion(info["year"]),info["syst"])
    return treePath
    
def getHistPath_module(logPath):       # get correct histPath from logPath for modules other than distributions
    info = getInfoFromLogPath(args.logPath)
    treePath = "/net/data_cms1b/user/dmeuser/top_analysis/{}/{}/output_framework/{}/{}/".format(info["year"],getNTupleVersion(info["year"]),info["module"],info["syst"])
    return treePath

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("logPath", type=str, help='<Required> Path to log directory of jobs to merge')
    parser.add_argument("--onlyTrees", action='store_true', help="Only minTrees are merged")
    parser.add_argument("--onlyHists", action='store_true', help="Only histograms are merged")
    parser.add_argument('--singleDataset', type=str, default="", help="Defines single dataset to merge (otherwise all samples are merged)")
    parser.add_argument("--mergeAllHists", action='store_true', default=True, help="Merge all hist of different datasets to one file")
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
    
    for dataset in datasetList:
            if isComplete(args.logPath,dataset):
                print "-------------Start merging", dataset, "----------------------"
                if args.onlyHists == False:
                    mergeTree(dataset,args.logPath,treePath)
                if args.onlyTrees == False:
                    mergeHist(dataset,args.logPath,histPath)
    
    if args.mergeAllHists and (isDistributions or isTriggerEff):
        print "-------------Start merging all datasets to one histFile----------------------"
        mergeAllHists(datasetList,args.logPath,histPath)
