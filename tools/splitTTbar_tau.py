#!/usr/bin/env python2
import argparse
import os
import ROOT
import multiprocessing
import tqdm
import re
import glob
from functools import partial
import subprocess as sp
import configparser


def getTargets(inputFile):
    datasetName = re.findall("TTTo2\w+_13TeV",inputFile)
    datasetName = datasetName[0].replace("_13TeV","")
    target = inputFile.replace("nTuple","nTuple/tau_splitted/"+datasetName)
    target = os.path.join(target.replace(".root", "_tau.root"))
    targetFolder = os.path.dirname(target)
    return target,targetFolder
    

def createTauOutput(inputFile):
    print "Started processing", inputFile
    
    outputName,directoryPath = getTargets(inputFile)
    
    f = ROOT.TFile(inputFile)
    hists = []
    for key in f.GetDirectory("TreeWriter").GetListOfKeys():
        h = key.ReadObj()
        if not isinstance(h, ROOT.TH1): continue
        h.SetDirectory(0)
        hists.append(h)

    ch = ROOT.TChain("TreeWriter/eventTree")
    ch.AddFile(inputFile)
    
    fout = ROOT.TFile(outputName, "recreate")
    fout.mkdir("TreeWriter")
    fout.cd("TreeWriter")
    new = ch.CopyTree("ttbarDecayMode > 3")
    new.Write("", ROOT.TObject.kWriteDelete)
    for h in hists: h.Write()
    fout.Close()
    print "Created", outputName

def getMergingScheme(localfiles, maxMergeSize=10):
    scheme = {}
    currentSize = 0
    currentList = []
    for filePath in localfiles:
        fileSize = os.path.getsize(filePath)
        currentSize += fileSize
        if((currentSize*1e-9)<maxMergeSize):
           currentList.append(filePath)
        else:
            scheme[len(scheme)+1]=currentList
            currentList = [filePath]
            currentSize = fileSize
    scheme[len(scheme)+1]=currentList
    return scheme
    
def getAllSamples(year):
    config = configparser.ConfigParser()
    config.read("../config"+year+".ini")
    dataBasePath = config["input"]["dataBasePath"]+config["input"]["version"]
    samples = []
    for sampleFile in glob.glob(dataBasePath+"/nTuple/*"):
        if sampleFile.find(".root")>0:
            continue
        if sampleFile.find("TTTo2L2Nu")>0:
            samples.append(sampleFile)
    return samples
    

def getInputFiles(sample):
    if sample.find(".root")>0:
        print "Error: Please enter sample name, not file name (without .root and _1)"
        exit(98)
    files = []
    for f in glob.glob(sample+"*root"):
        if f.find("tau") == -1:
            files.append(f)
    
    return files
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #  ~parser.add_argument("inputFiles")
    #  ~parser.add_argument("inputFiles", nargs='+', help='<Required> List of input files')
    #  ~parser.add_argument("inputSamples", nargs='+', help='<Required> List of input samples')
    parser.add_argument("year", type=str, help='<Required> Year for which samples should be processed')
    parser.add_argument("--singleSample", type=str, default="", help='Option to run on a single sample, instead of all samples of a year')
    args = parser.parse_args()
        
    if args.singleSample == "":
        inputSamples = getAllSamples(args.year)
    else:
        inputSamples = [args.singleSample]
    
    i = 0
    for sample in inputSamples:
        i+=1
        print "-----------------Starting with sample({}/{}) ".format(i,len(inputSamples))+sample.split("/")[-1]+"-------------------------"
        
        inputFiles = getInputFiles(sample)
        targetFolder = getTargets(inputFiles[0])[1]
        
        if os.path.exists(targetFolder):
            runAnyway = input("Already found tau samples connected to {}\n If you still want to continue, input 1 :\n".format(sample))
            if runAnyway!= "1":
                continue
            
        # Create target folder
        if not os.path.exists(targetFolder):
            print "create "+targetFolder
            os.makedirs(targetFolder)
            
        # Separate tau events
        p = multiprocessing.Pool(7)
        for _ in tqdm.tqdm(p.imap_unordered(partial(createTauOutput),inputFiles), total=len(inputFiles)):
            pass
        
        # Merge to given size
        localFiles = [f for f in glob.glob("{}/*root".format(targetFolder))]
        mergingScheme = getMergingScheme(localFiles)
        for outputNr in mergingScheme:
            outputFile = inputFiles[0].split("_1.root")[0]
            if sp.call(["hadd","-f",outputFile+"_tau_"+str(outputNr)+".root"]+mergingScheme[outputNr]):
                sys.exit(1)
    



