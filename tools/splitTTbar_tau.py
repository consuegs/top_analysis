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

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #  ~parser.add_argument("inputFiles")
    parser.add_argument("inputFiles", nargs='+', help='<Required> List of input files')
    args = parser.parse_args()
    
    targetFolder = getTargets(args.inputFiles[0])[1]
    
    # Create target folder
    if not os.path.exists(targetFolder):
        print "create "+targetFolder
        os.makedirs(targetFolder)
        
    # Separate tau events
    p = multiprocessing.Pool(7)
    for _ in tqdm.tqdm(p.imap_unordered(partial(createTauOutput),args.inputFiles), total=len(args.inputFiles)):
        pass
    
    # Merge to given size
    localFiles = [f for f in glob.glob("{}/*root".format(targetFolder))]
    mergingScheme = getMergingScheme(localFiles)
    for outputNr in mergingScheme:
        outputFile = args.inputFiles[0].split("_1.root")[0]
        if sp.call(["hadd","-f",outputFile+"_tau_"+str(outputNr)+".root"]+mergingScheme[outputNr]):
            sys.exit(1)
    



