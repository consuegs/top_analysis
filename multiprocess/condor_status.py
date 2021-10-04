#!/usr/bin/env python
import argparse
import subprocess
import time
import os
import glob
from termcolor import colored
import re

def getInfos():
    out = subprocess.check_output(["condor_q", "-long"])
    for jobStrings in out.split("\n\n"):
        for line in jobStrings.split("\n"):
            if line:
                break
                print [line.split(" = ")]
    return [ dict([line.replace("\"","").split(" = ") for line in jobStrings.split("\n") if " = " in line]) for jobStrings in out.split("\n\n") if jobStrings ]

def getSummary():
    out = subprocess.check_output(["condor_q"])
    return out.split("\n")[-3]

def getNameFromOutFile(outFile):
    return outFile.split("/")[-1].replace("distributions_","").replace(".out","")

def getNameFromFile(fname):
    return fname.split("dataset=")[-1]
    
def getMachineFromOut(outName):
    logName = outName.replace(".out",".log")
    with open(logName,"r") as f:
        next(f)
        for line in f:
                if line.find("alias") >= 0:
                    machine = line.split("alias=")[-1]
                    machine = machine.split("&")[0]
    return machine

def getProgressFromOut(outName,checkCompleted=False):
    processing = False
    if os.path.exists(outName):
        lines = len(open(outName).readlines(  ))
        if lines !=0:
            with open(outName,"r") as f:
                next(f)
                for line in f:
                    if line.find("Processing")==0:
                        processing = True
                        progress = line.count(".")
                        color = "yellow"
                        if progress==10:
                            color = "green"
                        print colored("--"+line.split(".")[0]+" "+str(progress*10)+"%",color)
        if(processing == False and checkCompleted and lines>2):
            color = "red"
            print colored(getNameFromOutFile(outName)+" failed",color)
        

def getStatusFromOut(outName,running):
    finished = False
    if os.path.exists(outName) == False or len(open(outName).readlines(  ))<=2:
        if running:
            color = "blue"
            print colored(getNameFromOutFile(outName)+" no output yet",color)
        else:
            color = "red"
            print colored(getNameFromOutFile(outName)+" failed",color) 
    else:
        datasetName = getNameFromOutFile(outName)
        if len(open(outName).readlines(  )) !=0:
            with open(outName,"r") as f:
                next(f)
                for line in f:
                    if line.find("took")>=0:
                        finished = True
                        color = "green"
                        print colored(datasetName+" is finished, "+line.split(")")[-1].replace("=","").replace("\n","")+" on "+getMachineFromOut(outName),color)
        if finished == False and running:
            getProgressFromOut(outName,True) 
        elif finished == False and running == False:
            color = "red"
            print colored(getNameFromOutFile(outName)+" failed"+" on "+getMachineFromOut(outName),color)
    
    return finished
    
def checkStatusFromLog(logPath,runningLogs):
    with open(args.checkCompleted+"finished.txt","w") as f:
        for log in glob.glob(args.checkCompleted+"*.log"):
            outFile = log.replace(".log",".out")
            if outFile in runningLogs:
                getStatusFromOut(outFile,True)
            else:
                if getStatusFromOut(outFile,False):
                    f.write(getNameFromOutFile(outFile)+"\n")
    return args.checkCompleted+"finished.txt"

def checkStatusFromQueue(printOutput=True):
    jobs = getInfos()
    jobs = sorted(jobs, key=lambda l: l["JobStatus"]+l["ClusterId"])
    susJobs = []
    runningLogs = []

    for job in jobs:
        if job["Cmd"].split("/")[-1]=="run.sh":
            name = getNameFromFile(job["Args"])
        else:
            name = (job["Args"].split(",")[0]).split("/")[2]
        name += " "*max([0,(40-len(name))])
        runningLogs.append(job["Out"])
        if(printOutput == False): continue     # show job status only in nominal mode
        jStatus = job["JobStatus"]
        if jStatus == "1":
            print colored("\033[1m"+name+"    "+job["ClusterId"]+"       idle"+"\033[0m","yellow")
        elif jStatus == "2":
            print colored("\033[1m"+name+"    "+job["RemoteHost"]+"       running"+"\033[0m","green")
            getProgressFromOut(job["Out"])
        elif jStatus == "5":
            print colored("\033[1m"+name+"    "+job["ClusterId"]+"       held"+"\033[0m","red")
        elif jStatus == "7":
            susTime = (time.time()-int(job["LastSuspensionTime"]))/60.
            print colored("\033[1m"+name+"    "+job["RemoteHost"]+"       suspended since {:.2f} min".format(susTime)+"         "+job["RemoteHost"]+"\033[0m","red")
            getProgressFromOut(job["Out"])
            if susTime > 10 :
                value = input("Job suspended for more than 10 Minutes, please enter 1 to kill job and start again or 0 to continue:\n")
                if value==1:
                    print "resubmitting job"
                    susJobs.append(job["ClusterId"])
        else:
            print "job status = ", jStatus
    
    if(printOutput):
        print getSummary()

        for job in susJobs:
            os.system("condor_hold "+job)
        if len(susJobs)!=0: os.system("condor_release dmeuser")
    
    return runningLogs
    
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--checkCompleted', type=str, default="", help="Checks status of all jobs in given log folder")
    args = parser.parse_args()
    
    runningLogs = checkStatusFromQueue(args.checkCompleted == "")
    
    if(args.checkCompleted != ""):
        if(os.path.exists(args.checkCompleted)):      # check status from log (does also inlcude finished jobs) and write finished names to list
            checkStatusFromLog(args.checkCompleted,runningLogs)
        else:
            print "Wrong Path"

