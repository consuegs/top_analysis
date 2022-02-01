#!/usr/bin/env python
import argparse
import subprocess
import time
import os
import glob
from termcolor import colored
import re
import subprocess as sp
import mergeOutputs


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
    moduleName = outFile.split("/")[-2]
    return outFile.split("/")[-1].replace(moduleName+"_","").replace(".out","")

def getSystFromOutFile(outFile):
    systName = outFile.split("/")[2]
    return systName

def getNameFromFile(fname):
    return fname.split(" ")[1]+" "+(fname.split("dataset=")[-1]).split("/")[0]
    
def getMachineFromOut(outName):
    logName = outName.replace(".out",".log")
    with open(logName,"r") as f:
        for line in f:
            if line.find("alias") >= 0:
                machine = line.split("alias=")[-1]
                machine = machine.split("&")[0]
    return machine
    
def getCopyTimeFromError(outName):
    logName = outName.replace(".out",".error")
    with open(logName,"r") as f:
        for line in f:
            if line.find("seconds") >= 0:
                time = line.split("in")[-1]
                time = time.split("seconds")[0]
    return time
    
def resubmitJob(outName,automaticResubmit=False):
    if automaticResubmit: value = 1
    else :
        value = input("For resubmitting "+getNameFromOutFile(outName)+" enter 1. For showing .out, .log and .error enter 2:\n")
    if value==1:
        print "Resubmitting..."
        sp.call(["condor_submit", outName.replace(".out",".submit")])
        sp.call(["rm", outName])
        sp.call(["rm", outName.replace(".out",".log")])
        sp.call(["rm", outName.replace(".out",".error")])
    if value==2:
        sp.call(["cat", outName])
        sp.call(["cat", outName.replace(".out",".log")])
        sp.call(["cat", outName.replace(".out",".error")])

def getProgressFromOut(outName,checkCompleted=False,printStatus=True,resubmit=False):
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
                        if printStatus: print colored("--"+line.split(".")[0]+" "+str(progress*10)+"%",color)
        if(processing == False and checkCompleted and lines>2):
            color = "red"
            if printStatus: 
                print colored(getNameFromOutFile(outName)+" failed",color)
                resubmitJob(outName,resubmit)
            if resubmit:
                resubmitJob(outName,resubmit)
        

def getStatusFromOut(outName,running,printStatus=True,resubmit=False):
    finished = False
    if os.path.exists(outName) == False or len(open(outName).readlines(  ))<=10:
        if running:
            color = "blue"
            if printStatus: print colored(getNameFromOutFile(outName)+" no output yet",color)
        else:
            color = "red"
            if printStatus:
                print colored(getNameFromOutFile(outName)+" failed",color)
                resubmitJob(outName,resubmit)
            if resubmit:
                resubmitJob(outName,resubmit)
    else:
        datasetName = getNameFromOutFile(outName)
        if len(open(outName).readlines(  )) !=0:
            with open(outName,"r") as f:
                next(f)
                for line in f:
                    if line.find("took")>=0:
                        finished = True
                        color = "green"
                        if printStatus: print colored(datasetName+" is finished, copy took "+getCopyTimeFromError(outName).strip()+"s, script"+line.split(")")[-1].replace("=","").replace("\n","")+" on "+getMachineFromOut(outName),color)
        if finished == False and running:
            getProgressFromOut(outName,True,printStatus,resubmit) 
        elif finished == False and running == False:
            color = "red"
            if printStatus:
                print colored(getNameFromOutFile(outName)+" failed"+" on "+getMachineFromOut(outName),color)
                resubmitJob(outName,resubmit)
            if resubmit:
                resubmitJob(outName,resubmit)
    
    return finished
    
def checkStatusFromLog(logPath,runningLogs,printSingleJobs=True,resubmit=False):
    allFinished = True
    nTotal = 0
    nRunning = 0
    nFinished = 0
    with open(logPath+"finished.txt","w") as f:
        for log in glob.glob(logPath+"*.submit"):
            nTotal+=1
            outFile = log.replace(".submit",".out")
            if outFile in runningLogs:
                nRunning+=1
                getStatusFromOut(outFile,True,printSingleJobs,resubmit)
                allFinished = False
            else:
                if getStatusFromOut(outFile,False,printSingleJobs,resubmit):
                    nFinished+=1
                    f.write(getNameFromOutFile(outFile)+"\n")
                else:
                    allFinished = False
    return allFinished,[nTotal,nRunning,nFinished]

def checkStatusFromQueue(printOutput=True,checkSuspended=False):
    jobs = getInfos()
    jobs = sorted(jobs, key=lambda l: l["JobStatus"]+l["ClusterId"])
    susJobs = []
    runningLogs = {}

    for job in jobs:
        if job["Args"]=="0":
            continue
        if job["Cmd"].split("/")[-1]=="run.sh":
            name = getNameFromFile(job["Args"])
        else:
            name = (job["Args"].split(",")[0]).split("/")[2]
        name += " "*max([0,(40-len(name))])
        jStatus = job["JobStatus"]
        runningLogs[job["Out"]] = (0 if jStatus=="1" else 1)
        if(printOutput == False): continue     # show job status only in nominal mode
        if jStatus == "1":
            print colored("\033[1m"+name+"    "+job["ClusterId"]+"       idle"+"\033[0m","yellow")
        elif jStatus == "2":
            print colored("\033[1m"+name+"    "+job["RemoteHost"]+"   "+job["ClusterId"]+"       running"+"\033[0m","green")
            getProgressFromOut(job["Out"])
        elif jStatus == "5":
            print colored("\033[1m"+name+"    "+job["ClusterId"]+"       held"+"\033[0m","red")
        elif jStatus == "7":
            susTime = (time.time()-int(job["LastSuspensionTime"]))/60.
            print colored("\033[1m"+name+"    "+job["RemoteHost"]+"       suspended since {:.2f} min".format(susTime)+"         "+job["RemoteHost"]+"\033[0m","red")
            getProgressFromOut(job["Out"])
            if susTime > 10 and checkSuspended :
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
    
def merge(distrLogPath,mergeAll=False,onlyHist=True):
    if mergeAll:
        value = 1
    else: value = input("All jobs finished succesfully. For merging the outputs enter 1:\n")
    if value==1:
        print "merging output for "+distrLogPath
        if onlyHist:
            print "merging only hists"
            sp.call(["python","mergeOutputs.py",distrLogPath,"--onlyHists"])
        else:
            sp.call(["python","mergeOutputs.py",distrLogPath])
    

def summaryDistributionJobs(year,runningLogs,resubmit,mergeAll,forceMergeAll,ignorePDF):
    for systPath in glob.glob("logs/"+year+"/*"):
        if ignorePDF:
            if systPath.find("PDF")>0 and systPath.find("PDF_ALPHAS")<0:
                continue
        distrLogPath = systPath+"/1.0/distributions/"
        if(os.path.exists(distrLogPath)):
            status = checkStatusFromLog(distrLogPath,runningLogs.keys(),False,resubmit)
            if status[0]:
                print colored(distrLogPath,"green",attrs=['bold'])
                if (mergeOutputs.mergeRequired(distrLogPath,True,False) or forceMergeAll) and resubmit==False:
                    merge(distrLogPath,(mergeAll or forceMergeAll))
            else:
                systName = getSystFromOutFile(distrLogPath)
                idleCheck = [value for key, value in runningLogs.iteritems() if key.find(systName)>0]
                nIdle = idleCheck.count(0)
                nFailed = status[1][0]-status[1][1]-status[1][2]
                print colored(distrLogPath,"cyan",attrs=['bold'])
                print colored("-----Idle:{}/{}".format(nIdle,status[1][0]),"yellow")
                print colored("-----Running:{}/{}".format(status[1][1]-nIdle,status[1][0]),"blue")
                print colored("-----Failed:{}/{}".format(nFailed,status[1][0]),"red")
                print colored("-----Finished:{}/{}".format(status[1][2],status[1][0]),"green")

def isDistributions(logPath):
    return logPath.find("distributions")>0
    
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",'--checkCompleted', type=str, default="", help="Checks status of all jobs in given log folder")
    parser.add_argument('--checkSuspended', action='store_true', help="Checks if jobs are suspended for more than 10 minutes, and offers resubmit")
    parser.add_argument('--distributions', action='store_true', help="Shows summary for distribution jobs per systematic")
    parser.add_argument('--resubmit', action='store_true', help="Automatic resubmit, avoids asking if job should be resubmitted for every single failed job")
    parser.add_argument('--mergeAll', action='store_true', help="Automatic merge, avoids asking if jobs should be merged")
    parser.add_argument('--forceMergeAll', action='store_true', help="Forces merge, even if not necessary due to updated input")
    parser.add_argument('--ignorePDF', action='store_true', help="Ignores all PDF jobs")
    args = parser.parse_args()
    
    if(args.checkCompleted == "" and args.distributions==False):  # print running jobs only on nominal mode
        printRunningJobs = True
    else:
        printRunningJobs = False
    
    runningLogs = checkStatusFromQueue(printRunningJobs,args.checkSuspended)
    
    if(args.distributions):
        summaryDistributionJobs("2018",runningLogs,args.resubmit,args.mergeAll,args.forceMergeAll,args.ignorePDF)
    
    if(args.checkCompleted != ""):
        if(os.path.exists(args.checkCompleted)):      # check status from log (does also inlcude finished jobs) and write finished names to list
            if checkStatusFromLog(args.checkCompleted,runningLogs.keys(),True,args.resubmit)[0]:
                if mergeOutputs.mergeRequired(args.checkCompleted,isDistributions(args.checkCompleted),False):
                    merge(args.checkCompleted)
                else:
                    print "No merge required"
        else:
            print "Wrong Path"

