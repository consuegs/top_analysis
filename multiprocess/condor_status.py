#!/usr/bin/env python2
import argparse
import subprocess
import time
import os
import glob
from termcolor import colored
import re
import subprocess as sp
import mergeOutputs
import sys
sys.path.append("../users")
from getPath import getPath

def getInfos():
    #  ~out = subprocess.check_output(["condor_q", "-long", "-grid"])
    out = subprocess.check_output(["condor_q", "-long"])
    #  ~for jobStrings in out.split("\n\n"):
        #  ~for line in jobStrings.split("\n"):
            #  ~if line:
                #  ~break
                #  ~print [line.split(" = ")]
    return [ dict([line.replace("\"","").split(" = ") for line in jobStrings.split("\n") if (" = " in line and "TransferOutputRemaps" not in line)]) for jobStrings in out.split("\n\n") if jobStrings ]

def getSummary():
    out = subprocess.check_output(["condor_q", "-constraint","JobUniverse==9"])
    #  ~return out.split("\n")[-3]
    return out.split("\n")[-5]

def getNameFromOutFile(outFile):
    moduleName = outFile.split("/")[-2]
    return outFile.split("/")[-1].replace(moduleName+"_","").replace(".out","")

def getSystFromOutFile(outFile):
    systName = outFile.split("/")[2]
    return systName

def getNameFromFile(fname):
    return fname.split(" ")[1]+" "+(fname.split("dataset=")[-1]).split("/")[0]
    
def getMachineFromOut(outName):
    #  ~logName = outName.replace(".out",".log")
    machine = "machine not known"
    #  ~with open(logName,"r") as f:
        #  ~for line in f:
            #  ~if line.find("alias") >= 0:
                #  ~machine = line.split("alias=")[-1]
                #  ~machine = machine.split("&")[0]
    with open(outName,"r") as f:
        for line in f:
            if line.find("grid-") >= 0:
                machine = line[:-1]
                break
    return machine
    
def getCopyTimeFromError(outName):
    logName = outName.replace(".out",".error")
    with open(logName,"r") as f:
        for line in f:
            if line.find("seconds") >= 0:
                line = line.split("in")[-1]
                time = line.split("seconds")[0]
                speed = line.split("seconds")[1]
    return time,speed
    
def checkErrorFile(outName):
    with open(outName,"r") as f:
        for line in f:
            if line.find("ERROR (70)")>=0:
                return False
    logName = outName.replace(".out",".error")
    if logName.find("TUnfold")>=0:
        return True
    size = os.path.getsize(logName)/(1024*1024)
    with open(logName,"r") as f:
        for line in f:
            if line.find("NaN")>=0:
                return False
            if line.find("SysError")>=0:
                return False
    return size < 1
    
    
def resubmitJob(outName,automaticResubmit=False):
    if automaticResubmit: value = 1
    else :
        value = input("For resubmitting "+outName+" enter 1. For showing .out, .log and .error enter 2:\n")
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

def getProgressFromOut(outName,checkCompleted=False,printStatus=True,resubmit=False,running=False):
    processing = running
    if (outName.find("TUnfold")<0):
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
        if checkErrorFile(outName) == False:
            finished = False
        elif len(open(outName).readlines(  )) !=0:
            with open(outName,"r") as f:
                next(f)
                for line in f:
                    if line.find("took")>=0:
                        finished = True
                        color = "green"
                        if printStatus: print colored(datasetName+" is finished, copy took "+getCopyTimeFromError(outName)[0].strip()+"s with "+getCopyTimeFromError(outName)[1].strip()+", script"+line.split(")")[-1].replace("=","").replace("\n","")+" on "+getMachineFromOut(outName),color)
        if finished == False and running:
            getProgressFromOut(outName,True,printStatus,resubmit,True) 
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
        
        if job.has_key("RemoteHost") == False:
            job["RemoteHost"] = "Grid"
        
        if job["Args"]=="0":
            continue
        #  ~if job["Cmd"].split("/")[-1]=="run.sh":
        if job["Cmd"].split("/")[-1].find("run")>=0:
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
            print colored("\033[1m"+name+"    "+job["RemoteHost"]+"       suspended since {:.2f} min".format(susTime)+"         "+job["RemoteHost"]+"       "+job["ClusterId"]+"\033[0m","red")
            getProgressFromOut(job["Out"])
            if susTime > 10 and checkSuspended :
                value = input("Job suspended for more than 10 Minutes, please enter 1 to kill job and start again or 0 to continue:\n")
                if value==1:
                    print "resubmitting job"
                    susJobs.append(job["ClusterId"])
        #  ~else:
            #  ~print "job status = ", jStatus
    
    if(printOutput):
        print getSummary()

        for job in susJobs:
            os.system("condor_hold "+job)
        if len(susJobs)!=0: os.system("condor_release {}".format(getPath("name")))
    
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

def dCacheTreeUpload(year):
    treePath = mergeOutputs.getTreePath("logs/{0}/Nominal/1.0/distributions/".format(str(year))).replace("Nominal/","")
    targetPath = "minTrees/"+str(year)+(treePath.split(str(year))[-1])
    command = "eval `scram unsetenv -sh`;","gfal-copy",treePath,"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/".format(getPath("dCacheBasePath"))+targetPath,"-t 36000","-f","-r","--dry-run"
    command = " ".join(command)
    sp.call(command,shell=True)

def summaryJobs(year,runningLogs,resubmit,mergeAll,forceMergeAll,ignorePDF,module="distributions"):
    allFinished = True;
    for systPath in glob.glob("logs/"+year+"/*"):
        if ignorePDF:
            if systPath.find("PDF")>0 and systPath.find("PDF_ALPHAS")<0:
                continue
        distrLogPath = systPath+"/1.0/"+module+"/"
        if(os.path.exists(distrLogPath)):
            status = checkStatusFromLog(distrLogPath,runningLogs.keys(),False,resubmit)
            if status[0]:
                print colored(distrLogPath,"green",attrs=['bold'])
                if (module != "TUnfold_binning"):
                    if (mergeOutputs.mergeRequired(distrLogPath,module=="distributions",False) or forceMergeAll) and resubmit==False:
                        merge(distrLogPath,(mergeAll or forceMergeAll))
            else:
                systName = getSystFromOutFile(distrLogPath)
                idleCheck = [value for key, value in runningLogs.iteritems() if (key.find(systName)>0 and key.find(year)>0)]
                nIdle = idleCheck.count(0)
                nFailed = status[1][0]-status[1][1]-status[1][2]
                print colored(distrLogPath,"cyan",attrs=['bold'])
                print colored("-----Idle:{}/{}".format(nIdle,status[1][0]),"yellow")
                print colored("-----Running:{}/{}".format(status[1][1]-nIdle,status[1][0]),"blue")
                print colored("-----Failed:{}/{}".format(nFailed,status[1][0]),"red")
                print colored("-----Finished:{}/{}".format(status[1][2],status[1][0]),"green")
                allFinished = False;
    #  ~if allFinished and module == "distributions":
        #  ~value = input("All jobs finished succesfully. For uploading the minTrees to dCache enter 1:\n")
        #  ~if value == 1:
            #  ~dCacheTreeUpload(year)

def isDistributions(logPath):
    return logPath.find("distributions")>0
    
    

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-c",'--checkCompleted', type=str, default="", help="Checks status of all jobs in given log folder")
    parser.add_argument('-y', type=str, help="year to be set as ANALYSIS_YEAR_CONFIG",default="2018")
    parser.add_argument('--checkSuspended', action='store_true', help="Checks if jobs are suspended for more than 10 minutes, and offers resubmit")
    parser.add_argument('--distributions', action='store_true', help="Shows summary for distribution jobs per systematic")
    parser.add_argument('--bTagEff', action='store_true', help="Shows summary for bTaggEff jobs per systematic")
    parser.add_argument('--triggerEff', action='store_true', help="Shows summary for triggerEff jobs per systematic")
    parser.add_argument('--TUnfold', action='store_true', help="Shows summary for TUnfold jobs per systematic")
    parser.add_argument('--resubmit', action='store_true', help="Automatic resubmit, avoids asking if job should be resubmitted for every single failed job")
    parser.add_argument('--mergeAll', action='store_true', help="Automatic merge, avoids asking if jobs should be merged")
    parser.add_argument('--forceMergeAll', action='store_true', help="Forces merge, even if not necessary due to updated input")
    parser.add_argument('--ignorePDF', action='store_true', help="Ignores all PDF jobs")
    parser.add_argument('--repeat', action='store_true', help="Repeats tasks after 15 minutes of sleeping")
    args = parser.parse_args()
    
    while True:
        if(args.checkCompleted == "" and args.distributions==False and args.bTagEff==False and args.TUnfold==False):  # print running jobs only on nominal mode
            printRunningJobs = True
        else:
            printRunningJobs = False
        
        runningLogs = checkStatusFromQueue(printRunningJobs,args.checkSuspended)
        
        if(args.distributions):
            summaryJobs(args.y,runningLogs,args.resubmit,args.mergeAll,args.forceMergeAll,args.ignorePDF)
        elif(args.bTagEff):
            summaryJobs(args.y,runningLogs,args.resubmit,args.mergeAll,args.forceMergeAll,args.ignorePDF,"bTagEff")
        elif(args.TUnfold):
            summaryJobs(args.y,runningLogs,args.resubmit,args.mergeAll,args.forceMergeAll,args.ignorePDF,"TUnfold_binning")
        elif(args.triggerEff):
            summaryJobs(args.y,runningLogs,args.resubmit,args.mergeAll,args.forceMergeAll,args.ignorePDF,"triggerEff")
        
        if(args.checkCompleted != ""):
            if(os.path.exists(args.checkCompleted)):      # check status from log (does also inlcude finished jobs) and write finished names to list
                if checkStatusFromLog(args.checkCompleted,runningLogs.keys(),True,args.resubmit)[0]:
                    if mergeOutputs.mergeRequired(args.checkCompleted,isDistributions(args.checkCompleted),False) or args.forceMergeAll:
                        merge(args.checkCompleted)
                    else:
                        print "No merge required"
            else:
                print "Wrong Path"
        
        if (args.repeat):
            print time.strftime("%H-%M-%S:"), "Sleeping for 15 minutes."
            time.sleep(900)
        else:
            break

