#!/usr/bin/env python

import subprocess
import time
import os
from termcolor import colored

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

def getNameFromFile(fname):
    return fname.split("=")[-1]

def getProgressFromLog(logName):
    if len(open(logName).readlines(  )) !=0:
        with open(logName,"r") as f:
            next(f)
            for line in f:
                if line.find("Processing")==0:
                    progress=line.count(".")
                    color="yellow"
                    if progress==10:
                        color="green"
                    print colored("--"+line.split(".")[0]+" "+str(progress*10)+"%",color)

jobs = getInfos()
jobs = sorted(jobs, key=lambda l: l["JobStatus"]+l["ClusterId"])
susJobs = []

for job in jobs:
    if job["Cmd"].split("/")[-1]=="run.sh":
        name = getNameFromFile(job["Args"])
    else:
        name = (job["Args"].split(",")[0]).split("/")[2]
    name += " "*max([0,(40-len(name))])
    jStatus = job["JobStatus"]
    if jStatus == "1":
        print colored("\033[1m"+name+job["ClusterId"]+"       idle"+"\033[0m","yellow")
    elif jStatus == "2":
        print colored("\033[1m"+name+job["RemoteHost"]+"       running"+"\033[0m","green")
        getProgressFromLog(job["Out"])
    elif jStatus == "5":
        print colored("\033[1m"+name+job["ClusterId"]+"       held"+"\033[0m","red")
    elif jStatus == "7":
        susTime = (time.time()-int(job["LastSuspensionTime"]))/60.
        print colored("\033[1m"+name+job["RemoteHost"]+"       suspended since {:.2f} min".format(susTime)+"         "+job["RemoteHost"]+"\033[0m","red")
        if susTime > 10 :
            value = input("Job suspended for more than 10 Minutes, please enter 1 to kill job and start again or 0 to continue:\n")
            if value==1:
                print "resubmitting job"
                susJobs.append(job["ClusterId"])
    else:
        print "job status = ", jStatus

print getSummary()

for job in susJobs:
    os.system("condor_hold "+job)
if len(susJobs)!=0: os.system("condor_release dmeuser")

