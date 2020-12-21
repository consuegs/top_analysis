#!/usr/bin/env python2
import argparse
import os
import multiprocessing
import glob
import subprocess

class Range(object):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def __eq__(self, other):
        return self.start <= other <= self.end

#  ~toProcess_mc=["TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","ttG"]
toProcess_mc=["TTbar_diLepton"]
#  ~toProcess_data=["DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"]
#  ~toProcess_data=["DoubleMuon"]
toProcess_data=[]
#  ~toProcess_signal=["T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"]
#  ~toProcess_signal=["T1tttt_1200_800"]
toProcess_signal=[]

#############################################
# Select datasets to process
#############################################


parser = argparse.ArgumentParser()
parser.add_argument('-m', type=str, default="distributions", help="module")
parser.add_argument('-f', type=float, choices=[Range(0.0, 1.0)], default=1.0, help="process fraction")
parser.add_argument('--condor', action="store_true", help="running with condor")
args = parser.parse_args()

print "Running "+args.m

print "Process "+str(args.f*100)+"% of events"

print toProcess_mc
print toProcess_data
print toProcess_signal

requ_mem=1500   #standard value, allocated if not defined

if args.condor:
    sampleStr = ""
    for sel in [[toProcess_mc,"--mc_dataset="],[toProcess_data,"--data_dataset="],[toProcess_signal,"--signal_dataset="]]:
        for x in sel[0]:
                sampleStr=sel[1]+x
                if x=="TTbar_diLepton": requ_mem=6000
                with open("submitCondor.txt","w") as f:
                    f.write("""
Universe        = vanilla
Executable      = run.sh
Arguments       = {0} {1} {2}
Log             = logs/{3}.log
Output          = logs/{3}.out
Error           = logs/{3}.error
Request_Memory  = {4} Mb
Queue
""".format("-f"+str(args.f),args.m,sampleStr,x,str(requ_mem)))
                subprocess.call(["condor_submit", "submitCondor.txt"])

#  ~else: # local processing
        #  ~subprocess.call(["./run.sh","-f"+str(args.f),args.m])
        #  ~subprocess.call(["./../build/run.x","-f"+str(args.f),args.m])
    #  ~if len(toProcess)==1:
            #  ~run.run(dir+toProcess[0])
    #  ~else:
            #  ~files = [dir+x for x in toProcess]
            #  ~files.sort(key=os.path.getsize, reverse=True)
            #  ~p = multiprocessing.Pool()
            #  ~p.map(run.run, files)
