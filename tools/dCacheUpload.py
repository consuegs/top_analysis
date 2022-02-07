#!/usr/bin/env python2
# need to be run without setting CMSSW environment!!
import argparse
import subprocess as sp
import os

def getFullPath(fileName):
    return "file://"+fileName

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #  ~parser.add_argument("inputFiles")
    parser.add_argument("inputFiles", nargs='+', help='<Required> List of input files')
    parser.add_argument('target', type=str, help='<Required> Target on dCache starting at "mergedNtuple" (e.g. 2018/v06)' )
    parser.add_argument('--dryRun', action='store_true' )
    args = parser.parse_args()
    
    with open("tempFiles.txt","w") as f:
        for fileName in args.inputFiles:
            f.write(getFullPath(fileName)+"\n")
    
    if args.dryRun:
        sp.call(["gfal-copy","--from-file",os.getcwd()+"/tempFiles.txt","srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/dmeuser/mergedNtuple/"+args.target,"-t 36000","-f","--dry-run"])
    else:
        sp.call(["gfal-copy","--from-file",os.getcwd()+"/tempFiles.txt","srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN=/pnfs/physik.rwth-aachen.de/cms/store/user/dmeuser/mergedNtuple/"+args.target,"-t 36000","-f"])
    



