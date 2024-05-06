#!/usr/bin/env python3
# need to be run without setting CMSSW environment!!
import argparse
import subprocess as sp
import os
import sys
import time
sys.path.append("../users")
from getPath import getPath

def getFullPath(pathName):
    return "file://"+pathName

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #  ~parser.add_argument("inputFiles")
    parser.add_argument("inputFiles", nargs='+', help='<Required> List of input files')
    parser.add_argument('target', type=str, help='<Required> Target on dCache starting at "mergedNtuple" (e.g. 2018/v06)')
    parser.add_argument('--dryRun', action='store_true', help='Run gfal-copy in dry-run mode (no copying, but only saying what would be done)')
    parser.add_argument('--srmcp', action='store_true', help='Use srmcp instead of gfal-copy')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing files (currently only for gfal-copy')
    args = parser.parse_args()
    
    sp.call(["gfal-mkdir","srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target])
    
    # no dry-run in srmcp
    if args.dryRun and args.srmcp:
        print("Error: srmcp has no dry-run option")
        exit()
    
    with open("tempFiles.txt","w") as f:
        i=0
        for pathName in args.inputFiles:
            i+=1
            f.write(getFullPath(pathName)+"\n")
            filename = pathName.split("/")[-1]
            print(getFullPath(pathName)+"  ({}/{})".format(i,len(args.inputFiles)))
            
            if args.srmcp:
                start_time = time.time()
                sp.call(["srmcp",getFullPath(pathName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target+"/"+filename])
                stop_time = time.time()
                print("Copy took {}s".format(round(stop_time-start_time,1)))
            else:
                if args.dryRun:
                    if args.overwrite:
                        sp.call(["gfal-copy",getFullPath(pathName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000","-f","--dry-run"])
                    else:
                        sp.call(["gfal-copy",getFullPath(pathName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000","--dry-run"])
                else:
                    if args.overwrite:
                        sp.call(["gfal-copy",getFullPath(pathName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000","-f"])
                    else:
                        sp.call(["gfal-copy",getFullPath(pathName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000"])
            



