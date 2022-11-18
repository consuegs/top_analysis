#!/usr/bin/env python2
# need to be run without setting CMSSW environment!!
import argparse
import subprocess as sp
import os
import sys
sys.path.append("../users")
from getPath import getPath

def getFullPath(fileName):
    return "file://"+fileName

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    #  ~parser.add_argument("inputFiles")
    parser.add_argument("inputFiles", nargs='+', help='<Required> List of input files')
    parser.add_argument('target', type=str, help='<Required> Target on dCache starting at "mergedNtuple" (e.g. 2018/v06)' )
    parser.add_argument('--dryRun', action='store_true' )
    args = parser.parse_args()
    
    sp.call(["gfal-mkdir","srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target])
    
    with open("tempFiles.txt","w") as f:
        for fileName in args.inputFiles:
            f.write(getFullPath(fileName)+"\n")
            print(getFullPath(fileName))
            
            if args.dryRun:
                sp.call(["gfal-copy",getFullPath(fileName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000","--dry-run"])
            else:
                sp.call(["gfal-copy",getFullPath(fileName),"srm://grid-srm.physik.rwth-aachen.de:8443/srm/managerv2?SFN={}/mergedNtuple/".format(getPath("dCacheBasePath"))+args.target,"-t 36000","-f"])
            



