#Tool to create datacards

import os, re, subprocess as sp, math
from distutils import spawn
from array import array
from operator import add,sub
import ROOT as rt
import numpy as np
from array import array

backgrounds=["GJ","TTcomb","Vg","diboson","efake"]

class Scan:
    T5gg,T5Wg,T6Wg,T6gg,GGM,TChiWg,TChiNg,GGM_M1_M2,GGM_M1_M3,TChiNg_gg,TChiNg_gz,TChiNg_zz,TChiNg_gg_C1N2,T5Wg_Wg,T5Wg_thirds=range(15)

#~ def prepareDatacard(obs,count,estat,esyst,eISR,correlated,mcUncertainties,pointName,sScan):
def prepareDatacard(obs,count,estat,esyst,eISR,ePileUp,eJES,correlated,mcUncertainties,pointName,sScan):
    #~ directory=outdir+"datacards/"+selection+"/"+sScan
    #~ directory=outdir+"datacards_newNui/"+selection+"/"+sScan
    directory=outdir+"datacards_final/"+selection+"/"+sScan
    if rho!=-0.0: directory=directory+"-"+str(rho)
    if not os.path.exists(directory):
		os.makedirs(directory)
    if sScan=="TChiNg_gg" or sScan=="TChiNg_gz" or sScan=="TChiNg_zz" or sScan=="TChiNg_gg_C1N2":
        dataCardPath=directory+"/datacard_%s.txt"%(sScan+"_"+(pointName.split("_"))[1])
    elif sScan=="TChiNg" or sScan=="TChiWg" :
        dataCardPath=directory+"/datacard_%s.txt"%pointName
    else :
        dataCardPath=directory+"/datacard_%s.txt"%pointName
    nBins=len(obs)
    card="# "+pointName
    card+="""
imax %d  number of channels
jmax *  number of backgrounds ('*' = automatic)
kmax *  number of nuisance parameters (sources of systematical uncertainties)
------------
bin                        bin1       bin2       bin3       bin4
observation          """%nBins
    for o in obs:
        card+="%10d "%o

    backgrounds=["GJ","rare","Vg","efake"]
    card+="\n------------\n"
    allProc=["sig"]+backgrounds
    card+="bin                  "
    for i in range(nBins):
        card+=("      bin%d "%(i+1))*(len(backgrounds)+1)
    card+="\nprocess              "
    for i in range(nBins):
        for bkg in allProc:
            card+="%10s "%bkg
    card+="\nprocess              "
    for i in range(nBins):
        for j in range(len(backgrounds)+1):
            card+="%10d "%j
    card+="\nrate                 "
    for i in range(nBins):
        for b in allProc:
            card+="%10.2f "%count[b][i]

    card+="\n------------\n"
    # statistical errors
    col=0
    maxCol=nBins*len(allProc)
    fill="%10s "%"-"
    for i in range(nBins):
        for b in allProc:
            name="statB%d%s"%(i+1,b)
            card+="%-16s lnN "%name
            card+=fill*col
            card+="%10.2f "%estat[b][i]
            card+=fill*(maxCol-col-1)
            card+="\n"
            col+=1
    # systematic uncertainties
    for b,errors in esyst.iteritems():
        index=allProc.index(b)
        name="syst_%s"%b
        if b=="rare": name="xs"
        elif b=="efake": name="e_to_pho_syst"
        if b=="sig":
            name="GenMet"
            card+="%-16s lnN "%name
        else :
            card+="%-16s lnN "%name
        # looping bins (all the same anyway...)
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    for b,errors in eISR.iteritems():
        index=allProc.index(b)
        name="ISR"
        card+="%-16s lnN "%name
        # looping bins
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    for b,errors in ePileUp.iteritems():
        index=allProc.index(b)
        name="PU"
        card+="%-16s lnN "%name
        # looping bins
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    for b,errors in eJES.iteritems():
        index=allProc.index(b)
        name="JES"
        card+="%-16s lnN "%name
        # looping bins
        for e in errors:
            card+=fill*index
            card+="%10.2f "%e
            card+=fill*(len(allProc)-index-1)
        card+="\n"
    # correlated part of systematic uncertainties
    corr={}
    name="cor_"
    for b in correlated:
        corr[allProc.index(b)]=correlated[b]
        name+=b
    card+="%-16s lnN "%name
    for i in range(nBins):
        lastIndex=-1
        for ind in corr:
            card+=fill*(ind-lastIndex-1)
            card+="%10.2f "%corr[ind][i]
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    # fully correlated uncertainties for all MC samples
    indices=[]
    mc_bkg=["rare","sig"]
    for s in mc_bkg:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    for name,val in mcUncertainties.iteritems():
        if name=="phSF": name="PhotonSF"
        card+="%-16s lnN "%name
        for i in range(nBins):
            lastIndex=-1
            for ind in indices:
                card+=fill*(ind-lastIndex-1)
                card+="%10.2f "%val
                lastIndex=ind
            card+=fill*(len(allProc)-lastIndex-1)
        card+="\n"

    indices=[]
    for s in ["GJ"]:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    card+="%-16s lnN "%"GJ_systematics"
    for i in range(nBins):
        lastIndex=-1
        for ind in indices:
            card+=fill*(ind-lastIndex-1)
            if i == 0:
               card+="%10.2f "%1.097
            elif i == 1:
               card+="%10.2f "%1.0797
            elif i == 2:
               card+="%10.2f "%1.0912
            elif i == 3:
               card+="%10.2f "%1.335
            else: print "error too many bins"
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    indices=[]
    for s in ["Vg"]:
        indices.append(allProc.index(s))
    indices=sorted(indices)
    card+="%-16s lnN "%"Vg_systematics"
    for i in range(nBins):
        lastIndex=-1
        for ind in indices:
            card+=fill*(ind-lastIndex-1)
            if i == 0:
               card+="%10.2f "%1.0722
            elif i == 1:
               card+="%10.2f "%1.08045
            elif i == 2:
               card+="%10.2f "%1.08787
            elif i == 3:
               card+="%10.2f "%1.110
            else: print "error too many bins"
            lastIndex=ind
        card+=fill*(len(allProc)-lastIndex-1)
    card+="\n"
    with open(dataCardPath, "w") as f:
        f.write(card)
    return dataCardPath

def getSignalYield_Unc(point,sScan):
    f=rt.TFile(outdir+signal_scan,"read")
    hist=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point)
    hist_gen=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"_gen")
    hist_JESu=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"_JESu")
    hist_JESd=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"_JESd")
    hist_ISR=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"SRErrISR")
    hist_PU=f.Get(sScan+"/PuUnc/"+point)
    if math.isnan(hist.Integral()):
        # input tree for model was bad->Ngen=0->weight=inf
        f.Close()
        return None
    sigYield=[]
    statErr=[]
    metErr=[]
    jseErr=[]
    ISRErr=[]
    puErr=[]
    for i in range(1,5):
        y=hist.GetBinContent(i)
        yg=hist_gen.GetBinContent(i)
        yJESu=hist_JESu.GetBinContent(i)
        yJESd=hist_JESd.GetBinContent(i)
        yISR=hist_ISR.GetBinContent(i)        
        sigYield.append(.5*(y+yg)) # use mean of reco and gen met yields
        if (y > 0):
           statErr.append(1+hist.GetBinError(i)/y)
           metErr.append(1+abs(y-yg)/(y+yg)) # uncertainty: half of the difference
           jseErr.append(1+abs(yJESu-yJESd)/(yJESu+yJESd)) # uncertainty: half of the difference
           ISRErr.append(1+(abs(y-yISR)/y)) # uncertainty: full difference
           puErr.append(1+hist_PU.GetBinContent(1))
        else:
           statErr.append(1+0.1)
           metErr.append(1+0.1) # uncertainty: half of the difference
           jseErr.append(1+0.1) # uncertainty: half of the difference
           ISRErr.append(1+0.1) # uncertainty: full difference
           puErr.append(1+0.1)
    f.Close()
    return sigYield,statErr,metErr,ISRErr,puErr,jseErr
    
def getSignalYield(point,sScan):
    f=rt.TFile(outdir+signal_scan,"read")
    hist=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point)
    hist_gen=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"_gen")
    hist_ISR=f.Get(sScan+"/pre_ph165/c_MET300/MT300/STg/"+point+"SRErrISR")
    if math.isnan(hist.Integral()):
        # input tree for model was bad->Ngen=0->weight=inf
        f.Close()
        return None
    sigYield=[]
    statErr=[]
    metErr=[]
    ISRErr=[]
    for i in range(1,5):
        y=hist.GetBinContent(i)
        yg=hist_gen.GetBinContent(i)
        yISR=hist_ISR.GetBinContent(i)        
        sigYield.append(.5*(y+yg)) # use mean of reco and gen met yields
        if (y > 0):
           statErr.append(1+hist.GetBinError(i)/y)
           metErr.append(1+abs(y-yg)/(y+yg)) # uncertainty: half of the difference
           ISRErr.append(1+(abs(y-yISR)/y)) # uncertainty: full difference
        else:
           statErr.append(1+0.1)
           metErr.append(1+0.1) # uncertainty: half of the difference
           ISRErr.append(1+0.1) # uncertainty: half of the difference
    f.Close()
    return sigYield,statErr,metErr,ISRErr

def getSignalContamination(point,sScan):
    f=rt.TFile(outdir+signal_scan,"read")
    hist=f.Get(sScan+"/pre_ph165/c_MET100/MT100/METl300vMTl300/absphiMETnJetPh/"+point)
    c=hist.Integral()
    if math.isnan(c):
        # input tree for model was bad->Ngen=0->weight=inf
        f.Close()
        return None
    f.Close()
    return c

def decomposeCorrelations(esyst,count):
    cor={"Vg":[],"GJ":[]}
    unc={"GJ":[]}
    for i in range(len(count["Vg"])):
        eX=(esyst["Vg"][i]-1)*count["Vg"][i]
        eY=(esyst["GJ"][i]-1)*count["GJ"][i]
        eXtoY=rho*eY
        eYunc=(1-rho**2)**.5 * eY
        # print eY,eXtoY,eYunc
        cor["Vg"].append((count["Vg"][i]+eX)/count["Vg"][i])
        cor["GJ"].append((count["GJ"][i]+eXtoY)/count["GJ"][i])
        unc["GJ"].append((count["GJ"][i]+eYunc)/count["GJ"][i])
    return cor,unc

def getMasses(point,scan):
    if scan==Scan.T5gg:   pattern="T5gg_(.*)_(.*)"
    elif scan==Scan.T5Wg:   pattern="T5Wg_(.*)_(.*)"
    elif scan==Scan.T5Wg_Wg:   pattern="T5Wg_(.*)_(.*)"
    elif scan==Scan.T5Wg_thirds:   pattern="T5Wg_(.*)_(.*)"
    elif scan==Scan.T6gg:   pattern="T6gg_(.*)_(.*)"   
    elif scan==Scan.T6Wg:   pattern="T6Wg_(.*)_(.*)"
    elif scan==Scan.GGM:    pattern=".*_M2_(.*)_M1_(.*)"
    elif scan==Scan.TChiWg: pattern="TChiWG_(.*)"
    elif scan==Scan.TChiNg: pattern="TChiNG_(.*)"
    elif scan==Scan.GGM_M1_M2: pattern="GGM_M1_M2_(.*)_(.*)"
    elif scan==Scan.GGM_M1_M3: pattern="GGM_M1_M3_(.*)_(.*)"
    elif scan==Scan.TChiNg_gg: pattern="TChiNG_(.*)"
    elif scan==Scan.TChiNg_gz: pattern="TChiNG_(.*)"
    elif scan==Scan.TChiNg_zz: pattern="TChiNG_(.*)"
    elif scan==Scan.TChiNg_gg_C1N2: pattern="TChiNG_(.*)"
      
    m=re.search(pattern,point)
    masses=[]
    if m and scan==Scan.TChiWg and (len(m.groups())==1):
        masses=[int(m.groups()[0]),0]
    elif m and (scan==Scan.TChiNg or scan==Scan.TChiNg_gg or scan==Scan.TChiNg_gz or scan==Scan.TChiNg_zz or scan==Scan.TChiNg_gg_C1N2) and (len(m.groups())==1):
        masses=[int(m.groups()[0]),0]        
    elif m and (len(m.groups())==2):
        for s in m.groups():
            masses.append(int(s))
    else:
        print "don't know what this is:",point
        exit(-1)
    return masses


def fillDatacards(scan,SignalTree):
    print "using",os.environ["CMSSW_VERSION"]
    print

    scanFile=basedir
    if scan==Scan.T5gg:   scanFile+="T5gg_scan.txt"
    elif scan==Scan.T5Wg:   scanFile+="T5Wg_scan.txt"
    elif scan==Scan.T5Wg_Wg:   scanFile+="T5Wg_scan.txt"
    elif scan==Scan.T5Wg_thirds:   scanFile+="T5Wg_scan.txt"
    elif scan==Scan.T6gg:   scanFile+="T6gg_scan.txt"
    elif scan==Scan.T6Wg:   scanFile+="T6Wg_scan.txt"
    elif scan==Scan.GGM:    scanFile+="GGM_WinoBino_scan.txt"
    elif scan==Scan.TChiWg: scanFile+="TChiWg_scan.txt"
    elif scan==Scan.TChiNg: scanFile+="TChiNg_scan.txt"
    elif scan==Scan.GGM_M1_M2: scanFile+="GGM_M1_M2_scan.txt"
    elif scan==Scan.GGM_M1_M3: scanFile+="GGM_M1_M3_scan.txt"
    elif scan==Scan.TChiNg_gg: scanFile+="TChiNg_scan.txt"
    elif scan==Scan.TChiNg_gz: scanFile+="TChiNg_scan.txt"
    elif scan==Scan.TChiNg_zz: scanFile+="TChiNg_scan.txt"
    elif scan==Scan.TChiNg_gg_C1N2: scanFile+="TChiNg_scan.txt"
    
    sScan="unkown_scan"
    if scan==Scan.T5gg:   sScan="T5gg"
    elif scan==Scan.T5Wg:   sScan="T5Wg"
    elif scan==Scan.T5Wg_Wg:   sScan="T5Wg_Wg"
    elif scan==Scan.T5Wg_thirds:   sScan="T5Wg_thirds"
    elif scan==Scan.T6gg:   sScan="T6gg"
    elif scan==Scan.T6Wg:   sScan="T6Wg"
    elif scan==Scan.GGM:    sScan="GGM"
    elif scan==Scan.TChiWg: sScan="TChiWg"
    elif scan==Scan.TChiNg: sScan="TChiNg"
    elif scan==Scan.GGM_M1_M2: sScan="GGM_M1_M2"
    elif scan==Scan.GGM_M1_M3: sScan="GGM_M1_M3"
    elif scan==Scan.TChiNg_gg: sScan="TChiNg_gg"
    elif scan==Scan.TChiNg_gz: sScan="TChiNg_gz"
    elif scan==Scan.TChiNg_zz: sScan="TChiNg_zz"
    elif scan==Scan.TChiNg_gg_C1N2: sScan="TChiNg_gg_C1N2"
    
    count={}
    estat={}
    esyst={}
    eISR={}
    ePileUp={}
    eJES={}
    f=rt.TFile(outdir+"yields.root","read")
    for bkg in backgrounds:
        #~ hist=f.Get("pre_ph165/c_MET300/MT300/STg/"+bkg)
        hist=f.Get("pre_ph165/c_MET300/MT300/"+selection+"/STg/"+bkg)
        count[bkg]=[hist.GetBinContent(i) for i in range(1,5)]
        # deal with zero bin content
        for i in range(1,5):
            if count[bkg][i-1]==0:
                count[bkg][i-1]=-1
        # subtract "additional" 1 from count index, because 1st bin is already left out
        estat[bkg]=[1+hist.GetBinError(i)/count[bkg][i-1] for i in range(1,5)]
        #~ hist=f.Get("pre_ph165/c_MET300/MT300/STg/"+bkg+"_esyst")
        hist=f.Get("pre_ph165/c_MET300/MT300/"+selection+"/STg/"+bkg+"_esyst")
        esyst[bkg]=[1+hist.GetBinContent(i)/count[bkg][i-1] for i in range(1,5)]
        # remove negative count used for proper error calculation
        for i in range(1,5):
            if count[bkg][i-1]==-1:
                count[bkg][i-1]=0
    #~ hist=f.Get("pre_ph165/c_MET300/MT300/STg/data")
    hist=f.Get("pre_ph165/c_MET300/MT300/"+selection+"/STg/data")
    obs=[int(hist.GetBinContent(i)) for i in range(1,5)]
    f.Close()

    # the same and 100% correlated for all MC
    mcUncertainties={
        "lumi"   : 1+ 2.6 /100,
        "phSF"   : 1+ 2. /100,
        "trigger": 1+ 0.43 /100,
    }
    
    #Merge diboson and ttgamma to rare background
    for i in range(0,4):
        count["TTcomb"][i]=count["diboson"][i]+count["TTcomb"][i]
        estat["TTcomb"][i]=np.sqrt((estat["TTcomb"][i]-1)**2+(estat["diboson"][i]-1)**2)+1
        esyst["TTcomb"][i]=1.50     #Change from 30% to 50%
    
    count["rare"] = count.pop("TTcomb")    
    estat["rare"] = estat.pop("TTcomb")    
    esyst["rare"] = esyst.pop("TTcomb")    
    del count["diboson"]
    del estat["diboson"]
    del esyst["diboson"]

    points=[]
    list_m1=[]
    list_m2=[]
    with open(scanFile) as f:
        for p in f.read().split():
            p=p.split(".")[0]
            if scan==Scan.T5gg or scan==Scan.T5Wg or scan==Scan.T5Wg_Wg or scan==Scan.T5Wg_thirds or scan==Scan.T6Wg or scan==Scan.T6gg:
                m2,m1=getMasses(p,scan)
                list_m1.append(m1)
                list_m2.append(m2)
                #if m2<1100: continue
            elif scan==Scan.GGM_M1_M2 or scan==Scan.GGM_M1_M3:
                m1,m2=getMasses(p,scan)
                list_m1.append(m1)
                list_m2.append(m2)
            points.append(p)
    
    if SignalTree:
        if scan==Scan.T5gg or scan==Scan.T5Wg or scan==Scan.T5Wg_Wg or scan==Scan.T5Wg_thirds or scan==Scan.T6Wg or scan==Scan.T6gg:
            y=sorted(list(set(list_m1)))
            x=sorted(list(set(list_m2)))
        elif scan==Scan.GGM_M1_M2 or scan==Scan.GGM_M1_M3:
            x=sorted(list(set(list_m1)))
            y=sorted(list(set(list_m2)))
            
        #2D Signal-Histograms for combination
        hist_SigRate_0=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigRate_1=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigRate_2=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigRate_3=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        
        hist_SigStat_0=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigStat_1=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigStat_2=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigStat_3=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        
        hist_SigSys_0=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigSys_1=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigSys_2=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigSys_3=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        
        hist_SigISR_0=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigISR_1=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigISR_2=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        hist_SigISR_3=rt.TH2F("","",len(x)-1,array("f",x),len(y)-1,array("f",y))
        
        hist_GJStat=rt.TH1F("","",4,-0.5,3.5)
        hist_TTcombStat=rt.TH1F("","",4,-0.5,3.5)
        hist_VgStat=rt.TH1F("","",4,-0.5,3.5)
        hist_DibosonStat=rt.TH1F("","",4,-0.5,3.5)
        hist_EfakeStat=rt.TH1F("","",4,-0.5,3.5)
        
        hist_TTcombSys=rt.TH1F("","",4,-0.5,3.5)
        hist_DibosonSys=rt.TH1F("","",4,-0.5,3.5)
        hist_GJSys=rt.TH1F("","",4,-0.5,3.5)
        hist_EfakeSys=rt.TH1F("","",4,-0.5,3.5)
        
        hist_corGJ=rt.TH1F("","",4,-0.5,3.5)
        hist_corVg=rt.TH1F("","",4,-0.5,3.5)
        hist_trigger=rt.TH1F("","",4,-0.5,3.5)
        hist_phSF=rt.TH1F("","",4,-0.5,3.5)
        hist_lumi=rt.TH1F("","",4,-0.5,3.5)
        hist_GJFit=rt.TH1F("","",4,-0.5,3.5)
        hist_VgFit=rt.TH1F("","",4,-0.5,3.5)
        
    graph_puUnc=rt.TGraph2D()
    graph_jesUnc=rt.TGraph2D()

    for i,point in enumerate(points):
        print point
        if scan==Scan.GGM_M1_M2 or scan==Scan.GGM_M1_M3:
            m1,m2=getMasses(point,scan)
            x_val,y_val=m1,m2
        else:
            m2,m1=getMasses(point,scan)
            x_val,y_val=m2,m1
        key=m2
        if scan==Scan.GGM: key=m2*100000+m1
        
        #~ sigYield=getSignalYield(point,sScan)
        sigYield=getSignalYield_Unc(point,sScan)
        
        contamin=getSignalContamination(point,sScan)
        if not sigYield:
            print " broken!"
            continue # broken point
        c,st,sy,sISR,sPileUp,sJES=dict(count),dict(estat),dict(esyst),dict(eISR),dict(ePileUp),dict(eJES)
        # decompose partially correlated parts:
        cor,unc=decomposeCorrelations(esyst,count)
        del sy["GJ"]
        del sy["Vg"]
        sy.update(unc) # re-add the uncorrelated part
        # add signal
        c['sig']=sigYield[0]
        st['sig']=sigYield[1]
        sy['sig']=sigYield[2]
        sISR['sig']=sigYield[3]
        sPileUp['sig']=sigYield[4]
        sJES['sig']=sigYield[5]
        # subtract bkg overestimation from signal contamination
        subtractGJ=[x*contamin for x in c["GJ"]]
        subtractVG=[x*contamin for x in c["Vg"]]
        subtract=map(add, subtractGJ, subtractVG)
        # actual subtraction done from S to not destroy the B-only hypothesis
        c["sig"]=map(sub,c["sig"],subtract)
        # avoid negative counts
        c["sig"]=[max(0,x) for x in c["sig"]]
        datacard=prepareDatacard(obs,c,st,sy,sISR,sPileUp,sJES,cor,mcUncertainties,point,sScan)
        
        #~ graph_puUnc.SetPoint(graph_puUnc.GetN()+1,x_val,y_val,np.mean(sigYield[4])-1)
        #~ graph_jesUnc.SetPoint(graph_puUnc.GetN()+1,x_val,y_val,np.mean(sigYield[5])-1)
        
        
        """
        if SignalTree:
            i=0
            for hist in [hist_SigRate_0,hist_SigRate_1,hist_SigRate_2,hist_SigRate_3]:
                hist.Fill(x_val,y_val,c["sig"][i])
                i+=1
            i=0
            for hist in [hist_SigStat_0,hist_SigStat_1,hist_SigStat_2,hist_SigStat_3]:
                hist.Fill(x_val,y_val,st["sig"][i])
                i+=1
            i=0
            for hist in [hist_SigSys_0,hist_SigSys_1,hist_SigSys_2,hist_SigSys_3]:
                hist.Fill(x_val,y_val,sy["sig"][i])
                i+=1
            i=0
            for hist in [hist_SigISR_0,hist_SigISR_1,hist_SigISR_2,hist_SigISR_3]:
                hist.Fill(x_val,y_val,sISR["sig"][i])
                i+=1
    
    if SignalTree:
        for i in xrange(4):
            hist_GJStat.Fill(i,st["GJ"][i])
            hist_TTcombStat.Fill(i,st["TTcomb"][i])
            hist_VgStat.Fill(i,st["Vg"][i])
            hist_DibosonStat.Fill(i,st["diboson"][i])
            hist_EfakeStat.Fill(i,st["efake"][i])
        
            hist_TTcombSys.Fill(i,sy["TTcomb"][i])
            hist_DibosonSys.Fill(i,sy["diboson"][i])
            hist_GJSys.Fill(i,sy["GJ"][i])
            hist_EfakeSys.Fill(i,sy["efake"][i])
        
            hist_corGJ.Fill(i,cor["GJ"][i])
            hist_corVg.Fill(i,cor["Vg"][i])
            hist_trigger.Fill(i,mcUncertainties["trigger"])
            hist_phSF.Fill(i,mcUncertainties["phSF"])
            hist_lumi.Fill(i,mcUncertainties["lumi"])
        
        hist_GJFit.Fill(0,1.097)
        hist_GJFit.Fill(1,1.0797)
        hist_GJFit.Fill(2,1.0912)
        hist_GJFit.Fill(3,1.335)
        hist_VgFit.Fill(0,1.0722)
        hist_VgFit.Fill(1,1.08045)
        hist_VgFit.Fill(2,1.08787)
        hist_VgFit.Fill(3,1.110)
            
            
        RootOut=rt.TFile(outdir+"SignalTrees/SignalTree_"+sScan+"_"+selection+".root","update")
        
        i=0
        for hist in [hist_SigRate_0,hist_SigRate_1,hist_SigRate_2,hist_SigRate_3]:
            hist.Write("SignalRate_bin"+str(i),rt.TObject.kOverwrite)
            i+=1
        i=0
        for hist in [hist_SigStat_0,hist_SigStat_1,hist_SigStat_2,hist_SigStat_3]:
            hist.Write("SignalStat_bin"+str(i),rt.TObject.kOverwrite)
            i+=1
        i=0
        for hist in [hist_SigSys_0,hist_SigSys_1,hist_SigSys_2,hist_SigSys_3]:
            hist.Write("SignalSys_bin"+str(i),rt.TObject.kOverwrite)
            i+=1
        i=0
        for hist in [hist_SigISR_0,hist_SigISR_1,hist_SigISR_2,hist_SigISR_3]:
            hist.Write("SignalISR_bin"+str(i),rt.TObject.kOverwrite)
            i+=1
        
        hist_GJStat.Write("Stat_GJ",rt.TObject.kOverwrite)
        hist_TTcombStat.Write("Stat_TTcomb",rt.TObject.kOverwrite)
        hist_VgStat.Write("Stat_Vg",rt.TObject.kOverwrite)
        hist_DibosonStat.Write("Stat_Diboson",rt.TObject.kOverwrite)
        hist_EfakeStat.Write("Stat_Efake",rt.TObject.kOverwrite)
        
        hist_TTcombSys.Write("Sys_TTcomb",rt.TObject.kOverwrite)
        hist_DibosonSys.Write("Sys_Diboson",rt.TObject.kOverwrite)
        hist_GJSys.Write("Sys_GJ",rt.TObject.kOverwrite)
        hist_EfakeSys.Write("Sys_Efake",rt.TObject.kOverwrite)
        
        hist_corGJ.Write("Corr_GJ",rt.TObject.kOverwrite)
        hist_corVg.Write("Corr_Vg",rt.TObject.kOverwrite)
        hist_trigger.Write("Trigger",rt.TObject.kOverwrite)
        hist_phSF.Write("PhSF",rt.TObject.kOverwrite)
        hist_lumi.Write("Lumi",rt.TObject.kOverwrite)
        hist_GJFit.Write("Fit_GJ",rt.TObject.kOverwrite)
        hist_VgFit.Write("Fit_Vg",rt.TObject.kOverwrite)
        RootOut.Close()
        """
        
    #~ RootOut=rt.TFile(outdir+"SignalUnc.root","update")
    #~ graph_puUnc.GetHistogram().Write(sScan+"_PuUnc",rt.TObject.kOverwrite)
    #~ graph_jesUnc.GetHistogram().Write(sScan+"_meanJesUnc",rt.TObject.kOverwrite)


if __name__ == '__main__':
    #~ selection="exclusiv"
    #~ selection="inclusiv"
    #~ selection="htgVeto"
    #~ selection="leptonVeto"
    #~ selection="diphotonVeto"
    #~ selection="htgHighVeto"
    #~ selection="htgHighLeptonVeto"
    selection="exclusiv_highHTG"
    #~ selection="leptonDiphotonVeto"
    basedir="../"
    outdir=basedir+"output/"
    signal_scan="signal_scan_"+selection+"_v03D.root"
    #~ rho=-0.0
    rho=-1.0
    #~ fillDatacards(Scan.T5gg,True)
    fillDatacards(Scan.T5Wg,True)
    #~ fillDatacards(Scan.T5Wg_Wg,True)
    #~ fillDatacards(Scan.T5Wg_thirds,True)
    #~ fillDatacards(Scan.T6Wg,True)
    #~ fillDatacards(Scan.T6gg,True)
    #~ fillDatacards(Scan.GGM,True)
    #~ fillDatacards(Scan.TChiWg,False)
    #~ fillDatacards(Scan.TChiNg,False)
    #~ fillDatacards(Scan.GGM_M1_M2,True)
    #~ fillDatacards(Scan.GGM_M1_M3,True)
    #~ fillDatacards(Scan.TChiNg_gg,False)
    #~ fillDatacards(Scan.TChiNg_gz,False)
    #~ fillDatacards(Scan.TChiNg_zz,False)
    #~ fillDatacards(Scan.TChiNg_gg_C1N2,False)


