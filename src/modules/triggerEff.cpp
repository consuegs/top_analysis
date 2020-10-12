#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>

Config const &cfg=Config::get();

void saveHistograms(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH1F> &hs, std::vector<TString> const &Samples)
{
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         for (TString sSample: Samples){
            saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
         }       
      }
   }
}
int saveRatios(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH1F> &hs, std::vector<TString> const &Samples)
{  
   if (Samples.size()!=2){
      std::cout<<"Number of samples != 2, no ratios are saved"<<std::endl;
      return 0;
   }
   else{
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            sVar=sPresel+sVar;
            TH1F* temp=hs.getHistogram(sVar,Samples[0]);
            temp->Divide(hs.getHistogram(sVar,Samples[1]));
            saver_hist.save(*temp,"ratios/"+sVar);    
         }
      }
      return 0;
   }
}
void saveHistograms2D(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH2F> &hs, std::vector<TString> const &Samples)
{
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         for (TString sSample: Samples){
            saver_hist.save(*hs.getHistogram(sVar,sSample),sVar+"/"+sSample);
         }       
      }
   }
}

bool matchLepton(TLorentzVector recoLep, TLorentzVector genLep) {
   // ~return (abs(recoLep.DeltaR(genLep))<0.5) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.5); //probably wrong numbers
   return (abs(recoLep.DeltaR(genLep))<0.05) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.1);
}

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   //Define histograms in the following
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   for(TString selection:{"baseline","Zpeak","Zpeak_noJetRequ"}){
      for(TString base:{"referenceTrigg","analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         hs.addHist(selection+"/"+base+"/ee/pTl1"   ,";%pTl1;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/emu/pTl1"   ,";%pTl1;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/mumu/pTl1"   ,";%pTl1;EventsBIN"           ,10,0,200);
         
         hs.addHist(selection+"/"+base+"/ee/pTl2"   ,";%pTl2;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/emu/pTl2"   ,";%pTl2;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/mumu/pTl2"   ,";%pTl2;EventsBIN"           ,10,0,200);
         
         hs.addHist(selection+"/"+base+"/ee/etal1"   ,";#eta_{l1};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/emu/etal1"   ,";#eta_{l1};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/mumu/etal1"   ,";#eta_{l1};EventsBIN"           ,15,-2.4,2.4);
         
         hs.addHist(selection+"/"+base+"/ee/etal2"   ,";#eta_{l2};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/emu/etal2"   ,";#eta_{l2};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/mumu/etal2"   ,";#eta_{l2};EventsBIN"           ,15,-2.4,2.4);
      }
   }

   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      
      hs.setCurrentSample(dss.name);
      hs2d.setCurrentSample(dss.name);
      
      bool const isData=dss.isData;
      
      //Check if current sample is Run2016H
      bool Run2016H=false;
      if (dss.datasetName.find("Run2016H")!=std::string::npos) Run2016H=true;
      
      //Check if current sample is DY MC
      bool DY_MC=false;
      if (dss.datasetName.find("DY")!=std::string::npos) DY_MC=true;

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<UInt_t> runNo(reader, "runNo");
      TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
      TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
      TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<bool> is_ee   (reader, "ee");
      TTreeReaderValue<bool> is_emu   (reader, "emu");
      TTreeReaderValue<bool> is_mumu   (reader, "mumu");
      TTreeReaderValue<float> mll   (reader, "mll");
      TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
      TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
      TTreeReaderValue<float> btagWeight(reader, "bTagWeight");
      TTreeReaderValue<bool> muonTrigg1(reader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v");
      TTreeReaderValue<bool> muonTrigg2(reader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v");
      TTreeReaderValue<bool> muonTrigg3(reader, "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v");
      TTreeReaderValue<bool> muonTrigg4(reader, "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v");
      TTreeReaderValue<bool> singleMuonTrigg1(reader, "HLT_IsoMu24_v");
      TTreeReaderValue<bool> singleMuonTrigg2(reader, "HLT_IsoTkMu24_v");
      TTreeReaderValue<bool> eleTrigg(reader, "HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TTreeReaderValue<bool> eleMuTrigg1(reader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TTreeReaderValue<bool> eleMuTrigg2(reader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v");
      TTreeReaderValue<bool> eleMuTrigg3(reader, "HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v");
      TTreeReaderValue<bool> eleMuTrigg4(reader, "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v");
      TTreeReaderValue<bool> singleEleTrigg(reader, "HLT_Ele27_WPTight_Gsf_v");
      TTreeReaderValue<bool> baselineTrigg1(reader, "HLT_PFHT300_PFMET110_v");
      TTreeReaderValue<bool> baselineTrigg2(reader, "HLT_PFMET120_PFMHT120_IDTight_v");
      TTreeReaderValue<bool> baselineTrigg3(reader, "HLT_PFMET170_HBHECleaned_v");
      TTreeReaderValue<bool> baselineTrigg4(reader, "HLT_PFMET300_v");
      TTreeReaderValue<bool> baselineTrigg5(reader, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v");
      TTreeReaderValue<bool> baselineTrigg6(reader, "HLT_MET200_v");
   
      int iEv=0;
      int processEvents=cfg.processFraction*dss.entries;
      // ~if(DY_MC) processEvents*=0.1;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }
         
         
         float fEventWeight=*w_pu * *w_mc;     //Set event weight 
         // ~float SFWeight=*sf_lep1 * *sf_lep2 * *btagWeight;     //Set combined SF weight
         float SFWeight=*sf_lep1 * *sf_lep2;     //Set combined SF weight
         // ~if(!isData) {
            // ~hs.setFillWeight(fEventWeight*SFWeight);
            // ~hs2d.setFillWeight(fEventWeight*SFWeight);
         // ~}
         // ~else {
            // ~hs.setFillWeight(1.);
            // ~hs2d.setFillWeight(1.);
         // ~}
         
         float met=MET->p.Pt();
         
         //Booleans for reco and zpeak selection
         bool rec_selection=true;
         bool mllRequ=false;
         bool metRequ=true;
         
         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         float p_l1_miniIso;
         float p_l2_miniIso;
         
         int charge_l1;
         int charge_l2;
         
         if (*is_ee){
            if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) rec_selection=false; //currently double check since trees only have tight leptons!!
            if(abs((*electrons)[0].etaSC)>2.4 || abs((*electrons)[1].etaSC>2.4)) rec_selection=false; //To use same region as for muons, cut on supercluster eta
            p_l1=(*electrons)[0].p*(*electrons)[0].corr;
            p_l2=(*electrons)[1].p*(*electrons)[1].corr;
            p_l1_miniIso=(*electrons)[0].PFminiIso;
            p_l2_miniIso=(*electrons)[1].PFminiIso;
            charge_l1=(*electrons)[0].charge;
            charge_l2=(*electrons)[1].charge;
         }
         else if (*is_mumu){
            if(!(*muons)[0].isTight || !(*muons)[1].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 || (*muons)[1].rIso>0.15) rec_selection=false;
            if(abs((*muons)[0].p.Eta())>2.4 || abs((*muons)[1].p.Eta())>2.4) rec_selection=false;
            p_l1=(*muons)[0].p*(*muons)[0].rochesterCorrection;
            p_l2=(*muons)[1].p*(*muons)[1].rochesterCorrection;
            p_l1_miniIso=(*muons)[0].PFminiIso;
            p_l2_miniIso=(*muons)[1].PFminiIso;
            charge_l1=(*muons)[0].charge;
            charge_l2=(*muons)[1].charge;
         }
         else if (*is_emu){
            if(!(*muons)[0].isTight || !(*electrons)[0].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 ) rec_selection=false;
            if(abs((*muons)[0].p.Eta())>2.4) rec_selection=false;
            if(abs((*electrons)[0].etaSC>2.4) ) rec_selection=false;
            if ((*muons)[0].p.Pt()>(*electrons)[0].p.Pt()){
               p_l1=(*muons)[0].p*(*muons)[0].rochesterCorrection;
               p_l2=(*electrons)[0].p*(*electrons)[0].corr;
               p_l1_miniIso=(*muons)[0].PFminiIso;
               p_l2_miniIso=(*electrons)[0].PFminiIso;
               charge_l1=(*muons)[0].charge;
               charge_l2=(*electrons)[0].charge;
            }
            else {
               p_l1=(*electrons)[0].p*(*electrons)[0].corr;
               p_l2=(*muons)[0].p*(*muons)[0].rochesterCorrection;
               p_l1_miniIso=(*electrons)[0].PFminiIso;
               p_l2_miniIso=(*muons)[0].PFminiIso;
               charge_l1=(*electrons)[0].charge;
               charge_l2=(*muons)[0].charge;
            }
         }
         
         if (p_l1.Pt()<25 || p_l2.Pt()<20) rec_selection=false; //eta cuts already done in TreeWriter
         if (*mll<20) rec_selection=false;
         if ((*is_ee || *is_mumu) && *mll<106 && *mll>76) mllRequ=true;
         if ((*is_ee || *is_mumu) && met<40) metRequ=false;
         
         
         bool jetRequirement=true;
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         if (cjets.size()<2) jetRequirement=false;
         
         std::vector<tree::Jet> lJets;
         for (tree::Jet const &jet : *jets) {
            if (jet.isLoose) lJets.push_back(jet);
         }
         
         bool bTag=false;
         std::vector<tree::Jet> BJets;
         for (tree::Jet const &jet : cjets) {
            // ~if (jet.bTagCSVv2>0.5426) {      //Loose working point for CSVv2 (Should be replaced in the future by deep CSV!!!)
            if (jet.bTagDeepCSV>0.2217) {      //Loose working point for deepCSV
               bTag=true;
               BJets.push_back(jet);
            }
         }
         
         //Define trigger selections
         bool baselineTriggers=*baselineTrigg1 || *baselineTrigg2 || *baselineTrigg3 || *baselineTrigg4 || *baselineTrigg5 || *baselineTrigg6;
         bool diElectronTriggers=*eleTrigg || *singleEleTrigg;
         bool diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2;
         bool electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         if(Run2016H){
            diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2;
            electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         }
         
         //Use only events if baseline triggers fired
         if(!baselineTriggers) continue;
         
         //Fill hists
         TString path_cat="ee";
         bool trigger=diElectronTriggers;
         bool doubleTrigger=*eleTrigg;
         bool doubleTrigger_DZ=*eleTrigg;
         bool singleTrigger=*singleEleTrigg;
         if (*is_emu) {
            path_cat="emu";
            trigger=electronMuonTriggers;
            doubleTrigger=*eleMuTrigg3 || *eleMuTrigg4;
            doubleTrigger_DZ=*eleMuTrigg1 || *eleMuTrigg2;
            singleTrigger=*singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         }
         else if (*is_mumu) {
            path_cat="mumu";
            trigger=diMuonTriggers;
            doubleTrigger=*muonTrigg3 || *muonTrigg4;
            doubleTrigger_DZ=*muonTrigg1 || *muonTrigg2;
            singleTrigger=*singleMuonTrigg1 || *singleMuonTrigg2;
         }
         
         std::vector<bool> triggerVec={true,trigger,doubleTrigger_DZ,doubleTrigger,singleTrigger};
         int i=0;
         for(TString base:{"referenceTrigg","analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
            if(rec_selection==true && mllRequ==false && jetRequirement && bTag && metRequ && triggerVec[i]){
               hs.fill("baseline/"+base+"/"+path_cat+"/pTl1",p_l1.Pt());
               hs.fill("baseline/"+base+"/"+path_cat+"/pTl2",p_l2.Pt());
               hs.fill("baseline/"+base+"/"+path_cat+"/etal1",p_l1.Eta());
               hs.fill("baseline/"+base+"/"+path_cat+"/etal2",p_l2.Eta());
            }
            
            if(rec_selection==true && mllRequ==true && jetRequirement && bTag && metRequ && triggerVec[i]){
               hs.fill("Zpeak/"+base+"/"+path_cat+"/pTl1",p_l1.Pt());
               hs.fill("Zpeak/"+base+"/"+path_cat+"/pTl2",p_l2.Pt());
               hs.fill("Zpeak/"+base+"/"+path_cat+"/etal1",p_l1.Eta());
               hs.fill("Zpeak/"+base+"/"+path_cat+"/etal2",p_l2.Eta());
            }
            
            if(rec_selection && mllRequ && triggerVec[i]){
               hs.fill("Zpeak_noJetRequ/"+base+"/"+path_cat+"/pTl1",p_l1.Pt());
               hs.fill("Zpeak_noJetRequ/"+base+"/"+path_cat+"/pTl2",p_l2.Pt());
               hs.fill("Zpeak_noJetRequ/"+base+"/"+path_cat+"/etal1",p_l1.Eta());
               hs.fill("Zpeak_noJetRequ/"+base+"/"+path_cat+"/etal2",p_l2.Eta());
            }
            i++;
         }
         
      }// evt loop
      io::log<<"";
      
      // ~hs.scaleLumi();
      // ~hs2d.scaleLumi();
      hs.mergeOverflow();
      hs2d.mergeOverflow();
      file.Close();
      
   } // dataset loop
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars;
   for(TString selection:{"baseline","Zpeak","Zpeak_noJetRequ"}){
      for(TString base:{"referenceTrigg","analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         for(TString channel:{"ee","mumu","emu"}){
            msPresel_vVars.insert(std::pair<TString,std::vector<TString>>(selection+"/"+base+"/"+channel+"/",{"pTl1","pTl2","etal1","etal2"}));
         }
      }
   }
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("triggerEff%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs,samplesToCombine);
   
}
