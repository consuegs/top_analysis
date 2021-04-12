#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"

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

bool matchLepton(TLorentzVector recoLep, TLorentzVector genLep) {
   // ~return (abs(recoLep.DeltaR(genLep))<0.5) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.5); //probably wrong numbers
   return (abs(recoLep.DeltaR(genLep))<0.05) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.1);
}

extern "C"
void run()
{
   std::vector<TString> baselineTriggerNames = {
         "HLT_PFHT300_PFMET110_v",
         "HLT_PFMET120_PFMHT120_IDTight_v",
         "HLT_PFMET170_HBHECleaned_v",
         "HLT_PFMET300_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
         "HLT_MET200_v",
         "HLT_MET200_v",   //Dummies to have same length as 2017_2018
         "HLT_MET200_v",
         "HLT_MET200_v",
         "HLT_MET200_v",
   };
   std::vector<TString> baselineTriggerNames2017_2018 = {
         "HLT_PFMET200_HBHECleaned_v",
         "HLT_PFMET200_HBHE_BeamHaloCleaned_v",
         "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
         "HLT_PFMET120_PFMHT120_IDTight_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
         "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
         "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v",
         "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v",
         "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v",
   };
         
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      if (dss.datasetName.find("MET")==std::string::npos && dss.datasetName.find("TTbar_diLepton")==std::string::npos) {
         std::cout<<dss.datasetName<<" found in running schedule, but module should ony be running with ttbar MC or MET PD"<<std::endl;
         return;
      }
   }
   
   //Define histograms in the following
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   for(TString selection:{"baseline"}){
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
         
         hs2d.addHist(selection+"/"+base+"/ee/pTl1_pTl2"   ,";%pTl1;%pTl2;EventsBIN"           ,40,0,200,40,0,200);
         hs2d.addHist(selection+"/"+base+"/emu/pTlmu_pTle"   ,";%pTlmu;%pTle;EventsBIN"           ,40,0,200,40,0,200);
         hs2d.addHist(selection+"/"+base+"/mumu/pTl1_pTl2"   ,";%pTl1;%pTl2;EventsBIN"           ,40,0,200,40,0,200);
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
      
      //Check if current sample is Run2017
      bool Run2017=(cfg.year=="2017");
      //Check if current sample is Run2017AB
      bool Run2017AB=false;
      if (dss.name.find("Run2017A")!=std::string::npos) Run2017AB=true;
      else if (dss.name.find("Run2017B")!=std::string::npos) Run2017AB=true;
      
      //Check if current sample is Run2018
      bool Run2018=(cfg.year=="2018");
      if (dss.name.find("Run2018")!=std::string::npos) Run2018=true;
      
      //Set correct baselineTriggers
      if (Run2018 || Run2017) baselineTriggerNames=baselineTriggerNames2017_2018;
      
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
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<bool> is_ee   (reader, "ee");
      TTreeReaderValue<bool> is_emu   (reader, "emu");
      TTreeReaderValue<bool> is_mumu   (reader, "mumu");
      TTreeReaderValue<float> mll   (reader, "mll");
      TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
      TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
      TTreeReaderValue<float> btagWeight(reader, "bTagWeight");
      
      TTreeReaderValue<bool> muonTrigg1(reader, cfg.muonTrigg1);
      TTreeReaderValue<bool> muonTrigg2(reader, cfg.muonTrigg2);
      TTreeReaderValue<bool> muonTrigg3(reader, cfg.muonTrigg3);
      TTreeReaderValue<bool> muonTrigg4(reader, cfg.muonTrigg4);
      TTreeReaderValue<bool> singleMuonTrigg1(reader, cfg.singleMuonTrigg1);
      TTreeReaderValue<bool> singleMuonTrigg2(reader, cfg.singleMuonTrigg2);
      TTreeReaderValue<bool> eleTrigg1(reader, cfg.eleTrigg1);
      TTreeReaderValue<bool> eleTrigg2(reader, cfg.eleTrigg2);
      TTreeReaderValue<bool> eleMuTrigg1(reader, cfg.eleMuTrigg1);
      TTreeReaderValue<bool> eleMuTrigg2(reader, cfg.eleMuTrigg2);
      TTreeReaderValue<bool> eleMuTrigg3(reader, cfg.eleMuTrigg3);
      TTreeReaderValue<bool> eleMuTrigg4(reader, cfg.eleMuTrigg4);
      TTreeReaderValue<bool> singleEleTrigg(reader, cfg.singleEleTrigg);
      
      TTreeReaderValue<bool> baselineTrigg1(reader, baselineTriggerNames[0]);
      TTreeReaderValue<bool> baselineTrigg2(reader, baselineTriggerNames[1]);
      TTreeReaderValue<bool> baselineTrigg3(reader, baselineTriggerNames[2]);
      TTreeReaderValue<bool> baselineTrigg4(reader, baselineTriggerNames[3]);
      TTreeReaderValue<bool> baselineTrigg5(reader, baselineTriggerNames[4]);
      TTreeReaderValue<bool> baselineTrigg6(reader, baselineTriggerNames[5]);
      TTreeReaderValue<bool> baselineTrigg7(reader, baselineTriggerNames[6]);
      TTreeReaderValue<bool> baselineTrigg8(reader, baselineTriggerNames[7]);
      TTreeReaderValue<bool> baselineTrigg9(reader, baselineTriggerNames[8]);
      TTreeReaderValue<bool> baselineTrigg10(reader, baselineTriggerNames[9]);
   
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
         float SFWeight=*sf_lep1 * *sf_lep2 * *btagWeight;     //Set combined SF weight
         // ~float SFWeight=*sf_lep1 * *sf_lep2;     //Set combined SF weight
         // ~if(!isData) {
            // ~hs.setFillWeight(fEventWeight*SFWeight);
            // ~hs2d.setFillWeight(fEventWeight*SFWeight);
         // ~}
         // ~else {
            // ~hs.setFillWeight(1.);
            // ~hs2d.setFillWeight(1.);
         // ~}
         
         float met=MET->p.Pt();
         float const met_puppi=MET_Puppi->p.Pt();
         
         //Booleans for reco and zpeak selection
         bool rec_selection=true;
         
         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter)
         std::vector<bool> channel={*is_ee,*is_mumu,*is_emu};
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         int flavor_l1=0;  //1 for electron and 2 for muon
         int flavor_l2=0;
         bool muonLead=true; //Boolean for emu channel
         TString cat="";
         
         rec_selection=selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,cat,muonLead);
               
         std::vector<tree::Jet> cjets;
         std::vector<tree::Jet> BJets;
         std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
         if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end(), [](bool v) { return v; })) rec_selection=false;
         
         //Define trigger selections
         bool baselineTriggers=*baselineTrigg1 || *baselineTrigg2 || *baselineTrigg3 || *baselineTrigg4 || *baselineTrigg5 || *baselineTrigg6 || *baselineTrigg7 || *baselineTrigg8|| *baselineTrigg9 || *baselineTrigg10;
         bool diElectronTriggers=*eleTrigg1 || *eleTrigg2 || *singleEleTrigg;
         bool diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2;
         bool electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         
         if(Run2017)diMuonTriggers=*muonTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2;
         if(Run2016H){
            diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2;
            electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         }
         else if(Run2017AB){
            diMuonTriggers=*muonTrigg1|| *singleMuonTrigg1 || *singleMuonTrigg2;
         }
         
         //Use only events if baseline triggers fired
         if(!baselineTriggers) continue;
         
         //Fill hists
         TString path_cat="ee";
         bool trigger=diElectronTriggers;
         bool doubleTrigger=*eleTrigg1 || *eleTrigg2;
         bool doubleTrigger_DZ=*eleTrigg1 || *eleTrigg2;
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
            if(rec_selection==true && triggerVec[i]){
               hs.fill("baseline/"+base+"/"+path_cat+"/pTl1",p_l1.Pt());
               hs.fill("baseline/"+base+"/"+path_cat+"/pTl2",p_l2.Pt());
               hs.fill("baseline/"+base+"/"+path_cat+"/etal1",p_l1.Eta());
               hs.fill("baseline/"+base+"/"+path_cat+"/etal2",p_l2.Eta());
               if(!*is_emu) hs2d.fill("baseline/"+base+"/"+path_cat+"/pTl1_pTl2",p_l1.Pt(),p_l2.Pt());
               else{
                  if(muonLead) hs2d.fill("baseline/"+base+"/"+path_cat+"/pTlmu_pTle",p_l1.Pt(),p_l2.Pt());
                  else hs2d.fill("baseline/"+base+"/"+path_cat+"/pTlmu_pTle",p_l2.Pt(),p_l1.Pt());
               }
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
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("triggerEff%.1f",cfg.processFraction*100),false);
   hs.saveHistograms(saver_hist,samplesToCombine);
   hs2d.saveHistograms(saver_hist,samplesToCombine);
   
}
