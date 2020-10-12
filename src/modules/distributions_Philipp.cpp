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
   hs.addHist("baseline/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("baseline/ee/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("baseline/emu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("baseline/mumu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   
   hs.addHist("baseline/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline_noTau/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau_tightRIsoMu/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("baseline_noTau_matchedLep/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau_matchedLep/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau_matchedLep/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("baseline_noTau/ee/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau/emu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("baseline_noTau/mumu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   
   hs.addHist("baseline_noTau/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline_noTau/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline_noTau/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline_noTau_trigg_lastBin/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_lastBin/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_lastBin/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   
   hs.addHist("baseline_noTau_trigg_secondlastBin/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_secondlastBin/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_secondlastBin/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   
   hs.addHist("baseline_noTau_trigg_firstBin/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_firstBin/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_firstBin/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   
   hs.addHist("baseline_noTau_trigg_firstFirstBin/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_firstFirstBin/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   hs.addHist("baseline_noTau_trigg_firstFirstBin/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,1,0,3.2);
   
   hs.addHist("Zpeak/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("Zpeak/ee/genHT"   ,";genHT;EventsBIN"           ,300,0,3000);
   hs.addHist("Zpeak/emu/genHT"   ,";genHT;EventsBIN"           ,300,0,3000);
   hs.addHist("Zpeak/mumu/genHT"   ,";genHT;EventsBIN"           ,300,0,3000);
   
   hs.addHist("Zpeak/ee/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak/emu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak/mumu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   
   hs.addHist("Zpeak/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("Zpeak/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("Zpeak/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
  
   hs.addHist("Zpeak/ee/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/emu/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/mumu/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak/ee/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/emu/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/mumu/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak/ee/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/emu/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/mumu/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak/ee/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/emu/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak/mumu/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak/ee/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak/emu/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak/mumu/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   
   hs.addHist("Zpeak/ee/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak/emu/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak/mumu/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   
   hs.addHist("Zpeak/ee/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak/emu/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak/mumu/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   
   hs.addHist("Zpeak/ee/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak/emu/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak/mumu/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   
   hs.addHist("Zpeak/ee/phil1"   ,";phil1;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak/emu/phil1"   ,";phil1;EventsBIN"          ,30,-3.14,3.14);
   hs.addHist("Zpeak/mumu/phil1"   ,";phil1;EventsBIN"           ,30,-3.14,3.14);
   
   hs.addHist("Zpeak/ee/phil2"   ,";phil2;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak/emu/phil2"   ,";phil2;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak/mumu/phil2"   ,";phil2;EventsBIN"          ,30,-3.14,3.14);
   
   hs.addHist("Zpeak/ee/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   hs.addHist("Zpeak/emu/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   hs.addHist("Zpeak/mumu/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   
   hs.addHist("Zpeak/ee/mll"   ,";mll;EventsBIN"           ,30,76,106);
   hs.addHist("Zpeak/emu/mll"   ,";mll;EventsBIN"           ,30,76,106);
   hs.addHist("Zpeak/mumu/mll"   ,";mll;EventsBIN"           ,30,76,106);
   
   hs.addHist("Zpeak/ee/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   hs.addHist("Zpeak/emu/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   hs.addHist("Zpeak/mumu/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   
   hs.addHist("Zpeak_noTrigg/ee/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   hs.addHist("Zpeak_noTrigg/emu/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   hs.addHist("Zpeak_noTrigg/mumu/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   
   hs.addHist("Zpeak/ee/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak/emu/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak/mumu/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   
   hs.addHist("Zpeak/ee/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak/emu/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak/mumu/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   
   hs.addHist("Zpeak/ee/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   hs.addHist("Zpeak/emu/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   hs.addHist("Zpeak/mumu/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   
   hs.addHist("Zpeak/ee/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak/emu/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak/mumu/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   
   hs.addHist("Zpeak/ee/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak/emu/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak/mumu/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   
   hs.addHist("Zpeak/ee/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak/emu/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak/mumu/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   
   hs.addHist("Zpeak/ee/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak/emu/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak/mumu/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   
   hs.addHist("Zpeak/ee/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/emu/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/mumu/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   
   hs.addHist("Zpeak/ee/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/emu/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/mumu/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   
   hs.addHist("Zpeak/ee/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/emu/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak/mumu/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/ee/met"   ,";%MET;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/met"   ,";%MET;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/met"   ,";%MET;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak_noJetRequ/ee/PuppiMet"   ,";PuppiMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,100,0,200);
  
   hs.addHist("Zpeak_noJetRequ/ee/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/CaloMet"   ,";CaloMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak_noJetRequ/ee/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/RawMet"   ,";RawMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak_noJetRequ/ee/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/DeepMet"   ,";DeepMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak_noJetRequ/ee/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/emu/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   hs.addHist("Zpeak_noJetRequ/mumu/NoHFMet"   ,";NoHFMet;EventsBIN"           ,100,0,200);
   
   hs.addHist("Zpeak_noJetRequ/ee/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak_noJetRequ/emu/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak_noJetRequ/mumu/pTl1"   ,";pTl1;EventsBIN"           ,60,0,300);
   
   hs.addHist("Zpeak_noJetRequ/ee/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak_noJetRequ/emu/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   hs.addHist("Zpeak_noJetRequ/mumu/pTl2"   ,";pTl2;EventsBIN"           ,60,0,300);
   
   hs.addHist("Zpeak_noJetRequ/ee/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak_noJetRequ/emu/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak_noJetRequ/mumu/etal1"   ,";etal1;EventsBIN"           ,30,-2.4,2.4);
   
   hs.addHist("Zpeak_noJetRequ/ee/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak_noJetRequ/emu/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   hs.addHist("Zpeak_noJetRequ/mumu/etal2"   ,";etal2;EventsBIN"           ,30,-2.4,2.4);
   
   hs.addHist("Zpeak_noJetRequ/ee/phil1"   ,";phil1;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak_noJetRequ/emu/phil1"   ,";phil1;EventsBIN"          ,30,-3.14,3.14);
   hs.addHist("Zpeak_noJetRequ/mumu/phil1"   ,";phil1;EventsBIN"           ,30,-3.14,3.14);
   
   hs.addHist("Zpeak_noJetRequ/ee/phil2"   ,";phil2;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak_noJetRequ/emu/phil2"   ,";phil2;EventsBIN"           ,30,-3.14,3.14);
   hs.addHist("Zpeak_noJetRequ/mumu/phil2"   ,";phil2;EventsBIN"          ,30,-3.14,3.14);
   
   hs.addHist("Zpeak_noJetRequ/ee/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   hs.addHist("Zpeak_noJetRequ/emu/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   hs.addHist("Zpeak_noJetRequ/mumu/dPhill"   ,";dPhill;EventsBIN"           ,30,0,3.14);
   
   hs.addHist("Zpeak_noJetRequ/ee/mll"   ,";mll;EventsBIN"           ,30,76,106);
   hs.addHist("Zpeak_noJetRequ/emu/mll"   ,";mll;EventsBIN"           ,30,76,106);
   hs.addHist("Zpeak_noJetRequ/mumu/mll"   ,";mll;EventsBIN"           ,30,76,106);
   
   hs.addHist("Zpeak_noJetRequ/ee/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   hs.addHist("Zpeak_noJetRequ/emu/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   hs.addHist("Zpeak_noJetRequ/mumu/ZpT"   ,";ZpT;EventsBIN"           ,50,0,100);
   
   hs.addHist("Zpeak_noJetRequ_noTrigg/ee/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   hs.addHist("Zpeak_noJetRequ_noTrigg/emu/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   hs.addHist("Zpeak_noJetRequ_noTrigg/mumu/pTl2"   ,";pTl2;EventsBIN"           ,100,0,100);
   
   hs.addHist("Zpeak_noJetRequ/ee/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak_noJetRequ/emu/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak_noJetRequ/mumu/nJets"   ,";nJets;EventsBIN"           ,7,0,7);
   
   hs.addHist("Zpeak_noJetRequ/ee/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak_noJetRequ/emu/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   hs.addHist("Zpeak_noJetRequ/mumu/nLooseJets"   ,";nLooseJets;EventsBIN"           ,7,0,7);
   
   hs.addHist("Zpeak_noJetRequ/ee/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   hs.addHist("Zpeak_noJetRequ/emu/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   hs.addHist("Zpeak_noJetRequ/mumu/looseJetpT1"   ,";looseJetpT1;EventsBIN"           ,100,0,600);
   
   hs.addHist("Zpeak_noJetRequ/ee/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak_noJetRequ/emu/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak_noJetRequ/mumu/nGoodVertices"   ,";nGoodVertices;EventsBIN"           ,100,0,60);
   
   hs.addHist("Zpeak_noJetRequ/ee/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak_noJetRequ/emu/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   hs.addHist("Zpeak_noJetRequ/mumu/nPV"   ,";nPV;EventsBIN"           ,100,0,60);
   
   hs.addHist("Zpeak_noJetRequ/ee/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak_noJetRequ/emu/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak_noJetRequ/mumu/PFminiIsoL1"   ,";PFminiIsoL1;EventsBIN"           ,100,0,0.2);
   
   hs.addHist("Zpeak_noJetRequ/ee/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak_noJetRequ/emu/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   hs.addHist("Zpeak_noJetRequ/mumu/PFminiIsoL2"   ,";PFminiIsoL2;EventsBIN"           ,100,0,0.2);
   
   hs.addHist("Zpeak_noJetRequ/ee/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/emu/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/mumu/ChargeL1"   ,";ChargeL1;EventsBIN"           ,3,-1,1);
   
   hs.addHist("Zpeak_noJetRequ/ee/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/emu/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/mumu/ChargeL2"   ,";ChargeL2;EventsBIN"           ,3,-1,1);
   
   hs.addHist("Zpeak_noJetRequ/ee/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/emu/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   hs.addHist("Zpeak_noJetRequ/mumu/ChargeProduct"   ,";ChargeProduct;EventsBIN"           ,3,-1,1);
   
   hs.addHist("Zpeak_noTau/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau_tightRIsoMu/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("Zpeak_noTau_matchedLep/ee/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau_matchedLep/emu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau_matchedLep/mumu/met"   ,";%MET;EventsBIN"           ,60,0,600);
   
   hs.addHist("Zpeak_noTau/ee/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau/emu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   hs.addHist("Zpeak_noTau/mumu/PuppiMet"   ,";PuppiMet;EventsBIN"           ,60,0,600);
   
   hs.addHist("Zpeak_noTau/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("Zpeak_noTau/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   hs.addHist("Zpeak_noTau/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
   
   //Define 2D histograms in the following
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   // ~hs2d.addHist("baseline_noTau_trigg/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   // ~hs2d.addHist("baseline_noTau_trigg/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   // ~hs2d.addHist("baseline_noTau_trigg/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline_noTau_trigg/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,1000,0,1000,320,0,3.2);
   hs2d.addHist("baseline_noTau_trigg/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,1000,0,1000,320,0,3.2);
   hs2d.addHist("baseline_noTau_trigg/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,1000,0,1000,320,0,3.2);

   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      
      hs.setCurrentSample(dss.name);
      hs2d.setCurrentSample(dss.name);
      
      bool const isData=dss.isData;
      
      //Check if current sample is SUSY scenario
      bool SUSY_scen=false;
      if (dss.name.find("SMS")!=std::string::npos) SUSY_scen=true;
      
      //Check if current sample is DY Incl scenario
      bool DYIncl=false;
      if (dss.datasetName=="DrellYan") DYIncl=true;
      
      //Check if current sample is Run2016H
      bool Run2016H=false;
      if (dss.datasetName.find("Run2016H")!=std::string::npos) Run2016H=true;
      
      //Check which PD is currently used
      bool SingleElectron=false;
      bool SingleMuon=false;
      bool DoubleEG=false;
      bool DoubleMuon=false;
      if (dss.datasetName.find("SingleElectron")!=std::string::npos) SingleElectron=true;
      else if (dss.datasetName.find("SingleMuon")!=std::string::npos) SingleMuon=true;
      else if (dss.datasetName.find("DoubleEG")!=std::string::npos) DoubleEG=true;
      else if (dss.datasetName.find("DoubleMuon")!=std::string::npos) DoubleMuon=true;

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
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<std::vector<tree::GenParticle>> DYgenLep(reader, "DYgenLep");
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
      TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<tree::MET> MET_Deep(reader, "metDeep");
      TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
      TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
      TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
      TTreeReaderValue<int> n_Interactions(reader, "true_nPV");
      TTreeReaderValue<float> HTgen(reader, "genHt");
      TTreeReaderValue<bool> is_ee   (reader, "ee");
      TTreeReaderValue<bool> is_emu   (reader, "emu");
      TTreeReaderValue<bool> is_mumu   (reader, "mumu");
      TTreeReaderValue<bool> addLepton   (reader, "addLepton");
      TTreeReaderValue<float> mll   (reader, "mll");
      TTreeReaderValue<float> mt2   (reader, "MT2");
      TTreeReaderValue<float> genMT2   (reader, "genMT2");
      TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
      TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
      TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
      TTreeReaderValue<float> btagWeight(reader, "bTagWeight");
      TTreeReaderValue<TLorentzVector> genTop(reader, "pseudoTop");
      TTreeReaderValue<TLorentzVector> genAntiTop(reader, "pseudoAntiTop");
      TTreeReaderValue<TLorentzVector> genLepton(reader, "pseudoLepton");
      TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "pseudoAntiLepton");
      TTreeReaderValue<TLorentzVector> genTau(reader, "pseudoTau");
      TTreeReaderValue<TLorentzVector> genAntiTau(reader, "pseudoAntiTau");
      TTreeReaderValue<int> genLeptonPdgId(reader, "pseudoLeptonPdgId");
      TTreeReaderValue<int> genAntiLeptonPdgId(reader, "pseudoAntiLeptonPdgId");
      TTreeReaderValue<TLorentzVector> genB(reader, "pseudoBJet");
      TTreeReaderValue<TLorentzVector> genAntiB(reader, "pseudoAntiBJet");
      TTreeReaderValue<TLorentzVector> genNeutrino(reader, "pseudoNeutrino");
      TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "pseudoAntiNeutrino");
      TTreeReaderValue<TLorentzVector> genWMinus(reader, "pseudoWMinus");
      TTreeReaderValue<TLorentzVector> genWPlus(reader, "pseudoWPlus");
      TTreeReaderValue<int> genDecayMode_pseudo(reader, "ttbarPseudoDecayMode");
      // ~TTreeReaderValue<TLorentzVector> genTop(reader, "genTop");    //Alternative genParticles, which are not based on RIVET
      // ~TTreeReaderValue<TLorentzVector> genAntiTop(reader, "genAntiTop");
      // ~TTreeReaderValue<TLorentzVector> genLepton(reader, "genLepton");
      // ~TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "genAntiLepton");
      // ~TTreeReaderValue<TLorentzVector> genTau(reader, "genTau");
      // ~TTreeReaderValue<TLorentzVector> genAntiTau(reader, "genAntiTau");
      // ~TTreeReaderValue<int> genLeptonPdgId(reader, "genLeptonPdgId");
      // ~TTreeReaderValue<int> genAntiLeptonPdgId(reader, "genAntiLeptonPdgId");
      // ~TTreeReaderValue<TLorentzVector> genB(reader, "genB");
      // ~TTreeReaderValue<TLorentzVector> genAntiB(reader, "genAntiB");
      // ~TTreeReaderValue<TLorentzVector> genNeutrino(reader, "genNeutrino");
      // ~TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "genAntiNeutrino");
      // ~TTreeReaderValue<TLorentzVector> genWMinus(reader, "genWMinus");
      // ~TTreeReaderValue<TLorentzVector> genWPlus(reader, "genWPlus");
      TTreeReaderValue<int> genDecayMode(reader, "ttbarDecayMode");
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
      TTreeReaderValue<int> nGoodVertices(reader, "nGoodVertices");
      TTreeReaderValue<int> nPV(reader, "nPV");
   
   
      int iEv=0;
      int processEvents=cfg.processFraction*dss.entries;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }
         
         // ~if(DYIncl && *HTgen>100) continue;
         
         float fEventWeight=*w_pu * *w_mc;     //Set event weight 
         // ~float SFWeight=*sf_lep1 * *sf_lep2 * *btagWeight;     //Set combined SF weight
         float SFWeight=*sf_lep1 * *sf_lep2;     //Set combined SF weight
         if(!isData) {
            hs.setFillWeight(fEventWeight*SFWeight);
            hs2d.setFillWeight(fEventWeight*SFWeight);
         }
         else {
            hs.setFillWeight(1.);
            hs2d.setFillWeight(1.);
         }
         
         float met=MET->p.Pt();
         
         //Booleans for reco and zpeak selection
         bool rec_selection=true;
         bool mllRequ=false;
         bool metRequ=true;
         bool tightRIso_mu=true;
         
         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         float p_l1_miniIso;
         float p_l2_miniIso;
         
         int charge_l1;
         int charge_l2;
         
         if (*is_ee){
            if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) rec_selection=false; //currently double check since trees only have tight leptons!!
            if(abs((*electrons)[0].etaSC)>2.4 || abs((*electrons)[1].etaSC)>2.4) rec_selection=false; //To use same region as for muons, cut on supercluster eta
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
            if(abs((*electrons)[0].etaSC)>2.4 ) rec_selection=false;
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
         // ~if (!bTag) rec_selection=false; // end reco baseline selection
         
         bool tauEvent = false;
         for (auto genLep : *DYgenLep){
            if (genLep.pdgId==15) tauEvent = true;
         }
         
         bool matchedLeptons = false;
         // ~if(!isData) {
            // ~if (matchLepton(p_l1,(*DYgenLep)[0].p) && matchLepton(p_l2,(*DYgenLep)[1].p)) matchedLeptons = true;
            // ~if (matchLepton(p_l2,(*DYgenLep)[0].p) && matchLepton(p_l1,(*DYgenLep)[1].p)) matchedLeptons = true;
         // ~}
         
         // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
         TLorentzVector neutrinoPair(0,0,0,0);
         neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
         
         //Calculate MET before parton shower and hadronization for DM scenario and SUSY scenarios if gen level selection is true
         for (auto const &genParticle : *genParticles){
            if((!SUSY_scen && abs(genParticle.pdgId)>100000) || (SUSY_scen && abs(genParticle.pdgId)==1000022)){
               neutrinoPair+=genParticle.p;
            }
         }
         
         //Get DeltaPhi between MET (or genMet or neutrino pT) and nearest Lepton
         float dPhiMETnearLep=4;
         float dPhiMETnearLepPuppi=4;
         float dPhiMETnearLepDeep=4;
         for (TLorentzVector const lep : {p_l1,p_l2}){
            const float dPhi=MET->p.DeltaPhi(lep);
            const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
            const float dPhi_Deep=MET_Deep->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
               dPhiMETnearLep=dPhi;
            }
            if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLepPuppi)) {
               dPhiMETnearLepPuppi=dPhi_Puppi;
            }
            if (std::abs(dPhi_Deep) < std::abs(dPhiMETnearLepDeep)) {
               dPhiMETnearLepDeep=dPhi_Deep;
            }
         }
         float dPhigenMETnearLep=4;
         for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
            const float dPhi=GENMET->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhigenMETnearLep)) dPhigenMETnearLep=dPhi;
         }
         float dPhiPtNunearLep=4;
         for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
            const float dPhi=neutrinoPair.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiPtNunearLep)) {
               dPhiPtNunearLep=dPhi;
            }
         }
         
         //Calculate genMET for DM scenario and SUSY scenarios (!!!Currently only neutrinos with pT above 10 GeV)
         TLorentzVector DMgenMET(0,0,0,0);
         for (auto const &genParticle : *genParticles){
            if (abs(genParticle.pdgId)==12 || abs(genParticle.pdgId)==14 || abs(genParticle.pdgId)==16){
               DMgenMET+=genParticle.p;
            }
            else if((!SUSY_scen && abs(genParticle.pdgId)>100000) || (SUSY_scen && abs(genParticle.pdgId)==1000022)){
               DMgenMET+=genParticle.p;
            }
         }
         
         //Define trigger selections and channel
         bool diElectronTriggers=*eleTrigg || *singleEleTrigg;
         bool diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2;
         bool electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         
         bool trigger=true;
         bool trigger_veto=false;
         TString path_cat="";
         if(!isData){
            path_cat="ee";
            trigger=diElectronTriggers;
            if (*is_emu) {
               path_cat="emu";
               trigger=electronMuonTriggers;
            }
            else if (*is_mumu) {
               path_cat="mumu";
               trigger=diMuonTriggers;
            }
         }
         else{    //Here no emu since actually only same sign is considered for Zpeak
            if(SingleElectron && *is_ee){
               trigger=*singleEleTrigg;
               trigger_veto=*eleTrigg;
               path_cat="ee";
            }
            else if(DoubleEG && *is_ee){
               trigger=*eleTrigg;
               trigger_veto=false;
               path_cat="ee";
            }
            else if(SingleMuon && *is_mumu){
               trigger=*singleMuonTrigg1 || *singleMuonTrigg2;
               if (Run2016H) trigger_veto=*muonTrigg3 || *muonTrigg4;
               else trigger_veto=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4;
               path_cat="mumu";
            }
            else if(DoubleMuon && *is_mumu){
               if (Run2016H) trigger=*muonTrigg3 || *muonTrigg4;
               else trigger=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4;
               trigger_veto=false;
               path_cat="mumu";
            }
            else continue;
         }
            
         
         if(rec_selection==true && mllRequ==false && jetRequirement && bTag && metRequ && trigger && !trigger_veto){
            hs.fill("baseline/"+path_cat+"/met",met);
            hs.fill("baseline/"+path_cat+"/PuppiMet",MET_Puppi->p.Pt());
            hs.fill("baseline/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
            hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
            if(!tauEvent) {
               hs.fill("baseline_noTau/"+path_cat+"/met",met);
               hs.fill("baseline_noTau/"+path_cat+"/PuppiMet",MET_Puppi->p.Pt());
               hs.fill("baseline_noTau/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
               if(matchedLeptons) hs.fill("baseline_noTau_matchedLep/"+path_cat+"/met",met);
               // ~else if (*is_mumu && abs(dPhiMETnearLep)<0.35 && met>140)
               if(tightRIso_mu && *is_mumu) hs.fill("baseline_noTau_tightRIsoMu/"+path_cat+"/met",met);
               if (trigger) {
                  hs2d.fill("baseline_noTau_trigg/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
                  
                  if (met>230 && abs(dPhiMETnearLep)>2.27) hs.fill("baseline_noTau_trigg_lastBin/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
                  if (met>=195 && met<=230 && abs(dPhiMETnearLep)>=2.27) hs.fill("baseline_noTau_trigg_secondlastBin/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
                  if (met <60 && abs(dPhiMETnearLep)>2.27) hs.fill("baseline_noTau_trigg_firstBin/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
                  if (met <60 && abs(dPhiMETnearLep)<0.35) hs.fill("baseline_noTau_trigg_firstFirstBin/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
               }
            }
         }
         
         if(rec_selection==true && mllRequ==true && jetRequirement && bTag && metRequ && trigger && !trigger_veto){
            hs.fill("Zpeak/"+path_cat+"/met",met);
            hs.fill("Zpeak/"+path_cat+"/PuppiMet",MET_Puppi->p.Pt());
            hs.fill("Zpeak/"+path_cat+"/genHT",*HTgen);
            hs.fill("Zpeak/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
            hs.fill("Zpeak/"+path_cat+"/CaloMet",MET_Calo->p.Pt());
            hs.fill("Zpeak/"+path_cat+"/DeepMet",MET_Deep->p.Pt());
            hs.fill("Zpeak/"+path_cat+"/RawMet",MET_Raw->p.Pt());
            hs.fill("Zpeak/"+path_cat+"/NoHFMet",MET_NoHF->p.Pt());
            hs.fill("Zpeak/"+path_cat+"/pTl1",p_l1.Pt());
            hs.fill("Zpeak/"+path_cat+"/pTl2",p_l2.Pt());
            hs.fill("Zpeak/"+path_cat+"/etal1",(*is_ee)? (*electrons)[0].etaSC : p_l1.Eta());
            hs.fill("Zpeak/"+path_cat+"/etal2",(*is_ee)? (*electrons)[1].etaSC : p_l2.Eta());
            hs.fill("Zpeak/"+path_cat+"/phil1",p_l1.Phi());
            hs.fill("Zpeak/"+path_cat+"/phil2",p_l2.Phi());
            hs.fill("Zpeak/"+path_cat+"/ZpT",(p_l1+p_l2).Pt());
            hs.fill("Zpeak/"+path_cat+"/mll",(p_l1+p_l2).M());
            hs.fill("Zpeak/"+path_cat+"/dPhill",abs(p_l1.DeltaPhi(p_l2)));
            hs.fill("Zpeak/"+path_cat+"/nJets",cjets.size());
            hs.fill("Zpeak/"+path_cat+"/nLooseJets",lJets.size());
            hs.fill("Zpeak/"+path_cat+"/nGoodVertices",*nGoodVertices);
            hs.fill("Zpeak/"+path_cat+"/nPV",*nPV);
            hs.fill("Zpeak/"+path_cat+"/PFminiIsoL1",p_l1_miniIso);
            hs.fill("Zpeak/"+path_cat+"/PFminiIsoL2",p_l2_miniIso);
            hs.fill("Zpeak/"+path_cat+"/ChargeL1",charge_l1);
            hs.fill("Zpeak/"+path_cat+"/ChargeL2",charge_l2);
            hs.fill("Zpeak/"+path_cat+"/ChargeProduct",charge_l1*charge_l2);
            if(!tauEvent){
               hs.fill("Zpeak_noTau/"+path_cat+"/met",met);
               hs.fill("Zpeak_noTau/"+path_cat+"/PuppiMet",MET_Puppi->p.Pt());
               hs.fill("Zpeak_noTau/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
               if(matchedLeptons) hs.fill("Zpeak_noTau_matchedLep/"+path_cat+"/met",met);
               if(tightRIso_mu && *is_mumu) hs.fill("Zpeak_noTau_tightRIsoMu/"+path_cat+"/met",met);
            }
         }
         
         
         if(rec_selection && mllRequ && trigger && !trigger_veto){
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/met",met);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/PuppiMet",MET_Puppi->p.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/CaloMet",MET_Calo->p.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/DeepMet",MET_Deep->p.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/RawMet",MET_Raw->p.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/NoHFMet",MET_NoHF->p.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/pTl1",p_l1.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/pTl2",p_l2.Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/etal1",(*is_ee)? (*electrons)[0].etaSC : p_l1.Eta());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/etal2",(*is_ee)? (*electrons)[1].etaSC : p_l2.Eta());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/phil1",p_l1.Phi());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/phil2",p_l2.Phi());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/ZpT",(p_l1+p_l2).Pt());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/mll",(p_l1+p_l2).M());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/dPhill",abs(p_l1.DeltaPhi(p_l2)));
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/nJets",cjets.size());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/nLooseJets",lJets.size());
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/nGoodVertices",*nGoodVertices);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/nPV",*nPV);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/PFminiIsoL1",p_l1_miniIso);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/PFminiIsoL2",p_l2_miniIso);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/ChargeL1",charge_l1);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/ChargeL2",charge_l2);
            hs.fill("Zpeak_noJetRequ/"+path_cat+"/ChargeProduct",charge_l1*charge_l2);
            if(lJets.size()>0) hs.fill("Zpeak_noJetRequ/"+path_cat+"/looseJetpT1",lJets[0].p.Pt());
         }
         
         if(rec_selection && mllRequ){
            hs.fill("Zpeak_noJetRequ_noTrigg/"+path_cat+"/pTl2",p_l2.Pt());
         }
         
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs2d.scaleLumi();
      hs.mergeOverflow();
      hs2d.mergeOverflow();
      file.Close();
      
   } // dataset loop
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   //Combine ttBar madGraph with high genMet sample
   // ~hs.combineSamples("DrellYanCOMB",{"DrellYan","DrellYan_HT"});
   // ~hs2d.combineSamples("DrellYanCOMB",{"DrellYan","DrellYan_HT"});
   // ~samplesToCombine.push_back("DrellYanCOMB");
   
   //Plotting part 1D
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("distributions_Philipp%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline/emu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline/mumu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline_noTau/ee/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline_noTau/emu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline_noTau/mumu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"baseline_noTau_tightRIsoMu/mumu/",{"met"}},
      {"baseline_noTau_matchedLep/ee/",{"met"}},
      {"baseline_noTau_matchedLep/emu/",{"met"}},
      {"baseline_noTau_matchedLep/mumu/",{"met"}},
      {"Zpeak/ee/",{"met","PuppiMet","dphi_metNearLep","genHT","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak/emu/",{"met","PuppiMet","dphi_metNearLep","genHT","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak/mumu/",{"met","PuppiMet","dphi_metNearLep","genHT","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ/ee/",{"met","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","PuppiMet","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ/emu/",{"met","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","PuppiMet","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ/mumu/",{"met","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","PuppiMet","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ_noTrigg/ee/",{"pTl2"}},
      {"Zpeak_noJetRequ_noTrigg/emu/",{"pTl2"}},
      {"Zpeak_noJetRequ_noTrigg/mumu/",{"pTl2"}},
      {"Zpeak_noTau/ee/",{"met","PuppiMet","dphi_metNearLep"}},
      {"Zpeak_noTau/emu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"Zpeak_noTau/mumu/",{"met","PuppiMet","dphi_metNearLep"}},
      {"Zpeak_noTau_tightRIsoMu/mumu/",{"met"}},
      {"Zpeak_noTau_matchedLep/ee/",{"met"}},
      {"Zpeak_noTau_matchedLep/emu/",{"met"}},
      {"Zpeak_noTau_matchedLep/mumu/",{"met"}},
      {"baseline_noTau_trigg_lastBin/ee/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_lastBin/emu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_lastBin/mumu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_secondlastBin/ee/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_secondlastBin/emu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_secondlastBin/mumu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstBin/ee/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstBin/emu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstBin/mumu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstFirstBin/ee/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstFirstBin/emu/",{"dphi_metNearLep"}},
      {"baseline_noTau_trigg_firstFirstBin/mumu/",{"dphi_metNearLep"}},
      };
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions_Philipp%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs,samplesToCombine);
   saveRatios(msPresel_vVars,saver_hist,hs,samplesToCombine);
   
   //Plotting part 2D
   std::map<TString,std::vector<TString>> msPresel_vVars2D={
      {"baseline/ee/",{"2d_MetVSdPhiMetNearLep"}},
      {"baseline/emu/",{"2d_MetVSdPhiMetNearLep"}},
      {"baseline/mumu/",{"2d_MetVSdPhiMetNearLep"}},
      {"baseline_noTau_trigg/ee/",{"2d_MetVSdPhiMetNearLep"}},
      {"baseline_noTau_trigg/emu/",{"2d_MetVSdPhiMetNearLep"}},
      {"baseline_noTau_trigg/mumu/",{"2d_MetVSdPhiMetNearLep"}},
      };
   
   //Save 2d histograms
   saveHistograms2D(msPresel_vVars2D,saver_hist,hs2d,samplesToCombine);
}
