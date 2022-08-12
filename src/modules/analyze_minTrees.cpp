//Script to analyze minTrees (currently used to get simple distributions from minTrees)

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/systematics.hpp"
#include "tools/util.hpp"
#include "tools/minTreeReader.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <iomanip>
#include <chrono>
#include <algorithm>
#include <TTree.h>

Config const &cfg=Config::get();

hist::Histograms<TH1F> defineHistograms(std::vector<TString> const &datasets){
   hist::Histograms<TH1F> hs(datasets);
      
   for(TString channel:{"/ee","/mumu","/emu","/all"}){
      hs.addHist(channel+"/Jet1_pt"   ,";p_{T}^{Jet1} (GeV);EventsBIN"           ,50,0,500);
   }
   
   return hs;
}

hist::Histograms<TH2F> defineHistograms2D(std::vector<TString> const &datasets){
   hist::Histograms<TH2F> hs2D(datasets);
      
   for(TString channel:{"/ee","/mumu","/emu","/all"}){
      hs2D.addHist(channel+"/GenMetDiffMETRel_dPhiMETLep"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
      hs2D.addHist(channel+"/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
      hs2D.addHist(channel+"/genMet_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];genmet",20,0,3.14,6000,0,1000);
      hs2D.addHist(channel+"/leadTopPT_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];p_{T}(lead. top)",20,0,3.14,6000,0,1000);
      
      hs2D.addHist("baseline_met120"+channel+"/GenMetDiffMETRel_dPhiMETLep"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
      hs2D.addHist("baseline_met120"+channel+"/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
      hs2D.addHist("baseline_met120"+channel+"/genMet_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];genmet",20,0,3.14,6000,0,1000);
      hs2D.addHist("baseline_met120"+channel+"/leadTopPT_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];p_{T}(lead. top)",20,0,3.14,6000,0,1000);
      
      hs2D.addHist("baseline_genmet120"+channel+"/GenMetDiffMETRel_dPhiMETLep"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
      hs2D.addHist("baseline_genmet120"+channel+"/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";min[#Delta#phi(p_{T}^{miss},l)];|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   }
   
   return hs2D;
}

class systHists_minTree
{
   public:
      systHists_minTree(TString const &systematicName, std::vector<TString> const &datasets):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(defineHistograms(datasets)),
      hists2D_(defineHistograms2D(datasets)),
      systematicName_(systematicName),
      datasets_(datasets)
      {
         if (std::find(Systematic::weightTypes.begin(), Systematic::weightTypes.end(), systematic_.type()) != Systematic::weightTypes.end()){
            weightType = true;
         }
         else {
            weightType = false;
         }
      }
      
      Systematic::Systematic systematic_;
      hist::Histograms<TH1F> hists_;
      hist::Histograms<TH2F> hists2D_;
      std::vector<TString> datasets_;
      TString const systematicName_;
      bool weightType;
};

std::vector<systHists_minTree*> buildSystHistVector(std::vector<TString> const &systematics, std::vector<TString> const &mcSamples, std::vector<TString> const &dataSamples, 
                                                      std::vector<TString> const &ttbarSamples, std::vector<TString> const &signalSamples) {
   std::vector<TString> datasets_temp;
   std::vector<systHists_minTree*> systHists_vec;
   for (TString syst : systematics){
      
      if (syst == "Nominal") datasets_temp = util::addVectors(mcSamples,dataSamples);
      else if (syst=="CR_ENVELOPE_UP" || syst=="CR_ENVELOPE_DOWN" || syst=="MTOP_DOWN" || syst=="MTOP_UP") datasets_temp = signalSamples;
      else datasets_temp = mcSamples;
      
      systHists_minTree* temp = new systHists_minTree(syst,datasets_temp);
      systHists_vec.push_back(temp);
   }
   return systHists_vec;
}


void ProgressBar(const uint &progress){
   std::string progressBar = "[";
   for(uint i = 0; i <= progress; ++i){
      if(i%1 == 0) progressBar += "#";
   }
   for(uint i = 0; i < 100 - progress; ++i){
      if(i%1 == 0) progressBar += " ";
   }
   progressBar = progressBar + "] " + std::to_string(progress) + "% of Events processed";
   std::cout << "\r" << progressBar << std::flush;
   if(progress == 100) std::cout << std::endl;
}

void fillHistograms(systHists_minTree* syst, minTreeReader const &reader, float const &weight){
   TString cat = "";
   if (reader.emu == 1)  cat = "/emu";
   else if (reader.mumu == 1)  cat = "/mumu";
   else cat = "/ee";
   
   for (const TString& catString : {cat,TString("/all")}){
      syst->hists_.fillweight(catString+"/Jet1_pt"                      ,reader.Jet1_pt,weight);
      
      syst->hists2D_.fillweight(catString+"/GenMetDiffMETRel_dPhiMETLep"                      ,reader.PhiRec_xy,(reader.genMet-reader.MET_xy)/reader.genMet,weight);
      syst->hists2D_.fillweight(catString+"/GenMetDiffMETRel_dPhiMETLep_Puppi"                ,reader.PhiRecPuppi_xy,(reader.genMet-reader.PuppiMET_xy)/reader.genMet,weight);
      
      syst->hists2D_.fillweight(catString+"/genMet_dPhiMETLep_Puppi"                ,reader.PhiRecPuppi_xy,reader.genMet,weight);
      syst->hists2D_.fillweight(catString+"/leadTopPT_dPhiMETLep_Puppi"                ,reader.PhiRecPuppi_xy,reader.leadTop_pT,weight);
            
      if (reader.MET_xy>120){
         syst->hists2D_.fillweight("baseline_met120"+catString+"/GenMetDiffMETRel_dPhiMETLep"        ,reader.PhiRec_xy,(reader.genMet-reader.MET_xy)/reader.genMet,weight);
      }
      if (reader.PuppiMET_xy>120){
         syst->hists2D_.fillweight("baseline_met120"+catString+"/GenMetDiffMETRel_dPhiMETLep_Puppi"  ,reader.PhiRecPuppi_xy,(reader.genMet-reader.PuppiMET_xy)/reader.genMet,weight);
         syst->hists2D_.fillweight("baseline_met120"+catString+"/genMet_dPhiMETLep_Puppi"                ,reader.PhiRecPuppi_xy,reader.genMet,weight);
         syst->hists2D_.fillweight("baseline_met120"+catString+"/leadTopPT_dPhiMETLep_Puppi"                ,reader.PhiRecPuppi_xy,reader.leadTop_pT,weight);
      }
      
      if (reader.genMet>120) {
         syst->hists2D_.fillweight("baseline_genmet120"+catString+"/GenMetDiffMETRel_dPhiMETLep"        ,reader.PhiRec_xy,(reader.genMet-reader.MET_xy)/reader.genMet,weight);
         syst->hists2D_.fillweight("baseline_genmet120"+catString+"/GenMetDiffMETRel_dPhiMETLep_Puppi"  ,reader.PhiRecPuppi_xy,(reader.genMet-reader.PuppiMET_xy)/reader.genMet,weight);
      }
      
   }
}

void analyzeMinTree(std::vector<systHists_minTree*> &systHists_vec){
   
   // Store weightType systematic separately
   std::vector<systHists_minTree*> systHists_weightType_vec;
   for(systHists_minTree* &syst : systHists_vec){
      if (syst->weightType) systHists_weightType_vec.push_back(syst);
   }
   
   for(systHists_minTree* syst : systHists_vec){
      if (syst->weightType) continue;     //weightType systematics have no minTrees, but only weights in nominal minTree
      for (TString const &dsName : syst->datasets_){
      // ~for (TString const &dsName : {"TTbar_diLepton"}){
         std::cout<<"Processing "<<dsName<<std::endl;
         TString fileName = TString::Format("%s/../minTrees/100.0/%s%s_merged.root",cfg.outputDirectory.Data(),(syst->systematicName_+"/").Data(),dsName.Data());
         TFile* minTreeFile = TFile::Open(fileName, "READ");
         TTree* eventTree = (TTree*)minTreeFile->Get(TString::Format("ttbar_res100.0/%s",dsName.Data()));
         minTreeReader reader(*eventTree,syst->systematic_);
         
         syst->hists_.setCurrentSample(dsName);
         syst->hists2D_.setCurrentSample(dsName);
         if (syst->systematic_.type() == Systematic::nominal) {   //weightType systematics
             for (systHists_minTree* syst_weightType : systHists_weightType_vec){
                syst_weightType->hists_.setCurrentSample(dsName);
                syst_weightType->hists2D_.setCurrentSample(dsName);
             }
         }
         
         ProgressBar(0);
         int processed = 0;
         int processEvents=cfg.processFraction*eventTree->GetEntries();
         //begin event loop
         for(Int_t i=0; i<processEvents; ++i) {
            eventTree->GetEntry(i);
            
            ++processed;
            if(processed % 1000 == 0) {
                uint progress = 100*(float)processed/processEvents;
                ProgressBar(progress);
            }
            
            if (reader.PuppiMET<0) continue;    // take only reco events into account
            
            fillHistograms(syst,reader,reader.totalWeight);
            
            // only used for systematics
            // ~if (syst->systematic_.type() == Systematic::nominal) {    //weightType systematics
               // ~for (systHists_minTree* syst_weightType : systHists_weightType_vec){
                  // ~if(syst_weightType->systematic_.variation() == Systematic::up){
                     // ~fillHistograms(syst_weightType,reader,reader.systTotalWeights[syst_weightType->systematic_.type()].first);
                  // ~}
                  // ~else{
                     // ~fillHistograms(syst_weightType,reader,reader.systTotalWeights[syst_weightType->systematic_.type()].second);
                  // ~}
               // ~}
            // ~}
            
         }
         std::cout<<std::endl;
      }
   }
}

extern "C"
void run()
{
   std::vector<TString> mcSamples={};
   std::vector<TString> dataSamples={};
   std::vector<TString> ttbarSamples={};
   std::vector<TString> signalSamples={};
   switch(cfg.year_int){
      case(3): //2018
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      // ~ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      // ~signalSamples = {"TTbar_diLepton"};
      mcSamples = {"TTbar_diLepton"};
      dataSamples = {};
      ttbarSamples = {};
      signalSamples = {};
      break;
      case(2): //2017
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      signalSamples = {"TTbar_diLepton"};
      break;
      case(1): //2016
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      signalSamples = {"TTbar_diLepton"};
      dataSamples = {};
      break;
   }
   
   std::vector<TString> systToPlot = {"Nominal"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP","PU_DOWN"};
   
   std::vector<systHists_minTree*> systHists_vec = buildSystHistVector(systToPlot,mcSamples,dataSamples,ttbarSamples,signalSamples);
   
   analyzeMinTree(systHists_vec);
   
   for(systHists_minTree* syst : systHists_vec){      
      TString loc=TString::Format("minTree_analyze/histograms_%s.root",syst->systematicName_.Data());
      io::RootFileSaver saver_hist(loc,TString::Format("analyze_minTrees%.1f",cfg.processFraction*100),false);
      syst->hists_.saveHistograms(saver_hist,mcSamples);
      syst->hists2D_.saveHistograms(saver_hist,mcSamples);
   }
}
