//Script to analyze minTrees (not finished since runtime is way to long for being kind of a duplicate of distributions.cpp)

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
      
   for(TString channel:{"/ee","/mumu","/emu"}){
         hs.addHist(channel+"/Jet1_pt"   ,";p_{T}^{Jet1} (GeV);EventsBIN"           ,50,0,500);
         
         hs.addHist(channel+"/PuppiMET*cos(PuppiMET_phi)"   ,";Puppi p_{x}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(channel+"/PuppiMET*sin(PuppiMET_phi)"   ,";Puppi p_{y}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(channel+"/METunc_Puppi"                 ,";Puppi #sigma(p_{t}^{miss}) (GeV);EventsBIN"           ,50,0,100);
         hs.addHist(channel+"/MET*cos(PFMET_phi)"           ,";PF p_{x}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(channel+"/MET*sind(PFMET_phi)"          ,";PF p_{y}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(channel+"/HT*cos(HT_phi)"               ,";H_{x} (GeV);EventsBIN"                                ,50,-250,250);
         hs.addHist(channel+"/HT*sin(HT_phi)"               ,";H_{y} (GeV);EventsBIN"                                ,50,-250,250);
         hs.addHist(channel+"/nJets"                        ,";N_{jets} ;EventsBIN"                                  ,10,1.5,12.5);
         hs.addHist(channel+"/n_Interactions"               ,";N_{interactions} ;EventsBIN"                          ,50,0,100);
         hs.addHist(channel+"/Lep1_flavor"                  ,";Flavor_{l1} ;EventsBIN"                               ,2,0.5,2.5);
         hs.addHist(channel+"/Lep2_flavor"                  ,";Flavor_{l2} ;EventsBIN"                               ,2,0.5,2.5);
         hs.addHist(channel+"/Lep1_pt*cos(Lep1_phi)"        ,";p_{x}^{lep1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Lep1_pt*sin(Lep1_phi)"        ,";p_{y}^{lep1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Lep1_eta"                     ,";#eta^{lep1} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(channel+"/Lep1_E"                       ,";E^{lep1} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(channel+"/Lep2_pt*cos(Lep2_phi)"        ,";p_{x}^{lep2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Lep2_pt*sin(Lep2_phi)"        ,";p_{y}^{lep2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Lep2_eta"                     ,";#eta^{lep2} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(channel+"/Lep2_E"                       ,";E^{lep2} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(channel+"/Jet1_pt*cos(Jet1_phi)"        ,";p_{x}^{jet1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Jet1_pt*sin(Jet1_phi)"        ,";p_{y}^{jet1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Jet1_eta"                     ,";#eta^{jet1} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(channel+"/Jet1_E"                       ,";E^{jet1} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(channel+"/Jet2_pt*cos(Jet2_phi)"        ,";p_{x}^{jet2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Jet2_pt*sin(Jet2_phi)"        ,";p_{y}^{jet2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(channel+"/Jet2_eta"                     ,";#eta^{jet2} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(channel+"/Jet2_E"                       ,";E^{jet2} (GeV) ;EventsBIN"                            ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETnearJet"               ,";|#delta#phi(PF MET, near. jet)| ;EventsBIN"           ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETfarJet"                ,";|#delta#phi(PF MET, far. jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETleadJet"               ,";|#delta#phi(PF MET, lead jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETlead2Jet"              ,";|#delta#phi(PF MET, sublead jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETbJet"                  ,";|#delta#phi(PF MET, lead b jet)| ;EventsBIN"          ,32,0,3.2);
         hs.addHist(channel+"/dPhiLep1Lep2"                 ,";|#delta#phi(lep1, lep2)| ;EventsBIN"                  ,32,0,3.2);
         hs.addHist(channel+"/dPhiJet1Jet2"                 ,";|#delta#phi(jet1, jet2)| ;EventsBIN"                  ,32,0,3.2);
         hs.addHist(channel+"/METsig"                       ,";PF MET sig. ;EventsBIN"                               ,20,0,200);
         hs.addHist(channel+"/MHT"                          ,";M_{HT} (GeV);EventsBIN"                               ,100,0,2000);
         hs.addHist(channel+"/MT"                           ,";M_{T} (GeV);EventsBIN"                                ,50,0,1000);
         hs.addHist(channel+"/looseLeptonVeto"              ,";looseLeptonVeto ;EventsBIN"                           ,2,-0.5,1.5);
         hs.addHist(channel+"/dPhiMETnearJet_Puppi"         ,";|#delta#phi(Puppi MET, near. jet)| ;EventsBIN"        ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETfarJet_Puppi"          ,";|#delta#phi(Puppi MET, far. jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETleadJet_Puppi"         ,";|#delta#phi(Puppi MET, lead jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETlead2Jet_Puppi"        ,";|#delta#phi(Puppi MET, sublead jet)| ;EventsBIN"      ,32,0,3.2);
         hs.addHist(channel+"/dPhiMETbJet_Puppi"            ,";|#delta#phi(Puppi MET, lead b jet)| ;EventsBIN"       ,32,0,3.2);
         hs.addHist(channel+"/dPhiLep1bJet"                 ,";|#delta#phi(lep1, lead b jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(channel+"/dPhiLep1Jet1"                 ,";|#delta#phi(lep1, ljet1)| ;EventsBIN"                 ,32,0,3.2);
         hs.addHist(channel+"/mLL"                          ,";m_{ll} (GeV) ;EventsBIN"                              ,32,0,3.2);
         hs.addHist(channel+"/CaloMET*cos(CaloMET)"         ,";Calo p_{x}^{miss} (GeV);EventsBIN"                    ,50,-250,250);
         hs.addHist(channel+"/CaloMET*sind(CaloMET)"        ,";Calo p_{y}^{miss} (GeV);EventsBIN"                    ,50,-250,250);
         hs.addHist(channel+"/MT2"                          ,";MT_{2} (GeV);EventsBIN"                               ,50,0,200);
         hs.addHist(channel+"/vecsum_pT_allJet"             ,";vecsum_pT_allJet (GeV);EventsBIN"                     ,50,0,500);
         hs.addHist(channel+"/vecsum_pT_l1l2_allJet"        ,";vecsum_pT_l1l2_allJet (GeV);EventsBIN"                ,50,0,500);
         hs.addHist(channel+"/mass_l1l2_allJet"             ,";mass_l1l2_allJet (GeV);EventsBIN"                     ,50,0,500);
         hs.addHist(channel+"/ratio_vecsumpTlep_vecsumpTjet",";ratio_vecsumpTlep_vecsumpTjet ;EventsBIN"             ,50,0,500);
         hs.addHist(channel+"/mjj"                          ,";m_{jj} (GeV);EventsBIN"                               ,100,0,2000);
      }
   
   return hs;
}

class systHists_minTree
{
   public:
      systHists_minTree(TString const &systematicName, std::vector<TString> const &datasets):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(defineHistograms(datasets)),
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
   
   syst->hists_.fillweight(cat+"/Jet1_pt"                      ,reader.Jet1_pt,weight);
   
   syst->hists_.fillweight(cat+"/PuppiMET*cos(PuppiMET_phi)"   ,reader.PuppiMET*cos(reader.PuppiMET_phi),weight);
   syst->hists_.fillweight(cat+"/PuppiMET*sin(PuppiMET_phi)"   ,reader.PuppiMET*sin(reader.PuppiMET_phi),weight);
   syst->hists_.fillweight(cat+"/METunc_Puppi"                 ,reader.METunc_Puppi,weight);
   syst->hists_.fillweight(cat+"/MET*cos(PFMET_phi)"           ,reader.MET*cos(reader.PFMET_phi),weight);
   syst->hists_.fillweight(cat+"/MET*sind(PFMET_phi)"          ,reader.MET*sin(reader.PFMET_phi),weight);
   syst->hists_.fillweight(cat+"/HT*cos(HT_phi)"               ,reader.HT*cos(reader.HT_phi),weight);
   syst->hists_.fillweight(cat+"/HT*sin(HT_phi)"               ,reader.HT*sin(reader.HT_phi),weight);
   syst->hists_.fillweight(cat+"/nJets"                        ,reader.nJets,weight);
   syst->hists_.fillweight(cat+"/n_Interactions"               ,reader.n_Interactions,weight);
   syst->hists_.fillweight(cat+"/Lep1_flavor"                  ,reader.Lep1_flavor,weight);
   syst->hists_.fillweight(cat+"/Lep2_flavor"                  ,reader.Lep2_flavor,weight);
   syst->hists_.fillweight(cat+"/Lep1_pt*cos(Lep1_phi)"        ,reader.Lep1_pt*cos(reader.Lep1_phi),weight);
   syst->hists_.fillweight(cat+"/Lep1_pt*sin(Lep1_phi)"        ,reader.Lep1_pt*sin(reader.Lep1_phi),weight);
   syst->hists_.fillweight(cat+"/Lep1_eta"                     ,reader.Lep1_eta,weight);
   syst->hists_.fillweight(cat+"/Lep1_E"                       ,reader.Lep1_E,weight);
   syst->hists_.fillweight(cat+"/Lep2_pt*cos(Lep2_phi)"        ,reader.Lep2_pt*cos(reader.Lep2_phi),weight);
   syst->hists_.fillweight(cat+"/Lep2_pt*sin(Lep2_phi)"        ,reader.Lep2_pt*sin(reader.Lep2_phi),weight);
   syst->hists_.fillweight(cat+"/Lep2_eta"                     ,reader.Lep2_eta,weight);
   syst->hists_.fillweight(cat+"/Lep2_E"                       ,reader.Lep2_E,weight);
   syst->hists_.fillweight(cat+"/Jet1_pt*cos(Jet1_phi)"        ,reader.Jet1_pt*cos(reader.Jet1_phi),weight);
   syst->hists_.fillweight(cat+"/Jet1_pt*sin(Jet1_phi)"        ,reader.Jet1_pt*sin(reader.Jet1_phi),weight);
   syst->hists_.fillweight(cat+"/Jet1_eta"                     ,reader.Jet1_eta,weight);
   syst->hists_.fillweight(cat+"/Jet1_E"                       ,reader.Jet1_E,weight);
   syst->hists_.fillweight(cat+"/Jet2_pt*cos(Jet2_phi)"        ,reader.Jet2_pt*cos(reader.Jet2_phi),weight);
   syst->hists_.fillweight(cat+"/Jet2_pt*sin(Jet2_phi)"        ,reader.Jet2_pt*sin(reader.Jet2_phi),weight);
   syst->hists_.fillweight(cat+"/Jet2_eta"                     ,reader.Jet2_eta,weight);
   syst->hists_.fillweight(cat+"/Jet2_E"                       ,reader.Jet2_E,weight);
   syst->hists_.fillweight(cat+"/dPhiMETnearJet"               ,reader.dPhiMETnearJet,weight);
   syst->hists_.fillweight(cat+"/dPhiMETfarJet"                ,reader.dPhiMETfarJet,weight);
   syst->hists_.fillweight(cat+"/dPhiMETleadJet"               ,reader.dPhiMETleadJet,weight);
   syst->hists_.fillweight(cat+"/dPhiMETlead2Jet"              ,reader.dPhiMETlead2Jet,weight);
   syst->hists_.fillweight(cat+"/dPhiMETbJet"                  ,reader.dPhiMETbJet,weight);
   syst->hists_.fillweight(cat+"/dPhiLep1Lep2"                 ,reader.dPhiLep1Lep2,weight);
   syst->hists_.fillweight(cat+"/dPhiJet1Jet2"                 ,reader.dPhiJet1Jet2,weight);
   syst->hists_.fillweight(cat+"/METsig"                       ,reader.METsig,weight);
   syst->hists_.fillweight(cat+"/MHT"                          ,reader.MHT,weight);
   syst->hists_.fillweight(cat+"/MT"                           ,reader.MT,weight);
   syst->hists_.fillweight(cat+"/looseLeptonVeto"              ,reader.looseLeptonVeto,weight);
   syst->hists_.fillweight(cat+"/dPhiMETnearJet_Puppi"         ,reader.dPhiMETnearJet_Puppi,weight);
   syst->hists_.fillweight(cat+"/dPhiMETfarJet_Puppi"          ,reader.dPhiMETfarJet_Puppi,weight);
   syst->hists_.fillweight(cat+"/dPhiMETleadJet_Puppi"         ,reader.dPhiMETleadJet_Puppi,weight);
   syst->hists_.fillweight(cat+"/dPhiMETlead2Jet_Puppi"        ,reader.dPhiMETlead2Jet_Puppi,weight);
   syst->hists_.fillweight(cat+"/dPhiMETbJet_Puppi"            ,reader.dPhiMETbJet_Puppi,weight);
   syst->hists_.fillweight(cat+"/dPhiLep1bJet"                 ,reader.dPhiLep1bJet,weight);
   syst->hists_.fillweight(cat+"/dPhiLep1Jet1"                 ,reader.dPhiLep1Jet1,weight);
   syst->hists_.fillweight(cat+"/mLL"                          ,reader.mLL,weight);
   syst->hists_.fillweight(cat+"/CaloMET*cos(CaloMET)"         ,reader.CaloMET*cos(reader.CaloMET),weight);
   syst->hists_.fillweight(cat+"/CaloMET*sind(CaloMET)"        ,reader.CaloMET*sin(reader.CaloMET),weight);
   syst->hists_.fillweight(cat+"/MT2"                          ,reader.MT2,weight);
   syst->hists_.fillweight(cat+"/vecsum_pT_allJet"             ,reader.vecsum_pT_allJet,weight);
   syst->hists_.fillweight(cat+"/vecsum_pT_l1l2_allJet"        ,reader.vecsum_pT_l1l2_allJet,weight);
   syst->hists_.fillweight(cat+"/mass_l1l2_allJet"             ,reader.mass_l1l2_allJet,weight);
   syst->hists_.fillweight(cat+"/ratio_vecsumpTlep_vecsumpTjet",reader.ratio_vecsumpTlep_vecsumpTjet,weight);
   syst->hists_.fillweight(cat+"/mjj"                          ,reader.mjj,weight);
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
         if (syst->systematic_.type() == Systematic::nominal) {   //weightType systematics
             for (systHists_minTree* syst_weightType : systHists_weightType_vec){
                syst_weightType->hists_.setCurrentSample(dsName);
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
            
            if (syst->systematic_.type() == Systematic::nominal) {    //weightType systematics
               for (systHists_minTree* syst_weightType : systHists_weightType_vec){
                  if(syst_weightType->systematic_.variation() == Systematic::up){
                     fillHistograms(syst_weightType,reader,reader.systTotalWeights[syst_weightType->systematic_.type()].first);
                  }
                  else{
                     fillHistograms(syst_weightType,reader,reader.systTotalWeights[syst_weightType->systematic_.type()].second);
                  }
               }
            }
            
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
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
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
   
   // ~std::vector<TString> systToPlot = {"Nominal"};
   std::vector<TString> systToPlot = {"Nominal","PU_UP","PU_DOWN"};
   
   std::vector<systHists_minTree*> systHists_vec = buildSystHistVector(systToPlot,mcSamples,dataSamples,ttbarSamples,signalSamples);
   
   analyzeMinTree(systHists_vec);
   
   for(systHists_minTree* syst : systHists_vec){
      TH1F* temp = syst->hists_.getHistogram("/ee/Jet1_pt","TTbar_diLepton");
      
      std::cout<<syst->systematicName_<<"   "<<temp->GetBinContent(20)<<std::endl;
   }
}
