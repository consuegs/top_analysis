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

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   float met_sf = cfg.sf.met_sf;
   TString met_sf_string="";
   if (met_sf!=1.0) met_sf_string="_SF"+std::to_string(met_sf);
   
   //Read TriggerSF hists
   io::RootFileReader triggerSF_ee(cfg.trigger_SF_ee,"");
   TH2F* triggerSF_ee_hist = (TH2F*)(triggerSF_ee.read<TH2F>("eff_histo"));
   io::RootFileReader triggerSF_mumu(cfg.trigger_SF_mumu,"");
   TH2F* triggerSF_mumu_hist = (TH2F*)(triggerSF_mumu.read<TH2F>("eff_histo"));
   io::RootFileReader triggerSF_emu(cfg.trigger_SF_emu,"");
   TH2F* triggerSF_emu_hist = (TH2F*)(triggerSF_emu.read<TH2F>("eff_histo"));
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   hist::Histograms<TH1F> hs_cutflow(vsDatasubsets);
   
   hs_cutflow.addHist("cutflow/ee"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   hs_cutflow.addHist("cutflow/emu"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   hs_cutflow.addHist("cutflow/mumu"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   
   for(TString selection:{"baseline","baseline_Met200"}){ //Reco 1D Histograms
      hs.addHist(selection+"/ee/met"   ,";%MET;EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/emu/met"   ,";%MET;EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/mumu/met"   ,";%MET;EventsBIN"           ,100,0,500);
      
      hs.addHist(selection+"/ee/met_puppi"   ,";%MET;EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/emu/met_puppi"   ,";%MET;EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/mumu/met_puppi"   ,";%MET;EventsBIN"           ,100,0,500);
      
      hs.addHist(selection+"/ee/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/mll"   ,";mll(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTsumlep"   ,";p_{T}^{ll};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTsumlep"   ,";p_{T}^{ll};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTsumlep"   ,";p_{T}^{ll};EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/sumpTlep"   ,";%pTl1+%pTl2;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/sumpTlep"   ,";%pTl1+%pTl2;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/sumpTlep"   ,";%pTl1+%pTl2;EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTbJet"   ,";p_{T}^{b};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTbJet"   ,";p_{T}^{b};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTbJet"   ,";p_{T}^{b};EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTJet1"   ,";p_{T}^{Jet1};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTJet1"   ,";p_{T}^{Jet1};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTJet1"   ,";p_{T}^{Jet1};EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTJet2"   ,";p_{T}^{Jet2};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTJet2"   ,";p_{T}^{Jet2};EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTJet2"   ,";p_{T}^{Jet2};EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_metLeadJet"   ,";|#Delta#phi|(p_{T}^{miss},leading jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metLeadJet"   ,";|#Delta#phi|(p_{T}^{miss},leading jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metLeadJet"   ,";|#Delta#phi|(p_{T}^{miss},leading jet);EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_metLead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd leading jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metLead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd leading jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metLead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd leading jet);EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
      hs.addHist(selection+"/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
     
      hs.addHist(selection+"/ee/dphi_metNearLep_puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
      hs.addHist(selection+"/emu/dphi_metNearLep_puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metNearLep_puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
      
      hs.addHist(selection+"/ee/COSdphi_metNearLep"   ,";cos(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,-1.,1);
      hs.addHist(selection+"/emu/COSdphi_metNearLep"   ,";cos(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,-1.,1);
      hs.addHist(selection+"/mumu/COSdphi_metNearLep"   ,";cos(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,-1.,1);
      
      hs.addHist(selection+"/ee/SINdphi_metNearLep"   ,";sin(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,0.,1);
      hs.addHist(selection+"/emu/SINdphi_metNearLep"   ,";sin(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,0.,1);
      hs.addHist(selection+"/mumu/SINdphi_metNearLep"   ,";sin(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,0.,1);
      
      hs.addHist(selection+"/ee/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},bjet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dR_bJetLep1"   ,";|#Delta R|(b Jet,l_{1});EventsBIN"           ,100,0,5);
      hs.addHist(selection+"/emu/dR_bJetLep1"   ,";|#Delta R|(b Jet,l_{1});EventsBIN"           ,100,0,5);
      hs.addHist(selection+"/mumu/dR_bJetLep1"   ,";|#Delta R|(b Jet,l_{1});EventsBIN"           ,100,0,5);
      
      hs.addHist(selection+"/ee/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_bJetnearLep"   ,";|#Delta#phi|(b Jet,next l);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_bJetnearLep"   ,";|#Delta#phi|(b Jet,next l);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_bJetnearLep"   ,";|#Delta#phi|(b Jet,next l);EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_b1b2"   ,";|#Delta#phi|(b Jet1,b Jet2);EventsBIN"           ,100,0,4);
      hs.addHist(selection+"/emu/dphi_b1b2"   ,";|#Delta#phi|(b Jet1,b Jet2);EventsBIN"           ,100,0,4);
      hs.addHist(selection+"/mumu/dphi_b1b2"   ,";|#Delta#phi|(b Jet1,b Jet2);EventsBIN"           ,100,0,4);
      
      hs.addHist(selection+"/ee/dR_b1b2"   ,";|#Delta R|(b Jet1,b Jet2);EventsBIN"           ,100,0,6);
      hs.addHist(selection+"/emu/dR_b1b2"   ,";|#Delta R|(b Jet1,b Jet2);EventsBIN"           ,100,0,6);
      hs.addHist(selection+"/mumu/dR_b1b2"   ,";|#Delta R|(b Jet1,b Jet2);EventsBIN"           ,100,0,6);
      
      hs.addHist(selection+"/ee/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_metLepsum"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}+l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metLepsum"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}+l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metLepsum"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}+l_{2});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2});EventsBIN"           ,100,0,3.2);
      
      hs.addHist(selection+"/ee/dR_Lep1Lep2"   ,";|#Delta R|(l_{1},l_{2});EventsBIN"           ,100,0,5);
      hs.addHist(selection+"/emu/dR_Lep1Lep2"   ,";|#Delta R|(l_{1},l_{2});EventsBIN"           ,100,0,5);
      hs.addHist(selection+"/mumu/dR_Lep1Lep2"   ,";|#Delta R|(l_{1},l_{2});EventsBIN"           ,100,0,5);
      
      hs.addHist(selection+"/ee/nJets"   ,";N_{Jets};EventsBIN"           ,11,-0.5,10.5);
      hs.addHist(selection+"/emu/nJets"   ,";N_{Jets};EventsBIN"           ,11,-0.5,10.5);
      hs.addHist(selection+"/mumu/nJets"   ,";N_{Jets};EventsBIN"           ,11,-0.5,10.5);
      
      hs.addHist(selection+"/ee/nBjets"   ,";N_{bJets};EventsBIN"           ,5,-0.5,4.5);
      hs.addHist(selection+"/emu/nBjets"   ,";N_{bJets};EventsBIN"           ,5,-0.5,4.5);
      hs.addHist(selection+"/mumu/nBjets"   ,";N_{bJets};EventsBIN"           ,5,-0.5,4.5);
      
      hs.addHist(selection+"/ee/mt2"   ,";MT2 (GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"/emu/mt2"   ,";MT2 (GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"/mumu/mt2"   ,";MT2 (GeV);EventsBIN"           ,100,0,600);
      
      hs.addHist(selection+"/ee/mt_MetLep1"   ,";M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/mt_MetLep1"   ,";M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/mt_MetLep1"   ,";M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/mt_MetLep2"   ,";M_{T}(p_{T}^{miss},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/mt_MetLep2"   ,";M_{T}(p_{T}^{miss},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/mt_MetLep2"   ,";M_{T}(p_{T}^{miss},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/mt_MetNextLep"   ,";M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/mt_MetNextLep"   ,";M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/mt_MetNextLep"   ,";M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/conMt_Lep1Lep2"   ,";conM_{T}(l_{1},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/conMt_Lep1Lep2"   ,";conM_{T}(l_{1},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/conMt_Lep1Lep2"   ,";conM_{T}(l_{1},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/ST"   ,";S_{T} (GeV);EventsBIN"           ,100,0,1500);
      hs.addHist(selection+"/emu/ST"   ,";S_{T} (GeV);EventsBIN"           ,100,0,1500);
      hs.addHist(selection+"/mumu/ST"   ,";S_{T} (GeV);EventsBIN"           ,100,0,1500);
      
      hs.addHist(selection+"/ee/HT"   ,";H_{T} (GeV);EventsBIN"           ,100,0,2500);
      hs.addHist(selection+"/emu/HT"   ,";H_{T} (GeV);EventsBIN"           ,100,0,2500);
      hs.addHist(selection+"/mumu/HT"   ,";H_{T} (GeV);EventsBIN"           ,100,0,2500);
      
      hs.addHist(selection+"/ee/sum_STHT"   ,";S_{T}+H_{T} (GeV);EventsBIN"           ,100,0,4000);
      hs.addHist(selection+"/emu/sum_STHT"   ,";S_{T}+H_{T} (GeV);EventsBIN"           ,100,0,4000);
      hs.addHist(selection+"/mumu/sum_STHT"   ,";S_{T}+H_{T} (GeV);EventsBIN"           ,100,0,4000);
      
      hs.addHist(selection+"/ee/sum_mlb"   ,";sum m_{lb} (GeV);EventsBIN"           ,100,0,3000);
      hs.addHist(selection+"/emu/sum_mlb"   ,";sum m_{lb} (GeV);EventsBIN"           ,100,0,3000);
      hs.addHist(selection+"/mumu/sum_mlb"   ,";sum m_{lb} (GeV);EventsBIN"           ,100,0,3000);
   }//End 1D reco histograms
   
   
   for(TString selection:{"genParticles"}){ //1D Gen Histograms
      hs.addHist(selection+"/ee/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/dphi_NeutrinoLep"   ,";|#Delta#phi|_{gen}(#nu,l);EventsBIN"           ,100,0,4);
      hs.addHist(selection+"/emu/dphi_NeutrinoLep"   ,";|#Delta#phi|_{gen}(#nu,l);EventsBIN"           ,100,0,4);
      hs.addHist(selection+"/mumu/dphi_NeutrinoLep"   ,";|#Delta#phi|_{gen}(#nu,l);EventsBIN"           ,100,0,4);
      
      hs.addHist(selection+"/ee/dR_NeutrinoLep"   ,";|#Delta R|_{gen}(#nu,l);EventsBIN"           ,100,0,6);
      hs.addHist(selection+"/emu/dR_NeutrinoLep"   ,";|#Delta R|_{gen}(#nu,l);EventsBIN"           ,100,0,6);
      hs.addHist(selection+"/mumu/dR_NeutrinoLep"   ,";|#Delta R|_{gen}(#nu,l);EventsBIN"           ,100,0,6);
      
      hs.addHist(selection+"/ee/pTtop1"   ,";p_{T}^{gen}(t_{1});EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTtop1"   ,";p_{T}^{gen}(t_{1});EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTtop1"   ,";%p_{T}^{gen}(t_{1});EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/pTtop2"   ,";p_{T}^{gen}(t_{2});EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pTtop2"   ,";p_{T}^{gen}(t_{2});EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pTtop2"   ,";%p_{T}^{gen}(t_{2});EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/genHT"   ,";genHT;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/genHT"   ,";genHT;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/genHT"   ,";genHT;EventsBIN"           ,100,0,1000);
   }//End 1D Gen Histograms
   
   //Define 2D histograms in the following
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   
   hs2d.addHist("genParticles/ee/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/emu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/mumu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   //Ntuple and file to save minimal ttbar tree used for binning studies
   float minTree_MET, minTree_PtNuNu, minTree_PhiRec, minTree_PhiGen, minTree_PhiNuNu, minTree_PhiMetNearJet, minTree_PhiMetFarJet, minTree_PhiMetLeadJet, minTree_PhiMetLead2Jet,
   minTree_PhiMetbJet, minTree_dPhiLep1Lep2, minTree_dPhiJet1Jet2, minTree_METsig, minTree_N, minTree_SF, minTree_genMet, minTree_PuppiMet, minTree_XYcorrMet, 
   minTree_HT, minTree_MHT, minTree_MT, minTree_genMT, minTree_MT_nextLep, minTree_genMT_nextLep,
   minTree_PhiPtnunuMet, minTree_leadTop, minTree_dPhiNuNu, minTree_PhiRecPuppi, minTree_PhiRecXYcorr, minTree_PhiMetNearJet_Puppi, minTree_PhiMetFarJet_Puppi,
   minTree_PhiMetLeadJet_Puppi, minTree_PhiMetLead2Jet_Puppi, minTree_PhiMetbJet_Puppi, minTree_dPhiLep1bJet, minTree_dPhiLep1Jet1, minTree_ratioMET_sqrtMETunc_Puppi,
   minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet, minTree_METunc_Puppi, minTree_METunc_PF, minTree_absmetres_PUPPI,
   minTree_Lep1_pt, minTree_Lep1_phi, minTree_Lep1_eta, minTree_Lep1_E, minTree_Lep1_flavor,
   minTree_Lep2_pt, minTree_Lep2_phi, minTree_Lep2_eta, minTree_Lep2_E, minTree_Lep2_flavor,
   minTree_Jet1_pt, minTree_Jet1_phi, minTree_Jet1_eta, minTree_Jet1_E, minTree_Jet1_bTagScore, minTree_Jet1_unc,
   minTree_Jet2_pt, minTree_Jet2_phi, minTree_Jet2_eta, minTree_Jet2_E, minTree_Jet2_bTagScore, minTree_Jet2_unc,
   minTree_PFMET_phi, minTree_PuppiMET_phi, minTree_CaloMET, minTree_CaloMET_phi, minTree_nJets, minTree_n_Interactions;
   UInt_t minTree_runNo, minTree_lumNo, minTree_genDecayMode, minTree_n_Interactions_gen, minTree_looseLeptonVeto, minTree_NpromptNeutrinos, minTree_NnonpromptNeutrinos;
   ULong64_t minTree_evtNo;
   // ~std::vector<float> minTree_v_bJet_muonFraction;
   // ~std::vector<float> minTree_v_bJet_electronFraction;
   // ~std::vector<float> minTree_v_Jet_muonFraction;
   // ~std::vector<float> minTree_v_Jet_electronFraction;
   // ~std::vector<float> minTree_v_bJet_muonEnergy;
   // ~std::vector<float> minTree_v_bJet_electronEnergy;
   // ~std::vector<float> minTree_v_Jet_muonEnergy;
   // ~std::vector<float> minTree_v_Jet_electronEnergy;
   io::RootFileSaver ttbar_res_saver(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res"+met_sf_string+"%.1f_new.root",cfg.processFraction*100),TString::Format("ttbar_res%.1f",cfg.processFraction*100),true,false);
   TTree ttbar_res("ttbar_res","ttbar_res");
   ttbar_res.Branch("MET",&minTree_MET,"MET/f");
   ttbar_res.Branch("PtNuNu",&minTree_PtNuNu,"PtNuNu/f");
   ttbar_res.Branch("Phi_rec",&minTree_PhiRec,"Phi_rec/f");
   ttbar_res.Branch("Phi_gen",&minTree_PhiGen,"Phi_gen/f");
   ttbar_res.Branch("Phi_NuNu",&minTree_PhiNuNu,"Phi_NuNu/f");
   ttbar_res.Branch("dPhiMETnearJet",&minTree_PhiMetNearJet,"dPhiMETnearJet/f");
   ttbar_res.Branch("dPhiMETfarJet",&minTree_PhiMetFarJet,"dPhiMETfarJet/f");
   ttbar_res.Branch("dPhiMETleadJet",&minTree_PhiMetLeadJet,"dPhiMETleadJet/f");
   ttbar_res.Branch("dPhiMETlead2Jet",&minTree_PhiMetLead2Jet,"dPhiMETlead2Jet/f");
   ttbar_res.Branch("dPhiMETbJet",&minTree_PhiMetbJet,"dPhiMETbJet/f");
   ttbar_res.Branch("dPhiLep1Lep2",&minTree_dPhiLep1Lep2,"dPhiLep1Lep2/f");
   ttbar_res.Branch("dPhiJet1Jet2",&minTree_dPhiJet1Jet2,"dPhiJet1Jet2/f");
   ttbar_res.Branch("METsig",&minTree_METsig,"METsig/f");
   ttbar_res.Branch("N",&minTree_N,"N/f");
   ttbar_res.Branch("SF",&minTree_SF,"SF/f");
   ttbar_res.Branch("runNo",&minTree_runNo,"runNo/i");
   ttbar_res.Branch("lumNo",&minTree_lumNo,"lumNo/i");
   ttbar_res.Branch("evtNo",&minTree_evtNo,"evtNo/l");
   ttbar_res.Branch("genDecayMode",&minTree_genDecayMode,"genDecayMode/i");
   ttbar_res.Branch("genMET",&minTree_genMet,"genMET/f");
   ttbar_res.Branch("PuppiMET",&minTree_PuppiMet,"PuppiMET/f");
   ttbar_res.Branch("XYcorrMET",&minTree_XYcorrMet,"XYcorrMET/f");
   ttbar_res.Branch("HT",&minTree_HT,"HT/f");
   ttbar_res.Branch("MHT",&minTree_MHT,"MHT/f");
   ttbar_res.Branch("MT",&minTree_MT,"MT/f");
   ttbar_res.Branch("genMT",&minTree_genMT,"genMT/f");
   ttbar_res.Branch("MT_nextLep",&minTree_MT_nextLep,"MT_nextLep/f");
   ttbar_res.Branch("genMT_nextLep",&minTree_genMT_nextLep,"genMT_nextLep/f");
   ttbar_res.Branch("n_Interactions",&minTree_n_Interactions,"n_Interactions/f");
   ttbar_res.Branch("n_Interactions_gen",&minTree_n_Interactions_gen,"n_Interactions_gen/i");
   ttbar_res.Branch("dPhiPtnunuMet",&minTree_PhiPtnunuMet,"dPhiPtnunuMet/f");
   ttbar_res.Branch("leadTop_pT",&minTree_leadTop,"leadTop_pT/f");
   ttbar_res.Branch("dPhiNuNu",&minTree_dPhiNuNu,"dPhiNuNu/f");
   ttbar_res.Branch("Phi_recPuppi",&minTree_PhiRecPuppi,"Phi_recPuppi/f");
   ttbar_res.Branch("Phi_recXYcorr",&minTree_PhiRecXYcorr,"Phi_recXYcorr/f");
   ttbar_res.Branch("looseLeptonVeto",&minTree_looseLeptonVeto,"looseLeptonVeto/i");
   ttbar_res.Branch("nJets",&minTree_nJets,"nJets/f");
   ttbar_res.Branch("dPhiMETnearJet_Puppi",&minTree_PhiMetNearJet_Puppi,"dPhiMETnearJet_Puppi/f");
   ttbar_res.Branch("dPhiMETfarJet_Puppi",&minTree_PhiMetFarJet_Puppi,"dPhiMETfarJet_Puppi/f");
   ttbar_res.Branch("dPhiMETleadJet_Puppi",&minTree_PhiMetLeadJet_Puppi,"dPhiMETleadJet_Puppi/f");
   ttbar_res.Branch("dPhiMETlead2Jet_Puppi",&minTree_PhiMetLead2Jet_Puppi,"dPhiMETlead2Jet_Puppi/f");
   ttbar_res.Branch("dPhiMETbJet_Puppi",&minTree_PhiMetbJet_Puppi,"dPhiMETbJet_Puppi/f");
   ttbar_res.Branch("dPhiLep1bJet",&minTree_dPhiLep1bJet,"dPhiLep1bJet/f");
   ttbar_res.Branch("dPhiLep1Jet1",&minTree_dPhiLep1Jet1,"dPhiLep1Jet1/f");
   ttbar_res.Branch("ratioMET_sqrtMETunc_Puppi",&minTree_ratioMET_sqrtMETunc_Puppi,"ratioMET_sqrtMETunc_Puppi/f");
   ttbar_res.Branch("ratio_pTj1_vecsum_pT_l1_l2_bjet",&minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet,"ratio_pTj1_vecsum_pT_l1_l2_bjet/f");
   ttbar_res.Branch("METunc_Puppi",&minTree_METunc_Puppi,"METunc_Puppi/f");
   ttbar_res.Branch("METunc_PF",&minTree_METunc_PF,"METunc_PF/f");
   ttbar_res.Branch("absmetres_PUPPI",&minTree_absmetres_PUPPI,"absmetres_PUPPI/f");
   ttbar_res.Branch("Lep1_pt",&minTree_Lep1_pt,"Lep1_pt/f");
   ttbar_res.Branch("Lep1_phi",&minTree_Lep1_phi,"Lep1_phi/f");
   ttbar_res.Branch("Lep1_eta",&minTree_Lep1_eta,"Lep1_eta/f");
   ttbar_res.Branch("Lep1_E",&minTree_Lep1_E,"Lep1_E/f");
   ttbar_res.Branch("Lep1_flavor",&minTree_Lep1_flavor,"Lep1_flavor/f");
   ttbar_res.Branch("Lep2_pt",&minTree_Lep2_pt,"Lep2_pt/f");
   ttbar_res.Branch("Lep2_phi",&minTree_Lep2_phi,"Lep2_phi/f");
   ttbar_res.Branch("Lep2_eta",&minTree_Lep2_eta,"Lep2_eta/f");
   ttbar_res.Branch("Lep2_E",&minTree_Lep2_E,"Lep2_E/f");
   ttbar_res.Branch("Lep2_flavor",&minTree_Lep2_flavor,"Lep2_flavor/f");
   ttbar_res.Branch("Jet1_pt",&minTree_Jet1_pt,"Jet1_pt/f");
   ttbar_res.Branch("Jet1_phi",&minTree_Jet1_phi,"Jet1_phi/f");
   ttbar_res.Branch("Jet1_eta",&minTree_Jet1_eta,"Jet1_eta/f");
   ttbar_res.Branch("Jet1_E",&minTree_Jet1_E,"Jet1_E/f");
   ttbar_res.Branch("Jet1_bTagScore",&minTree_Jet1_bTagScore,"Jet1_bTagScore/f");
   ttbar_res.Branch("Jet1_unc",&minTree_Jet1_unc,"Jet1_unc/f");
   ttbar_res.Branch("Jet2_pt",&minTree_Jet2_pt,"Jet2_pt/f");
   ttbar_res.Branch("Jet2_phi",&minTree_Jet2_phi,"Jet2_phi/f");
   ttbar_res.Branch("Jet2_eta",&minTree_Jet2_eta,"Jet2_eta/f");
   ttbar_res.Branch("Jet2_E",&minTree_Jet2_E,"Jet2_E/f");
   ttbar_res.Branch("Jet2_bTagScore",&minTree_Jet2_bTagScore,"Jet2_bTagScore/f");
   ttbar_res.Branch("Jet2_unc",&minTree_Jet2_unc,"Jet2_unc/f");
   ttbar_res.Branch("PFMET_phi",&minTree_PFMET_phi,"PFMET_phi/f");
   ttbar_res.Branch("PuppiMET_phi",&minTree_PuppiMET_phi,"PuppiMET_phi/f");
   ttbar_res.Branch("CaloMET",&minTree_CaloMET,"CaloMET/f");
   ttbar_res.Branch("CaloMET_phi",&minTree_CaloMET_phi,"CaloMET_phi/f");
   ttbar_res.Branch("NpromptNeutrinos",&minTree_NpromptNeutrinos,"NpromptNeutrinos/i");
   ttbar_res.Branch("NnonpromptNeutrinos",&minTree_NnonpromptNeutrinos,"NnonpromptNeutrinos/i");
   // ~ttbar_res.Branch("bJet_muonFraction",&minTree_v_bJet_muonFraction);
   // ~ttbar_res.Branch("bJet_electronFraction",&minTree_v_bJet_electronFraction);
   // ~ttbar_res.Branch("Jet_muonFraction",&minTree_v_Jet_muonFraction);
   // ~ttbar_res.Branch("Jet_electronFraction",&minTree_v_Jet_electronFraction);
   // ~ttbar_res.Branch("bJet_muonEnergy",&minTree_v_bJet_muonEnergy);
   // ~ttbar_res.Branch("bJet_electronEnergy",&minTree_v_bJet_electronEnergy);
   // ~ttbar_res.Branch("Jet_muonEnergy",&minTree_v_Jet_muonEnergy);
   // ~ttbar_res.Branch("Jet_electronEnergy",&minTree_v_Jet_electronEnergy);
   
   
   //Additional map to calculate signal efficiencies
   std::map<TString,float> count;
   std::map<TString,float> Ngen;
   
   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");

      bool const isData=dss.isData;
      bool const isSignal=dss.isSignal;
      
      hs.setCurrentSample(dss.name);
      hs_cutflow.setCurrentSample(dss.name);
      hs2d.setCurrentSample(dss.name);
      
      //Lumi weight for current sample
      float lumi_weight=dss.xsec/float(dss.Ngen)*cfg.lumi;
      
      //Save number of gen events for efficiency
      Ngen[dss.datasetName]=dss.Ngen;
      
      //Check if current sample is TTbar powheg (standard ttbar MC)
      bool ttBar_standard=false;
      if (dss.datasetName=="TTbar") ttBar_standard=true;
      
      //Check if current sample is TTbar powheg dilepton
      bool ttBar_dilepton=false;
      if (dss.datasetName=="TTbar_diLepton") ttBar_dilepton=true;
      
      //Check if current sample is TTbar powheg dilepton
      bool ttBar_dilepton_CP5=false;
      if (dss.datasetName=="TTbar_diLepton_CP5") ttBar_dilepton_CP5=true;
      
      //Check if current sample is TTbar madGraph (has extra genMet>150 part)
      bool ttBar_madGraph=false;
      if (dss.datasetName=="TTbar_madGraph") ttBar_madGraph=true;
      
      //Check if current sample is TTbar madGraph high MET
      bool ttBar_madGraph150=false;
      if (dss.datasetName=="TTbar_madGraph150") ttBar_madGraph150=true;
      
      //Check if current sample is TTbar madGraph high MET
      bool ttBar_amcatnlo=false;
      if (dss.datasetName=="TTbar_amcatnlo") ttBar_amcatnlo=true;
      
      //Check if current sample is selectedSusy scenario 
      bool SUSY_T2tt_650_350=false;
      if (dss.datasetName=="T2tt_650_350") SUSY_T2tt_650_350=true;
      
      //Check if current sample is selectedDM scenario 
      bool DM_scalar_1_200=false;
      if (dss.datasetName=="DM_scalar_1_200") DM_scalar_1_200=true;
      
      //Check if current sample is SUSY scenario
      bool SUSY_scen=false;
      if (dss.name.find("SMS")!=std::string::npos) SUSY_scen=true;
      
      //Check which PD is currently used
      bool SingleElectron=false;
      bool SingleMuon=false;
      bool DoubleEG=false;
      bool DoubleMuon=false;
      bool MuonEG=false;
      if (dss.datasetName.find("SingleElectron")!=std::string::npos) SingleElectron=true;
      else if (dss.datasetName.find("SingleMuon")!=std::string::npos) SingleMuon=true;
      else if (dss.datasetName.find("DoubleEG")!=std::string::npos) DoubleEG=true;
      else if (dss.datasetName.find("DoubleMuon")!=std::string::npos) DoubleMuon=true;
      else if (dss.datasetName.find("MuonEG")!=std::string::npos) MuonEG=true;
      
      //Check if current sample is Run2016H
      bool Run2016H=false;
      if (dss.datasetName.find("Run2016H")!=std::string::npos) Run2016H=true;

      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<UInt_t> runNo(reader, "runNo");
      TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
      TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
      TTreeReaderValue<float> w_topPT(reader, "topPTweight");
      TTreeReaderValue<float> w_bTag(reader, "bTagWeight");
      TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
      TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<tree::MET> MET_XYcorr(reader, "metXYcorr");
      TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
      TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
      TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
      TTreeReaderValue<int> n_Interactions(reader, "nPV");
      TTreeReaderValue<int> n_Interactions_gen(reader, "true_nPV");
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
   
   
      int iEv=0;
      int processEvents=cfg.processFraction*dss.entries;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }
         
         float met=MET->p.Pt();
         float const met_puppi=MET_Puppi->p.Pt();
         float const genMet=GENMET->p.Pt();
         
         //Booleans for reco and pseudo selection
         // ~bool rec_selection=false;
         bool rec_selection=false;
         bool pseudo_selection=true;
         
         //For ttBar madGraph sample, only use genMet<150 due to extension
         if (ttBar_madGraph && genMet>150) continue;
         
         //Do not use tau events in signal sample
         if (ttBar_dilepton && *genDecayMode>3) continue;
         
         //Trigger selection
         bool diElectronTriggers=*eleTrigg || *singleEleTrigg;
         bool diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2;
         bool electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
         bool triggerData=false;
         bool triggerData_veto=true;
         bool triggerMC=true;
         
         if (!isData){
            if(!(*is_ee && diElectronTriggers) && !(*is_mumu && diMuonTriggers) && !(*is_emu && electronMuonTriggers)) triggerMC=false;
         }
         else{
            //ee
            if(SingleElectron && *is_ee){
               triggerData=*singleEleTrigg;
               triggerData_veto=*eleTrigg;
            }
            else if(DoubleEG && *is_ee){
               triggerData=*eleTrigg;
               triggerData_veto=false;
            }
            //mumu
            else if(SingleMuon && *is_mumu){
               triggerData=*singleMuonTrigg1 || *singleMuonTrigg2;
               if (Run2016H) triggerData_veto=*muonTrigg3 || *muonTrigg4;
               else triggerData_veto=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4;
            }
            else if(DoubleMuon && *is_mumu){
               if (Run2016H) triggerData=*muonTrigg3 || *muonTrigg4;
               else triggerData=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4;
               triggerData_veto=false;
            }
            //emu
            else if(SingleMuon && *is_emu){
               triggerData=*singleMuonTrigg1 || *singleMuonTrigg2;
               if (Run2016H) triggerData_veto=*eleMuTrigg3 || *eleMuTrigg4;
               else triggerData_veto=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4;
            }
            else if(SingleElectron && *is_emu){
               triggerData=*singleEleTrigg;
               if (Run2016H) triggerData_veto=*eleMuTrigg3 || *eleMuTrigg4;
               else triggerData_veto=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4;
            }
            else if(MuonEG && *is_emu){
               if (Run2016H) triggerData=*eleMuTrigg3 || *eleMuTrigg4;
               else triggerData=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4;
               triggerData_veto=false;
            }
            else continue;
            
            if(!(triggerData && !triggerData_veto)) continue;
         }
         
         
         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter)
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         int flavor_l1=0;  //1 for electron and 2 for muon
         int flavor_l2=0;
         
         int leadLepton=0;
         int subleadLepton=1;
         bool muonLead=true; //Boolean for emu channel
         TString cat="";
         
         
         if (*is_ee){
            rec_selection=true;
            if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) rec_selection=false; //currently double check since trees only have tight leptons!!
            if(abs((*electrons)[0].etaSC)>2.4 || abs((*electrons)[1].etaSC>2.4)) rec_selection=false; //To use same region as for muons, cut on supercluster eta
            if((*electrons)[0].p.Pt()*(*electrons)[0].corr<(*electrons)[1].p.Pt()*(*electrons)[1].corr) leadLepton=1,subleadLepton=0;
            p_l1=(*electrons)[leadLepton].p*(*electrons)[leadLepton].corr;
            p_l2=(*electrons)[subleadLepton].p*(*electrons)[subleadLepton].corr;
            p_l1=(*electrons)[leadLepton].p*(*electrons)[leadLepton].corr;
            p_l2=(*electrons)[subleadLepton].p*(*electrons)[subleadLepton].corr;
            flavor_l1=1;
            flavor_l2=1;
            cat="ee";
         }
         else if (*is_mumu){
            rec_selection=true;
            if(!(*muons)[0].isTight || !(*muons)[1].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 || (*muons)[1].rIso>0.15) rec_selection=false;
            if(abs((*muons)[0].p.Eta())>2.4 || abs((*muons)[1].p.Eta())>2.4) rec_selection=false;
            if((*muons)[0].p.Pt()*(*muons)[0].rochesterCorrection<(*muons)[1].p.Pt()*(*muons)[1].rochesterCorrection) leadLepton=1,subleadLepton=0;
            p_l1=(*muons)[leadLepton].p*(*muons)[leadLepton].rochesterCorrection;
            p_l2=(*muons)[1].p*(*muons)[subleadLepton].rochesterCorrection;
            flavor_l1=2;
            flavor_l2=2;
            cat="mumu";
         }
         else if (*is_emu){
            rec_selection=true;
            if(!(*muons)[0].isTight || !(*electrons)[0].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 ) rec_selection=false;
            if(abs((*muons)[0].p.Eta())>2.4) rec_selection=false;
            if(abs((*electrons)[0].etaSC>2.4) ) rec_selection=false;
            if ((*muons)[0].p.Pt()*(*muons)[0].rochesterCorrection>(*electrons)[0].p.Pt()*(*electrons)[0].corr){
               p_l1=(*muons)[0].p*(*muons)[0].rochesterCorrection;
               p_l2=(*electrons)[0].p*(*electrons)[0].corr;
               flavor_l1=2;
               flavor_l2=1;
            }
            else {
               p_l1=(*electrons)[0].p*(*electrons)[0].corr;
               p_l2=(*muons)[0].p*(*muons)[0].rochesterCorrection;
               flavor_l1=1;
               flavor_l2=2;
               muonLead=false;
            }
            cat="emu";
         }
         
         if (triggerMC==false) rec_selection=false;
         
         if(triggerMC){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,1,1);
            else hs_cutflow.fillweight("cutflow/"+cat,1,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
         }
         
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,2,1);
            else hs_cutflow.fillweight("cutflow/"+cat,2,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
         }
         
         if (p_l1.Pt()<25 || p_l2.Pt()<20) rec_selection=false; //eta cuts already done in TreeWriter
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,3,1);
            else hs_cutflow.fillweight("cutflow/"+cat,3,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
         }
         
         // ~if (*mll<20 || ((*is_ee || *is_mumu) && *mll<106 && *mll>76)) rec_selection=false;
         float mll_corr=(p_l1+p_l2).M();
         if (mll_corr<20 || ((*is_ee || *is_mumu) && mll_corr<106 && mll_corr>76)) rec_selection=false;
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,4,1);
            else hs_cutflow.fillweight("cutflow/"+cat,4,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
         }
         
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         if (cjets.size()<2) rec_selection=false;
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,5,1);
            else hs_cutflow.fillweight("cutflow/"+cat,5,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
         }
         
         // ~if ((*is_ee || *is_mumu) && met<40) rec_selection=false;
         if ((*is_ee || *is_mumu) && met_puppi<40) rec_selection=false;
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,6,1);
            else hs_cutflow.fillweight("cutflow/"+cat,6,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_topPT);
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
         if (!bTag) rec_selection=false;
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,7,1);
            else hs_cutflow.fillweight("cutflow/"+cat,7,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_bTag * *w_topPT);
         }
                           
         // end reco baseline selection
         
         if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
         
         if(rec_selection==false && pseudo_selection==false) continue;  // only proceed with events selected by one of the baseline selection (reco or pseudo)
         
         //Set Event weights
         float fEventWeight = 1.;
         float SFWeight = 1.;
         float triggerSF = 1.;
         float p_l1_trigg = std::min(199.,p_l1.Pt());
         float p_l2_trigg = std::min(199.,p_l2.Pt());
         if(!isData){
            if (*is_ee) triggerSF = triggerSF_ee_hist->GetBinContent(triggerSF_ee_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_ee_hist->GetYaxis()->FindBin(p_l2_trigg));
            else if (*is_mumu) triggerSF = triggerSF_mumu_hist->GetBinContent(triggerSF_mumu_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_mumu_hist->GetYaxis()->FindBin(p_l2_trigg));
            else if (*is_emu){
               if(muonLead)triggerSF = triggerSF_emu_hist->GetBinContent(triggerSF_emu_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_emu_hist->GetYaxis()->FindBin(p_l2_trigg));
               else triggerSF = triggerSF_emu_hist->GetBinContent(triggerSF_emu_hist->GetXaxis()->FindBin(p_l2_trigg), triggerSF_emu_hist->GetYaxis()->FindBin(p_l1_trigg));
            }
            fEventWeight=*w_pu * *w_mc;     //Set event weight 
            SFWeight=*sf_lep1 * *sf_lep2 * *w_topPT * *w_bTag * triggerSF;     //Set combined SF weight
            hs.setFillWeight(fEventWeight*SFWeight);
            hs2d.setFillWeight(fEventWeight*SFWeight);
         }
         else{
            hs.setFillWeight(1);
            hs2d.setFillWeight(1);
            hs_cutflow.setFillWeight(1);
         }
         if(rec_selection){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,8,1);
            else hs_cutflow.fillweight("cutflow/"+cat,8,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_bTag * *w_topPT *triggerSF);
         }
         if(rec_selection && !*addLepton){
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,9,1);
            else hs_cutflow.fillweight("cutflow/"+cat,9,*w_pu * *w_mc * *sf_lep2 * *sf_lep1 * *w_bTag * *w_topPT *triggerSF);
         }
                  
         //Muon and Electron Fraction for bJets and cleaned jets
         // ~minTree_v_bJet_muonFraction.clear();
         // ~minTree_v_bJet_electronFraction.clear();
         // ~minTree_v_bJet_muonEnergy.clear();
         // ~minTree_v_bJet_electronEnergy.clear();
         // ~minTree_v_bJet_muonFraction.reserve(BJets.size());
         // ~minTree_v_bJet_electronFraction.reserve(BJets.size());
         // ~minTree_v_bJet_muonEnergy.reserve(BJets.size());
         // ~minTree_v_bJet_electronEnergy.reserve(BJets.size());
         // ~for (tree::Jet const &jet : BJets) {
            // ~minTree_v_bJet_muonFraction.push_back(jet.muonf);
            // ~minTree_v_bJet_electronFraction.push_back(jet.electronf);
            // ~minTree_v_bJet_muonEnergy.push_back(jet.muonf*jet.p.Pt());
            // ~minTree_v_bJet_electronEnergy.push_back(jet.electronf*jet.p.Pt());
         // ~}
         // ~minTree_v_Jet_muonFraction.clear();
         // ~minTree_v_Jet_electronFraction.clear();
         // ~minTree_v_Jet_muonEnergy.clear();
         // ~minTree_v_Jet_electronEnergy.clear();
         // ~minTree_v_Jet_muonFraction.reserve(cjets.size());
         // ~minTree_v_Jet_electronFraction.reserve(cjets.size());
         // ~minTree_v_Jet_muonEnergy.reserve(cjets.size());
         // ~minTree_v_Jet_electronEnergy.reserve(cjets.size());
         // ~for (tree::Jet const &jet : cjets) {
            // ~minTree_v_Jet_muonFraction.push_back(jet.muonf);
            // ~minTree_v_Jet_electronFraction.push_back(jet.electronf);
            // ~minTree_v_Jet_muonEnergy.push_back(jet.muonf*jet.p.Pt());
            // ~minTree_v_Jet_electronEnergy.push_back(jet.electronf*jet.p.Pt());
         // ~}
         
         //Scale measured met by met_sf from config
         met=met_sf*met;
         
         // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
         TLorentzVector neutrinoPair(0,0,0,0);
         neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
         
         //Calculate MET before parton shower and hadronization for DM scenario and SUSY scenarios
         //And get number of (non)prompt neutrinos in event
         int NpromptNeutrinos=0;
         int NnonpromptNeutrinos=0;
         for (auto const &genParticle : *genParticles){
            if((!SUSY_scen && abs(genParticle.pdgId)>100000) || (SUSY_scen && abs(genParticle.pdgId)==1000022)){
               neutrinoPair+=genParticle.p;
            }
            if(abs(genParticle.pdgId)==12 || abs(genParticle.pdgId)==14 || abs(genParticle.pdgId)==16){
               if(genParticle.isPrompt) NpromptNeutrinos++;
               else {
                  NnonpromptNeutrinos++;
               }
            }
         }
         
         //Get DeltaPhi between MET (or genMet or neutrino pT) and nearest Lepton
         float dPhiMETnearLep=4;
         float dPhiMETnearLepPuppi=4;
         float dPhiMETnearLepXYcorr=4;
         float mt_MetNextLep=0; 
         float mt_NuNuNextLep=0; 
         for (TLorentzVector const lep : {p_l1,p_l2}){
            const float dPhi=MET->p.DeltaPhi(lep);
            const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
            const float dPhi_XYcorr=MET_XYcorr->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
               dPhiMETnearLep=dPhi;
               mt_MetNextLep=phys::M_T(MET->p,lep);
            }
            if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLepPuppi)) {
               dPhiMETnearLepPuppi=dPhi_Puppi;
            }
            if (std::abs(dPhi_XYcorr) < std::abs(dPhiMETnearLepXYcorr)) {
               dPhiMETnearLepXYcorr=dPhi_XYcorr;
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
               mt_NuNuNextLep=phys::M_T(neutrinoPair,lep);
            }
         }
         
         //Get DeltaPhi between MET nearest/leading/farest (b) Jet and HT
         float dPhiMETnearJet=4;
         float dPhiMETfarJet=0;
         float dPhiMETleadJet=4;
         float dPhiMETlead2Jet=4;
         float dPhiMetBJet=4;
         float dPhiMETnearJet_Puppi=4;
         float dPhiMETfarJet_Puppi=0;
         float dPhiMETleadJet_Puppi=4;
         float dPhiMETlead2Jet_Puppi=4;
         float dPhiMetBJet_Puppi=4;
         float dPhiLep1Lep2=4;
         float dPhiJet1Jet2=4;
         float dPhiLep1bJet=4;
         float dPhiLep1Jet1=4;
         TLorentzVector MHT(0.,0.,0.,0.);
         float HT=0;
         float ratio_pTj1_vecsum_pT_l1_l2_bjet=0;
         if (rec_selection){ 
            for (tree::Jet const &jet : cjets) {
               const float dPhi=MET->p.DeltaPhi(jet.p);
               const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(jet.p);
               if (std::abs(dPhi) < std::abs(dPhiMETnearJet)) dPhiMETnearJet=dPhi;
               if (std::abs(dPhi) > std::abs(dPhiMETfarJet)) dPhiMETfarJet=dPhi;
               if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearJet_Puppi)) dPhiMETnearJet_Puppi=dPhi_Puppi;
               if (std::abs(dPhi_Puppi) > std::abs(dPhiMETfarJet_Puppi)) dPhiMETfarJet_Puppi=dPhi_Puppi;
               HT+=jet.p.Pt();
               MHT+=jet.p;
            }
            dPhiMETleadJet=MET->p.DeltaPhi(cjets[0].p);
            dPhiMETlead2Jet=MET->p.DeltaPhi(cjets[1].p);
            dPhiMetBJet=MET->p.DeltaPhi(BJets[0].p);
            dPhiMETleadJet_Puppi=MET_Puppi->p.DeltaPhi(cjets[0].p);
            dPhiMETlead2Jet_Puppi=MET_Puppi->p.DeltaPhi(cjets[1].p);
            dPhiMetBJet_Puppi=MET_Puppi->p.DeltaPhi(BJets[0].p);
            dPhiLep1Lep2=p_l1.DeltaPhi(p_l2);
            dPhiJet1Jet2=cjets[0].p.DeltaPhi(cjets[1].p);
            dPhiLep1bJet=p_l1.DeltaPhi(BJets[0].p);
            dPhiLep1Jet1=p_l1.DeltaPhi(cjets[0].p);
            ratio_pTj1_vecsum_pT_l1_l2_bjet=cjets[0].p.Pt()/(p_l1+p_l2+BJets[0].p).Pt();
         }
         
         //Further variables
         float mt_MetLep1=phys::M_T(MET->p,p_l1);
         float mt_MetLep2=phys::M_T(MET->p,p_l2);
         float sum_mlb=phys::sumMlb(p_l1,p_l2,cjets,BJets);
         float conMt_Lep1Lep2=phys::conM_T(p_l1,p_l2);
         float dPhiNuNu=genNeutrino->DeltaPhi(*genAntiNeutrino);
         
         TLorentzVector leadGenLepton;
         if (genLepton->Pt()>genAntiLepton->Pt()) leadGenLepton=*genLepton;
         else leadGenLepton=*genAntiLepton;
         
         float mt_genMetLep1=phys::M_T(GENMET->p,leadGenLepton);
         float mt_genNeutrinosLep1=phys::M_T(neutrinoPair,leadGenLepton);
         
         //Sort genTop pT
         std::vector<TLorentzVector> gen_tops;
         gen_tops.push_back(*genTop);
         gen_tops.push_back(*genAntiTop);
         sort(gen_tops.begin(), gen_tops.end(), tree::PtGreaterLorentz);
         float pT_top1=0;
         float pT_top2=0;
         pT_top1=gen_tops[0].Pt();
         pT_top2=gen_tops[1].Pt();
         
         //Fill minimal tree for TTbar resolution used in binning studies
         if (ttBar_dilepton || ttBar_dilepton_CP5 || ttBar_madGraph || ttBar_madGraph150 || SUSY_T2tt_650_350 || DM_scalar_1_200 || ttBar_standard || ttBar_amcatnlo){
            minTree_MET=met;
            minTree_PtNuNu=neutrinoPair.Pt();
            minTree_PhiRec=abs(dPhiMETnearLep);
            minTree_PhiGen=abs(dPhigenMETnearLep);
            minTree_PhiNuNu=abs(dPhiPtNunearLep);
            minTree_PhiMetNearJet=abs(dPhiMETnearJet);
            minTree_PhiMetFarJet=abs(dPhiMETfarJet);
            minTree_PhiMetLeadJet=abs(dPhiMETleadJet);
            minTree_PhiMetLead2Jet=abs(dPhiMETlead2Jet);
            minTree_PhiMetbJet=abs(dPhiMetBJet);
            minTree_dPhiLep1Lep2=abs(dPhiLep1Lep2);
            minTree_dPhiJet1Jet2=abs(dPhiJet1Jet2);
            minTree_METsig=MET->sig;
            minTree_N=lumi_weight*fEventWeight;
            minTree_SF=SFWeight;
            minTree_runNo=*runNo;
            minTree_lumNo=*lumNo;
            minTree_evtNo=*evtNo;
            minTree_genDecayMode=*genDecayMode;
            minTree_genMet=genMet;
            minTree_PuppiMet=met_puppi;
            minTree_XYcorrMet=MET_XYcorr->p.Pt();
            minTree_HT=HT;
            minTree_MHT=MHT.M();
            minTree_MT=mt_MetLep1;
            minTree_genMT=mt_genNeutrinosLep1;
            minTree_n_Interactions=*n_Interactions;
            minTree_n_Interactions_gen=*n_Interactions_gen;
            minTree_MT_nextLep=mt_MetNextLep;
            minTree_genMT_nextLep=mt_NuNuNextLep;
            minTree_leadTop=pT_top1;
            minTree_dPhiNuNu=abs(dPhiNuNu);
            minTree_PhiRecPuppi=abs(dPhiMETnearLepPuppi);
            minTree_PhiRecXYcorr=abs(dPhiMETnearLepXYcorr);
            minTree_looseLeptonVeto= *addLepton? 1: 0;
            minTree_nJets=cjets.size();
            minTree_PhiMetNearJet_Puppi=abs(dPhiMETnearJet_Puppi);
            minTree_PhiMetFarJet_Puppi=abs(dPhiMETfarJet_Puppi);
            minTree_PhiMetLeadJet_Puppi=abs(dPhiMETleadJet_Puppi);
            minTree_PhiMetLead2Jet_Puppi=abs(dPhiMETlead2Jet_Puppi);
            minTree_PhiMetbJet_Puppi=abs(dPhiMetBJet_Puppi);
            minTree_dPhiLep1bJet=abs(dPhiLep1bJet);
            minTree_dPhiLep1Jet1=abs(dPhiLep1Jet1);
            minTree_ratioMET_sqrtMETunc_Puppi=MET_Puppi->uncertainty/sqrt(met_puppi);
            minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet=ratio_pTj1_vecsum_pT_l1_l2_bjet;
            minTree_METunc_Puppi=MET_Puppi->uncertainty;
            minTree_METunc_PF=MET->uncertainty;
            minTree_absmetres_PUPPI=abs(met_puppi-genMet);
            minTree_PFMET_phi=MET->p.Phi();
            minTree_PuppiMET_phi=MET_Puppi->p.Phi();
            minTree_CaloMET=MET_Calo->p.Pt();
            minTree_CaloMET_phi=MET_Calo->p.Phi();
            minTree_NpromptNeutrinos=NpromptNeutrinos;
            minTree_NnonpromptNeutrinos=NnonpromptNeutrinos;
            if (rec_selection){
               minTree_Lep1_pt=p_l1.Pt();
               minTree_Lep1_phi=p_l1.Phi();
               minTree_Lep1_eta=p_l1.Eta();
               minTree_Lep1_E=p_l1.E();
               minTree_Lep1_flavor=flavor_l1;
               minTree_Lep2_pt=p_l2.Pt();
               minTree_Lep2_phi=p_l2.Phi();
               minTree_Lep2_eta=p_l2.Eta();
               minTree_Lep2_E=p_l2.E();
               minTree_Lep2_flavor=flavor_l2;
               minTree_Jet1_pt=cjets[0].p.Pt();
               minTree_Jet1_phi=cjets[0].p.Phi();
               minTree_Jet1_eta=cjets[0].p.Eta();
               minTree_Jet1_E=cjets[0].p.E();
               minTree_Jet1_bTagScore=cjets[0].bTagDeepCSV;
               minTree_Jet1_unc=cjets[0].uncert;
               minTree_Jet2_pt=cjets[1].p.Pt();
               minTree_Jet2_phi=cjets[1].p.Phi();
               minTree_Jet2_eta=cjets[1].p.Eta();
               minTree_Jet2_E=cjets[1].p.E();
               minTree_Jet2_bTagScore=cjets[1].bTagDeepCSV;
               minTree_Jet2_unc=cjets[1].uncert;
            }
            if (rec_selection==false) {
               minTree_MET=-1.;
               minTree_PhiRec=-1.;
               minTree_PhiMetNearJet=-1.;
               minTree_PhiMetFarJet=-1.;
               minTree_PhiMetLeadJet=-1.;
               minTree_PhiMetLead2Jet=-1.;
               minTree_PhiMetbJet=-1.;
               minTree_PhiMetbJet=-1.;
               minTree_METsig=-1.;
               minTree_PuppiMet=-1.;
               minTree_XYcorrMet=-1.;
               minTree_MHT=-1.;
               minTree_HT=-1.;
               minTree_MT=-1.;
               minTree_MT_nextLep=-1.;
               minTree_SF=0.;
               minTree_PhiPtnunuMet=-1;
               minTree_PhiRecPuppi=-1;
               minTree_PhiRecXYcorr=-1;
               minTree_n_Interactions=0;
               minTree_looseLeptonVeto=5;
               minTree_nJets=-1;
               minTree_ratioMET_sqrtMETunc_Puppi=-20;
               minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet=-1;
               minTree_METunc_PF=-1;
               minTree_METunc_Puppi=-1;
               minTree_absmetres_PUPPI=-1;
               minTree_Lep1_pt=-1;
               minTree_Lep1_phi=-5;
               minTree_Lep1_eta=-3;
               minTree_Lep1_E=-1;
               minTree_Lep1_flavor=-1;
               minTree_Lep2_pt=-1;
               minTree_Lep2_phi=-5;
               minTree_Lep2_eta=-3;
               minTree_Lep2_E=-1;
               minTree_Lep2_flavor=-1;
               minTree_Jet1_pt=-1;
               minTree_Jet1_phi=-5;
               minTree_Jet1_eta=-3;
               minTree_Jet1_E=-1;
               minTree_Jet1_bTagScore=-5;
               minTree_Jet1_unc=-1;
               minTree_Jet2_pt=-1;
               minTree_Jet2_phi=-5;
               minTree_Jet2_eta=-3;
               minTree_Jet2_E=-1;
               minTree_Jet2_unc=-1;
               minTree_CaloMET=-1;
               minTree_CaloMET_phi=-1;
               minTree_PuppiMET_phi=-1;
               minTree_PFMET_phi=-1;
            }
            else if (pseudo_selection==false) {
               minTree_PtNuNu=-1.;
               minTree_PhiGen=-1.;
               minTree_PhiNuNu=-1.;
               minTree_genMT=-1;
               minTree_genMT_nextLep=-1;
               minTree_PhiPtnunuMet=-1;
               minTree_leadTop=-1;
               minTree_dPhiNuNu=-1;
            }
            else {
               minTree_PhiPtnunuMet=abs(neutrinoPair.DeltaPhi(MET->p));
            }
            ttbar_res.Fill();
         }
         
         if(rec_selection==false) continue;  //fill the following histograms only with events selected by the reco baseline selection
         
         // Bjet and angular variables
         int nBjets=BJets.size();
         float dPhiLep1BJet=p_l1.DeltaPhi(BJets[0].p);
         float dRLep1BJet=p_l1.DeltaR(BJets[0].p);
         float dPhiLep2BJet=p_l2.DeltaPhi(BJets[0].p);
         float dphi_bJetnearLep=std::min(abs(dPhiLep1BJet),abs(dPhiLep2BJet));
         float dPhiLep1MET=p_l1.DeltaPhi(MET->p);
         float dPhiLep2MET=p_l2.DeltaPhi(MET->p);
         float dR_Lep1Lep2=p_l1.DeltaR(p_l2);
         float dPhiMetLepSum=MET->p.DeltaPhi(p_l1+p_l2);
         
         float dPhiBjets=4;
         float dRBjets=6;
         if (nBjets>1) {
            dPhiBjets=BJets[0].p.DeltaPhi(BJets[1].p);
            dRBjets=BJets[0].p.DeltaR(BJets[1].p);
         }
         
         //Get distance between neutrino and lepton from same W
         float dPhiNeutrinoLep1=4;
         float dPhiNeutrinoLep2=4;
         float dRNeutrinoLep1=6;
         float dRNeutrinoLep2=6;
         dPhiNeutrinoLep1=genNeutrino->DeltaPhi(*genAntiLepton);
         dPhiNeutrinoLep2=genAntiNeutrino->DeltaPhi(*genLepton);
         dRNeutrinoLep1=genNeutrino->DeltaR(*genAntiLepton);
         dRNeutrinoLep2=genAntiNeutrino->DeltaR(*genLepton);
         
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
         
         //Calculate ST
         float ST=met+p_l1.Pt()+p_l2.Pt();
         
         //Count individual events for efficiency
         count[dss.datasetName]++;
         
         //Fill hists
         TString path_cat="ee";
         if (*is_emu) path_cat="emu";
         else if (*is_mumu) path_cat="mumu";
         
         hs.fill("baseline/"+path_cat+"/met",met);
         hs.fill("baseline/"+path_cat+"/met_puppi",met_puppi);
         hs.fill("baseline/"+path_cat+"/met1000",met);
         hs.fill("baseline/"+path_cat+"/mll",*mll);
         hs.fill("baseline/"+path_cat+"/pTlep1",p_l1.Pt());
         hs.fill("baseline/"+path_cat+"/pTlep2",p_l2.Pt());
         hs.fill("baseline/"+path_cat+"/pTsumlep",(p_l1+p_l2).Pt());
         hs.fill("baseline/"+path_cat+"/sumpTlep",p_l1.Pt()+p_l2.Pt());
         hs.fill("baseline/"+path_cat+"/pTbJet",BJets[0].p.Pt());
         hs.fill("baseline/"+path_cat+"/pTJet1",cjets[0].p.Pt());
         hs.fill("baseline/"+path_cat+"/pTJet2",cjets[1].p.Pt());
         hs.fill("baseline/"+path_cat+"/dphi_metJet",abs(dPhiMETnearJet));
         hs.fill("baseline/"+path_cat+"/dphi_metLeadJet",abs(dPhiMETleadJet));
         hs.fill("baseline/"+path_cat+"/dphi_metLead2Jet",abs(dPhiMETlead2Jet));
         hs.fill("baseline/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
         hs.fill("baseline/"+path_cat+"/dphi_metNearLep_puppi",abs(dPhiMETnearLepPuppi));
         hs.fill("baseline/"+path_cat+"/COSdphi_metNearLep",TMath::Cos(abs(dPhiMETnearLep)));
         hs.fill("baseline/"+path_cat+"/SINdphi_metNearLep",TMath::Sin(abs(dPhiMETnearLep)));
         hs.fill("baseline/"+path_cat+"/dphi_metBJet",abs(dPhiMetBJet));
         hs.fill("baseline/"+path_cat+"/dphi_bJetLep1",abs(dPhiLep1BJet));
         hs.fill("baseline/"+path_cat+"/dR_bJetLep1",abs(dRLep1BJet));
         hs.fill("baseline/"+path_cat+"/dphi_bJetLep2",abs(dPhiLep2BJet));
         hs.fill("baseline/"+path_cat+"/dphi_bJetnearLep",dphi_bJetnearLep);
         hs.fill("baseline/"+path_cat+"/dphi_b1b2",abs(dPhiBjets));
         hs.fill("baseline/"+path_cat+"/dR_b1b2",abs(dRBjets));
         hs.fill("baseline/"+path_cat+"/dphi_metLep1",abs(dPhiLep1MET));
         hs.fill("baseline/"+path_cat+"/dphi_metLep2",abs(dPhiLep2MET));
         hs.fill("baseline/"+path_cat+"/dphi_metLepsum",abs(dPhiMetLepSum));
         hs.fill("baseline/"+path_cat+"/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
         hs.fill("baseline/"+path_cat+"/dR_Lep1Lep2",abs(dR_Lep1Lep2));
         hs.fill("baseline/"+path_cat+"/nJets",cjets.size());
         hs.fill("baseline/"+path_cat+"/nBjets",nBjets);
         hs.fill("baseline/"+path_cat+"/mt2",*mt2);
         hs.fill("baseline/"+path_cat+"/mt_MetLep1",mt_MetLep1);
         hs.fill("baseline/"+path_cat+"/mt_MetLep2",mt_MetLep2);
         hs.fill("baseline/"+path_cat+"/mt_MetNextLep",mt_MetNextLep);
         hs.fill("baseline/"+path_cat+"/conMt_Lep1Lep2",conMt_Lep1Lep2);
         hs.fill("baseline/"+path_cat+"/ST",ST);
         hs.fill("baseline/"+path_cat+"/HT",HT);
         hs.fill("baseline/"+path_cat+"/sum_STHT",ST+HT);
         hs.fill("baseline/"+path_cat+"/sum_mlb",sum_mlb);
         hs.fill("genParticles/"+path_cat+"/pT_nunu",neutrinoPair.Pt());
         hs.fill("genParticles/"+path_cat+"/genMet",genMet);
         hs.fill("genParticles/"+path_cat+"/DMgenMet",DMgenMET.Pt());
         hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep1));
         hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep2));
         hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep1));
         hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep2));
         hs.fill("genParticles/"+path_cat+"/pTtop1",pT_top1);
         hs.fill("genParticles/"+path_cat+"/pTtop2",pT_top2);
         hs.fill("genParticles/"+path_cat+"/genHT",*HTgen);
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep_Puppi",met_puppi,abs(dPhiMETnearLepPuppi));
         hs2d.fill("genParticles/"+path_cat+"/2d_PtNuNuVSdPhiNuNuNearLep",neutrinoPair.Pt(),abs(dPhiPtNunearLep));
         
         if (met>200){
            hs.fill("baseline_Met200/"+path_cat+"/met",met);
            hs.fill("baseline_Met200/"+path_cat+"/met_puppi",met_puppi);
            hs.fill("baseline_Met200/"+path_cat+"/met1000",met);
            hs.fill("baseline_Met200/"+path_cat+"/mll",*mll);
            hs.fill("baseline_Met200/"+path_cat+"/pTlep1",p_l1.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/pTlep2",p_l2.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/pTsumlep",(p_l1+p_l2).Pt());
            hs.fill("baseline_Met200/"+path_cat+"/sumpTlep",p_l1.Pt()+p_l2.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/pTbJet",BJets[0].p.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/pTJet1",cjets[0].p.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/pTJet2",cjets[1].p.Pt());
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metJet",abs(dPhiMETnearJet));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metNearLep_puppi",abs(dPhiMETnearLepPuppi));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metLeadJet",abs(dPhiMETleadJet));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metLead2Jet",abs(dPhiMETlead2Jet));
            hs.fill("baseline_Met200/"+path_cat+"/COSdphi_metNearLep",TMath::Cos(abs(dPhiMETnearLep)));
            hs.fill("baseline_Met200/"+path_cat+"/SINdphi_metNearLep",TMath::Sin(abs(dPhiMETnearLep)));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metBJet",abs(dPhiMetBJet));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_bJetLep1",abs(dPhiLep1BJet));
            hs.fill("baseline_Met200/"+path_cat+"/dR_bJetLep1",abs(dRLep1BJet));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_bJetLep2",abs(dPhiLep2BJet));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_bJetnearLep",dphi_bJetnearLep);
            hs.fill("baseline_Met200/"+path_cat+"/dphi_b1b2",abs(dPhiBjets));
            hs.fill("baseline_Met200/"+path_cat+"/dR_b1b2",abs(dRBjets));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metLep1",abs(dPhiLep1MET));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metLep2",abs(dPhiLep2MET));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_metLepsum",abs(dPhiMetLepSum));
            hs.fill("baseline_Met200/"+path_cat+"/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
            hs.fill("baseline_Met200/"+path_cat+"/dR_Lep1Lep2",abs(dR_Lep1Lep2));
            hs.fill("baseline_Met200/"+path_cat+"/nJets",cjets.size());
            hs.fill("baseline_Met200/"+path_cat+"/nBjets",nBjets);
            hs.fill("baseline_Met200/"+path_cat+"/mt2",*mt2);
            hs.fill("baseline_Met200/"+path_cat+"/mt_MetLep1",mt_MetLep1);
            hs.fill("baseline_Met200/"+path_cat+"/mt_MetLep2",mt_MetLep2);
            hs.fill("baseline_Met200/"+path_cat+"/mt_MetNextLep",mt_MetNextLep);
            hs.fill("baseline_Met200/"+path_cat+"/conMt_Lep1Lep2",conMt_Lep1Lep2);
            hs.fill("baseline_Met200/"+path_cat+"/ST",ST);
            hs.fill("baseline_Met200/"+path_cat+"/HT",HT);
            hs.fill("baseline_Met200/"+path_cat+"/sum_STHT",ST+HT);
            hs.fill("baseline_Met200/"+path_cat+"/sum_mlb",sum_mlb);
         }
               
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs_cutflow.scaleLumi();
      hs2d.scaleLumi();
      hs.mergeOverflow();
      hs_cutflow.mergeOverflow();
      hs2d.mergeOverflow();
      file.Close();
      
      //Save ntuple for TTbar resolution used in binning studies
      if (ttBar_standard) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res");
         ttbar_res.Reset();
      }
      else if (ttBar_madGraph150) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_MadGraph");
         ttbar_res.Reset();
      }
      else if (SUSY_T2tt_650_350) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_T2tt_650_350");
         ttbar_res.Reset();
      }
      else if (DM_scalar_1_200) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_DM_scalar_1_200");
         ttbar_res.Reset();
      }
      else if (ttBar_dilepton) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_dilepton");
         ttbar_res.Reset();
      }
      else if (ttBar_dilepton_CP5) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_dilepton_CP5");
         ttbar_res.Reset();
      }
      else if (ttBar_amcatnlo) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_amcatnlo");
         ttbar_res.Reset();
      }
      
   } // dataset loop
   // ~ttbar_res_saver.closeFile();
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   // ~std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar_diLepton","TTbar_madGraph","TTbar_madGraph150","TTbar_singleLepton",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM","ttH_SM"};
   // ~std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar_madGraph","TTbar_madGraph150",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM","ttH_SM"};
   // ~std::vector<TString> samplesToCombine={"TTbar","TTbar_diLepton","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM","ttH_SM"};
   // ~std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM","ttH_SM"};
   // ~std::vector<TString> samplesToCombine={"TTbar",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM","ttH_SM"};
   hs.combineFromSubsamples(samplesToCombine);
   hs_cutflow.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   //Combine ttBar madGraph with high genMet sample
   // ~hs.combineSamples("TTbar_madGraphCOMB",{"TTbar_madGraph","TTbar_madGraph150"});
   // ~hs2d.combineSamples("TTbar_madGraphCOMB",{"TTbar_madGraph","TTbar_madGraph150"});
   // ~samplesToCombine.push_back("TTbar_madGraphCOMB");
   
   //Plotting part 1D
   io::RootFileSaver saver(TString::Format("plots"+met_sf_string+"%.1f.root",cfg.processFraction*100),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/emu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metBJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/mumu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metBJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/ee/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/emu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/mumu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","dphi_metNearLep_puppi","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"genParticles/ee/",{"pT_nunu","genMet","DMgenMet","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2","genHT"}},
      {"genParticles/emu/",{"pT_nunu","genMet","DMgenMet","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2","genHT"}},
      {"genParticles/mumu/",{"pT_nunu","genMet","DMgenMet","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2","genHT"}},
      };
   
   // The following can be used for plotting, which is currently done with an extra module   
   // ~for (auto const &sPresel_vVars:msPresel_vVars){
      // ~TString const &sPresel=sPresel_vVars.first;
      // ~for (TString sVar:sPresel_vVars.second){
         // ~TString loc;
         // ~loc=sPresel+sVar;
         // ~TH1F* hist=hs.getHistogram(loc,"TTbar");
         // ~gfx::LegendEntries le=hs.getLegendEntries();
         // ~TString cat;
         // ~if (sPresel.Contains("ee/")) cat="ee";
         // ~else if (sPresel.Contains("emu/")) cat="e#mu";
         // ~else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         // ~TLatex label=gfx::cornerLabel(cat,2);
         // ~hist->SetStats(0);
         // ~hist->SetMarkerSize(0);
         // ~if (sVar=="diff_ptNuNu_genMET") can.SetLogy(0);
         // ~hist->Draw("histE");
         // ~TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         // ~leg.Draw();
         // ~label.Draw();
         // ~saver.save(can,"tt_only/"+loc);
         
         // ~THStack st_mc=hs.getStack(loc,{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar"});
         // ~le=hs.getLegendEntries();
         // ~st_mc.Draw();
         // ~TLegend leg2=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         // ~leg2.Draw();
         // ~label.Draw();
         // ~saver.save(can,"all/"+loc);
         
         // ~//Plot Stack also with Signal
         // ~le=hs.getLegendEntries();
         // ~st_mc.Draw();
         // ~auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM"});
         // ~for (auto const &h: hists) h->Draw("same hist");
         // ~le+=hs.getLegendEntries();
         // ~TLegend leg3=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         // ~leg3.Draw();
         // ~label.Draw();
         // ~saver.save(can,"all_withSignal/"+loc);
      // ~}
   // ~}
      
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms"+met_sf_string+"_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs,samplesToCombine);
   saveHistograms({{"cutflow/",{"ee","emu","mumu"}}},saver_hist,hs_cutflow,samplesToCombine);
   
   //Plotting part 2D
   std::map<TString,std::vector<TString>> msPresel_vVars2D={
      {"genParticles/ee/",{"2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"genParticles/emu/",{"2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"genParticles/mumu/",{"2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"baseline/ee/",{"2d_MetVSdPhiMetNearLep","2d_MetVSdPhiMetNearLep_Puppi"}},
      {"baseline/emu/",{"2d_MetVSdPhiMetNearLep","2d_MetVSdPhiMetNearLep_Puppi"}},
      {"baseline/mumu/",{"2d_MetVSdPhiMetNearLep","2d_MetVSdPhiMetNearLep_Puppi"}},
      };
   
   // The following can be used for plotting, which is currently done with an extra module 
   // ~TCanvas can_2d;
   // ~for (auto const &sPresel_vVars:msPresel_vVars2D){
      // ~TString const &sPresel=sPresel_vVars.first;
      // ~for (TString sVar:sPresel_vVars.second){
         // ~can_2d.cd();
         // ~can_2d.SetLogz();
         // ~gPad->SetRightMargin(0.2);
         // ~gPad->SetLeftMargin(0.13);
         // ~gPad->SetBottomMargin(0.10);
         // ~TString loc=sPresel+sVar;
         // ~TH2F *hist=hs2d.getHistogram(loc,"TTbar");
         
         // ~hist->GetYaxis()->SetTitleOffset(1.3);
         // ~hist->GetXaxis()->SetTitleOffset(0.9);
         // ~hist->GetZaxis()->SetTitleOffset(1.3);
         // ~hist->GetYaxis()->SetTitleSize(0.05);
         // ~hist->GetXaxis()->SetTitleSize(0.05);
         // ~hist->GetZaxis()->SetTitleSize(0.05);
         // ~hist->GetYaxis()->SetLabelSize(0.04);
         // ~hist->GetXaxis()->SetLabelSize(0.04);
         // ~hist->GetZaxis()->SetLabelSize(0.04);
                  
         // ~hist->SetStats(0);
         // ~hist->Draw("colz");
         // ~TString cat;
         // ~if (sPresel.Contains("ee/")) cat="ee";
         // ~else if (sPresel.Contains("emu/")) cat="e#mu";
         // ~else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         // ~TLatex label=gfx::cornerLabel(cat,2);
         // ~label.Draw();
         // ~saver.save(can_2d,"tt_only/"+loc);
      // ~}
   // ~}
   
   //Save 2d histograms
   saveHistograms2D(msPresel_vVars2D,saver_hist,hs2d,samplesToCombine);
   
   
   // ~//Print efficiencies
   // ~std::ofstream out;
   // ~out.open(TString::Format("../output/txt/efficiencies%.1f.txt",cfg.processFraction*100));
   // ~out.precision(4);
   // ~for (TString sample : samplesToCombine){
      // ~out<<sample<<"   "<<count[sample]/(1.*Ngen[sample])<<std::endl;
   // ~}
   // ~out.close();
   
}
