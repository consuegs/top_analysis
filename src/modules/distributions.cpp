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
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   
   for(TString selection:{"baseline","baseline_Met200","baseline_Met400"}){ //Reco 1D Histograms
      hs.addHist(selection+"/ee/met"   ,";%MET;EventsBIN"           ,100,0,600);
      hs.addHist(selection+"/emu/met"   ,";%MET;EventsBIN"           ,100,0,600);
      hs.addHist(selection+"/mumu/met"   ,";%MET;EventsBIN"           ,100,0,600);
      
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
      
      hs.addHist(selection+"/ee/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/emu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"/mumu/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
      
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
   
   
   for(TString selection:{"genParticles","genParticles_Met200"}){ //1D Gen Histograms
      hs.addHist(selection+"/ee/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,1000);
      
      hs.addHist(selection+"/ee/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu(+BSM)}-genMET(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu(+BSM)}-genMET(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu(+BSM)}-genMET(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu(+BSM)}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu(+BSM)}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu(+BSM)}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
      hs.addHist(selection+"/emu/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
      hs.addHist(selection+"/mumu/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
      
      hs.addHist(selection+"/ee/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
      hs.addHist(selection+"/emu/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
      hs.addHist(selection+"/mumu/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
      
      hs.addHist(selection+"/ee/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-5,5);
      hs.addHist(selection+"/emu/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-5,5);
      hs.addHist(selection+"/mumu/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-5,5);
      
      hs.addHist(selection+"/ee/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-5,5);
      hs.addHist(selection+"/emu/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-5,5);
      hs.addHist(selection+"/mumu/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-5,5);
      
      hs.addHist(selection+"/ee/diff_Met_genMET_normSUM"   ,";#frac{p_{T}^{miss}-genMET}{p_{T}^{miss}+genMET};EventsBIN"           ,200,-2,2);
      hs.addHist(selection+"/emu/diff_Met_genMET_normSUM"   ,";#frac{p_{T}^{miss}-genMET}{p_{T}^{miss}+genMET};EventsBIN"           ,200,-2,2);
      hs.addHist(selection+"/mumu/diff_Met_genMET_normSUM"   ,";#frac{p_{T}^{miss}-genMET}{p_{T}^{miss}+genMET};EventsBIN"           ,200,-2,2);
      
      hs.addHist(selection+"/ee/diff_Met_DMgenMET_normSUM"   ,";#frac{p_{T}^{miss}-DMgenMET}{p_{T}^{miss}+DMgenMET};EventsBIN"           ,200,-2,2);
      hs.addHist(selection+"/emu/diff_Met_DMgenMET_normSUM"   ,";#frac{p_{T}^{miss}-DMgenMET}{p_{T}^{miss}+DMgenMET};EventsBIN"           ,200,-2,2);
      hs.addHist(selection+"/mumu/diff_Met_DMgenMET_normSUM"   ,";#frac{p_{T}^{miss}-DMgenMET}{p_{T}^{miss}+DMgenMET};EventsBIN"           ,200,-2,2);
      
      hs.addHist(selection+"/ee/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_genMT2_MT2"   ,";genMT2-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_genMT2_MT2"   ,";genMT2-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_genMT2_MT2"   ,";genMT2-MT2(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_genMT2neutrino_MT2"   ,";genMT2neutrino-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_genMT2neutrino_MT2"   ,";genMT2neutrino-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_genMT2neutrino_MT2"   ,";genMT2neutrino-MT2(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_genMT_MT"   ,";gen M_{T}(genMET,l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_genMT_MT"   ,";gen M_{T}(genMET,l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_genMT_MT"   ,";gen M_{T}(genMET,l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_genMTneutrino_MT"   ,";gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/emu/diff_genMTneutrino_MT"   ,";gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"/mumu/diff_genMTneutrino_MT"   ,";gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1})-M_{T}(p_{T}^{miss},l_{1})(GeV);EventsBIN"           ,400,-100,100);
      
      hs.addHist(selection+"/ee/diff_dPhiMetNearLep_gen"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(genMet,nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"/emu/diff_dPhiMetNearLep_gen"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(genMet,nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"/mumu/diff_dPhiMetNearLep_gen"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(genMet,nearest l);;EventsBIN" ,400,-3.2,3.2);
      
      hs.addHist(selection+"/ee/diff_dPhiMetNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"/emu/diff_dPhiMetNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"/mumu/diff_dPhiMetNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);;EventsBIN" ,400,-3.2,3.2);
      
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
   }//End 1D Gen Histograms
   
   //Define 2D histograms in the following
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   hs2d.addHist("genParticles/ee/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_nunuVSMet1000", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("genParticles/emu/2d_nunuVSMet1000", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("genParticles/mumu/2d_nunuVSMet1000", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,1000,100,0,1000);
   
   hs2d.addHist("genParticles/ee/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles/emu/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles/mumu/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   
   hs2d.addHist("genParticles/ee/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles/emu/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles/mumu/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   hs2d.addHist("baseline/emu/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   hs2d.addHist("baseline/emu/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   hs2d.addHist("baseline/mumu/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   
   hs2d.addHist("baseline/ee/2d_MetVSMT", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("baseline/emu/2d_MetVSMT", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("baseline/mumu/2d_MetVSMT", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   
   hs2d.addHist("baseline/ee/2d_MetVSMT_nextLep", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("baseline/emu/2d_MetVSMT_nextLep", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   hs2d.addHist("baseline/mumu/2d_MetVSMT_nextLep", ";p_{T}^{miss} (GeV);M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN" ,100,0,1000,100,0,1000);
   
   hs2d.addHist("genParticles/ee/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/emu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/mumu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   // Only MET>200 GeV
   hs2d.addHist("genParticles_Met200/ee/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_MT2VSgenMT2", ";MT2 (GeV);genMT2(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_MT2VSgenMT2neutrino", ";MT2 (GeV);genMT2neutrino(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_MT_l1METVSgenMT", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(genMET,l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/emu/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles_Met200/mumu/2d_MT_l1METVSgenMTneutrino", ";M_{T}(p_{T}^{miss},l_{1}) (GeV);gen M_{T}(p_{T}^{#nu#nu(+BSM)},l_{1}) (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles_Met200/ee/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles_Met200/emu/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles_Met200/mumu/2d_dPhiMetNearLep_Response", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   
   hs2d.addHist("genParticles_Met200/ee/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles_Met200/emu/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   hs2d.addHist("genParticles_Met200/mumu/2d_dPhiMetNearLep_genResponse", ";|#Delta#phi|(p_{T}^{miss},nearest l);|#Delta#phi|(genMet,nearest l);EventsBIN" ,100,0,3.2,100,0,3.2);
   
   hs2d.addHist("baseline_Met200/ee/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   hs2d.addHist("baseline_Met200/emu/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   hs2d.addHist("baseline_Met200/mumu/2d_MT2VSdPhiMetNearLep", ";MT2 (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,600,100,0,3.2);
   
   hs2d.addHist("baseline_Met200/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline_Met200/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("baseline_Met200/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   hs2d.addHist("baseline_Met200/ee/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   hs2d.addHist("baseline_Met200/emu/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   hs2d.addHist("baseline_Met200/mumu/2d_MetVSCosdPhiMetNearLep", ";p_{T}^{miss} (GeV);cos(|#Delta#phi|(p_{T}^{miss}),nearest l);EventsBIN" ,100,0,1000,100,-1,1);
   
   //Ntuple and file to save minimal ttbar tree used for binning studies
   float minTree_MET, minTree_PtNuNu, minTree_PhiRec, minTree_PhiGen, minTree_PhiNuNu, minTree_PhiMetNearJet, minTree_PhiMetFarJet, minTree_PhiMetLeadJet, minTree_PhiMetLead2Jet,
   minTree_PhiMetbJet, minTree_PhiLep1Lep2, minTree_METsig, minTree_N, minTree_SF, minTree_genMet, minTree_PuppiMet, minTree_HT, minTree_MT, minTree_genMT, minTree_MT_nextLep, minTree_genMT_nextLep;
   UInt_t minTree_runNo, minTree_lumNo, minTree_genDecayMode, minTree_n_Interactions;
   ULong64_t minTree_evtNo;
   io::RootFileSaver ttbar_res_saver(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res"+met_sf_string+"%.1f.root",cfg.processFraction*100),TString::Format("ttbar_res%.1f",cfg.processFraction*100),true,false);
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
   ttbar_res.Branch("dPhiLep1Lep2",&minTree_PhiLep1Lep2,"dPhiLep1Lep2/f");
   ttbar_res.Branch("METsig",&minTree_METsig,"METsig/f");
   ttbar_res.Branch("N",&minTree_N,"N/f");
   ttbar_res.Branch("SF",&minTree_SF,"SF/f");
   ttbar_res.Branch("runNo",&minTree_runNo,"runNo/i");
   ttbar_res.Branch("lumNo",&minTree_lumNo,"lumNo/i");
   ttbar_res.Branch("evtNo",&minTree_evtNo,"evtNo/l");
   ttbar_res.Branch("genDecayMode",&minTree_genDecayMode,"genDecayMode/i");
   ttbar_res.Branch("genMET",&minTree_genMet,"genMET/f");
   ttbar_res.Branch("PuppiMET",&minTree_PuppiMet,"PuppiMET/f");
   ttbar_res.Branch("HT",&minTree_HT,"HT/f");
   ttbar_res.Branch("MT",&minTree_MT,"MT/f");
   ttbar_res.Branch("genMT",&minTree_genMT,"genMT/f");
   ttbar_res.Branch("MT_nextLep",&minTree_MT_nextLep,"MT_nextLep/f");
   ttbar_res.Branch("genMT_nextLep",&minTree_genMT_nextLep,"genMT_nextLep/f");
   ttbar_res.Branch("n_Interactions",&minTree_n_Interactions,"n_Interactions/i");
   
   //Additional map to calculate signal efficiencies
   std::map<TString,float> count;
   std::map<TString,float> Ngen;

   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");

      //~ bool const isData=dss.isData;
      //~ bool const isSignal=dss.isSignal;
      
      hs.setCurrentSample(dss.name);
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
      
      //Check if current sample is TTbar madGraph (has extra genMet>150 part)
      bool ttBar_madGraph=false;
      if (dss.datasetName=="TTbar_madGraph") ttBar_madGraph=true;
      
      //Check if current sample is TTbar madGraph high MET
      bool ttBar_madGraph150=false;
      if (dss.datasetName=="TTbar_madGraph150") ttBar_madGraph150=true;
      
      //Check if current sample is selectedSusy scenario 
      bool SUSY_T2tt_650_350=false;
      if (dss.datasetName=="T2tt_650_350") SUSY_T2tt_650_350=true;

      
      //Check if current sample is SUSY scenario
      bool SUSY_scen=false;
      if (dss.name.find("SMS")!=std::string::npos) SUSY_scen=true;

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
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
      TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
      TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
      TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
      TTreeReaderValue<int> n_Interactions(reader, "true_nPV");
      TTreeReaderValue<float> HTgen(reader, "genHt");
      TTreeReaderValue<bool> is_ee   (reader, "ee");
      TTreeReaderValue<bool> is_emu   (reader, "emu");
      TTreeReaderValue<bool> is_mumu   (reader, "mumu");
      TTreeReaderValue<float> mll   (reader, "mll");
      TTreeReaderValue<float> mt2   (reader, "MT2");
      TTreeReaderValue<float> genMT2   (reader, "genMT2");
      TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
      TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
      TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
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
         
         float fEventWeight=*w_pu * *w_mc;     //Set event weight 
         float SFWeight=*sf_lep1 * *sf_lep2;     //Set combined SF weight
         hs.setFillWeight(fEventWeight*SFWeight);
         hs2d.setFillWeight(fEventWeight*SFWeight);
         
         float met=MET->p.Pt();
         float const genMet=GENMET->p.Pt();
         
         //Booleans for reco and pseudo selection
         bool rec_selection=true;
         bool pseudo_selection=true;
         
         //For ttBar madGraph sample, only use genMet<150 due to extension
         if (ttBar_madGraph and genMet>150) continue;
         
         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         
         if (*is_ee){
            if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) rec_selection=false; //currently double check since trees only have tight leptons!!
            if((*electrons)[0].etaSC>2.4 || (*electrons)[1].etaSC>2.4) rec_selection=false; //To use same region as for muons, cut on supercluster eta
            p_l1=(*electrons)[0].p;
            p_l2=(*electrons)[1].p;
         }
         else if (*is_mumu){
            if(!(*muons)[0].isTight || !(*muons)[1].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 || (*muons)[1].rIso>0.15) rec_selection=false;
            p_l1=(*muons)[0].p;
            p_l2=(*muons)[1].p;
         }
         else if (*is_emu){
            if(!(*muons)[0].isTight || !(*electrons)[0].isTight) rec_selection=false;
            if((*muons)[0].rIso>0.15 ) rec_selection=false;
            if((*electrons)[0].etaSC>2.4 ) rec_selection=false;
            if ((*muons)[0].p.Pt()>(*electrons)[0].p.Pt()){
               p_l1=(*muons)[0].p;
               p_l2=(*electrons)[0].p;
            }
            else {
               p_l1=(*electrons)[0].p;
               p_l2=(*muons)[0].p;
            }
         }
         
         if (p_l1.Pt()<25 || p_l2.Pt()<20) rec_selection=false; //eta cuts already done in TreeWriter
         if (*mll<20 || ((*is_ee || *is_mumu) && *mll<106 && *mll>76)) rec_selection=false;
         if ((*is_ee || *is_mumu) && met<40) rec_selection=false;
         
         std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
         if (cjets.size()<2) rec_selection=false;
         
         bool bTag=false;
         std::vector<tree::Jet> BJets;
         for (tree::Jet const &jet : cjets) {
            // ~if (jet.bTagCSVv2>0.5426) {      //Loose working point for CSVv2 (Should be replaced in the future by deep CSV!!!)
            if (jet.bTagDeepCSV>0.2217) {      //Loose working point for deepCSV
               bTag=true;
               BJets.push_back(jet);
            }
         }
         if (!bTag) rec_selection=false; // end reco baseline selection
         
         if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
         
         if(rec_selection==false && pseudo_selection==false) continue;  // only proceed with events selected by one of the baseline selection (reco or pseudo)
         
         //Scale measured met by met_sf from config
         met=met_sf*met;
         
         // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
         TLorentzVector neutrinoPair(0,0,0,0);
         neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
         
         //Calculate MET before parton shower and hadronization for DM scenario and SUSY scenarios if gen level selection is true
         if (pseudo_selection==true) {
            for (auto const &genParticle : *genParticles){
               if((!SUSY_scen && abs(genParticle.pdgId)>100000) || (SUSY_scen && abs(genParticle.pdgId)==1000022)){
                  neutrinoPair+=genParticle.p;
               }
            }
         }
         
         //Get DeltaPhi between MET (or genMet or neutrino pT) and nearest Lepton
         float dPhiMETnearLep=4;
         float mt_MetNextLep=0; 
         float mt_NuNuNextLep=0; 
         for (TLorentzVector const lep : {p_l1,p_l2}){
            const float dPhi=MET->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
               dPhiMETnearLep=dPhi;
               mt_MetNextLep=phys::M_T(MET->p,lep);
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
         float dPhiLep1Lep2=4;
         float HT=0;
         if (rec_selection){ 
            for (tree::Jet const &jet : cjets) {
               const float dPhi=MET->p.DeltaPhi(jet.p);
               if (std::abs(dPhi) < std::abs(dPhiMETnearJet)) dPhiMETnearJet=dPhi;
               if (std::abs(dPhi) > std::abs(dPhiMETfarJet)) dPhiMETfarJet=dPhi;
               HT+=jet.p.Pt();
            }
            dPhiMETleadJet=MET->p.DeltaPhi(cjets[0].p);
            dPhiMETlead2Jet=MET->p.DeltaPhi(cjets[1].p);
            dPhiMetBJet=MET->p.DeltaPhi(BJets[0].p);
            dPhiLep1Lep2=p_l1.DeltaPhi(p_l2);
         }
         
         //Further variables
         float mt_MetLep1=phys::M_T(MET->p,p_l1);
         float mt_MetLep2=phys::M_T(MET->p,p_l2);
         float sum_mlb=phys::sumMlb(p_l1,p_l2,cjets,BJets);
         float conMt_Lep1Lep2=phys::conM_T(p_l1,p_l2);
         
         TLorentzVector leadGenLepton;
         if (genLepton->Pt()>genAntiLepton->Pt()) leadGenLepton=*genLepton;
         else leadGenLepton=*genAntiLepton;
         
         float mt_genMetLep1=phys::M_T(GENMET->p,leadGenLepton);
         float mt_genNeutrinosLep1=phys::M_T(neutrinoPair,leadGenLepton);
         
         //Fill minimal tree for TTbar resolution used in binning studies
         if (ttBar_dilepton || ttBar_madGraph || ttBar_madGraph150 || SUSY_T2tt_650_350 || ttBar_standard){
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
            minTree_PhiLep1Lep2=abs(dPhiLep1Lep2);
            minTree_METsig=MET->sig;
            minTree_N=lumi_weight*fEventWeight;
            minTree_SF=SFWeight;
            minTree_runNo=*runNo;
            minTree_lumNo=*lumNo;
            minTree_evtNo=*evtNo;
            minTree_genDecayMode=*genDecayMode;
            minTree_genMet=genMet;
            minTree_PuppiMet=MET_Puppi->p.Pt();
            minTree_HT=HT;
            minTree_MT=mt_MetLep1;
            minTree_genMT=mt_genNeutrinosLep1;
            minTree_n_Interactions=*n_Interactions;
            minTree_MT_nextLep=mt_MetNextLep;
            minTree_genMT_nextLep=mt_NuNuNextLep;
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
               minTree_HT=-1.;
               minTree_MT=-1.;
               minTree_MT_nextLep=-1.;
               minTree_SF=0.;
            }
            else if (pseudo_selection==false) {
               minTree_PtNuNu=-1.;
               minTree_PhiGen=-1.;
               minTree_PhiNuNu=-1.;
               minTree_genDecayMode=0.;
               minTree_genMet=-1.;
               minTree_n_Interactions=0;
               minTree_genMT=-1;
               minTree_genMT_nextLep=-1;
            }
            // ~std::cout<<minTree_lumNo<<std::endl;
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
         
         //Sort genTop pT
         std::vector<TLorentzVector> gen_tops;
         gen_tops.push_back(*genTop);
         gen_tops.push_back(*genAntiTop);
         sort(gen_tops.begin(), gen_tops.end(), tree::PtGreaterLorentz);
         float pT_top1=0;
         float pT_top2=0;
         pT_top1=gen_tops[0].Pt();
         pT_top2=gen_tops[1].Pt();
         
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
         hs.fill("genParticles/"+path_cat+"/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
         hs.fill("genParticles/"+path_cat+"/diff_ptNuNu_DMgenMET",neutrinoPair.Pt()-DMgenMET.Pt());
         hs.fill("genParticles/"+path_cat+"/diff_Met_genMET",met-genMet);
         hs.fill("genParticles/"+path_cat+"/diff_Met_DMgenMET",met-DMgenMET.Pt());
         hs.fill("genParticles/"+path_cat+"/diff_Met_genMET_norm",(met-genMet)/genMet);
         hs.fill("genParticles/"+path_cat+"/diff_Met_DMgenMET_norm",(met-DMgenMET.Pt())/DMgenMET.Pt());
         hs.fill("genParticles/"+path_cat+"/diff_Met_genMET_normSUM",(met-genMet)/(met+genMet));
         hs.fill("genParticles/"+path_cat+"/diff_Met_DMgenMET_normSUM",(met-DMgenMET.Pt())/(met+DMgenMET.Pt()));
         hs.fill("genParticles/"+path_cat+"/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
         hs.fill("genParticles/"+path_cat+"/diff_genMT2_MT2",*genMT2-*mt2);
         hs.fill("genParticles/"+path_cat+"/diff_genMT2neutrino_MT2",*genMT2neutrino-*mt2);
         hs.fill("genParticles/"+path_cat+"/diff_genMT_MT",mt_genMetLep1-mt_MetLep1);
         hs.fill("genParticles/"+path_cat+"/diff_genMTneutrino_MT",mt_genNeutrinosLep1-mt_MetLep1);
         hs.fill("genParticles/"+path_cat+"/diff_dPhiMetNearLep_gen",abs(dPhigenMETnearLep)-abs(dPhiMETnearLep));
         hs.fill("genParticles/"+path_cat+"/diff_dPhiMetNearLep",abs(dPhiPtNunearLep)-abs(dPhiMETnearLep));
         hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep1));
         hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep2));
         hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep1));
         hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep2));
         hs.fill("genParticles/"+path_cat+"/pTtop1",pT_top1);
         hs.fill("genParticles/"+path_cat+"/pTtop2",pT_top2);
         hs2d.fill("genParticles/"+path_cat+"/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
         hs2d.fill("genParticles/"+path_cat+"/2d_nunuVSDMgenMet",DMgenMET.Pt(),neutrinoPair.Pt());
         hs2d.fill("genParticles/"+path_cat+"/MetVSgenMet",met,genMet);
         hs2d.fill("genParticles/"+path_cat+"/MetVSDMgenMet",met,DMgenMET.Pt());
         hs2d.fill("genParticles/"+path_cat+"/2d_nunuVSMet",met,neutrinoPair.Pt());
         hs2d.fill("genParticles/"+path_cat+"/2d_nunuVSMet1000",met,neutrinoPair.Pt());
         hs2d.fill("genParticles/"+path_cat+"/2d_MT2VSgenMT2",*mt2,*genMT2);
         hs2d.fill("genParticles/"+path_cat+"/2d_MT2VSgenMT2neutrino",*mt2,*genMT2neutrino);
         hs2d.fill("genParticles/"+path_cat+"/2d_MT_l1METVSgenMT",mt_MetLep1,mt_genMetLep1);
         hs2d.fill("genParticles/"+path_cat+"/2d_MT_l1METVSgenMTneutrino",mt_MetLep1,mt_genNeutrinosLep1);
         hs2d.fill("genParticles/"+path_cat+"/2d_dPhiMetNearLep_genResponse",abs(dPhiMETnearLep),abs(dPhigenMETnearLep));
         hs2d.fill("genParticles/"+path_cat+"/2d_dPhiMetNearLep_Response",abs(dPhiMETnearLep),abs(dPhiPtNunearLep));
         hs2d.fill("baseline/"+path_cat+"/2d_MT2VSdPhiMetNearLep",*mt2,abs(dPhiMETnearLep));
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSCosdPhiMetNearLep",met,TMath::Cos(abs(dPhiMETnearLep)));
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSMT",met,mt_MetLep1);
         hs2d.fill("baseline/"+path_cat+"/2d_MetVSMT_nextLep",met,mt_MetNextLep);
         hs2d.fill("genParticles/"+path_cat+"/2d_PtNuNuVSdPhiNuNuNearLep",neutrinoPair.Pt(),abs(dPhiPtNunearLep));
         
         if (met>200){
            hs.fill("baseline_Met200/"+path_cat+"/met",met);
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
            hs.fill("genParticles_Met200/"+path_cat+"/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_ptNuNu_DMgenMET",neutrinoPair.Pt()-DMgenMET.Pt());
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_genMET",met-genMet);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_DMgenMET",met-DMgenMET.Pt());
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_genMET_norm",(met-genMet)/genMet);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_DMgenMET_norm",(met-DMgenMET.Pt())/DMgenMET.Pt());
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_genMET_normSUM",(met-genMet)/(met+genMet));
            hs.fill("genParticles_Met200/"+path_cat+"/diff_Met_DMgenMET_normSUM",(met-DMgenMET.Pt())/(met+DMgenMET.Pt()));
            hs.fill("genParticles_Met200/"+path_cat+"/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_genMT2_MT2",*genMT2-*mt2);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_genMT2neutrino_MT2",*genMT2neutrino-*mt2);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_genMT_MT",mt_genMetLep1-mt_MetLep1);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_genMTneutrino_MT",mt_genNeutrinosLep1-mt_MetLep1);
            hs.fill("genParticles_Met200/"+path_cat+"/diff_dPhiMetNearLep_gen",abs(dPhigenMETnearLep)-abs(dPhiMETnearLep));
            hs.fill("genParticles_Met200/"+path_cat+"/diff_dPhiMetNearLep",abs(dPhiPtNunearLep)-abs(dPhiMETnearLep));
            hs.fill("genParticles_Met200/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep1));
            hs.fill("genParticles_Met200/"+path_cat+"/dphi_NeutrinoLep",abs(dPhiNeutrinoLep2));
            hs.fill("genParticles_Met200/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep1));
            hs.fill("genParticles_Met200/"+path_cat+"/dR_NeutrinoLep",abs(dRNeutrinoLep2));
            hs.fill("genParticles_Met200/"+path_cat+"/pTtop1",pT_top1);
            hs.fill("genParticles_Met200/"+path_cat+"/pTtop2",pT_top2);
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_nunuVSDMgenMet",DMgenMET.Pt(),neutrinoPair.Pt());
            hs2d.fill("genParticles_Met200/"+path_cat+"/MetVSgenMet",met,genMet);
            hs2d.fill("genParticles_Met200/"+path_cat+"/MetVSDMgenMet",met,DMgenMET.Pt());
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_nunuVSMet",met,neutrinoPair.Pt());
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_MT2VSgenMT2",*mt2,*genMT2);
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_MT2VSgenMT2neutrino",*mt2,*genMT2neutrino);
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_MT_l1METVSgenMT",mt_MetLep1,mt_genMetLep1);
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_MT_l1METVSgenMTneutrino",mt_MetLep1,mt_genNeutrinosLep1);
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_dPhiMetNearLep_genResponse",abs(dPhiMETnearLep),abs(dPhigenMETnearLep));
            hs2d.fill("genParticles_Met200/"+path_cat+"/2d_dPhiMetNearLep_Response",abs(dPhiMETnearLep),abs(dPhiPtNunearLep));
            hs2d.fill("baseline_Met200/"+path_cat+"/2d_MT2VSdPhiMetNearLep",*mt2,abs(dPhiMETnearLep));
            hs2d.fill("baseline_Met200/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
            hs2d.fill("baseline_Met200/"+path_cat+"/2d_MetVSCosdPhiMetNearLep",met,TMath::Cos(abs(dPhiMETnearLep)));
            
            if (met>400){
               hs.fill("baseline_Met400/"+path_cat+"/met",met);
               hs.fill("baseline_Met400/"+path_cat+"/met1000",met);
               hs.fill("baseline_Met400/"+path_cat+"/mll",*mll);
               hs.fill("baseline_Met400/"+path_cat+"/pTlep1",p_l1.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/pTlep2",p_l2.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/pTsumlep",(p_l1+p_l2).Pt());
               hs.fill("baseline_Met400/"+path_cat+"/sumpTlep",p_l1.Pt()+p_l2.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/pTbJet",BJets[0].p.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/pTJet1",cjets[0].p.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/pTJet2",cjets[1].p.Pt());
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metJet",abs(dPhiMETnearJet));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metNearLep",abs(dPhiMETnearLep));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metLeadJet",abs(dPhiMETleadJet));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metLead2Jet",abs(dPhiMETlead2Jet));
               hs.fill("baseline_Met400/"+path_cat+"/COSdphi_metNearLep",TMath::Cos(abs(dPhiMETnearLep)));
               hs.fill("baseline_Met400/"+path_cat+"/SINdphi_metNearLep",TMath::Sin(abs(dPhiMETnearLep)));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metBJet",abs(dPhiMetBJet));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_bJetLep1",abs(dPhiLep1BJet));
               hs.fill("baseline_Met400/"+path_cat+"/dR_bJetLep1",abs(dRLep1BJet));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_bJetLep2",abs(dPhiLep2BJet));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_bJetnearLep",dphi_bJetnearLep);
               hs.fill("baseline_Met400/"+path_cat+"/dphi_b1b2",abs(dPhiBjets));
               hs.fill("baseline_Met400/"+path_cat+"/dR_b1b2",abs(dRBjets));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metLep1",abs(dPhiLep1MET));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metLep2",abs(dPhiLep2MET));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_metLepsum",abs(dPhiMetLepSum));
               hs.fill("baseline_Met400/"+path_cat+"/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
               hs.fill("baseline_Met400/"+path_cat+"/dR_Lep1Lep2",abs(dR_Lep1Lep2));
               hs.fill("baseline_Met400/"+path_cat+"/nJets",cjets.size());
               hs.fill("baseline_Met400/"+path_cat+"/nBjets",nBjets);
               hs.fill("baseline_Met400/"+path_cat+"/mt2",*mt2);
               hs.fill("baseline_Met400/"+path_cat+"/mt_MetLep1",mt_MetLep1);
               hs.fill("baseline_Met400/"+path_cat+"/mt_MetLep2",mt_MetLep2);
               hs.fill("baseline_Met400/"+path_cat+"/mt_MetNextLep",mt_MetNextLep);
               hs.fill("baseline_Met400/"+path_cat+"/conMt_Lep1Lep2",conMt_Lep1Lep2);
               hs.fill("baseline_Met400/"+path_cat+"/ST",ST);
               hs.fill("baseline_Met400/"+path_cat+"/HT",HT);
               hs.fill("baseline_Met400/"+path_cat+"/sum_STHT",ST+HT);
               hs.fill("baseline_Met400/"+path_cat+"/sum_mlb",sum_mlb);
            }
         }
               
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs2d.scaleLumi();
      hs.mergeOverflow();
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
      else if (ttBar_dilepton) {
         ttbar_res_saver.save(ttbar_res,"ttbar_res_dilepton");
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
   hs2d.combineFromSubsamples(samplesToCombine);
   
   //Combine ttBar madGraph with high genMet sample
   hs.combineSamples("TTbar_madGraphCOMB",{"TTbar_madGraph","TTbar_madGraph150"});
   hs2d.combineSamples("TTbar_madGraphCOMB",{"TTbar_madGraph","TTbar_madGraph150"});
   samplesToCombine.push_back("TTbar_madGraphCOMB");
   
   //Plotting part 1D
   io::RootFileSaver saver(TString::Format("plots"+met_sf_string+"%.1f.root",cfg.processFraction*100),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/emu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metBJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/mumu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metBJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/ee/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/emu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/mumu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met400/ee/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met400/emu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met400/mumu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTbJet","pTJet1","pTJet2","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"genParticles/ee/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles/emu/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles/mumu/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/ee/",{"diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_dPhiMetNearLep_gen","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/emu/",{"diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_dPhiMetNearLep_gen","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/mumu/",{"diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM","diff_Met_DMgenMET_normSUM","diff_ptNuNu_Met","diff_dPhiMetNearLep_gen","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_genMT_MT","diff_genMTneutrino_MT","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
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
   
   //Plotting part 2D
   std::map<TString,std::vector<TString>> msPresel_vVars2D={
      {"genParticles/ee/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSgenMet","MetVSDMgenMet","2d_nunuVSMet","2d_nunuVSMet1000","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response","2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"genParticles/emu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet","2d_nunuVSMet1000","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response","2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"genParticles/mumu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet","2d_nunuVSMet1000","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response","2d_PtNuNuVSdPhiNuNuNearLep"}},
      {"genParticles_Met200/ee/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSgenMet","MetVSDMgenMet","2d_nunuVSMet","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response"}},
      {"genParticles_Met200/emu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response"}},
      {"genParticles_Met200/mumu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino","2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response"}},
      {"baseline/ee/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep","2d_MetVSMT","2d_MetVSMT_nextLep"}},
      {"baseline/emu/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep","2d_MetVSMT","2d_MetVSMT_nextLep"}},
      {"baseline/mumu/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep","2d_MetVSMT","2d_MetVSMT_nextLep"}},
      {"baseline_Met200/ee/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep"}},
      {"baseline_Met200/emu/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep"}},
      {"baseline_Met200/mumu/",{"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep","2d_MetVSCosdPhiMetNearLep"}},
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
