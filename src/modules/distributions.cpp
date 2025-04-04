#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"
#include "tools/jetCorrections.hpp"
#include "tools/bTagWeights.hpp"
#include "tools/leptonSF.hpp"
#include "tools/leptonCorrections.hpp"
#include "tools/triggerSF.hpp"
#include "tools/mcWeights.hpp"
#include "tools/dnnRegression.hpp"
#include "tools/jetPileupIDSF.hpp"
#include "tools/jetVetoMaps.hpp"

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
#include <chrono>
#include <numeric>

using namespace std::chrono;

Config const &cfg=Config::get();

struct dim2D {
   TString name;
   TString axisTitle;
   float xMin;
   float xMax;
   int nBins;
   float value = 0;
};

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   TString dssName_multi="";
   
   //Read systematic from command line
   Systematic::Systematic currentSystematic(cfg.systematic);
   //Set boolean for running nominal (also used when testing different selections)
   bool isNominal = (std::find(Systematic::nominalTypes.begin(), Systematic::nominalTypes.end(), currentSystematic.type()) != Systematic::nominalTypes.end());
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   hist::Histograms<TH1F> hs_cutflow(vsDatasubsets);
   
   hs_cutflow.addHist("cutflow/ee"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   hs_cutflow.addHist("cutflow/emu"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   hs_cutflow.addHist("cutflow/mumu"   ,";cut;EventsBIN"           ,9,0.5,9.5);
   
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee","/mumu","/emu"}){
         hs.addHist(selection+channel+"/MET"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/PuppiMET"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/PuppiMET_xy"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/MET_xy"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/DNN_MET_pT"   ,";DNN_MET_pT (GeV);EventsBIN"           ,500,0,500);
         hs.addHist(selection+channel+"/DNN_MET_dPhi_nextLep"   ,";DNN_MET_dPhi_nextLep;EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/DeepMET_reso"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/DeepMET_resp"   ,";%MET;EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/Lep1_pt"   ,";%pTl1;EventsBIN"           ,100,0,600);
         hs.addHist(selection+channel+"/Lep2_pt"   ,";%pTl2;EventsBIN"           ,100,0,600);
         hs.addHist(selection+channel+"/Lep_e_pt"   ,";%pTle;EventsBIN"           ,100,0,600);
         hs.addHist(selection+channel+"/Lep_mu_pt"   ,";%pTlmu;EventsBIN"           ,100,0,600);
         hs.addHist(selection+channel+"/pTsumlep"   ,";p_{T}^{ll} (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/sumpTlep"   ,";%pTl1+%pTl2 (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/pTbJet"   ,";p_{T}^{b} (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/bJet_eta"   ,";#eta^{b};EventsBIN"           ,50,-2.5,2.5);
         hs.addHist(selection+channel+"/Jet1_pt"   ,";p_{T}^{Jet1} (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/Jet2_pt"   ,";p_{T}^{Jet2} (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dphi_metNearLep_puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dphi_metNearLep_puppi_xy"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dphi_metNearLep_xy"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/COSdphi_metNearLep"   ,";cos(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,-1.,1);
         hs.addHist(selection+channel+"/SINdphi_metNearLep"   ,";sin(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,0.,1);
         hs.addHist(selection+channel+"/dR_bJetLep1"   ,";|#Delta R|(b Jet,l_{1});EventsBIN"           ,100,0,5);
         hs.addHist(selection+channel+"/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2});EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dphi_bJetnearLep"   ,";|#Delta#phi|(b Jet,next l);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dphi_b1b2"   ,";|#Delta#phi|(b Jet1,b Jet2);EventsBIN"           ,320,0,3.2);
         hs.addHist(selection+channel+"/dR_b1b2"   ,";|#Delta R|(b Jet1,b Jet2);EventsBIN"           ,100,0,5);
         hs.addHist(selection+channel+"/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1});EventsBIN"           ,100,0,3.2);
         hs.addHist(selection+channel+"/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2});EventsBIN"           ,100,0,3.2);
         hs.addHist(selection+channel+"/dphi_metLepsum"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}+l_{2});EventsBIN"           ,100,0,3.2);
         hs.addHist(selection+channel+"/dR_Lep1Lep2"   ,";|#Delta R|(l_{1},l_{2});EventsBIN"           ,100,0,5);
         hs.addHist(selection+channel+"/nBjets"   ,";N_{bJets};EventsBIN"           ,5,-0.5,4.5);
         hs.addHist(selection+channel+"/C_em_W_p"   ,";C_{em,W,+} (GeV);EventsBIN"           ,100,0,300);
         hs.addHist(selection+channel+"/C_em_W_m"   ,";C_{em,W,-} (GeV);EventsBIN"           ,100,0,300);
         hs.addHist(selection+channel+"/mt_MetLep2"   ,";M_{T}(p_{T}^{miss},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/mt_MetNextLep"   ,";M_{T}(p_{T}^{miss},nearest l) (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/conMt_Lep1Lep2"   ,";conM_{T}(l_{1},l_{2}) (GeV);EventsBIN"           ,100,0,1000);
         hs.addHist(selection+channel+"/ST"   ,";S_{T} (GeV);EventsBIN"           ,100,0,1500);
         hs.addHist(selection+channel+"/HT"   ,";H_{T} (GeV);EventsBIN"           ,100,0,2500);
         hs.addHist(selection+channel+"/sum_STHT"   ,";S_{T}+H_{T} (GeV);EventsBIN"           ,100,0,4000);
         hs.addHist(selection+channel+"/sum_mlb"   ,";sum m_{lb} (GeV);EventsBIN"           ,100,0,3000);
         hs.addHist(selection+channel+"/Lep1_phi"   ,";Lep1_phi;EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/Lep2_phi"   ,";Lep2_phi;EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/Jet1_phi"   ,";Jet1_phi;EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/Jet2_phi"   ,";Jet2_phi;EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/PFMET_phi"   ,";#phi (PFMET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/PuppiMET_phi"   ,";#phi (PuppiMET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/PuppiMET_xy_phi"   ,";#phi (PuppiMET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/MET_xy_phi"   ,";#phi (MET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/CaloMET"   ,";CaloMET (GeV);EventsBIN"           ,100,0,500);
         hs.addHist(selection+channel+"/CaloMET_phi"   ,";#phi (CaloMET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/DeepMET_reso_phi"   ,";#phi (DeepMET);EventsBIN"           ,100,-3.2,3.2);
         hs.addHist(selection+channel+"/DeepMET_resp_phi"   ,";#phi (DeepMET);EventsBIN"           ,100,-3.2,3.2);
         
         //DNN Inputs
         hs.addHist(selection+channel+"/PuppiMET*cos(PuppiMET_phi)"   ,";Puppi p_{x}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(selection+channel+"/PuppiMET*sin(PuppiMET_phi)"   ,";Puppi p_{y}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(selection+channel+"/PuppiMET_xy*cos(PuppiMET_xy_phi)"   ,";Puppi p_{x}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(selection+channel+"/PuppiMET_xy*sin(PuppiMET_xy_phi)"   ,";Puppi p_{y}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(selection+channel+"/METunc_Puppi"                 ,";Puppi #sigma(p_{t}^{miss}) (GeV);EventsBIN"           ,50,0,100);
         hs.addHist(selection+channel+"/MET*cos(PFMET_phi)"           ,";PF p_{x}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/MET*sin(PFMET_phi)"           ,";PF p_{y}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/MET_xy*cos(MET_xy_phi)"           ,";PF p_{x}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/MET_xy*sin(MET_xy_phi)"           ,";PF p_{y}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/DeepMET_reso*cos(DeepMET_reso_phi)",";Deep p_{x}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/DeepMET_reso*sin(DeepMET_reso_phi)",";Deep p_{y}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/DeepMET_resp*cos(DeepMET_resp_phi)",";Deep p_{x}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/DeepMET_resp*sin(DeepMET_resp_phi)",";Deep p_{y}^{miss} (GeV);EventsBIN"                      ,50,-250,250);
         hs.addHist(selection+channel+"/vecsum_pT_allJet*cos(HT_phi)" ,";vecsum_pT_allJet_{x} (GeV);EventsBIN"                 ,80,-400,400);
         hs.addHist(selection+channel+"/vecsum_pT_allJet*sin(HT_phi)" ,";vecsum_pT_allJet_{y} (GeV);EventsBIN"                 ,80,-400,400);
         hs.addHist(selection+channel+"/nJets"                        ,";N_{jets} ;EventsBIN"                                  ,11,1.5,12.5);
         hs.addHist(selection+channel+"/n_Interactions"               ,";N_{interactions} ;EventsBIN"                          ,50,0,100);
         hs.addHist(selection+channel+"/Lep1_flavor"                  ,";Flavor_{l1} ;EventsBIN"                               ,2,0.5,2.5);
         hs.addHist(selection+channel+"/Lep2_flavor"                  ,";Flavor_{l2} ;EventsBIN"                               ,2,0.5,2.5);
         hs.addHist(selection+channel+"/Lep1_pt*cos(Lep1_phi)"        ,";p_{x}^{lep1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Lep1_pt*sin(Lep1_phi)"        ,";p_{y}^{lep1} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Lep1_eta"                     ,";#eta^{lep1} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(selection+channel+"/Lep1_E"                       ,";E^{lep1} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(selection+channel+"/Lep2_pt*cos(Lep2_phi)"        ,";p_{x}^{lep2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Lep2_pt*sin(Lep2_phi)"        ,";p_{y}^{lep2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Lep2_eta"                     ,";#eta^{lep2} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(selection+channel+"/Lep2_E"                       ,";E^{lep2} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(selection+channel+"/Jet1_pt*cos(Jet1_phi)"        ,";p_{x}^{jet1} (GeV) ;EventsBIN"                        ,80,-400,400);
         hs.addHist(selection+channel+"/Jet1_pt*sin(Jet1_phi)"        ,";p_{y}^{jet1} (GeV) ;EventsBIN"                        ,80,-400,400);
         hs.addHist(selection+channel+"/Jet1_eta"                     ,";#eta^{jet1} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(selection+channel+"/Jet1_E"                       ,";E^{jet1} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(selection+channel+"/Jet2_pt*cos(Jet2_phi)"        ,";p_{x}^{jet2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Jet2_pt*sin(Jet2_phi)"        ,";p_{y}^{jet2} (GeV) ;EventsBIN"                        ,50,-250,250);
         hs.addHist(selection+channel+"/Jet2_eta"                     ,";#eta^{jet2} ;EventsBIN"                               ,50,-2.5,2.5);
         hs.addHist(selection+channel+"/Jet2_E"                       ,";E^{jet2} (GeV) ;EventsBIN"                            ,50,0,500);
         hs.addHist(selection+channel+"/dPhiMETnearJet"               ,";|#delta#phi(PF MET, near. jet)| ;EventsBIN"           ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETfarJet"                ,";|#delta#phi(PF MET, far. jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETleadJet"               ,";|#delta#phi(PF MET, lead jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETlead2Jet"              ,";|#delta#phi(PF MET, sublead jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETbJet"                  ,";|#delta#phi(PF MET, lead b jet)| ;EventsBIN"          ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiLep1Lep2"                 ,";|#delta#phi(lep1, lep2)| ;EventsBIN"                  ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiJet1Jet2"                 ,";|#delta#phi(jet1, jet2)| ;EventsBIN"                  ,32,0,3.2);
         hs.addHist(selection+channel+"/METsig"                       ,";PF MET sig. ;EventsBIN"                               ,20,0,200);
         hs.addHist(selection+channel+"/MHT"                          ,";M_{HT} (GeV);EventsBIN"                               ,100,0,2000);
         hs.addHist(selection+channel+"/MT"                           ,";M_{T} (GeV);EventsBIN"                                ,50,0,1000);
         hs.addHist(selection+channel+"/looseLeptonVeto"              ,";looseLeptonVeto ;EventsBIN"                           ,2,-0.5,1.5);
         hs.addHist(selection+channel+"/dPhiMETnearJet_Puppi"         ,";|#delta#phi(Puppi MET, near. jet)| ;EventsBIN"        ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETfarJet_Puppi"          ,";|#delta#phi(Puppi MET, far. jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETleadJet_Puppi"         ,";|#delta#phi(Puppi MET, lead jet)| ;EventsBIN"         ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETlead2Jet_Puppi"        ,";|#delta#phi(Puppi MET, sublead jet)| ;EventsBIN"      ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiMETbJet_Puppi"            ,";|#delta#phi(Puppi MET, lead b jet)| ;EventsBIN"       ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiLep1bJet"                 ,";|#delta#phi(lep1, lead b jet)| ;EventsBIN"            ,32,0,3.2);
         hs.addHist(selection+channel+"/dPhiLep1Jet1"                 ,";|#delta#phi(lep1, ljet1)| ;EventsBIN"                 ,32,0,3.2);
         hs.addHist(selection+channel+"/mLL"                          ,";m_{ll} (GeV) ;EventsBIN"                              ,50,0,200);
         hs.addHist(selection+channel+"/CaloMET*cos(CaloMET_phi)"     ,";Calo p_{x}^{miss} (GeV);EventsBIN"                    ,50,-250,250);
         hs.addHist(selection+channel+"/CaloMET*sin(CaloMET_phi)"     ,";Calo p_{y}^{miss} (GeV);EventsBIN"                    ,50,-250,250);
         hs.addHist(selection+channel+"/MT2"                          ,";MT_{2} (GeV);EventsBIN"                               ,50,0,200);
         hs.addHist(selection+channel+"/vecsum_pT_allJet"             ,";vecsum_pT_allJet (GeV);EventsBIN"                     ,50,0,500);
         hs.addHist(selection+channel+"/vecsum_pT_l1l2_allJet"        ,";vecsum_pT_l1l2_allJet (GeV);EventsBIN"                ,50,0,500);
         hs.addHist(selection+channel+"/mass_l1l2_allJet"             ,";mass_l1l2_allJet (GeV);EventsBIN"                     ,100,0,2000);
         hs.addHist(selection+channel+"/ratio_vecsumpTlep_vecsumpTjet",";ratio_vecsumpTlep_vecsumpTjet ;EventsBIN"             ,50,0,10);
         hs.addHist(selection+channel+"/mjj"                          ,";m_{jj} (GeV);EventsBIN"                               ,100,0,2000);
         
         hs.addHist(selection+channel+"/DNNMET*cos(DNNMET_phi)"   ,";DNN p_{x}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
         hs.addHist(selection+channel+"/DNNMET*sin(DNNMET_phi)"   ,";DNN p_{y}^{miss} (GeV);EventsBIN"                   ,50,-250,250);
      }
   }//End 1D reco histograms
   
   
   for(TString selection:{"genParticles"}){ //1D Gen Histograms
      hs.addHist(selection+"/ee/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/emu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/mumu/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,500);
      
      hs.addHist(selection+"/ee/genMET"   ,";genMET(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/emu/genMET"   ,";genMET(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/mumu/genMET"   ,";genMET(GeV);EventsBIN"           ,100,0,500);
      
      hs.addHist(selection+"/ee/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/emu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"/mumu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,500);
      
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
      
      hs.addHist(selection+"/ee/n_Interactions_gen"   ,";n_Interactions_gen;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/emu/n_Interactions_gen"   ,";n_Interactions_gen;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"/mumu/n_Interactions_gen"   ,";n_Interactions_gen;EventsBIN"           ,100,0,1000);
   }//End 1D Gen Histograms
   
   //Define 2D histograms in the following
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep_Puppi", ";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   
   // ~hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   // ~hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   // ~hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,100,0,500,320,0,3.2);
   hs2d.addHist("baseline/ee/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,800,0,400,160,0,3.2);
   hs2d.addHist("baseline/emu/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,800,0,400,160,0,3.2);
   hs2d.addHist("baseline/mumu/2d_MetVSdPhiMetNearLep_DNN", ";p_{T,DNN}^{miss} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,800,0,400,160,0,3.2);
   
   hs2d.addHist("genParticles/ee/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/emu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   hs2d.addHist("genParticles/mumu/2d_PtNuNuVSdPhiNuNuNearLep", ";p_{T}^{#nu#nu(+BSM)} (GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);EventsBIN" ,100,0,1000,100,0,3.2);
   
   hs2d.addHist("baseline/ee/2d_CemVSdPhiMetNearLep", ";C_{em,W,-} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,20,0,300,20,0,3.2);
   hs2d.addHist("baseline/emu/2d_CemVSdPhiMetNearLep", ";C_{em,W,-} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,20,0,300,20,0,3.2);
   hs2d.addHist("baseline/mumu/2d_CemVSdPhiMetNearLep", ";C_{em,W,-} (GeV);|#Delta#phi|(p_{T,DNN}^{miss},nearest l);EventsBIN" ,20,0,300,20,0,3.2);
   
   //Define 2D histograms for GOF studies
   std::vector<dim2D> dim2D_vec;
   hist::Histograms<TH2F> hs2d_GOF(vsDatasubsets);
   dim2D_vec.push_back({"PuppiMET_xy*cos(PuppiMET_xy_phi)","Puppi p_{x}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"PuppiMET_xy*sin(PuppiMET_xy_phi)","Puppi p_{y}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"MET_xy*cos(MET_xy_phi)","PF p_{x}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"MET_xy*sin(MET_xy_phi)","PF p_{y}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"vecsum_pT_allJet*cos(HT_phi)","p_{x}^{sum all jets} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"vecsum_pT_allJet*sin(HT_phi)","p_{y}^{sum all jets} (GeV) ",-150, 150, 6});
   dim2D_vec.push_back({"mass_l1l2_allJet", "mass_l1l2_allJet (GeV) ",0, 1800, 6});
   dim2D_vec.push_back({"Jet1_pt*sin(Jet1_phi)","p_{y}^{jet1} (GeV) ",-400, 400, 5});
   dim2D_vec.push_back({"MHT", "M_{HT} (GeV) ",0, 1800, 6});
   dim2D_vec.push_back({"Lep1_pt*cos(Lep1_phi)","p_{x}^{lep1} (GeV) ",-250, 250, 5});
   dim2D_vec.push_back({"Lep1_pt*sin(Lep1_phi)","p_{y}^{lep1} (GeV) ",-250, 250, 5});
   dim2D_vec.push_back({"Jet1_pt*cos(Jet1_phi)","p_{x}^{jet1} (GeV) ",-400, 400, 5});
   dim2D_vec.push_back({"CaloMET","CaloMET (GeV) ",0, 200, 5});
   dim2D_vec.push_back({"MT2","MT_{2} (GeV) ",0, 160, 5});
   dim2D_vec.push_back({"mjj", "m_{jj} (GeV) ",0, 1800, 6});
   dim2D_vec.push_back({"nJets","n_{jets}",1.5, 7.5, 6});
   dim2D_vec.push_back({"Jet1_E","E^{jet1} (GeV) ",0, 500, 5});
   dim2D_vec.push_back({"HT","H_{T} (GeV) ",0, 1200, 6});
   dim2D_vec.push_back({"Jet2_pt*cos(Jet2_phi)","p_{x}^{jet2} (GeV) ",-250, 250, 5});
   dim2D_vec.push_back({"Jet2_pt*sin(Jet2_phi)","p_{y}^{jet2} (GeV) ",-250, 250, 5});
   dim2D_vec.push_back({"DeepMET_reso*cos(DeepMET_reso_phi)","Deep(reso) p_{x}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"DeepMET_reso*sin(DeepMET_reso_phi)","Deep(reso) p_{y}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"DeepMET_resp*cos(DeepMET_resp_phi)","Deep(resp) p_{x}^{miss} (GeV)",-150, 150, 6});
   dim2D_vec.push_back({"DeepMET_resp*sin(DeepMET_resp_phi)","Deep(resp) p_{y}^{miss} (GeV)",-150, 150, 6});
   
   for (int i=0; i<dim2D_vec.size(); i++){
      for (int j=(i+1); j<dim2D_vec.size(); j++){
         for(TString channel:{"/ee","/mumu","/emu"}){
            TString path = "baseline_GOF2D"+channel+"/"+dim2D_vec[i].name+"_VS_"+dim2D_vec[j].name;
            TString title = ";"+dim2D_vec[i].axisTitle+";"+dim2D_vec[j].axisTitle+";EventsBIN";
            hs2d_GOF.addHist(path,title,dim2D_vec[i].nBins,dim2D_vec[i].xMin,dim2D_vec[i].xMax,dim2D_vec[j].nBins,dim2D_vec[j].xMin,dim2D_vec[j].xMax);
         }
      }
   }
   

   //Additional map to calculate signal efficiencies
   std::map<TString,float> count;
   std::map<TString,float> Ngen;
   
   //Define variable used for event runtime measurements
   double totalTime;
   double timePerEvent;
   
   for (TString ds_name: cfg.datasets.getDatasetNames()){
      auto ds=cfg.datasets.getDataset(ds_name);
      
      bool const isData=ds.isData;
      //Break if systematic shift is selected for data
      if(!isNominal && isData){
         std::cerr<<"Systematic shift should not be applied to data!"<<std::endl;
         exit(98);
      }
      
      // Check if current systematic and sample match
      Systematic::checkAlternativeSample(currentSystematic,ds.systName,ds.name);
      
      // Define Saver for minTrees
      io::RootFileSaver ttbar_res_saver(TString::Format("../minTrees/%.1f/%s%s.root",cfg.processFraction*100,("/"+currentSystematic.name()+"/").Data(),TString(ds_name+((cfg.fileNR==0)?"":"_"+std::to_string(cfg.fileNR))).Data()),TString::Format("ttbar_res%.1f",cfg.processFraction*100),true,true,true);
      TTree ttbar_res("ttbar_res","ttbar_res");
      
      int runEra=cfg.fileNR;     //int to store run Era in minimal Trees
                  
      for (auto dss: cfg.datasets.getDatasubsets({ds.name})){  // dataset loop (only used if multiple datasets are collected, usually not for condor jobs)
         TFile* file = TFile::Open(dss.getPath(),"read");
         if (file->IsZombie()) {
            return;
         }

         bool const isSignal=dss.isSignal;
         int year_int=cfg.year_int;
         
         if(cfg.fileNR==0)runEra++;
         if(!isData) runEra=0;
         
         // Configure JES/JER Corrections
         jesCorrections jesCorrector = jesCorrections(cfg.getJESPath(runEra,false).Data(),cfg.jes_UNC_mc_regrouped,currentSystematic,cfg.year);
         jesCorrections jesCorrector_puppi = jesCorrections(cfg.getJESPath(runEra,true).Data(),cfg.jes_UNC_mc_puppi_regrouped,currentSystematic,cfg.year);
         jerCorrections jerCorrector = jerCorrections(isData? cfg.jer_SF_data.Data() : cfg.jer_SF_mc.Data(),isData? cfg.jer_RES_data.Data() : cfg.jer_RES_mc.Data(),currentSystematic);
         
         // Configure pileupID application (applied in nominal selection)
         // ~bool applyPileupID = (currentSystematic.type() == Systematic::jetPileupIDapplied);
         bool applyPileupID = true;
         JetPileupIDWeights jetPileupIDWeighter = JetPileupIDWeights(cfg.jetPileupID_file.Data(),cfg.jetPileupID_sfHist.Data(),cfg.jetPileupID_effHist.Data(),currentSystematic);
         if(applyPileupID) std::cout<<"!!!!!!!!!!!!!!!JetPileupID is applied!!!!!!!!!!!!!!!!!!!"<<std::endl;
         
         //Configure jet veto maps (applied for UL18 only using HEM1516 region)
         // ~bool applyJetVetoMaps = (currentSystematic.type() == Systematic::applyJetVetoMaps || currentSystematic.type() == Systematic::applyJetVetoMaps_loose || currentSystematic.type() == Systematic::applyJetVetoMaps_HEM1516);
         bool applyJetVetoMaps_leading = (currentSystematic.type() == Systematic::applyJetVetoMaps_leading);
         bool applyJetVetoMaps_subleading = (currentSystematic.type() == Systematic::applyJetVetoMaps_subleading);
         bool applyJetVetoMaps_cleanedJets = (currentSystematic.type() == Systematic::applyJetVetoMaps_cleanedJets);
         bool applyJetVetoMaps = (year_int==3);    // apply for UL18 only
         JetVetoMaps jetVetoMaps = JetVetoMaps(cfg.jetVetoMap_file.Data(),(currentSystematic.type() == Systematic::applyJetVetoMaps_loose)? cfg.jetVetoMap_vetoMapHist_loose.Data():cfg.jetVetoMap_vetoMapHist.Data(),cfg.jetVetoMap_vetoMapHist_MC16.is_initialized()? cfg.jetVetoMap_vetoMapHist_MC16.get() : "",currentSystematic,year_int);
         if(applyJetVetoMaps || applyJetVetoMaps_leading || applyJetVetoMaps_subleading || applyJetVetoMaps_cleanedJets) std::cout<<"!!!!!!!!!!!!!!!JetVetoMaps are applied!!!!!!!!!!!!!!!!!!!"<<std::endl;
         
         // Configure lepton Correction
         leptonCorrections leptonCorretor = leptonCorrections(currentSystematic);
         
         // Configure bTag Weights
         BTagWeights bTagWeighter = BTagWeights(cfg.bTagSF_file.Data(),cfg.bTagEffPath.Data(),cfg.bTagger.Data(),BTagEntry::OperatingPoint(cfg.bTagWP),cfg.bTagWPcut,currentSystematic);
         std::vector<BTagWeights> bTagWeighter_vec;   //define vector for all bTagWeights related syst. (for nominal minTree)
         for (const Systematic::Type type : Systematic::weightTypes){
            if (std::find(Systematic::btagTypes.begin(), Systematic::btagTypes.end(), type) != Systematic::btagTypes.end()){
               std::cout.setstate(std::ios_base::failbit);     //avoid spaming cout
               bTagWeighter_vec.emplace_back(cfg.bTagSF_file.Data(),cfg.bTagEffPath.Data(),cfg.bTagger.Data(),BTagEntry::OperatingPoint(cfg.bTagWP),cfg.bTagWPcut,Systematic::Systematic(Systematic::convertType(type)+"_UP"));
               bTagWeighter_vec.emplace_back(cfg.bTagSF_file.Data(),cfg.bTagEffPath.Data(),cfg.bTagger.Data(),BTagEntry::OperatingPoint(cfg.bTagWP),cfg.bTagWPcut,Systematic::Systematic(Systematic::convertType(type)+"_DOWN"));
               std::cout.clear();
            }
         }
         
         // Configure leptonSF
         LeptonScaleFactors leptonSF = LeptonScaleFactors(cfg.electronID_file,cfg.electronID_hist,cfg.electronRECO_file,cfg.electronRECO_hist,
                                                         cfg.muonID_file,cfg.muonID_hist,cfg.muonISO_file,cfg.muonISO_hist,currentSystematic);
         leptonSF.setDYExtrapolationUncFactors(cfg.muonDYunc,cfg.electronDYunc);
         std::vector<LeptonScaleFactors> leptonSF_vec;   //define vector for all leptonSF related syst. (for nominal minTree)
         for (const Systematic::Type type : Systematic::weightTypes){
            if (std::find(Systematic::leptonsfTypes.begin(), Systematic::leptonsfTypes.end(), type) != Systematic::leptonsfTypes.end()){
               std::cout.setstate(std::ios_base::failbit);     //avoid spaming cout
               leptonSF_vec.push_back(LeptonScaleFactors(cfg.electronID_file,cfg.electronID_hist,cfg.electronRECO_file,cfg.electronRECO_hist,
                                                         cfg.muonID_file,cfg.muonID_hist,cfg.muonISO_file,cfg.muonISO_hist,Systematic::Systematic(Systematic::convertType(type)+"_UP")));
               leptonSF_vec.back().setDYExtrapolationUncFactors(cfg.muonDYunc,cfg.electronDYunc);
               leptonSF_vec.push_back(LeptonScaleFactors(cfg.electronID_file,cfg.electronID_hist,cfg.electronRECO_file,cfg.electronRECO_hist,
                                                         cfg.muonID_file,cfg.muonID_hist,cfg.muonISO_file,cfg.muonISO_hist,Systematic::Systematic(Systematic::convertType(type)+"_DOWN")));
               leptonSF_vec.back().setDYExtrapolationUncFactors(cfg.muonDYunc,cfg.electronDYunc);
               std::cout.clear();
            }
         }
         
         // Configure triggerSF
         TriggerScaleFactors triggerSFcalc = TriggerScaleFactors(cfg.trigger_SF.Data(),currentSystematic);
         std::vector<TriggerScaleFactors> triggerSF_vec;   //define vector for all trigger related syst. (for nominal minTree)
         std::cout.setstate(std::ios_base::failbit);     //avoid spaming cout
         triggerSF_vec.push_back(TriggerScaleFactors(cfg.trigger_SF.Data(),Systematic::Systematic("TRIG_UP")));
         triggerSF_vec.push_back(TriggerScaleFactors(cfg.trigger_SF.Data(),Systematic::Systematic("TRIG_DOWN")));
         std::cout.clear();
         
         //Configure unclustered energy shift
         TString metAddition = Systematic::metNameAddition(currentSystematic);
         
         //Configure mcWeights
         mcWeights mcWeighter = mcWeights(currentSystematic,dss.isMadGraph,dss.isPythiaOnly);
         std::vector<mcWeights> mcWeighter_vec;   //define vector for all leptonSF related syst. (for nominal minTree)
         for (const Systematic::Type type : Systematic::weightTypes){
            if (std::find(Systematic::mcWeightTypes.begin(), Systematic::mcWeightTypes.end(), type) != Systematic::mcWeightTypes.end()){
               std::cout.setstate(std::ios_base::failbit);     //avoid spaming cout
               if (type == Systematic::pdf){       // include all PDF variations
                  for (int i=1; i<=50; i++){
                     mcWeighter_vec.push_back(mcWeights(Systematic::Systematic("PDF_"+std::to_string(i)+"_UP"),dss.isMadGraph,dss.isPythiaOnly));
                     mcWeighter_vec.back().prepareLumiWeight(dss,cfg.lumi);
                     mcWeighter_vec.push_back(mcWeights(Systematic::Systematic("PDF_"+std::to_string(i)+"_DOWN"),dss.isMadGraph,dss.isPythiaOnly));
                     mcWeighter_vec.back().prepareLumiWeight(dss,cfg.lumi);
                  }
               }
               else{
                  mcWeighter_vec.push_back(mcWeights(Systematic::Systematic(Systematic::convertType(type)+"_UP"),dss.isMadGraph,dss.isPythiaOnly));
                  mcWeighter_vec.back().prepareLumiWeight(dss,cfg.lumi);
                  mcWeighter_vec.push_back(mcWeights(Systematic::Systematic(Systematic::convertType(type)+"_DOWN"),dss.isMadGraph,dss.isPythiaOnly));
                  mcWeighter_vec.back().prepareLumiWeight(dss,cfg.lumi);
               }
               std::cout.clear();
            }
         }
         
         //Set current sample for histograms
         hs.setCurrentSample(dss.name);
         hs_cutflow.setCurrentSample(dss.name);
         hs2d.setCurrentSample(dss.name);
         hs2d_GOF.setCurrentSample(dss.name);
         
         //Lumi weight for current sample
         float lumi_weight=dss.xsec/dss.getNgen_syst(currentSystematic)*cfg.lumi;
         float lumi_weight_puUp=dss.xsec/dss.getNgen_syst(Systematic::Systematic("PU_UP"))*cfg.lumi;
         float lumi_weight_puDown=dss.xsec/dss.getNgen_syst(Systematic::Systematic("PU_DOWN"))*cfg.lumi;
         float lumi_weight_topPT=dss.xsec/dss.getNgen_syst(Systematic::Systematic("TOP_PT"))*cfg.lumi;
         float lumi_unc = cfg.systUncFactor.at("LUMI").first;
                  
         //Save number of gen events for efficiency
         Ngen[dss.datasetName]=dss.Ngen_woWeight;
         
         //Check if current sample is TTbar 2L sample (later used to veto tau events)
         bool ttBar_dilepton=dss.isTTbar2L;
                  
         //Check if current sample is TTbar powheg dilepton tau
         bool ttBar_dilepton_tau=false;
         if (dss.datasetName=="TTbar_diLepton_tau") ttBar_dilepton_tau=true;
         
         //Check if current sample is TTbar powheg semileptonic
         bool ttBar_singleLepton=false;
         if (dss.datasetName=="TTbar_singleLepton") ttBar_singleLepton=true;
         
         //Check if current sample is TTbar powheg hadronic
         bool ttBar_hadronic=false;
         if (dss.datasetName=="TTbar_hadronic") ttBar_hadronic=true;
         
         //Check if current sample is TTbar amc@NLO
         bool ttBar_amc=false;
         if (dss.datasetName=="TTbar_amcatnlo") ttBar_amc=true;
         
         //Check if current sample is TTbar herwig
         bool ttBar_herwig=false;
         if (dss.datasetName=="TTbar_herwig") ttBar_herwig=true;
         
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
         bool EGamma=false;
         if (dss.datasetName.find("SingleElectron")!=std::string::npos) SingleElectron=true;
         else if (dss.datasetName.find("SingleMuon")!=std::string::npos) SingleMuon=true;
         else if (dss.datasetName.find("DoubleEG")!=std::string::npos) DoubleEG=true;
         else if (dss.datasetName.find("DoubleMuon")!=std::string::npos) DoubleMuon=true;
         else if (dss.datasetName.find("MuonEG")!=std::string::npos) MuonEG=true;
         else if (dss.datasetName.find("EGamma")!=std::string::npos) EGamma=true;
         
         //Check if current sample is Run2016H
         bool Run2016H=false;
         if (dss.name.find("Run2016H")!=std::string::npos) Run2016H=true;
         
         //Check if current sample is Run2017AB
         bool Run2017AB=false;
         if (dss.name.find("Run2017A")!=std::string::npos) Run2017AB=true;
         else if (dss.name.find("Run2017B")!=std::string::npos) Run2017AB=true;
         
         //Check if gen level deltaR cut should be applied (should be always applied in nominal)
         // ~bool applyDeltaRgen = (currentSystematic.type() == Systematic::applyGenLevel_DeltaRcut);
         bool applyDeltaRgen = true;
         
         //Configure topPT reweighting
         bool applytopPTreweighting = (checkTopPTreweighting(currentSystematic) && (ttBar_dilepton || ttBar_dilepton_tau || ttBar_amc || ttBar_herwig));
         // ~bool applytopPTreweighting = checkTopPTreweighting(currentSystematic);
         
         //variables used for storing in minimal trees
         float minTree_MET, minTree_PtNuNu, minTree_PhiRec, minTree_PhiGen, minTree_PhiNuNu, minTree_PhiMetNearJet, minTree_PhiMetFarJet, minTree_PhiMetLeadJet, minTree_PhiMetLead2Jet,
         minTree_PhiMetbJet, minTree_dPhiLep1Lep2, minTree_dPhiJet1Jet2, minTree_METsig, minTree_N, minTree_SF, minTree_totalWeight, minTree_genMet, minTree_PuppiMet,minTree_PuppiMet_xy,minTree_Met_xy,minTree_DeepMet_reso,minTree_DeepMet_resp,
         minTree_HT, minTree_HT_phi, minTree_MHT, minTree_MT, minTree_genMT, minTree_MT_nextLep, minTree_genMT_nextLep,
         minTree_PhiPtnunuMet, minTree_leadTop, minTree_dPhiNuNu, minTree_PhiRecPuppi,minTree_PhiRecPuppi_xy,minTree_PhiRec_xy, minTree_PhiMetNearJet_Puppi, minTree_PhiMetFarJet_Puppi,
         minTree_PhiMetLeadJet_Puppi, minTree_PhiMetLead2Jet_Puppi, minTree_PhiMetbJet_Puppi, minTree_dPhiLep1bJet, minTree_dPhiLep1Jet1, minTree_ratioMET_sqrtMETunc_Puppi,
         minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet, minTree_METunc_Puppi, minTree_METunc_PF, minTree_absmetres_PUPPI,
         minTree_Lep1_pt, minTree_Lep1_phi, minTree_Lep1_eta, minTree_Lep1_E, minTree_Lep1_flavor,
         minTree_Lep2_pt, minTree_Lep2_phi, minTree_Lep2_eta, minTree_Lep2_E, minTree_Lep2_flavor,
         minTree_Jet1_pt, minTree_Jet1_phi, minTree_Jet1_eta, minTree_Jet1_E, minTree_Jet1_bTagScore, minTree_Jet1_unc,
         minTree_Jet2_pt, minTree_Jet2_phi, minTree_Jet2_eta, minTree_Jet2_E, minTree_Jet2_bTagScore, minTree_Jet2_unc,
         minTree_PFMET_phi, minTree_PuppiMET_phi, minTree_PuppiMET_xy_phi, minTree_MET_xy_phi, minTree_CaloMET, minTree_CaloMET_phi, minTree_PFMETxy, minTree_PFMETxy_phi,minTree_DeepMet_reso_phi, minTree_DeepMet_resp_phi, minTree_genMET_phi, minTree_PtNuNu_phi, minTree_nJets, minTree_nGenJets, minTree_n_Interactions,minTree_DNN_MET_pT, minTree_DNN_MET_phi,minTree_DNN_MET_dPhi_nextLep,
         minTree_mLL, minTree_PtLL, minTree_PtLLgen, minTree_MT2, minTree_vecsum_pT_allJet, minTree_vecsum_pT_l1l2_allJet, minTree_mass_l1l2_allJet, minTree_ratio_vecsumpTlep_vecsumpTjet, minTree_mjj;
         UInt_t minTree_runNo, minTree_lumNo, minTree_genDecayMode, minTree_n_Interactions_gen, minTree_looseLeptonVeto, minTree_NpromptNeutrinos, minTree_NnonpromptNeutrinos, minTree_ee, minTree_mumu, minTree_emu;
         ULong64_t minTree_evtNo;
         std::vector<float> minTree_systWeights(Systematic::numberOfWeightTypes());
         
         //Define branches of mininmal tree
         ttbar_res.Branch("ee",&minTree_ee,"ee/i");
         ttbar_res.Branch("mumu",&minTree_mumu,"mumu/i");
         ttbar_res.Branch("emu",&minTree_emu,"emu/i");
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
         ttbar_res.Branch("totalWeight",&minTree_totalWeight,"totalWeight/f");
         ttbar_res.Branch("runNo",&minTree_runNo,"runNo/i");
         ttbar_res.Branch("lumNo",&minTree_lumNo,"lumNo/i");
         ttbar_res.Branch("evtNo",&minTree_evtNo,"evtNo/l");
         ttbar_res.Branch("runEra",&runEra,"runEra/i");
         ttbar_res.Branch("genDecayMode",&minTree_genDecayMode,"genDecayMode/i");
         ttbar_res.Branch("genMET",&minTree_genMet,"genMET/f");
         ttbar_res.Branch("PuppiMET",&minTree_PuppiMet,"PuppiMET/f");
         ttbar_res.Branch("PuppiMET_xy",&minTree_PuppiMet_xy,"PuppiMET_xy/f");
         ttbar_res.Branch("MET_xy",&minTree_Met_xy,"MET_xy/f");
         ttbar_res.Branch("DeepMET_reso",&minTree_DeepMet_reso,"DeepMET_reso/f");
         ttbar_res.Branch("DeepMET_resp",&minTree_DeepMet_resp,"DeepMET_resp/f");
         ttbar_res.Branch("HT",&minTree_HT,"HT/f");
         ttbar_res.Branch("HT_phi",&minTree_HT_phi,"HT_phi/f");
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
         ttbar_res.Branch("Phi_recPuppi_xy",&minTree_PhiRecPuppi_xy,"Phi_recPuppi_xy/f");
         ttbar_res.Branch("Phi_rec_xy",&minTree_PhiRec_xy,"Phi_rec_xy/f");
         ttbar_res.Branch("looseLeptonVeto",&minTree_looseLeptonVeto,"looseLeptonVeto/i");
         ttbar_res.Branch("nJets",&minTree_nJets,"nJets/f");
         ttbar_res.Branch("nGenJets",&minTree_nGenJets,"nGenJets/f");
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
         ttbar_res.Branch("mLL",&minTree_mLL,"mLL/f");
         ttbar_res.Branch("pTll",&minTree_PtLL,"pTll/f");
         ttbar_res.Branch("pTllGen",&minTree_PtLLgen,"pTllGen/f");
         ttbar_res.Branch("PFMET_phi",&minTree_PFMET_phi,"PFMET_phi/f");
         ttbar_res.Branch("PuppiMET_phi",&minTree_PuppiMET_phi,"PuppiMET_phi/f");
         ttbar_res.Branch("PuppiMET_xy_phi",&minTree_PuppiMET_xy_phi,"PuppiMET_xy_phi/f");
         ttbar_res.Branch("MET_xy_phi",&minTree_MET_xy_phi,"MET_xy_phi/f");
         ttbar_res.Branch("DeepMET_reso_phi",&minTree_DeepMet_reso_phi,"DeepMET_reso_phi/f");
         ttbar_res.Branch("DeepMET_resp_phi",&minTree_DeepMet_resp_phi,"DeepMET_resp_phi/f");
         ttbar_res.Branch("CaloMET",&minTree_CaloMET,"CaloMET/f");
         ttbar_res.Branch("CaloMET_phi",&minTree_CaloMET_phi,"CaloMET_phi/f");
         ttbar_res.Branch("genMET_phi",&minTree_genMET_phi,"genMET_phi/f");
         ttbar_res.Branch("PtNuNu_phi",&minTree_PtNuNu_phi,"PtNuNu_phi/f");
         ttbar_res.Branch("NpromptNeutrinos",&minTree_NpromptNeutrinos,"NpromptNeutrinos/i");
         ttbar_res.Branch("NnonpromptNeutrinos",&minTree_NnonpromptNeutrinos,"NnonpromptNeutrinos/i");
         ttbar_res.Branch("DNN_MET_pT",&minTree_DNN_MET_pT,"DNN_MET_pT/f");
         ttbar_res.Branch("DNN_MET_phi",&minTree_DNN_MET_phi,"DNN_MET_phi/f");
         ttbar_res.Branch("DNN_MET_dPhi_nextLep",&minTree_DNN_MET_dPhi_nextLep,"DNN_MET_dPhi_nextLep/f");
         ttbar_res.Branch("MT2",&minTree_MT2,"MT2/f");
         ttbar_res.Branch("vecsum_pT_allJet",&minTree_vecsum_pT_allJet,"vecsum_pT_allJet/f");
         ttbar_res.Branch("vecsum_pT_l1l2_allJet",&minTree_vecsum_pT_l1l2_allJet,"vecsum_pT_l1l2_allJet/f");
         ttbar_res.Branch("mass_l1l2_allJet",&minTree_mass_l1l2_allJet,"mass_l1l2_allJet/f");
         ttbar_res.Branch("ratio_vecsumpTlep_vecsumpTjet",&minTree_ratio_vecsumpTlep_vecsumpTjet,"ratio_vecsumpTlep_vecsumpTjet/f");
         ttbar_res.Branch("mjj",&minTree_mjj,"mjj/f");
         //add branches for syst. weight only to nominal min Tree
         if (isNominal){
            int weightCounter = 0;
            for (const Systematic::Type type : Systematic::weightTypes){
               TString systName = convertType(type);
               if (std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), type) != Systematic::upDownTypes.end()){
                  if (type == Systematic::pdf){    //Store all pdf shifts
                     for (int i=1; i<=50; i++){
                        ttbar_res.Branch("weight_"+systName+"_"+std::to_string(i)+"_UP",&(minTree_systWeights[weightCounter]),"weight_"+systName+"_"+std::to_string(i)+"_UP/f");
                        weightCounter++;
                        ttbar_res.Branch("weight_"+systName+"_"+std::to_string(i)+"_DOWN",&(minTree_systWeights[weightCounter]),"weight_"+systName+"_"+std::to_string(i)+"_DOWN/f");
                        weightCounter++;
                     }
                  }
                  else{
                     ttbar_res.Branch("weight_"+systName+"_UP",&(minTree_systWeights[weightCounter]),"weight_"+systName+"_UP/f");
                     weightCounter++;
                     ttbar_res.Branch("weight_"+systName+"_DOWN",&(minTree_systWeights[weightCounter]),"weight_"+systName+"_DOWN/f");
                     weightCounter++;
                  }
               }
               else{
                  ttbar_res.Branch("weight_"+systName,&(minTree_systWeights[weightCounter]),"weight_"+systName+"/f");
                  weightCounter++;
               }
            }
         } 
         
         //Initialize DNN regression
         TString dnnPath = cfg.DNN_Path;
         if(currentSystematic.type() == Systematic::useDNNnoMETcut){
            dnnPath = "../data/DNN/Inlusive_noMetCut_amcatnlo_xyComponent_JetLepXY_50EP__diff_xy_2018_20230111-1004genMETweighted";
         }
         else if (currentSystematic.type() == Systematic::useDNNmumu){
            dnnPath = "../data/DNN/Inlusive_amcatnlo_xyComponent_JetLepXY_50EP_mumuOnly__diff_xy_2018_20230112-1124genMETweighted";
         }
         else if (currentSystematic.type() == Systematic::useDNNnoMetCutDY){
            dnnPath = "../data/DNN/Inlusive_noMetCut_amcatnlo_DrellYan_xyComponent_JetLepXY_50EP__diff_xy_2018_20230113-1059genMETweighted";
         }
         DNNregression dnnRegression(dnnPath.Data());
         std::vector<float> input_vec(17);
         std::vector<float> output_vec(2);
         
         
         //Set Tree Input variables
         TTreeReader reader(cfg.treeName, file);
         TTreeReaderValue<float> w_pu(reader, Systematic::puWeightName(currentSystematic));
         TTreeReaderValue<float> w_pu_up(reader, "pu_weight_up" );
         TTreeReaderValue<float> w_pu_down(reader, "pu_weight_down" );
         TTreeReaderValue<UInt_t> runNo(reader, "runNo");
         TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
         TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
         TTreeReaderValue<float> w_mc(reader, "mc_weight");
         TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
         TTreeReaderValue<std::vector<float>> w_ps(reader, "ps_weights");
         // ~TTreeReaderValue<float> w_topPT(reader, "topPTweight");
         TTreeReaderValue<float> w_topPT(reader, "topPTweightNNLO");
         TTreeReaderValue<float> fragUpWeight(reader, "weightFragUp");
         TTreeReaderValue<float> fragCentralWeight(reader, "weightFragCentral");
         TTreeReaderValue<float> fragDownWeight(reader, "weightFragDown");
         TTreeReaderValue<float> fragPetersonWeight(reader, "weightFragPeterson");
         TTreeReaderValue<float> semilepbrUpWeight(reader, "weightSemilepbrUp");
         TTreeReaderValue<float> semilepbrDownWeight(reader, "weightSemilepbrDown");
         TTreeReaderValue<double> w_prefiring(reader, Systematic::prefiringWeightName(currentSystematic));
         TTreeReaderValue<double> w_prefiring_up(reader, "prefiring_weight_up");
         TTreeReaderValue<double> w_prefiring_down(reader, "prefiring_weight_down");
         TTreeReaderValue<float> matchWeight_up(reader, Systematic::matchDCTRWeightName(ds.name,true));
         TTreeReaderValue<float> matchWeight_down(reader, Systematic::matchDCTRWeightName(ds.name,false));
         TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
         TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
         TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
         TTreeReaderValue<std::vector<tree::Jet>>      jets_puppi     (reader, "jetsPuppi");
         TTreeReaderValue<std::vector<tree::Particle>>      genJets     (reader, "genJets");
         TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
         TTreeReaderValue<tree::MET> MET(reader, "met"+metAddition);
         TTreeReaderValue<tree::MET> MET_xy(reader, "metXYcorr"+metAddition);
         TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
         TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi"+metAddition);
         TTreeReaderValue<tree::MET> MET_Puppi_xy(reader, "metPuppiXYcorr"+metAddition);
         TTreeReaderValue<tree::MET> DeepMET_reso(reader, "metDeepResolution"+metAddition);
         TTreeReaderValue<tree::MET> DeepMET_resp(reader, "metDeepResponse"+metAddition);
         TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
         TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
         TTreeReaderValue<int> n_Interactions(reader, "nPV");
         TTreeReaderValue<int> n_Interactions_gen(reader, "true_nPV");
         TTreeReaderValue<float> HTgen(reader, "genHt");
         TTreeReaderValue<float> rho(reader, "rho");
         TTreeReaderValue<bool> addLepton   (reader, "addLepton");
         TTreeReaderValue<float> genMT2   (reader, "genMT2");
         TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
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
      
         // Log current sample
         io::log * ("Processing '"+dss.name+"' ");
         
         int iEv=0;
         int processEvents=cfg.processFraction*dss.entries;
         auto startTime_evtLoop = high_resolution_clock::now();
         while (reader.Next()){     // Event loop
            iEv++;
            if (iEv>processEvents) break;
            if (iEv%(std::max(processEvents/10,1))==0){
               io::log*".";
               io::log.flush();
            }
                                                
            //Booleans for reco and pseudo selection
            bool rec_selection=false;
            bool pseudo_selection=true;
            
            //Do not use tau events in signal sample
            if (ttBar_dilepton && *genDecayMode>3) continue;
            
            //Do only use ee,emu,mumu in amc ttbar or herwig ttbar
            if ((ttBar_amc || ttBar_herwig) && (*genDecayMode>3 || *genDecayMode==0)) continue;
            
            // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
            TLorentzVector neutrinoPair(0,0,0,0);
            neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
            
            if (*genDecayMode_pseudo==0 || (*genDecayMode_pseudo!=3 && neutrinoPair.Pt()<40)) pseudo_selection=false; //pseudo baseline selection
            
            // Apply gen level delta Rcut between bJets and leptons
            if (applyDeltaRgen && pseudo_selection){
               if(genB->DeltaR(*genLepton)<0.4 || genB->DeltaR(*genAntiLepton)<0.4 || genAntiB->DeltaR(*genLepton)<0.4 || genAntiB->DeltaR(*genAntiLepton)<0.4) pseudo_selection=false;
            }
            
            // Construct vector of different METs for correction
            std::vector<tree::MET*> PFMETs = {&(*MET),&(*MET_xy),&(*MET_Calo),&(*DeepMET_reso),&(*DeepMET_resp)};
            std::vector<tree::MET*> PuppiMETs = {&(*MET_Puppi),&(*MET_Puppi_xy)};
            std::vector<tree::MET*> allMETs = {&(*MET_Puppi),&(*MET_Puppi_xy),&(*MET),&(*MET_xy),&(*MET_Calo),&(*DeepMET_reso),&(*DeepMET_resp)};
            
            // Correct and select leptons
            *muons = leptonCorretor.correctMuons(*muons,PFMETs,PuppiMETs);
            *electrons = leptonCorretor.correctElectrons(*electrons,PFMETs,PuppiMETs);
            
            //Baseline selection (including separation into ee, emu, mumu)
            TLorentzVector p_l1;
            TLorentzVector p_l2;
            int flavor_l1=0;  //1 for electron and 2 for muon
            int flavor_l2=0;
            double etaSC_l1=0;  // etaSC needed for electron ID SF
            double etaSC_l2=0;
            bool muonLead=true; //Boolean for emu channel
            std::vector<bool> channel={false,false,false};     //ee, mumu, emu
            TString cat="";
            
            rec_selection=selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,etaSC_l1,etaSC_l2,cat,muonLead);
            if(rec_selection==false && pseudo_selection==false) continue;
                        
            //Trigger selection
            std::vector<bool> diElectronTriggers={*eleTrigg1,*eleTrigg2,*singleEleTrigg};
            std::vector<bool> diMuonTriggers={*muonTrigg1,*muonTrigg2,*muonTrigg3,*muonTrigg4,*singleMuonTrigg1,*singleMuonTrigg2};
            std::vector<bool> electronMuonTriggers={*eleMuTrigg1,*eleMuTrigg2,*eleMuTrigg3,*eleMuTrigg4,*singleMuonTrigg1,*singleMuonTrigg2,*singleEleTrigg};
            std::vector<bool> PD={SingleElectron,DoubleEG,SingleMuon,DoubleMuon,MuonEG,EGamma};
            bool triggerMC=true;
            
            if (!isData){
               triggerMC=selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,false,year_int);
            }
            else{
               if(!selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,true,year_int,PD,Run2016H,Run2017AB)) continue;
            }
                        
            if (triggerMC==false) rec_selection=false;
            if(rec_selection==false && pseudo_selection==false) continue;
            
            
            //Apply JES
            if(!isData){
               jesCorrector.applySystematics(*jets,PFMETs);
               jesCorrector_puppi.applySystematics(*jets_puppi,PuppiMETs);    // Needed for correction of Puppi MET
            }
            
            //Apply jeto veto map if selected (for all jets)
            if(applyJetVetoMaps) {
               if(jetVetoMaps.checkVetoMap(*jets) == false) rec_selection=false;
            } 
                        
            // Get leptonSF weight
            float leptonSFweight =(isData)? 1. : leptonSF.getSFDilepton(p_l1,p_l2,flavor_l1,flavor_l2,etaSC_l1,etaSC_l2);
            
            // Get mcWeight
            std::vector<float> bFrag_vec = {*fragUpWeight,*fragCentralWeight,*fragDownWeight,*fragPetersonWeight,*semilepbrUpWeight,*semilepbrDownWeight};
            std::vector<float> match_vec = {*matchWeight_up,*matchWeight_down};
            float mcWeight = (isData)? 1. : mcWeighter.getMCweight(*w_mc,*w_pdf,*w_ps,bFrag_vec,match_vec);
                        
            if(isSignal || !applytopPTreweighting) *w_topPT=1.;  //Set top weight for signals or if topPTsystematic to 1 
            float cutFlow_weight=(isData)? 1: *w_pu * mcWeight * *w_topPT * leptonSFweight * *w_prefiring;
            
            if(rec_selection) hs_cutflow.fillweight("cutflow/"+cat,1,cutFlow_weight);
            
            // Get cleaned and smeared jets
            std::vector<tree::Jet> cjets = phys::getCleanedJets(*jets, p_l1, p_l2,jerCorrector,*rho,allMETs,applyPileupID);
            
            //Apply jeto veto map if selected (for leading or leading+subleading jets)
            if(applyJetVetoMaps_subleading) {
               if(cjets.size()>1){
                  if(jetVetoMaps.checkVetoMap({cjets.begin(), cjets.begin()+2}) == false) rec_selection=false;
               }
            }
            else if(applyJetVetoMaps_leading){
               if(cjets.size()>1){
                  if(jetVetoMaps.checkVetoMap({cjets.begin(),cjets.begin()+1}) == false) rec_selection=false;
               }  
            }
            else if(applyJetVetoMaps_cleanedJets){
               if(cjets.size()>1){
                  if(jetVetoMaps.checkVetoMap(cjets) == false) rec_selection=false;
               }  
            }
            
            //Set some met variables
            float met=MET->p.Pt();
            float const met_puppi=MET_Puppi->p.Pt();
            float const met_puppi_xy=MET_Puppi_xy->p.Pt();
            float const genMet=GENMET->p.Pt();
            float mt2 = phys::MT2(p_l1,p_l2,MET_Puppi->p);
            
            // Apply ttBar selection
            std::vector<tree::Jet> BJets;
            std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi_xy,channel,*jets,cjets,BJets);
            
            if(rec_selection && ttbarSelection[0]){
               hs_cutflow.fillweight("cutflow/"+cat,2,cutFlow_weight);
               if(ttbarSelection[1]){
                  hs_cutflow.fillweight("cutflow/"+cat,3,cutFlow_weight);
               }
            }
            
            //Check booleans from ttBar selection (Apply relaxed selection with only two lepton) 
            if(currentSystematic.type()==Systematic::removeMLLcut){
               if(!std::all_of(ttbarSelection.begin()+1, ttbarSelection.end(), [](bool v) { return v; })) rec_selection=false;
            }
            else{ // no MET cut on Puppi because MET cut is performed for DNN MET further below
               if(!ttbarSelection[0] || !ttbarSelection[1] || !ttbarSelection[3]) rec_selection=false;
            }
            
            // Ignore events which do not fulfill reco or gen selection
            if(rec_selection==false && pseudo_selection==false) continue;
            
            // Cross check for add. MET cut in emu
            if(currentSystematic.type()==Systematic::met40Cut && met_puppi<40) rec_selection=false;
            if(rec_selection==false && pseudo_selection==false) continue;
            
            //Get bTag weight
            float bTagWeight = 1.;
            int channelID = 0;
            if(rec_selection) {
               if (channel[0]) channelID = 1;
               else if (channel[1]) channelID = 2;
               else if (channel[2]) channelID = 3;
               if (!isData) bTagWeight = bTagWeighter.getEventWeight(cjets,channelID);
            }
            
            //Get jet pileup ID weight
            float jetPileupIDWeight = 1.;
            if (!isData && rec_selection && applyPileupID) jetPileupIDWeight = jetPileupIDWeighter.getEventWeight(cjets);
 
            // end reco baseline selection
                        
            if(rec_selection==false && pseudo_selection==false) continue;  // only proceed with events selected by one of the baseline selection (reco or pseudo)
            
            // Continue if data events and not selected
            if(rec_selection==false && isData) continue;
            
            
            //Set Event weights
            float fEventWeight = 1.;
            float SFWeight = 1.;
            float triggerSF = 1.;
            if(!isData){
               if(rec_selection) triggerSF = triggerSFcalc.getTriggerSF(p_l1.Pt(),p_l2.Pt(),channelID,muonLead);
               fEventWeight=*w_pu * mcWeight * *w_topPT;     //Set event weight 
               SFWeight=leptonSFweight * bTagWeight * jetPileupIDWeight * triggerSF * *w_prefiring;     //Set combined SF weight
               hs.setFillWeight(fEventWeight*SFWeight);
               hs2d.setFillWeight(fEventWeight*SFWeight);
               hs2d_GOF.setFillWeight(fEventWeight*SFWeight);
            }
            else{
               hs.setFillWeight(1);
               hs2d.setFillWeight(1);
               hs2d_GOF.setFillWeight(1);
               hs_cutflow.setFillWeight(1);
            }
            if(rec_selection){
               if(isData) hs_cutflow.fillweight("cutflow/"+cat,4,1);
               else hs_cutflow.fillweight("cutflow/"+cat,4,*w_pu * mcWeight * leptonSFweight * *w_prefiring * bTagWeight * jetPileupIDWeight * *w_topPT *triggerSF);
            }
            
            //Fill syst. weights to minimal tree if nominal
            if(isNominal && !isData){
               if (rec_selection){
                  int currentPos = 0;
                  for (int i=0; i<leptonSF_vec.size(); i++){
                     minTree_systWeights[i] = lumi_weight * fEventWeight * bTagWeight * jetPileupIDWeight * triggerSF * *w_prefiring * leptonSF_vec[i].getSFDilepton(p_l1,p_l2,flavor_l1,flavor_l2,etaSC_l1,etaSC_l2);
                  }
                  currentPos += leptonSF_vec.size(); 
                  for (int i=0; i<bTagWeighter_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = lumi_weight * fEventWeight * triggerSF * leptonSFweight * *w_prefiring * jetPileupIDWeight * bTagWeighter_vec[i].getEventWeight(cjets,channelID);
                  }
                  currentPos += bTagWeighter_vec.size();
                  for (int i=0; i<triggerSF_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = lumi_weight * fEventWeight * bTagWeight * jetPileupIDWeight * leptonSFweight * *w_prefiring * triggerSF_vec[i].getTriggerSF(p_l1.Pt(),p_l2.Pt(),channelID,muonLead);
                  }
                  currentPos += triggerSF_vec.size();
                  for (int i=0; i<mcWeighter_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = SFWeight * *w_pu * *w_topPT * mcWeighter_vec[i].getMCweight_lumiWeighted(*w_mc,*w_pdf,*w_ps,bFrag_vec,match_vec);
                  }
                  // ~std::cout<<minTree_systWeights[currentPos]<<std::endl;
                  // ~exit(1);
                  currentPos += mcWeighter_vec.size();
                  minTree_systWeights[currentPos] = lumi_weight_puUp * SFWeight * *w_pu_up * *w_topPT * mcWeight;   //PU_UP
                  minTree_systWeights[currentPos+1] = lumi_weight_puDown * SFWeight * *w_pu_down * *w_topPT * mcWeight;     //PU_DOWN
                  minTree_systWeights[currentPos+2] = lumi_weight_topPT * SFWeight * *w_pu * mcWeight;   //TOP_PT
                  minTree_systWeights[currentPos+3] = (lumi_weight*(1.+lumi_unc)) * SFWeight * *w_pu * *w_topPT * mcWeight;   //LUMI_UP
                  minTree_systWeights[currentPos+4] = (lumi_weight*(1.-lumi_unc)) * SFWeight * *w_pu * *w_topPT * mcWeight;   //LUMI_DOWN
               }
               else{    // do not include SF weights for non-reco events
                  int currentPos = 0;
                  for (int i=0; i<leptonSF_vec.size(); i++){
                     minTree_systWeights[i] = lumi_weight * fEventWeight;
                  }
                  currentPos += leptonSF_vec.size();
                  for (int i=0; i<bTagWeighter_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = lumi_weight * fEventWeight;
                  }
                  currentPos += bTagWeighter_vec.size();
                  for (int i=0; i<triggerSF_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = lumi_weight * fEventWeight;
                  }
                  currentPos += triggerSF_vec.size();
                  for (int i=0; i<mcWeighter_vec.size(); i++){
                     minTree_systWeights[currentPos+i] = *w_pu * *w_topPT * mcWeighter_vec[i].getMCweight_lumiWeighted(*w_mc,*w_pdf,*w_ps,bFrag_vec,match_vec);
                  }
                  currentPos += mcWeighter_vec.size();
                  minTree_systWeights[currentPos] = lumi_weight_puUp * *w_pu_up * *w_topPT * mcWeight;   //PU_UP
                  minTree_systWeights[currentPos+1] = lumi_weight_puDown * *w_pu_down * *w_topPT * mcWeight;     //PU_DOWN
                  minTree_systWeights[currentPos+2] = lumi_weight_topPT * *w_pu * mcWeight;   //TOP_PT
                  minTree_systWeights[currentPos+3] = (lumi_weight*(1.+lumi_unc)) * *w_pu * *w_topPT * mcWeight;   //LUMI_UP
                  minTree_systWeights[currentPos+4] = (lumi_weight*(1.-lumi_unc)) * *w_pu * *w_topPT * mcWeight;   //LUMI_DOWN
               }
            }
            
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
            float dPhiMETnearLepPuppi_xy=4;
            float dPhiMETnearLep_xy=4;
            float mt_MetNextLep=0; 
            float mt_NuNuNextLep=0; 
            for (TLorentzVector const lep : {p_l1,p_l2}){
               const float dPhi=MET->p.DeltaPhi(lep);
               const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
               const float dPhi_Puppi_xy=MET_Puppi_xy->p.DeltaPhi(lep);
               const float dPhi_xy=MET_xy->p.DeltaPhi(lep);
               if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
                  dPhiMETnearLep=dPhi;
                  mt_MetNextLep=phys::M_T(MET->p,lep);
               }
               if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLepPuppi)) {
                  dPhiMETnearLepPuppi=dPhi_Puppi;
               }
               if (std::abs(dPhi_Puppi_xy) < std::abs(dPhiMETnearLepPuppi_xy)) {
                  dPhiMETnearLepPuppi_xy=dPhi_Puppi_xy;
               }
               if (std::abs(dPhi_xy) < std::abs(dPhiMETnearLep_xy)) {
                  dPhiMETnearLep_xy=dPhi_xy;
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
                  const float dPhi_Puppi=MET_Puppi_xy->p.DeltaPhi(jet.p);
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
               dPhiMETleadJet_Puppi=MET_Puppi_xy->p.DeltaPhi(cjets[0].p);
               dPhiMETlead2Jet_Puppi=MET_Puppi_xy->p.DeltaPhi(cjets[1].p);
               dPhiMetBJet_Puppi=MET_Puppi_xy->p.DeltaPhi(BJets[0].p);
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
            
            // Gen nGenJets with pT>30
            std::vector<tree::Particle> genJets30;
            for (const tree::Particle &j: *genJets){
               if (j.p.Pt()>30){
                  genJets30.push_back(j);
               }
            }
            
            //Fill minimal tree for TTbar resolution used in binning/unfolding studies
            minTree_ee=channel[0];
            minTree_mumu=channel[1];
            minTree_emu=channel[2];
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
            minTree_genDecayMode=*genDecayMode_pseudo;
            minTree_genMet=genMet;
            minTree_PuppiMet=met_puppi;
            minTree_PuppiMet_xy=MET_Puppi_xy->p.Pt();
            minTree_Met_xy=MET_xy->p.Pt();
            minTree_DeepMet_reso=DeepMET_reso->p.Pt();
            minTree_DeepMet_resp=DeepMET_resp->p.Pt();
            minTree_HT=HT;
            minTree_HT_phi=MHT.Phi();
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
            minTree_PhiRecPuppi_xy=abs(dPhiMETnearLepPuppi_xy);
            minTree_PhiRec_xy=abs(dPhiMETnearLep_xy);
            minTree_looseLeptonVeto= *addLepton? 1: 0;
            minTree_nJets=cjets.size();
            minTree_nGenJets=genJets30.size();
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
            minTree_PuppiMET_xy_phi=MET_Puppi_xy->p.Phi();
            minTree_MET_xy_phi=MET_xy->p.Phi();
            minTree_DeepMet_reso_phi=DeepMET_reso->p.Phi();
            minTree_DeepMet_resp_phi=DeepMET_resp->p.Phi();
            minTree_CaloMET=MET_Calo->p.Pt();
            minTree_CaloMET_phi=MET_Calo->p.Phi();
            minTree_genMET_phi=GENMET->p.Phi();
            minTree_PtNuNu_phi=neutrinoPair.Phi();
            minTree_NpromptNeutrinos=NpromptNeutrinos;
            minTree_NnonpromptNeutrinos=NnonpromptNeutrinos;
            minTree_MT2=mt2;
            minTree_vecsum_pT_allJet=MHT.Pt();
            minTree_vecsum_pT_l1l2_allJet=(MHT+p_l1+p_l2).Pt();
            minTree_mass_l1l2_allJet=(MHT+p_l1+p_l2).M();
            minTree_ratio_vecsumpTlep_vecsumpTjet=(p_l1+p_l2).Pt()/MHT.Pt();
            minTree_PtLLgen=(*genLepton+*genAntiLepton).Pt();
            
            minTree_DNN_MET_pT=-1;
            minTree_DNN_MET_phi=-4;
            minTree_DNN_MET_dPhi_nextLep=-1;
            
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
               minTree_Jet2_pt=cjets[1].p.Pt();
               minTree_Jet2_phi=cjets[1].p.Phi();
               minTree_Jet2_eta=cjets[1].p.Eta();
               minTree_Jet2_E=cjets[1].p.E();
               minTree_Jet2_bTagScore=cjets[1].bTagDeepCSV;
               minTree_mLL=(p_l1+p_l2).M();
               minTree_PtLL=(p_l1+p_l2).Pt();
               minTree_mjj=(cjets[0].p+cjets[1].p).M();
               minTree_totalWeight=lumi_weight*fEventWeight*SFWeight;
            }
            else{
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
               minTree_PuppiMet_xy=-1.;
               minTree_Met_xy=-1.;
               minTree_DeepMet_reso=-1.;
               minTree_DeepMet_resp=-1.;
               minTree_MHT=-1.;
               minTree_HT=-1.;
               minTree_HT_phi=-4.;
               minTree_MT=-1.;
               minTree_MT_nextLep=-1.;
               minTree_SF=0.;
               minTree_totalWeight=minTree_N;
               minTree_PhiPtnunuMet=-1;
               minTree_PhiRecPuppi=-1;
               minTree_PhiRecPuppi_xy=-1;
               minTree_PhiRec_xy=-1;
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
               minTree_CaloMET_phi=-5;
               minTree_PFMETxy=-1;
               minTree_PFMETxy_phi=-5;
               minTree_PuppiMET_phi=-5;
               minTree_PuppiMET_xy_phi=-5;
               minTree_MET_xy_phi=-5;
               minTree_PFMET_phi=-5;
               minTree_DeepMet_reso_phi=-5;
               minTree_DeepMet_resp_phi=-5;
               minTree_MT2=-1;
               minTree_vecsum_pT_allJet=-1;
               minTree_vecsum_pT_l1l2_allJet=-1;
               minTree_mass_l1l2_allJet=-1;
               minTree_ratio_vecsumpTlep_vecsumpTjet=-1;
               minTree_mjj=-1;
               minTree_PtLL=-1;

            }
            if (pseudo_selection==false) {
               minTree_PtNuNu=-1.;
               minTree_PhiGen=-1.;
               minTree_PhiNuNu=-1.;
               minTree_genMT=-1;
               minTree_genMT_nextLep=-1;
               minTree_PhiPtnunuMet=-1;
               minTree_leadTop=-1;
               minTree_dPhiNuNu=-1;
               minTree_PtNuNu_phi=-5;
               minTree_nGenJets=-1;
               minTree_PtLLgen=-1;
            }
            else {
               minTree_PhiPtnunuMet=abs(neutrinoPair.DeltaPhi(MET->p));
            }
            
            //Evaluate DNN Regression
            TLorentzVector DNN_MET;
            float DNN_MET_x;
            float DNN_MET_y;
            float DNN_MET_pT=-1.;
            float DNN_MET_dPhi_nextLep=5.;
            if (rec_selection){
               if(cfg.applyDNN){
                  //ONNX
                  input_vec[0]=minTree_PuppiMet_xy*sin(minTree_PuppiMET_xy_phi);
                  input_vec[1]=minTree_PuppiMet_xy*cos(minTree_PuppiMET_xy_phi);
                  input_vec[2]=minTree_Met_xy*sin(minTree_MET_xy_phi);
                  input_vec[3]=minTree_Met_xy*cos(minTree_MET_xy_phi);
                  input_vec[4]=minTree_vecsum_pT_allJet*sin(minTree_HT_phi);
                  input_vec[5]=minTree_vecsum_pT_allJet*cos(minTree_HT_phi);
                  input_vec[6]=minTree_mass_l1l2_allJet;
                  input_vec[7]=minTree_Jet1_pt*sin(minTree_Jet1_phi);
                  input_vec[8]=minTree_MHT;
                  input_vec[9]=minTree_Lep1_pt*cos(minTree_Lep1_phi);
                  input_vec[10]=minTree_Lep1_pt*sin(minTree_Lep1_phi);
                  input_vec[11]=minTree_Jet1_pt*cos(minTree_Jet1_phi);
                  input_vec[12]=minTree_mjj;
                  input_vec[13]=minTree_Jet1_E;
                  input_vec[14]=minTree_HT;
                  input_vec[15]=minTree_Jet2_pt*sin(minTree_Jet2_phi);
                  input_vec[16]=minTree_Jet2_pt*cos(minTree_Jet2_phi);
                                    
                  dnnRegression.evaluate(input_vec,output_vec);
                  DNN_MET_x = minTree_PuppiMet_xy*cos(minTree_PuppiMET_phi)-output_vec[0];
                  DNN_MET_y = minTree_PuppiMet_xy*sin(minTree_PuppiMET_phi)-output_vec[1];
                  
                  DNN_MET.SetXYZM(DNN_MET_x,DNN_MET_y,0.,0.);
                  DNN_MET_pT=DNN_MET.Pt();
                  
                  
                  for (TLorentzVector const lep : {p_l1,p_l2}){
                     const float dPhi=DNN_MET.DeltaPhi(lep);
                     if (std::abs(dPhi) < std::abs(DNN_MET_dPhi_nextLep)) {
                        DNN_MET_dPhi_nextLep=abs(dPhi);
                     }
                  }
                  minTree_DNN_MET_pT=DNN_MET_pT;
                  minTree_DNN_MET_phi=DNN_MET.Phi();
                  minTree_DNN_MET_dPhi_nextLep=DNN_MET_dPhi_nextLep;
               }
            }
            
            //apply MET cut for SF channels based on DNN MET (not applied if syst removeMetCut, which is used for DNN training)
            if (currentSystematic.type()!=Systematic::removeMetCut){
               if ((channel[0] || channel[1]) && DNN_MET_pT<40) rec_selection=false;
            }
            
            //fill event in minimal tree and set reco values to non physical values if rec selection is not fulfilled
            if(rec_selection==false && pseudo_selection==false) continue; //fill minTrees only with reco or gen selected events
            if(!rec_selection){
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
               minTree_PuppiMet_xy=-1.;
               minTree_Met_xy=-1.;
               minTree_DeepMet_reso=-1.;
               minTree_DeepMet_resp=-1.;
               minTree_MHT=-1.;
               minTree_HT=-1.;
               minTree_HT_phi=-4.;
               minTree_MT=-1.;
               minTree_MT_nextLep=-1.;
               minTree_SF=0.;
               minTree_totalWeight=minTree_N;
               minTree_PhiPtnunuMet=-1;
               minTree_PhiRecPuppi=-1;
               minTree_PhiRecPuppi_xy=-1;
               minTree_PhiRec_xy=-1;
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
               minTree_CaloMET_phi=-5;
               minTree_PFMETxy=-1;
               minTree_PFMETxy_phi=-5;
               minTree_PuppiMET_phi=-5;
               minTree_PuppiMET_xy_phi=-5;
               minTree_MET_xy_phi=-5;
               minTree_PFMET_phi=-5;
               minTree_DeepMet_reso_phi=-5;
               minTree_DeepMet_resp_phi=-5;
               minTree_MT2=-1;
               minTree_vecsum_pT_allJet=-1;
               minTree_vecsum_pT_l1l2_allJet=-1;
               minTree_mass_l1l2_allJet=-1;
               minTree_ratio_vecsumpTlep_vecsumpTjet=-1;
               minTree_mjj=-1;
               minTree_PtLL=-1;
            }
            ttbar_res.Fill();
            
            if(rec_selection==false) continue;  //fill the following histograms only with events selected by the reco baseline selection
            
            if(isData) hs_cutflow.fillweight("cutflow/"+cat,5,1);
            else hs_cutflow.fillweight("cutflow/"+cat,5,*w_pu * mcWeight * leptonSFweight * *w_prefiring * bTagWeight * jetPileupIDWeight * *w_topPT *triggerSF);
            
            
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
            TString path_cat=cat;
                        
            hs.fill("baseline/"+path_cat+"/MET"                         ,met);
            hs.fill("baseline/"+path_cat+"/PuppiMET"                    ,met_puppi);
            hs.fill("baseline/"+path_cat+"/PuppiMET_xy"                 ,MET_Puppi_xy->p.Pt());
            hs.fill("baseline/"+path_cat+"/MET_xy"                      ,MET_xy->p.Pt());
            hs.fill("baseline/"+path_cat+"/met1000"                     ,met);
            hs.fill("baseline/"+path_cat+"/DNN_MET_pT"                  ,DNN_MET_pT);
            hs.fill("baseline/"+path_cat+"/DNN_MET_dPhi_nextLep"        ,DNN_MET_dPhi_nextLep);
            hs.fill("baseline/"+path_cat+"/DeepMET_reso"                ,DeepMET_reso->p.Pt());
            hs.fill("baseline/"+path_cat+"/DeepMET_resp"                ,DeepMET_resp->p.Pt());
            hs.fill("baseline/"+path_cat+"/Lep1_pt"                     ,p_l1.Pt());
            hs.fill("baseline/"+path_cat+"/Lep2_pt"                     ,p_l2.Pt());
            hs.fill("baseline/"+path_cat+"/pTsumlep"                    ,(p_l1+p_l2).Pt());
            hs.fill("baseline/"+path_cat+"/sumpTlep"                    ,p_l1.Pt()+p_l2.Pt());
            hs.fill("baseline/"+path_cat+"/pTbJet"                      ,BJets[0].p.Pt());
            hs.fill("baseline/"+path_cat+"/bJet_eta"                    ,BJets[0].p.Eta());
            hs.fill("baseline/"+path_cat+"/Jet1_pt"                     ,cjets[0].p.Pt());
            hs.fill("baseline/"+path_cat+"/Jet2_pt"                     ,cjets[1].p.Pt());
            hs.fill("baseline/"+path_cat+"/dphi_metNearLep"             ,abs(dPhiMETnearLep));
            hs.fill("baseline/"+path_cat+"/dphi_metNearLep_puppi"       ,abs(dPhiMETnearLepPuppi));
            hs.fill("baseline/"+path_cat+"/dphi_metNearLep_puppi_xy"    ,abs(dPhiMETnearLepPuppi_xy));
            hs.fill("baseline/"+path_cat+"/dphi_metNearLep_xy"          ,abs(dPhiMETnearLep_xy));
            hs.fill("baseline/"+path_cat+"/COSdphi_metNearLep"          ,TMath::Cos(abs(dPhiMETnearLep)));
            hs.fill("baseline/"+path_cat+"/SINdphi_metNearLep"          ,TMath::Sin(abs(dPhiMETnearLep)));
            hs.fill("baseline/"+path_cat+"/dPhiLep1bJet"                ,abs(dPhiLep1BJet));
            hs.fill("baseline/"+path_cat+"/dR_bJetLep1"                 ,abs(dRLep1BJet));
            hs.fill("baseline/"+path_cat+"/dphi_bJetLep2"               ,abs(dPhiLep2BJet));
            hs.fill("baseline/"+path_cat+"/dphi_bJetnearLep"            ,dphi_bJetnearLep);
            hs.fill("baseline/"+path_cat+"/dphi_b1b2"                   ,abs(dPhiBjets));
            hs.fill("baseline/"+path_cat+"/dR_b1b2"                     ,abs(dRBjets));
            hs.fill("baseline/"+path_cat+"/dphi_metLep1"                ,abs(dPhiLep1MET));
            hs.fill("baseline/"+path_cat+"/dphi_metLep2"                ,abs(dPhiLep2MET));
            hs.fill("baseline/"+path_cat+"/dphi_metLepsum"              ,abs(dPhiMetLepSum));
            hs.fill("baseline/"+path_cat+"/dR_Lep1Lep2"                 ,abs(dR_Lep1Lep2));
            hs.fill("baseline/"+path_cat+"/nBjets"                      ,nBjets);
            hs.fill("baseline/"+path_cat+"/C_em_W_p"                    ,mt2+0.2*(200-met_puppi));
            hs.fill("baseline/"+path_cat+"/C_em_W_m"                    ,mt2-0.2*(200-met_puppi));
            hs.fill("baseline/"+path_cat+"/mt_MetLep2"                  ,mt_MetLep2);
            hs.fill("baseline/"+path_cat+"/mt_MetNextLep"               ,mt_MetNextLep);
            hs.fill("baseline/"+path_cat+"/conMt_Lep1Lep2"              ,conMt_Lep1Lep2);
            hs.fill("baseline/"+path_cat+"/ST"                          ,ST);
            hs.fill("baseline/"+path_cat+"/HT"                          ,HT);
            hs.fill("baseline/"+path_cat+"/sum_STHT"                    ,ST+HT);
            hs.fill("baseline/"+path_cat+"/sum_mlb"                     ,sum_mlb);
            hs.fill("baseline/"+path_cat+"/Lep1_phi"                    ,p_l1.Phi());
            hs.fill("baseline/"+path_cat+"/Lep2_phi"                    ,p_l2.Phi());
            hs.fill("baseline/"+path_cat+"/Jet1_phi"                    ,cjets[0].p.Phi());
            hs.fill("baseline/"+path_cat+"/Jet2_phi"                    ,cjets[1].p.Phi());
            hs.fill("baseline/"+path_cat+"/PFMET_phi"                   ,MET->p.Phi());
            hs.fill("baseline/"+path_cat+"/PuppiMET_phi"                ,MET_Puppi->p.Phi());
            hs.fill("baseline/"+path_cat+"/PuppiMET_xy_phi"             ,MET_Puppi_xy->p.Phi());
            hs.fill("baseline/"+path_cat+"/MET_xy_phi"                  ,MET_xy->p.Phi());
            hs.fill("baseline/"+path_cat+"/DeepMET_reso_phi"            ,DeepMET_reso->p.Phi());
            hs.fill("baseline/"+path_cat+"/DeepMET_resp_phi"            ,DeepMET_resp->p.Phi());
            hs.fill("baseline/"+path_cat+"/CaloMET"                     ,MET_Calo->p.Pt());
            hs.fill("baseline/"+path_cat+"/CaloMET_phi"                 ,MET_Calo->p.Phi());
            
            hs.fill("baseline/"+path_cat+"/PuppiMET*cos(PuppiMET_phi)"                 ,met_puppi*cos(MET_Puppi->p.Phi()));
            hs.fill("baseline/"+path_cat+"/PuppiMET*sin(PuppiMET_phi)"                 ,met_puppi*sin(MET_Puppi->p.Phi()));
            hs.fill("baseline/"+path_cat+"/METunc_Puppi"                               ,MET_Puppi->uncertainty);
            hs.fill("baseline/"+path_cat+"/MET*cos(PFMET_phi)"                         ,met*cos(MET->p.Phi()));
            hs.fill("baseline/"+path_cat+"/MET*sin(PFMET_phi)"                         ,met*sin(MET->p.Phi()));
            hs.fill("baseline/"+path_cat+"/vecsum_pT_allJet*cos(HT_phi)"               ,MHT.Pt()*cos(MHT.Phi()));
            hs.fill("baseline/"+path_cat+"/vecsum_pT_allJet*sin(HT_phi)"               ,MHT.Pt()*sin(MHT.Phi()));
            hs.fill("baseline/"+path_cat+"/nJets"                                      ,cjets.size());
            hs.fill("baseline/"+path_cat+"/n_Interactions"                             ,*n_Interactions);
            hs.fill("baseline/"+path_cat+"/Lep1_flavor"                                ,flavor_l1);
            hs.fill("baseline/"+path_cat+"/Lep2_flavor"                                ,flavor_l2);
            hs.fill("baseline/"+path_cat+"/Lep1_pt*cos(Lep1_phi)"                      ,p_l1.Pt()*cos(p_l1.Phi()));
            hs.fill("baseline/"+path_cat+"/Lep1_pt*sin(Lep1_phi)"                      ,p_l1.Pt()*sin(p_l1.Phi()));
            hs.fill("baseline/"+path_cat+"/Lep1_eta"                                   ,p_l1.Eta());
            hs.fill("baseline/"+path_cat+"/Lep1_E"                                     ,p_l1.E());
            hs.fill("baseline/"+path_cat+"/Lep2_pt*cos(Lep2_phi)"                      ,p_l2.Pt()*cos(p_l2.Phi()));
            hs.fill("baseline/"+path_cat+"/Lep2_pt*sin(Lep2_phi)"                      ,p_l2.Pt()*sin(p_l2.Phi()));
            hs.fill("baseline/"+path_cat+"/Lep2_eta"                                   ,p_l2.Eta());
            hs.fill("baseline/"+path_cat+"/Lep2_E"                                     ,p_l2.E());
            hs.fill("baseline/"+path_cat+"/Jet1_pt*cos(Jet1_phi)"                      ,cjets[0].p.Pt()*cos(cjets[0].p.Phi()));
            hs.fill("baseline/"+path_cat+"/Jet1_pt*sin(Jet1_phi)"                      ,cjets[0].p.Pt()*sin(cjets[0].p.Phi()));
            hs.fill("baseline/"+path_cat+"/Jet1_eta"                                   ,cjets[0].p.Eta());
            hs.fill("baseline/"+path_cat+"/Jet1_E"                                     ,cjets[0].p.E());
            hs.fill("baseline/"+path_cat+"/Jet2_pt*cos(Jet2_phi)"                      ,cjets[1].p.Pt()*cos(cjets[1].p.Phi()));
            hs.fill("baseline/"+path_cat+"/Jet2_pt*sin(Jet2_phi)"                      ,cjets[1].p.Pt()*sin(cjets[1].p.Phi()));
            hs.fill("baseline/"+path_cat+"/Jet2_eta"                                   ,cjets[1].p.Eta());
            hs.fill("baseline/"+path_cat+"/Jet2_E"                                     ,cjets[1].p.E());
            hs.fill("baseline/"+path_cat+"/dPhiMETnearJet"                             ,abs(dPhiMETnearJet));
            hs.fill("baseline/"+path_cat+"/dPhiMETfarJet"                              ,abs(dPhiMETfarJet));
            hs.fill("baseline/"+path_cat+"/dPhiMETleadJet"                             ,abs(dPhiMETleadJet));
            hs.fill("baseline/"+path_cat+"/dPhiMETlead2Jet"                            ,abs(dPhiMETlead2Jet));
            hs.fill("baseline/"+path_cat+"/dPhiMETbJet"                                ,abs(dPhiMetBJet));
            hs.fill("baseline/"+path_cat+"/dPhiLep1Lep2"                               ,abs(dPhiLep1Lep2));
            hs.fill("baseline/"+path_cat+"/dPhiJet1Jet2"                               ,abs(dPhiJet1Jet2));
            hs.fill("baseline/"+path_cat+"/METsig"                                     ,MET->sig);
            hs.fill("baseline/"+path_cat+"/MHT"                                        ,MHT.M());
            hs.fill("baseline/"+path_cat+"/MT"                                         ,mt_MetLep1);
            hs.fill("baseline/"+path_cat+"/looseLeptonVeto"                            ,*addLepton? 1: 0);
            hs.fill("baseline/"+path_cat+"/dPhiMETnearJet_Puppi"                       ,abs(dPhiMETnearJet_Puppi));
            hs.fill("baseline/"+path_cat+"/dPhiMETfarJet_Puppi"                        ,abs(dPhiMETfarJet_Puppi));
            hs.fill("baseline/"+path_cat+"/dPhiMETleadJet_Puppi"                       ,abs(dPhiMETleadJet_Puppi));
            hs.fill("baseline/"+path_cat+"/dPhiMETlead2Jet_Puppi"                      ,abs(dPhiMETlead2Jet_Puppi));
            hs.fill("baseline/"+path_cat+"/dPhiMETbJet_Puppi"                          ,abs(dPhiMetBJet_Puppi));
            hs.fill("baseline/"+path_cat+"/dPhiLep1bJet"                               ,abs(dPhiLep1BJet));
            hs.fill("baseline/"+path_cat+"/dPhiLep1Jet1"                               ,abs(dPhiLep1Jet1));
            hs.fill("baseline/"+path_cat+"/mLL"                                        ,minTree_mLL);
            hs.fill("baseline/"+path_cat+"/CaloMET*cos(CaloMET_phi)"                   ,MET_Calo->p.Pt()*cos(MET_Calo->p.Phi()));
            hs.fill("baseline/"+path_cat+"/CaloMET*sin(CaloMET_phi)"                   ,MET_Calo->p.Pt()*sin(MET_Calo->p.Phi()));
            hs.fill("baseline/"+path_cat+"/MET_xy*cos(MET_xy_phi)"                     ,MET_xy->p.Pt()*cos(MET_xy->p.Phi()));
            hs.fill("baseline/"+path_cat+"/MET_xy*sin(MET_xy_phi)"                     ,MET_xy->p.Pt()*sin(MET_xy->p.Phi()));
            hs.fill("baseline/"+path_cat+"/PuppiMET_xy*cos(PuppiMET_xy_phi)"           ,MET_Puppi_xy->p.Pt()*cos(MET_Puppi_xy->p.Phi()));
            hs.fill("baseline/"+path_cat+"/PuppiMET_xy*sin(PuppiMET_xy_phi)"           ,MET_Puppi_xy->p.Pt()*sin(MET_Puppi_xy->p.Phi()));
            hs.fill("baseline/"+path_cat+"/DeepMET_reso*cos(DeepMET_reso_phi)"         ,DeepMET_reso->p.Pt()*cos(DeepMET_reso->p.Phi()));
            hs.fill("baseline/"+path_cat+"/DeepMET_reso*sin(DeepMET_reso_phi)"         ,DeepMET_reso->p.Pt()*sin(DeepMET_reso->p.Phi()));
            hs.fill("baseline/"+path_cat+"/DeepMET_resp*cos(DeepMET_resp_phi)"         ,DeepMET_resp->p.Pt()*cos(DeepMET_resp->p.Phi()));
            hs.fill("baseline/"+path_cat+"/DeepMET_resp*sin(DeepMET_resp_phi)"         ,DeepMET_resp->p.Pt()*sin(DeepMET_resp->p.Phi()));
            hs.fill("baseline/"+path_cat+"/MT2"                                        ,mt2);
            hs.fill("baseline/"+path_cat+"/vecsum_pT_allJet"                           ,minTree_vecsum_pT_allJet);
            hs.fill("baseline/"+path_cat+"/vecsum_pT_l1l2_allJet"                      ,minTree_vecsum_pT_l1l2_allJet);
            hs.fill("baseline/"+path_cat+"/mass_l1l2_allJet"                           ,minTree_mass_l1l2_allJet);
            hs.fill("baseline/"+path_cat+"/ratio_vecsumpTlep_vecsumpTjet"              ,minTree_ratio_vecsumpTlep_vecsumpTjet);
            hs.fill("baseline/"+path_cat+"/mjj"                                        ,minTree_mjj);
            
            hs.fill("baseline/"+path_cat+"/DNNMET*cos(DNNMET_phi)"                     ,DNN_MET_x);
            hs.fill("baseline/"+path_cat+"/DNNMET*sin(DNNMET_phi)"                     ,DNN_MET_y);
            
            hs.fill("genParticles/"+path_cat+"/pT_nunu"                 ,neutrinoPair.Pt());
            hs.fill("genParticles/"+path_cat+"/genMET"                  ,genMet);
            hs.fill("genParticles/"+path_cat+"/DMgenMet"                ,DMgenMET.Pt());
            hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep"        ,abs(dPhiNeutrinoLep1));
            hs.fill("genParticles/"+path_cat+"/dphi_NeutrinoLep"        ,abs(dPhiNeutrinoLep2));
            hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep"          ,abs(dRNeutrinoLep1));
            hs.fill("genParticles/"+path_cat+"/dR_NeutrinoLep"          ,abs(dRNeutrinoLep2));
            hs.fill("genParticles/"+path_cat+"/pTtop1"                  ,pT_top1);
            hs.fill("genParticles/"+path_cat+"/pTtop2"                  ,pT_top2);
            hs.fill("genParticles/"+path_cat+"/genHT"                   ,*HTgen);
            hs.fill("genParticles/"+path_cat+"/n_Interactions_gen"      ,*n_Interactions_gen);
            
            
            hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep",met,abs(dPhiMETnearLep));
            hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep_Puppi",met_puppi,abs(dPhiMETnearLepPuppi));
            hs2d.fill("baseline/"+path_cat+"/2d_MetVSdPhiMetNearLep_DNN",DNN_MET_pT,abs(DNN_MET_dPhi_nextLep));
            hs2d.fill("genParticles/"+path_cat+"/2d_PtNuNuVSdPhiNuNuNearLep",neutrinoPair.Pt(),abs(dPhiPtNunearLep));
            hs2d.fill("baseline/"+path_cat+"/2d_CemVSdPhiMetNearLep",mt2-0.2*(200-met_puppi),abs(dPhiMETnearLepPuppi));
            
            if(channel[2]){
               if (muonLead){
                  hs.fill("baseline/"+path_cat+"/Lep_e_pt",p_l2.Pt());
                  hs.fill("baseline/"+path_cat+"/Lep_mu_pt",p_l1.Pt());
               }
               else{
                  hs.fill("baseline/"+path_cat+"/Lep_e_pt",p_l1.Pt());
                  hs.fill("baseline/"+path_cat+"/Lep_mu_pt",p_l2.Pt());
               }
            }
            
            //Fill hists for 2D GOF studies
            // First set values
            dim2D_vec[0].value = MET_Puppi_xy->p.Pt()*cos(MET_Puppi_xy->p.Phi());
            dim2D_vec[1].value = MET_Puppi_xy->p.Pt()*sin(MET_Puppi_xy->p.Phi());
            dim2D_vec[2].value = MET_xy->p.Pt()*cos(MET_xy->p.Phi());
            dim2D_vec[3].value = MET_xy->p.Pt()*sin(MET_xy->p.Phi());
            dim2D_vec[4].value = MHT.Pt()*cos(MHT.Phi());
            dim2D_vec[5].value = MHT.Pt()*sin(MHT.Phi());
            dim2D_vec[6].value = minTree_mass_l1l2_allJet;
            dim2D_vec[7].value = cjets[0].p.Pt()*sin(cjets[0].p.Phi());
            dim2D_vec[8].value = MHT.M();
            dim2D_vec[9].value = p_l1.Pt()*cos(p_l1.Phi());
            dim2D_vec[10].value = p_l1.Pt()*sin(p_l1.Phi());
            dim2D_vec[11].value = cjets[0].p.Pt()*cos(cjets[0].p.Phi());
            dim2D_vec[12].value = MET_Calo->p.Pt();
            dim2D_vec[13].value = mt2;
            dim2D_vec[14].value = minTree_mjj;
            dim2D_vec[15].value = cjets.size();
            dim2D_vec[16].value = cjets[0].p.E();
            dim2D_vec[17].value = HT;
            dim2D_vec[18].value = cjets[1].p.Pt()*cos(cjets[1].p.Phi());
            dim2D_vec[19].value = cjets[1].p.Pt()*sin(cjets[1].p.Phi());
            
            dim2D_vec[20].value = DeepMET_reso->p.Pt()*cos(DeepMET_reso->p.Phi());
            dim2D_vec[21].value = DeepMET_reso->p.Pt()*sin(DeepMET_reso->p.Phi());
            dim2D_vec[22].value = DeepMET_resp->p.Pt()*cos(DeepMET_resp->p.Phi());
            dim2D_vec[23].value = DeepMET_resp->p.Pt()*sin(DeepMET_resp->p.Phi());
            
            // Now loop over 2D hists an fill
            for (int i=0; i<dim2D_vec.size(); i++){
               for (int j=(i+1); j<dim2D_vec.size(); j++){
                  hs2d_GOF.fill("baseline_GOF2D/"+path_cat+"/"+dim2D_vec[i].name+"_VS_"+dim2D_vec[j].name,dim2D_vec[i].value,dim2D_vec[j].value);
               }
            }
            
                  
         }// evt loop
         io::log<<"";
         
         //Derive runtime
         auto stopTime_evtLoop = high_resolution_clock::now();
         auto duration = duration_cast<milliseconds>(stopTime_evtLoop - startTime_evtLoop);
         totalTime = duration.count();
         timePerEvent = totalTime/(1.0*iEv);
         
         hs.scaleLumi(currentSystematic);
         hs_cutflow.scaleLumi(currentSystematic);
         hs2d.scaleLumi(currentSystematic);
         hs2d_GOF.scaleLumi(currentSystematic);
         hs.mergeOverflow();
         hs_cutflow.mergeOverflow();
         hs2d.mergeOverflow();
         hs2d_GOF.mergeOverflow();
         file->Close();
         
         //Save ntuple for TTbar resolution used in binning studies (only if last file of sample e.g. for SingleTop)
         if (cfg.datasets.getDatasubsets({dss.datasetName}).back().name==dss.name) {
            ttbar_res_saver.save(ttbar_res,dss.datasetName);
            ttbar_res.Reset();
         }
         
         //For multi save dss name
         dssName_multi=TString(dss.datasetName);
      
      }  // datasubset loop
      
   } // dataset loop
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs.combineFromSubsamples(samplesToCombine);
   hs_cutflow.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   hs2d_GOF.combineFromSubsamples(samplesToCombine);
   
   // Save histograms
   TString loc=TString::Format("hists/histograms_%s%s.root",cfg.treeVersion.Data(),("_"+currentSystematic.name()).Data());
   if(cfg.multi && cfg.processFraction<1) loc=TString::Format("multiHists/test/%shistograms_%s%s_%s.root",(currentSystematic.name()+"/").Data(),dssName_multi.Data(),(cfg.fileNR==0)?TString("").Data():TString("_"+std::to_string(cfg.fileNR)).Data(),cfg.treeVersion.Data());
   else if (cfg.multi) loc=TString::Format("multiHists/%shistograms_%s%s_%s.root",(currentSystematic.name()+"/").Data(),dssName_multi.Data(),(cfg.fileNR==0)?TString("").Data():TString("_"+std::to_string(cfg.fileNR)).Data(),cfg.treeVersion.Data());
   io::RootFileSaver saver_hist(loc,TString::Format("distributions%.1f",cfg.processFraction*100),false);
   hs.saveHistograms(saver_hist,samplesToCombine);
   hs_cutflow.saveHistograms(saver_hist,samplesToCombine);
   hs2d.saveHistograms(saver_hist,samplesToCombine);
   
   // Save histograms for 2D GOF studies (transform 2D to 1D first)
   hs2d_GOF.saveHistograms2D_as1D(saver_hist,samplesToCombine);
   
   // Print runtime for event loop
   std::cout<<"Event loop had a total runtime of "<<totalTime<<"ms correspond to "<<timePerEvent<<"ms per event"<<std::endl;
   
}
