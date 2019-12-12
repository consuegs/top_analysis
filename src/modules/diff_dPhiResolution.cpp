//Script to extract distributions separetely for events with good and bad dPhi(MET,next l) resolution
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

extern "C"
void run()
{
   // ~std::vector<TString> vsDatasubsets({"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_merged"});
   std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8"});
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   // ~for (TString selection:{"badPhiRes/","goodPhiRes/"}){
   for (TString selection:{"badMETRes/","goodMETRes/"}){
      hs.addHist(selection+"baseline/met"   ,";%MET;EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
      hs.addHist(selection+"baseline/metSig"   ,";METsignificance;EventsBIN"           ,100,0,30);
      hs.addHist(selection+"baseline/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/etalep1"   ,";|#eta_{l1}|;EventsBIN"           ,100,0,2.5);
      hs.addHist(selection+"baseline/etalep2"   ,";|#eta_{l2}|;EventsBIN"           ,100,0,2.5);
      hs.addHist(selection+"baseline/pTbJet"   ,";p_{T}^{b}(GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/etabJet"   ,";|#eta_{b}|;EventsBIN"           ,100,0,2.5);
      hs.addHist(selection+"baseline/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metLeadJet"   ,";|#Delta#phi|(p_{T}^{miss},lead jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metLead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd lead jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metFarJet"   ,";|#Delta#phi|(p_{T}^{miss},far jet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/COSdphi_metNearLep"   ,";cos(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,-1.,1);
      hs.addHist(selection+"baseline/SINdphi_metNearLep"   ,";sin(|#Delta#phi|(p_{T}^{miss},nearest l));EventsBIN"           ,100,0.,1);
      hs.addHist(selection+"baseline/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},bjet);EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_metLepsum"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}+l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2});EventsBIN"           ,100,0,3.2);
      hs.addHist(selection+"baseline/dR_Lep1Lep2"   ,";|#Delta R|(l_{1},l_{2});EventsBIN"           ,100,0,5);
      hs.addHist(selection+"baseline/nBjets"   ,";N_{bJets};EventsBIN"           ,5,-0.5,4.5);
      hs.addHist(selection+"baseline/njets"   ,";N_{Jets};EventsBIN"           ,9,-0.5,8.5);
      hs.addHist(selection+"baseline/mt2"   ,";MT2 (GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"baseline/mt_MetLep1"   ,";M_{T}(p_{T}^{miss},l_{1}) (GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"baseline/mt_MetLep2"   ,";M_{T}(p_{T}^{miss},l_{2}) (GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"baseline/conMt_Lep1Lep2"   ,";conM_{T}(l_{1},l_{2}) (GeV);EventsBIN"           ,100,0,500);
      hs.addHist(selection+"baseline/ST"   ,";S_{T} (GeV);EventsBIN"           ,100,0,1500);
      hs.addHist(selection+"baseline/HT"   ,";H_{T} (GeV);EventsBIN"           ,100,0,2500);
      hs.addHist(selection+"baseline/METoverHT"   ,";MET/#sqrt{H_{T}p_{T}} ;EventsBIN"           ,100,0,50);
      hs.addHist(selection+"baseline/METoverSUMpt"   ,";MET/#sqrt{#Sum } ;EventsBIN"           ,100,0,30);
      hs.addHist(selection+"baseline/sum_STHT"   ,";S_{T}+H_{T} (GeV);EventsBIN"           ,100,0,4000);
      hs.addHist(selection+"baseline/sum_mlb"   ,";sum m_{lb} (GeV);EventsBIN"           ,100,0,3000);
      hs.addHist(selection+"baseline/diff_met_Puppi"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|(GeV);EventsBIN"           ,100,0,300);
      hs.addHist(selection+"baseline/diff_met_NoHF"   ,";|PFp_{T}^{miss}-NHFp_{T}^{miss}|(GeV);EventsBIN"           ,100,0,300);
      hs.addHist(selection+"baseline/diff_met_Calo"   ,";|PFp_{T}^{miss}-CALOp_{T}^{miss}|(GeV);EventsBIN"           ,100,0,300);
      hs.addHist(selection+"baseline/diff_met_Raw"   ,";|PFp_{T}^{miss}-RAWp_{T}^{miss}|(GeV);EventsBIN"           ,100,0,300);
      hs.addHist(selection+"baseline/HTnormed_diff_met_Puppi"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,1.5);
      hs.addHist(selection+"baseline/HTnormed_diff_met_NoHF"   ,";|PFp_{T}^{miss}-NHFp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,1.5);
      hs.addHist(selection+"baseline/HTnormed_diff_met_Calo"   ,";|PFp_{T}^{miss}-CALOp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,1.5);
      hs.addHist(selection+"baseline/HTnormed_diff_met_Raw"   ,";|PFp_{T}^{miss}-RAWp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,1.5);
      hs.addHist(selection+"baseline/n_Interactions"   ,";n_Interactions;EventsBIN"           ,100,0,70);
      hs.addHist(selection+"genParticles/pT_nunu"   ,";p_{T}^{#nu#nu(+BSM)}(GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"genParticles/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
      hs.addHist(selection+"genParticles/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu(+BSM)}-genMET(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"genParticles/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
      hs.addHist(selection+"genParticles/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-5,5);
      hs.addHist(selection+"genParticles/diff_Met_genMET_normSUM"   ,";#frac{p_{T}^{miss}-genMET}{p_{T}^{miss}+genMET};EventsBIN"           ,200,-2,2);
      hs.addHist(selection+"genParticles/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"genParticles/diff_genMT2_MT2"   ,";genMT2-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"genParticles/diff_genMT2neutrino_MT2"   ,";genMT2neutrino-MT2(GeV);EventsBIN"           ,400,-100,100);
      hs.addHist(selection+"genParticles/diff_dPhiMetNearLep_gen"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(genMet,nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"genParticles/diff_dPhiMetNearLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l)-|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);;EventsBIN" ,400,-3.2,3.2);
      hs.addHist(selection+"genParticles/dphi_NeutrinoLep"   ,";|#Delta#phi|_{gen}(#nu,l);EventsBIN"           ,100,0,4);
      hs.addHist(selection+"genParticles/dR_NeutrinoLep"   ,";|#Delta R|_{gen}(#nu,l);EventsBIN"           ,100,0,6);
      hs.addHist(selection+"genParticles/pTtop1"   ,";p_{T}^{gen}(t_{1});EventsBIN"           ,100,0,600);   
      hs.addHist(selection+"genParticles/pTtop2"   ,";p_{T}^{gen}(t_{2});EventsBIN"           ,100,0,600);
   }

   
   auto const dss = cfg.datasets.getDatasubset(vsDatasubsets[0]);
   TFile file(dss.getPath(),"read");
   if (file.IsZombie()) {
      return;
   }
   io::log * ("Processing '"+dss.datasetName+"' ");
;
   
   hs.setCurrentSample(dss.name);
   
   //Lumi weight for current sample
   float lumi_weight=dss.xsec/float(dss.Ngen)*cfg.lumi;

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
   // ~TTreeReaderValue<tree::MET> MET(reader, "met");
   TTreeReaderValue<tree::MET> MET(reader, "metPuppi");
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
      
      float fEventWeight=*w_pu * *w_mc * *sf_lep1 * *sf_lep2;     //Set event weight also taking lepton scale factors into account
      // ~hs.setFillWeight(fEventWeight);
      hs.setFillWeight(1.);
      
      float const met=MET->p.Pt();
      float const genMet=GENMET->p.Pt();
      
      //Booleans for reco and pseudo selection
      bool rec_selection=true;
      bool pseudo_selection=true;
      
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
         if (jet.bTagCSVv2>0.5426) {      //Loose working point for CSVv2 (Should be replaced in the future by deep CSV!!!)
            bTag=true;
            BJets.push_back(jet);
         }
      }
      if (!bTag) rec_selection=false; // end reco baseline selection
      
      if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
      
      if(rec_selection==false || pseudo_selection==false) continue;  // only proceed with events selected by both of the baseline selection (reco or pseudo)
      
      // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
      TLorentzVector neutrinoPair(0,0,0,0);
      neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
      
      //Get DeltaPhi between MET (or genMet or neutrino pT) and nearest Lepton
      float dPhiMETnearLep=4; 
      for (TLorentzVector const lep : {p_l1,p_l2}){
         const float dPhi=MET->p.DeltaPhi(lep);
         if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) dPhiMETnearLep=dPhi;
      }
      float dPhigenMETnearLep=4;
      for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
         const float dPhi=GENMET->p.DeltaPhi(lep);
         if (std::abs(dPhi) < std::abs(dPhigenMETnearLep)) dPhigenMETnearLep=dPhi;
      }
      float dPhiPtNunearLep=4;
      for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
         const float dPhi=neutrinoPair.DeltaPhi(lep);
         if (std::abs(dPhi) < std::abs(dPhiPtNunearLep)) dPhiPtNunearLep=dPhi;
      }
      
      //Get DeltaPhi between MET nearest/leading Jet and HT
      float dPhiMETnearJet=4;
      float dPhiMETfarJet=0;
      float dPhiMETleadJet=4;
      float dPhiMETlead2Jet=4;
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
      }
      
      if(rec_selection==false) continue;  //fill the following histograms only with events selected by the reco baseline selection
      
      // Bjet and angular variables
      int nBjets=BJets.size();
      float dPhiMetBJet=MET->p.DeltaPhi(BJets[0].p);
      float dPhiLep1BJet=p_l1.DeltaPhi(BJets[0].p);
      float dPhiLep2BJet=p_l2.DeltaPhi(BJets[0].p);
      float dPhiLep1MET=p_l1.DeltaPhi(MET->p);
      float dPhiLep2MET=p_l2.DeltaPhi(MET->p);
      float dPhiLep1Lep2=p_l1.DeltaPhi(p_l2);
      float dR_Lep1Lep2=p_l1.DeltaR(p_l2);
      float dPhiMetLepSum=MET->p.DeltaPhi(p_l1+p_l2);
      
      //Further variables
      float mt_MetLep1=phys::M_T(MET->p,p_l1);
      float mt_MetLep2=phys::M_T(MET->p,p_l2);
      float sum_mlb=phys::sumMlb(p_l1,p_l2,cjets,BJets);
      float conMt_Lep1Lep2=phys::conM_T(p_l1,p_l2);
      
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
      
      //Calculate ST
      float ST=met+p_l1.Pt()+p_l2.Pt();
      float METoverSUMpt=met/sqrt(HT+p_l1.Pt()+p_l2.Pt());
            
      // ~if (met<120 || abs(dPhiMETnearLep)<1.4) continue; //only events in last two bins
      if (met<230 || abs(dPhiMETnearLep)<1.4) continue; //only events in last bin
      if (*genDecayMode>3) continue;  //remove tau events
      // ~if (abs(dPhiMETlead2Jet)<0.8) continue;
      // ~if (abs(dPhiMETnearJet)<0.8) continue;
      // ~if ((abs(met-MET_Puppi->p.Pt())/HT)>0.05) continue;
      // ~TString selection="goodPhiRes";
      // ~if (abs(abs(dPhiMETnearLep)-abs(dPhiPtNunearLep))>0.3) selection="badPhiRes";
      TString selection="goodMETRes";
      // ~if (abs(neutrinoPair.Pt()-met)>50) selection="badMETRes";
      if (neutrinoPair.Pt()<230 || abs(dPhiPtNunearLep)<1.4) selection="badMETRes";
      // ~if (abs(abs(dPhiMETnearLep)-abs(dPhiPtNunearLep))>2) {
         // ~std::cout<<"------------------------------------"<<std::endl;
         // ~std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
         // ~std::cout<<neutrinoPair.Pt()<<"   "<<met<<std::endl;
         // ~std::cout<<neutrinoPair.Phi()<<"   "<<MET->p.Phi()<<std::endl;
         // ~std::cout<<genLepton->Pt()<<"   "<<p_l2.Pt()<<std::endl;
         // ~std::cout<<genLepton->Phi()<<"   "<<p_l2.Phi()<<std::endl;
         // ~std::cout<<genAntiLepton->Pt()<<"   "<<p_l1.Pt()<<std::endl;
         // ~std::cout<<genAntiLepton->Phi()<<"   "<<p_l1.Phi()<<std::endl;
         // ~std::cout<<abs(dPhiPtNunearLep)<<"   "<<abs(dPhiMETnearLep)<<std::endl;
         // ~std::cout<<dPhiMETleadJet<<std::endl;
         // ~break;
      // ~}
      
      hs.fill(selection+"/baseline/met",met);
      hs.fill(selection+"/baseline/met1000",met);
      hs.fill(selection+"/baseline/metSig",MET->sig);
      hs.fill(selection+"/baseline/mll",*mll);
      hs.fill(selection+"/baseline/pTlep1",p_l1.Pt());
      hs.fill(selection+"/baseline/pTlep2",p_l2.Pt());
      hs.fill(selection+"/baseline/etalep1",abs(p_l1.Eta()));
      hs.fill(selection+"/baseline/etalep2",abs(p_l2.Eta()));
      hs.fill(selection+"/baseline/pTbJet",BJets[0].p.Pt());
      hs.fill(selection+"/baseline/etabJet",abs(BJets[0].p.Eta()));
      hs.fill(selection+"/baseline/dphi_metJet",abs(dPhiMETnearJet));
      hs.fill(selection+"/baseline/dphi_metLeadJet",abs(dPhiMETleadJet));
      hs.fill(selection+"/baseline/dphi_metLead2Jet",abs(dPhiMETlead2Jet));
      hs.fill(selection+"/baseline/dphi_metFarJet",abs(dPhiMETfarJet));
      hs.fill(selection+"/baseline/dphi_metNearLep",abs(dPhiMETnearLep));
      hs.fill(selection+"/baseline/COSdphi_metNearLep",TMath::Cos(abs(dPhiMETnearLep)));
      hs.fill(selection+"/baseline/SINdphi_metNearLep",TMath::Sin(abs(dPhiMETnearLep)));
      hs.fill(selection+"/baseline/dphi_metBJet",abs(dPhiMetBJet));
      hs.fill(selection+"/baseline/dphi_bJetLep1",abs(dPhiLep1BJet));
      hs.fill(selection+"/baseline/dphi_bJetLep2",abs(dPhiLep2BJet));
      hs.fill(selection+"/baseline/dphi_metLep1",abs(dPhiLep1MET));
      hs.fill(selection+"/baseline/dphi_metLep2",abs(dPhiLep2MET));
      hs.fill(selection+"/baseline/dphi_metLepsum",abs(dPhiMetLepSum));
      hs.fill(selection+"/baseline/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
      hs.fill(selection+"/baseline/dR_Lep1Lep2",abs(dR_Lep1Lep2));
      hs.fill(selection+"/baseline/nBjets",nBjets);
      hs.fill(selection+"/baseline/njets",cjets.size());
      hs.fill(selection+"/baseline/mt2",*mt2);
      hs.fill(selection+"/baseline/mt_MetLep1",mt_MetLep1);
      hs.fill(selection+"/baseline/mt_MetLep2",mt_MetLep2);
      hs.fill(selection+"/baseline/conMt_Lep1Lep2",conMt_Lep1Lep2);
      hs.fill(selection+"/baseline/ST",ST);
      hs.fill(selection+"/baseline/HT",HT);
      hs.fill(selection+"/baseline/METoverHT",met/sqrt(HT));
      hs.fill(selection+"/baseline/METoverSUMpt",METoverSUMpt);
      hs.fill(selection+"/baseline/sum_STHT",ST+HT);
      hs.fill(selection+"/baseline/sum_mlb",sum_mlb);
      hs.fill(selection+"/baseline/diff_met_Puppi",abs(met-MET_Puppi->p.Pt()));
      hs.fill(selection+"/baseline/diff_met_NoHF",abs(met-MET_NoHF->p.Pt()));
      hs.fill(selection+"/baseline/diff_met_Calo",abs(met-MET_Calo->p.Pt()));
      hs.fill(selection+"/baseline/diff_met_Raw",abs(met-MET_Raw->p.Pt()));
      hs.fill(selection+"/baseline/HTnormed_diff_met_Puppi",abs(met-MET_Puppi->p.Pt())/HT);
      hs.fill(selection+"/baseline/HTnormed_diff_met_NoHF",abs(met-MET_NoHF->p.Pt())/HT);
      hs.fill(selection+"/baseline/HTnormed_diff_met_Calo",abs(met-MET_Calo->p.Pt())/HT);
      hs.fill(selection+"/baseline/HTnormed_diff_met_Raw",abs(met-MET_Raw->p.Pt())/HT);
      hs.fill(selection+"/baseline/n_Interactions",*n_Interactions);
      hs.fill(selection+"/genParticles/pT_nunu",neutrinoPair.Pt());
      hs.fill(selection+"/genParticles/genMet",genMet);
      hs.fill(selection+"/genParticles/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
      hs.fill(selection+"/genParticles/diff_Met_genMET",met-genMet);
      hs.fill(selection+"/genParticles/diff_Met_genMET_norm",(met-genMet)/genMet);
      hs.fill(selection+"/genParticles/diff_Met_genMET_normSUM",(met-genMet)/(met+genMet));
      hs.fill(selection+"/genParticles/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
      hs.fill(selection+"/genParticles/diff_genMT2_MT2",*genMT2-*mt2);
      hs.fill(selection+"/genParticles/diff_genMT2neutrino_MT2",*genMT2neutrino-*mt2);
      hs.fill(selection+"/genParticles/diff_dPhiMetNearLep_gen",abs(dPhigenMETnearLep)-abs(dPhiMETnearLep));
      hs.fill(selection+"/genParticles/diff_dPhiMetNearLep",abs(dPhiPtNunearLep)-abs(dPhiMETnearLep));
      hs.fill(selection+"/genParticles/dphi_NeutrinoLep",abs(dPhiNeutrinoLep1));
      hs.fill(selection+"/genParticles/dphi_NeutrinoLep",abs(dPhiNeutrinoLep2));
      hs.fill(selection+"/genParticles/dR_NeutrinoLep",abs(dRNeutrinoLep1));
      hs.fill(selection+"/genParticles/dR_NeutrinoLep",abs(dRNeutrinoLep2));
      hs.fill(selection+"/genParticles/pTtop1",pT_top1);
      hs.fill(selection+"/genParticles/pTtop2",pT_top2);
            
   }// evt loop
   io::log<<"";
   
   // ~hs.scaleLumi();
   hs.mergeOverflow();
   file.Close();
   
   // ~std::vector<TString> samplesToCombine={"TTbar"};
   std::vector<TString> samplesToCombine={"TTbar_diLepton"};
   hs.combineFromSubsamples(samplesToCombine);
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   // ~std::map<TString,std::vector<TString>> msPresel_vVars={
      // ~{"goodPhiRes/baseline/",{"met","met1000","metSig","mll","pTlep1","pTlep2","pTbJet","dphi_metJet","dphi_metLeadJet","dphi_metFarJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","METoverHT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"goodPhiRes/genParticles/",{"pT_nunu","genMet","diff_ptNuNu_genMET","diff_Met_genMET","diff_Met_genMET_norm","diff_Met_genMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      // ~{"badPhiRes/baseline/",{"met","met1000","metSig","mll","pTlep1","pTlep2","pTbJet","dphi_metJet","dphi_metLeadJet","dphi_metFarJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","METoverHT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"badPhiRes/genParticles/",{"pT_nunu","genMet","diff_ptNuNu_genMET","diff_Met_genMET","diff_Met_genMET_norm","diff_Met_genMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      // ~};
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"goodMETRes/baseline/",{"met","met1000","metSig","mll","pTlep1","pTlep2","etalep1","etalep2","pTbJet","etabJet","dphi_metJet","dphi_metLeadJet",
         "dphi_metLead2Jet","dphi_metFarJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","njets","mt2",
         "dR_Lep1Lep2","ST","HT","METoverHT","METoverSUMpt","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep",
         "COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum","diff_met_Puppi","diff_met_NoHF","diff_met_Calo","diff_met_Raw",
         "HTnormed_diff_met_Puppi","HTnormed_diff_met_NoHF","HTnormed_diff_met_Calo","HTnormed_diff_met_Raw","n_Interactions"}},
      {"goodMETRes/genParticles/",{"pT_nunu","genMet","diff_ptNuNu_genMET","diff_Met_genMET","diff_Met_genMET_norm","diff_Met_genMET_normSUM","diff_ptNuNu_Met",
         "diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"badMETRes/baseline/",{"met","met1000","metSig","mll","pTlep1","pTlep2","etalep1","etalep2","pTbJet","etabJet","dphi_metJet","dphi_metLeadJet",
         "dphi_metLead2Jet","dphi_metFarJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","njets","mt2",
         "dR_Lep1Lep2","ST","HT","METoverHT","METoverSUMpt","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep",
         "COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum","diff_met_Puppi","diff_met_NoHF","diff_met_Calo","diff_met_Raw",
         "HTnormed_diff_met_Puppi","HTnormed_diff_met_NoHF","HTnormed_diff_met_Calo","HTnormed_diff_met_Raw","n_Interactions"}},
      {"badMETRes/genParticles/",{"pT_nunu","genMet","diff_ptNuNu_genMET","diff_Met_genMET","diff_Met_genMET_norm","diff_Met_genMET_normSUM","diff_ptNuNu_Met",
         "diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      };
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("diff_dPhiResolution%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs,samplesToCombine);
   
}
