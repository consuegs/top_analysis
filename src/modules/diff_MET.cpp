//Script to extract special distributions to study met resoultion 
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

void saveHistograms(std::map<TString,std::vector<TString>> const &msPresel_vVars, io::RootFileSaver const &saver_hist,hist::Histograms<TH2F> &hs, std::vector<TString> const &Samples)
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
   std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8"});
   // ~std::vector<TString> vsDatasubsets({"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4","ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4"});
   // ~std::vector<TString> vsDatasubsets({"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_merged"});
   // ~std::vector<TString> vsDatasubsets({"WW_TuneCUETP8M1_13TeV-pythia8_merged"});
   // ~std::vector<TString> vsDatasubsets({"SMS-T2tt_mStop-650_mLSP-350_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"});
   
   std::vector<TString> samplesToCombine={"TTbar_diLepton"};
   // ~std::vector<TString> samplesToCombine={"SingleTop"};
   // ~std::vector<TString> samplesToCombine={"DrellYan"};
   // ~std::vector<TString> samplesToCombine={"WW"};
   // ~std::vector<TString> samplesToCombine={"T2tt_650_350"};
   
   hist::Histograms<TH2F> hs2D(vsDatasubsets);    //Define histograms in the following
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met230/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met230/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_230/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met230/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_230/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/MetSig_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met230/MetSig_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_230/MetSig_dPhiMETLep_Deep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);

   
   for(auto const set : vsDatasubsets) {
      auto const dss = cfg.datasets.getDatasubset(set);
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      
      hs2D.setCurrentSample(dss.name);
      
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

      int recoEvents=0;
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
         hs2D.setFillWeight(fEventWeight);
         
         float const met=MET->p.Pt();
         float const genMet=GENMET->p.Pt();
         
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
            // ~if (jet.bTagCSVv2>0.5426) {      //Loose working point for CSVv2 (Should be replaced in the future by deep CSV!!!)
            if (jet.bTagDeepCSV>0.2217) {      //Loose working point for deepCSV
               bTag=true;
               BJets.push_back(jet);
            }
         }
         if (!bTag) rec_selection=false; // end reco baseline selection
         
         if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
               
         if(rec_selection==false) continue;  //fill the following histograms only with events selected by the reco baseline selection
         
         float dPhiMETnearLep=4;
         float dPhiMETnearLep_Puppi=4;
         float dPhiMETnearLep_Deep=4;
         for (TLorentzVector const lep : {p_l1,p_l2}){
            const float dPhi=MET->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
               dPhiMETnearLep=dPhi;
            }
            const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
            if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLep_Puppi)) {
               dPhiMETnearLep_Puppi=dPhi_Puppi;
            }
            const float dPhi_Deep=MET_Deep->p.DeltaPhi(lep);
            if (std::abs(dPhi_Deep) < std::abs(dPhiMETnearLep_Deep)) {
               dPhiMETnearLep_Deep=dPhi_Deep;
            }
         }
         
         if (met<120) continue;
         hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
         hs2D.fill("baseline_met120/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
         hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
         hs2D.fill("baseline_met120/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
         hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),(genMet-MET_Deep->p.Pt())/genMet);
         hs2D.fill("baseline_met120/MetSig_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),MET_Deep->sig);
         
         if (met<230) {
            hs2D.fill("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met120_230/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met120_230/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            hs2D.fill("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),(genMet-MET_Deep->p.Pt())/genMet);
            hs2D.fill("baseline_met120_230/MetSig_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),MET_Deep->sig);
         }
         else {
            hs2D.fill("baseline_met230/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met230/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met230/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            hs2D.fill("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),(genMet-MET_Deep->p.Pt())/genMet);
            hs2D.fill("baseline_met230/MetSig_dPhiMETLep_Deep",abs(dPhiMETnearLep_Deep),MET_Deep->sig);
         }
               
      }// evt loop
      io::log<<"";
      
      hs2D.scaleLumi();
      hs2D.mergeOverflow();
      file.Close();
   }
   
   hs2D.combineFromSubsamples(samplesToCombine);
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline_met120/",{"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep","GenMetDiffMETRel_dPhiMETLep_Puppi","MetSig_dPhiMETLep_Puppi","GenMetDiffMETRel_dPhiMETLep_Deep","MetSig_dPhiMETLep_Deep"}},
      {"baseline_met120_230/",{"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep","GenMetDiffMETRel_dPhiMETLep_Puppi","MetSig_dPhiMETLep_Puppi","GenMetDiffMETRel_dPhiMETLep_Deep","MetSig_dPhiMETLep_Deep"}},
      {"baseline_met230/",{"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep","GenMetDiffMETRel_dPhiMETLep_Puppi","MetSig_dPhiMETLep_Puppi","GenMetDiffMETRel_dPhiMETLep_Deep","MetSig_dPhiMETLep_Deep"}},
      };
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("diff_MET%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs2D,samplesToCombine);
   
}
