//Script to extract special distributions to study met resoultion 
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
   std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8"});
   // ~std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8"});
   // ~std::vector<TString> vsDatasubsets({"ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4","ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M2T4"});
   // ~std::vector<TString> vsDatasubsets({"DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8","DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_merged"});
   // ~std::vector<TString> vsDatasubsets({"WW_TuneCUETP8M1_13TeV-pythia8_merged"});
   // ~std::vector<TString> vsDatasubsets({"SMS-T2tt_mStop-650_mLSP-350_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"});
   
   std::vector<TString> samplesToCombine={"TTbar_diLepton"};
   // ~std::vector<TString> samplesToCombine={"TTbar_diLepton_CUETP8M2"};
   // ~std::vector<TString> samplesToCombine={"SingleTop"};
   // ~std::vector<TString> samplesToCombine={"DrellYan"};
   // ~std::vector<TString> samplesToCombine={"WW"};
   // ~std::vector<TString> samplesToCombine={"T2tt_650_350"};
   
   hist::Histograms<TH2F> hs2D(vsDatasubsets);    //Define histograms in the following
   hs2D.addHist("baseline/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhigenMETLep"   ,";|#Delta#phi|(GenMet,nearest gen l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120/Met_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);met",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120/genMet_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);genmet",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120/GenMetDiffMETRelReco_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/|p_{T}^{miss}",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_matchedLep_met120/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_matchedLep_met120/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_genmet120/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_genmet120/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met230/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met230/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met200/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met200/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_230/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_ee/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_ee/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_emu/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_emu/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_mumu/GenMetDiffMETRel_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_mumu/MetSig_dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   
   hs2D.addHist("baseline/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/GenMetDiffMETRel_dPhigenMETLep_Puppi"   ,";|#Delta#phi|(GenMet,nearest gen l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120/Met_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);met",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120/GenMetDiffMETRelReco_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/|p_{T}^{miss}",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_matchedLep_met120/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_matchedLep_met120/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_genmet120/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_genmet120/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met230/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met200/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met200/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_230/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_ee/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_ee/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_emu/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_emu/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   hs2D.addHist("baseline_met120_mumu/GenMetDiffMETRel_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",20,0,3.14,6000,-5,5);
   hs2D.addHist("baseline_met120_mumu/MetSig_dPhiMETLep_Puppi"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",20,0,3.14,6000,0,1000);
   
   hs2D.addHist("baseline_met120/nVertex"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);nVertices",20,0,3.14,6000,0,60);
   
   hs2D.addHist("baseline/diffMET_BJetLBRegr_diffGenMET_PtNuNu"   ,";MET-MET(BJetRegrLB);GenMET-p_{T}^{#nu#nu}",1000,-100,100,1000,-100,100);
   hs2D.addHist("baseline/diffMET_BJetLBRegrMan_diffGenMET_PtNuNu"   ,";MET-MET(BJetRegrLBman);GenMET-p_{T}^{#nu#nu}",1000,-100,100,1000,-100,100);
   
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hs.addHist("baseline/dPhiMETLep"   ,";|#Delta#phi|(p_{T}^{miss},nearest l);EventsBIN",200,0,3.14);
   hs.addHist("baseline/dPhiMETLep_gen"   ,";|#Delta#phi|(genMET,nearest gen l);EventsBIN",200,0,3.14);
   
   hs.addHist("baseline/METres"   ,";MET-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline/METresPuppi"   ,";MET(Puppi)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline/METresBJetRegr"   ,";MET(BJetRegr)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline/METresBJetLBRegr"   ,";MET(BJetRegrLB)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline/METresBJetLBRegr_manual"   ,";MET(BJetRegrLB_man)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   
   hs.addHist("baseline_genmet120/METres"   ,";MET-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_genmet120/METresPuppi"   ,";MET(Puppi)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_genmet120/METresBJetRegr"   ,";MET(BJetRegr)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_genmet120/METresBJetLBRegr"   ,";MET(BJetRegrLB)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   
   hs.addHist("baseline_met120/METres"   ,";MET-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_met120/METresPuppi"   ,";MET(Puppi)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_met120/METresBJetRegr"   ,";MET(BJetRegr)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   hs.addHist("baseline_met120/METresBJetLBRegr"   ,";MET(BJetRegrLB)-p_{T}^{#nu#nu};EventsBIN",1000,-150,150);
   
   for(auto const set : vsDatasubsets) {
      auto const dss = cfg.datasets.getDatasubset(set);
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      
      hs2D.setCurrentSample(dss.name);
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
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      // ~TTreeReaderValue<tree::MET> MET_BReg(reader, "metBJetRegression");
      // ~TTreeReaderValue<tree::MET> MET_BRegLB(reader, "metBJetRegressionLoose");
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
      TTreeReaderValue<TLorentzVector> genLepton_noPseudo(reader, "genLepton");
      TTreeReaderValue<TLorentzVector> genAntiLepton_noPseudo(reader, "genAntiLepton");
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
         hs.setFillWeight(fEventWeight);
         
         float const met=MET->p.Pt();
         float const met_puppi=MET_Puppi->p.Pt();
         float const genMet=GENMET->p.Pt();
         
         bool rec_selection=true;
         bool pseudo_selection=true;
         
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
         // ~std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
         std::vector<bool> ttbarSelection=selection::ttbarSelection_looseJetID(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
         if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end(), [](bool v) { return v; })) rec_selection=false;
         
         // end reco baseline selection
         
         if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
               
         if(rec_selection==false) continue;  //fill the following histograms only with events selected by the reco baseline selection
         
         if (isnan(MET_Puppi->uncertainty)) {
            std::cout<<MET_Puppi->uncertainty<<std::endl;
            std::cout<<*lumNo<<":"<<*evtNo<<std::endl;
         }
         
         float dPhiMETnearLep=4;
         float dPhiMETnearLep_Puppi=4;
         for (TLorentzVector const lep : {p_l1,p_l2}){
            const float dPhi=MET->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
               dPhiMETnearLep=dPhi;
            }
            const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
            if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLep_Puppi)) {
               dPhiMETnearLep_Puppi=dPhi_Puppi;
            }
         }
         
         float dPhiMETnearLep_gen=4;
         for (TLorentzVector const lep : {*genLepton_noPseudo,*genAntiLepton_noPseudo}){
            const float dPhi=GENMET->p.DeltaPhi(lep);
            if (std::abs(dPhi) < std::abs(dPhiMETnearLep_gen)) {
               dPhiMETnearLep_gen=dPhi;
            }
         }
         
         //manually correct met for bJetRegression
         TLorentzVector oldBJets(0,0,0,0);
         TLorentzVector newBJets(0,0,0,0);
         for(auto jet : *jets){
         // ~for(auto jet : BJets){
            if(jet.bTagDeepCSV>0.2217){
               oldBJets+=jet.p;
               TLorentzVector temp=jet.p;
               temp.SetPtEtaPhiE(jet.p.Pt()*jet.bJetRegressionCorr, jet.p.Eta(), jet.p.Phi(), jet.p.Energy()*jet.bJetRegressionCorr);
               newBJets+=temp;
            }
         }
         
         TLorentzVector met_manual_BJetRegr=MET->p-newBJets+oldBJets;
         
         
         TLorentzVector neutrinoPair(0,0,0,0);
         neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
         
         bool matchedLeptons=(matchLepton(p_l1,*genLepton_noPseudo) || matchLepton(p_l1,*genAntiLepton_noPseudo)) && (matchLepton(p_l2,*genLepton_noPseudo) || matchLepton(p_l2,*genAntiLepton_noPseudo));
         
         hs.fill("baseline/dPhiMETLep",abs(dPhiMETnearLep));
         hs.fill("baseline/dPhiMETLep_gen",abs(dPhiMETnearLep_gen));
         
         if(pseudo_selection){
            hs.fill("baseline/METres",met-neutrinoPair.Pt());
            hs.fill("baseline/METresPuppi",MET_Puppi->p.Pt()-neutrinoPair.Pt());
            // ~hs.fill("baseline/METresBJetRegr",MET_BReg->p.Pt()-neutrinoPair.Pt());
            // ~hs.fill("baseline/METresBJetLBRegr",MET_BRegLB->p.Pt()-neutrinoPair.Pt());
            hs.fill("baseline/METresBJetLBRegr_manual",met_manual_BJetRegr.Pt()-neutrinoPair.Pt());
            // ~hs2D.fill("baseline/diffMET_BJetLBRegr_diffGenMET_PtNuNu",met-MET_BRegLB->p.Pt(),genMet-neutrinoPair.Pt());
            hs2D.fill("baseline/diffMET_BJetLBRegrMan_diffGenMET_PtNuNu",met-met_manual_BJetRegr.Pt(),genMet-neutrinoPair.Pt());
         }
         
         hs2D.fill("baseline/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep_Puppi),(genMet-met)/genMet);
         hs2D.fill("baseline/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
         
         if (genMet>120) {
            hs2D.fill("baseline_genmet120/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_genmet120/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_genmet120/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_genmet120/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            
            if(pseudo_selection){
               hs.fill("baseline_genmet120/METres",met-neutrinoPair.Pt());
               hs.fill("baseline_genmet120/METresPuppi",MET_Puppi->p.Pt()-neutrinoPair.Pt());
               // ~hs.fill("baseline_genmet120/METresBJetRegr",MET_BReg->p.Pt()-neutrinoPair.Pt());
               // ~hs.fill("baseline_genmet120/METresBJetLBRegr",MET_BRegLB->p.Pt()-neutrinoPair.Pt());
            }
         }
         
         bool pfMET120 = (met>120);
         bool puppi120 = (MET_Puppi->p.Pt()>120);
         if (!pfMET120 && !puppi120) continue;
         
         if (pfMET120) {
            hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhigenMETLep",abs(dPhiMETnearLep_gen),(genMet-met)/genMet);
            hs2D.fill("baseline_met120/GenMetDiffMETRelReco_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/met);
            hs2D.fill("baseline_met120/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met120/Met_dPhiMETLep",abs(dPhiMETnearLep),met);
            hs2D.fill("baseline_met120/genMet_dPhiMETLep",abs(dPhiMETnearLep),genMet);
            hs2D.fill("baseline_met120/nVertex",abs(dPhiMETnearLep),*n_Interactions);
         }
         if (puppi120) {
            hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met120/GenMetDiffMETRel_dPhigenMETLep_Puppi",abs(dPhiMETnearLep_gen),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met120/GenMetDiffMETRelReco_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/MET_Puppi->p.Pt());
            hs2D.fill("baseline_met120/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            hs2D.fill("baseline_met120/Met_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->p.Pt());
            
            
            if(pseudo_selection){
               hs.fill("baseline_met120/METres",met-neutrinoPair.Pt());
               hs.fill("baseline_met120/METresPuppi",MET_Puppi->p.Pt()-neutrinoPair.Pt());
               // ~hs.fill("baseline_met120/METresBJetRegr",MET_BReg->p.Pt()-neutrinoPair.Pt());
               // ~hs.fill("baseline_met120/METresBJetLBRegr",MET_BRegLB->p.Pt()-neutrinoPair.Pt());
            }
         }
         
         if (*is_ee){
            if (pfMET120) {
               hs2D.fill("baseline_met120_ee/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
               hs2D.fill("baseline_met120_ee/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            }
            if (puppi120) {
               hs2D.fill("baseline_met120_ee/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
               hs2D.fill("baseline_met120_ee/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            }
         }
         else if (*is_emu){
            if (pfMET120) {
               hs2D.fill("baseline_met120_emu/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
               hs2D.fill("baseline_met120_emu/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            }
            if (puppi120) {
               hs2D.fill("baseline_met120_emu/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
               hs2D.fill("baseline_met120_emu/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            }
         }
         else {
            if (pfMET120) {
               hs2D.fill("baseline_met120_mumu/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
               hs2D.fill("baseline_met120_mumu/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            }
            if (puppi120) {
               hs2D.fill("baseline_met120_mumu/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
               hs2D.fill("baseline_met120_mumu/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            }
         }
         
         if (matchedLeptons) {
            if (pfMET120) {
               hs2D.fill("baseline_matchedLep_met120/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
               hs2D.fill("baseline_matchedLep_met120/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            }
            if (puppi120) {
               hs2D.fill("baseline_matchedLep_met120/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
               hs2D.fill("baseline_matchedLep_met120/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
            }
         }
         
         if (met>200) {
            hs2D.fill("baseline_met200/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met200/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met200/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met200/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
         }
         
         if (met<230) {
            hs2D.fill("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met120_230/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met120_230/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met120_230/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
         }
         else {
            hs2D.fill("baseline_met230/GenMetDiffMETRel_dPhiMETLep",abs(dPhiMETnearLep),(genMet-met)/genMet);
            hs2D.fill("baseline_met230/MetSig_dPhiMETLep",abs(dPhiMETnearLep),MET->sig);
            hs2D.fill("baseline_met230/GenMetDiffMETRel_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),(genMet-MET_Puppi->p.Pt())/genMet);
            hs2D.fill("baseline_met230/MetSig_dPhiMETLep_Puppi",abs(dPhiMETnearLep_Puppi),MET_Puppi->sig);
         }
               
      }// evt loop
      io::log<<"";
      
      hs2D.scaleLumi();
      hs2D.mergeOverflow();
      hs.scaleLumi();
      hs.mergeOverflow();
      file.Close();
   }
   
   hs2D.combineFromSubsamples(samplesToCombine);
   hs.combineFromSubsamples(samplesToCombine);
   
   // Save histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("diff_MET%.1f",cfg.processFraction*100),false);
   hs2D.saveHistograms(saver_hist,samplesToCombine);
   hs.saveHistograms(saver_hist,samplesToCombine);
   
}
