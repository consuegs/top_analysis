//Script to test performance of BJet regression 
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

extern "C"
void run()
{
   std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8"});
   
   std::vector<TString> samplesToCombine={"TTbar_diLepton_CUETP8M2"};
   
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hs.addHist("baseline/pt_genJet"   ,";%pT;EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet"   ,";%pT;EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet_BReg"   ,";%pT;EventsBIN"           ,1000,0,1000);
   
   hs.addHist("baseline/pt_genJet1"   ,";%pT(Jet1);EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet1"   ,";%pT(Jet1);EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet1_BReg"   ,";%pT(Jet1);EventsBIN"           ,1000,0,1000);
   
   hs.addHist("baseline/diff_pt_recoJet1"   ,";%pT(Jet1)-%pT(genJet1);EventsBIN"           ,50,-100,100);
   hs.addHist("baseline/diff_pt_recoJet1_BReg"   ,";%pT(Jet1)-%pT(genJet1);EventsBIN"           ,50,-100,100);
   hs.addHist("baseline/ratio_pt_recoJet1"   ,";%pT(genJet1)/%pT(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline/ratio_pt_recoJet1_BReg"   ,";%pT(genJet1)/%pT(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline_onlyBTag/ratio_pt_recoJet1"   ,";%pT(genJet1)/%pT(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline_onlyBTag/ratio_pt_recoJet1_BReg"   ,";%pT(genJet1)/%pT(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline/ratio_E_recoJet1"   ,";E(genJet1)/E(Jet1);EventsBIN"           ,50,0,2);
   hs.addHist("baseline/ratio_E_recoJet1_BReg"   ,";E(genJet1)/E(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline_onlyBTag/ratio_E_recoJet1"   ,";E(genJet1)/E(Jet1)-;EventsBIN"           ,50,0,2);
   hs.addHist("baseline_onlyBTag/ratio_E_recoJet1_BReg"   ,";E(genJet1)/E(Jet1)-;EventsBIN"           ,50,0,2);
   
   hs.addHist("baseline_vetoJ1noCSVscore/diff_pt_recoJet1"   ,";%pT(Jet1)-%pT(genJet1);EventsBIN"           ,500,-100,100);
   hs.addHist("baseline_vetoJ1noCSVscore/diff_pt_recoJet1_BReg"   ,";%pT(Jet1)-%pT(genJet1);EventsBIN"           ,500,-100,100);
   
   hs.addHist("baseline/pt_genJet2"   ,";%pT(Jet2);EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet2"   ,";%pT(Jet2);EventsBIN"           ,1000,0,1000);
   hs.addHist("baseline/pt_recoJet2_BReg"   ,";%pT(Jet2);EventsBIN"           ,1000,0,1000);
   
   hs.addHist("baseline/diff_MET_PF"   ,";genMET-PFMET;EventsBIN"           ,50,-100,100);
   hs.addHist("baseline/diff_MET_bJetRegression"   ,";genMET-MET_bJetRegression;EventsBIN"           ,50,-100,100);
   hs.addHist("baseline/diff_MET_bJetRegressionLoose"   ,";genMET-MET_bJetRegressionLoose;EventsBIN"           ,50,-100,100);

   
   for(auto const set : vsDatasubsets) {
      auto const dss = cfg.datasets.getDatasubset(set);
      TFile file(dss.getPath(),"read");
      // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v13/TTTo2L2Nu_TuneCUETP8M2_ttHtranche3_13TeV-powheg-pythia8.root","read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");
      
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
      TTreeReaderValue<std::vector<tree::Particle>>      genJets     (reader, "genJets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
      TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      // ~TTreeReaderValue<tree::MET> MET_Deep(reader, "metDeep");
      TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
      TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
      TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
      TTreeReaderValue<tree::MET> MET_BJetRegression(reader, "metBJetRegression");
      TTreeReaderValue<tree::MET> MET_BJetRegressionLoose(reader, "metBJetRegressionLoose");
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
         hs.setFillWeight(fEventWeight);
         
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
         
         std::vector<tree::Jet> cjets=phys::getCleanedJets_looseID(*jets);
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
         
         
         for (tree::Particle const &jet : *genJets){
            hs.fill("baseline/pt_genJet",jet.p.Pt());
         }
         for (tree::Jet const &jet : *jets){
            hs.fill("baseline/pt_recoJet",jet.p.Pt());
            hs.fill("baseline/pt_recoJet_BReg",jet.p.Pt()*jet.bJetRegressionCorr);
         }
         
         if((cjets)[0].p.DeltaPhi((*genJets)[0].p)>0.5 || (cjets)[0].p.DeltaR((*genJets)[0].p)>0.5) continue;
         
         hs.fill("baseline/pt_genJet1",(*genJets)[0].p.Pt());
         hs.fill("baseline/pt_recoJet1",(cjets)[0].p.Pt());
         hs.fill("baseline/pt_recoJet1_BReg",(cjets)[0].p.Pt()*(cjets)[0].bJetRegressionCorr);
         hs.fill("baseline/diff_pt_recoJet1",(cjets)[0].p.Pt()-(*genJets)[0].p.Pt());
         hs.fill("baseline/diff_pt_recoJet1_BReg",(cjets)[0].p.Pt()*(cjets)[0].bJetRegressionCorr-(*genJets)[0].p.Pt());
         hs.fill("baseline/ratio_pt_recoJet1",(*genJets)[0].p.Pt()/((cjets)[0].p.Pt()));
         hs.fill("baseline/ratio_pt_recoJet1_BReg",(*genJets)[0].p.Pt()/((cjets)[0].p.Pt()*(cjets)[0].bJetRegressionCorr));
         hs.fill("baseline/ratio_E_recoJet1",(*genJets)[0].p.E()/((cjets)[0].p.E()));
         hs.fill("baseline/ratio_E_recoJet1_BReg",(*genJets)[0].p.E()/((cjets)[0].p.E()*(cjets)[0].bJetRegressionCorr));
         hs.fill("baseline/diff_MET_PF",genMet-met);
         hs.fill("baseline/diff_MET_bJetRegression",genMet-MET_BJetRegression->p.Pt());
         hs.fill("baseline/diff_MET_bJetRegressionLoose",genMet-MET_BJetRegressionLoose->p.Pt());
         if((cjets)[0].bTagDeepCSV>0) {
            hs.fill("baseline_vetoJ1noCSVscore/diff_pt_recoJet1",(cjets)[0].p.Pt()-(*genJets)[0].p.Pt());
            hs.fill("baseline_vetoJ1noCSVscore/diff_pt_recoJet1_BReg",(cjets)[0].p.Pt()*(cjets)[0].bJetRegressionCorr-(*genJets)[0].p.Pt());
         }
         if((cjets)[0].bTagDeepCSV>0.2213) {
            hs.fill("baseline_onlyBTag/ratio_pt_recoJet1",(*genJets)[0].p.Pt()/((cjets)[0].p.Pt()));
            hs.fill("baseline_onlyBTag/ratio_pt_recoJet1_BReg",(*genJets)[0].p.Pt()/((cjets)[0].p.Pt()*(cjets)[0].bJetRegressionCorr));
            hs.fill("baseline_onlyBTag/ratio_E_recoJet1",(*genJets)[0].p.E()/((cjets)[0].p.E()));
            hs.fill("baseline_onlyBTag/ratio_E_recoJet1_BReg",(*genJets)[0].p.E()/((cjets)[0].p.E()*(cjets)[0].bJetRegressionCorr));
         }
         
         hs.fill("baseline/pt_genJet2",(*genJets)[1].p.Pt());
         hs.fill("baseline/pt_recoJet2",(cjets)[1].p.Pt());
         hs.fill("baseline/pt_recoJet2_BReg",(cjets)[1].p.Pt()*(cjets)[1].bJetRegressionCorr);
               
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs.mergeOverflow();
      file.Close();
   }
   
   hs.combineFromSubsamples(samplesToCombine);
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("test_BJetRegression_v13%.1f",cfg.processFraction*100),false);
   hs.saveHistograms(saver_hist,samplesToCombine);
   
}
