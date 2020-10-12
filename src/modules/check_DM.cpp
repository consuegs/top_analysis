//Script to check genParticle content of samples and other stuff (mainly for debugging)
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
#include <iostream>
#include <fstream>

Config const &cfg=Config::get();

extern "C"
void run()
{

   
   // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v03/TTbarDMJets_DiLept_pseudoscalar_Mchi-50_Mphi-50_TuneCUETP8M1_v2_13TeV-madgraphMLM-pythia8.root","read");
   // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v05/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","read");
   // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v05/SMS-T2tt_mStop-850_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","read");
   // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v05/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root","read");
   TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v08/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8.root","read");
   // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v05/TTZToLLNuNu_M-10_TuneCUETP8M1_13TeV-amcatnlo-pythia8.root","read");

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
   TTreeReaderValue<float> HTgen(reader, "genHt");
   TTreeReaderValue<bool> is_ee   (reader, "ee");
   TTreeReaderValue<bool> is_emu   (reader, "emu");
   TTreeReaderValue<bool> is_mumu   (reader, "mumu");
   TTreeReaderValue<float> mll   (reader, "mll");
   TTreeReaderValue<TLorentzVector> genNeutrino(reader, "genNeutrino");
   TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "genAntiNeutrino");
   TTreeReaderValue<TLorentzVector> pseudoNeutrino(reader, "pseudoNeutrino");
   TTreeReaderValue<TLorentzVector> pseudoAntiNeutrino(reader, "pseudoAntiNeutrino");
   TTreeReaderValue<TLorentzVector> genLepton(reader, "genLepton");
   TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "genAntiLepton");
   TTreeReaderValue<TLorentzVector> pseudoLepton(reader, "pseudoLepton");
   TTreeReaderValue<TLorentzVector> pseudoAntiLepton(reader, "pseudoAntiLepton");
   TTreeReaderValue<int> genDecayMode_pseudo(reader, "ttbarPseudoDecayMode");
   
   TH1F hist_neutrinos("","",6,-0.5,5.5);
   TH1F hist_neutrinosPt("","",20,0,150);


   int iEv=0;
   //~ int events =0;
   int processEvents=1000;
   while (reader.Next()){
      iEv++;
      if (iEv>processEvents) break;
      
      // ~if (*lumNo!=2681) continue;
      
      float const met=MET->p.Pt();
      float const genMet=GENMET->p.Pt();
      
      //Baseline selection
      TLorentzVector p_l1;
      TLorentzVector p_l2;
      
      if (*is_ee){
         if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) continue; //currently double check since trees only have tight leptons!!
         if(abs((*electrons)[0].etaSC)>2.4 || abs((*electrons)[1].etaSC)>2.4) continue; //To use same region as for muons, cut on supercluster eta
         p_l1=(*electrons)[0].p;
         p_l2=(*electrons)[1].p;
      }
      else if (*is_mumu){
         if(!(*muons)[0].isTight || !(*muons)[1].isTight) continue;
         if((*muons)[0].rIso>0.15 || (*muons)[1].rIso>0.15) continue;
         if(abs((*muons)[0].p.Eta())>2.4 || abs((*muons)[1].p.Eta())>2.4) continue;
         p_l1=(*muons)[0].p;
         p_l2=(*muons)[1].p;
      }
      else if (*is_emu){
         if(!(*muons)[0].isTight || !(*electrons)[0].isTight) continue;
         if((*muons)[0].rIso>0.15 ) continue;
         if(abs((*muons)[0].p.Eta())>2.4) continue;
         if(abs((*electrons)[0].etaSC)>2.4 ) continue;
         if ((*muons)[0].p.Pt()>(*electrons)[0].p.Pt()){
            p_l1=(*muons)[0].p;
            p_l2=(*electrons)[0].p;
         }
         else {
            p_l1=(*electrons)[0].p;
            p_l2=(*muons)[0].p;
         }
      }
      
      if (p_l1.Pt()<25 || p_l2.Pt()<20) continue; //eta cuts already done in TreeWriter
      if (*mll<20 || ((*is_ee || *is_mumu) && *mll<106 && *mll>76)) continue;
      if ((*is_ee || *is_mumu) && met<40) continue;
      
      std::vector<tree::Jet> cjets=phys::getCleanedJets(*jets);
      if (cjets.size()<2) continue;
      
      bool bTag=false;
      int nBjets=0;
      for (tree::Jet const &jet : cjets) {
         if (jet.bTagCSVv2>0.5426) {
            bTag=true;
            nBjets++;
         }
      }
      if (!bTag) continue; // end baseline selection
      
      // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
      TLorentzVector neutrinoPair(0,0,0,0);
      // ~neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
      neutrinoPair=(*pseudoNeutrino)+(*pseudoAntiNeutrino);
      
      //Calculate MET before parton shower and hadronization for DM scenario and SUSY scenarios
      for (auto const &genParticle : *genParticles){
         if((abs(genParticle.pdgId)==1000022)){
            neutrinoPair+=genParticle.p;
         }
      }
      
      if(*genDecayMode_pseudo==1){
         std::cout<<neutrinoPair.Pt()<<std::endl;
      }
      
      // ~std::cout<<"------------------------------"<<std::endl;
      // ~std::cout<<genNeutrino->Pt()<<"   "<<pseudoNeutrino->Pt()<<std::endl;
      // ~std::cout<<genAntiNeutrino->Pt()<<"   "<<pseudoAntiNeutrino->Pt()<<std::endl;
      // ~std::cout<<"------------------------------"<<std::endl;
      // ~std::cout<<genLepton->Pt()<<"   "<<pseudoLepton->Pt()<<std::endl;
      // ~std::cout<<genAntiLepton->Pt()<<"   "<<pseudoAntiLepton->Pt()<<std::endl;
      
      // ~std::cout<<"---------------------"<<std::endl;
      // ~std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
      // ~std::cout<<genNeutrino->Pt()<<"   "<<genAntiNeutrino->Pt()<<std::endl;
      // ~std::cout<<neutrinoPair.Pt()<<std::endl;
      // ~std::cout<<genMet<<std::endl;

      
      // ~int n_neutrinos=0;
      // ~for (auto const &genParticle : *genParticles){
         // ~if (abs(genParticle.pdgId)==12 || abs(genParticle.pdgId)==14 || abs(genParticle.pdgId)==16){
            // ~n_neutrinos++;
            // ~if (n_neutrinos>2)neutrinoPair+=genParticle.p;
            // ~std::cout<<"neutrino"<<"   "<<genParticle.p.Pt()<<"   "<<genParticle.p.Phi()<<"   "<<genParticle.p.Eta()<<std::endl;
            // ~if(n_neutrinos>2)hist_neutrinosPt.Fill(genParticle.p.Pt());
         // ~}
         // ~if((abs(genParticle.pdgId)==1000022)){
            // ~std::cout<<"BSM"<<"   "<<genParticle.p.Pt()<<"   "<<genParticle.p.Phi()<<"   "<<genParticle.p.Eta()<<std::endl;
         // ~}
      // ~}
      // ~std::cout<<neutrinoPair.Pt()<<std::endl;
      
      // ~TLorentzVector temp;
      // ~temp.SetPtEtaPhiM(2.780,-1.285,2.419,0.);
      // ~neutrinoPair+=temp;
      
      //Get DeltaPhi between MET and nearest Jet
      // ~float dPhiMETnearJet=4; // nearest jet or photon
      // ~for (tree::Jet const &jet : cjets){
         // ~const float dPhi=MET->p.DeltaPhi(jet.p);
         // ~if (std::abs(dPhi) < std::abs(dPhiMETnearJet))
            // ~dPhiMETnearJet=dPhi;
      // ~}
      
      // ~std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
      // ~TLorentzVector DMgenMET=neutrinoPair;
      // ~for (auto const &genParticle : *genParticles){
         // ~//~ std::cout<<genParticle.pdgId<<std::endl;
         // ~//~ if(abs(genParticle.pdgId)>100000){
         // ~if(abs(genParticle.pdgId)==1000022){
            // ~DMgenMET+=genParticle.p;
            // ~std::cout<<abs(genParticle.pdgId)<<"   "<<genParticle.p.Pt()<<std::endl;
         // ~}
      // ~}
      //~ std::cout<<met<<std::endl;
      // ~std::cout<<genMet<<std::endl;
      // ~std::cout<<DMgenMET.Pt()<<std::endl;
      
      // ~hist_neutrinos.Fill(n_neutrinos);
      // ~std::cout<<n_neutrinos<<std::endl;
      // ~std::cout<<"---------------------------------------------"<<std::endl;
   }// evt loop
   io::log<<"";
   
   file.Close();
   
   hist_neutrinos.SaveAs("test.root");
   // ~hist_neutrinosPt.SaveAs("test.root");
   
}
