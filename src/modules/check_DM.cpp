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

   
   TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v02/TTbarDMJets_DiLept_pseudoscalar_Mchi-50_Mphi-50_TuneCUETP8M1_v2_13TeV-madgraphMLM-pythia8.root","read");
   //~ TFile file("/net/data_cms1b/user/dmeuser/top_analysis/2016/v02/SMS-T1tttt_mGluino-1200_mLSP-800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","read");

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


   int iEv=0;
   //~ int events =0;
   int processEvents=10;
   while (reader.Next()){
      iEv++;
      if (iEv>processEvents) break;
      
      float const met=MET->p.Pt();
      float const genMet=GENMET->p.Pt();
      
      //Baseline selection
      TLorentzVector p_l1;
      TLorentzVector p_l2;
      
      if (*is_ee){
         if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) continue; //currently double check since trees only have tight leptons!!
         p_l1=(*electrons)[0].p;
         p_l2=(*electrons)[1].p;
      }
      else if (*is_mumu){
         if(!(*muons)[0].isTight || !(*muons)[1].isTight) continue;
         p_l1=(*muons)[0].p;
         p_l2=(*muons)[1].p;
      }
      else if (*is_emu){
         if(!(*muons)[0].isTight || !(*electrons)[0].isTight) continue;
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
      
      // Get pT of Neutrino Pair
      TLorentzVector neutrinoPair(0,0,0,0);
      for (tree::IntermediateGenParticle const &inter: *intermediateGenParticles){
         if (abs(inter.pdgId)==24) {
            for (tree::GenParticle const &daughter: inter.daughters){
               if (abs(daughter.pdgId)==12 || abs(daughter.pdgId)==14 || abs(daughter.pdgId)==16){
                  neutrinoPair+=daughter.p;
               }
            }
         }
      }
      
      //Get DeltaPhi between MET and nearest Jet
      float dPhiMETnearJet=4; // nearest jet or photon
      for (tree::Jet const &jet : cjets){
         const float dPhi=MET->p.DeltaPhi(jet.p);
         if (std::abs(dPhi) < std::abs(dPhiMETnearJet))
            dPhiMETnearJet=dPhi;
      }
      
      std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
      TLorentzVector DMgenMET=neutrinoPair;
      for (auto const &genParticle : *genParticles){
         std::cout<<genParticle.pdgId<<std::endl;
         if(abs(genParticle.pdgId)>100000){
            DMgenMET+=genParticle.p;
         }
      }
      std::cout<<met<<std::endl;
      std::cout<<genMet<<std::endl;
      std::cout<<DMgenMET.Pt()<<std::endl;
      
      

      std::cout<<"---------------------------------------------"<<std::endl;
   }// evt loop
   io::log<<"";
   
   file.Close();
   
}
