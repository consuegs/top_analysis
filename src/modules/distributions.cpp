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
   
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hs.addHist("baseline/ee/met"   ,";%MET;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/met"   ,";%MET;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/met"   ,";%MET;EventsBIN"           ,100,0,600);
   
   hs.addHist("baseline/ee/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
   hs.addHist("baseline/emu/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
   hs.addHist("baseline/mumu/met1000"   ,";%MET;EventsBIN"           ,100,0,1000);
   
   hs.addHist("baseline/ee/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   
   hs.addHist("baseline/ee/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/pTlep1"   ,";%pTl1;EventsBIN"           ,100,0,600);
   
   hs.addHist("baseline/ee/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/pTlep2"   ,";%pTl2;EventsBIN"           ,100,0,600);
   
   hs.addHist("baseline/ee/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_metJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},bjet),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_metBJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_bJetLep1"   ,";|#Delta#phi|(b Jet,l_{1}),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_bJetLep2"   ,";|#Delta#phi|(b Jet,l_{2}),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_metLep1"   ,";|#Delta#phi|(p_{T}^{miss},l_{1}),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_metLep2"   ,";|#Delta#phi|(p_{T}^{miss},l_{2}),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/emu/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2}),;EventsBIN"           ,100,0,3.2);
   hs.addHist("baseline/mumu/dphi_Lep1Lep2"   ,";|#Delta#phi|(l_{1},l_{2}),;EventsBIN"           ,100,0,3.2);
   
   hs.addHist("baseline/ee/nBjets"   ,";N_{bJets},;EventsBIN"           ,5,-0.5,4.5);
   hs.addHist("baseline/emu/nBjets"   ,";N_{bJets},;EventsBIN"           ,5,-0.5,4.5);
   hs.addHist("baseline/mumu/nBjets"   ,";N_{bJets},;EventsBIN"           ,5,-0.5,4.5);
   
   hs.addHist("genParticles/ee/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/emu/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/mumu/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   
   hs.addHist("genParticles/ee/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/emu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/mumu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   
   hs.addHist("genParticles/ee/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/emu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/mumu/DMgenMet"   ,";DMgenMET(GeV);EventsBIN"           ,100,0,600);
   
   hs.addHist("genParticles/ee/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu}-genMET(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/emu/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu}-genMET(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/mumu/diff_ptNuNu_genMET"   ,";p_{T}^{#nu#nu}-genMET(GeV);EventsBIN"           ,400,-100,100);
   
   hs.addHist("genParticles/ee/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/emu/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/mumu/diff_ptNuNu_DMgenMET"   ,";p_{T}^{#nu#nu}-DMgenMET(GeV);EventsBIN"           ,400,-100,100);
   
   hs.addHist("genParticles/ee/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
   hs.addHist("genParticles/emu/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
   hs.addHist("genParticles/mumu/diff_Met_genMET"   ,";p_{T}^{miss}-genMET(GeV);EventsBIN"           ,200,-100,100);
   
   hs.addHist("genParticles/ee/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
   hs.addHist("genParticles/emu/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
   hs.addHist("genParticles/mumu/diff_Met_DMgenMET"   ,";p_{T}^{miss}-DMgenMET(GeV);EventsBIN"           ,200,-100,100);
   
   hs.addHist("genParticles/ee/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-10,10);
   hs.addHist("genParticles/emu/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-10,10);
   hs.addHist("genParticles/mumu/diff_Met_genMET_norm"   ,";#frac{p_{T}^{miss}-genMET}{genMET};EventsBIN"           ,200,-10,10);
   
   hs.addHist("genParticles/ee/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-10,10);
   hs.addHist("genParticles/emu/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-10,10);
   hs.addHist("genParticles/mumu/diff_Met_DMgenMET_norm"   ,";#frac{p_{T}^{miss}-DMgenMET}{DMgenMET};EventsBIN"           ,200,-10,10);
   
   hs.addHist("genParticles/ee/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/emu/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
   hs.addHist("genParticles/mumu/diff_ptNuNu_Met"   ,";p_{T}^{#nu#nu}-p_{T}^{miss}(GeV);EventsBIN"           ,400,-100,100);
   
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   hs2d.addHist("genParticles/ee/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSDMgenMet", ";DMgenMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/MetVSgenMet", ";p_{T}^{miss} (GeV);genMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/MetVSDMgenMet", ";p_{T}^{miss} (GeV);DMgenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
   
   hs2d.addHist("genParticles/ee/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSMet", ";p_{T}^{miss} (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);


   
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
      int processEvents=cfg.processFraction*dss.entries;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }
         
         float fEventWeight=*w_pu * *w_mc;
         hs.setFillWeight(fEventWeight);
         hs2d.setFillWeight(fEventWeight);
         
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
         std::vector<tree::Jet> BJets;
         for (tree::Jet const &jet : cjets) {
            if (jet.bTagCSVv2>0.5426) {
               bTag=true;
               BJets.push_back(jet);
            }
         }
         if (!bTag) continue; // end baseline selection
         
         // Bjet and angular variables
         int nBjets=BJets.size();
         float dPhiMetBJet=MET->p.DeltaPhi(BJets[0].p);
         float dPhiLep1BJet=p_l1.DeltaPhi(BJets[0].p);
         float dPhiLep2BJet=p_l2.DeltaPhi(BJets[0].p);
         float dPhiLep1MET=p_l1.DeltaPhi(MET->p);
         float dPhiLep2MET=p_l2.DeltaPhi(MET->p);
         float dPhiLep1Lep2=p_l1.DeltaPhi(p_l2);
         
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
         
         //Calculate genMET for DM scenario
         TLorentzVector DMgenMET=neutrinoPair;
         for (auto const &genParticle : *genParticles){
            if(abs(genParticle.pdgId)>100000){
               DMgenMET+=genParticle.p;
            }
         }
         
         //Fill hists
         if (*is_ee){
            hs.fill("baseline/ee/met",met);
            hs.fill("baseline/ee/met1000",met);
            hs.fill("baseline/ee/mll",*mll);
            hs.fill("baseline/ee/pTlep1",p_l1.Pt());
            hs.fill("baseline/ee/pTlep2",p_l2.Pt());
            hs.fill("baseline/ee/dphi_metJet",abs(dPhiMETnearJet));
            hs.fill("baseline/ee/dphi_metBJet",abs(dPhiMetBJet));
            hs.fill("baseline/ee/dphi_bJetLep1",abs(dPhiLep1BJet));
            hs.fill("baseline/ee/dphi_bJetLep2",abs(dPhiLep2BJet));
            hs.fill("baseline/ee/dphi_metLep1",abs(dPhiLep1MET));
            hs.fill("baseline/ee/dphi_metLep2",abs(dPhiLep2MET));
            hs.fill("baseline/ee/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
            hs.fill("baseline/ee/nBjets",nBjets);
            hs.fill("genParticles/ee/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/ee/genMet",genMet);
            hs.fill("genParticles/ee/DMgenMet",DMgenMET.Pt());
            hs.fill("genParticles/ee/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
            hs.fill("genParticles/ee/diff_ptNuNu_DMgenMET",neutrinoPair.Pt()-DMgenMET.Pt());
            hs.fill("genParticles/ee/diff_Met_genMET",met-genMet);
            hs.fill("genParticles/ee/diff_Met_DMgenMET",met-DMgenMET.Pt());
            hs.fill("genParticles/ee/diff_Met_genMET_norm",(met-genMet)/genMet);
            hs.fill("genParticles/ee/diff_Met_DMgenMET_norm",(met-DMgenMET.Pt())/DMgenMET.Pt());
            hs.fill("genParticles/ee/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
            hs2d.fill("genParticles/ee/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
            hs2d.fill("genParticles/ee/2d_nunuVSDMgenMet",genMet,DMgenMET.Pt());
            hs2d.fill("genParticles/ee/MetVSgenMet",met,genMet);
            hs2d.fill("genParticles/ee/MetVSDMgenMet",met,DMgenMET.Pt());
            hs2d.fill("genParticles/ee/2d_nunuVSMet",met,neutrinoPair.Pt());
         }
         else if (*is_emu){
            hs.fill("baseline/emu/met",met);
            hs.fill("baseline/emu/met1000",met);
            hs.fill("baseline/emu/mll",*mll);
            hs.fill("baseline/emu/pTlep1",p_l1.Pt());
            hs.fill("baseline/emu/pTlep2",p_l2.Pt());
            hs.fill("baseline/emu/dphi_metJet",abs(dPhiMETnearJet));
            hs.fill("baseline/emu/dphi_metBJet",abs(dPhiMetBJet));
            hs.fill("baseline/emu/dphi_bJetLep1",abs(dPhiLep1BJet));
            hs.fill("baseline/emu/dphi_bJetLep2",abs(dPhiLep2BJet));
            hs.fill("baseline/emu/dphi_metLep1",abs(dPhiLep1MET));
            hs.fill("baseline/emu/dphi_metLep2",abs(dPhiLep2MET));
            hs.fill("baseline/emu/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
            hs.fill("baseline/emu/nBjets",nBjets);
            hs.fill("genParticles/emu/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/emu/genMet",genMet);
            hs.fill("genParticles/emu/DMgenMet",DMgenMET.Pt());
            hs.fill("genParticles/emu/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
            hs.fill("genParticles/emu/diff_ptNuNu_DMgenMET",neutrinoPair.Pt()-DMgenMET.Pt());
            hs.fill("genParticles/emu/diff_Met_genMET",met-genMet);
            hs.fill("genParticles/emu/diff_Met_DMgenMET",met-DMgenMET.Pt());
            hs.fill("genParticles/emu/diff_Met_genMET_norm",(met-genMet)/genMet);
            hs.fill("genParticles/emu/diff_Met_DMgenMET_norm",(met-DMgenMET.Pt())/DMgenMET.Pt());
            hs.fill("genParticles/emu/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
            hs2d.fill("genParticles/emu/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
            hs2d.fill("genParticles/emu/2d_nunuVSDMgenMet",genMet,DMgenMET.Pt());
            hs2d.fill("genParticles/emu/MetVSgenMet",met,genMet);
            hs2d.fill("genParticles/emu/MetVSDMgenMet",met,DMgenMET.Pt());
            hs2d.fill("genParticles/emu/2d_nunuVSMet",met,neutrinoPair.Pt());
         }
         else if (*is_mumu){
            hs.fill("baseline/mumu/met",met);
            hs.fill("baseline/mumu/met1000",met);
            hs.fill("baseline/mumu/mll",*mll);
            hs.fill("baseline/mumu/pTlep1",p_l1.Pt());
            hs.fill("baseline/mumu/pTlep2",p_l2.Pt());
            hs.fill("baseline/mumu/dphi_metJet",abs(dPhiMETnearJet));
            hs.fill("baseline/mumu/dphi_metBJet",abs(dPhiMetBJet));
            hs.fill("baseline/mumu/dphi_bJetLep1",abs(dPhiLep1BJet));
            hs.fill("baseline/mumu/dphi_bJetLep2",abs(dPhiLep2BJet));
            hs.fill("baseline/mumu/dphi_metLep1",abs(dPhiLep1MET));
            hs.fill("baseline/mumu/dphi_metLep2",abs(dPhiLep2MET));
            hs.fill("baseline/mumu/dphi_Lep1Lep2",abs(dPhiLep1Lep2));
            hs.fill("baseline/mumu/nBjets",nBjets);
            hs.fill("genParticles/mumu/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/mumu/genMet",genMet);
            hs.fill("genParticles/mumu/DMgenMet",DMgenMET.Pt());
            hs.fill("genParticles/mumu/diff_ptNuNu_genMET",neutrinoPair.Pt()-genMet);
            hs.fill("genParticles/mumu/diff_ptNuNu_DMgenMET",neutrinoPair.Pt()-DMgenMET.Pt());
            hs.fill("genParticles/mumu/diff_Met_genMET",met-genMet);
            hs.fill("genParticles/mumu/diff_Met_DMgenMET",met-DMgenMET.Pt());
            hs.fill("genParticles/mumu/diff_Met_genMET_norm",(met-genMet)/genMet);
            hs.fill("genParticles/mumu/diff_Met_DMgenMET_norm",(met-DMgenMET.Pt())/DMgenMET.Pt());
            hs.fill("genParticles/mumu/diff_ptNuNu_Met",neutrinoPair.Pt()-met);
            hs2d.fill("genParticles/mumu/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
            hs2d.fill("genParticles/mumu/2d_nunuVSDMgenMet",genMet,DMgenMET.Pt());
            hs2d.fill("genParticles/mumu/MetVSgenMet",met,genMet);
            hs2d.fill("genParticles/mumu/MetVSDMgenMet",met,DMgenMET.Pt());
            hs2d.fill("genParticles/mumu/2d_nunuVSMet",met,neutrinoPair.Pt());
         }
         
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs2d.scaleLumi();
      hs.mergeOverflow();
      hs2d.mergeOverflow();
      file.Close();
   } // dataset loop
   
   std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      "T1tttt_1200_800","T2tt_650_250","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"};
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   //Plotting part 1D
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
      {"baseline/emu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
      {"baseline/mumu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
      {"genParticles/ee/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_ptNuNu_Met"}},
      {"genParticles/emu/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_ptNuNu_Met"}},
      {"genParticles/mumu/",{"pT_nunu","genMet","DMgenMet","diff_ptNuNu_genMET","diff_ptNuNu_DMgenMET","diff_Met_genMET","diff_Met_DMgenMET","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_ptNuNu_Met"}},
      };
      
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc;
         loc=sPresel+sVar;
         TH1F* hist=hs.getHistogram(loc,"TTbar");
         gfx::LegendEntries le=hs.getLegendEntries();
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         TLatex label=gfx::cornerLabel(cat,2);
         hist->SetStats(0);
         hist->SetMarkerSize(0);
         if (sVar=="diff_ptNuNu_genMET") can.SetLogy(0);
         hist->Draw("histE");
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         label.Draw();
         saver.save(can,"tt_only/"+loc);
         
         THStack st_mc=hs.getStack(loc,{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar"});
         le=hs.getLegendEntries();
         st_mc.Draw();
         TLegend leg2=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg2.Draw();
         label.Draw();
         saver.save(can,"all/"+loc);
      }
   }
   
   // Save 1d histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),false);
   saveHistograms(msPresel_vVars,saver_hist,hs,samplesToCombine);
   
   //Plotting part 2D
   std::map<TString,std::vector<TString>> msPresel_vVars2D={
      {"genParticles/ee/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSgenMet","MetVSDMgenMet","2d_nunuVSMet"}},
      {"genParticles/emu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet"}},
      {"genParticles/mumu/",{"2d_nunuVSgenMet","2d_nunuVSDMgenMet","MetVSDMgenMet","MetVSgenMet","2d_nunuVSMet"}},
      };
   TCanvas can_2d;
   for (auto const &sPresel_vVars:msPresel_vVars2D){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         can_2d.cd();
         can_2d.SetLogz();
         gPad->SetRightMargin(0.2);
         gPad->SetLeftMargin(0.13);
         gPad->SetBottomMargin(0.10);
         TString loc=sPresel+sVar;
         TH2F *hist=hs2d.getHistogram(loc,"TTbar");
         
         hist->GetYaxis()->SetTitleOffset(1.3);
         hist->GetXaxis()->SetTitleOffset(0.9);
         hist->GetZaxis()->SetTitleOffset(1.3);
         hist->GetYaxis()->SetTitleSize(0.05);
         hist->GetXaxis()->SetTitleSize(0.05);
         hist->GetZaxis()->SetTitleSize(0.05);
         hist->GetYaxis()->SetLabelSize(0.04);
         hist->GetXaxis()->SetLabelSize(0.04);
         hist->GetZaxis()->SetLabelSize(0.04);
                  
         hist->SetStats(0);
         hist->Draw("colz");
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         TLatex label=gfx::cornerLabel(cat,2);
         label.Draw();
         saver.save(can_2d,"tt_only/"+loc);
      }
   }
   
   //Save 2d histograms
   saveHistograms2D(msPresel_vVars2D,saver_hist,hs2d,samplesToCombine);
}
