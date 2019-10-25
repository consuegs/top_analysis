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
   std::vector<TString> vsDatasubsets({"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_merged"});
   
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   hs.addHist("baseline/ee/met"   ,";%MET;EventsBIN"           ,100,0,600);
   
   hist::Histograms<TH2F> hs2d(vsDatasubsets);     //Define 2D histograms in the following
   hs2d.addHist("genTopVSpseudoTop/bQuarks_pT", ";p_{T}^{b} (genTop) (GeV);p_{T}^{b} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genTopVSpseudoTop/leptons_pMag", ";|p^{l}| (genTop) (GeV);|p^{l}| (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genTopVSpseudoTop/leptons_pT", ";p_{T}^{l} (genTop) (GeV);p_{T}^{l} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genTopVSpseudoTop/leptons_phi", ";#phi^{l} (genTop) ;phi^{l} (pseudo);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genTopVSpseudoTop/neutrinos_pT", ";p_{T}^{#nu} (genTop) (GeV);p_{T}^{#nu} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genTopVSpseudoTop/neutrinos_phi", ";#phi^{#nu} (genTop) ;phi^{#nu} (pseudo);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genTopVSpseudoTop/sumNeutrinos_pT", ";p_{T}^{#nu#nu} (genTop) (GeV);p_{T}^{#nu#nu} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   
   hs2d.addHist("genVSpseudoTop/leptons_pT", ";p_{T}^{l} (gen) (GeV);p_{T}^{l} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSpseudoTop/leptons_pMag", ";|p^{l}| (gen) (GeV);|p^{l}| (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSpseudoTop/leptons_phi", ";#phi^{l} (gen) ;phi^{l} (pseudo);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genVSpseudoTop/neutrinos_pT", ";p_{T}^{#nu} (gen) (GeV);p_{T}^{#nu} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSpseudoTop/neutrinos_phi", ";#phi^{#nu} (gen) ;phi^{#nu} (pseudo);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genVSpseudoTop/sumNeutrinos_pT", ";p_{T}^{#nu#nu} (gen) (GeV);p_{T}^{#nu#nu} (pseudo)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   
   hs2d.addHist("genVSgenTop/leptons_pT", ";p_{T}^{l} (gen) (GeV);p_{T}^{l} (genTop)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSgenTop/leptons_pMag", ";|p^{l}| (gen) (GeV);|p^{l}| (genTop)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSgenTop/leptons_phi", ";#phi^{l} (gen) ;phi^{l} (genTop);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genVSgenTop/neutrinos_pT", ";p_{T}^{#nu} (gen) (GeV);p_{T}^{#nu} (genTop)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   hs2d.addHist("genVSgenTop/neutrinos_phi", ";#phi^{#nu} (gen) ;phi^{#nu} (genTop);Particles/BIN" ,100,-3.14,3.14,100,-3.14,-3.14);
   hs2d.addHist("genVSgenTop/sumNeutrinos_pT", ";p_{T}^{#nu#nu} (gen) (GeV);p_{T}^{#nu#nu} (genTop)(GeV);Particles/BIN" ,100,0,400,100,0,400);
   

   auto const dss = cfg.datasets.getDatasubset(vsDatasubsets[0]);
   TFile file(dss.getPath(),"read");
   if (file.IsZombie()) {
      return;
   }
   io::log * ("Processing '"+dss.datasetName+"' ");
   
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
   TTreeReaderValue<float> mt2   (reader, "MT2");
   TTreeReaderValue<float> genMT2   (reader, "genMT2");
   TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
   TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
   TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
   TTreeReaderValue<TLorentzVector> pseudoTop(reader, "pseudoTop");
   TTreeReaderValue<TLorentzVector> pseudoAntiTop(reader, "pseudoAntiTop");
   TTreeReaderValue<TLorentzVector> pseudoLepton(reader, "pseudoLepton");
   TTreeReaderValue<TLorentzVector> pseudoAntiLepton(reader, "pseudoAntiLepton");
   TTreeReaderValue<TLorentzVector> pseudoTau(reader, "pseudoTau");
   TTreeReaderValue<TLorentzVector> pseudoAntiTau(reader, "pseudoAntiTau");
   TTreeReaderValue<int> pseudoLeptonPdgId(reader, "pseudoLeptonPdgId");
   TTreeReaderValue<int> pseudoAntiLeptonPdgId(reader, "pseudoAntiLeptonPdgId");
   TTreeReaderValue<TLorentzVector> pseudoB(reader, "pseudoBJet");
   TTreeReaderValue<TLorentzVector> pseudoAntiB(reader, "pseudoAntiBJet");
   TTreeReaderValue<TLorentzVector> pseudoNeutrino(reader, "pseudoNeutrino");
   TTreeReaderValue<TLorentzVector> pseudoAntiNeutrino(reader, "pseudoAntiNeutrino");
   TTreeReaderValue<TLorentzVector> pseudoWMinus(reader, "pseudoWMinus");
   TTreeReaderValue<TLorentzVector> pseudoWPlus(reader, "pseudoWPlus");
   TTreeReaderValue<int> pseudoDecayMode(reader, "ttbarPseudoDecayMode");
   TTreeReaderValue<TLorentzVector> genTop(reader, "genTop");    //Alternative genParticles, which are not based on RIVET
   TTreeReaderValue<TLorentzVector> genAntiTop(reader, "genAntiTop");
   TTreeReaderValue<TLorentzVector> genLepton(reader, "genLepton");
   TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "genAntiLepton");
   TTreeReaderValue<TLorentzVector> genTau(reader, "genTau");
   TTreeReaderValue<TLorentzVector> genAntiTau(reader, "genAntiTau");
   TTreeReaderValue<int> genLeptonPdgId(reader, "genLeptonPdgId");
   TTreeReaderValue<int> genAntiLeptonPdgId(reader, "genAntiLeptonPdgId");
   TTreeReaderValue<TLorentzVector> genB(reader, "genB");
   TTreeReaderValue<TLorentzVector> genAntiB(reader, "genAntiB");
   TTreeReaderValue<TLorentzVector> genNeutrino(reader, "genNeutrino");
   TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "genAntiNeutrino");
   TTreeReaderValue<TLorentzVector> genWMinus(reader, "genWMinus");
   TTreeReaderValue<TLorentzVector> genWPlus(reader, "genWPlus");
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
      
      // ~float fEventWeight=*w_pu * *w_mc * *sf_lep1 * *sf_lep2;     //Set event weight also taking lepton scale factors into account
      hs.setFillWeight(1.);
      hs2d.setFillWeight(1.);
      
      float const met=MET->p.Pt();
      
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
      
      if (*pseudoDecayMode==0) pseudo_selection=false; //pseudo baseline selection
      
      if(pseudo_selection==false) continue;  // only proceed with events selected by pseudo baseline selection
      
      if(*genDecayMode>3) continue;  // decays including taus
      if(genNeutrino->Pt()<1 || genAntiNeutrino->Pt()<1) continue;  // neutrinos not found
      
      TLorentzVector pureNeutrino(0,0,0,0);
      TLorentzVector pureAntiNeutrino(0,0,0,0);
      TLorentzVector pureLepton(0,0,0,0);
      TLorentzVector pureAntiLepton(0,0,0,0);
      for (auto const &genParticle : *genParticles){
         if ((genParticle.pdgId==12 || genParticle.pdgId==14) && pureNeutrino.Pt()==0){
            pureNeutrino=genParticle.p;
         }
         if ((genParticle.pdgId==-12 || genParticle.pdgId==-14) && pureAntiNeutrino.Pt()==0){
            pureAntiNeutrino=genParticle.p;
         }
         if ((genParticle.pdgId==11 || genParticle.pdgId==13) && pureLepton.Pt()==0){
            pureLepton=genParticle.p;
         }
         if ((genParticle.pdgId==-11 || genParticle.pdgId==-13) && pureAntiLepton.Pt()==0){
            pureAntiLepton=genParticle.p;
         }
      }
      
      hs2d.fill("genTopVSpseudoTop/bQuarks_pT",genB->Pt(),pseudoB->Pt());
      hs2d.fill("genTopVSpseudoTop/bQuarks_pT",genAntiB->Pt(),pseudoAntiB->Pt());
      hs2d.fill("genTopVSpseudoTop/leptons_pT",genLepton->Pt(),pseudoLepton->Pt());
      hs2d.fill("genTopVSpseudoTop/leptons_pT",genAntiLepton->Pt(),pseudoAntiLepton->Pt());
      hs2d.fill("genTopVSpseudoTop/leptons_pMag",genLepton->P(),pseudoLepton->P());
      hs2d.fill("genTopVSpseudoTop/leptons_pMag",genAntiLepton->P(),pseudoAntiLepton->P());
      hs2d.fill("genTopVSpseudoTop/leptons_phi",genLepton->Phi(),pseudoLepton->Phi());
      hs2d.fill("genTopVSpseudoTop/leptons_phi",genAntiLepton->Phi(),pseudoAntiLepton->Phi());
      hs2d.fill("genTopVSpseudoTop/neutrinos_pT",genNeutrino->Pt(),pseudoNeutrino->Pt());
      hs2d.fill("genTopVSpseudoTop/neutrinos_pT",genAntiNeutrino->Pt(),pseudoAntiNeutrino->Pt());
      hs2d.fill("genTopVSpseudoTop/neutrinos_phi",genNeutrino->Phi(),pseudoNeutrino->Phi());
      hs2d.fill("genTopVSpseudoTop/neutrinos_phi",genAntiNeutrino->Phi(),pseudoAntiNeutrino->Phi());
      hs2d.fill("genTopVSpseudoTop/sumNeutrinos_pT",(*(genAntiNeutrino)+*(genNeutrino)).Pt(),(*(pseudoAntiNeutrino)+*(pseudoNeutrino)).Pt());
      
      if(abs(pureLepton.Phi()-pseudoLepton->Phi())>0.05) continue;
      if(abs(pureAntiLepton.Phi()-pseudoAntiLepton->Phi())>0.05) continue;
      if(abs(pureLepton.Phi()-genLepton->Phi())>0.05) continue;
      if(abs(pureAntiLepton.Phi()-genAntiLepton->Phi())>0.05) continue;
      
      hs2d.fill("genVSpseudoTop/leptons_pT",pureLepton.Pt(),pseudoLepton->Pt());
      hs2d.fill("genVSpseudoTop/leptons_pT",pureAntiLepton.Pt(),pseudoAntiLepton->Pt());
      hs2d.fill("genVSpseudoTop/leptons_pMag",pureLepton.P(),pseudoLepton->P());
      hs2d.fill("genVSpseudoTop/leptons_pMag",pureAntiLepton.P(),pseudoAntiLepton->P());
      hs2d.fill("genVSpseudoTop/leptons_phi",pureLepton.Phi(),pseudoLepton->Phi());
      hs2d.fill("genVSpseudoTop/leptons_phi",pureAntiLepton.Phi(),pseudoAntiLepton->Phi());
      hs2d.fill("genVSpseudoTop/neutrinos_pT",pureNeutrino.Pt(),pseudoNeutrino->Pt());
      hs2d.fill("genVSpseudoTop/neutrinos_pT",pureAntiNeutrino.Pt(),pseudoAntiNeutrino->Pt());
      hs2d.fill("genVSpseudoTop/neutrinos_phi",pureNeutrino.Phi(),pseudoNeutrino->Phi());
      hs2d.fill("genVSpseudoTop/neutrinos_phi",pureAntiNeutrino.Phi(),pseudoAntiNeutrino->Phi());
      hs2d.fill("genVSpseudoTop/sumNeutrinos_pT",(pureAntiNeutrino+pureNeutrino).Pt(),(*(pseudoAntiNeutrino)+*(pseudoNeutrino)).Pt());
      
      hs2d.fill("genVSgenTop/leptons_pT",pureLepton.Pt(),genLepton->Pt());
      hs2d.fill("genVSgenTop/leptons_pT",pureAntiLepton.Pt(),genAntiLepton->Pt());
      hs2d.fill("genVSgenTop/leptons_pMag",pureLepton.P(),genLepton->P());
      hs2d.fill("genVSgenTop/leptons_pMag",pureAntiLepton.P(),genAntiLepton->P());
      hs2d.fill("genVSgenTop/leptons_phi",pureLepton.Phi(),genLepton->Phi());
      hs2d.fill("genVSgenTop/leptons_phi",pureAntiLepton.Phi(),genAntiLepton->Phi());
      hs2d.fill("genVSgenTop/neutrinos_pT",pureNeutrino.Pt(),genNeutrino->Pt());
      hs2d.fill("genVSgenTop/neutrinos_pT",pureAntiNeutrino.Pt(),genAntiNeutrino->Pt());
      hs2d.fill("genVSgenTop/neutrinos_phi",pureNeutrino.Phi(),genNeutrino->Phi());
      hs2d.fill("genVSgenTop/neutrinos_phi",pureAntiNeutrino.Phi(),genAntiNeutrino->Phi());
      hs2d.fill("genVSgenTop/sumNeutrinos_pT",(pureAntiNeutrino+pureNeutrino).Pt(),(*(genAntiNeutrino)+*(genNeutrino)).Pt());
      
   }// evt loop
   io::log<<"";
   
   // ~hs.scaleLumi();
   // ~hs2d.scaleLumi();
   hs.mergeOverflow();
   hs2d.mergeOverflow();
   file.Close();
   
   std::vector<TString> samplesToCombine={"TTbar"};
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   
   
   //Plotting part 2D
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("compare_pseudoGen%.1f",cfg.processFraction*100));
   std::map<TString,std::vector<TString>> msPresel_vVars2D={
      {"genTopVSpseudoTop/",{"bQuarks_pT","leptons_pT","leptons_pMag","leptons_phi","neutrinos_pT","neutrinos_phi","sumNeutrinos_pT"}},
      {"genVSpseudoTop/",{"leptons_pT","leptons_pMag","leptons_phi","neutrinos_pT","neutrinos_phi","sumNeutrinos_pT"}},
      {"genVSgenTop/",{"leptons_pT","leptons_pMag","leptons_phi","neutrinos_pT","neutrinos_phi","sumNeutrinos_pT"}},
      };
   
   // The following can be used for plotting, which is currently done with an extra module 
   TCanvas can_2d;
   for (auto const &sPresel_vVars:msPresel_vVars2D){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         can_2d.cd();
         // ~can_2d.SetLogz();
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
         saver.save(can_2d,loc);
      }
   }
   
   // Save histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("compare_pseudoGen%.1f",cfg.processFraction*100),false);   
   saveHistograms2D(msPresel_vVars2D,saver_hist,hs2d,samplesToCombine);
}
