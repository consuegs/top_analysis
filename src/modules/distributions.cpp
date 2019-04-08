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
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hs.addHist("baseline/ee/met"   ,";%MET;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/met"   ,";%MET;EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/met"   ,";%MET;EventsBIN"           ,100,0,600);
   
   hs.addHist("baseline/ee/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("baseline/emu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("baseline/mumu/mll"   ,";mll(GeV);EventsBIN"           ,100,0,600);
   
   hs.addHist("genParticles/ee/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/emu/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/mumu/pT_nunu"   ,";p_{T}^{#nu#nu}(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/ee/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/emu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   hs.addHist("genParticles/mumu/genMet"   ,";genMET(GeV);EventsBIN"           ,100,0,600);
   
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   hs2d.addHist("genParticles/ee/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/emu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);
   hs2d.addHist("genParticles/mumu/2d_nunuVSgenMet", ";genMET (GeV);p_{T}^{#nu#nu}(GeV);EventsBIN" ,100,0,600,100,0,600);


   
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
         for (tree::Jet const &jet : cjets) {
            if (jet.bTagCSVv2>0.5426) bTag=true;
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
         
         //Fill hists
         if (*is_ee){
            hs.fill("baseline/ee/met",met);
            hs.fill("baseline/ee/mll",*mll);
            hs.fill("genParticles/ee/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/ee/genMet",genMet);
            hs2d.fill("genParticles/ee/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
         }
         else if (*is_emu){
            hs.fill("baseline/emu/met",met);
            hs.fill("baseline/emu/mll",*mll);
            hs.fill("genParticles/emu/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/emu/genMet",genMet);
            hs2d.fill("genParticles/emu/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
         }
         else if (*is_mumu){
            hs.fill("baseline/mumu/met",met);
            hs.fill("baseline/mumu/mll",*mll);
            hs.fill("genParticles/mumu/pT_nunu",neutrinoPair.Pt());
            hs.fill("genParticles/mumu/genMet",genMet);
            hs2d.fill("genParticles/mumu/2d_nunuVSgenMet",genMet,neutrinoPair.Pt());
         }
         
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs2d.scaleLumi();
      hs.mergeOverflow();
      hs2d.mergeOverflow();
      file.Close();
   } // dataset loop
   
   std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ"};
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   //Plotting part 1D
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","mll"}},
      {"baseline/emu/",{"met","mll"}},
      {"baseline/mumu/",{"met","mll"}},
      {"genParticles/ee/",{"pT_nunu","genMet"}},
      {"genParticles/emu/",{"pT_nunu","genMet"}},
      {"genParticles/mumu/",{"pT_nunu","genMet"}},
      };
      
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         TH1F* hist=hs.getHistogram(sVar,"TTbar");
         gfx::LegendEntries le=hs.getLegendEntries();
         hist->SetStats(0);
         hist->SetMarkerSize(0);
         hist->Draw("histE");
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         saver.save(can,"tt_only/"+sVar);
         
         THStack st_mc=hs.getStack(sVar,{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar"});
         le=hs.getLegendEntries();
         st_mc.Draw();
         TLegend leg2=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg2.Draw();
         saver.save(can,"all/"+sVar);
      }
   }
   
   //Plotting part 2D
   TCanvas can_2d;
   for (TString selection:{"ee","emu","mumu"}){
      can_2d.cd();
      TString loc="genParticles/"+selection+"/2d_nunuVSgenMet";
      TH2F *hist=hs2d.getHistogram(loc,"TTbar");
      hist->SetStats(0);
      hist->Draw("colz");
      saver.save(can_2d,"tt_only/"+loc);
   }
}
