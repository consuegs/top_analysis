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
   
   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.datasetName+"' ");

      bool const isData=dss.isData;
      bool const isSignal=dss.isSignal;
      
      hs.setCurrentSample(dss.name);

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
         
         float const met=MET->p.Pt();
         
         //Baseline selection
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         
         if (*is_ee){
            if(!(*electrons)[0].isTight || !(*electrons)[1].isTight) continue;
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
         for (tree::Jet jet : cjets) {
            if (jet.bTagCSVv2>0.5426) bTag=true;
         }
         if (!bTag) continue;
         
         
         //Fill hists
         if (*is_ee){
            hs.fill("baseline/ee/met",met);
            hs.fill("baseline/ee/mll",*mll);
         }
         else if (*is_emu){
            hs.fill("baseline/emu/met",met);
            hs.fill("baseline/emu/mll",*mll);
         }
         else if (*is_mumu){
            hs.fill("baseline/mumu/met",met);
            hs.fill("baseline/mumu/mll",*mll);
         }
         
      }// evt loop
      io::log<<"";
      
      hs.scaleLumi();
      hs.mergeOverflow();
      file.Close();
   } // dataset loop
   
   std::vector<TString> samplesToCombine={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ"};
   hs.combineFromSubsamples(samplesToCombine);
   
   //Plotting part
   io::RootFileSaver saver("plots.root",TString::Format("distributions%.1f",cfg.processFraction*100));
   
   TCanvas can;
   can.SetLogy();
   // what to plot in which preselection
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","mll"}},
      {"baseline/emu/",{"met","mll"}},
      {"baseline/mumu/",{"met","mll"}},};
      
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sVar=sPresel+sVar;
         THStack st_mc=hs.getStack(sVar,{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","TTbar"});
         gfx::LegendEntries le=hs.getLegendEntries();
         st_mc.Draw();
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         saver.save(can,sVar);
      }
   }
}
