//Script to (re-)plot distributions from distributions.cpp 

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   std::vector<TString> samplesToPlot={"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      "T1tttt_1200_800","T2tt_650_250","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
      {"baseline/emu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
      {"baseline/mumu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets"}},
   };
   
   hist::Histograms<TH1F> hs(samplesToPlot);
   
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            TString loc;
            loc=sPresel+sVar;
            TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
            if (sVar!="nBjets") tempHist->Rebin(4);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
      }
   }
   
   hs.combineSamples("SM bkg",{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ"});
   
    
   TCanvas can;
   can.SetLogy();
   for (auto const &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs.getStack(loc,{"SM bkg","TTbar"});
         gfx::LegendEntries le=hs.getLegendEntries();
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         TLatex label=gfx::cornerLabel(cat,2);
         st_mc.Draw();
         
         auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T2tt_650_250","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"});
         for (auto const &h: hists) h->Draw("same hist");
         le+=hs.getLegendEntries();
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         label.Draw();
         saver.save(can,loc);
      }
   }
}
