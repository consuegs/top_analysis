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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions_Philipp");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions_Philipp%.1f",cfg.processFraction*100));
   
   std::vector<TString> samplesToPlot=cfg.datasets.getDatasetNames();
   
   // Define distributions to be plotted
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","dphi_metNearLep"}},
      {"baseline/emu/",{"met","dphi_metNearLep"}},
      {"baseline/mumu/",{"met","dphi_metNearLep"}},
   };
   
   hist::Histograms<TH1F> hs(samplesToPlot);
   
   //Perform rebinning of histograms
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            TString loc;
            loc=sPresel+sVar;
            TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
            if (sVar!="nBjets" and sVar!="nJets") tempHist->Rebin(5);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
      }
   }
   
   
   //Plot stacked distributions
   TCanvas can;
   can.SetLogy();
   for (auto const &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs.getStack(loc,samplesToPlot);
         gfx::LegendEntries le=hs.getLegendEntries();
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         if (sPresel.Contains("Met200")) cat+="  p_{T}^{miss}>200 GeV";
         if (sPresel.Contains("Met400")) cat+="  p_{T}^{miss}>400 GeV";
         TLatex label=gfx::cornerLabel(cat,1);
         if (sVar.Contains("phi")){
            st_mc.SetMinimum(0.01);
            st_mc.SetMaximum(1e6);
         }
         st_mc.SetMinimum(0.01);
         st_mc.SetMaximum(1e6);
         st_mc.Draw();
         
         TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         leg.Draw();
         label.Draw();
         saver.save(can,loc);
         
      }
   }
}
