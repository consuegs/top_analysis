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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_ttMC");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   // ~std::vector<TString> samplesToImport={"TTbar","TTbar_diLepton","TTbar_singleLepton","TTbar_madGraphCOMB"};
   std::vector<TString> samplesToImport={"TTbar","TTbar_diLepton","TTbar_madGraphCOMB"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met1000","pTlep2"}},
      {"baseline/emu/",{"met1000","pTlep2"}},
      {"baseline/mumu/",{"met1000","pTlep2",}},
   };
   
   hist::Histograms<TH1F> hs(samplesToImport);
   
   for (TString sSample : samplesToImport){
      hs.setCurrentSample(sSample);
      
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            TString loc;
            loc=sPresel+sVar;
            TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
            // ~tempHist->Rebin(5);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
      }
   }
   
   // ~hs.combineSamples("ll_sl_comb",{"TTbar_diLepton","TTbar_singleLepton"});
   
   // ~std::vector<TString> samplesToPlot={"TTbar","TTbar_diLepton","ll_sl_comb","TTbar_madGraphCOMB"};
   std::vector<TString> samplesToPlot={"TTbar","TTbar_diLepton","TTbar_madGraphCOMB"};
    
   gfx::SplitCan spcan;
   spcan.cdUp();
   spcan.pU_.SetLogy();
   for (auto const &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc;
         loc=sPresel+sVar;
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         TLatex label=gfx::cornerLabel(cat,1);
         
         auto hist_tt=hs.getHistogram(loc,"TTbar");
         auto hist_ttDilep=hs.getHistogram(loc,"TTbar_diLepton");
         // ~auto hist_ttComb=hs.getHistogram(loc,"ll_sl_comb");
         auto hist_ttMad=hs.getHistogram(loc,"TTbar_madGraphCOMB");
         hist_tt->SetStats(0);
         hist_tt->Draw("axis");
         // ~hist_ttComb->SetLineColor(kCyan);
         hist_ttMad->SetLineColor(kGreen);
         // ~for (auto const &h: {hist_tt,hist_ttDilep,hist_ttComb,hist_ttMad}) {
         for (auto const &h: {hist_tt,hist_ttDilep,hist_ttMad}) {
            h->SetStats(0);
            h->SetMarkerSize(0);
            h->Draw("same e");
         }
         gfx::LegendEntries le;
         le.append(*hist_tt,"Powheg incl.","l");
         le.append(*hist_ttDilep,"Powheg dilep.","l");
         // ~le.append(*hist_ttComb,"Powheg dilep.+slep","l");
         le.append(*hist_ttMad,"Madgraph dilep.","l");
         TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
         leg.Draw();
         label.Draw();
         
         spcan.cdLow();
         TH1F hRatio=hist::getRatio(*hist_tt,*hist_tt,"Ratio.",hist::ONLY1);
         TH1F hRatio_Dilep=hist::getRatio(*hist_ttDilep,*hist_tt,"Ratio.",hist::ONLY1);
         // ~TH1F hRatio_comb=hist::getRatio(*hist_ttComb,*hist_tt,"Ratio.",hist::ONLY1);
         TH1F hRatio_Mad=hist::getRatio(*hist_ttMad,*hist_tt,"Ratio.",hist::ONLY1);
         hRatio.SetStats(false);
         hRatio.SetMarkerSize(0);
         hRatio_Dilep.SetMarkerSize(0);
         // ~hRatio_comb.SetMarkerSize(0);
         // ~hRatio_comb.SetLineColor(hist_ttComb->GetLineColor());
         hRatio_Mad.SetMarkerSize(0);
         hRatio.Draw("axis e0");
         hRatio.Draw("same pe0");
         hRatio_Dilep.Draw("same pe0");
         // ~hRatio_comb.Draw("same pe0");
         hRatio_Mad.Draw("same pe0");
         saver.save(spcan,loc);
      }
   }
}
