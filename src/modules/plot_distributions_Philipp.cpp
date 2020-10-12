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
   // ~io::RootFileSaver saver(TString::Format("plots%.1f._oldTriggSelec.root",cfg.processFraction*100),"plot_distributions_Philipp");
   // ~io::RootFileReader histReader(TString::Format("histograms_%s_oldTriggSelec.root",cfg.treeVersion.Data()),TString::Format("distributions_Philipp%.1f",cfg.processFraction*100));
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions_Philipp");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions_Philipp%.1f",cfg.processFraction*100));
   
   std::vector<TString> samplesToPlot=cfg.datasets.getDatasetNames();
   std::reverse(samplesToPlot.begin(),samplesToPlot.end());
   // ~std::vector<TString> samplesToPlot={"DrellYan_NLO","DoubleMuon"};
   
   // Define distributions to be plotted
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"Zpeak/ee/",{"met","PuppiMet","dphi_metNearLep","genHT","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak/mumu/",{"met","PuppiMet","dphi_metNearLep","genHT","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ/ee/",{"met","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","PuppiMet","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
      {"Zpeak_noJetRequ/mumu/",{"met","pTl1","pTl2","etal1","etal2","phil1","phil2","ZpT","mll","dPhill","nJets","looseJetpT1","nGoodVertices","nPV","PuppiMet","CaloMet","DeepMet","RawMet","NoHFMet","PFminiIsoL1","PFminiIsoL2","nLooseJets","ChargeL1","ChargeL2","ChargeProduct"}},
   };
   
   hist::Histograms<TH1F> hs(samplesToPlot);
   
   //Read histogram into HS container
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            TString loc;
            loc=sPresel+sVar;
            TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
            // ~if (sVar!="nBjets" and sVar!="nJets") tempHist->Rebin(5);
            if(sSample=="DrellYan_NLO")tempHist->Scale(5931.9/6225.4);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
      }
   }
   
   //Plot stacked distributions
   gfx::SplitCan can;
   can.pU_.SetLogy();
   for (auto const &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         can.pU_.cd();
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs.getStack(loc,samplesToPlot);
         // ~TH1F* MC=hs.getHistogram(loc,"DrellYan_NLO");
         // ~TH1F* MC=hs.getHistogram(loc,"DrellYan");
         gfx::LegendEntries le=hs.getLegendEntries();
         TH1F* data;
         TH1F* data_single;
         TString cat;
         if (sPresel.Contains("ee/")) {
            cat="ee";
            data=hs.getHistogram(loc,"DoubleEG");
            data_single=hs.getHistogram(loc,"SingleElectron");
            }
         else if (sPresel.Contains("mumu/")) {
            cat="#mu#mu";
            data=hs.getHistogram(loc,"DoubleMuon");
            data_single=hs.getHistogram(loc,"SingleMuon");
         }
         data->Add(data_single);
         le.prepend(*data,"data","pe1");
         TLatex label=gfx::cornerLabel(cat,1);
         TString selection="ZPeak";
         if(sPresel.Contains("noJetRequ")) selection+=", no jet cuts";
         TLatex label2=gfx::cornerLabel(selection,3);
         // ~if (sVar.Contains("phi")){
            // ~st_mc.SetMinimum(0.01);
            // ~st_mc.SetMaximum(1e6);
         // ~}
         data->SetMinimum(0.11);
         data->SetMaximum(5*st_mc.GetMaximum());
         data->SetStats(0);
         data->GetYaxis()->SetTitleOffset(0.85);
         data->Draw("axis");
         st_mc.Draw("same");
         data->Draw("pe1 same");
         
         TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         leg.Draw();
         label.Draw();
         label2.Draw();
         can.pL_.cd();
         TH1F ratio=hist::getRatio(*data,st_mc,"data/MC",hist::ONLY1);
         TH1F ratio_mc=hist::getRatio(st_mc,st_mc,"data/MC",hist::ONLY1);
         ratio_mc.GetYaxis()->SetTitleOffset(0.3);
         ratio_mc.SetStats(0);
         ratio.SetLineColor(kBlack);
         ratio_mc.SetMaximum(1.25);
         ratio_mc.SetMinimum(0.75);
         // ~ratio_mc.SetMaximum(1.4);
         // ~ratio_mc.SetMinimum(0.6);
         ratio_mc.SetFillColor(kGray);
         ratio_mc.SetMarkerSize(0);
         ratio_mc.Draw("e2");
         ratio.Draw("pe1 same");
         if(std::find(samplesToPlot.begin(),samplesToPlot.end(),"DrellYan_NLO") != samplesToPlot.end()){
            saver.save(can,"DYnlo/"+loc,false);
         }
         else saver.save(can,loc,false);
         
      }
   }
}
