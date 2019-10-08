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
      "T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","TTbar_diLepton","ttZ_SM"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/emu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/mumu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/ee/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/emu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline_Met200/mumu/",{"met","met1000","mll","pTlep1","pTlep2","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"genParticles/ee/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles/emu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles/mumu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/ee/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/emu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"genParticles_Met200/mumu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
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
         if (sPresel.Contains("Met200")) cat+="  p_{T}^{miss}>200 GeV";
         TLatex label=gfx::cornerLabel(cat,1);
         if (sVar.Contains("phi")){
            st_mc.SetMinimum(0.01);
            st_mc.SetMaximum(1e6);
         }
         st_mc.SetMinimum(0.01);
         st_mc.SetMaximum(1e6);
         st_mc.Draw();
         
         auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","ttZ_SM"});
         for (auto const &h: hists) h->Draw("same hist");
         le+=hs.getLegendEntries();
         TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         TH1F axis=*(hists[0]);
         axis.SetStats(0);
         axis.Draw("same axis");
         leg.Draw();
         label.Draw();
         saver.save(can,loc);
      }
   }
   
   for (TString cat:{"ee","emu","mumu"}){    //Get the number of events per category
      TH1F* temp_hist=hs.getHistogram("baseline/"+cat+"/met","TTbar");
      //~ TH1F* temp_hist=hs.getHistogram("baseline/"+cat+"/met","TTbar_diLepton");
      std::cout<<cat<<"   "<<temp_hist->Integral()<<std::endl;
   }
}
