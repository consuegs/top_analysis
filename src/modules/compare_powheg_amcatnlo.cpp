//Script to compare DNN inputs for POWHEG and amcatnlo 

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
#include <TObjString.h>

Config const &cfg=Config::get();

void compare_tt_powheg_madgraph(){
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_ttMC");
   io::RootFileReader histReader(TString::Format("multiHists/Nominal/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   std::vector<TString> samplesToImport={"TTbar_diLepton","TTbar_amcatnlo"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={};
   for(TString channel:{"ee/","emu/","mumu/","combined/"}){
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("baseline/"+channel,
      {"PuppiMET","METunc_Puppi","MET","HT","nJets","n_Interactions","Lep1_flavor","Lep2_flavor","Lep1_pt","Lep1_phi","Lep1_eta","Lep1_E","Lep2_pt","Lep2_phi","Lep2_eta","Lep2_E","Jet1_pt","Jet1_phi","Jet1_eta","Jet1_E","Jet2_pt","Jet2_phi","Jet2_eta","Jet2_E","dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","dPhiJet1Jet2","METsig","MHT","MT","looseLeptonVeto","dPhiMETnearJet_Puppi","dPhiMETfarJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiMETbJet_Puppi","dPhiLep1bJet","dPhiLep1Jet1","mLL","PFMET_phi","PuppiMET_phi","CaloMET","CaloMET_phi","MT2","vecsum_pT_allJet","vecsum_pT_l1l2_allJet","mass_l1l2_allJet","ratio_vecsumpTlep_vecsumpTjet","mjj"}));
      // ~msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("genParticles/"+channel,{"genMet"}));
   }
   
   hist::Histograms<TH1F> hs(samplesToImport);
   
   for (TString sSample : samplesToImport){
      hs.setCurrentSample(sSample);
      
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            if(sPresel.Contains("combined/")){
               TString loc;
               loc=sPresel+sVar;
               TString sel=((TObjString*)sPresel.Tokenize("/")->First())->String();
               // ~TH1F* tempHist=histReader.read<TH1F>("baseline/ee/"+sVar+"/"+sSample);
               // ~tempHist->Add(histReader.read<TH1F>("baseline/emu/"+sVar+"/"+sSample));
               // ~tempHist->Add(histReader.read<TH1F>("baseline/mumu/"+sVar+"/"+sSample));
               TH1F* tempHist=histReader.read<TH1F>(sel+"/ee/"+sVar+"/"+sSample);
               tempHist->Add(histReader.read<TH1F>(sel+"/emu/"+sVar+"/"+sSample));
               tempHist->Add(histReader.read<TH1F>(sel+"/mumu/"+sVar+"/"+sSample));
               hs.addFilledHist(loc,sSample,*(tempHist));
            }
            else{
               TString loc;
               loc=sPresel+sVar;
               TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
               hs.addFilledHist(loc,sSample,*(tempHist));
            }
         }
      }
      hs.normHists();
   }
    
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
         else if (sPresel.Contains("combined/")) cat="all";
         TLatex label=gfx::cornerLabel(cat,1);
         
         auto hist_ttDilep=hs.getHistogram(loc,"TTbar_diLepton");
         auto hist_ttAmc=hs.getHistogram(loc,"TTbar_amcatnlo");
         hist_ttDilep->SetStats(0);
         hist_ttDilep->Draw("axis");
         hist_ttAmc->SetLineColor(kBlue);
         for (auto const &h: {hist_ttDilep,hist_ttAmc}) {
            // ~if (hist_ttDilep->GetNbinsX()>20) h->Rebin(2);
            h->SetStats(0);
            h->SetMarkerSize(0);
            h->Draw("same e");
         }
         gfx::LegendEntries le;
         le.append(*hist_ttDilep,"POWHEG dilep.","l");
         le.append(*hist_ttAmc,"amc@NLO dilep.","l");
         TLegend leg=le.buildLegend(.6,.8,1-gPad->GetRightMargin(),-1,1);
         leg.Draw();
         label.Draw();
         
         spcan.cdLow();
         TH1F hRatio=hist::getRatio(*hist_ttDilep,*hist_ttDilep,"Ratio.",hist::ONLY1);
         TH1F hRatio_amc=hist::getRatio(*hist_ttAmc,*hist_ttDilep,"Ratio.",hist::ONLY1);
         hRatio.SetStats(false);
         hRatio.SetMarkerSize(0);
         hRatio.SetMaximum(1.25);
         hRatio.SetMinimum(0.75);
         hRatio_amc.SetMarkerSize(0);
         hRatio.Draw("axis e0");
         hRatio.Draw("same pe0");
         hRatio_amc.Draw("same pe0");
         saver.save(spcan,loc,true,true);
      }
   }
}

void compare_tW_ds_dr(const bool add_tt=false){
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_tW_ds_dr");
   
   std::vector<TString> samplesToImport={"SingleTop","SingleTop_TWDS"};
   if(add_tt) samplesToImport.push_back("TTbar_diLepton");
   
   std::map<TString,std::vector<TString>> msPresel_vVars={};
   for(TString channel:{"ee/","emu/","mumu/","combined/"}){
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("baseline/"+channel,
      {"PuppiMET","METunc_Puppi","MET","HT","nJets","n_Interactions","Lep1_flavor","Lep2_flavor","Lep1_pt","Lep1_phi","Lep1_eta","Lep1_E","Lep2_pt","Lep2_phi","Lep2_eta","Lep2_E","Jet1_pt","Jet1_phi","Jet1_eta","Jet1_E","Jet2_pt","Jet2_phi","Jet2_eta","Jet2_E","dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","dPhiJet1Jet2","METsig","MHT","MT","looseLeptonVeto","dPhiMETnearJet_Puppi","dPhiMETfarJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiMETbJet_Puppi","dPhiLep1bJet","dPhiLep1Jet1","mLL","PFMET_phi","PuppiMET_phi","CaloMET","CaloMET_phi","MT2","vecsum_pT_allJet","vecsum_pT_l1l2_allJet","mass_l1l2_allJet","ratio_vecsumpTlep_vecsumpTjet","mjj","DNN_MET_pT"}));
      // ~msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("genParticles/"+channel,{"genMet"}));
   }
   
   hist::Histograms<TH1F> hs(samplesToImport);
   
   for (TString sSample : samplesToImport){
      hs.setCurrentSample(sSample);
      
      io::RootFileReader histReader(TString::Format("multiHists/%s/histograms_merged_%s.root",(sSample != "SingleTop_TWDS")? "Nominal":"TWDS",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
            
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            if(sPresel.Contains("combined/")){
               TString loc;
               loc=sPresel+sVar;
               TString sel=((TObjString*)sPresel.Tokenize("/")->First())->String();
               // ~TH1F* tempHist=histReader.read<TH1F>("baseline/ee/"+sVar+"/"+sSample);
               // ~tempHist->Add(histReader.read<TH1F>("baseline/emu/"+sVar+"/"+sSample));
               // ~tempHist->Add(histReader.read<TH1F>("baseline/mumu/"+sVar+"/"+sSample));
               TH1F* tempHist=histReader.read<TH1F>(sel+"/ee/"+sVar+"/"+sSample);
               tempHist->Add(histReader.read<TH1F>(sel+"/emu/"+sVar+"/"+sSample));
               tempHist->Add(histReader.read<TH1F>(sel+"/mumu/"+sVar+"/"+sSample));
               hs.addFilledHist(loc,sSample,*(tempHist));
            }
            else{
               TString loc;
               loc=sPresel+sVar;
               TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
               hs.addFilledHist(loc,sSample,*(tempHist));
            }
         }
      }
      if(!add_tt)hs.normHists();    //only normalize here if ttbar is not added, if ttbar added normalize later
   }
    
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
         else if (sPresel.Contains("combined/")) cat="all";
         TLatex label=gfx::cornerLabel(cat,1);
         
         auto hist_tW=hs.getHistogram(loc,"SingleTop");
         auto hist_tWDS=hs.getHistogram(loc,"SingleTop_TWDS");
         auto hist_tt=hs.getHistogram(loc,"TTbar_diLepton");
         
         if(add_tt){
            hist_tW->Add(hist_tt);
            hist_tWDS->Add(hist_tt);
            hist_tW->Scale(1./hist_tW->Integral());
            hist_tWDS->Scale(1./hist_tWDS->Integral());
         }
         
         hist_tW->SetStats(0);
         hist_tW->Draw("axis");
         hist_tW->SetLineColor(kGray);
         hist_tWDS->SetLineColor(kBlue);
         for (auto const &h: {hist_tW,hist_tWDS}) {
            if (hist_tW->GetNbinsX()>20) h->Rebin(2);
            h->SetStats(0);
            h->SetMarkerSize(0);
            h->Draw("same e");
         }
         gfx::LegendEntries le;
         if(add_tt){
            le.append(*hist_tW,"tt+SingleTop","l");
            le.append(*hist_tWDS,"tt+SingleTop DS","l");
         }
         else{
            le.append(*hist_tW,"SingleTop","l");
            le.append(*hist_tWDS,"SingleTop DS","l");
         }
         TLegend leg=le.buildLegend(.6,.8,1-gPad->GetRightMargin(),-1,1);
         leg.Draw();
         label.Draw();
         
         spcan.cdLow();
         TH1F hRatio=hist::getRatio(*hist_tW,*hist_tW,"Ratio.",hist::ONLY1);
         TH1F hRatio_amc=hist::getRatio(*hist_tWDS,*hist_tW,"Ratio.",hist::ONLY1);
         hRatio.SetStats(false);
         hRatio.SetMarkerSize(0);
         if(add_tt){
            hRatio.SetMaximum(1.15);
            hRatio.SetMinimum(0.85);
         }
         else{
            hRatio.SetMaximum(1.45);
            hRatio.SetMinimum(0.55);
         }
         hRatio_amc.SetMarkerSize(0);
         hRatio.Draw("axis e0");
         hRatio.Draw("same pe0");
         hRatio_amc.Draw("same pe0");
         if(add_tt){
            saver.save(spcan,"add_tt/"+loc,true,true);
         }
         else{
            saver.save(spcan,loc,true,true);
         }
      }
   }
}

extern "C"
void run()
{
   // ~compare_tt_powheg_madgraph();
   
   // ~compare_tW_ds_dr();
   compare_tW_ds_dr(true);
}
