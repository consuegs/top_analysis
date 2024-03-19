//Script to compare DNN inputs for POWHEG and amcatnlo 

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/tunfoldPlottingHelper.hpp"

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

void compare_tW_ds_dr(const bool add_tt=false, const bool compare_bb4l=false){
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_tW_ds_dr");
   
   std::vector<TString> samplesToImport={"SingleTop","SingleTop_DS"};
   samplesToImport.push_back("TTbar_diLepton");
   samplesToImport.push_back("TTbar_diLepton_tau");
   // ~samplesToImport.push_back("bb4l");
   samplesToImport.push_back("bb4l_new");
   
   std::map<TString,std::vector<TString>> msPresel_vVars={};
   for(TString channel:{"ee/","emu/","mumu/","combined/"}){
      // ~msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("baseline/"+channel,
      // ~{"PuppiMET","METunc_Puppi","MET","HT","nJets","n_Interactions","Lep1_flavor","Lep2_flavor","Lep1_pt","Lep1_phi","Lep1_eta","Lep1_E","Lep2_pt","Lep2_phi","Lep2_eta","Lep2_E","Jet1_pt","Jet1_phi","Jet1_eta","Jet1_E","Jet2_pt","Jet2_phi","Jet2_eta","Jet2_E","dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","dPhiJet1Jet2","METsig","MHT","MT","looseLeptonVeto","dPhiMETnearJet_Puppi","dPhiMETfarJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiMETbJet_Puppi","dPhiLep1bJet","dPhiLep1Jet1","mLL","PFMET_phi","PuppiMET_phi","CaloMET","CaloMET_phi","MT2","vecsum_pT_allJet","vecsum_pT_l1l2_allJet","mass_l1l2_allJet","ratio_vecsumpTlep_vecsumpTjet","mjj","DNN_MET_pT"}));
      // ~{"DNN_MET_pT","DNN_MET_dPhi_nextLep"}));
      // ~{"sum_mlb"}));
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("genParticles/"+channel,{"genMET"}));
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
      if(!add_tt && !compare_bb4l)hs.normHists();    //only normalize here if ttbar is not added, if ttbar added normalize later
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
         auto hist_tWDS=hs.getHistogram(loc,"SingleTop_DS");
         auto hist_tt=hs.getHistogram(loc,"TTbar_diLepton");
         auto hist_tt_tau=hs.getHistogram(loc,"TTbar_diLepton_tau");
         // ~auto hist_bb4l=hs.getHistogram(loc,"bb4l");
         auto hist_bb4l=hs.getHistogram(loc,"bb4l_new");
         
         if(add_tt){
            hist_tW->Add(hist_tt);
            hist_tWDS->Add(hist_tt);
            hist_tW->Scale(1./hist_tW->Integral());
            hist_tWDS->Scale(1./hist_tWDS->Integral());
         }
         else if(compare_bb4l){
            hist_tW->Add(hist_tt);
            hist_tWDS->Add(hist_tt);
            // ~hist_tW->Add(hist_tt_tau);
            // ~hist_tWDS->Add(hist_tt_tau);
            hist_tW->Scale(1./hist_tW->Integral());
            hist_tWDS->Scale(1./hist_tWDS->Integral());
            hist_bb4l->Scale(1./hist_bb4l->Integral());
         }
         
         if(sVar=="genMET")hist_tWDS->GetXaxis()->SetTitle("gen. %MET");
         hist_tWDS->GetYaxis()->SetTitle("normalized distribution");
         hist_tWDS->SetStats(0);
         hist_tWDS->Draw("axis");
         hist_tW->SetLineColor(kRed);
         hist_tWDS->SetLineColor(kGreen+2);
         hist_bb4l->SetLineColor(kBlue);
         if(compare_bb4l) {
            for (auto const &h: {hist_tW,hist_tWDS,hist_bb4l}) {
               // ~if (h->GetNbinsX()>20) h->Rebin(2);
               h->Rebin(4);
               h->SetStats(0);
               h->SetMarkerSize(0);
               h->Draw("same e");
            }
         }
         else{
            for (auto const &h: {hist_tW,hist_tWDS}) {
               // ~if (h->GetNbinsX()>20) h->Rebin(2);
               h->Rebin(4);
               h->SetStats(0);
               h->SetMarkerSize(0);
               h->Draw("same e");
            }
         }
         gfx::LegendEntries le;
         if(add_tt){
            le.append(*hist_tW,"tt+SingleTop","l");
            le.append(*hist_tWDS,"tt+SingleTop DS","l");
         }
         else if(compare_bb4l){
            le.append(*hist_tW,"tt+SingleTop","l");
            le.append(*hist_tWDS,"tt+SingleTop DS","l");
            le.append(*hist_bb4l,"bb4l","l");
         }
         else{
            le.append(*hist_tW,"Single top DR","le");
            le.append(*hist_tWDS,"Single top DS","le");
         }
         TLegend leg=le.buildLegend(.6,.8,1-gPad->GetRightMargin(),-1,1);
         leg.Draw();
         label.Draw();
         
         spcan.cdLow();
         TH1F hRatio=hist::getRatio(*hist_tW,*hist_tW,"Ratio",hist::ONLY1);
         TH1F hRatio_amc=hist::getRatio(*hist_tWDS,*hist_tW,"Ratio",hist::ONLY1);
         TH1F hRatio_bb4l=hist::getRatio(*hist_bb4l,*hist_tW,"Ratio",hist::ONLY1);
         if(compare_bb4l){
            hRatio=hist::getRatio(*hist_tW,*hist_tW,"Ratio",hist::ONLY1);
            hRatio_amc=hist::getRatio(*hist_tWDS,*hist_tW,"Ratio",hist::ONLY1);
            hRatio_bb4l=hist::getRatio(*hist_bb4l,*hist_tW,"Ratio",hist::ONLY1);
         }
         hRatio.SetStats(false);
         hRatio.SetMarkerSize(0);
         if(add_tt){
            hRatio.SetMaximum(1.15);
            hRatio.SetMinimum(0.85);
         }
         else{
            hRatio.SetMaximum(1.25);
            hRatio.SetMinimum(0.75);
         }
         hRatio_amc.SetMarkerSize(0);
         hRatio_amc.GetYaxis()->SetTitleOffset(0.45);
         hRatio_amc.Draw("axis e0");
         hRatio.Draw("same pe0");
         hRatio_amc.Draw("same pe0");
         if(!compare_bb4l && !add_tt) hRatio_amc.GetYaxis()->SetTitle("DS/DR");
         if(compare_bb4l) hRatio_bb4l.Draw("same pe0");
         if(add_tt){
            saver.save(spcan,"add_tt/"+loc,true,true);
         }
         else if(compare_bb4l){
            saver.save(spcan,"bb4l/"+loc,true,true);
         }
         else{
            saver.save(spcan,loc,true,true);
         }
         
         //Print chi2 comparison
         std::pair<float,int> chi_DR = tunfoldplotting::getChi2NDF(hist_tW,hist_bb4l,true);
         std::pair<float,int> chi_DS = tunfoldplotting::getChi2NDF(hist_tWDS,hist_bb4l,true);
         
         std::cout<<"chi2/ndf(DR vs. bb4l): "<<chi_DR.first<<"/"<<chi_DR.second<<std::endl;
         std::cout<<"chi2/ndf(DS vs. bb4l): "<<chi_DS.first<<"/"<<chi_DS.second<<std::endl;
      }
   }
}


// compare hDamp shifted samples between alternative MC and DCTR approach
void compare_hDamp(){
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_hdamp_DCTR");
   
   io::RootFileReader histReader_DCTR_UP(TString::Format("multiHists/MATCH_DCTR_UP/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader_DCTR_DOWN(TString::Format("multiHists/MATCH_DCTR_DOWN/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader_UP(TString::Format("multiHists/MATCH_UP/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader_DOWN(TString::Format("multiHists/MATCH_DOWN/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader_nominal(TString::Format("multiHists/Nominal/histograms_merged_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   std::map<TString,std::vector<TString>> msPresel_vVars={};
   // ~for(TString channel:{"ee/","emu/","mumu/","combined/"}){
   for(TString channel:{"combined/"}){
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("baseline/"+channel,
      // ~{"PuppiMET","METunc_Puppi","MET","HT","nJets","n_Interactions","Lep1_flavor","Lep2_flavor","Lep1_pt","Lep1_phi","Lep1_eta","Lep1_E","Lep2_pt","Lep2_phi","Lep2_eta","Lep2_E","Jet1_pt","Jet1_phi","Jet1_eta","Jet1_E","Jet2_pt","Jet2_phi","Jet2_eta","Jet2_E","dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","dPhiJet1Jet2","METsig","MHT","MT","looseLeptonVeto","dPhiMETnearJet_Puppi","dPhiMETfarJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiMETbJet_Puppi","dPhiLep1bJet","dPhiLep1Jet1","mLL","PFMET_phi","PuppiMET_phi","CaloMET","CaloMET_phi","MT2","vecsum_pT_allJet","vecsum_pT_l1l2_allJet","mass_l1l2_allJet","ratio_vecsumpTlep_vecsumpTjet","mjj","DNN_MET_pT"}));
      {"DNN_MET_pT","DNN_MET_dPhi_nextLep"}));
      // ~{"sum_mlb"}));
      // ~msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("genParticles/"+channel,{"genMET","pT_nunu"}));
   }
   
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc;
         loc=sPresel+sVar;
         
         //Define empty histograms
         TH1F* hist_DCTR_UP;
         TH1F* hist_DCTR_DOWN;
         TH1F* hist_UP;
         TH1F* hist_DOWN;
         TH1F* hist_nominal;
         
         if(sPresel.Contains("combined/")){
            TString sel=((TObjString*)sPresel.Tokenize("/")->First())->String();
            
            hist_DCTR_UP=histReader_DCTR_UP.read<TH1F>(sel+"/ee/"+sVar+"/TTbar_diLepton");
            hist_DCTR_UP->Add(histReader_DCTR_UP.read<TH1F>(sel+"/emu/"+sVar+"/TTbar_diLepton"));
            hist_DCTR_UP->Add(histReader_DCTR_UP.read<TH1F>(sel+"/mumu/"+sVar+"/TTbar_diLepton"));
            
            hist_DCTR_DOWN=histReader_DCTR_DOWN.read<TH1F>(sel+"/ee/"+sVar+"/TTbar_diLepton");
            hist_DCTR_DOWN->Add(histReader_DCTR_DOWN.read<TH1F>(sel+"/emu/"+sVar+"/TTbar_diLepton"));
            hist_DCTR_DOWN->Add(histReader_DCTR_DOWN.read<TH1F>(sel+"/mumu/"+sVar+"/TTbar_diLepton"));
            
            hist_UP=histReader_UP.read<TH1F>(sel+"/ee/"+sVar+"/TTbar_diLepton_MATCH_UP");
            hist_UP->Add(histReader_UP.read<TH1F>(sel+"/emu/"+sVar+"/TTbar_diLepton_MATCH_UP"));
            hist_UP->Add(histReader_UP.read<TH1F>(sel+"/mumu/"+sVar+"/TTbar_diLepton_MATCH_UP"));
            
            hist_DOWN=histReader_DOWN.read<TH1F>(sel+"/ee/"+sVar+"/TTbar_diLepton_MATCH_DOWN");
            hist_DOWN->Add(histReader_DOWN.read<TH1F>(sel+"/emu/"+sVar+"/TTbar_diLepton_MATCH_DOWN"));
            hist_DOWN->Add(histReader_DOWN.read<TH1F>(sel+"/mumu/"+sVar+"/TTbar_diLepton_MATCH_DOWN"));
            
            hist_nominal=histReader_nominal.read<TH1F>(sel+"/ee/"+sVar+"/TTbar_diLepton");
            hist_nominal->Add(histReader_nominal.read<TH1F>(sel+"/emu/"+sVar+"/TTbar_diLepton"));
            hist_nominal->Add(histReader_nominal.read<TH1F>(sel+"/mumu/"+sVar+"/TTbar_diLepton"));
         }
         else{
            hist_DCTR_UP=histReader_DCTR_UP.read<TH1F>(loc+"/TTbar_diLepton");
            hist_DCTR_DOWN=histReader_DCTR_DOWN.read<TH1F>(loc+"/TTbar_diLepton");
            hist_UP=histReader_UP.read<TH1F>(loc+"/TTbar_diLepton_MATCH_UP");
            hist_DOWN=histReader_DOWN.read<TH1F>(loc+"/TTbar_diLepton_MATCH_DOWN");
            hist_nominal=histReader_nominal.read<TH1F>(loc+"/TTbar_diLepton");
         }
         
         int rebinFactor = 1;
         if (sVar == "genMET" || sVar == "pT_nunu") rebinFactor = 4;
         else if (sVar == "DNN_MET_pT") rebinFactor = 25;
         else if (sVar == "DNN_MET_dPhi_nextLep") rebinFactor = 16;
         
         hist_DCTR_UP->Rebin(rebinFactor);
         hist_DCTR_DOWN->Rebin(rebinFactor);
         hist_UP->Rebin(rebinFactor);
         hist_DOWN->Rebin(rebinFactor);
         hist_nominal->Rebin(rebinFactor);
   
         gfx::SplitCan spcan;
         spcan.cdUp();
         spcan.pU_.SetLogy();
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         else if (sPresel.Contains("combined/")) cat="all";
         TLatex label=gfx::cornerLabel(cat,1);
         
         hist_DCTR_UP->SetStats(0);
         hist_nominal->Draw("axis");
         
         hist_DCTR_UP->SetLineColor(kRed);
         hist_DCTR_DOWN->SetLineColor(kGreen+2);
         hist_UP->SetLineColor(kRed);
         hist_DOWN->SetLineColor(kGreen+2);
         
         hist_nominal->SetMarkerSize(0);
         hist_DCTR_UP->SetMarkerSize(0);
         hist_DCTR_DOWN->SetMarkerSize(0);
         hist_UP->SetMarkerSize(0);
         hist_DOWN->SetMarkerSize(0);
         
         hist_DCTR_UP->SetLineStyle(7);
         hist_DCTR_DOWN->SetLineStyle(7);
         
         hist_nominal->Draw("same");
         hist_DCTR_UP->Draw("same");
         hist_DCTR_DOWN->Draw("same");
         hist_UP->Draw("same");
         hist_DOWN->Draw("same");

         gfx::LegendEntries le;
         le.append(*hist_nominal,"Nominal","l");
         le.append(*hist_UP,"hDamp up","l");
         le.append(*hist_DOWN,"hDamp down","l");
         le.append(*hist_DCTR_UP,"hDamp up (DCTR)","l");
         le.append(*hist_DCTR_DOWN,"hDamp down (DCTR)","l");
         
         TLegend leg=le.buildLegend(.7,.8,1-gPad->GetRightMargin(),-1,1);
         leg.Draw();
         label.Draw();
         
         //Get chi2 and KS comparison
         std::pair<float,int> chi_UP = tunfoldplotting::getChi2NDF(hist_DCTR_UP,hist_UP,true);
         std::pair<float,int> chi_DOWN = tunfoldplotting::getChi2NDF(hist_DCTR_DOWN,hist_DOWN,true);
         double KS_UP = hist_DCTR_UP->KolmogorovTest(hist_UP);
         double KS_DOWN = hist_DCTR_DOWN->KolmogorovTest(hist_DOWN);
         
         std::cout<<"chi2/ndf(UP): "<<chi_UP.first<<"/"<<chi_UP.second<<std::endl;
         std::cout<<"chi2/ndf(DOWN): "<<chi_DOWN.first<<"/"<<chi_DOWN.second<<std::endl;
         std::cout<<"KS-Test(UP): "<<hist_DCTR_UP->KolmogorovTest(hist_UP)<<std::endl;
         std::cout<<"KS-Test(DOWN): "<<hist_DCTR_DOWN->KolmogorovTest(hist_DOWN)<<std::endl;
         
         TLatex testResults;
         testResults.SetTextSize(0.03);
         testResults.SetNDC();
         testResults.DrawLatex(0.2,0.6,TString::Format("#splitline{#chi^{2}/ndf (UP) = %.1f/%i}{#chi^{2}/ndf (DOWN) = %.1f/%i}",chi_UP.first,chi_UP.second,chi_DOWN.first,chi_DOWN.second));
         testResults.DrawLatex(0.2,0.5,TString::Format("#splitline{KS-Test (UP) = %.2f}{KS-Test (DOWN) = %.2f}",KS_UP,KS_DOWN));
         
         spcan.cdLow();
         TH1F hRatio_nominal=hist::getRatio(*hist_nominal,*hist_nominal,"Ratio",hist::ONLY1);
         TH1F hRatio_DCTR_UP=hist::getRatio(*hist_DCTR_UP,*hist_nominal,"Ratio",hist::ONLY1);
         TH1F hRatio_DCTR_DOWN=hist::getRatio(*hist_DCTR_DOWN,*hist_nominal,"Ratio",hist::ONLY1);
         TH1F hRatio_UP=hist::getRatio(*hist_UP,*hist_nominal,"Ratio",hist::ONLY1);
         TH1F hRatio_DOWN=hist::getRatio(*hist_DOWN,*hist_nominal,"Ratio",hist::ONLY1);
         hRatio_nominal.SetStats(false);
         hRatio_nominal.SetMarkerSize(0);
         
         hRatio_nominal.SetMaximum(1.15);
         hRatio_nominal.SetMinimum(0.85);
         
         hRatio_nominal.Draw();
         hRatio_DCTR_UP.Draw("same");
         hRatio_DCTR_DOWN.Draw("same");
         hRatio_UP.Draw("same");
         hRatio_DOWN.Draw("same");
         
         saver.save(spcan,loc,true,true);
      }
   }
}

extern "C"
void run()
{
   // ~compare_tt_powheg_madgraph();
   
   // ~compare_tW_ds_dr();
   // ~compare_tW_ds_dr(true);
   // ~compare_tW_ds_dr(false,false);
   // ~compare_tW_ds_dr(false,true);
   
   compare_hDamp();
}
