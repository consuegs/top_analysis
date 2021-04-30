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

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"compare_ttMC");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   std::vector<TString> samplesToImport={"TTbar_diLepton","TTbar_amcatnlo"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={};
   for(TString channel:{"ee/","emu/","mumu/","combined/"}){
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("baseline/"+channel,
      {"PuppiMET","METunc_Puppi","MET","HT","nJets","n_Interactions","Lep1_flavor","Lep2_flavor","Lep1_pt","Lep1_phi","Lep1_eta","Lep1_E","Lep2_pt","Lep2_phi","Lep2_eta","Lep2_E","Jet1_pt","Jet1_phi","Jet1_eta","Jet1_E","Jet2_pt","Jet2_phi","Jet2_eta","Jet2_E","dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","dPhiJet1Jet2","METsig","MHT","MT","looseLeptonVeto","dPhiMETnearJet_Puppi","dPhiMETfarJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiMETbJet_Puppi","dPhiLep1bJet","dPhiLep1Jet1","mLL","PFMET_phi","PuppiMET_phi","CaloMET","CaloMET_phi","MT2","vecsum_pT_allJet","vecsum_pT_l1l2_allJet","mass_l1l2_allJet","ratio_vecsumpTlep_vecsumpTjet","mjj"}));
      msPresel_vVars.insert(std::pair<TString,std::vector<TString>>("genParticles/"+channel,{"genMet"}));
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
         saver.save(spcan,loc);
      }
   }
}
