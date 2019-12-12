//Script to plot distributions for events with good and bad dPhi(MET,next l) resolution

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
   // ~io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_dPhiResolution");
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_dPhiResolution_dilepton");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("diff_dPhiResolution%.1f",cfg.processFraction*100));
   TCanvas can;
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/",{"met","met1000","metSig","mll","pTlep1","pTlep2","etalep1","etalep2","pTbJet","etabJet","dphi_metJet","dphi_metLeadJet","dphi_metLead2Jet",
         "dphi_metFarJet","dphi_metBJet","dphi_bJetLep1","dphi_bJetLep2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","njets","mt2","dR_Lep1Lep2","ST","HT",
         "METoverHT","METoverSUMpt","sum_STHT","mt_MetLep1","mt_MetLep2","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep",
         "dphi_metLepsum","diff_met_Puppi","diff_met_NoHF","diff_met_Calo","diff_met_Raw","HTnormed_diff_met_Puppi","HTnormed_diff_met_NoHF",
         "HTnormed_diff_met_Calo","HTnormed_diff_met_Raw","n_Interactions"}},
      {"genParticles/",{"pT_nunu","genMet","diff_ptNuNu_genMET","diff_Met_genMET","diff_Met_genMET_norm","diff_Met_genMET_normSUM","diff_ptNuNu_Met","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep","dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      };
   
   for (auto const &sPresel_vVars:msPresel_vVars){
      TString const &sPresel=sPresel_vVars.first;
      
      for (TString sVar:sPresel_vVars.second){
            
         TH1F* hist_good;
         TH1F* hist_bad;
         
         // ~hist_good=(TH1F*) histReader.read<TH1F>("goodPhiRes/"+sPresel+sVar+"/TTbar");
         // ~hist_bad=(TH1F*) histReader.read<TH1F>("badPhiRes/"+sPresel+sVar+"/TTbar");
         // ~hist_good=(TH1F*) histReader.read<TH1F>("goodMETRes/"+sPresel+sVar+"/TTbar");
         // ~hist_bad=(TH1F*) histReader.read<TH1F>("badMETRes/"+sPresel+sVar+"/TTbar");
         hist_good=(TH1F*) histReader.read<TH1F>("goodMETRes/"+sPresel+sVar+"/TTbar_diLepton");
         hist_bad=(TH1F*) histReader.read<TH1F>("badMETRes/"+sPresel+sVar+"/TTbar_diLepton");
         
         if (sVar!="njets" and sVar!="nBjets"){
            // ~hist_bad->Rebin(4);
            // ~hist_good->Rebin(4);
            hist_bad->Rebin(10);
            hist_good->Rebin(10);
         }
            
         hist_good->Scale(1.0/(hist_good->Integral()));   //Normalize the hist to the integral
         hist_bad->Scale(1.0/(hist_bad->Integral()));
         
         can.cd();
         can.SetLogz();
         
         hist_bad->GetYaxis()->SetTitle("normalized distribution");
         hist_good->SetLineColor(kBlue);
         hist_bad->SetLineColor(kRed);
         
         hist_bad->SetStats(0);
         hist_good->SetStats(0);
         
         hist_bad->Draw("hist");
         hist_good->Draw("same hist");
         
         hist_bad->SetMaximum(1.3*std::max(hist_bad->GetMaximum(),hist_good->GetMaximum()));
         
         gfx::LegendEntries legE;
         // ~legE.append(*hist_good,"good resolution (<0.3)","l");
         // ~legE.append(*hist_bad,"bad resolution (>0.3)","l");
         // ~legE.append(*hist_good,"good resolution (<50)","l");
         // ~legE.append(*hist_bad,"bad resolution (>50)","l");
         legE.append(*hist_good,"true","l");
         legE.append(*hist_bad,"migrated","l");
         TLegend leg=legE.buildLegend(.5,.65,0.75,.9,1);
         leg.SetTextSize(0.035);
         leg.Draw();
         
         // ~TLatex label=gfx::cornerLabel("MET>120 GeV, dPhi>1.4",1);
         TLatex label=gfx::cornerLabel("MET>230 GeV, dPhi>1.4",1);
         // ~TLatex label=gfx::cornerLabel("All bins",1);
         // ~TLatex label2=gfx::cornerLabel("|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T}<0.05",2);
         label.Draw();
         // ~label2.Draw();
         
         can.RedrawAxis();
         // ~TString plotLoc=sPresel+"/"+sVar;
         TString plotLoc="METresolution/"+sPresel+"/"+sVar;
         saver.save(can,plotLoc,true,true);
         // ~if(sPresel=="baseline/")can.SaveAs(sVar+".pdf");
         can.Clear();
      }
   }
}
