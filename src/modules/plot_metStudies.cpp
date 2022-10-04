//Script to plot studies from diff_MET script

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
   // ~TString sampleName = "TTbar_diLepton";
   TString sampleName = "T2tt_525_350";
   
   io::RootFileSaver saver(TString::Format("minTree_analyze/plots%.1f.root",cfg.processFraction*100),"plot_metStudies/"+sampleName);
   
   TCanvas can;
   can.cd();
   io::RootFileReader histReader("minTree_analyze/histograms_Nominal.root");
   
   //Plotting comparison between PF and Puppi for different selections
   // ~for (TString sel : {"baseline_met120","baseline_met120_ee","baseline_met120_emu","baseline_met120_mumu","baseline_genmet120","baseline_met200","baseline_met120_230","baseline_met230","baseline_matchedLep_met120","baseline"}){
   for (TString sel : {"baseline_genmet120/all"}){
         // ~std::vector<TString> var_vec = {"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep"};
         std::vector<TString> var_vec = {"GenMetDiffMETRel_dPhiMETLep"};
         if (sel=="baseline_met120") {
            var_vec.push_back("GenMetDiffMETRelReco_dPhiMETLep");
            var_vec.push_back("GenMetDiffMETRel_dPhigenMETLep");
            var_vec.push_back("Met_dPhiMETLep");
         }
         if(sel=="baseline") var_vec = {"GenMetDiffMETRel_dPhiMETLep"};
         
         for (TString var :var_vec){
            
            TH2F* TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+sel+"/"+var+"/"+sampleName);
            TH2F* TTbar_2D_Puppi = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+sel+"/"+var+"_Puppi/"+sampleName);
            
            TProfile TTbar_profile;
            TProfile TTbar_profile_Puppi;
            TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
            TTbar_profile_Puppi=*(TTbar_2D_Puppi->ProfileX("Profile"));
            
            TTbar_profile.SetLineColor(kBlue);
            TTbar_profile_Puppi.SetLineColor(kRed);
            TTbar_profile.SetMarkerColor(kBlue);
            TTbar_profile_Puppi.SetMarkerColor(kRed);
            TTbar_profile.SetMarkerSize(1);
            TTbar_profile_Puppi.SetMarkerSize(1);
            
            TTbar_profile.SetStats(0);
            TTbar_profile.SetMaximum(0.1);
            TTbar_profile.SetMinimum(-0.7);
            TTbar_profile.GetYaxis()->SetTitle("mean[(genMET-p_{T}^{miss})/genMET]");
            if (var=="GenMetDiffMETRelReco_dPhiMETLep") {
               TTbar_profile.SetMaximum(0.1);
               TTbar_profile.SetMinimum(-0.7);
               TTbar_profile.GetYaxis()->SetTitle("mean[(genMET-p_{T}^{miss})/p_{T}^{miss}]");
            }
            else if (var=="Met_dPhiMETLep") {
               TTbar_profile.SetMaximum(190);
               TTbar_profile.SetMinimum(120);
               TTbar_profile.GetYaxis()->SetTitle("mean(p_{T}^{miss}) (GeV)");
            }
            else if (var=="MetSig_dPhiMETLep") {
               TTbar_profile.SetMaximum(400);
               TTbar_profile.SetMinimum(20);
               TTbar_profile.GetYaxis()->SetTitle("mean(metSig)");
            }
            
            gPad->SetLeftMargin(0.15);
            TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
            TTbar_profile.Draw("e1");
            TTbar_profile_Puppi.Draw("same e1");
            
            gfx::LegendEntries le;
            le.append(TTbar_profile,"PF","pl");
            le.append(TTbar_profile_Puppi,"Puppi","pl");
            TLegend leg=le.buildLegend(.7,.2,1-1.5*gPad->GetRightMargin(),0.3,1);
            leg.Draw();
            
            TString stringLabel ="p_{T}^{miss}>120 GeV";
            if (sel=="baseline_genmet120/all") stringLabel="genMET>120 GeV";
            else if (sel=="baseline_met120_ee") stringLabel="ee, p_{T}^{miss}>120 GeV";
            else if (sel=="baseline_met120_emu") stringLabel="e#mu, p_{T}^{miss}>120 GeV";
            else if (sel=="baseline_met120_mumu") stringLabel="#mu#mu, p_{T}^{miss}>120 GeV";
            else if (sel=="baseline_met200") stringLabel="p_{T}^{miss}>200 GeV";
            else if (sel=="baseline_met230") stringLabel="p_{T}^{miss}>230 GeV";
            else if (sel=="baseline_met120_230") stringLabel="120 GeV<p_{T}^{miss}<230 GeV";
            else if (sel=="baseline_matchedLep_met120") stringLabel="p_{T}^{miss}>120 GeV, matched genLeptons";
            TLatex label=gfx::cornerLabel(stringLabel,1);
            label.Draw();
            saver.save(can,"pfVSpuppi/"+sel+"/"+var,true,true,true);
         }
   }
   
   //Plotting comparison between PF and Puppi comparing full and met>120
   for (auto dummy: {1}){
      TH2F* TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"baseline_met120/all/GenMetDiffMETRel_dPhiMETLep/"+sampleName);
      TH2F* TTbar_2D_Puppi = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"baseline_met120/all/GenMetDiffMETRel_dPhiMETLep_Puppi/"+sampleName);
      TH2F* TTbar_2D_full = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"all/GenMetDiffMETRel_dPhiMETLep/"+sampleName);
      TH2F* TTbar_2D_full_Puppi = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"all/GenMetDiffMETRel_dPhiMETLep_Puppi/"+sampleName);
      
      TProfile TTbar_profile;
      TProfile TTbar_profile_Puppi;
      TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
      TTbar_profile_Puppi=*(TTbar_2D_Puppi->ProfileX("Profile"));
      
      TProfile TTbar_profile_full;
      TProfile TTbar_profile_full_Puppi;
      TTbar_profile_full=*(TTbar_2D_full->ProfileX("Profile"));
      TTbar_profile_full_Puppi=*(TTbar_2D_full_Puppi->ProfileX("Profile"));
      
      TTbar_profile.SetLineColor(kBlue);
      TTbar_profile_Puppi.SetLineColor(kRed);
      TTbar_profile.SetMarkerColor(kBlue);
      TTbar_profile_Puppi.SetMarkerColor(kRed);
      TTbar_profile.SetMarkerSize(1);
      TTbar_profile_Puppi.SetMarkerSize(1);
      
      TTbar_profile_full.SetLineColor(kBlue);
      TTbar_profile_full_Puppi.SetLineColor(kRed);
      TTbar_profile_full.SetMarkerColor(kBlue);
      TTbar_profile_full_Puppi.SetMarkerColor(kRed);
      TTbar_profile_full.SetMarkerSize(1);
      TTbar_profile_full_Puppi.SetMarkerSize(1);
      TTbar_profile_full.SetMarkerStyle(4);
      TTbar_profile_full_Puppi.SetMarkerStyle(4);
      
      TTbar_profile.SetStats(0);
      TTbar_profile.SetMaximum(0.1);
      TTbar_profile.SetMinimum(-0.7);
      TTbar_profile.GetYaxis()->SetTitle("mean[(genMET-p_{T}^{miss})/genMET]");
      
      gPad->SetLeftMargin(0.15);
      TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
      TTbar_profile.Draw("e1");
      TTbar_profile_Puppi.Draw("same e1");
      TTbar_profile_full.Draw("same e1");
      TTbar_profile_full_Puppi.Draw("same e1");
      
      gfx::LegendEntries le;
      le.append(TTbar_profile_full,"PF","pl");
      le.append(TTbar_profile_full_Puppi,"Puppi","pl");
      le.append(TTbar_profile,"PF>120 GeV","pl");
      le.append(TTbar_profile_Puppi,"Puppi> 120 GeV","pl");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      saver.save(can,"pfVSpuppi/baseline_vs_MET120/GenMetDiffMETRel_dPhiMETLep",true,true,true);
   }
   
   //Plotting genMet as a function of dPhi (met120)
   TH2F* TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"baseline_met120/all/genMet_dPhiMETLep_Puppi/"+sampleName);
   TProfile TTbar_profile;
   TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   gPad->SetLeftMargin(0.15);
   TTbar_profile.SetStats(0);
   TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
   TTbar_profile.GetYaxis()->SetTitle("mean(genMET) (GeV)");
   TTbar_profile.Draw("e1");
   TString stringLabel ="p_{T}^{miss}>120 GeV";
   TLatex label=gfx::cornerLabel(stringLabel,2);
   label.Draw();
   saver.save(can,"genMetVSdPhi/met120",true,true,true); 
   
   //Plotting genMet as a function of dPhi
   TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"all/genMet_dPhiMETLep_Puppi/"+sampleName);
   TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   gPad->SetLeftMargin(0.15);
   TTbar_profile.SetStats(0);
   TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
   TTbar_profile.GetYaxis()->SetTitle("mean(genMET) (GeV)");
   TTbar_profile.Draw("e1");
   saver.save(can,"genMetVSdPhi/baseline",true,true,true); 
   
   //Plotting lead top pT as a function of dPhi (met120)
   TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"baseline_met120/all/leadTopPT_dPhiMETLep_Puppi/"+sampleName);
   TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   gPad->SetLeftMargin(0.15);
   TTbar_profile.SetStats(0);
   TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
   TTbar_profile.GetYaxis()->SetTitle("mean(p_{T}(lead. top)) (GeV)");
   TTbar_profile.Draw("e1");
   label.Draw();
   saver.save(can,"leadTopPtVSdPhi/met120",true,true,true); 
   
   //Plotting lead top pT as a function of dPhi 
   TTbar_2D = histReader.read<TH2F>(TString::Format("analyze_minTrees%.1f/",cfg.processFraction*100)+"all/leadTopPT_dPhiMETLep_Puppi/"+sampleName);
   TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   gPad->SetLeftMargin(0.15);
   TTbar_profile.SetStats(0);
   TTbar_profile.GetYaxis()->SetTitleOffset(1.2);
   TTbar_profile.GetYaxis()->SetTitle("mean(p_{T}(lead. top)) (GeV)");
   TTbar_profile.Draw("e1");
   saver.save(can,"leadTopPtVSdPhi/baseline",true,true,true); 
   
   // ~//Plotting nVertices as a function of dPhi
   // ~TTbar_2D = histReader.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline_met120/nVertex/"+sampleName);
   // ~TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   // ~gPad->SetLeftMargin(0.13);
   // ~TTbar_profile.SetStats(0);
   // ~TTbar_profile.GetYaxis()->SetTitleOffset(1.1);
   // ~TTbar_profile.GetYaxis()->SetTitle("mean(nVertices)");
   // ~TTbar_profile.Draw("e1");  
   // ~label.Draw();
   // ~saver.save(can,"nVerticesVSdPhi/met120",true,true);  
   
   /*
   //Plotting dPhi gen and reco distr in same plot
   gfx::SplitCan spcan;
   spcan.pU_.cd();
   gPad->SetLeftMargin(0.13);
   TH1F* TTbar_reco = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/dPhiMETLep/"+sampleName); 
   TH1F* TTbar_gen = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/dPhiMETLep_gen/"+sampleName);
   
   TTbar_reco->Rebin(5);
   TTbar_gen->Rebin(5);
   
   TTbar_reco->GetYaxis()->SetTitleOffset(1.1);
   TTbar_reco->SetLineColor(kBlue);
   TTbar_gen->SetLineColor(kRed);
   TTbar_reco->SetMarkerSize(0);
   TTbar_gen->SetMarkerSize(0);
   TTbar_reco->SetStats(0);
   TTbar_reco->Draw("e1");
   TTbar_gen->Draw("same e1");
   
   gfx::LegendEntries le;
   le.append(*TTbar_reco,"reco","pl");
   le.append(*TTbar_gen,"gen","pl");;
   TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
   leg.Draw();
   
   spcan.pL_.cd();
   gPad->SetLeftMargin(0.13);
   TH1F ratio=hist::getRatio(*TTbar_reco,*TTbar_gen,hist::ONLY1);
   ratio.GetYaxis()->SetTitleOffset(0.2);
   ratio.SetLineColor(kRed);
   ratio.GetYaxis()->SetTitle("reco/gen");
   ratio.SetMaximum(1.2);
   ratio.SetMinimum(0.8);
   ratio.Draw("e1");
   saver.save(spcan,"dPhiComparison/fullSelection");
   
   //Plotting MET resolution as a function of nInteractions
   can.Clear();
   TTbar_2D = histReader.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/nVertex_vs_MetRes/"+sampleName);
   TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
   TProfile TTbar_profile_RMS=*(TTbar_2D->ProfileX("ProfileRMS",1,-1,"s"));
   TH1D RMS("","",100,0,100);
   for (int i=1; i<=TTbar_profile_RMS.GetNbinsX(); i++){
      RMS.SetBinContent(i,TTbar_profile_RMS.GetBinError(i));
      RMS.SetBinError(i,0);
   }
   gPad->SetLeftMargin(0.13);
   TTbar_profile.SetStats(0);
   TTbar_profile.GetYaxis()->SetTitleOffset(1.1);
   TTbar_profile.GetYaxis()->SetTitle("mean(genMET-p_{T}^{miss}) (GeV)");
   TTbar_profile.Draw("e1");
   RMS.Draw("same");
   saver.save(can,"nVertex_vs_MetRes",true,true);
   */
   
   
   //Plotting different METS and also comparision to SUSY
   /* 
   for (TString var :{"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep"}){
      
      TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"/"+sampleName);
      TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"/T2tt_650_350");
      TH2F* TTbar_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"_Puppi/"+sampleName);
      TH2F* OtherSample_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"_Puppi/T2tt_650_350");
      TH2F* TTbar_2D_Deep = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"_Deep/"+sampleName);
      TH2F* OtherSample_2D_Deep = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"_Deep/T2tt_650_350");
      // ~TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"/"+sampleName);
      // ~TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"/T2tt_650_350");
      // ~TH2F* TTbar_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"_Puppi/"+sampleName);
      // ~TH2F* OtherSample_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"_Puppi/T2tt_650_350");
      // ~TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"/"+sampleName);
      // ~TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"/T2tt_650_350");
      // ~TH2F* TTbar_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Puppi/"+sampleName);
      // ~TH2F* OtherSample_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Puppi/T2tt_650_350");
      // ~TH2F* TTbar_2D_Deep = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Deep/"+sampleName);
      // ~TH2F* OtherSample_2D_Deep = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Deep/T2tt_650_350");
      
      TProfile TTbar_profile;
      TProfile OtherSample_profile;
      TProfile TTbar_profile_Puppi;
      TProfile OtherSample_profile_Puppi;
      TProfile TTbar_profile_Deep;
      TProfile OtherSample_profile_Deep;
      TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
      OtherSample_profile=*(OtherSample_2D->ProfileX("Profile"));
      TTbar_profile_Puppi=*(TTbar_2D_Puppi->ProfileX("Profile"));
      OtherSample_profile_Puppi=*(OtherSample_2D_Puppi->ProfileX("Profile"));
      TTbar_profile_Deep=*(TTbar_2D_Deep->ProfileX("Profile"));
      OtherSample_profile_Deep=*(OtherSample_2D_Deep->ProfileX("Profile"));
      
      TTbar_profile.SetLineColor(kRed-6);
      OtherSample_profile.SetLineColor(kBlue);
      // ~TTbar_profile.SetMarkerSize(0);
      // ~OtherSample_profile.SetMarkerSize(0);
      TTbar_profile_Puppi.SetLineColor(kRed-6);
      OtherSample_profile_Puppi.SetLineColor(kBlue);
      TTbar_profile_Puppi.SetMarkerStyle(24);
      OtherSample_profile_Puppi.SetMarkerStyle(24);
      TTbar_profile_Deep.SetLineColor(kRed-6);
      OtherSample_profile_Deep.SetLineColor(kBlue);
      TTbar_profile_Deep.SetMarkerStyle(22);
      OtherSample_profile_Deep.SetMarkerStyle(22);
      
      TTbar_profile.SetStats(0);
      if (var=="GenMetDiffMETRel_dPhiMETLep") {
            TTbar_profile.SetMaximum(0.1);
            TTbar_profile.SetMinimum(-0.5);
            TTbar_profile.GetYaxis()->SetTitle("mean[genMET-p_{T}^{miss}]/genMET]");
      }
      else if (var=="MetSig_dPhiMETLep") {
            // ~TTbar_profile.SetMaximum(200);
            TTbar_profile.SetMaximum(400);
            TTbar_profile.SetMinimum(20);
            TTbar_profile.GetYaxis()->SetTitle("mean[metSig]");
      }
      
      gPad->SetLeftMargin(0.13);
      TTbar_profile.GetYaxis()->SetTitleOffset(1.0);
      TTbar_profile.Draw("e1");
      OtherSample_profile.Draw("same e1");
      TTbar_profile_Puppi.Draw("same e1");
      OtherSample_profile_Puppi.Draw("same e1");
      TTbar_profile_Deep.Draw("same e1");
      OtherSample_profile_Deep.Draw("same e1");
      
      gfx::LegendEntries le;
      le.append(TTbar_profile,"TTbar","pl");
      le.append(OtherSample_profile,"T2tt_650_350","pl");
      le.append(TTbar_profile_Puppi,"TTbar(Puppi)","pl");
      le.append(OtherSample_profile_Puppi,"T2tt_650_350(Puppi)","pl");
      le.append(TTbar_profile_Deep,"TTbar(Deep)","pl");
      le.append(OtherSample_profile_Deep,"T2tt_650_350(Deep)","pl");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>120GeV",1);
      // ~TLatex label=gfx::cornerLabel("p_{T}^{miss}>230GeV",1);
      // ~TLatex label=gfx::cornerLabel("120GeV<p_{T}^{miss}<230GeV",1);
      label.Draw();
      saver.save(can,var,true,true);
   }
   */
   
   /*
   //Plotting MET Res Comparison including BJet Regression
   for (TString selection:{"baseline","baseline_genmet120"}){
      can.Clear();
      can.cd();
      // ~gPad->SetLeftMargin(0.13);
      // ~TH1F* PFMET_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METres/"+sampleName); 
      // ~TH1F* Puppi_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresPuppi/"+sampleName);
      // ~TH1F* BReg_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresBJetRegr/"+sampleName);
      // ~TH1F* BRegLB_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresBJetLBRegr/"+sampleName);
      TH1F* PFMET_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METres/TTbar_diLepton_CUETP8M2"); 
      TH1F* Puppi_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresPuppi/TTbar_diLepton_CUETP8M2");
      TH1F* BReg_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresBJetRegr/TTbar_diLepton_CUETP8M2");
      TH1F* BRegLB_res = histReader.read<TH1F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+selection+"/METresBJetLBRegr/TTbar_diLepton_CUETP8M2");
      
      Puppi_res->SetStats(0);
      Puppi_res->Draw("axis");
      
      int color=1;
      gfx::LegendEntries le_res;
      std::vector<TString> names={"PFMET","PUPPI","BReg","BRegLB"};
      for (auto hist:{PFMET_res,Puppi_res,BReg_res,BRegLB_res}){
         hist->Rebin(2);
         hist->SetMarkerColor(color);
         hist->SetMarkerSize(0.8);
         hist->Draw("same");
         le_res.append(*hist,TString::Format(names[color-1]+", #mu=%.1f #sigma=%.1f",hist->GetMean(),hist->GetRMS()),"p");
         color++;
      }
      
      TLegend leg_res=le_res.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg_res.Draw();
      if (selection=="baseline") saver.save(can,"BRegrComparisons/CompareMetRes",true,true);
      else saver.save(can,"BRegrComparisons/CompareMetRes_genMet120",true,true);
   }
   */
}
