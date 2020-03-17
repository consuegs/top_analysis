//Script to plot different MET scale factors

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_metSF");
   
   //Plot ptnunu/met and ptnunu/genmet for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      for (TString phiRegion : {"","_phi1","_phi2","_phi3"}){
         TCanvas can;
         can.cd();
         
         // ~TGraph* SF_profile = histReader.read<TGraph>("binningUnfolding/METdiff_vs_MET_profile");
         // ~TGraph* SFgen_profile = histReader.read<TGraph>("binningUnfolding/METdiffgen_vs_MET_profile");
         TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/METdiff"+phiRegion+"_vs_MET_profile");
         TH1D* SFgen_profile = histReader.read<TH1D>("binningUnfolding/METdiffgen"+phiRegion+"_vs_MET_profile");
         
         SF_profile->SetLineColor(kBlue);
         SFgen_profile->SetLineColor(kRed);
         SF_profile->SetMarkerColor(kBlue);
         SFgen_profile->SetMarkerColor(kRed);
         
         SF_profile->SetStats(0);
         SF_profile->GetYaxis()->SetRangeUser(0.5,SF_profile->GetMaximum());
         SF_profile->Draw("hist");
         SFgen_profile->Draw("same hist");
         
         gfx::LegendEntries le;
         le.append(*SF_profile,"p_{T}^{#nu#nu}/p_{T}^{miss}","l");
         le.append(*SFgen_profile,"p_{T}^{#nu#nu}/genMET","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(phiRegion,1);
         label.Draw();
         saver.save(can,"ScaleFactor_"+sSample+phiRegion,true,true);
      }
   }
   
   //Plot ptnunu/met and ptnunu/genmet for each sample in one canvas
   
	for (TString phiRegion : {"","_phi1","_phi2","_phi3"}){
      TCanvas can;
      can.cd();
      io::RootFileReader histReader_incl(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100));
      io::RootFileReader histReader_dilepton(TString::Format("binningUnfolding_dilepton%.1f.root",cfg.processFraction*100));
      io::RootFileReader histReader_MadGraph(TString::Format("binningUnfolding_MadGraph%.1f.root",cfg.processFraction*100));
         
      TH1D* SF_profile = histReader_incl.read<TH1D>("binningUnfolding/METdiff"+phiRegion+"_vs_MET_profile");
      TH1D* SFgen_profile = histReader_incl.read<TH1D>("binningUnfolding/METdiffgen"+phiRegion+"_vs_MET_profile");
      TH1D* SF_profile_dilepton = histReader_dilepton.read<TH1D>("binningUnfolding/METdiff"+phiRegion+"_vs_MET_profile");
      TH1D* SFgen_profile_dilepton = histReader_dilepton.read<TH1D>("binningUnfolding/METdiffgen"+phiRegion+"_vs_MET_profile");
      TH1D* SF_profile_MadGraph = histReader_MadGraph.read<TH1D>("binningUnfolding/METdiff"+phiRegion+"_vs_MET_profile");
      TH1D* SFgen_profile_MadGraph = histReader_MadGraph.read<TH1D>("binningUnfolding/METdiffgen"+phiRegion+"_vs_MET_profile");
      
      SF_profile->SetLineColor(kBlue);
      SFgen_profile->SetLineColor(kBlue);
      SF_profile_dilepton->SetLineColor(kRed);
      SFgen_profile_dilepton->SetLineColor(kRed);
      SF_profile_MadGraph->SetLineColor(kGreen);
      SFgen_profile_MadGraph->SetLineColor(kGreen);
      
      SF_profile->SetStats(0);
      SF_profile->GetYaxis()->SetRangeUser(0.5,SF_profile->GetMaximum());
      SF_profile->Draw("hist");
      SFgen_profile->Draw("same hist");
      SF_profile_dilepton->Draw("same hist");
      SFgen_profile_dilepton->Draw("same hist");
      SF_profile_MadGraph->Draw("same hist");
      SFgen_profile_MadGraph->Draw("same hist");
      
      // ~gfx::LegendEntries le;
      // ~le.append(*SF_profile,"p_{T}^{#nu#nu}/p_{T}^{miss}","l");
      // ~le.append(*SFgen_profile,"p_{T}^{#nu#nu}/genMET","l");
      // ~TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      // ~leg.Draw();
      
      // ~TLatex label=gfx::cornerLabel(phiRegion,1);
      // ~label.Draw();
      saver.save(can,"Compare_Samples/ScaleFactor_"+phiRegion,true,true);
   }
   
   //Plot ptnunu/met and ptnunu/puppimet for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      TCanvas can;
      can.cd();
      
      TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/METdiff_vs_MET_profile");
      TH1D* SFgen_profile = histReader.read<TH1D>("binningUnfolding/METdiff_vs_PuppiMet_profile");
      
      SF_profile->SetLineColor(kBlue);
      SFgen_profile->SetLineColor(kRed);
      SF_profile->SetMarkerColor(kBlue);
      SFgen_profile->SetMarkerColor(kRed);
      
      SF_profile->SetStats(0);
      SF_profile->GetYaxis()->SetRangeUser(0.5,SF_profile->GetMaximum());
      SF_profile->GetXaxis()->SetTitle("p_{T}^{miss} or PuppiMet (GeV)");
      SF_profile->Draw("hist");
      SFgen_profile->Draw("same hist");
      
      gfx::LegendEntries le;
      le.append(*SF_profile,"p_{T}^{#nu#nu}/p_{T}^{miss}","l");
      le.append(*SFgen_profile,"p_{T}^{#nu#nu}/PuppiMet","l");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel(sSample,1);
      label.Draw();
      saver.save(can,"PuppiMet/ScaleFactor_"+sSample,true,true);
   }
   
   //Plot ptnunu/met as a function of ptnuunu each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      TCanvas can;
      can.cd();
      
      TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/METdiff_vs_PtNuNu_profile");
      
      SF_profile->SetLineColor(kBlue);
      
      SF_profile->SetStats(0);
      SF_profile->SetMaximum(1.3*SF_profile->GetMaximum());
      SF_profile->GetYaxis()->SetRangeUser(0.5,SF_profile->GetMaximum());
      SF_profile->GetXaxis()->SetTitle("p_{T}^{#nu#nu}(GeV)");
      SF_profile->Draw("hist");
      
      gfx::LegendEntries le;
      le.append(*SF_profile,"p_{T}^{#nu#nu}/p_{T}^{miss}","l");
      TLegend leg=le.buildLegend(.7,.75,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel(sSample,1);
      label.Draw();
      saver.save(can,"SF_vs_pTnunu/ScaleFactor_"+sSample,true,true);
   }
   
   //Plot leading top pt as a function of dPhi(met,next lep) for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      TCanvas can;
      can.cd();
      
      TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_profile");
      
      SF_profile->SetLineColor(kBlue);
      
      SF_profile->SetStats(0);
      SF_profile->SetMaximum(1.3*SF_profile->GetMaximum());
      SF_profile->GetYaxis()->SetRangeUser(0.5,SF_profile->GetMaximum());
      SF_profile->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{#nu#nu},nearest l)");
      SF_profile->Draw("hist");
      
      gfx::LegendEntries le;
      le.append(*SF_profile,"p_{T}^{t1} (GeV)","l");
      TLegend leg=le.buildLegend(.7,.75,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel(sSample,1);
      label.Draw();
      saver.save(can,"leadTopPt_vs_dPhiMetLep/trend_"+sSample,true,true);
   }
   
   //Plot 2D ptnunu/met vs met each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      for (TString phiRegion : {"","_phi1","_phi2","_phi3"}){
         TCanvas can;
         can.cd();
         can.SetLogz();
         
         TH2D* SF_profile = histReader.read<TH2D>("binningUnfolding/METdiff"+phiRegion+"_vs_MET");
         gPad->SetRightMargin(0.15);
         SF_profile->SetStats(0);
         SF_profile->GetYaxis()->SetRangeUser(0,3);
         SF_profile->GetYaxis()->SetTitle("p_{T}^{miss}/p_{T}^{#nu#nu}");
         SF_profile->GetYaxis()->SetTitleOffset(0.6);
         SF_profile->SetMinimum(0.0001);
         SF_profile->Draw("hcolz");
         
         TLatex label=gfx::cornerLabel(phiRegion,2);
         label.Draw();
         saver.save(can,"ScatterPlots/ScaleFactor_"+sSample+phiRegion,true,true);
      }
   }
   
   //Plot 2D dPhi(ptnunu,met) vs met for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      for (TString phiRegion : {"_phi1","_phi2","_phi3"}){
         TCanvas can;
         can.cd();
         can.SetLogz();
         
         TH2D* SF_profile = histReader.read<TH2D>("binningUnfolding/dPhiMETpTnunu"+phiRegion+"_vs_MET");
         gPad->SetRightMargin(0.15);
         SF_profile->SetStats(0);
         SF_profile->GetYaxis()->SetRangeUser(0,3);
         SF_profile->GetYaxis()->SetTitle("|#Delta#phi|(p_{T}^{miss},p_{T}^{#nu#nu})");
         SF_profile->GetYaxis()->SetTitleOffset(0.6);
         SF_profile->SetMinimum(0.0001);
         SF_profile->Draw("hcolz");
         
         TLatex label=gfx::cornerLabel(phiRegion,2);
         label.Draw();
         saver.save(can,"ScatterPlots/dPhiMETpTnunu"+sSample+phiRegion,true,true);
      }
   }
   
   //Plot dPhi(ptnunu,met) for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      for (TString phiRegion : {"_phi1","_phi2","_phi3"}){
         TCanvas can;
         can.cd();
         
         TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/dPhiMETpTnunu"+phiRegion+"_vs_MET_profile");
         
         SF_profile->SetLineColor(kBlue);
         SF_profile->SetMarkerColor(kBlue);
         
         SF_profile->SetStats(0);
         SF_profile->Draw("hist");
         
         gfx::LegendEntries le;
         le.append(*SF_profile,"|#Delta#phi|(p_{T}^{miss},p_{T}^{#nu#nu})","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(phiRegion,1);
         label.Draw();
         saver.save(can,"dPhiMETpTnunu_"+sSample+phiRegion,true,true);
      }
   }
   
}
