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
   
   //Plot leading top pt as a function of dPhi(met,next lep) for each sample in met bins
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

      TCanvas can;
      can.cd();
      
      TH1D* SF_profile = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_profile");
      TH1D* SF_profile_met1 = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_met1_profile");
      TH1D* SF_profile_met2 = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_met2_profile");
      TH1D* SF_profile_met3 = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_met3_profile");
      TH1D* SF_profile_met4 = histReader.read<TH1D>("binningUnfolding/TopPt_vs_dPhi_met4_profile");
      
      SF_profile->SetLineColor(kBlue);
      SF_profile_met1->SetLineColor(kRed);
      SF_profile_met2->SetLineColor(kMagenta);
      SF_profile_met3->SetLineColor(kGreen+3);
      SF_profile_met4->SetLineColor(kOrange-3);
      
      SF_profile->SetMarkerSize(0);
      SF_profile_met1->SetMarkerSize(0);
      SF_profile_met2->SetMarkerSize(0);
      SF_profile_met3->SetMarkerSize(0);
      SF_profile_met4->SetMarkerSize(0);
      
      SF_profile_met4->SetStats(0);
      SF_profile_met4->SetMaximum(1.3*SF_profile_met4->GetMaximum());
      SF_profile_met4->GetYaxis()->SetRangeUser(0.5,SF_profile_met4->GetMaximum());
      SF_profile_met4->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{#nu#nu},nearest l)");
      SF_profile_met4->GetYaxis()->SetTitle("p_{T}^{t1} (GeV)");
      
      SF_profile_met4->Draw("e1");
      SF_profile_met1->Draw("e1 same");
      SF_profile_met2->Draw("e1 same");
      SF_profile_met3->Draw("e1 same");
      SF_profile->Draw("e1 same");
      
      gfx::LegendEntries le;
      le.append(*SF_profile,"all","l");
      le.append(*SF_profile_met1,"met1","l");
      le.append(*SF_profile_met2,"met2","l");
      le.append(*SF_profile_met3,"met3","l");
      le.append(*SF_profile_met4,"met4","l");
      TLegend leg=le.buildLegend(.6,.65,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      saver.save(can,"leadTopPt_vs_dPhiMetLep/trend_"+sSample+"_combined",true,true);
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
         SF_profile->GetYaxis()->SetTitle("p_{T}^{#nu#nu}/p_{T}^{miss}");
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
   
   //Plot Profile2D lead top pt for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

         TCanvas can;
         can.cd();
         
         TH2D* SF_profile = histReader.read<TH2D>("binningUnfolding/TopPt_profile");
         gPad->SetRightMargin(0.2);
         SF_profile->SetStats(0);
         SF_profile->GetYaxis()->SetRangeUser(0,3);
         SF_profile->GetYaxis()->SetTitleOffset(0.6);
         SF_profile->GetZaxis()->SetTitleOffset(1.2);
         SF_profile->GetZaxis()->SetLabelOffset(0.01);
         SF_profile->GetZaxis()->SetTitle("mean p_{T}^{t1} (GeV)");
         SF_profile->SetMarkerColor(kRed);
         SF_profile->Draw("hcolz text");
         
         TLatex label=gfx::cornerLabel(sSample,2);
         label.Draw();
         saver.save(can,"Profile_leadTopPt/Profile_"+sSample,true,true);
   }
   
   //Plot SF vs. dPhiLepLep or dPhiNuNu for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      for (TString var:{"dPhiLep","dPhiNuNu"}) {
         io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));

         TCanvas can;
         can.cd();
         
         TH1D* SF_profile_phi1 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_phi1_profile");
         TH1D* SF_profile_phi2 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_phi2_profile");
         TH1D* SF_profile_phi3 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_phi3_profile");
         
         SF_profile_phi1->SetLineColor(kBlue);
         SF_profile_phi2->SetLineColor(kRed);
         SF_profile_phi3->SetLineColor(kGreen+3);
         
         SF_profile_phi1->SetMarkerSize(0);
         SF_profile_phi2->SetMarkerSize(0);
         SF_profile_phi3->SetMarkerSize(0);
         
         if (var=="dPhiLep"){
            SF_profile_phi1->GetXaxis()->SetTitle("|#Delta#phi|(l1,l2)");
            SF_profile_phi1->SetMaximum(1.0);
            SF_profile_phi1->SetMinimum(0.8);
         }
         else {
            SF_profile_phi1->GetXaxis()->SetTitle("|#Delta#phi|(#nu,#nu)");
            SF_profile_phi1->SetMaximum(1.2);
            SF_profile_phi1->SetMinimum(0.6);
         }
         gPad->SetLeftMargin(0.13);
         SF_profile_phi1->GetYaxis()->SetTitleOffset(1.0);
         SF_profile_phi1->GetYaxis()->SetTitle("p_{T}^{#nu#nu}/p_{T}^{miss}");
         
         SF_profile_phi1->SetStats(0);
         SF_profile_phi1->Draw("e1");
         SF_profile_phi2->Draw("e1 same");
         SF_profile_phi3->Draw("e1 same");
         
         gfx::LegendEntries le;
         le.append(*SF_profile_phi1,"phi1","l");
         le.append(*SF_profile_phi2,"phi2","l");
         le.append(*SF_profile_phi3,"phi3","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(sSample,1);
         label.Draw();
         saver.save(can,"SF_"+var+"/combined_"+sSample,true,true);
      }
   }
   
   //Plot SF vs. dPhiLepLep or dPhiNuNu for each sample met bins
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      for (TString var:{"dPhiLep","dPhiNuNu"}) {
         TCanvas can;
         can.cd();
         
         TH1D* SF_profile_met1 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_met1_profile");
         TH1D* SF_profile_met2 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_met2_profile");
         TH1D* SF_profile_met3 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_met3_profile");
         TH1D* SF_profile_met4 = histReader.read<TH1D>("binningUnfolding/METdiff_"+var+"_met4_profile");
         
         SF_profile_met1->SetLineColor(kBlue);
         SF_profile_met2->SetLineColor(kRed);
         SF_profile_met3->SetLineColor(kGreen+3);
         SF_profile_met4->SetLineColor(kOrange-3);
         
         SF_profile_met1->SetMarkerSize(0);
         SF_profile_met2->SetMarkerSize(0);
         SF_profile_met3->SetMarkerSize(0);
         SF_profile_met4->SetMarkerSize(0);
         
         if (var=="dPhiLep"){
            SF_profile_met1->GetXaxis()->SetTitle("|#Delta#phi|(l1,l2)");
            SF_profile_met1->SetMaximum(1.5);
            SF_profile_met1->SetMinimum(0.5);
         }
         else {
            SF_profile_met1->GetXaxis()->SetTitle("|#Delta#phi|(#nu,#nu)");
            // ~SF_profile_met1->SetMaximum(1.2);
            SF_profile_met1->SetMaximum(2.);
            SF_profile_met1->SetMinimum(0.6);
         }
         
         gPad->SetLeftMargin(0.13);
         SF_profile_met1->GetYaxis()->SetTitleOffset(1.0);
         SF_profile_met1->GetYaxis()->SetTitle("p_{T}^{#nu#nu}/p_{T}^{miss}");
         
         SF_profile_met1->SetStats(0);
         SF_profile_met1->Draw("e1");
         SF_profile_met2->Draw("e1 same");
         SF_profile_met3->Draw("e1 same");
         SF_profile_met4->Draw("e1 same");
         
         gfx::LegendEntries le;
         le.append(*SF_profile_met1,"met1","l");
         le.append(*SF_profile_met2,"met2","l");
         le.append(*SF_profile_met3,"met3","l");
         le.append(*SF_profile_met4,"met4","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(sSample,1);
         label.Draw();
         saver.save(can,"SF_"+var+"/combined_metBins_"+sSample,true,true);
      }
   }
   
   //Plot relDiffMet vs. dPhiLepNu for met bins
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      for (TString var:{"PtNuNuDiffMETRel_dPhiMETLep","PtNuNuDiffGenMETRel_dPhiMETLep","GenMetDiffMETRel_dPhiMETLep"}) {
         TCanvas can;
         can.cd();
         
         TH1D* SF_profile_met1 = histReader.read<TH1D>("binningUnfolding/"+var+"_met1_profile");
         TH1D* SF_profile_met2 = histReader.read<TH1D>("binningUnfolding/"+var+"_met2_profile");
         TH1D* SF_profile_met3 = histReader.read<TH1D>("binningUnfolding/"+var+"_met3_profile");
         TH1D* SF_profile_met4 = histReader.read<TH1D>("binningUnfolding/"+var+"_met4_profile");
         
         SF_profile_met1->SetLineColor(kBlue);
         SF_profile_met2->SetLineColor(kRed);
         SF_profile_met3->SetLineColor(kGreen+3);
         SF_profile_met4->SetLineColor(kOrange-3);
         
         SF_profile_met1->SetMarkerSize(0);
         SF_profile_met2->SetMarkerSize(0);
         SF_profile_met3->SetMarkerSize(0);
         SF_profile_met4->SetMarkerSize(0);
         
         if (var=="PtNuNuDiffMETRel_dPhiMETLep"){
            SF_profile_met1->GetYaxis()->SetTitle("(p_{T}^{#nu#nu}-p_{T}^{miss})/genMET");
            SF_profile_met1->SetMaximum(0.5);
            SF_profile_met1->SetMinimum(-0.7);
         }
         else if (var=="PtNuNuDiffGenMETRel_dPhiMETLep"){
            SF_profile_met1->GetYaxis()->SetTitle("(p_{T}^{#nu#nu}-genMet)/genMET");
            SF_profile_met1->SetMaximum(0.5);
            SF_profile_met1->SetMinimum(-0.5);
         }
         else {
            SF_profile_met1->GetYaxis()->SetTitle("(genMet-p_{T}^{miss})/genMET");
            SF_profile_met1->SetMaximum(0.5);
            SF_profile_met1->SetMinimum(-0.5);
         }

         
         gPad->SetLeftMargin(0.13);
         SF_profile_met1->GetYaxis()->SetTitleOffset(1.0);
         SF_profile_met1->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{miss},nearest l)");
         
         SF_profile_met1->SetStats(0);
         SF_profile_met1->Draw("e1");
         SF_profile_met2->Draw("e1 same");
         SF_profile_met3->Draw("e1 same");
         SF_profile_met4->Draw("e1 same");
         
         gfx::LegendEntries le;
         le.append(*SF_profile_met1,"met1","l");
         le.append(*SF_profile_met2,"met2","l");
         le.append(*SF_profile_met3,"met3","l");
         le.append(*SF_profile_met4,"met4","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(sSample,1);
         label.Draw();
         saver.save(can,"RelDiff/"+var+"_"+sSample,true,true);
      }
   }
   
   //Plot relDiffMet vs. dPhiLepNu for met>120
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      TCanvas can;
      can.cd();
      
      TH1D* SF_profile_met1 = histReader.read<TH1D>("binningUnfolding/GenMetDiffMETRel_dPhiMETLep_met120_profile");
      
      SF_profile_met1->SetLineColor(kBlue);
      
      SF_profile_met1->SetMarkerSize(0);
      
      SF_profile_met1->GetYaxis()->SetTitle("mean[(genMET-p_{T}^{miss})/genMET]");

      gPad->SetLeftMargin(0.13);
      SF_profile_met1->GetYaxis()->SetTitleOffset(1.0);
      SF_profile_met1->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{miss},nearest l)");
      
      SF_profile_met1->SetStats(0);
      SF_profile_met1->Draw("e1");
      
      // ~gfx::LegendEntries le;
      // ~le.append(*SF_profile_met1,"met1","l");
      // ~le.append(*SF_profile_met2,"met2","l");
      // ~le.append(*SF_profile_met3,"met3","l");
      // ~le.append(*SF_profile_met4,"met4","l");
      // ~TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      // ~leg.Draw();
      
      // ~TLatex label=gfx::cornerLabel(sSample,1);
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>120 GeV",2);
      label.Draw();
      saver.save(can,"RelDiff_met120/metRes_"+sSample,true,true);
   }
   
   //Plot relDiffMet vs. dPhiLepNu for met>120 compare puppi and pfmet
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      TCanvas can;
      can.cd();
      
      TH1D* SF_profile_pfmet = histReader.read<TH1D>("binningUnfolding/GenMetDiffMETRel_dPhiMETLep_met120_profile");
      TH1D* SF_profile_puppi = histReader.read<TH1D>("binningUnfolding/GenMetDiffPuppiMETRel_dPhiMETLep_met120_profile");
      
      SF_profile_pfmet->SetLineColor(kBlue);
      SF_profile_puppi->SetLineColor(kRed);
      
      SF_profile_pfmet->SetMarkerSize(0);
      SF_profile_puppi->SetMarkerSize(0);
      
      SF_profile_pfmet->GetYaxis()->SetTitle("mean[(genMET-p_{T}^{miss})/genMET]");

      gPad->SetLeftMargin(0.13);
      SF_profile_pfmet->GetYaxis()->SetTitleOffset(1.0);
      SF_profile_pfmet->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{miss},nearest l)");
      
      SF_profile_pfmet->SetStats(0);
      SF_profile_pfmet->SetMaximum(-0.05);
      SF_profile_pfmet->Draw("e1");
      SF_profile_puppi->Draw("same e1");
      
      gfx::LegendEntries le;
      le.append(*SF_profile_pfmet,"pfMET","l");
      le.append(*SF_profile_puppi,"Puppi","l");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      // ~TLatex label=gfx::cornerLabel("p_{T}^{miss}>120 GeV",2);
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>180 GeV",2);
      label.Draw();
      saver.save(can,"RelDiff_met120/compare_PuppiPF"+sSample,true,true);
   }
   
   //Plot metSig vs. dPhiLepNu for met>120
   for (TString sSample :{"","dilepton","MadGraph"}){
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      TCanvas can;
      can.cd();
      
      TH1D* SF_profile_met1 = histReader.read<TH1D>("binningUnfolding/METsig_dPhiMETLep_met120_profile");
      
      SF_profile_met1->SetLineColor(kBlue);
      
      SF_profile_met1->SetMarkerSize(0);
      
      SF_profile_met1->GetYaxis()->SetTitle("metSig");

      gPad->SetLeftMargin(0.13);
      SF_profile_met1->GetYaxis()->SetTitleOffset(1.0);
      SF_profile_met1->GetXaxis()->SetTitle("|#Delta#phi|(p_{T}^{miss},nearest l)");
      
      SF_profile_met1->SetStats(0);
      SF_profile_met1->Draw("e1");
      
      // ~gfx::LegendEntries le;
      // ~le.append(*SF_profile_met1,"met1","l");
      // ~le.append(*SF_profile_met2,"met2","l");
      // ~le.append(*SF_profile_met3,"met3","l");
      // ~le.append(*SF_profile_met4,"met4","l");
      // ~TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      // ~leg.Draw();
      
      // ~TLatex label=gfx::cornerLabel(sSample,1);
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>120 GeV",2);
      label.Draw();
      saver.save(can,"MetSig_met120/MetSig_"+sSample,true,true);
   }
   
}
