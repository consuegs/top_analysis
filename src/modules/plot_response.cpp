//Script to plot response Matrix

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
#include <TProfile.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_reponse");
   io::RootFileReader histReader(TString::Format("TUnfold%.1f.root",cfg.processFraction*100),"");
   TCanvas can;
   
   TH2F *hist;
   TH2F *hist_alt;
   
   // ~hist = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton/histMCGenRec");
   // ~hist_alt = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_MadGraph/histMCGenRec");
   // ~hist = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton/histMCGenRec_sameBins");
   // ~hist_alt = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_MadGraph/histMCGenRec_sameBins");
   // ~hist = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton_PTreweight/histMCGenRec_sameBins");
   // ~hist_alt = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton_PTreweight_sanity/histMCGenRec_sameBins");
   // ~hist_alt = (TH2F*) histReader.read<TH2F>("TUnfold_binning__/histMCGenRec_sameBins");
   hist = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton_Puppi/histMCGenRec_sameBins");
   hist_alt = (TH2F*) histReader.read<TH2F>("TUnfold_binning_dilepton_dilepton_Puppi/histMCGenRec_sameBins");
   
   for (TString norm:{"","LineNorm"}){    //Plot nominal response and response normalized in each row
      
      if (norm==""){
         hist->Scale(1.0/(hist->Integral()));
      }
      else{
         //Normalize each individual line of diagram
         float sum;
         for (int y=1; y<=hist->GetNbinsY(); y++){
            sum=hist->Integral(1,hist->GetNbinsX(),y,y);
            if (sum==0) continue;
            for (int x=1; x<=hist->GetNbinsY(); x++){
               if (hist->GetBinContent(x,y)!=0)hist->SetBinContent(x,y,hist->GetBinContent(x,y)/sum);
               else hist->SetBinContent(x,y,0.000002);
            }
         }
      }
      
      // ~for (int y=1; y<=hist->GetNbinsY(); y++){
         // ~for (int x=1; x<=hist->GetNbinsX(); x++){
            // ~if (hist->GetBinContent(x,y)==0) hist->SetBinContent(x,y,0.000001);
         // ~}
      // ~}
      
      // Plotting stuff in the following
      can.cd();
      // ~can.SetLogz();
      
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.11);
      
      hist->GetYaxis()->SetTitleOffset(1.3);
      hist->GetXaxis()->SetTitleOffset(0.9);
      hist->GetZaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetTitleSize(0.05);
      hist->GetZaxis()->SetTitleSize(0.05);
      hist->GetYaxis()->SetLabelSize(0.04);
      hist->GetXaxis()->SetLabelSize(0.04);
      hist->GetZaxis()->SetLabelSize(0.04);
      
      hist->SetMarkerColor(kRed);
      hist->SetMarkerSize(1.5);
      hist->SetTitle("");
      
      if (norm==""){
         hist->GetZaxis()->SetTitle("normalized distribution");
         
         hist->SetMinimum(0.000001);
         hist->SetMaximum(0.1);
      }
      else{
         hist->GetZaxis()->SetTitle("line normalized distribution");
         
         hist->SetMinimum(0.000001);
         hist->SetMaximum(1.);
      }

      hist->SetStats(false);
      hist->Draw("colz text");
      // ~hist->Draw("colz");
      
      hist->GetYaxis()->SetTitle("reco binNumber");
      hist->GetXaxis()->SetTitle("gen binNumber");
      // ~hist->GetYaxis()->SetRangeUser(2,13);
      // ~hist->GetXaxis()->SetRangeUser(2,13);
      
      can.RedrawAxis();
      TString plotLoc="dilepton_dilepton/reponse_"+norm;
      saver.save(can,plotLoc,true,true);
      can.Clear();
      
      if(norm=="LineNorm"){
         float sum;
         for (int y=1; y<=hist_alt->GetNbinsY(); y++){
            sum=hist_alt->Integral(1,hist_alt->GetNbinsX(),y,y);
            if (sum==0) continue;
            for (int x=1; x<=hist_alt->GetNbinsY(); x++){
               if (hist_alt->GetBinContent(x,y)!=0)hist_alt->SetBinContent(x,y,hist_alt->GetBinContent(x,y)/sum);
               else hist_alt->SetBinContent(x,y,0.000002);
            }
         }
         saver.save(*hist_alt,"dilepton_MadGraph/response_"+norm);
         hist->Add(hist_alt,-1.);
         for (int i=0; i<=hist->GetNbinsX()+1; i++){
            for (int j=0; j<=hist->GetNbinsY()+1; j++){
               hist->SetBinContent(i,j,hist->GetBinContent(i,j)*100);
            }
         }
         hist->GetZaxis()->SetTitle("difference (*100)");
         TH1D* projectionX=hist_alt->ProjectionX();
         saver.save(*hist,"diff_dilepton_MadGraph");
         saver.save(*projectionX,"diff_dilepton_MadGraph_project");
         hist->Draw("colz text");
         saver.save(can,"diff_dilepton_MadGraph_plot",true,true);
      }
   }
}
