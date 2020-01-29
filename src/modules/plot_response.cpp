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

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_reponse");
   io::RootFileReader histReader(TString::Format("binningUnfolding_dilepton%.1f.root",cfg.processFraction*100),"binningUnfolding");
   TCanvas can;
   
   TH2F *hist;
   
   //~hist = (TH2F*) histReader.read<TH2F>("Binning_Met1_Phi1/Unfolding/response_sameBins");
   hist = (TH2F*) histReader.read<TH2F>("Binning_Met1_Phi1/Unfolding/response");
   
   for (TString norm:{"","LineNorm"}){    //Plot nominal response and response normalized in each row
      
      if (norm==""){
         hist->Scale(1.0/(hist->Integral()));
      }
      else{
         //Normalize each individual row of diagram
         float sum;
         for (int y=1; y<=hist->GetNbinsY(); y++){
            sum=hist->Integral(1,hist->GetNbinsX(),y,y);
            if (sum==0) continue;
            for (int x=1; x<=hist->GetNbinsX(); x++){
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
      can.SetLogz();
      
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
      hist->Draw("colz");
      
      hist->GetXaxis()->SetTitle("reco binNumber");
      hist->GetYaxis()->SetTitle("gen binNumber");
      
      can.RedrawAxis();
      TString plotLoc="Binning_Met1_Phi1/reponse_"+norm;
      saver.save(can,plotLoc,true,true);
      can.Clear();
   }
}
