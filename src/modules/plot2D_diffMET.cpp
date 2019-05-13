//Script to plot 2D distributions for met, genMet and pT of both neutrinos

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_diffMET");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   TCanvas can;
   
   
   for (TString sSample :{"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
   "T1tttt_1200_800","T2tt_650_250","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
   
      for (TString sVar :{"2d_nunuVSgenMet","MetVSgenMet","2d_nunuVSMet","2d_nunuVSDMgenMet","MetVSDMgenMet"}) {
         
         TH2F hist_add;
         TH2F *hist;
         
         for (TString cat:{"ee","emu","mumu","all"}){
            
            if (cat!="all"){
               hist = (TH2F*) histReader.read<TH2F>("genParticles/"+cat+"/"+sVar+"/"+sSample);
               if (cat=="ee") hist_add=(TH2F) *(histReader.read<TH2F>("genParticles/"+cat+"/"+sVar+"/"+sSample));
               else hist_add.Add(hist);
            }
            else {
               hist = &hist_add;
            }
               
            
            hist->RebinX(2);
            hist->RebinY(2);
            
            for (TString norm:{"","LineNorm"}){
               
               if (norm==""){
                  hist->Scale(1.0/(hist->Integral()));
               }
               else{
                  //Normalize each individual line of diagram
                  float sum;
                  for (int y=1; y<=hist->GetNbinsY(); y++){
                     sum=hist->Integral(1,hist->GetNbinsX(),y,y);
                     if (sum==0) continue;
                     for (int x=1; x<=hist->GetNbinsX(); x++){
                        hist->SetBinContent(x,y,hist->GetBinContent(x,y)/sum);
                     }
                  }
               }
               
               TGraph diag;
               diag.SetPoint(1,0,0);
               diag.SetPoint(2,600,600);
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
                  
                  hist->SetMinimum(0.00001);
                  hist->SetMaximum(1.);
               }

               hist->SetStats(false);
               hist->Draw("colz");
               diag.Draw("same");
               
               TString cat_label=cat;
               if (cat=="emu") cat_label="e#mu";
               else if (cat=="mumu") cat_label="#mu#mu";
               else if (cat=="all") cat_label="all";
               TLatex label=gfx::cornerLabel(cat_label+"    "+sSample,1);
               label.Draw();
               
               can.RedrawAxis();
               TString plotLoc=sSample+"/"+sVar+norm+"/"+cat;
               saver.save(can,plotLoc,true,true);
               can.Clear();
            }
         }
         
      }
   }
}
