//Script to plot difference between met, genMet and pT of both neutrinos

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_diffMET");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   TCanvas can;
   
   
   for (TString sSample :{"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
   "T1tttt_1200_800","T2tt_650_250","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
   
      for (TString sVar :{"diff_ptNuNu_genMET","diff_Met_genMET","diff_ptNuNu_Met","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_ptNuNu_DMgenMET","diff_Met_DMgenMET"}) {
         
         TH1F hist_add;
         TH1F *hist;
         
         for (TString cat:{"ee","emu","mumu","all"}){
            
            if (cat!="all"){
               hist=(TH1F*) histReader.read<TH2F>("genParticles/"+cat+"/"+sVar+"/"+sSample);
               if (cat=="ee") hist_add=(TH1F) *(histReader.read<TH1F>("genParticles/"+cat+"/"+sVar+"/"+sSample));
               else hist_add.Add(hist);
            }
            else {
               hist = &hist_add;
            }
            
            if (sVar!="diff_ptNuNu_genMET" and sVar!="diff_Met_genMET_norm" and sVar!="diff_Met_DMgenMET_norm") hist->RebinX(4);
            //~ else hist->RebinX(2);
            
            hist->Scale(1.0/(hist->Integral()));
            
            can.cd();
            can.SetLogz();
            
            hist->GetYaxis()->SetTitle("normalized distribution");

            hist->Fit("gaus","Q");
            
            if (sVar=="diff_Met_genMET_norm"){
               hist->GetXaxis()->SetTitle("(p_{T}^{miss}-genMET)/genMET");
               hist->GetXaxis()->SetRangeUser(-5,5);
            }
            else if (sVar=="diff_Met_DMgenMET_norm"){
               hist->GetXaxis()->SetTitle("(p_{T}^{miss}-DMgenMET)/DMgenMET");
               hist->GetXaxis()->SetRangeUser(-5,5);
            }
            gStyle->SetStatY(0.9);
            hist->Draw();
            
            TString cat_label=cat;
            if (cat=="emu") cat_label="e#mu";
            else if (cat=="mumu") cat_label="#mu#mu";
            else if (cat=="all") cat_label="all";
            TLatex label=gfx::cornerLabel(cat_label+"    "+sSample,1);
            label.Draw();
            
            can.RedrawAxis();
            TString plotLoc=sSample+"/"+sVar+"/"+cat;
            saver.save(can,plotLoc,true,true);
            can.Clear();
         }
      }
   }
}
