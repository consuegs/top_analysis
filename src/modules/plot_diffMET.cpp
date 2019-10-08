//Script to plot difference between met, genMet and pT of both neutrinos(+BSM)

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
   
   for (TString sSelection : {"genParticles","genParticles_Met200"}){
      for (TString sSample :{"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      "T1tttt_1200_800","T2tt_650_350","T1tttt_1500_100","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","TTbar_diLepton"}){
      
         for (TString sVar :{"diff_ptNuNu_genMET","diff_Met_genMET","diff_ptNuNu_Met","diff_Met_genMET_norm","diff_Met_DMgenMET_norm","diff_Met_genMET_normSUM",
            "diff_Met_DMgenMET_normSUM","diff_ptNuNu_DMgenMET","diff_Met_DMgenMET","diff_genMT2_MT2","diff_genMT2neutrino_MT2","diff_dPhiMetNearLep_gen","diff_dPhiMetNearLep"}) {
            
            TH1F hist_add;
            TH1F *hist;
            
            for (TString cat:{"ee","emu","mumu","all"}){
               
               if (cat!="all"){     //Add categories
                  hist=(TH1F*) histReader.read<TH2F>(sSelection+"/"+cat+"/"+sVar+"/"+sSample);
                  if (cat=="ee") hist_add=(TH1F) *(histReader.read<TH1F>(sSelection+"/"+cat+"/"+sVar+"/"+sSample));
                  else hist_add.Add(hist);
               }
               else {
                  hist = &hist_add;
               }
               
               //Rebin in case of poor statistics
               if (sVar!="diff_ptNuNu_genMET" and sVar!="diff_Met_genMET_norm" and sVar!="diff_Met_DMgenMET_norm" and sVar!="diff_Met_genMET_normSUM" and sVar!="diff_Met_DMgenMET_normSUM") hist->RebinX(4);
               //~ else hist->RebinX(2);
               
               hist->Scale(1.0/(hist->Integral()));   //Normalize the hist to the integral
               
               can.cd();
               can.SetLogz();
               
               hist->GetYaxis()->SetTitle("normalized distribution");

               hist->Fit("gaus","Q");     //Try to fit a gaussian to the response
               
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
               TString labelString=cat_label+"    "+sSample;
               if (sSelection=="genParticles_Met200") labelString=cat_label+"    "+sSample+"    p_{T}^{miss}>200 GeV";
               TLatex label=gfx::cornerLabel(labelString,1);
               label.Draw();
               
               can.RedrawAxis();
               TString plotLoc=sSample+"/"+sVar+"/"+cat;
               if (sSelection=="genParticles_Met200") plotLoc="Met200/"+sSample+"/"+sVar+"/"+cat;
               saver.save(can,plotLoc,true,true);
               can.Clear();
            }
         }
      }
   }
}
