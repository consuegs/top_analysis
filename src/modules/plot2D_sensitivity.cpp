//Script to plot 2D sensitivity distributions

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
#include <TPie.h>
#include <TPieSlice.h>

Config const &cfg=Config::get();

void add_Categories(TString const path, io::RootFileReader const &reader_hist, TH2F &out_hist) {   //Function to add the three different categories
   TH2F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH2F*) reader_hist.read<TH2F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH2F) *(reader_hist.read<TH2F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

std::vector<float> get_Ratio(TH2F bkg, TH2F const &total) {    //Function to get the ratio of the 2D histograms
   std::vector<float> ratios;
   bkg.Divide(&total);
   for (int j=1; j<=bkg.GetNbinsY(); j++){
      for (int i=1; i<=bkg.GetNbinsX(); i++){
         ratios.push_back(bkg.GetBinContent(i,j));
      }
   }
   return ratios;
}

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100),"plot2D_sensitvity");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   TCanvas can;
   
   TH2F signal;
   TH2F ttbar;
   TH2F SMbkg;
   TH2F sensitivity;
   TH2F temp_hist;
   
   // Define different binning schemes
   std::vector<float> met_bins1={0,100,200,300,400};
   std::vector<float> phi_bins1={0,0.8,1.6,2.4,3.2};
   
   std::vector<float> met_bins2={0,50,100,200,400};
   std::vector<float> phi_bins2={0,0.4,0.8,1.6,3.2};
   
   std::vector<float> met_bins3={0,70,140,250,400};
   std::vector<float> phi_bins3={0,0.4,0.8,1.2,3.14};
   
   std::vector<float> met_bins4={0,80,160,260,400};
   std::vector<float> phi_bins4={0,0.7,1.4,3.14};
   
   //Used for names in the savin process
   int numberBinningSchemeMet=1;
   int numberBinningSchemePhi=1;
   
   for(std::vector<float> met_bins : {met_bins1,met_bins2,met_bins3,met_bins4}){
      for(std::vector<float> phi_bins : {phi_bins1,phi_bins2,phi_bins3,phi_bins4}){
               
         for (TString bkgSample :{"SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ_SM","ttH_SM"}) {     //Add all SM backgrounds except ttbar
            add_Categories("2d_MetVSdPhiMetNearLep/"+bkgSample,histReader,temp_hist);
            if (bkgSample=="SingleTop") add_Categories("2d_MetVSdPhiMetNearLep/"+bkgSample,histReader,SMbkg);
            else SMbkg.Add(&temp_hist);
         }
         
         add_Categories("2d_MetVSdPhiMetNearLep/TTbar",histReader,ttbar);
         
         SMbkg=hist::rebinned(SMbkg,met_bins,phi_bins);
         ttbar=hist::rebinned(ttbar,met_bins,phi_bins);
         
         
         for (TString signalSample :{"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
            
            add_Categories("2d_MetVSdPhiMetNearLep/"+signalSample,histReader,signal);
            signal=hist::rebinned(signal,met_bins,phi_bins);
            
            temp_hist=ttbar;
            temp_hist.Add(&SMbkg);
            temp_hist.Add(&signal);
            
            sensitivity=signal;
            sensitivity.Divide(&temp_hist);     //Get S/(S+B)
            
            //Plotting stuff in the following
            can.cd();
            can.SetLogz();
            
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.11);
            
            sensitivity.GetYaxis()->SetTitleOffset(1.3);
            sensitivity.GetXaxis()->SetTitleOffset(0.9);
            sensitivity.GetZaxis()->SetTitleOffset(1.3);
            sensitivity.GetYaxis()->SetTitleSize(0.05);
            sensitivity.GetXaxis()->SetTitleSize(0.05);
            sensitivity.GetZaxis()->SetTitleSize(0.05);
            sensitivity.GetYaxis()->SetLabelSize(0.04);
            sensitivity.GetXaxis()->SetLabelSize(0.04);
            sensitivity.GetZaxis()->SetLabelSize(0.04);
            sensitivity.GetZaxis()->SetTitle("S/(S+B)");
            
            sensitivity.SetMaximum(1.0);
            sensitivity.SetStats(false);
            sensitivity.Draw("colz");
            
            signal.SetMarkerColor(kRed);
            signal.SetMarkerSize(1.5);
            signal.Draw("same text");
            
            TString cat_label="all";
            TString labelString=cat_label+"    "+signalSample;
            TLatex label=gfx::cornerLabel(labelString,1);
            label.Draw();
            
            can.RedrawAxis();
            TString plotLoc="Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+"/"+signalSample;
            saver.save(can,plotLoc,true,true);
            can.Clear();
         }
         
         // BKG composition for individual bins
         int nTotBins=(met_bins.size()-1)*(phi_bins.size()-1);
         Float_t bin_ratios[nTotBins][7];
         int j=0;
         for (TString bkgSample :{"SingleTop","DrellYan","ttZ_SM","ttH_SM","WJetsToLNu","Diboson","TTbar"}) {
            TH2F temp;
            if (bkgSample=="Diboson"){    //If diboson add samples
               TH2F temp2;
               for (TString dibosonSample :{"WW","WZ","ZZ"}) {
                  if (dibosonSample=="WW") {
                     add_Categories("2d_MetVSdPhiMetNearLep/"+dibosonSample,histReader,temp);
                  }
                  else {
                     add_Categories("2d_MetVSdPhiMetNearLep/"+dibosonSample,histReader,temp2);
                     temp.Add(&temp2);
                  }
               }
            }
            else add_Categories("2d_MetVSdPhiMetNearLep/"+bkgSample,histReader,temp);
            
            temp=hist::rebinned(temp,met_bins,phi_bins);    //Rebin following the current binning scheme
            
            
            std::vector<float> ratios=get_Ratio(temp,SMbkg);      //Get ratios needed for the pie charts
            for (int i=0;i<int(ratios.size());i++){
               bin_ratios[i][j]=ratios[i];
            }
            j++;
         }
         
         //Plot pie charts with and without ttbar
         TCanvas can;
         TCanvas can_noTTbar;
         can.Divide(met_bins.size()-1,phi_bins.size()-1);
         can_noTTbar.Divide(met_bins.size()-1,phi_bins.size()-1);
         Int_t colors[] = {2,3,4,5,6,7,8};
         Int_t fix[] ={-3,3,1,-1};
         for (int i=0;i<nTotBins;i++){
            Float_t vals[7];
            for (int j=0;j<7;j++){
               vals[j]= bin_ratios[i][j];
            }
            TString binName="Bin"+std::to_string(i+1);
            TPie *pie1 = new TPie("pie1",binName,7,vals,colors);
            
            pie1->GetSlice(0)->SetTitle("SingleTop");
            pie1->GetSlice(1)->SetTitle("DrellYan");
            pie1->GetSlice(2)->SetTitle("ttZ");
            pie1->GetSlice(3)->SetTitle("ttH");
            pie1->GetSlice(4)->SetTitle("W+Jets");
            pie1->GetSlice(5)->SetTitle("Diboson");
            pie1->GetSlice(6)->SetTitle("TTbar");
            pie1->SetLabelsOffset(-.2);
            
            can.cd(nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)]);
            
            pie1->Draw("nol R <");
            if ((nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)])==1) pie1->MakeLegend();
            
            can_noTTbar.cd(nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)]);
            TPie *pie2=(TPie*) pie1->Clone("pie2");
            pie2->GetSlice(6)->SetTitle("");
            pie2->GetSlice(6)->SetValue(0);
            pie2->Draw("nol R <");
            if ((nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)])==1) pie2->MakeLegend();
         }
         // Save pie charts
         TString plotLoc="../output/BKGratios/BinningMet_"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+".pdf";
         can.SaveAs(plotLoc);    //Save as pdf right away since there is not proper vizualisation in TBrwoser for PieCharts
         plotLoc="../output/BKGratios/BinningMet_"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+"_noTTbar.pdf";
         can_noTTbar.SaveAs(plotLoc);
         
         numberBinningSchemePhi++;
      }
      numberBinningSchemeMet++;
      numberBinningSchemePhi=1;
   }
      
}
