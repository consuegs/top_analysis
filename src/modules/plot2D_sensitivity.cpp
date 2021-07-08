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

TH2F get_sqrt(TH2F inital) {    //Function to get the square-root of the 2D histogram
   TH2F sqrts=inital;
   for (int j=1; j<=inital.GetNbinsY(); j++){
      for (int i=1; i<=inital.GetNbinsX(); i++){
         sqrts.SetBinContent(i,j,sqrt(inital.GetBinContent(i,j)));
      }
   }
   return sqrts;
}

extern "C"
void run()
{
   std::cout<<"---------------------------------------"<<std::endl;
   std::cout<<"Yields are manually scaled to ful Run2 lumi (check code!!!)"<<std::endl;
   std::cout<<"---------------------------------------"<<std::endl;
   
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_sensitivity");
   // ~io::RootFileSaver saver(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100),"plot2D_sensitivity_0.91");
   // ~io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader(TString::Format("multiHists/histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   // ~io::RootFileReader histReader("histograms_v10.root",TString::Format("distributions%.1f",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("histograms_SF0.910000_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   TCanvas can;
   
   TH2F signal;
   TH2F ttbar;
   TH2F SMbkg;
   TH2F sensitivity;
   TH2F sensitivity_sqrtB;
   TH2F temp_hist;
   
   // Define different binning schemes
   std::vector<float> met_bins1={0,100,200,300,400};
   std::vector<float> phi_bins1={0,0.8,1.6,2.4,3.2};
   
   std::vector<float> met_bins2={0,50,100,200,400};
   std::vector<float> phi_bins2={0,0.4,0.8,1.6,3.2};
   
   std::vector<float> met_bins3={0,70,140,250,400};
   std::vector<float> phi_bins3={0,0.4,0.8,1.2,3.14};
   
   // ~std::vector<float> met_bins4={0,40,80,120,160,230,400};
   // ~std::vector<float> phi_bins4={0,0.7,1.4,3.141};
   std::vector<float> met_bins4={0,50,70,100,130,160,200,400,410};
   std::vector<float> phi_bins4={0,0.64,1.2,3.2};
   // ~std::vector<float> met_bins4={0,400};
   // ~std::vector<float> phi_bins4={0,3.14};
   
   //Used for names in the savin process
   int numberBinningSchemeMet=1;
   int numberBinningSchemePhi=1;
   
   //Select base 2D hist
   TString hist_2d="2d_MetVSdPhiMetNearLep_DNN/";
   // ~TString hist_2d="2d_MetVSdPhiMetNearLep_Puppi/";
   
   // ~for(std::vector<float> met_bins : {met_bins1,met_bins2,met_bins3,met_bins4}){
      // ~for(std::vector<float> phi_bins : {phi_bins1,phi_bins2,phi_bins3,phi_bins4}){
   for(std::vector<float> met_bins : {met_bins4}){
      for(std::vector<float> phi_bins : {phi_bins4}){
               
         // ~for (TString bkgSample :{"SingleTop","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW"}) {     //Add all SM backgrounds except ttbar
         for (TString bkgSample :{"SingleTop","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","ttG"}) {     //Add all SM backgrounds except ttbar
            add_Categories(hist_2d+bkgSample,histReader,temp_hist);
            if (bkgSample=="SingleTop") add_Categories(hist_2d+bkgSample,histReader,SMbkg);
            else SMbkg.Add(&temp_hist);
         }
         
         // ~add_Categories("2d_MetVSdPhiMetNearLep/TTbar",histReader,ttbar);
         // ~add_Categories("2d_MetVSdPhiMetNearLep/TTbar_diLepton",histReader,ttbar);
         add_Categories(hist_2d+"TTbar_diLepton",histReader,ttbar);
         
         std::cout<<"Correlation in ttbar: "<<ttbar.GetCorrelationFactor()<<std::endl;
         
         // ~SMbkg=hist::rebinned(SMbkg,met_bins,phi_bins);
         // ~ttbar=hist::rebinned(ttbar,met_bins,phi_bins);
         SMbkg=hist::rebinned(SMbkg,met_bins,phi_bins,false);
         ttbar=hist::rebinned(ttbar,met_bins,phi_bins,false);
         
         // ~SMbkg.Scale(137191.0/35867.05998);
         // ~ttbar.Scale(137191.0/35867.05998);
         
         
         for (TString signalSample :{"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
            
            add_Categories(hist_2d+signalSample,histReader,signal);
            signal=hist::rebinned(signal,met_bins,phi_bins);
            // ~signal.Scale(137191.0/35867.05998);
            
            for (TString sensitivity_type : {"SplusB","sqrtB"}){
               temp_hist=ttbar;
               temp_hist.Add(&SMbkg);
               if (sensitivity_type=="SplusB") temp_hist.Add(&signal);
               else if (sensitivity_type=="sqrtB") temp_hist=get_sqrt(temp_hist);
               
               sensitivity=signal;
               sensitivity.Divide(&temp_hist);
               
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
               sensitivity.GetZaxis()->SetLabelOffset(0.02);
               if (sensitivity_type=="SplusB") sensitivity.GetZaxis()->SetTitle("S/(S+B)");
               else if (sensitivity_type=="sqrtB") sensitivity.GetZaxis()->SetTitle("S/#sqrt{B}");
               
               sensitivity.SetMaximum(2.0);
               sensitivity.SetStats(false);
               sensitivity.SetMarkerColor(kRed);
               sensitivity.SetMarkerSize(2.5);
               sensitivity.Draw("colz text");
               // ~sensitivity.Draw("colz");
               
               // ~signal.SetMarkerColor(kRed);
               // ~signal.SetMarkerSize(2.5);
               // ~signal.Draw("same text");
               
               TString cat_label="all";
               TString labelString=cat_label+"    "+signalSample;
               TLatex label=gfx::cornerLabel(labelString,1);
               label.Draw();
               
               can.RedrawAxis();
               TString plotLoc="Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+"/"+sensitivity_type+"/"+signalSample;
               saver.save(can,plotLoc,true,true);
               can.Clear();
            }
         }
         
         // BKG composition for individual bins
         int nTotBins=(met_bins.size()-1)*(phi_bins.size()-1);
         Float_t bin_ratios[nTotBins][7];
         int j=0;
         SMbkg.Add(&ttbar);
         for (TString bkgSample :{"SingleTop","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","TTbar_diLepton"}) {
            TH2F temp;
            if (bkgSample=="Diboson"){    //If diboson add samples
               TH2F temp2;
               for (TString dibosonSample :{"WW","WZ","ZZ"}) {
                  if (dibosonSample=="WW") {
                     add_Categories(hist_2d+dibosonSample,histReader,temp);
                  }
                  else {
                     add_Categories(hist_2d+dibosonSample,histReader,temp2);
                     temp.Add(&temp2);
                  }
               }
            }
            else add_Categories(hist_2d+bkgSample,histReader,temp);
            
            temp=hist::rebinned(temp,met_bins,phi_bins);    //Rebin following the current binning scheme
            
            std::vector<float> ratios=get_Ratio(temp,SMbkg);      //Get ratios needed for the pie charts
            for (int i=0;i<int(ratios.size());i++){
               bin_ratios[i][j]=ratios[i];
               std::cout<<ratios[i]<<std::endl;
            }
            j++;
         }
         
         //Plot pie charts with and without ttbar
         TCanvas can("","",300,300);
         TCanvas can_noTTbar;
         can.Divide(met_bins.size()-1,phi_bins.size()-1);
         can_noTTbar.Divide(met_bins.size()-1,phi_bins.size()-1);
         Int_t colors[] = {2,3,4,5,6,7,8};
         Int_t fix[] ={-3,3,1,-1};
         for (int i=0;i<nTotBins;i++){
            Float_t vals[7];
            TString binName="Bin"+std::to_string(i+1);
            // ~TPie *pie1 = new TPie("pie1",binName,7,vals,colors);
            TPie *pie1 = new TPie("pie1",binName,7,bin_ratios[i],colors);
            
            pie1->GetSlice(0)->SetTitle("SingleTop");
            pie1->GetSlice(1)->SetTitle("DrellYan");
            // ~pie1->GetSlice(2)->SetTitle("ttZ");
            // ~pie1->GetSlice(3)->SetTitle("ttH");
            // ~pie1->GetSlice(4)->SetTitle("W+Jets");
            // ~pie1->GetSlice(5)->SetTitle("Diboson");
            pie1->GetSlice(2)->SetTitle("");
            pie1->GetSlice(3)->SetTitle("");
            pie1->GetSlice(4)->SetTitle("");
            pie1->GetSlice(5)->SetTitle("");
            pie1->GetSlice(6)->SetTitle("TTbar");
            pie1->SetLabelsOffset(-.2);
            pie1->SetLabelFormat("%txt(%perc)");
            pie1->SetPercentFormat("%2.1f");
            pie1->SetRadius(.3);
            
            can.cd(nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)]);
            
            pie1->Draw("nol R <");
            if ((nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)])==1) pie1->MakeLegend();
            
            // ~can_noTTbar.cd(nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)]);
            // ~TPie *pie2=(TPie*) pie1->Clone("pie2");
            // ~pie2->GetSlice(6)->SetTitle("");
            // ~pie2->GetSlice(6)->SetValue(0);
            // ~pie2->Draw("nol R <");
            // ~if ((nTotBins-i+fix[(nTotBins-i)%(met_bins.size()-1)])==1) pie2->MakeLegend();
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
