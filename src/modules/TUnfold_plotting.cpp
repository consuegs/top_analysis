//Script to plot the result of TUnfold_unfolding
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include "TUnfoldDensity.h"
#include <TRandom3.h>
#include <TProfile.h>
#include <TVectorD.h>


#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;

Config const &cfg=Config::get();

std::pair<float,int> getChi2NDF(TH1F* hist_res, TH1F* hist_true) {
   if (hist_res->GetNbinsX()!=hist_true->GetNbinsX()){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=1; i<=hist_res->GetNbinsX(); i++) {
         float diff = hist_res->GetBinContent(i)-hist_true->GetBinContent(i);
         float err = hist_res->GetBinError(i)*1.;
         chi2+= diff*diff/(err*err);
      }
      std::pair<float,int>result(chi2,hist_res->GetNbinsX());
      return result;
   }
}

std::pair<float,int> getChi2NDF_withCorr(TH1F* hist_res, TH1F* hist_true, TH2F* corr_res) {
   if (hist_res->GetNbinsX()!=hist_true->GetNbinsX()){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      
      TMatrixD diff(hist_res->GetNbinsX(),1);
      for (int i=1; i<=hist_res->GetNbinsX();i++){
         diff[i-1][0]=hist_res->GetBinContent(i)-hist_true->GetBinContent(i);
      }
      
      TMatrixD corr(hist_res->GetNbinsX(),hist_res->GetNbinsX());
      for (int i=1; i<=corr_res->GetNbinsX();i++){
         for (int j=1; j<=corr_res->GetNbinsY();j++){
            corr[i-1][j-1]=corr_res->GetBinContent(i,j);
         }
      }
      corr.Invert();
      
      TMatrixD resultMatrix(diff, TMatrixD::kTransposeMult,corr*diff);
      
      float chi2 = resultMatrix[0][0];
      std::pair<float,int>result(chi2,hist_res->GetNbinsX());
      return result;
   }
}

void plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver) {
   
   TH2F* tempHist;
   for (TString norm :{"","column"}){
      float sum=0;
      if (norm=="column"){ //Normalize each individual column of diagram
         tempHist=(TH2F*)responseHist->Clone();
         for (int x=1; x<=tempHist->GetNbinsX(); x++){
            sum=tempHist->Integral(x,x,1,tempHist->GetNbinsY());
            if (sum==0) continue;
            for (int y=1; y<=tempHist->GetNbinsY(); y++){
               if (tempHist->GetBinContent(x,y)!=0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
               else tempHist->SetBinContent(x,y,0.000002);
            }
         }
      }
      else { //Normalize each individual line of diagram
         tempHist=(TH2F*)responseHist->Clone();
         for (int y=1; y<=tempHist->GetNbinsY(); y++){
            sum=tempHist->Integral(1,tempHist->GetNbinsX(),y,y);
            if (sum==0) continue;
            for (int x=1; x<=tempHist->GetNbinsY(); x++){
               if (tempHist->GetBinContent(x,y)!=0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
               else tempHist->SetBinContent(x,y,0.000002);
            }
         }
      }
            
      TCanvas can;
      can.cd();
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.11);
      
      tempHist->GetYaxis()->SetTitleOffset(1.3);
      tempHist->GetXaxis()->SetTitleOffset(0.9);
      tempHist->GetZaxis()->SetTitleOffset(1.3);
      tempHist->GetYaxis()->SetTitleSize(0.05);
      tempHist->GetXaxis()->SetTitleSize(0.05);
      tempHist->GetZaxis()->SetTitleSize(0.05);
      tempHist->GetYaxis()->SetLabelSize(0.04);
      tempHist->GetXaxis()->SetLabelSize(0.04);
      tempHist->GetZaxis()->SetLabelSize(0.04);
      
      tempHist->SetMarkerColor(kRed);
      tempHist->SetMarkerSize(1.5);
      tempHist->SetTitle("");
      tempHist->GetZaxis()->SetTitle((norm=="column")?"column normalized distribution":"line normalized distribution");
      tempHist->SetMinimum(0.000001);
      tempHist->SetMaximum(1.);
      
      tempHist->SetStats(false);
      tempHist->Draw("colz text");
      
      tempHist->GetYaxis()->SetTitle("reco binNumber");
      tempHist->GetXaxis()->SetTitle("gen binNumber");
      tempHist->GetZaxis()->SetTitleOffset(0.55);
      tempHist->GetZaxis()->SetLabelOffset(0.0015);
      
      can.RedrawAxis();
      saver->save(can,"response"+norm+"/"+name,true,true);
   }
}

void plot_correlation(TH2F* corrMatrix, TString name, io::RootFileSaver* saver){
   TCanvas can2D;
   can2D.cd();
   
   gPad->SetRightMargin(0.2);
   corrMatrix->SetStats(0);
   corrMatrix->SetTitle("");
   corrMatrix->GetYaxis()->SetTitleOffset(0.6);
   corrMatrix->GetYaxis()->SetTitle("BinNo");
   corrMatrix->GetXaxis()->SetTitle("BinNo");
   corrMatrix->GetZaxis()->SetTitle("Correlation");
   corrMatrix->GetZaxis()->SetTitleOffset(1.2);
   corrMatrix->GetZaxis()->SetLabelOffset(0.01);
   corrMatrix->SetMarkerColor(kRed);
   corrMatrix->Draw("hcolz text");
   saver->save(can2D,name,true,true);
}
   

extern "C"
void run()
{
   // unfolded sample
   // ~TString sample="MadGraph";
   TString sample="TTbar_diLepton";
   // ~TString sample="";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="TTbar_diLepton";
   // ~TString sample_response="";
   
   // Use pT reweighted
   bool withPTreweight = false;
   // ~bool withPTreweight = true;
   TString scale="0.001";
   
   // ~// Use DNN instead of pfMET
   // ~bool withDNN = false;
   bool withDNN = true;
   
   // ~// Use puppi instead of pfMET
   bool withPuppi = false;
   // ~bool withPuppi = true;
   
   // Use same bin numbers for gen/true
   bool withSameBins = false;
   // ~bool withSameBins = true;
   
   // include signal to pseudo data
   // ~bool withBSM = true;
   bool withBSM = false;
   
   //Use scale factor
   bool withScaleFactor = false;
   // ~bool withScaleFactor = true;
   
   //Plot comparison
   // ~bool plotComparison = false;
   bool plotComparison = true;
   
   //Plot toy studies
   bool plotToyStudies = false;
   // ~bool plotToyStudies = true;
   
   //==============================================
   // step 1 : open output file
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format(!withScaleFactor ? "TUnfold_plotting%.1f" : "TUnfold_plotting_SF_%.1f",cfg.processFraction*100));

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
   TString input_loc_result="TUnfold_results_"+sample+"_"+sample_response;
   if (withBSM) {
      input_loc+="_BSM";
      input_loc_result+="_BSM";
   }
   if (withPuppi) {
      input_loc+="_Puppi";
      input_loc_result+="_Puppi";
   }
   if (withDNN) {
      input_loc+="_DNN";
      input_loc_result+="_DNN";
   }
   if (withSameBins) {
      input_loc+="_SameBins";
      input_loc_result+="_SameBins";
   }
   if (withPTreweight) {
      input_loc+="_PTreweight"+scale;
      input_loc_result+="_PTreweight"+scale;
   }

   TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
   TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

   if((!detectorBinning)||(!generatorBinning)) {
      cout<<"problem to read binning schemes\n";
   }

   // read histograms
   TH1F *realDis=histReader.read<TH1F>(input_loc+"/histDataTruth");
   TH1F *realDis_response=histReader.read<TH1F>(input_loc+"/histMCGenRec_projX");
   TH1F *unfolded=histReader.read<TH1F>(input_loc_result+"/hist_unfoldedResult");
   TH1F *unfolded_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_unfoldedResult");
   TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
   TH2F *corrMatrix=histReader.read<TH2F>(input_loc_result+"/corr_matrix");
   TH2F *corrMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/corr_matrix");
   TH2F *covMatrix=histReader.read<TH2F>(input_loc_result+"/cov_output");
   TH2F *covMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/cov_output");
   TH2F *covMatrix_bbb=histReader.read<TH2F>(input_loc_result+"/BBB/cov_output");
   TH2F *responseMatrix=histReader.read<TH2F>(input_loc+"/histMCGenRec_sameBins");

   if((!realDis)||(!realDis_response)||(!unfolded)) {
      cout<<"problem to read input histograms\n";
   }

   //========================
   // Step 3: plotting
   
   gfx::SplitCan can;
   can.can_.SetWindowSize(1800,600);
   can.pU_.SetLogy();
   can.pU_.cd();  //Use upper part of canvas
   
   
   //Initialize proper binning for plotting
   TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
   TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning(1)));
   
   int num_bins = binning_met.GetNoElements()*(binning_phi.GetNoElements()-1);
   int num_bins_met = binning_met.GetNoElements();
   
   binning_met.ResizeTo(binning_met.GetNoElements()+1);
   binning_met[binning_met.GetNoElements()-1] = 400;  //Plotting end for overflow bin
   
   Double_t xbins[num_bins+1];
   Double_t xbins_minus[num_bins+1];
   Double_t xbins_plus[num_bins+1];
   xbins[0] = 0;   //Plotting start
   xbins_minus[0] = -10;
   xbins_plus[0] = 10;
   
   
   int phi_bin = 0;
   for (int i=0; i<(num_bins); i++)   {
      xbins[i+1] = binning_met[i%num_bins_met+1]+phi_bin*400;
      xbins_minus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*400-10;
      xbins_plus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*400+10;
      if (i%num_bins_met==num_bins_met-1) phi_bin++;
   }
   
   unfolded->SetBins(num_bins,xbins);
   unfolded_reg->SetBins(num_bins,xbins_minus);
   unfolded_bbb->SetBins(num_bins,xbins_plus);
   realDis->SetBins(num_bins,xbins);
   realDis_response->SetBins(num_bins,xbins);
   for (int i=1; i<=num_bins; i++) {  //Set proper label for x axis
      int bin_label_no = (i-1)%num_bins_met+1;
      TString label;
      if (bin_label_no == num_bins_met) label = ">"+std::to_string((int)binning_met[bin_label_no-1]);
      else label = std::to_string((int)binning_met[bin_label_no-1])+"-"+std::to_string((int)binning_met[bin_label_no]);
      unfolded->GetXaxis()->SetBinLabel(i,label);
      realDis->GetXaxis()->SetBinLabel(i,label);
   }
   unfolded->GetXaxis()->SetTickLength(0.);
   unfolded->GetYaxis()->SetTickLength(0.008);
   unfolded->GetXaxis()->SetTitleOffset(1.5);
   unfolded->GetYaxis()->SetTitleOffset(0.8);
   unfolded->GetXaxis()->CenterLabels(false);
   
   unfolded->LabelsOption("v");
   realDis->LabelsOption("v");
   unfolded->SetMaximum(4*unfolded->GetMaximum());
   unfolded->SetMinimum(2);
   unfolded->SetLineColor(kBlack);
   unfolded->SetTitle(";p_{T}^{#nu#nu} (GeV);arbitrary unit");
   unfolded->SetStats(false);
   realDis->SetTitle(";p_{T}^{#nu#nu} (GeV);arbitrary unit");
   realDis->SetStats(false);
   
   
   if (plotComparison) {
      unfolded_reg->SetLineColor(kGreen+2);
      unfolded_bbb->SetLineColor(kViolet);
      unfolded_reg->SetMarkerColor(kGreen+2);
      unfolded_bbb->SetMarkerColor(kViolet);
      unfolded->Draw("pe1x0");  //Draw into current canvas
      unfolded_reg->Draw("pe1x0 same");
      unfolded_bbb->Draw("pe1x0 same");
   }
   else unfolded->Draw("pe1");  //Draw into current canvas
   
   realDis->SetLineColor(kRed-6);
   realDis->SetFillColor(kRed-6);
   realDis->Draw("hist same");   //Draw into current canvas (axis are not drawn again due to "same")
   
   realDis_response->SetLineColor(kBlue);
   realDis_response->SetFillColor(kBlue);
   realDis_response->Draw("hist same");
   
   //Draw vertical lines and binning ranges for deltaPhi
   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.03);
   aline->SetLineWidth(2);
   aline->DrawLine(800,2,800,unfolded->GetMaximum());
   aline->DrawLine(400,2,400,unfolded->GetMaximum());
   aline->DrawLine(800,2,800,unfolded->GetMaximum());
   atext->DrawLatex(75,0.5*unfolded->GetMaximum(),"0<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<0.7");
   atext->DrawLatex(475,0.5*unfolded->GetMaximum(),"0.7<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<1.4");
   atext->DrawLatex(875,0.5*unfolded->GetMaximum(),"1.4<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<3.14");
   
   //Get Chi2 and NDF
   auto Chi2Pair = getChi2NDF(unfolded,realDis);
   auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
   
   //Draw legend
   gfx::LegendEntries legE;
   if (plotComparison) {
      auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
      auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
      auto Chi2Pair_bbb = getChi2NDF(unfolded_bbb,realDis);
      auto Chi2Pair_corr_bbb = getChi2NDF_withCorr(unfolded_bbb,realDis,covMatrix_bbb);
      legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair.first,Chi2Pair.second,Chi2Pair_corr.first,Chi2Pair_corr.second),"pe");
      legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
      legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_bbb.first,Chi2Pair_bbb.second,Chi2Pair_corr_bbb.first,Chi2Pair_corr_bbb.second),"pe");
   }
   else {
      legE.append(*unfolded,"Unfolded","pe");
      atext->DrawLatex(30,10,TString::Format("#chi^{2}/NDF=%.1f/%i",Chi2Pair.first,Chi2Pair.second));
      atext->DrawLatex(30,3,TString::Format("#chi^{2}/NDF(corr.)=%.1f/%i",Chi2Pair_corr.first,Chi2Pair_corr.second));
   }
   legE.append(*realDis,"MC true ttbar","l");
   legE.append(*realDis_response,"MC true ttbar (response)","l");
   TLegend leg=legE.buildLegend(.10,.45,0.25,.65,1);
   leg.SetTextSize(0.03);
   leg.Draw();
   
   unfolded->Draw("axis same");
   
   if (withScaleFactor) {
      atext->DrawLatex(75,10,"With ScaleFactor");
   }
   
   //Change to lower part of canvas
   can.pL_.cd();
   can.pL_.SetBottomMargin(0.45);
   can.pL_.SetTickx(0);
   TH1F ratio;
   TH1F ratio_response;
   TH1F ratio_unfolded;
   TH1F ratio_unfolded_reg;
   TH1F ratio_unfolded_bbb;
   if (plotComparison) {
      ratio=hist::getRatio(*realDis,*realDis,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
      ratio_response=hist::getRatio(*realDis_response,*realDis,"ratio",hist::NOERR);
      ratio_unfolded=hist::getRatio(*unfolded,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_reg=hist::getRatio(*unfolded_reg,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_bbb=hist::getRatio(*unfolded_bbb,*realDis,"ratio",hist::ONLY1);
   }
   else {
      ratio=hist::getRatio(*realDis,*unfolded,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
      ratio_response=hist::getRatio(*realDis_response,*unfolded,"ratio",hist::NOERR);
      ratio_unfolded=hist::getRatio(*unfolded,*unfolded,"ratio",hist::ONLY1);
   }
      
   // ~ratio.SetMaximum(1.12);
   // ~ratio.SetMinimum(0.9);
   ratio.SetMaximum(1.2);
   ratio.SetMinimum(0.8);
   ratio.SetLineColor(kRed-6);
   ratio.SetMarkerColor(kRed-6);
   ratio.GetYaxis()->SetTitleOffset(0.3);
   ratio.GetXaxis()->SetTitleOffset(1.3);
   ratio.GetXaxis()->SetTickLength(0.);
   if (plotComparison) ratio.Draw("hist");
   else ratio.Draw();
   
   ratio_response.SetLineColor(kBlue);
   ratio_response.SetMarkerColor(kBlue);
   ratio_response.Draw("same");
   
   if (plotComparison) {
      ratio_unfolded.Draw("pe1x0 same");
      ratio_unfolded_reg.Draw("pe1x0 same");
      ratio_unfolded_bbb.Draw("pe1x0 same");
   }
   else {
      ratio_unfolded.SetMarkerSize(0);
      ratio_unfolded.Draw("same");
   }
   
   // ~aline->SetLineStyle(2);
   aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
   aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   
   //Print rel. uncertainties:
   for (int i=1; i<=ratio_unfolded.GetNbinsX(); i++){
      std::cout<<roundf(ratio_unfolded.GetBinError(i)*100 * 100) / 100<<std::endl;
   }
   std::cout<<"-----------------------"<<std::endl;
   for (int i=1; i<=ratio_unfolded_reg.GetNbinsX(); i++){
      std::cout<<roundf(ratio_unfolded_reg.GetBinError(i)*100 * 100) / 100<<std::endl;
   }
   
   //===========================
   // Step 4: save plot
   TString saveName=sample+"_"+sample_response;
   TString saveName2D="correlations/"+sample+"_"+sample_response;
   if (withBSM) {
      saveName+="_BSM";
      saveName2D+="_BSM";
   }
   if (withPuppi)  {
      saveName+="_Puppi";
      saveName2D+="_Puppi";
   }
   if (withDNN)  {
      saveName+="_DNN";
      saveName2D+="_DNN";
   }
   if (withSameBins) {
      saveName+="_SameBins";
      saveName2D+="_SameBins";
   }
   if (withPTreweight) {
      saveName+="_PTreweight"+scale;
      saveName2D+="_PTreweight"+scale;
   }
   if (plotComparison) saveName="compareMethods/"+saveName;
   saver.save(can,saveName);
   
   //Plot response matrix
   plot_response(responseMatrix,saveName,&saver);
   plot_correlation(corrMatrix,saveName2D,&saver);
   plot_correlation(corrMatrix_reg,saveName2D+"_reg",&saver);
   
   //Plot toy studies
   if (plotComparison) saveName=sample+"_"+sample_response;
   if (plotToyStudies) {
      TProfile *profPull=histReader.read<TProfile>(input_loc_result+"/prof_pull");
      TProfile *profPull_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_pull");
      TProfile *profPull_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_pull");
      TProfile *profRes=histReader.read<TProfile>(input_loc_result+"/prof_res");
      TProfile *profRes_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_res");
      TProfile *profRes_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_res");
      TH1F *histPull=histReader.read<TH1F>(input_loc_result+"/hist_pull");
      TH1F *histPull_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_pull");
      TH1F *histPull_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_pull");
      TH1F *histRes=histReader.read<TH1F>(input_loc_result+"/hist_res");
      TH1F *histRes_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_res");
      TH1F *histRes_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_res");
      TH1F *histChi=histReader.read<TH1F>(input_loc_result+"/hist_chi");
      TH1F *histChi_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_chi");
      TH1F *histChi_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_chi");
      for (TString type : {"Pull","Residual"}) {   //Plot bin-wise pull and residual
         TCanvas canToy;
         canToy.cd();
         
         //Set temp profiles/hists
         TProfile* temp=0;
         TProfile* temp_reg=0;
         TProfile* temp_bbb=0;
         if (type=="Pull") {
            temp=profPull;
            temp_reg=profPull_reg;
            temp_bbb=profPull_bbb;
         }
         else {
            temp=profRes;
            temp_reg=profRes_reg;
            temp_bbb=profRes_bbb;
         }
         
         //Get RMS from profile
         TH1D tempRMS=*(temp->ProjectionX());
         TH1D tempRMS_reg=*(TH1D*)tempRMS.Clone();
         TH1D tempRMS_bbb=*(TH1D*)tempRMS.Clone();
         for (int i=1; i<=profPull->GetNbinsX(); i++){
            tempRMS.SetBinContent(i,temp->GetBinError(i)*sqrt(temp->GetBinEntries(i)));
            tempRMS_reg.SetBinContent(i,temp_reg->GetBinError(i)*sqrt(temp_reg->GetBinEntries(i)));
            tempRMS_bbb.SetBinContent(i,temp_bbb->GetBinError(i)*sqrt(temp_bbb->GetBinEntries(i)));
         }
         
         temp->SetStats(0);
         temp->GetXaxis()->SetTitle("Bin");
         if (type=="Pull") temp->SetMaximum(1.2);
         else temp->SetMaximum(1.5);
         
         temp->SetMarkerColor(kBlue);
         temp_reg->SetMarkerColor(kRed);
         temp_bbb->SetMarkerColor(kBlack);
         tempRMS.SetMarkerColor(kBlue);
         tempRMS_reg.SetMarkerColor(kRed);
         tempRMS_bbb.SetMarkerColor(kBlack);
         
         temp->SetLineWidth(0);
         temp_reg->SetLineWidth(0);
         temp_bbb->SetLineWidth(0);
         tempRMS.SetLineWidth(0);
         tempRMS_reg.SetLineWidth(0);
         tempRMS_bbb.SetLineWidth(0);
         
         temp->SetMarkerSize(1.3);
         temp_reg->SetMarkerSize(1.3);
         temp_bbb->SetMarkerSize(1.3);
         tempRMS.SetMarkerSize(1.3);
         tempRMS_reg.SetMarkerSize(1.3);
         tempRMS_bbb.SetMarkerSize(1.3);
         
         tempRMS.SetMarkerStyle(4);
         tempRMS_reg.SetMarkerStyle(4);
         tempRMS_bbb.SetMarkerStyle(4);
         
         temp->Draw("");
         temp_reg->Draw("same");
         temp_bbb->Draw("same");
         tempRMS.Draw("same");
         tempRMS_reg.Draw("same");
         tempRMS_bbb.Draw("same");
         
         aline->DrawLine(0.5,0,num_bins+0.5,0);
         if (type=="Pull") aline->DrawLine(0.5,1,num_bins+0.5,1);
         
         gfx::LegendEntries legE_pull;
         legE_pull.append(*temp,"NoReg "+type+" Mean","p");
         legE_pull.append(tempRMS,"NoReg "+type+" RMS","p");
         legE_pull.append(*temp_reg,"Reg "+type+" Mean","p");
         legE_pull.append(tempRMS_reg,"Reg "+type+" RMS","p");
         legE_pull.append(*temp_bbb,"BBB "+type+" Mean","p");
         legE_pull.append(tempRMS_bbb,"BBB "+type+" RMS","p");
         TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.8,.65,2);
         leg_pull.SetTextSize(0.04);
         leg_pull.Draw();
         
         saver.save(canToy,"toyStudies/"+saveName+"/prof"+type);
      }
      
      for (TString type : {"Pull","Residual","Chi"}) {
         TCanvas canToy;
         canToy.cd();
         
         //Set temp profiles/hists
         TH1* temp=0;
         TH1* temp_reg=0;
         TH1* temp_bbb=0;
         if (type=="Pull") {
            temp=histPull;
            temp_reg=histPull_reg;
            temp_bbb=histPull_bbb;
         }
         else if (type=="Residual") {
            temp=histRes;
            temp_reg=histRes_reg;
            temp_bbb=histRes_bbb;
         }
         else {
            temp=histChi;
            temp_reg=histChi_reg;
            temp_bbb=histChi_bbb;
         }
         
         temp->SetLineColor(kBlue);
         temp_reg->SetLineColor(kRed);
         temp_bbb->SetLineColor(kBlack);
         temp->SetLineWidth(2);
         temp_reg->SetLineWidth(2);
         temp_bbb->SetLineWidth(2);
         
         temp_bbb->SetStats(0);
         temp_bbb->Draw("hist");
         temp->Draw("hist same");
         temp_reg->Draw("hist same");
         gfx::LegendEntries legE_pull;
         legE_pull.append(*temp,TString::Format("NoReg [#mu=%.3f #sigma=%.3f]",temp->GetMean(),temp->GetRMS()),"l");
         legE_pull.append(*temp_reg,TString::Format("Reg [#mu=%.3f #sigma=%.3f]",temp_reg->GetMean(),temp_reg->GetRMS()),"l");
         legE_pull.append(*temp_bbb,TString::Format("BBB [#mu=%.3f #sigma=%.3f]",temp_bbb->GetMean(),temp_bbb->GetRMS()),"l");
         TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.4,.65,1);
         leg_pull.SetTextSize(0.04);
         leg_pull.Draw();
         
         saver.save(canToy,"toyStudies/"+saveName+"/hist"+type);
      }
         
   }
      
}
