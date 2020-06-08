//Script to plot the result of TUnfold_unfolding
#include <iostream>
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
   
   //Normalize each individual line of diagram
   float sum;
   for (int y=1; y<=responseHist->GetNbinsY(); y++){
      sum=responseHist->Integral(1,responseHist->GetNbinsX(),y,y);
      if (sum==0) continue;
      for (int x=1; x<=responseHist->GetNbinsY(); x++){
         if (responseHist->GetBinContent(x,y)!=0)responseHist->SetBinContent(x,y,responseHist->GetBinContent(x,y)/sum);
         else responseHist->SetBinContent(x,y,0.000002);
      }
   }
         
   TCanvas can;
   can.cd();
   gPad->SetRightMargin(0.2);
   gPad->SetLeftMargin(0.13);
   gPad->SetBottomMargin(0.11);
   
   responseHist->GetYaxis()->SetTitleOffset(1.3);
   responseHist->GetXaxis()->SetTitleOffset(0.9);
   responseHist->GetZaxis()->SetTitleOffset(1.3);
   responseHist->GetYaxis()->SetTitleSize(0.05);
   responseHist->GetXaxis()->SetTitleSize(0.05);
   responseHist->GetZaxis()->SetTitleSize(0.05);
   responseHist->GetYaxis()->SetLabelSize(0.04);
   responseHist->GetXaxis()->SetLabelSize(0.04);
   responseHist->GetZaxis()->SetLabelSize(0.04);
   
   responseHist->SetMarkerColor(kRed);
   responseHist->SetMarkerSize(1.5);
   responseHist->SetTitle("");
   responseHist->GetZaxis()->SetTitle("line normalized distribution");
   responseHist->SetMinimum(0.000001);
   responseHist->SetMaximum(1.);
   
   responseHist->SetStats(false);
   responseHist->Draw("colz text");
   
   responseHist->GetYaxis()->SetTitle("reco binNumber");
   responseHist->GetXaxis()->SetTitle("gen binNumber");
   
   can.RedrawAxis();
   saver->save(can,"response/"+name,true,true);
}
   

extern "C"
void run()
{
   // unfolded sample
   TString sample="MadGraph";
   // ~TString sample="dilepton";
   // ~TString sample="";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="dilepton";
   // ~TString sample_response="";
   
   // Use pT reweighted
   bool withPTreweight = false;
   // ~bool withPTreweight = true;
   
   // ~// Use deep instead of pfMET
   bool withDeep = false;
   // ~bool withDeep = true;
   
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
   if (withDeep) {
      input_loc+="_Deep";
      input_loc_result+="_Deep";
   }
   if (withSameBins) {
      input_loc+="_SameBins";
      input_loc_result+="_SameBins";
   }
   if (withPTreweight) {
      input_loc+="_PTreweight";
      input_loc_result+="_PTreweight";
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
   TH2F *corrMatrix=histReader.read<TH2F>(input_loc_result+"/corr_matrix");
   TH2F *covMatrix=histReader.read<TH2F>(input_loc_result+"/cov_output");
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
   xbins[0] = 0;   //Plotting start
   
   int phi_bin = 0;
   for (int i=0; i<(num_bins); i++)   {
      xbins[i+1] = binning_met[i%num_bins_met+1]+phi_bin*400;
      if (i%num_bins_met==num_bins_met-1) phi_bin++;
   }
   
   unfolded->SetBins(num_bins,xbins);
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
   
   unfolded->Draw("pe1");  //Draw into current canvas 
   realDis->SetLineColor(kRed-6);
   realDis->SetFillColor(kRed-6);
   realDis->Draw("hist same");   //Dra into current canvas (axis are not drawn again due to "same")
   
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
   
   //Draw legend
   gfx::LegendEntries legE;
   legE.append(*unfolded,"Unfolded","pe");
   legE.append(*realDis,"MC true ttbar","l");
   legE.append(*realDis_response,"MC true ttbar (response)","l");
   TLegend leg=legE.buildLegend(.15,.45,0.3,.6,1);
   leg.SetTextSize(0.035);
   leg.Draw();
   
   unfolded->Draw("axis same");
   
   if (withScaleFactor) {
      atext->DrawLatex(75,10,"With ScaleFactor");
   }
   
   //Get Chi2 and NDF
   auto Chi2Pair = getChi2NDF(unfolded,realDis);
   auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
   atext->DrawLatex(75,10,TString::Format("#chi^{2}/NDF=%.1f/%i",Chi2Pair.first,Chi2Pair.second));
   atext->DrawLatex(75,3,TString::Format("#chi^{2}/NDF(corr.)=%.1f/%i",Chi2Pair_corr.first,Chi2Pair_corr.second));
   
   //Change to lower part of canvas
   can.pL_.cd();
   can.pL_.SetBottomMargin(0.45);
   can.pL_.SetTickx(0);
   TH1F ratio=hist::getRatio(*realDis,*unfolded,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
   // ~ratio.SetMaximum(1.12);
   // ~ratio.SetMinimum(0.9);
   ratio.SetMaximum(1.2);
   ratio.SetMinimum(0.8);
   ratio.SetLineColor(kRed-6);
   ratio.SetMarkerColor(kRed-6);
   ratio.GetYaxis()->SetTitleOffset(0.3);
   ratio.GetXaxis()->SetTitleOffset(1.3);
   ratio.GetXaxis()->SetTickLength(0.);
   ratio.Draw();
   
   TH1F ratio_response=hist::getRatio(*realDis_response,*unfolded,"ratio",hist::NOERR);
   ratio_response.SetLineColor(kBlue);
   ratio_response.SetMarkerColor(kBlue);
   ratio_response.Draw("same");
   
   TH1F uncertainty_unfolded=hist::getRatio(*unfolded,*unfolded,"ratio",hist::ONLY1);  //Get Ratio between unfolded and unfolded to draw unfolded uncertainty around 1
   uncertainty_unfolded.SetMarkerSize(0);
   uncertainty_unfolded.Draw("same");
   
   // ~aline->SetLineStyle(2);
   aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
   aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   
   //Plot correlation matrix
   TCanvas can2D;
   can2D.cd();
   
   gPad->SetRightMargin(0.2);
   corrMatrix->SetStats(0);
   corrMatrix->SetTitle("");
   corrMatrix->GetYaxis()->SetTitleOffset(0.6);
   corrMatrix->GetYaxis()->SetTitle("BinNo");
   corrMatrix->GetXaxis()->SetTitle("BinNo");
   corrMatrix->GetZaxis()->SetTitleOffset(1.2);
   corrMatrix->GetZaxis()->SetLabelOffset(0.01);
   corrMatrix->SetMarkerColor(kRed);
   corrMatrix->Draw("hcolz text");

   
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
   if (withDeep)  {
      saveName+="_Deep";
      saveName2D+="_Deep";
   }
   if (withSameBins) {
      saveName+="_SameBins";
      saveName2D+="_SameBins";
   }
   if (withPTreweight) {
      saveName+="_PTreweight";
      saveName2D+="_PTreweight";
   }
   saver.save(can,saveName);
   saver.save(can2D,saveName2D);
   
   //Plot response matrix
   plot_response(responseMatrix,saveName,&saver);

}
