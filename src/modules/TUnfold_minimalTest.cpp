//Script to perform a minimal Test using TUnfold
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


#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;

Config const &cfg=Config::get();

std::pair<float,int> getChi2NDF_withCorr(TH1* hist_res, TH1* hist_true, TH2* corr_res) {
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

extern "C"
void run()
{
   // read histograms
   TH1D *histDataReco=new TH1D("Reco","",4,0.5,4.5);
   TH1D *histDataTruth=new TH1D("Gen","",2,0.5,2.5);
   TH2D *histMCGenRec=new TH2D("Response","",2,0,2,4,0,4);
   
   // ~histDataReco->Fill(0.5,4.);
   // ~histDataReco->Fill(1.5,3.);
   histDataReco->Fill(0.5,3.);
   histDataReco->Fill(1.5,4.);
   histDataReco->Fill(2.5,2.);
   histDataReco->Fill(3.5,1.);
   
   histMCGenRec->Fill(0.,0.,3.);
   histMCGenRec->Fill(0.,1.,2.);
   histMCGenRec->Fill(0.,2.,1.);
   histMCGenRec->Fill(0.,3.,0.);
   histMCGenRec->Fill(1.,0.,1.);
   histMCGenRec->Fill(1.,1.,1.);
   histMCGenRec->Fill(1.,2.,1.);
   histMCGenRec->Fill(1.,3.,1.);
   
   histDataTruth->Fill(0.5,6.);
   histDataTruth->Fill(1.5,4.);

   TUnfoldDensity unfold(histMCGenRec,TUnfold::kHistMapOutputHoriz,TUnfold::kRegModeSize,TUnfold::kEConstraintNone);

   // define the input vector (the measured data distribution)
   unfold.SetInput(histDataReco,1.0);
   // ~unfold.SetInput(histDataReco);

   // run the unfolding
   //
   std::list<float> chi_list;
   for (float tau = 0.001; tau<0.1; tau+=0.001) {
      
      // ~unfold.DoUnfold(0.);//Without regularization
      unfold.DoUnfold(tau);//Without regularization
      
      // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult");
      TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,true);
      TH1 *hist_folded=unfold.GetFoldedOutput("hist_foldedResult",";bin",0,0,true);
      TH2 *cov_input=unfold.GetEmatrixInput("cov_input",";bin",0,0,true);
      TH2 *cov_total=unfold.GetEmatrixTotal("cov_total",";bin",0,0,true);
      TH2 *cov_sys=unfold.GetEmatrixSysUncorr("cov_sys",";bin",0,0,true);
      TH2 *prob_matrix=unfold.GetProbabilityMatrix("prob_matrix",0,true);
      
      TH2 *cov_inputInverse=new TH2D("cov_inputInverse","",4,0,4,4,0,4);
      unfold.GetInputInverseEmatrix(cov_inputInverse);
      
      io::RootFileSaver saver(TString::Format("TUnfold_minimalTest%.1f.root",cfg.processFraction*100),"test/");
      saver.save(*histDataReco,"reco");
      saver.save(*hist_unfolded,"unfolded");
      saver.save(*hist_folded,"folded");
      saver.save(*histDataTruth,"true");
      saver.save(*cov_input,"cov_input");
      saver.save(*cov_inputInverse,"cov_inputInverse");
      saver.save(*cov_total,"cov_total");
      saver.save(*cov_sys,"cov_sys");
      saver.save(*prob_matrix,"prob_matrix");
      
      auto chi_pair = getChi2NDF_withCorr(hist_unfolded,histDataTruth,cov_total);
      
      chi_list.push_back(chi_pair.first);
      // ~std::cout<<tau<<"   "<<chi_pair.first<<std::endl;
   }
   
   for (auto chi : chi_list) {
      std::cout<<chi<<std::endl;
   }

}
