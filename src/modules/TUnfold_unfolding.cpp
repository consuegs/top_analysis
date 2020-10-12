//Script to perform the unfolding using the output of TUnfold_binning.cpp based on TUnfolding classes (Following TUnfold examples)
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
#include <TDecompSVD.h>
#include <TMatrixDEigen.h>

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "RooUnfold/RooUnfoldResponse.h"
#include "RooUnfold/RooUnfoldBinByBin.h"

using namespace std;


std::pair<float,int> getChi2NDF_withCorr(TH1 const *hist_res, TH1 const *hist_true, TH2 const *corr_res) {
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

void analyzeToy(TH1 const *hist_toy,
                       TH1 const *hist_truth, TH2 const *cov,
                       TProfile *&prof_pull,TH1 *&hist_coverage,TH1 *&hist_pull,TProfile *&prof_res,TH1 *&hist_res, TH1 *&hist_chi, int nToyTotal) {
   if(!prof_pull) {
      TString namePull(hist_truth->GetName()+TString("_toyPull"));
      TString nameRes(hist_truth->GetName()+TString("_toyRes"));
      TString nameCoverage(hist_truth->GetName()+TString("_toyCoverage"));
      TString namePullHist(hist_truth->GetName()+TString("_toyPullHist"));
      TString nameResHist(hist_truth->GetName()+TString("_toyResHist"));
      TString nameChiHist(hist_truth->GetName()+TString("_toyChiHist"));
      TString title=hist_truth->GetTitle();
      TArrayD const *xBins=hist_toy->GetXaxis()->GetXbins();
      if(xBins && (xBins->GetSize()>1)) {
         prof_pull=new TProfile(namePull,title,xBins->GetSize()-1,xBins->GetArray());
         prof_res=new TProfile(nameRes,title,xBins->GetSize()-1,xBins->GetArray());
         hist_coverage=new TH1D(nameCoverage,title,xBins->GetSize()-1,xBins->GetArray());
         hist_pull=new TH1D(namePullHist,";Pull;Entries",50,-3,3);
         hist_res=new TH1D(nameResHist,";Residual;Entries",50,-0.1,0.1);
         hist_chi=new TH1D(nameResHist,";#chi^{2}/ndf;Entries",50,0,3);
      } else {
         int nBins=hist_toy->GetNbinsX();
         double x0=hist_toy->GetXaxis()->GetXmin();
         double x1=hist_toy->GetXaxis()->GetXmax();
         prof_pull=new TProfile(namePull,title,nBins,x0,x1);
         prof_res=new TProfile(namePull,title,nBins,x0,x1);
         hist_coverage=new TH1D(nameCoverage,title,nBins,x0,x1);
         hist_pull=new TH1D(namePullHist,";Pull;Entries",50,-3,3);
         hist_res=new TH1D(nameResHist,";Residual;Entries",50,-0.1,0.1);
         hist_chi=new TH1D(nameResHist,";#chi^{2}/ndf;Entries",50,0,3);
      }
   }
   
   for(int i=1;i<hist_toy->GetNbinsX()+1;i++) {
      double pullI=(hist_truth->GetBinContent(i)-hist_toy->GetBinContent(i))/hist_truth->GetBinError(i);
      double resI=(hist_truth->GetBinContent(i)-hist_toy->GetBinContent(i))/hist_truth->GetBinContent(i);
      prof_pull->Fill(hist_toy->GetBinCenter(i),pullI);
      hist_pull->Fill(pullI);
      prof_res->Fill(hist_toy->GetBinCenter(i),resI);
      hist_res->Fill(resI);
      
      // ~std::cout<<hist_truth->GetBinContent(i)<<"   "<<hist_toy->GetBinContent(i)<<"   "<<hist_toy->GetBinCenter(i)<<std::endl;
      
      if(TMath::Abs(pullI)<1.) {
         hist_coverage->Fill(hist_toy->GetBinCenter(i),1./nToyTotal);
      }
   }
   
   auto chi2_pair=getChi2NDF_withCorr(hist_toy,hist_truth,cov);
   hist_chi->Fill(chi2_pair.first/(1.0*chi2_pair.second));
   
}

TH1 *generatePoissonToy(TH1 *base,int ntoy) {
   static TRandom *rnd=0;
   if(!rnd) rnd=new TRandom3();
   TH1 *r=(TH1 *)base->Clone(base->GetName()+TString::Format("_toy%d",ntoy));
   for(int ibin=0;ibin<=r->GetNbinsX()+1;ibin++) {
      double mu=r->GetBinContent(ibin);
      double c=0.;
      if(mu>0.) {
         c=rnd->Poisson(mu);
      }
      r->SetBinContent(ibin,c);
      r->SetBinError(ibin,TMath::Sqrt(c));
   }
   return r;
}

Config const &cfg=Config::get();
 
extern "C"
void run()
{
   // unfolded sample
   // ~TString sample="MadGraph";
   TString sample="dilepton";
   // ~TString sample="";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="dilepton";
   // ~TString sample_response="";
   
   // Use pT reweighted
   bool withPTreweight = false;
   // ~bool withPTreweight = true;
   TString scale="0.001";
   
   // Use Deep instead of pfMET
   bool withDeep = false;
   // ~bool withDeep = true;
   
   // Use puppi instead of pfMET
   // ~bool withPuppi = false;
   bool withPuppi = true;
   
   // Use same bin numbers for gen/true
   bool withSameBins = false;
   // ~bool withSameBins = true;
   
   // include signal to pseudo data
   // ~bool withBSM = true;
   bool withBSM = false;
   
   //Use scale factor
   bool withScaleFactor = false;
   // ~bool withScaleFactor = true;
   
   // perform toys studies?
   bool toy_studies=false;
   // ~bool toy_studies=true;
   int MAXTOY=3000;
   
   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //==============================================
   // step 1 : open output file
   TString save_path = "TUnfold_results_"+sample+"_"+sample_response;
   if (withBSM) save_path+="_BSM";
   if (withPuppi) save_path+="_Puppi";
   if (withDeep) save_path+="_Deep";
   if (withSameBins) save_path+="_SameBins";
   if (withPTreweight) save_path+="_PTreweight"+scale;
   io::RootFileSaver saver(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);
   // ~io::RootFileSaver saver(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
   if (withBSM) input_loc+="_BSM";
   if (withPuppi) input_loc+="_Puppi";
   if (withDeep) input_loc+="_Deep";
   if (withSameBins) input_loc+="_SameBins";
   if (withPTreweight) input_loc+="_PTreweight"+scale;

   TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
   TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

   if((!detectorBinning)||(!generatorBinning)) {
      cout<<"problem to read binning schemes\n";
   }

   // read histograms
   TH1 *histDataReco=histReader.read<TH1>(input_loc+"/histDataReco");
   TH1 *histDataReco_unscaled=(TH1*)histDataReco->Clone();
   TH1 *histDataReco_coarse=histReader.read<TH1>(input_loc+"/histDataReco_coarse");
   TH1 *histDataReco_coarse_unscaled=(TH1*)histDataReco_coarse->Clone();
   TH1 *histDataTruth=histReader.read<TH1>(input_loc+"/histDataTruth");
   TH1 *histSignalFraction=histReader.read<TH1>(input_loc+"/hist_SignalFraction");
   TH1 *histSignalFraction_coarse=histReader.read<TH1>(input_loc+"/hist_SignalFraction_coarse");
   // ~TH1 *histSignalFraction=histReader.read<TH1>("TUnfold_binning_MadGraph_MadGraph/hist_SignalFraction");
   // ~TH1 *histSignalFraction=histReader.read<TH1>("TUnfold_binning_dilepton_/hist_SignalFraction");
   TH2 *histMCGenRec=histReader.read<TH2>(input_loc+"/histMCGenRec");
   TH2 *histMCGenRec_sameBins=histReader.read<TH2>(input_loc+"/histMCGenRec_sameBins");

   if((!histDataReco)||(!histDataTruth)||(!histMCGenRec)) {
      cout<<"problem to read input histograms\n";
   }
   
   // remove fakes in reco by applying signal fraction
   histDataReco->Multiply(histSignalFraction);
   histDataReco_coarse->Multiply(histSignalFraction_coarse);
   
   // ~// test with flat reco distribution
   // ~for(int i=1; i<=histDataReco->GetNbinsX(); i++) {
      // ~histDataReco->SetBinContent(i,10000);
      // ~histDataReco->SetBinError(i,100);
   // ~}
   
   //========================
   // Step 3: unfolding
   
   
   for (bool regularisation : {false,true}) {
      
      // define saving folder
      TString saveFolder = (regularisation)? "reg/" : "";
      
      // preserve the area
      // ~TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
      TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;

      // basic choice of regularisation scheme:
      // ~TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
      TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
      // ~TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;

      // density flags
      TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
      // ~TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeNone;

      // detailed steering for regularisation
      const char *REGULARISATION_DISTRIBUTION=0;
      const char *REGULARISATION_AXISSTEERING="*[B]";

      // set up matrix of migrations
      TUnfoldDensity unfold(histMCGenRec,TUnfold::kHistMapOutputHoriz,
                           regMode,constraintMode,densityFlags,
                           generatorBinning,detectorBinning);

      // define the input vector (the measured data distribution)
      unfold.SetInput(histDataReco,1.0);
      // ~unfold.SetInput(histDataReco,0.);

      
      // ~unfold.ScanTau(1000,0,0,0);
      // ~unfold.ScanLcurve(1000,0,0,0);
      // ~std::cout<<unfold.GetTau()<<std::endl;
      // run the unfolding
      float tau=0;
      if (regularisation) {
         unfold.ScanTau(1000,0.0001,0.1,0); //With regularization
         tau=unfold.GetTau();
      }
      else unfold.DoUnfold(0.);//Without regularization
      
      // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult");
      // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
      TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
      TH1 *hist_folded=unfold.GetFoldedOutput("hist_foldedResult",";bin",0,0,false);
      TH2 *cov_input=unfold.GetEmatrixInput("cov_input",";bin",0,0,true);
      TH2 *cov_statResponse=unfold.GetEmatrixSysUncorr("cov_statResponse",";bin",0,0,true);
      TH2 *cov_Output=unfold.GetEmatrixTotal("cov_total",";bin",0,0,true);
      TH2 *corr_matrix=unfold.GetRhoIJtotal("Rho2D");
      
      // run the unfolding for toys
      if(toy_studies){
         //Define histograms for toy studies
         TProfile *prof_pull=0;
         TProfile *prof_res=0;
         TH1 *hist_coverage=0;
         TH1 *hist_pull=0;
         TH1 *hist_res=0;
         TH1 *hist_chi=0;
         TH1 *hist_unfolded_firstToy=0;
         for(int itoy=0;itoy<MAXTOY;itoy++) {
            std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
            // ~unfold.SetInput(generatePoissonToy(hist_folded,itoy),1.0);
            TH1 *reco_toy=generatePoissonToy(histDataReco_unscaled,itoy);
            reco_toy->Multiply(histSignalFraction);
            unfold.SetInput(reco_toy,1.0);
            if (regularisation) unfold.DoUnfold(tau); //With regularization
            else unfold.DoUnfold(0.);//Without regularization
            if (itoy==0)hist_unfolded_firstToy=unfold.GetOutput("",";bin",0,0,false);
            analyzeToy(unfold.GetOutput("","",0,0,false),
                          hist_unfolded,
                          // ~histDataTruth,
                          // ~cov_Output,
                          cov_input,
                          prof_pull,
                          hist_coverage,hist_pull,prof_res,hist_res,hist_chi,MAXTOY);
         }
         saver.save(*hist_unfolded_firstToy,saveFolder+"hist_unfolded_firstToy");
         saver.save(*prof_pull,saveFolder+"prof_pull");
         saver.save(*prof_res,saveFolder+"prof_res");
         saver.save(*hist_coverage,saveFolder+"hist_coverage");
         saver.save(*hist_pull,saveFolder+"hist_pull");
         saver.save(*hist_res,saveFolder+"hist_res");
         saver.save(*hist_chi,saveFolder+"hist_chi");
      }
      
      //Set stat. error to combined stat. error from input and response matrix
      for (int i=1; i<hist_unfolded->GetNbinsX();i++){
         hist_unfolded->SetBinError(i,sqrt(cov_input->GetBinContent(i,i)+cov_statResponse->GetBinContent(i,i)));
      }
      //===========================
      // Step 4: retreive and plot unfolding results
      // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult","P_{T}^{#nu#nu} [GeV]","signal");
      saver.save(*hist_unfolded,saveFolder+"hist_unfoldedResult");
      saver.save(*hist_folded,saveFolder+"hist_foldedResult");
      saver.save(*cov_input,saveFolder+"cov_input");
      saver.save(*cov_statResponse,saveFolder+"cov_statResponse");
      saver.save(*corr_matrix,saveFolder+"corr_matrix");
      saver.save(*cov_Output,saveFolder+"cov_output");
   }//end unfolding loop
   
   /*
   //SVD decomposition
   // ~TMatrixD test(histMCGenRec->GetNbinsY(),histMCGenRec->GetNbinsX());
   // ~for (int i=1; i<=histMCGenRec->GetNbinsX();i++){
      // ~for (int j=1; j<=histMCGenRec->GetNbinsY();j++){
         // ~test[j-1][i-1]=histMCGenRec->GetBinContent(i,j);
      // ~}
   // ~}
   
   TH2 *normResponse=unfold.GetProbabilityMatrix("prob_matrix",0,false);
   TMatrixD test(normResponse->GetNbinsY(),normResponse->GetNbinsX());
   for (int i=1; i<=normResponse->GetNbinsX();i++){
      for (int j=1; j<=normResponse->GetNbinsY();j++){
         test[j-1][i-1]=normResponse->GetBinContent(i,j);
      }
   }
   // ~test.Print();
   TDecompSVD testSVD(test);
   testSVD.Decompose();
   // ~testSVD.Print();
   std::cout<<testSVD.Condition()<<std::endl;
   
   // ~TMatrixDEigen testEigen(test);
   // ~TMatrixD eigenvalues=testEigen.GetEigenValues();
   // ~eigenvalues.Print();
   */
   
   //RooUnfold for Bin by Bin
   TH2F* histMCGenRec_sameBins_inv = (TH2F*)histMCGenRec_sameBins->Clone();
   for (int i=0; i<=histMCGenRec_sameBins->GetNbinsX()+1; i++){
      for (int j=0; j<=histMCGenRec_sameBins->GetNbinsY()+1; j++){
         histMCGenRec_sameBins_inv->SetBinContent(j,i,histMCGenRec_sameBins->GetBinContent(i,j));
         histMCGenRec_sameBins_inv->SetBinError(j,i,histMCGenRec_sameBins->GetBinError(i,j));
      }
   }
   RooUnfoldResponse responseRooUnfold(histMCGenRec_sameBins_inv->ProjectionX(), histMCGenRec_sameBins_inv->ProjectionY(), histMCGenRec_sameBins_inv);
   RooUnfoldBinByBin bbb(&responseRooUnfold, histDataReco_coarse);
   TH1* histBBB=bbb.Hreco(RooUnfold::kCovariance);
   TMatrixD covBBB_matrix=bbb.Ereco();
   TH2D covBBB(covBBB_matrix);
       
   saver.save(*histBBB,"BBB/hist_unfoldedResult");
   saver.save(covBBB,"BBB/cov_output");
   
   if(toy_studies){
      //Define histograms for toy studies
      TProfile *prof_pull=0;
      TProfile *prof_res=0;
      TH1 *hist_coverage=0;
      TH1 *hist_pull=0;
      TH1 *hist_res=0;
      TH1 *hist_chi=0;
      TH1 *hist_unfolded_firstToy=0;
      for(int itoy=0;itoy<MAXTOY;itoy++) {
         std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
         TH1 *reco_toy=generatePoissonToy(histDataReco_coarse_unscaled,itoy);
         reco_toy->Multiply(histSignalFraction_coarse);
         RooUnfoldBinByBin bbb(&responseRooUnfold, reco_toy);
         if (itoy==0)hist_unfolded_firstToy=bbb.Hreco(RooUnfold::kCovariance);
         analyzeToy(bbb.Hreco(RooUnfold::kCovariance),
                       histBBB,
                       &covBBB,
                       prof_pull,
                       hist_coverage,hist_pull,prof_res,hist_res,hist_chi,MAXTOY);
      }
      saver.save(*hist_unfolded_firstToy,"BBB/hist_unfolded_firstToy");
      saver.save(*prof_pull,"BBB/prof_pull");
      saver.save(*prof_res,"BBB/prof_res");
      saver.save(*hist_coverage,"BBB/hist_coverage");
      saver.save(*hist_pull,"BBB/hist_pull");
      saver.save(*hist_res,"BBB/hist_res");
      saver.save(*hist_chi,"BBB/hist_chi");
   }
}
