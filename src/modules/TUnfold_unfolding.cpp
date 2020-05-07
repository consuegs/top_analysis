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

using namespace std;

void analyzeToy(TH1 const *hist_toy,
                       TH1 const *hist_truth,
                       TProfile *&prof_pull,TH1 *&hist_coverage,TH1 *&hist_pull,int nToyTotal) {
   if(!prof_pull) {
      TString namePull(hist_truth->GetName()+TString("_toyPull"));
      TString nameCoverage(hist_truth->GetName()+TString("_toyCoverage"));
      TString namePullHist(hist_truth->GetName()+TString("_toyPullHist"));
      TString title=hist_truth->GetTitle();
      TArrayD const *xBins=hist_toy->GetXaxis()->GetXbins();
      if(xBins && (xBins->GetSize()>1)) {
         prof_pull=new TProfile(namePull,title,xBins->GetSize()-1,xBins->GetArray());
         hist_coverage=new TH1D(nameCoverage,title,xBins->GetSize()-1,xBins->GetArray());
         hist_pull=new TH1D(namePullHist,title,20,-3,3);
      } else {
         int nBins=hist_toy->GetNbinsX();
         double x0=hist_toy->GetXaxis()->GetXmin();
         double x1=hist_toy->GetXaxis()->GetXmax();
         prof_pull=new TProfile(namePull,title,nBins,x0,x1);
         hist_coverage=new TH1D(nameCoverage,title,nBins,x0,x1);
         hist_pull=new TH1D(namePullHist,title,20,-3,3);
      }
   }
   
   for(int i=0;i<=hist_toy->GetNbinsX()+1;i++) {
      double pullI=(hist_toy->GetBinContent(i)-hist_truth->GetBinContent(i))/hist_truth->GetBinError(i);
      prof_pull->Fill(hist_toy->GetBinCenter(i),pullI);
      hist_pull->Fill(pullI);
      
      if(TMath::Abs(pullI)<1.) {
         hist_coverage->Fill(hist_toy->GetBinCenter(i),1./nToyTotal);
      }
   }
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
   // ~bool withPTreweight = false;
   bool withPTreweight = true;
   
   // Use puppi instead of pfMET
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
   
   // perform toys studies?
   bool toy_studies=false;
   // ~bool toy_studies=true;
   
   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //==============================================
   // step 1 : open output file
   TString save_path = "TUnfold_results_"+sample+"_"+sample_response;
   if (withBSM) save_path+="_BSM";
   if (withPuppi) save_path+="_Puppi";
   if (withSameBins) save_path+="_SameBins";
   if (withPTreweight) save_path+="_PTreweight";
   io::RootFileSaver saver(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);
   // ~io::RootFileSaver saver(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
   if (withBSM) input_loc+="_BSM";
   if (withPuppi) input_loc+="_Puppi";
   if (withSameBins) input_loc+="_SameBins";
   if (withPTreweight) input_loc+="_PTreweight";

   TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
   TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

   if((!detectorBinning)||(!generatorBinning)) {
      cout<<"problem to read binning schemes\n";
   }

   // read histograms
   TH1 *histDataReco=histReader.read<TH1>(input_loc+"/histDataReco");
   TH1 *histDataTruth=histReader.read<TH1>(input_loc+"/histDataTruth");
   TH1 *histSignalFraction=histReader.read<TH1>(input_loc+"/hist_SignalFraction");
   // ~TH1 *histSignalFraction=histReader.read<TH1>("TUnfold_binning_MadGraph_MadGraph/hist_SignalFraction");
   // ~TH1 *histSignalFraction=histReader.read<TH1>("TUnfold_binning_dilepton_/hist_SignalFraction");
   TH2 *histMCGenRec=histReader.read<TH2>(input_loc+"/histMCGenRec");

   if((!histDataReco)||(!histDataTruth)||(!histMCGenRec)) {
      cout<<"problem to read input histograms\n";
   }
   
   // remove fakes in reco by applying signal fraction
   histDataReco->Multiply(histSignalFraction);
   
   // ~// test with flat reco distribution
   // ~for(int i=1; i<=histDataReco->GetNbinsX(); i++) {
      // ~histDataReco->SetBinContent(i,10000);
      // ~histDataReco->SetBinError(i,100);
   // ~}
   
   //========================
   // Step 3: unfolding

   // preserve the area
   // ~TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
   TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;

   // basic choice of regularisation scheme:
   TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
   // ~TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
   // ~TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;

   // density flags
   TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
   // ~TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeNone;

   // detailed steering for regularisation
   const char *REGULARISATION_DISTRIBUTION=0;
   const char *REGULARISATION_AXISSTEERING="*[B]";

   // set up matrix of migrations
   // ~TUnfoldDensity unfold(histMCGenRec,TUnfold::kHistMapOutputHoriz,
                        // ~regMode,constraintMode,densityFlags,
                        // ~generatorBinning,detectorBinning,
			// ~REGULARISATION_DISTRIBUTION,
			// ~REGULARISATION_AXISSTEERING);
   TUnfoldDensity unfold(histMCGenRec,TUnfold::kHistMapOutputHoriz,
                        regMode,constraintMode,densityFlags,
                        generatorBinning,detectorBinning);

   // define the input vector (the measured data distribution)
   unfold.SetInput(histDataReco,1.0);
   // ~unfold.SetInput(histDataReco,0.);

   
   // ~unfold.ScanTau(1000,0,0,0);
   // ~std::cout<<unfold.GetTau()<<std::endl;
   // run the unfolding
   //
   unfold.DoUnfold(0.);//Without regularization
   // ~unfold.DoUnfold(0.0006);//Without regularization
   
   // run the unfolding for toys
   TProfile *prof_pull_noRegularisation=0;
   TH1 *hist_coverage_noRegularisation=0;
   TH1 *hist_pull_noRegularisation=0;
   int MAXTOY=1000;
   
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult");
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
   TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
   TH1 *hist_folded=unfold.GetFoldedOutput("hist_foldedResult",";bin",0,0,false);
   TH2 *cov_input=unfold.GetEmatrixInput("cov_input",";bin",0,0,true);
   TH2 *cov_statResponse=unfold.GetEmatrixSysUncorr("cov_statResponse",";bin",0,0,true);
   TH2 *cov_Output=unfold.GetEmatrixTotal("cov_total",";bin",0,0,true);
   TH2 *corr_matrix=unfold.GetRhoIJtotal("Rho2D");
   TH1 *hist_unfolded_firstToy=0;
   
   std::cout<<unfold.GetRhoAvg()<<std::endl;
   
   if(toy_studies){
      for(int itoy=0;itoy<=MAXTOY;itoy++) {
         std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
         TH1 *hist_toybase_noRegularisation=unfold.GetFoldedOutput("",";bin",0,0,false,true);
         // ~unfold.SetInput(generatePoissonToy(hist_folded,itoy),1.0);
         unfold.SetInput(generatePoissonToy(histDataReco,itoy),1.0);
         unfold.DoUnfold(0.);//Without regularization
         if (itoy==0)hist_unfolded_firstToy=unfold.GetOutput("");
         analyzeToy(unfold.GetOutput("","",0,0,false),
                       hist_unfolded,
                       // ~histDataTruth,
                       prof_pull_noRegularisation,
                       hist_coverage_noRegularisation,hist_pull_noRegularisation,MAXTOY);
      }
      saver.save(*hist_unfolded_firstToy,"hist_unfolded_firstToy");
      saver.save(*prof_pull_noRegularisation,"prof_pull_noRegularisation");
      saver.save(*hist_coverage_noRegularisation,"hist_coverage_noRegularisation");
      saver.save(*hist_pull_noRegularisation,"hist_pull_noRegularisation");
   }
   
   //Set stat. error to combined stat. error from input and response matrix
   for (int i=1; i<hist_unfolded->GetNbinsX();i++){
      hist_unfolded->SetBinError(i,sqrt(cov_input->GetBinContent(i,i)+cov_statResponse->GetBinContent(i,i)));
   }
   //===========================
   // Step 4: retreive and plot unfolding results
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult","P_{T}^{#nu#nu} [GeV]","signal");
   saver.save(*hist_unfolded,"hist_unfoldedResult");
   saver.save(*hist_folded,"hist_foldedResult");
   saver.save(*cov_input,"cov_input");
   saver.save(*cov_statResponse,"cov_statResponse");
   saver.save(*corr_matrix,"corr_matrix");
   saver.save(*cov_Output,"cov_output");
   
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

}
