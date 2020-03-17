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
   TString sample="MadGraph";
   // ~TString sample="dilepton";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="dilepton";
   
   // include signal to pseudo data
   // ~bool withBSM = true;
   bool withBSM = false;
   
   //Use scale factor
   bool withScaleFactor = false;
   // ~bool withScaleFactor = true;
   
   // perform toys studies?
   // ~bool toy_studies=false;
   bool toy_studies=true;
   
   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //==============================================
   // step 1 : open output file
   TString save_path = "TUnfold_results_"+sample+"_"+sample_response;
   if (withBSM) save_path+="_BSM";
   io::RootFileSaver saver(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);
   // ~io::RootFileSaver saver(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("TUnfold_SF91_%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
   if (withBSM) input_loc+="_BSM";

   TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
   TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

   if((!detectorBinning)||(!generatorBinning)) {
      cout<<"problem to read binning schemes\n";
   }

   // read histograms
   TH1 *histDataReco=histReader.read<TH1>(input_loc+"/histDataReco");
   TH1 *histDataTruth=histReader.read<TH1>(input_loc+"/histDataTruth");
   TH2 *histMCGenRec=histReader.read<TH2>(input_loc+"/histMCGenRec");

   if((!histDataReco)||(!histDataTruth)||(!histMCGenRec)) {
      cout<<"problem to read input histograms\n";
   }

   //========================
   // Step 3: unfolding

   // preserve the area
   // ~TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
   TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;

   // basic choice of regularisation scheme:
   //    curvature (second derivative)
   TUnfold::ERegMode regMode = TUnfold::kRegModeSize;

   // density flags
   // ~TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;
   TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeNone;

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

   // run the unfolding
   //
   unfold.DoUnfold(0.);//Without regularization
   
   // run the unfolding for toys
   TProfile *prof_pull_noRegularisation=0;
   TH1 *hist_coverage_noRegularisation=0;
   TH1 *hist_pull_noRegularisation=0;
   int MAXTOY=300;
   
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult");
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
   TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult",";bin",0,0,false);
   TH1 *hist_folded=unfold.GetFoldedOutput("hist_foldedResult",";bin",0,0,false);
   TH2 *cov_input=unfold.GetEmatrixInput("cov_input",";bin",0,0,true);
   TH1 *hist_unfolded_firstToy=0;
   
   if(toy_studies){
      for(int itoy=0;itoy<=MAXTOY;itoy++) {
         std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
         TH1 *hist_toybase_noRegularisation=unfold.GetFoldedOutput("",";bin",0,0,false,true);
         // ~unfold.SetInput(generatePoissonToy(hist_folded,itoy),1.0);
         unfold.SetInput(generatePoissonToy(histDataReco,itoy),1.0);
         unfold.DoUnfold(0.);//Without regularization
         if (itoy==0)hist_unfolded_firstToy=unfold.GetOutput("");
         analyzeToy(unfold.GetOutput(""),
                       hist_unfolded,
                       prof_pull_noRegularisation,
                       hist_coverage_noRegularisation,hist_pull_noRegularisation,MAXTOY);
      }
      saver.save(*hist_unfolded_firstToy,"hist_unfolded_firstToy");
      saver.save(*prof_pull_noRegularisation,"prof_pull_noRegularisation");
      saver.save(*hist_coverage_noRegularisation,"hist_coverage_noRegularisation");
      saver.save(*hist_pull_noRegularisation,"hist_pull_noRegularisation");
   }

   //===========================
   // Step 4: retreive and plot unfolding results
   // ~TH1 *hist_unfolded=unfold.GetOutput("hist_unfoldedResult","P_{T}^{#nu#nu} [GeV]","signal");
   saver.save(*hist_unfolded,"hist_unfoldedResult");
   saver.save(*hist_folded,"hist_foldedResult");
   saver.save(*cov_input,"cov_input");

}
