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
#include <TParameter.h>

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
         hist_pull=new TH1D(namePullHist,";(x^{unf}_{i}-x_{i}^{toy})/#sigma_{i}^{unf};Number of toys",50,-3,3);
         hist_res=new TH1D(nameResHist,";x_{i}^{unf}-x_{i}^{toy};Number of toys",50,-0.1,0.1);
         hist_chi=new TH1D(nameChiHist,";#chi^{2};Number of toys",50,0,30);
      } else {
         int nBins=hist_toy->GetNbinsX();
         double x0=hist_toy->GetXaxis()->GetXmin();
         double x1=hist_toy->GetXaxis()->GetXmax();
         prof_pull=new TProfile(namePull,title,nBins,x0,x1);
         prof_res=new TProfile(nameRes,title,nBins,x0,x1);
         hist_coverage=new TH1D(nameCoverage,title,nBins,x0,x1);
         hist_pull=new TH1D(namePullHist,";(x^{unf}_{i}-x_{i}^{toy})/#sigma_{i}^{unf};Number of toys",50,-3,3);
         hist_res=new TH1D(nameResHist,";x_{i}^{unf}-x_{i}^{toy};Number of toys",50,-0.1,0.1);
         hist_chi=new TH1D(nameChiHist,";#chi^{2};Number of toys",50,0,30);
      }
   }
   
   for(int i=1;i<hist_toy->GetNbinsX()+1;i++) {
      double pullI=(hist_truth->GetBinContent(i)-hist_toy->GetBinContent(i))/hist_truth->GetBinError(i);
      double resI=(hist_truth->GetBinContent(i)-hist_toy->GetBinContent(i))/hist_truth->GetBinContent(i);
      prof_pull->Fill(hist_toy->GetBinCenter(i),pullI);
      hist_pull->Fill(pullI);
      prof_res->Fill(hist_toy->GetBinCenter(i),resI);
      hist_res->Fill(resI);
                  
      if(TMath::Abs(pullI)<1.) {
         hist_coverage->Fill(hist_toy->GetBinCenter(i),1./nToyTotal);
      }
   }
   
   auto chi2_pair=getChi2NDF_withCorr(hist_toy,hist_truth,cov);
   // ~hist_chi->Fill(chi2_pair.first/(1.0*chi2_pair.second));
   hist_chi->Fill(chi2_pair.first);
   
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
      if(c>0) r->SetBinError(ibin,TMath::Sqrt(c));
   }
   
   return r;
}

Config const &cfg=Config::get();
 
extern "C"
void run()
{
   
   gErrorIgnoreLevel = kWarning;    // Only print warnings. Otherwise scanTau method spams a lot!
   
   //////////////////////////////////
   // Set Systematic Uncertainties //
   //////////////////////////////////
   
   // Nominal Set
   
   // ~std::vector<TString> systVec = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // JES complete Set
   
   // ~std::vector<TString> systVec = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN","JESUserDefinedHEM1516_DOWN","LUMI_UP","LUMI_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MTOP169p5_UP","MTOP175p5_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Only JES complete
   // ~std::vector<TString> systVec = {"JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN"};
   
   
   // ~std::vector<TString> systVec = {"Nominal","JESTotal_UP","JESTotal_DOWN"};
   // ~std::vector<TString> systVec = {"Nominal","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systVec = {"Nominal","MATCH_UP","MATCH_DOWN","MTOP169p5","MTOP175p5"};
   // ~std::vector<TString> systVec = {"MEFACSCALE_UP"};
   // ~std::vector<TString> systVec = {"L1PREFIRING_UP","L1PREFIRING_DOWN"};
   // ~std::vector<TString> systVec = {"JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN"};
   // ~std::vector<TString> systVec = {"MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_ST_UP"};
   // ~std::vector<TString> systVec = {"TWDS"};
   // ~std::vector<TString> systVec = {"Nominal","TWDS"};
   // ~std::vector<TString> systVec = {"Nominal","XSEC_DY_UP"};
   // ~std::vector<TString> systVec = {"Nominal","MESCALE_UP"};
   // ~std::vector<TString> systVec = {"applyGenLevel_DeltaRcut"};
   // ~std::vector<TString> systVec = {"JETPILEUPID_DOWN"};
   // ~std::vector<TString> systVec = {"Nominal","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN"};
   std::vector<TString> systVec = {"Nominal","CR1","CR2","ERDON"};
   // ~std::vector<TString> systVec = {"Nominal"};
   
   //Remove HEM unc. for all year except 2018
   auto itr =std::find(systVec.begin(), systVec.end(), "JESUserDefinedHEM1516_DOWN");
   if (itr != systVec.end() && cfg.year_int != 3){
      systVec.erase(itr);
   }
   
   /////////////////////////
   // Add PDF Systematics //
   /////////////////////////
   
   // add pdf unc
   // ~for(int i=1; i<=50; i++){
      // ~systVec.push_back(TString::Format("PDF_%i_UP",i));
      // ~systVec.push_back(TString::Format("PDF_%i_DOWN",i));
   // ~}
   
   //////////////////////////
   // Define Distributions //
   //////////////////////////
   
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu","2D_dPhi_pTnunu_new","2D_dPhi_pTnunu_new40","2D_dPhi_pTnunu_DNN","2D_dPhi_pTnunu_new_DNN","2D_dPhi_pTnunu_new40_DNN","pTnunu","dPhi","pTnunu_DNN","pTnunu_new_DNN","dPhi_DNN","pTll","inclusive"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu","2D_dPhi_pTnunu_new","2D_dPhi_pTnunu_new40","pTnunu","dPhi","pTll","inclusive"};
   // ~std::vector<TString> distributions = {"pTll"};
   // ~std::vector<TString> distributions = {"inclusive"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_DNN"};
   // ~std::vector<TString> distributions = {"pTnunu_new"};
   // ~std::vector<TString> distributions = {"pTnunu_new_DNN"};
   // ~std::vector<TString> distributions = {"dPhi_DNN"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_30StabPur_DNN"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_mergedBins_DNN"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN"};
   // ~std::vector<TString> distributions = {"inclusive","dPhi_new_DNN","pTnunu_new_DNN","2D_dPhi_pTnunu_new_DNN","2D_dPhi_pTnunu_new_30StabPur_DNN","2D_dPhi_pTnunu_new_30StabPur12Bins_DNN","2D_dPhi_pTnunu_new_mergedBins_DNN","2D_dPhi_pTnunu_new_30StabPur2PhiBins16_DNN","2D_dPhi_pTnunu_new_30StabPur2PhiBins14_DNN"};
   // ~std::vector<TString> distributions = {"pTnunu_new_DNN","pTnunu_DNN","2D_dPhi_pTnunu_new_DNN","2D_dPhi_pTnunu_new40_DNN"};
   // ~std::vector<TString> distributions = {"pTnunu_new_DNN","2D_dPhi_pTnunu_new_30StabPur12Bins_DNN"};
   
   // ~std::vector<TString> distributions = {"inclusive","dPhi_new_DNN","pTnunu_new_DNN","2D_dPhi_pTnunu_new_30StabPur12Bins_DNN"};
   // ~std::vector<TString> distributions = {"pTnunu_new","dPhi","2D_dPhi_pTnunu_new_30StabPur12Bins"};
   // ~std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_30StabPur12Bins"};
   
   std::vector<TString> distributions = {"2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN","2D_dPhi_pTnunu_new_30StabPur9Bins_diffDet_DNN","2D_dPhi_pTnunu_new_30StabPur12Bins_DNN"};
   
   
   // ~bool verbose = false;
   bool verbose = true;
   
   for (auto syst : systVec){
      
      std::cout<<"Running "<<syst<<std::endl;
      
      // avoid spamming from TUnfold
      if (!verbose) {
         gErrorIgnoreLevel = kFatal;
         std::cout.setstate(std::ios_base::failbit);
      }
      
      Systematic::Systematic currentSystematic(syst);
      bool isNominal = (std::find(Systematic::nominalTypes.begin(), Systematic::nominalTypes.end(), currentSystematic.type()) != Systematic::nominalTypes.end());
      
      for (auto varName : distributions){
         
         std::cout<<"--------------------------------"<<currentSystematic.name()<<"-------------------------------------------------------"<<std::endl;
      
         TString sample = cfg.tunfold_InputSamples[0];
         TString sample_response = cfg.tunfold_ResponseSample;
         
         // Use pT reweighted
         bool withPTreweight = cfg.tunfold_withPTreweight;
         TString scale = cfg.tunfold_scalePTreweight;
         
         // Use phi reweighted
         bool withPHIreweight = cfg.tunfold_withPHIreweight;
         TString scale_phi = cfg.tunfold_scalePHIreweight;
         
         // Use DNN instead of pfMET
         bool withDNN = cfg.tunfold_withDNN;
         
         // Use pf instead of PuppiMET
         bool withPF = cfg.tunfold_withPF;
         
         // Use puppi instead of pfMET
         bool withPuppi = !withDNN && !withPF;
         
         // Use same bin numbers for gen/true
         bool withSameBins = cfg.tunfold_withSameBins;
         
         // include signal to pseudo data
         bool withBSM = cfg.tunfold_withBSM;
         
         //Use scale factor
         bool withScaleFactor = cfg.tunfold_withScaleFactor;
         
         //Use alternative pseudo data (amcAtNLO)
         // ~bool useAltReco = true;
         bool useAltReco = false;
         
         //Use real data
         // ~bool useRealData = false;
         bool useRealData = true;
         
         //Use Single Top DS
         // ~bool useSingleTopDS = false;
         bool useSingleTopDS = true;
         
         //Cannot use realData and altPseudo data simultaneously
         if (useRealData && useAltReco){
            std::cerr<<"Error: Cannot use realData and altPseudo data simultaneously"<<std::endl;
            exit(5);
         }
         
         // perform toys studies?
         bool toy_studies = cfg.tunfold_plotToyStudies;
         int MAXTOY=5000;
         // ~int MAXTOY=2;
         
         // switch on histogram errors
         TH1::SetDefaultSumw2();

         //==============================================
         // step 1 : read binning schemes and input histograms
         TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
         if (withBSM) input_loc+="_BSM";
         if (withPuppi) input_loc+="_Puppi";
         if (withPF) input_loc+="_PF";
         if (withDNN) input_loc+="_DNN";
         if (withSameBins) input_loc+="_SameBins";
         if (withPTreweight) input_loc+="_PTreweight"+scale;
         else if (withPHIreweight) input_loc+="_PHIreweight"+scale_phi;
         io::RootFileReader histReader(TString::Format(!cfg.tunfold_withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",input_loc.Data(),(isNominal)? currentSystematic.name().Data() : "Nominal",cfg.processFraction*100));
         io::RootFileReader histReader_syst(TString::Format(!cfg.tunfold_withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",input_loc.Data(),currentSystematic.name().Data(),cfg.processFraction*100));
         
         // ~TString input_loc_old = (useAltReco)? "TUnfold_binning_TTbar_amcatnlo_"+sample_response : input_loc;
         TString input_loc_old = input_loc;
         input_loc += "/"+varName;

         TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
         TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

         if((!detectorBinning)||(!generatorBinning)) {
            cout<<"problem to read binning schemes\n";
         }

         // read histograms
         TString recoName = (useAltReco)? "/histDataRecoAlt" : ((useRealData)? "/histDataReal": ((useSingleTopDS)? "/histDataRecoSTDS" : "/histDataReco"));
         TH1 *histDataReco=histReader.read<TH1>(input_loc+recoName);
         TH1 *histDataReco_unscaled=(TH1*)histDataReco->Clone();
         TH1 *histDataReco_coarse=histReader.read<TH1>(input_loc+recoName+"_coarse");
         TH1 *histDataReco_coarse_unscaled=(TH1*)histDataReco_coarse->Clone();
         TH1 *histDataTruth=histReader.read<TH1>(input_loc+"/histDataTruth");
         TH1 *histSignalFraction=histReader_syst.read<TH1>(input_loc+"/hist_SignalFraction");
         TH1 *histMCRec_bkg=histReader_syst.read<TH1>(input_loc+((useSingleTopDS)? "/histMCRecSTDS_bkg" : "/histMCRec_bkg"));
         TH1 *histMCRec_bkg_coarse=histReader_syst.read<TH1>(input_loc+((useSingleTopDS)? "/histMCRecSTDS_bkg_coarse" : "/histMCRec_bkg_coarse"));
         TH1 *histSignalFraction_coarse=histReader_syst.read<TH1>(input_loc+"/hist_SignalFraction_coarse");
         TH2 *histMCGenRec=histReader_syst.read<TH2>(input_loc+"/histMCGenRec");
         TH2 *histMCGenRec_sameBins=histReader_syst.read<TH2>(input_loc+"/histMCGenRec_sameBins");

         if((!histDataReco)||(!histDataTruth)||(!histMCGenRec)) {
            cout<<"problem to read input histograms\n";
         }
         
         if(!(withPTreweight || withPHIreweight)){    //No background and signal fraction considered in reweighting study
            // substract non-ttbar background
            histDataReco->Add(histMCRec_bkg,-1);
            histDataReco_coarse->Add(histMCRec_bkg_coarse,-1);
            
            // remove fakes in reco by applying signal fraction
            histDataReco->Multiply(histSignalFraction);
            histDataReco_coarse->Multiply(histSignalFraction_coarse);
         }
         
         // ~// test with flat reco distribution
         // ~for(int i=1; i<=histDataReco->GetNbinsX(); i++) {
            // ~histDataReco->SetBinContent(i,10000);
            // ~histDataReco->SetBinError(i,100);
         // ~}
      
         //==============================================
         // step 2 : open output file
         TString save_path = (useAltReco)? "TUnfold_results_TTbar_amcatnlo_"+sample_response : "TUnfold_results_"+sample+"_"+sample_response;
         if (useRealData) save_path = "TUnfold_results_realData_"+sample_response;
         if (useSingleTopDS) save_path+="_SingleTopDS";
         if (withBSM) save_path+="_BSM";
         if (withPuppi) save_path+="_Puppi";
         if (withPF) save_path+="_PF";
         if (withDNN) save_path+="_DNN";
         if (withSameBins) save_path+="_SameBins";
         if (withPTreweight) save_path+="_PTreweight"+scale;
         else if (withPHIreweight) save_path+="_PHIreweight"+scale_phi;
         save_path += "/"+varName;
         io::RootFileSaver saver(TString::Format(!cfg.tunfold_withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",input_loc_old.Data(),currentSystematic.name().Data(),cfg.processFraction*100),save_path);
         
         //========================
         // Step 3: unfolding
         
         if (varName != "inclusive") {    // use only BBB for inclusive xsec
            for (bool regularisation : {false,true}) {
               
               // define saving folder
               TString saveFolder = (regularisation)? "reg/" : "";
               
               // preserve the area
               TUnfold::EConstraint constraintMode= TUnfold::kEConstraintArea;
               // ~TUnfold::EConstraint constraintMode= TUnfold::kEConstraintNone;

               // basic choice of regularisation scheme:
               // ~TUnfold::ERegMode regMode = TUnfold::kRegModeSize;
               TUnfold::ERegMode regMode = TUnfold::kRegModeCurvature;
               // ~TUnfold::ERegMode regMode = TUnfold::kRegModeDerivative;

               // density flags
               TUnfoldDensity::EDensityMode densityFlags=TUnfoldDensity::kDensityModeBinWidth;

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
               TH2 *cov_input=unfold.GetEmatrixInput("cov_input",";bin",0,0,false);
               TH2 *cov_statResponse=unfold.GetEmatrixSysUncorr("cov_statResponse",";bin",0,0,true);
               // ~TH2 *cov_Output=unfold.GetEmatrixTotal("cov_total",";bin",0,0,true);
               TH2 *cov_Output=unfold.GetEmatrixTotal("cov_total",";bin",0,0,false);
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
                  TH2 *cov_input_firstToy=0;
                  for(int itoy=0;itoy<MAXTOY;itoy++) {
                     std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
                     TH1 *reco_toy=generatePoissonToy(histDataReco_unscaled,itoy);
                     
                     if (itoy==0)unfold.SetInput(histDataReco_unscaled,1.0);
                     else unfold.SetInput(reco_toy,1.0);
                     
                     if (regularisation) unfold.DoUnfold(tau); //With regularization
                     else unfold.DoUnfold(0.);//Without regularization
                     
                     if (itoy==0){
                        hist_unfolded_firstToy=unfold.GetOutput("",";bin",0,0,false);
                        cov_input_firstToy=unfold.GetEmatrixInput("",";bin",0,0,false);
                     }
                     else analyzeToy(unfold.GetOutput("","",0,0,false),
                                   // ~hist_unfolded,
                                   hist_unfolded_firstToy,
                                   // ~histDataTruth,
                                   // ~cov_Output,
                                   // ~cov_input,
                                   cov_input_firstToy,
                                   prof_pull,
                                   hist_coverage,hist_pull,prof_res,hist_res,hist_chi,MAXTOY-1);
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
               for (int i=1; i<=hist_unfolded->GetNbinsX();i++){
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
               if (regularisation) {
                  TParameter<float> tau_par("tau",tau);
                  saver.save(tau_par,saveFolder+"tau");
               }
            }//end unfolding loop
         }
         
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
         
         TH1D* BBB_factors = histMCGenRec_sameBins->ProjectionX();
         BBB_factors->Divide(histMCGenRec_sameBins->ProjectionY());
         saver.save(*BBB_factors,"BBB/correctionFactors");
         BBB_factors->Multiply(histSignalFraction_coarse);
         saver.save(*BBB_factors,"BBB/correctionFactorsTimesSF");
         
         // ~TH1 *histSignalFractionAlt_coarse=histReader.read<TH1>(input_loc+"/hist_SignalFractionAlt_coarse");
         // ~TH2* histMCGenRecAlt_sameBins = histReader_syst.read<TH2>(input_loc+"/histMCGenRecAlt_sameBins");
         // ~TH1D* BBB_factorsAlt = histMCGenRecAlt_sameBins->ProjectionX();
         // ~BBB_factorsAlt->Divide(histMCGenRecAlt_sameBins->ProjectionY());
         // ~saver.save(*BBB_factorsAlt,"BBB/correctionFactorsAlt");
         // ~BBB_factorsAlt->Multiply(histSignalFractionAlt_coarse);
         // ~saver.save(*BBB_factorsAlt,"BBB/correctionFactorsTimesSFAlt");
         
         bbb.PrintTable(std::cout);
             
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
            TH2 *covBBB_firstToy;
            for(int itoy=0;itoy<MAXTOY;itoy++) {
               // ~std::cout<<"================== itoy="<<itoy<<" =========================="<<std::endl;
               TH1 *reco_toy=generatePoissonToy(histDataReco_coarse_unscaled,itoy);
               
               if (itoy==0){
                  RooUnfoldBinByBin bbb(&responseRooUnfold, histDataReco_coarse_unscaled);
                  hist_unfolded_firstToy=bbb.Hreco(RooUnfold::kCovariance);
                  TMatrixD covBBB_matrix=bbb.Ereco();
                  covBBB_firstToy = new TH2D(covBBB_matrix);
               }
               else {
                  RooUnfoldBinByBin bbb(&responseRooUnfold, reco_toy);
                  analyzeToy(bbb.Hreco(RooUnfold::kCovariance),
                             // ~histBBB,
                             hist_unfolded_firstToy,
                             // ~covBBB,
                             covBBB_firstToy,
                             prof_pull,
                             hist_coverage,hist_pull,prof_res,hist_res,hist_chi,MAXTOY);
               }
            }
            saver.save(*hist_unfolded_firstToy,"BBB/hist_unfolded_firstToy");
            saver.save(*prof_pull,"BBB/prof_pull");
            saver.save(*prof_res,"BBB/prof_res");
            saver.save(*hist_coverage,"BBB/hist_coverage");
            saver.save(*hist_pull,"BBB/hist_pull");
            saver.save(*hist_res,"BBB/hist_res");
            saver.save(*hist_chi,"BBB/hist_chi");
         }
         
         // Condition number of response matrix
         
         TH2F* tempHist=(TH2F*)histMCGenRec->Clone();
         // ~TH2F* tempHist=(TH2F*)histMCGenRec_sameBins->Clone();
         
         float sum = 0;
         for (int x=1; x<=tempHist->GetNbinsX(); x++){   //Normalize to gen bin
            sum=tempHist->Integral(x,x,0,tempHist->GetNbinsY());
            if (sum==0) continue;
            for (int y=1; y<=tempHist->GetNbinsY(); y++){
               if (tempHist->GetBinContent(x,y)!=0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
               else tempHist->SetBinContent(x,y,0.0);
            }
         }
         
         TMatrixD responseMatrix(tempHist->GetNbinsY(),tempHist->GetNbinsX());
         // ~TMatrixD responseMatrix(28,7);
         for(int i=1; i<=tempHist->GetNbinsX(); i++){
            for(int j=1; j<=tempHist->GetNbinsY(); j++){
         // ~for(int i=15; i<=tempHist->GetNbinsX(); i++){
            // ~for(int j=57; j<=tempHist->GetNbinsY(); j++){
               responseMatrix(j-1,i-1) = tempHist->GetBinContent(i,j);
               // ~responseMatrix(j-57,i-15) = tempHist->GetBinContent(i,j);
            }
         }             
                    
         TDecompSVD svdMatrix(responseMatrix);
         std::cout<<"======================================"<<std::endl;
         std::cout<<"Condition Number of Matrix: "<<svdMatrix.Condition()<<std::endl;
         std::cout<<"======================================"<<std::endl;
      }
      if (!verbose) std::cout.clear();;
   }
}
