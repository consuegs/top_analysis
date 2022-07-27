//Script to define unfolding binning and generate distributions used for the unfolding based on TUnfolding classes (Following TUnfold examples)
//This Script takes fakes into account by correcting the signal fraction before unfolding
#include <iostream>
#include <map>
#include <cmath>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile.h>
#include <tuple>
#include "TUnfoldBinning.h"

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;
Config const &cfg=Config::get();

std::tuple<TUnfoldBinning*,TUnfoldBinning*> defineBinnings2D(std::vector<double> &metBinsCoarse_vector, std::vector<double> &phiBinsCoarse_vector, 
                                                            std::vector<double> &metBinsFine_vector, std::vector<double> &phiBinsFine_vector){
   
   // Get nBins and convert to arrays needed for TUnfoldBinning
   int NBIN_MET_COARSE = metBinsCoarse_vector.size()-1;
   int NBIN_PHI_COARSE = phiBinsCoarse_vector.size()-1;
   Double_t* metBinsCoarse = &metBinsCoarse_vector[0];
   Double_t* phiBinsCoarse = &phiBinsCoarse_vector[0];
   
   if(cfg.tunfold_withSameBins){
      metBinsFine_vector.resize(NBIN_MET_COARSE+1);
      phiBinsFine_vector.resize(NBIN_PHI_COARSE+1);
      metBinsFine_vector = metBinsCoarse_vector;
      phiBinsFine_vector = phiBinsCoarse_vector;
   }
   int NBIN_MET_FINE = metBinsFine_vector.size()-1;
   int NBIN_PHI_FINE = phiBinsFine_vector.size()-1;
   Double_t* metBinsFine=&metBinsFine_vector[0];
   Double_t* phiBinsFine=&phiBinsFine_vector[0];
   
   //=======================================================================
   // detector level binning scheme

   TUnfoldBinning *detectorDistribution=new TUnfoldBinning("detector");
   TUnfoldBinning *detectorBinning=detectorDistribution->AddBinning("detectordistribution");
   
   detectorBinning->AddAxis("met",NBIN_MET_FINE,metBinsFine,
                                 false, // no underflow bin (not reconstructed)
                                 true // overflow bin
                                 );
   detectorBinning->AddAxis("phi",NBIN_PHI_FINE,phiBinsFine,
                                 false, // no underflow bin (not reconstructed)
                                 false // no overflow bin (not reconstructed)
                                 );

   //=======================================================================
   // generator level binning
   TUnfoldBinning *generatorBinning=new TUnfoldBinning("generator");

   // signal distribution is measured with coarse binning
   // fake bin is used to check
   // what happens outside the phase-space
   TUnfoldBinning *signalBinning = generatorBinning->AddBinning("signal");
   signalBinning->AddAxis("metnunugen",NBIN_MET_COARSE,metBinsCoarse,
                           false, // underflow bin
                           true // overflow bin
                           );
   signalBinning->AddAxis("phigen",NBIN_PHI_COARSE,phiBinsCoarse,
                           false, // underflow bin
                           false // overflow bin
                           );
   
   return {detectorBinning,signalBinning};
}

std::tuple<TUnfoldBinning*,TUnfoldBinning*> defineBinnings1D(std::vector<double> &BinsCoarse_vector, std::vector<double> &BinsFine_vector){
   
   // Get nBins and convert to arrays needed for TUnfoldBinning
   int NBIN_COARSE = BinsCoarse_vector.size()-1;
   Double_t* BinsCoarse = &BinsCoarse_vector[0];
   
   if(cfg.tunfold_withSameBins){
      BinsFine_vector.resize(NBIN_COARSE+1);
      BinsFine_vector = BinsCoarse_vector;
   }
   int NBIN_FINE = BinsFine_vector.size()-1;
   Double_t* BinsFine=&BinsFine_vector[0];
   
   //=======================================================================
   // detector level binning scheme
   
   TUnfoldBinning *detectorDistribution=new TUnfoldBinning("detector");
   TUnfoldBinning *detectorBinning=detectorDistribution->AddBinning("detectordistribution");
   detectorBinning->AddAxis("reco",NBIN_FINE,BinsFine,
                                 false, // no underflow bin (not reconstructed)
                                 // ~true // overflow bin
                                 !(NBIN_COARSE == 1) // overflow bin
                                 );

   //=======================================================================
   // generator level binning
   TUnfoldBinning *generatorBinning=new TUnfoldBinning("generator");

   // signal distribution is measured with coarse binning
   // fake bin is used to check
   // what happens outside the phase-space
   TUnfoldBinning *signalBinning = generatorBinning->AddBinning("signal");
   signalBinning->AddAxis("gen",NBIN_COARSE,BinsCoarse,
                           false, // underflow bin
                           // ~true // overflow bin
                           !(NBIN_COARSE == 1) // overflow bin
                           );
   
   return {detectorBinning,signalBinning};
}

TH1* deriveSignalFraction(TH2* const histMCGenRec, TH1* const histMCRec_fakes, TUnfoldBinning* const &detectorBinning, io::RootFileSaver const &saver,
                           TString const &histName, TString const &varName){
   TH1 *hist_SignalFraction=detectorBinning->CreateHistogram(histName);
   TH1 *histMCRec=histMCGenRec->ProjectionY("",1,-1);
   saver.save(*histMCRec,varName+"/histMCRec");
   
   for(int i=1; i<=histMCRec->GetNbinsX(); i++){  //Is this correct or not???? Idea: mc stat. of ttbar signal is already taken into account in response matrix
      histMCRec->SetBinError(i,0);
   }
   
   hist_SignalFraction->Add(histMCRec);
   histMCRec->Add(histMCRec_fakes);
   hist_SignalFraction->Divide(histMCRec);
   
   return hist_SignalFraction;
}

TH1* deriveEfficiency(TH2* const histMCGenRec, TUnfoldBinning* const &signalBinning){
   TH1 *hist_efficiency=signalBinning->CreateHistogram("efficiency");
   for(int binGen=0;binGen<=histMCGenRec->GetNbinsX()+1;binGen++) {
      double sum0=0.;
      double sum1=0.;
      for(int binRec=0;binRec<=histMCGenRec->GetNbinsY()+1;binRec++) {
         double c=histMCGenRec->GetBinContent(binGen,binRec);
         sum0+=c;
         if((binRec>0)&&(binRec<=histMCGenRec->GetNbinsY()+1)) sum1+=c;
      } 
      if(sum0>0.0) hist_efficiency->SetBinContent(binGen,sum1/sum0);
   }
   return hist_efficiency;
}

void derivePurityStability(TH2* const histMCGenRec_purity, TUnfoldBinning* const &signalBinning, io::RootFileSaver const &saver, TString const &varName){
   TH1 *hist_purity=signalBinning->CreateHistogram("purity");
   TH1 *hist_stability=signalBinning->CreateHistogram("stability");
   TH1 *hist_N_GenRec=signalBinning->CreateHistogram("N_GenRec");
   TH1 *hist_N_Rec=signalBinning->CreateHistogram("N_Rec");
   for(int binRec=0;binRec<=hist_purity->GetNbinsX()+1;binRec++) {
      double sum=0.;
      for(int binGen=1;binGen<=hist_purity->GetNbinsX()+1;binGen++) {
         sum += histMCGenRec_purity->GetBinContent(binGen,binRec);
      }
      double p=0.;
      if(sum>0.0) {
         p=histMCGenRec_purity->GetBinContent(binRec,binRec)/sum;
      }
      hist_purity->SetBinContent(binRec,p);
      hist_N_GenRec->SetBinContent(binRec,histMCGenRec_purity->GetBinContent(binRec,binRec));
      hist_N_Rec->SetBinContent(binRec,sum);
   }
   
   saver.save(*hist_purity,varName+"/hist_purity");
   saver.save(*hist_N_GenRec,varName+"/hist_N_GenRec");
   saver.save(*hist_N_Rec,varName+"/hist_N_Rec");
   
   for(int binGen=0;binGen<=hist_purity->GetNbinsX()+1;binGen++) {
      double sum=0.;
      for(int binRec=1;binRec<=hist_purity->GetNbinsX()+1;binRec++) {
         sum += histMCGenRec_purity->GetBinContent(binGen,binRec);
      }
      double p=0.;
      if(sum>0.0) {
         p=histMCGenRec_purity->GetBinContent(binGen,binGen)/sum;
      }
      hist_stability->SetBinContent(binGen,p);
   }
   
   saver.save(*hist_stability,varName+"/hist_stability");
}

//small class containing all TUnfold related hist for one distribution
class Distribution
{
   public:
      // 2D constructor
      Distribution(TString const &varName,std::vector<double> xBinsCoarse_vector, std::vector<double> yBinsCoarse_vector, 
                  std::vector<double> xBinsFine_vector, std::vector<double> yBinsFine_vector):
      varName_(varName)
      {
         is2D = true;
         std::tie(detectorBinning_,signalBinning_) = defineBinnings2D(xBinsCoarse_vector,yBinsCoarse_vector,xBinsFine_vector,yBinsFine_vector);
      }
      
      // 1D constructor
      Distribution(TString const &varName,std::vector<double> xBinsCoarse_vector, std::vector<double> xBinsFine_vector):
      varName_(varName)
      {
         is2D = false;
         std::tie(detectorBinning_,signalBinning_) = defineBinnings1D(xBinsCoarse_vector,xBinsFine_vector);
      }
      
      void setupDataHists(){
         histDataReco=detectorBinning_->CreateHistogram("histDataReco_"+varName_);
         histDataReco_coarse=signalBinning_->CreateHistogram("histDataReco_coarse_"+varName_);
         histDataTruth=signalBinning_->CreateHistogram("histDataTruth_"+varName_);
         histDataTruth_fakes=signalBinning_->CreateHistogram("histDataTruth_fakes_"+varName_);
         histDataRecoAlt=detectorBinning_->CreateHistogram("histDataRecoAlt_"+varName_);
         histDataRecoAlt_coarse=signalBinning_->CreateHistogram("histDataRecoAlt_coarse_"+varName_);
         histDataTruthAlt=signalBinning_->CreateHistogram("histDataTruthAlt_"+varName_);
         histDataTruthAlt_fakes=signalBinning_->CreateHistogram("histDataTruthAlt_fakes_"+varName_);
      }
      
      void setupMCHists(){
         histMCGenRec=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning_,detectorBinning_,"histMCGenRec-"+varName_);
         histMCGenRec_purity=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning_,signalBinning_,"histMCGenRec_purity-"+varName_);
         histMCGenRecAlt=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning_,detectorBinning_,"histMCGenRecAlt-"+varName_);
         histMCGenRecAlt_purity=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning_,signalBinning_,"histMCGenRecAlt_purity-"+varName_);
         histMCRec_fakes=detectorBinning_->CreateHistogram("histMCRec_fakes-"+varName_);
         histMCRec_fakes_coarse=signalBinning_->CreateHistogram("histMCRec_fakes_coarse-"+varName_);
         histMCRecAlt_fakes=detectorBinning_->CreateHistogram("histMCRecAlt_fakes-"+varName_);
         histMCRecAlt_fakes_coarse=signalBinning_->CreateHistogram("histMCRecAlt_fakes_coarse-"+varName_);
         histMCRec_bkg=detectorBinning_->CreateHistogram("histMCRec_bkg-"+varName_);
         histMCRec_bkg_coarse=signalBinning_->CreateHistogram("histMCRec_bkg_coarse-"+varName_);
      }
      
      //2D
      void setVariables(Float_t &xReco, Float_t &yReco, Float_t &xGen, Float_t &yGen){
         xReco_ = &xReco;
         yReco_ = &yReco;
         xGen_ = &xGen;
         yGen_ = &yGen;
      }
      
      //1D
      void setVariables(Float_t &xReco, Float_t &xGen){
         xReco_ = &xReco;
         xGen_ = &xGen;
      }
      
      void fillDataTruth_fakes(float const &weight){
         Int_t genbinNumber = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         histDataTruth_fakes->Fill(genbinNumber,weight);
      }
      
      void fillDataTruth(float const &weight){
         Int_t genbinNumber = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         histDataTruth->Fill(genbinNumber,weight);
      }
      
      void fillDataTruthAlt_fakes(float const &weight){
         Int_t genbinNumber = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         histDataTruthAlt_fakes->Fill(genbinNumber,weight);
      }
      
      void fillDataTruthAlt(float const &weight){
         Int_t genbinNumber = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         histDataTruthAlt->Fill(genbinNumber,weight);
      }
      
      void fillDataReco(float const &weight){
         Int_t binNumber = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         histDataReco->Fill(binNumber,weight);
         
         Int_t binNumber_coarse = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         histDataReco_coarse->Fill(binNumber_coarse,weight);
      }
      
      void fillDataRecoAlt(float const &weight){
         Int_t binNumber = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         histDataRecoAlt->Fill(binNumber,weight);
         
         Int_t binNumber_coarse = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         histDataRecoAlt_coarse->Fill(binNumber_coarse,weight);
      }
      
      void fillMCbkg(float const &weight){
         Int_t recBin = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         Int_t recBin_purity = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         histMCRec_bkg->Fill(recBin,weight);
         histMCRec_bkg_coarse->Fill(recBin_purity,weight);
      }
      
      void fillMCfakes(float const &weight){
         Int_t recBin = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         Int_t recBin_purity = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         histMCRec_fakes->Fill(recBin,weight);
         histMCRec_fakes_coarse->Fill(recBin_purity,weight);
      }
      
      void fillMCfakesAlt(float const &weight){
         Int_t recBin = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         Int_t recBin_purity = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         histMCRecAlt_fakes->Fill(recBin,weight);
         histMCRecAlt_fakes_coarse->Fill(recBin_purity,weight);
      }
      
      void fillMCGenRec(float const &mcWeight,float const &recoWeight){
         Int_t recBin = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         Int_t recBin_purity = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         Int_t genBin = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         
         histMCGenRec->Fill(genBin,recBin,mcWeight*recoWeight);
         histMCGenRec_purity->Fill(genBin,recBin_purity,mcWeight*recoWeight);
         
         // count fraction of events which have a different weight on truth and reco (in reco underflow bin)
         // this is required for TUnfold to function properly
         histMCGenRec->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
         histMCGenRec_purity->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
      }
      
      void fillMCGenRecAlt(float const &mcWeight,float const &recoWeight){
         Int_t recBin = (is2D)? detectorBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : detectorBinning_->GetGlobalBinNumber(*xReco_);
         Int_t recBin_purity = (is2D)? signalBinning_->GetGlobalBinNumber(*xReco_,*yReco_) : signalBinning_->GetGlobalBinNumber(*xReco_);
         Int_t genBin = (is2D)? signalBinning_->GetGlobalBinNumber(*xGen_,*yGen_) : signalBinning_->GetGlobalBinNumber(*xGen_);
         
         histMCGenRecAlt->Fill(genBin,recBin,mcWeight*recoWeight);
         histMCGenRecAlt_purity->Fill(genBin,recBin_purity,mcWeight*recoWeight);
         
         // count fraction of events which have a different weight on truth and reco (in reco underflow bin)
         // this is required for TUnfold to function properly
         histMCGenRecAlt->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
         histMCGenRecAlt_purity->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
      }
      
      void setExpStatUnc(){
         // set reco bin error to error expected in data
         for(int i=0;i<=histDataReco->GetNbinsX()+1;i++) {
            histDataReco->SetBinError(i,sqrt(histDataReco->GetBinContent(i)));
            histDataRecoAlt->SetBinError(i,sqrt(histDataRecoAlt->GetBinContent(i)));
         }
         for(int i=0;i<=histDataReco_coarse->GetNbinsX()+1;i++) {
            histDataReco_coarse->SetBinError(i,sqrt(histDataReco_coarse->GetBinContent(i)));
            histDataRecoAlt_coarse->SetBinError(i,sqrt(histDataRecoAlt_coarse->GetBinContent(i)));
         }
      }
      
      void saveDataHists(io::RootFileSaver const &saver){
         saver.save(*detectorBinning_,varName_+"/detector_binning");
         saver.save(*signalBinning_,varName_+"/generator_binning");
         saver.save(*histDataReco,varName_+"/histDataReco");
         saver.save(*histDataReco_coarse,varName_+"/histDataReco_coarse");
         saver.save(*histDataRecoAlt,varName_+"/histDataRecoAlt");
         saver.save(*histDataRecoAlt_coarse,varName_+"/histDataRecoAlt_coarse");
         saver.save(*histDataTruth_fakes,varName_+"/histDataTruth_fakes");
         saver.save(*histDataTruth,varName_+"/histDataTruth");
         saver.save(*histDataTruthAlt_fakes,varName_+"/histDataTruthAlt_fakes");
         saver.save(*histDataTruthAlt,varName_+"/histDataTruthAlt");
      }
      
      void saveMCHists(io::RootFileSaver const &saver){
         saver.save(*histMCGenRec,varName_+"/histMCGenRec");
         saver.save(*histMCGenRec_purity,varName_+"/histMCGenRec_sameBins");
         saver.save(*histMCGenRecAlt,varName_+"/histMCGenRecAlt");
         saver.save(*histMCGenRecAlt_purity,varName_+"/histMCGenRecAlt_sameBins");
         saver.save(*histMCRec_fakes,varName_+"/histMCRec_fakes");
         saver.save(*histMCRecAlt_fakes,varName_+"/histMCRecAlt_fakes");
         saver.save(*histMCRec_bkg,varName_+"/histMCRec_bkg");
         saver.save(*histMCRec_bkg_coarse,varName_+"/histMCRec_bkg_coarse");
         saver.save(*(histMCGenRec->ProjectionX()),varName_+"/histMCGenRec_projX");
         saver.save(*(histMCGenRec->ProjectionY()),varName_+"/histMCGenRec_projY");
         saver.save(*(histMCGenRecAlt->ProjectionX()),varName_+"/histMCGenRecAlt_projX");
         saver.save(*(histMCGenRecAlt->ProjectionY()),varName_+"/histMCGenRecAlt_projY");
         
         //Calculate Signal fration
         TH1* hist_SignalFraction = deriveSignalFraction(histMCGenRec,histMCRec_fakes,detectorBinning_,saver,"hist_SignalFraction",varName_);
         saver.save(*hist_SignalFraction,varName_+"/hist_SignalFraction");
         
         //Calculate Signal fration for coarse binning
         TH1* hist_SignalFraction_coarse = deriveSignalFraction(histMCGenRec_purity,histMCRec_fakes_coarse,signalBinning_,saver,"hist_SignalFraction_coarse",varName_);
         saver.save(*hist_SignalFraction_coarse,varName_+"/hist_SignalFraction_coarse");
         
         //Calculate Signal fration for alt sample
         TH1* hist_SignalFractionAlt = deriveSignalFraction(histMCGenRecAlt,histMCRecAlt_fakes,detectorBinning_,saver,"hist_SignalFractionAlt",varName_);
         saver.save(*hist_SignalFractionAlt,varName_+"/hist_SignalFractionAlt");
         
         //Calculate Signal fration for coarse binning for alt sample
         TH1* hist_SignalFractionAlt_coarse = deriveSignalFraction(histMCGenRecAlt_purity,histMCRecAlt_fakes_coarse,signalBinning_,saver,"hist_SignalFractionAlt_coarse",varName_);
         saver.save(*hist_SignalFractionAlt_coarse,varName_+"/hist_SignalFractionAlt_coarse");
         
         //Save normalized reco distribution for fake and non fake events
         TH1D histGen_fakes=*(histMCGenRec_purity->ProjectionY("Fakes",1,1));
         TH1D histGen_noFakes=*(histMCGenRec_purity->ProjectionY("Fakes",2,-1));
         saver.save(histGen_fakes,varName_+"/histGen_fakes");
         saver.save(histGen_noFakes,varName_+"/histGen_noFakes");
         
         //efficiency
         TH1* hist_efficiency = deriveEfficiency(histMCGenRec,signalBinning_);
         saver.save(*hist_efficiency,varName_+"/hist_efficiency");
         
         // calculate bin purities and stabilities
         derivePurityStability(histMCGenRec_purity,signalBinning_,saver,varName_);
      }
      
      TUnfoldBinning* detectorBinning_;
      TUnfoldBinning* signalBinning_;
      TString varName_;
      bool is2D;
      
      Float_t* xGen_;
      Float_t* yGen_;
      Float_t* xReco_;
      Float_t* yReco_;
      
      TH1* histDataReco;
      TH1* histDataReco_coarse;
      TH1* histDataRecoAlt;
      TH1* histDataRecoAlt_coarse;
      TH1* histDataTruth;
      TH1* histDataTruth_fakes;
      TH1* histDataTruthAlt;
      TH1* histDataTruthAlt_fakes;
      
      TH2* histMCGenRec;
      TH2* histMCGenRec_purity;
      TH2* histMCGenRecAlt;
      TH2* histMCGenRecAlt_purity;
      TH1* histMCRec_fakes;
      TH1* histMCRec_fakes_coarse;
      TH1* histMCRecAlt_fakes;
      TH1* histMCRecAlt_fakes_coarse;
      TH1* histMCRec_bkg;
      TH1* histMCRec_bkg_coarse;
};

TString getSavePath(){
   TString save_path = "TUnfold_binning_"+cfg.tunfold_InputSamples[0]+"_"+cfg.tunfold_ResponseSample;
   if (cfg.tunfold_withBSM) save_path+="_BSM";
   if (!cfg.tunfold_withDNN && cfg.tunfold_withPF) save_path+="_PF";
   else if (!cfg.tunfold_withDNN) save_path+="_Puppi";
   if (cfg.tunfold_withDNN) save_path+="_DNN";
   if (cfg.tunfold_withSameBins) save_path+="_SameBins";
   if (cfg.tunfold_withPTreweight) {
      save_path+="_PTreweight"+cfg.tunfold_scalePTreweight;
   }
   return save_path;
}

void loopDataEvents(std::vector<Distribution> &distribution_vec, io::RootFileSaver const &saver){
   
   TString sample = cfg.tunfold_InputSamples[0];
   
   if (cfg.tunfold_withPTreweight) {
      sample+="_PTreweight"+cfg.tunfold_scalePTreweight;
   }
   
   // define variables for reading tree
   Float_t phiRec,metRec,phiRec_DNN,metRec_DNN,phiGen,metGen,pTllRec,pTllGen,mcWeight,recoWeight,genMET,reweight;
   UInt_t genDecayMode;
   
   // setup distributions
   for(Distribution& dist : distribution_vec) {
      dist.setupDataHists();
      if (dist.varName_ == "2D_dPhi_pTnunu") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new40") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new40_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "pTnunu") dist.setVariables(metRec,metGen);
      else if (dist.varName_ == "dPhi") dist.setVariables(phiRec,phiGen);
      else if (dist.varName_ == "pTnunu_DNN") dist.setVariables(metRec_DNN,metGen);
      else if (dist.varName_ == "pTnunu_new_DNN") dist.setVariables(metRec_DNN,metGen);
      else if (dist.varName_ == "dPhi_DNN") dist.setVariables(phiRec_DNN,phiGen);
      else if (dist.varName_ == "pTll") dist.setVariables(pTllRec,pTllGen);      
      else if (dist.varName_ == "inclusive") dist.setVariables(phiRec,phiGen);      
   }
   
   TString minTreePath_Nominal = cfg.minTreePath+TString::Format("/100.0/%s/","Nominal");
   
   for (TString currentSample : cfg.tunfold_InputSamples){     // loop over samples (data or pseudo data)
      bool isSignal = (currentSample == "TTbar_diLepton");
      bool isSignalAlt = (currentSample == "TTbar_amcatnlo");
      for (int i=1; i<=cfg.getTotalFileNR(currentSample.Data()); i++){     // loop over each minTree of sample
         TFile *dataFile=new TFile(minTreePath_Nominal+currentSample+"_"+std::to_string(i)+".root","read");
         TTree *dataTree=(TTree *) dataFile->Get("ttbar_res100.0/"+currentSample);

         if(!dataTree) {
            cout<<"could not read 'data' tree\n";
         }
         
         // set tree variables
         dataTree->ResetBranchAddresses();
         if(cfg.tunfold_withPF) {
            dataTree->SetBranchAddress("Phi_rec",&phiRec);
            dataTree->SetBranchAddress("MET",&metRec);
         }
         else if(cfg.tunfold_withDNN) {
            dataTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
            dataTree->SetBranchAddress("DNN_MET_pT",&metRec);
         }
         else{
            dataTree->SetBranchAddress("Phi_recPuppi",&phiRec);
            dataTree->SetBranchAddress("PuppiMET",&metRec);
         }
         dataTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec_DNN);
         dataTree->SetBranchAddress("DNN_MET_pT",&metRec_DNN);
         dataTree->SetBranchAddress("Phi_NuNu",&phiGen);
         dataTree->SetBranchAddress("PtNuNu",&metGen);
         dataTree->SetBranchAddress("pTll",&pTllRec);
         dataTree->SetBranchAddress("pTllGen",&pTllGen);
         dataTree->SetBranchAddress("genDecayMode",&genDecayMode);
         dataTree->SetBranchAddress("N",&mcWeight);
         dataTree->SetBranchAddress("SF",&recoWeight);
         if(cfg.tunfold_withPTreweight) dataTree->SetBranchAddress("reweight_PTnunu",&reweight);
         else reweight=1.0;
         

         std::cout<<"loop over "<<currentSample<<i<<std::endl;
         
         int entriesToRun=dataTree->GetEntriesFast()*cfg.processFraction;
         int iEv=0;
         for(Int_t ievent=0;ievent<entriesToRun;ievent++) {    // loop over events
            iEv++;
            if (iEv%(std::max(entriesToRun/10,1))==0){
               io::log*".";
               io::log.flush();
            }
            
            if(dataTree->GetEntry(ievent)<=0) break;
            
            // ~mcWeight=mcWeight*(137191.0/35867.05998);   //scale to full Run2 Lumi

            //only bin to bin migration
            // ~if(metRec<0 || genDecayMode>3 || metGen<0) continue;
            
            //ignore acceptance
            // ~if(metRec<0) continue;
            
            //ignore fakes
            // ~if(genDecayMode>3 || metGen<0) continue;
            
            // ~if (metGen>0) mcWeight*=sqrt(sqrt(sqrt(sqrt(metGen))));
            
            //reweight
            mcWeight*=reweight;
            
            if (isSignal){    // fill truth distributions (used for pseudo data)
               if (metGen<0 || genDecayMode>3 || (genDecayMode!=3 && metGen<40)) {
                  for(Distribution& dist : distribution_vec) dist.fillDataTruth_fakes(mcWeight);
               }
               else for(Distribution& dist : distribution_vec) dist.fillDataTruth(mcWeight);
            }
            else if (isSignalAlt){ // fill truth distributions with alternative MC sample (used for pseudo data)
               if (metGen<0 || genDecayMode>3 || (genDecayMode!=3 && metGen<40)) {
                  for(Distribution& dist : distribution_vec) dist.fillDataTruthAlt_fakes(mcWeight);
               }
               else for(Distribution& dist : distribution_vec) dist.fillDataTruthAlt(mcWeight);
            }

            // fill histogram with reconstructed quantities
            if (metRec<0) continue;   //events that are not reconstructed
            if (!isSignalAlt){
               for(Distribution& dist : distribution_vec) dist.fillDataReco(mcWeight*recoWeight);
            }
            if (!isSignal){
               for(Distribution& dist : distribution_vec) dist.fillDataRecoAlt(mcWeight*recoWeight);
            }
         }
         io::log<<"";
      }
   }
   
   /*
   TFile *bsmFile=new TFile(minTreePath_Nominal+"T2tt_650_350.root","read");
   TTree *BSMTree=(TTree *) bsmFile->Get("ttbar_res100.0/T2tt_650_350");
   if(!BSMTree) {
      cout<<"could not read 'BSM' tree\n";
   }
   
   if (cfg.tunfold_withBSM) {
      BSMTree->ResetBranchAddresses();
      BSMTree->SetBranchAddress("Phi_recPuppi",&phiRec);
      BSMTree->SetBranchAddress("PuppiMET",&metRec);
      if(cfg.tunfold_withPF) {
         BSMTree->SetBranchAddress("Phi_rec",&phiRec);
         BSMTree->SetBranchAddress("MET",&metRec);
      }
      else if(cfg.tunfold_withDNN) {
         BSMTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
         BSMTree->SetBranchAddress("DNN_MET_pT",&metRec);
      }
      else{
         BSMTree->SetBranchAddress("Phi_recPuppi",&phiRec);
         BSMTree->SetBranchAddress("PuppiMET",&metRec);
      }
      BSMTree->SetBranchAddress("Phi_NuNu",&phiGen);
      BSMTree->SetBranchAddress("PtNuNu",&metGen);
      BSMTree->SetBranchAddress("genDecayMode",&genDecayMode);
      BSMTree->SetBranchAddress("N",&mcWeight);
      BSMTree->SetBranchAddress("SF",&recoWeight);
      
      int entriesToRun=BSMTree->GetEntriesFast()*cfg.processFraction;
      
      for(Int_t ievent=0;ievent<entriesToRun;ievent++) {
         if(BSMTree->GetEntry(ievent)<=0) break;
         
         // ~mcWeight=mcWeight*(137191.0/35867.05998);   //scale to full Run2 Lumi

         // fill histogram with reconstructed quantities
         if (metRec<0) continue;   //events that are not reconstructed
         Int_t binNumber=detectorBinning->GetGlobalBinNumber(metRec,phiRec);
         histDataReco->Fill(binNumber,mcWeight);
         
         Int_t binNumber_coarse=signalBinning->GetGlobalBinNumber(metRec,phiRec);
         histDataReco_coarse->Fill(binNumber_coarse,mcWeight);
      }
   }
   delete BSMTree;
   */
   
   // correct stat unc. (if pseudo data) and save data hists
   for(Distribution& dist : distribution_vec){
      dist.setExpStatUnc();
      dist.saveDataHists(saver);
   }
}

std::tuple<TString,TString,float> getPath_SampleName_SF(TString const &sample, TString const &minTreePath_nominal, TString const &minTreePath_syst, bool const isSignal, bool const                                                isBKGother, Systematic::Systematic const &syst){
   // Change sample name in case systematic has own sample. Use nominal for non ttbar samples
   TString minTreePath_current = minTreePath_syst;
   TString currentSample = sample;
   float SF = 1.0;
   if(std::find(Systematic::altSampleTypes.begin(),Systematic::altSampleTypes.end(), syst.type()) != Systematic::altSampleTypes.end()){
      if(!isBKGother){
         currentSample+="_"+syst.name();
      }
      else{
         minTreePath_current = minTreePath_nominal;
      }
   }
   // Use nominal for all except ttbar in case of pdf unc.
   else if(syst.type()==Systematic::pdf){
      if(!isSignal) minTreePath_current = minTreePath_nominal;
   }
   // Use nominal for other bkg in case of bFrag bSemi and MEScale (currently not working for madgraph samples)
   else if(syst.type() == Systematic::bFrag || syst.type() == Systematic::bSemilep || syst.type() == Systematic::meFacScale || syst.type() == Systematic::meRenScale || syst.type() == Systematic::meScale){
      if(isBKGother) minTreePath_current = minTreePath_nominal;
   }
   // Change SF for lumi and xsec uncertainties
   else if(cfg.systUncFactor.find(syst.type_str()) != cfg.systUncFactor.end()){
      minTreePath_current = minTreePath_nominal;
      std::vector<std::string> samplesSF = cfg.systUncFactor.at(syst.type_str()).second;
      if (std::find(samplesSF.begin(),samplesSF.end(),sample) != samplesSF.end()){
         float unc = cfg.systUncFactor.at(syst.type_str()).first;
         SF = (syst.variation() == Systematic::up)? 1.0+unc : 1.0-unc;
      }
   }
   
   return {minTreePath_current,currentSample,SF};
}

void loopMCEvents(std::vector<Distribution> &distribution_vec, io::RootFileSaver const &saver, Systematic::Systematic const &syst){
   //set correct paths for nominal and current syst
   TString minTreePath_nominal = cfg.minTreePath+TString::Format("/100.0/%s/","Nominal");
   TString minTreePath_syst = cfg.minTreePath+TString::Format("/100.0/%s/",syst.name().Data());
   
   // define variables for reading tree
   Float_t phiRec,metRec,phiRec_DNN,metRec_DNN,phiGen,metGen,pTllRec,pTllGen,mcWeight,recoWeight,genMET,reweight;
   UInt_t genDecayMode;
   
   // setup distributions
   for(Distribution& dist : distribution_vec) {
      dist.setupMCHists();
      if (dist.varName_ == "2D_dPhi_pTnunu") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new40") dist.setVariables(metRec,phiRec,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "2D_dPhi_pTnunu_new40_DNN") dist.setVariables(metRec_DNN,phiRec_DNN,metGen,phiGen);
      else if (dist.varName_ == "pTnunu") dist.setVariables(metRec,metGen);
      else if (dist.varName_ == "dPhi_DNN") dist.setVariables(phiRec_DNN,phiGen);      
      else if (dist.varName_ == "pTnunu_DNN") dist.setVariables(metRec_DNN,metGen);
      else if (dist.varName_ == "pTnunu_new_DNN") dist.setVariables(metRec_DNN,metGen);
      else if (dist.varName_ == "dPhi") dist.setVariables(phiRec,phiGen);      
      else if (dist.varName_ == "pTll") dist.setVariables(pTllRec,pTllGen);      
      else if (dist.varName_ == "inclusive") dist.setVariables(phiRec,phiGen);      
   }
   
   // create vector with all input samples
   std::vector<std::string> samples = cfg.tunfold_bkgSamples_ttbar;
   samples.insert(samples.begin(),cfg.tunfold_ResponseSample.Data());
   if (syst.name().Data() == "Nominal") samples.insert(samples.begin(),cfg.tunfold_ResponseSampleAlt.Data());
   samples.insert(samples.end(),cfg.tunfold_bkgSamples_other.begin(),cfg.tunfold_bkgSamples_other.end());
   
   for (TString sample : samples){     // loop over all MC samples
      bool isSignal = (sample == cfg.tunfold_ResponseSample.Data());
      bool isSignalAlt = (sample == cfg.tunfold_ResponseSampleAlt.Data());
      bool isBKGother = (std::find(cfg.tunfold_bkgSamples_other.begin(),cfg.tunfold_bkgSamples_other.end(), sample) != cfg.tunfold_bkgSamples_other.end());
            
      auto [minTreePath_current, currentSample, sfUnc] = getPath_SampleName_SF(sample,minTreePath_nominal,minTreePath_syst,isSignal,isBKGother,syst);
      
      for (int i=1; i<=cfg.getTotalFileNR(currentSample.Data()); i++){     // loop over each minTree of sample
   
         TFile *signalFile=new TFile(minTreePath_current+currentSample+"_"+std::to_string(i)+".root","read");
         TTree *signalTree=(TTree *) signalFile->Get("ttbar_res100.0/"+currentSample);
            
         if(!signalTree) {
            cout<<"could not read 'signal' tree\n";
         }
        
         signalTree->ResetBranchAddresses();
         if(cfg.tunfold_withPF) {
            signalTree->SetBranchAddress("Phi_rec",&phiRec);
            signalTree->SetBranchAddress("MET",&metRec);
         }
         else if(cfg.tunfold_withDNN) {
            signalTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
            signalTree->SetBranchAddress("DNN_MET_pT",&metRec);
         }
         else{
            signalTree->SetBranchAddress("Phi_recPuppi",&phiRec);
            signalTree->SetBranchAddress("PuppiMET",&metRec);
         }
         signalTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec_DNN);
         signalTree->SetBranchAddress("DNN_MET_pT",&metRec_DNN);
         signalTree->SetBranchAddress("Phi_NuNu",&phiGen);
         signalTree->SetBranchAddress("PtNuNu",&metGen);
         signalTree->SetBranchAddress("pTll",&pTllRec);
         signalTree->SetBranchAddress("pTllGen",&pTllGen);
         signalTree->SetBranchAddress("genDecayMode",&genDecayMode);
         signalTree->SetBranchAddress("N",&mcWeight);
         signalTree->SetBranchAddress("SF",&recoWeight);
         
         std::cout<<"loop over "<<currentSample<<"_"<<i<<std::endl;
         
         int entriesToRun=signalTree->GetEntriesFast()*cfg.processFraction;
         int iEv=0;
         for(Int_t ievent=0;ievent<entriesToRun;ievent++) {
            iEv++;
            if (iEv%(std::max(entriesToRun/10,1))==0){
               io::log*".";
               io::log.flush();
            }
            if(signalTree->GetEntry(ievent)<=0) break;
            
            // change reco weight (used for xsec unc.)
            recoWeight*=sfUnc;
            
            // ~mcWeight=mcWeight*(137191.0/35867.05998);   //scale to full Run2 Lumi
            
            //only bin to bin migration
            // ~if(metRec<0 || genDecayMode>3 || metGen<0) continue;
            
            //ignore acceptance
            // ~if(metRec<0) continue;
            
            //ignore fakes
            // ~if(genDecayMode>3 || metGen<0) continue;
            
            // fill background distributions
            if(isBKGother){
               for(Distribution& dist : distribution_vec) dist.fillMCbkg(mcWeight*recoWeight);
               continue;
            }
            
            // store fakes and tt_other backgrounds for signal fraction
            if (metGen<0 || genDecayMode>3 || (genDecayMode!=3 && metGen<40) || (!isSignal && !isSignalAlt)) {
               if (!isSignalAlt) for(Distribution& dist : distribution_vec) dist.fillMCfakes(mcWeight*recoWeight);
               if (!isSignal) for(Distribution& dist : distribution_vec) dist.fillMCfakesAlt(mcWeight*recoWeight);
               continue;
            }
            
            // fill response matrix
            for(Distribution& dist : distribution_vec) {
               if (isSignal) dist.fillMCGenRec(mcWeight,recoWeight);
               else if (isSignalAlt) dist.fillMCGenRecAlt(mcWeight,recoWeight);
            }
         }
         io::log<<"";
      }
   }
   
   // derive additional quantities (e.g. purity and stability) and store histograms
   for(Distribution& dist : distribution_vec) dist.saveMCHists(saver);
   
}

extern "C"
void run()
{
   //Read systematic from command line
   Systematic::Systematic currentSystematic(cfg.systematic);
   bool isNominal = (currentSystematic.type()==Systematic::nominal);
      
   //=======================================================
   // Step 1: open file to save histograms and binning schemes
   TString save_path = getSavePath();
   io::RootFileSaver saver(TString::Format(!cfg.tunfold_withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",save_path.Data(),currentSystematic.name().Data(),cfg.processFraction*100),save_path,true,true,true);
   
   //=======================================================
   // Step 2: define binnings and distributions
   
   std::vector<Distribution> distribution_vec;
   
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu",
                                          {0,40,80,120,160,230},
                                          {0,0.7,1.4,3.141},
                                          {0,20,40,60,80,100,120,140,160,195,230},
                                          {0,0.35,0.7,1.05,1.4,2.27,3.141}
                                          ));
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu_new",     // Fabians optimized binning for DNN
                                          {0,50,70,100,130,160,200},
                                          {0,0.64,1.2,3.141},
                                          {0,25,50,60,70,85,100,115,130,145,160,180,200,300},
                                          {0,0.32,0.64,0.92,1.2,2.2,3.141}
                                          ));
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu_new40",     // Fabians optimized binning for DNN with 0-40
                                          {0,40,70,100,130,160,200},
                                          {0,0.64,1.2,3.141},
                                          {0,20,40,55,70,85,100,115,130,145,160,180,200,300},
                                          {0,0.32,0.64,0.92,1.2,2.2,3.141}
                                          ));
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu_DNN",
                                          {0,40,80,120,160,230},
                                          {0,0.7,1.4,3.141},
                                          {0,20,40,60,80,100,120,140,160,195,230},
                                          {0,0.35,0.7,1.05,1.4,2.27,3.141}
                                          ));
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu_new_DNN",     // Fabians optimized binning for DNN
                                          {0,50,70,100,130,160,200},
                                          {0,0.64,1.2,3.141},
                                          {0,25,50,60,70,85,100,115,130,145,160,180,200,300},
                                          {0,0.32,0.64,0.92,1.2,2.2,3.141}
                                          ));
   distribution_vec.push_back(Distribution("2D_dPhi_pTnunu_new40_DNN",     // Fabians optimized binning for DNN with 0-40
                                          {0,40,70,100,130,160,200},
                                          {0,0.64,1.2,3.141},
                                          {0,20,40,55,70,85,100,115,130,145,160,180,200,300},
                                          {0,0.32,0.64,0.92,1.2,2.2,3.141}
                                          ));
   distribution_vec.push_back(Distribution("pTnunu",
                                          {0,40,70,110,170,260,370},
                                          {0,20,40,55,70,90,110,140,170,215,260,315,370,435}
                                          ));
   distribution_vec.push_back(Distribution("dPhi",
                                          {0.,0.4,0.8,1.2,1.6,2.,2.4,2.8},
                                          {0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.}
                                          ));
   distribution_vec.push_back(Distribution("pTnunu_DNN",
                                          {0,40,70,110,170,260,370},
                                          {0,20,40,55,70,90,110,140,170,215,260,315,370,435}
                                          ));
   distribution_vec.push_back(Distribution("pTnunu_new_DNN",            // Fabians optimzed binning for DNN
                                          {0,40,65,95,135,195,255,315,380},
                                          {0,20,40,52.5,65,80,95,115,135,165,195,225,255,285,315,347.5,380,440}
                                          ));
   distribution_vec.push_back(Distribution("dPhi_DNN",
                                          {0.,0.4,0.8,1.2,1.6,2.,2.4,2.8},
                                          {0.,0.2,0.4,0.6,0.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,3.}
                                          ));
   distribution_vec.push_back(Distribution("pTll",
                                          {0.,40,80,100,120,150,200,250,300},
                                          {0.,20,40,60,80,90,100,110,120,135,150,175,200,225,250,275,300,350}
                                          ));
   distribution_vec.push_back(Distribution("inclusive",
                                          {0.,3.2},
                                          {0.,3.2}
                                          ));

   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //=======================================================
   // Step 3: book and fill data histograms (only for nominal case)
   if (isNominal) {
      loopDataEvents(distribution_vec,saver);
   }
   
   //=======================================================
   // Step 4: book and fill histogram of migrations
   //         it receives events from both signal MC and background MC (depends on choosen systematic)
   loopMCEvents(distribution_vec,saver,currentSystematic);
}
