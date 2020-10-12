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
#include "TUnfoldBinning.h"

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;

Config const &cfg=Config::get();

extern "C"
void run()
{
   std::cout<<"---------------------------------------"<<std::endl;
   std::cout<<"Still have to take weights into account (especially how to handle SF weight, which gets important using real data)"<<std::endl;
   std::cout<<"---------------------------------------"<<std::endl;
   
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
   
   // Use deep instead of pfMET
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
   
   // number of met and phi bins and binning
   // ~int NBIN_MET_FINE=6;
   // ~int NBIN_MET_FINE=8;
   int NBIN_MET_FINE=10;
   // ~int NBIN_MET_FINE=12;
   int NBIN_PHI_FINE=6;
   // ~std::vector<double> metBinsFine_vector={0,20,40,80,120,175,230};
   // ~std::vector<double> metBinsFine_vector={0,20,40,60,80,100,120,175,230};
   std::vector<double> metBinsFine_vector={0,20,40,60,80,100,120,140,160,195,230};
   std::vector<double> phiBinsFine_vector={0,0.35,0.7,1.05,1.4,2.27,3.141};
   // ~std::vector<double> phiBinsFine_vector={0,0.35,0.7,1.05,1.4,1.8,3.141};
   if(withSameBins){
      // ~NBIN_MET_FINE=3;
      NBIN_MET_FINE=3;
      NBIN_PHI_FINE=3;
      metBinsFine_vector.resize(NBIN_MET_FINE+1);
      phiBinsFine_vector.resize(NBIN_PHI_FINE+1);
      metBinsFine_vector={0,40,120,230};
      phiBinsFine_vector={0,0.7,1.4,3.141};
   }
   Double_t* metBinsFine=&metBinsFine_vector[0];
   Double_t* phiBinsFine=&phiBinsFine_vector[0];

   // ~int NBIN_MET_COARSE=3;
   // ~int NBIN_MET_COARSE=4;
   int NBIN_MET_COARSE=5;
   // ~int NBIN_MET_COARSE=6;
   int NBIN_PHI_COARSE=3;
   // ~Double_t metBinsCoarse[NBIN_MET_COARSE+1]={0,40,120,230};
   // ~Double_t metBinsCoarse[NBIN_MET_COARSE+1]={0,40,80,120,230};
   Double_t metBinsCoarse[NBIN_MET_COARSE+1]={0,40,80,120,160,230};
   Double_t phiBinsCoarse[NBIN_PHI_COARSE+1]={0,0.7,1.4,3.141};
   
   //=======================================================================
   // detector level binning scheme

   TUnfoldBinning *detectorBinning=new TUnfoldBinning("detector");
   TUnfoldBinning *detectorDistribution=detectorBinning->AddBinning("detectordistribution");
   detectorDistribution->AddAxis("met",NBIN_MET_FINE,metBinsFine,
                                 false, // no underflow bin (not reconstructed)
                                 true // overflow bin
                                 // ~false // overflow bin
                                 );
   detectorDistribution->AddAxis("phi",NBIN_PHI_FINE,phiBinsFine,
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
                           // ~false // overflow bin
                           );
   signalBinning->AddAxis("phigen",NBIN_PHI_COARSE,phiBinsCoarse,
                           false, // underflow bin
                           false // overflow bin
                           );
  
  

   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //=======================================================
   // Step 1: open file to save histograms and binning schemes
   TString save_path = "TUnfold_binning_"+sample+"_"+sample_response;
   if (withBSM) save_path+="_BSM";
   if (withPuppi) save_path+="_Puppi";
   if (withDeep) save_path+="_Deep";
   if (withSameBins) save_path+="_SameBins";
   if (withPTreweight) {
      save_path+="_PTreweight"+scale;
      sample+="_PTreweight"+scale;
   }
   
   io::RootFileSaver saver(TString::Format(!withScaleFactor ? "TUnfold%.1f.root" : "TUnfold_SF91_%.1f.root",cfg.processFraction*100),save_path);

   //=======================================================
   // Step 2: save binning to output file

   saver.save(*detectorBinning,"detector_binning");
   saver.save(*generatorBinning,"generator_binning");

   //=======================================================
   // Step 3: book and fill data histograms

   Float_t phiRec,metRec,phiGen,metGen,mcWeight,recoWeight,genMET,reweight;
   UInt_t genDecayMode;
   // ~Int_t istriggered,issignal;

   TH1 *histDataReco=detectorBinning->CreateHistogram("histDataReco");
   TH1 *histDataReco_coarse=generatorBinning->CreateHistogram("histDataReco_coarse");
   TH1 *histDataTruth=generatorBinning->CreateHistogram("histDataTruth");
   TH1 *histDataTruth_fakes=generatorBinning->CreateHistogram("histDataTruth_fakes");

   TString inputFile = TString::Format("ttbar_res%.1f.root",cfg.processFraction*100);
   // ~TString inputFile = TString::Format("ttbar_res%.1f_reweightAlsoFakes.root",cfg.processFraction*100);
   // ~TString inputFile = TString::Format("ttbar_res%.1f_extremeReweight.root",cfg.processFraction*100);
   TString inputString = TString::Format("ttbar_res%.1f",cfg.processFraction*100);
   if (withScaleFactor) inputFile = "ttbar_res_SF0.910000100.0.root";
   TFile *dataFile=new TFile("/net/data_cms1b/user/dmeuser/top_analysis/output/"+inputFile);
   TTree *dataTree=(TTree *) dataFile->Get(sample!=""? inputString+"/ttbar_res_"+sample: inputString+"/ttbar_res");
   // ~TFile *dataFile_alt=new TFile("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0_new.root");
   // ~TTree *dataTree=(TTree *) dataFile_alt->Get(sample!=""? inputString+"/ttbar_res_"+sample: inputString+"/ttbar_res");
   TTree *BSMTree=(TTree *) dataFile->Get(inputString+"/ttbar_res_T2tt_650_350");
   

   if(!dataTree) {
      cout<<"could not read 'data' tree\n";
   }
   if(!BSMTree) {
      cout<<"could not read 'BSM' tree\n";
   }
  
   dataTree->ResetBranchAddresses();
   dataTree->SetBranchAddress("Phi_rec",&phiRec);
   dataTree->SetBranchAddress("MET",&metRec);
   // ~dataTree->SetBranchAddress("genMET",&metRec);
   // ~dataTree->SetBranchAddress("genMET",&genMET);
   if(withPuppi) {
      dataTree->SetBranchAddress("Phi_recPuppi",&phiRec);
      dataTree->SetBranchAddress("PuppiMET",&metRec);
   }
   if(withDeep) {
      dataTree->SetBranchAddress("Phi_recDeep",&phiRec);
      dataTree->SetBranchAddress("DeepMET",&metRec);
   }
   dataTree->SetBranchAddress("Phi_NuNu",&phiGen);
   dataTree->SetBranchAddress("PtNuNu",&metGen);
   dataTree->SetBranchAddress("genDecayMode",&genDecayMode);
   dataTree->SetBranchAddress("N",&mcWeight);
   dataTree->SetBranchAddress("SF",&recoWeight);
   // ~dataTree->SetBranchAddress("issignal",&issignal);
   if(withPTreweight) dataTree->SetBranchAddress("reweight_PTnunu",&reweight);
   else reweight=1.0;

   cout<<"loop over data events\n";
   

   for(Int_t ievent=0;ievent<dataTree->GetEntriesFast();ievent++) {
   // ~for(Int_t ievent=0;ievent<dataTree->GetEntriesFast()/2;ievent++) {
      if(dataTree->GetEntry(ievent)<=0) break;

      //only bin to bin migration
      if(metRec<0 || genDecayMode>3 || metGen<0) continue;
      
      //ignore acceptance
      // ~if(metRec<0) continue;
      
      //ignore fakes
      // ~if(genDecayMode>3 || metGen<0) continue;
      
      // ~if (metGen>0) mcWeight*=sqrt(sqrt(sqrt(sqrt(metGen))));
      
      //remove tau events, which are not reconstructed (no gen match, no reco match, should be included into distributions.cpp!!!!)
      if(metRec<0 && genDecayMode>3) continue;

      // fill histogram with data truth parameters
      Int_t genbinNumber=signalBinning->GetGlobalBinNumber(metGen,phiGen);
      
      //reweight
      mcWeight*=reweight;
      
      if (metGen<0 || genDecayMode>3) {
         histDataTruth_fakes->Fill(genbinNumber,mcWeight);
      }
      else histDataTruth->Fill(genbinNumber,mcWeight);

      // fill histogram with reconstructed quantities
      if (metRec<0) continue;   //events that are not reconstructed
      Int_t binNumber=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
      histDataReco->Fill(binNumber,mcWeight);
      
      Int_t binNumber_coarse=signalBinning->GetGlobalBinNumber(metRec,phiRec);
      histDataReco_coarse->Fill(binNumber_coarse,mcWeight);
   }
   
   if (withBSM) {
      BSMTree->ResetBranchAddresses();
      BSMTree->SetBranchAddress("Phi_rec",&phiRec);
      BSMTree->SetBranchAddress("MET",&metRec);
      if(withPuppi) {
         BSMTree->SetBranchAddress("Phi_recPuppi",&phiRec);
         BSMTree->SetBranchAddress("PuppiMET",&metRec);
      }
      if(withDeep) {
         BSMTree->SetBranchAddress("Phi_recDeep",&phiRec);
         BSMTree->SetBranchAddress("DeepMET",&metRec);
      }
      BSMTree->SetBranchAddress("Phi_NuNu",&phiGen);
      BSMTree->SetBranchAddress("PtNuNu",&metGen);
      BSMTree->SetBranchAddress("genDecayMode",&genDecayMode);
      BSMTree->SetBranchAddress("N",&mcWeight);
      BSMTree->SetBranchAddress("SF",&recoWeight);
      
      for(Int_t ievent=0;ievent<BSMTree->GetEntriesFast();ievent++) {
         if(BSMTree->GetEntry(ievent)<=0) break;

         // fill histogram with reconstructed quantities
         if (metRec<0) continue;   //events that are not reconstructed
         Int_t binNumber=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
         histDataReco->Fill(binNumber,mcWeight);
      }
   }
   
   // set reco bin error to error expected in data
   for(int i=0;i<=histDataReco->GetNbinsX()+1;i++) {
      histDataReco->SetBinError(i,sqrt(histDataReco->GetBinContent(i)));
   }
   for(int i=0;i<=histDataReco_coarse->GetNbinsX()+1;i++) {
      histDataReco_coarse->SetBinError(i,sqrt(histDataReco_coarse->GetBinContent(i)));
   }
   
   saver.save(*histDataReco,"histDataReco");
   saver.save(*histDataReco_coarse,"histDataReco_coarse");
   saver.save(*histDataTruth_fakes,"histDataTruth_fakes");
   saver.save(*histDataTruth,"histDataTruth");

   delete dataTree;
   delete BSMTree;
   //=======================================================
   // Step 4: book and fill histogram of migrations
   //         it receives events from both signal MC and background MC (right now only signal MC)

   TH2 *histMCGenRec=TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning,detectorBinning,"histMCGenRec");
   TH2 *histMCGenRec_purity=TUnfoldBinning::CreateHistogramOfMigrations(generatorBinning,generatorBinning,"histMCGenRec_purity");
   TH1 *histMCRec_fakes=detectorBinning->CreateHistogram("histMCRec_fakes");
   TH1 *histMCRec_fakes_coarse=generatorBinning->CreateHistogram("histMCRec_fakes_coarse");

   TTree *signalTree=(TTree *) dataFile->Get(sample_response!=""? inputString+"/ttbar_res_"+sample_response: inputString+"/ttbar_res");
   
   if(!signalTree) {
      cout<<"could not read 'signal' tree\n";
   }
  
   signalTree->ResetBranchAddresses();
   signalTree->SetBranchAddress("Phi_rec",&phiRec);
   signalTree->SetBranchAddress("MET",&metRec);
   // ~signalTree->SetBranchAddress("genMET",&metRec);
   // ~signalTree->SetBranchAddress("genMET",&genMET);
   if(withPuppi) {
      signalTree->SetBranchAddress("Phi_recPuppi",&phiRec);
      signalTree->SetBranchAddress("PuppiMET",&metRec);
   }
   if(withDeep) {
      signalTree->SetBranchAddress("Phi_recDeep",&phiRec);
      signalTree->SetBranchAddress("DeepMET",&metRec);
   }
   // ~signalTree->SetBranchAddress("istriggered",&istriggered);
   signalTree->SetBranchAddress("Phi_NuNu",&phiGen);
   signalTree->SetBranchAddress("PtNuNu",&metGen);
   signalTree->SetBranchAddress("genDecayMode",&genDecayMode);
   signalTree->SetBranchAddress("N",&mcWeight);
   signalTree->SetBranchAddress("SF",&recoWeight);

   cout<<"loop over MC signal events\n";
   
  
   for(Int_t ievent=0;ievent<signalTree->GetEntriesFast();ievent++) {
   // ~for(Int_t ievent=signalTree->GetEntriesFast()/2;ievent<signalTree->GetEntriesFast();ievent++) {
      if(signalTree->GetEntry(ievent)<=0) break;
      
      //only bin to bin migration
      if(metRec<0 || genDecayMode>3 || metGen<0) continue;
      
      //ignore acceptance
      // ~if(metRec<0) continue;
      
      //ignore fakes
      // ~if(genDecayMode>3 || metGen<0) continue;

      //remove tau events, which are not reconstructed (no gen match, no reco match, should be included into distributions.cpp!!!!)
      if(metRec<0 && genDecayMode>3) continue;

      // bin number on generator level for signal
      Int_t genBin=signalBinning->GetGlobalBinNumber(metGen,phiGen);

      // bin number on reconstructed level
      // bin number 0 corresponds to non-reconstructed events
      Int_t recBin=0;
      Int_t recBin_purity=0;
      recBin=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
      recBin_purity=signalBinning->GetGlobalBinNumber(metRec,phiRec);
      
      // ~if (metGen<0 || genDecayMode>3) continue;
      if (metGen<0 || genDecayMode>3) {
         histMCRec_fakes->Fill(recBin,mcWeight);
         histMCRec_fakes_coarse->Fill(recBin_purity,mcWeight);
         continue;
      }

      histMCGenRec->Fill(genBin,recBin,mcWeight);
      histMCGenRec_purity->Fill(genBin,recBin_purity,mcWeight);
      
      // ~if(genBin==1 && recBin==0) std::cout<<metGen<<"   "<<phiGen<<"   "<<metRec<<"   "<<phiRec<<std::endl;

      /* Still have to fill with weights
      hist_migration_MC[k]->Fill(genBin,recoBin,weightRec);
      // count fraction of events which have a different weight on truth and reco (in reco underflow bin)
      // this is required for TUnfold to function properly
      hist_migration_MC[k]->Fill(genBin,0.,weightGen-weightRec);
      */
   }
  
   saver.save(*histMCGenRec,"histMCGenRec");
   saver.save(*histMCGenRec_purity,"histMCGenRec_sameBins");
   saver.save(*histMCRec_fakes,"histMCRec_fakes");
   saver.save(*(histMCGenRec->ProjectionX()),"histMCGenRec_projX");
   saver.save(*(histMCGenRec->ProjectionY()),"histMCGenRec_projY");
   
   //Calculate Signal fration
   TH1 *hist_SignalFraction=detectorBinning->CreateHistogram("hist_SignalFraction");
   TH1 *histMCRec=histMCGenRec->ProjectionY("",1,-1);
   saver.save(*histMCRec,"histMCRec");
   hist_SignalFraction->Add(histMCRec);
   histMCRec->Add(histMCRec_fakes);
   hist_SignalFraction->Divide(histMCRec);
   
   for(int i=1; i<=hist_SignalFraction->GetNbinsX(); i++){  //Is this correct or not????
      hist_SignalFraction->SetBinError(i,0);
   }
   saver.save(*hist_SignalFraction,"hist_SignalFraction");
   
   //Calculate Signal fration for coarse binning
   TH1 *hist_SignalFraction_coarse=generatorBinning->CreateHistogram("hist_SignalFraction_coarse");
   TH1 *histMCRec_coarse=histMCGenRec_purity->ProjectionY("",1,-1);
   saver.save(*histMCRec_coarse,"histMCRec_coarse");
   hist_SignalFraction_coarse->Add(histMCRec_coarse);
   histMCRec_coarse->Add(histMCRec_fakes_coarse);
   hist_SignalFraction_coarse->Divide(histMCRec_coarse);
   
   for(int i=1; i<=hist_SignalFraction_coarse->GetNbinsX(); i++){  //Is this correct or not????
      hist_SignalFraction_coarse->SetBinError(i,0);
   }
   saver.save(*hist_SignalFraction_coarse,"hist_SignalFraction_coarse");
   
   //Save normalized reco distribution for fake and non fake events
   TH1D histGen_fakes=*(histMCGenRec_purity->ProjectionY("Fakes",1,1));
   TH1D histGen_noFakes=*(histMCGenRec_purity->ProjectionY("Fakes",2,-1));
   // ~histRec_fakes.Scale(1/histRec_fakes.Integral());
   // ~histRec_noFakes.Scale(1/histRec_noFakes.Integral());
   saver.save(histGen_fakes,"histGen_fakes");
   saver.save(histGen_noFakes,"histGen_noFakes");

   delete signalTree;

   /*
   TFile *bgrFile=new TFile("testUnfold5_background.root");
   TTree *bgrTree=(TTree *) bgrFile->Get("background");

   if(!bgrTree) {
     cout<<"could not read 'background' tree\n";
   }

   bgrTree->ResetBranchAddresses();
   bgrTree->SetBranchAddress("phirec",&phiRec);
   bgrTree->SetBranchAddress("metrec",&metRec);
   bgrTree->SetBranchAddress("discr",&discr);
   bgrTree->SetBranchAddress("istriggered",&istriggered);
   bgrTree->SetBranchStatus("*",1);

   cout<<"loop over MC background events\n";

   for(Int_t ievent=0;ievent<bgrTree->GetEntriesFast();ievent++) {
      if(bgrTree->GetEntry(ievent)<=0) break;

      // here, for background only reconstructed quantities are known
      // and only the reconstructed events are relevant
      if(istriggered) {
         // bin number on generator level for background
         Int_t genBin=bgrBinning->GetGlobalBinNumber(metRec,phiRec);
         // bin number on reconstructed level
         Int_t recBin=detectordistribution->GetGlobalBinNumber
            (metRec,phiRec,discr);
         histMCGenRec->Fill(genBin,recBin);
      }
   }

   delete bgrTree;
   delete bgrFile;
   */
   //efficiency
   TH1 *hist_efficiency=generatorBinning->CreateHistogram("efficiency");
   for(int binGen=0;binGen<=histMCGenRec->GetNbinsX()+1;binGen++) {
      double sum0=0.;
      double sum1=0.;
      for(int binRec=0;binRec<=histMCGenRec->GetNbinsY()+1;binRec++) {
         double c=histMCGenRec->GetBinContent(binGen,binRec);
         sum0+=c;
         // ~if((binRec>0)&&(binRec<=histMCGenRec->GetNbinsY())) sum1+=c;
         if((binRec>0)&&(binRec<=histMCGenRec->GetNbinsY()+1)) sum1+=c;
      } 
      if(sum0>0.0) hist_efficiency->SetBinContent(binGen,sum1/sum0);
   }
   saver.save(*hist_efficiency,"hist_efficiency");
   
   // calculate bin purities and stabilities
   TH1 *hist_purity=generatorBinning->CreateHistogram("purity");
   TH1 *hist_stability=generatorBinning->CreateHistogram("stability");
   TH1 *hist_N_GenRec=generatorBinning->CreateHistogram("N_GenRec");
   TH1 *hist_N_Rec=generatorBinning->CreateHistogram("N_Rec");
   for(int binRec=0;binRec<=hist_purity->GetNbinsX()+1;binRec++) {
      double sum=0.;
      // ~for(int binGen=0;binGen<=hist_purity->GetNbinsX()+1;binGen++) {
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
   
   saver.save(*hist_purity,"hist_purity");
   saver.save(*hist_N_GenRec,"hist_N_GenRec");
   saver.save(*hist_N_Rec,"hist_N_Rec");
   
   for(int binGen=0;binGen<=hist_purity->GetNbinsX()+1;binGen++) {
      double sum=0.;
      // ~for(int binGen=0;binGen<=hist_purity->GetNbinsX()+1;binGen++) {
      for(int binRec=1;binRec<=hist_purity->GetNbinsX()+1;binRec++) {
         sum += histMCGenRec_purity->GetBinContent(binGen,binRec);
      }
      double p=0.;
      if(sum>0.0) {
         p=histMCGenRec_purity->GetBinContent(binGen,binGen)/sum;
      }
      hist_stability->SetBinContent(binGen,p);
   }
   
   saver.save(*hist_stability,"hist_stability");
   
   // ~// Bin by bin correction factors
   // ~TH1* correction = histMCGenRec_purity->ProjectionX();
   // ~correction->Divide(histMCGenRec_purity->ProjectionY());
   
   // ~TH1* result_binBYbin = (TH1*)histDataReco_coarse->Clone("");
   // ~result_binBYbin->Multiply(correction);
   
   // ~for (int i=1;i<=result_binBYbin->GetNbinsX();i++){
      // ~std::cout<<result_binBYbin->GetBinContent(i)<<"   "<<histDataTruth->GetBinContent(i)<<"   "<<correction->GetBinContent(i)<<std::endl;
   // ~}

}
