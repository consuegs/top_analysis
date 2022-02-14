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

TString getSavePath(){
   TString save_path = "TUnfold_binning_"+cfg.tunfold_InputSamples[0]+"_"+cfg.tunfold_ResponseSample;
   if (cfg.tunfold_withBSM) save_path+="_BSM";
   if (!cfg.tunfold_withDNN) save_path+="_Puppi";
   if (cfg.tunfold_withDNN) save_path+="_DNN";
   if (cfg.tunfold_withSameBins) save_path+="_SameBins";
   if (cfg.tunfold_withPTreweight) {
      save_path+="_PTreweight"+cfg.tunfold_scalePTreweight;
   }
   return save_path;
}

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

   TUnfoldBinning *detectorBinning=new TUnfoldBinning("detector");
   TUnfoldBinning *detectorDistribution=detectorBinning->AddBinning("detectordistribution");
   detectorDistribution->AddAxis("met",NBIN_MET_FINE,metBinsFine,
                                 false, // no underflow bin (not reconstructed)
                                 true // overflow bin
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
                           );
   signalBinning->AddAxis("phigen",NBIN_PHI_COARSE,phiBinsCoarse,
                           false, // underflow bin
                           false // overflow bin
                           );
   
   return {detectorDistribution,signalBinning};
}

void fillDataHists(TUnfoldBinning* const &detectorDistribution, TUnfoldBinning* const &signalBinning, io::RootFileSaver const &saver){
   
   TString sample = cfg.tunfold_InputSamples[0];
   
   if (cfg.tunfold_withPTreweight) {
      sample+="_PTreweight"+cfg.tunfold_scalePTreweight;
   }
   
   // define variables for reading tree
   Float_t phiRec,metRec,phiGen,metGen,mcWeight,recoWeight,genMET,reweight;
   UInt_t genDecayMode;
   
   // define histograms to be filles
   TH1 *histDataReco=detectorDistribution->CreateHistogram("histDataReco");
   TH1 *histDataReco_coarse=signalBinning->CreateHistogram("histDataReco_coarse");
   TH1 *histDataTruth=signalBinning->CreateHistogram("histDataTruth");
   TH1 *histDataTruth_fakes=signalBinning->CreateHistogram("histDataTruth_fakes");
   
   TString minTreePath_Nominal = TString::Format("%s/100.0/Nominal/",cfg.minTreePath.Data(),cfg.treeVersion.Data());
   
   for (TString currentSample : cfg.tunfold_InputSamples){     // loop over samples (data or pseudo data)
      bool isSignal = (currentSample == "TTbar_diLepton");
      for (int i=1; i<=cfg.getTotalFileNR(currentSample.Data()); i++){     // loop over each minTree of sample
         TFile *dataFile=new TFile(minTreePath_Nominal+currentSample+"_"+std::to_string(i)+".root","read");
         TTree *dataTree=(TTree *) dataFile->Get("ttbar_res100.0/"+currentSample);

         if(!dataTree) {
            cout<<"could not read 'data' tree\n";
         }
         
         // set tree variables
         dataTree->ResetBranchAddresses();
         dataTree->SetBranchAddress("Phi_recPuppi",&phiRec);
         dataTree->SetBranchAddress("PuppiMET",&metRec);
         if(cfg.tunfold_withDNN) {
            dataTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
            dataTree->SetBranchAddress("DNN_MET_pT",&metRec);
         }
         dataTree->SetBranchAddress("Phi_NuNu",&phiGen);
         dataTree->SetBranchAddress("PtNuNu",&metGen);
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

            // fill histogram with data truth parameters
            Int_t genbinNumber=signalBinning->GetGlobalBinNumber(metGen,phiGen);
            
            //reweight
            mcWeight*=reweight;
            
            if (isSignal){
               if (metGen<0 || genDecayMode>3 || (genDecayMode!=3 && metGen<40)) {
                  histDataTruth_fakes->Fill(genbinNumber,mcWeight);
               }
               else histDataTruth->Fill(genbinNumber,mcWeight);
            }

            // fill histogram with reconstructed quantities
            if (metRec<0) continue;   //events that are not reconstructed
            
            Int_t binNumber=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
            histDataReco->Fill(binNumber,mcWeight*recoWeight);
            
            Int_t binNumber_coarse=signalBinning->GetGlobalBinNumber(metRec,phiRec);
            histDataReco_coarse->Fill(binNumber_coarse,mcWeight*recoWeight);
         }
         io::log<<"";
      }
   }
   
   TFile *bsmFile=new TFile(minTreePath_Nominal+"T2tt_650_350.root","read");
   TTree *BSMTree=(TTree *) bsmFile->Get("ttbar_res100.0/T2tt_650_350");
   if(!BSMTree) {
      cout<<"could not read 'BSM' tree\n";
   }
   
   if (cfg.tunfold_withBSM) {
      BSMTree->ResetBranchAddresses();
      BSMTree->SetBranchAddress("Phi_recPuppi",&phiRec);
      BSMTree->SetBranchAddress("PuppiMET",&metRec);
      if(cfg.tunfold_withDNN) {
         BSMTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
         BSMTree->SetBranchAddress("DNN_MET_pT",&metRec);
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
         Int_t binNumber=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
         histDataReco->Fill(binNumber,mcWeight);
         
         Int_t binNumber_coarse=signalBinning->GetGlobalBinNumber(metRec,phiRec);
         histDataReco_coarse->Fill(binNumber_coarse,mcWeight);
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

   delete BSMTree;
}

TH1* deriveSignalFraction(TH2* const histMCGenRec, TH1* const histMCRec_fakes, TUnfoldBinning* const &detectorDistribution, io::RootFileSaver const &saver,
                           TString const &histName){
   TH1 *hist_SignalFraction=detectorDistribution->CreateHistogram(histName);
   TH1 *histMCRec=histMCGenRec->ProjectionY("",1,-1);
   saver.save(*histMCRec,"histMCRec");
   
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

void derivePurityStability(TH2* const histMCGenRec_purity, TUnfoldBinning* const &signalBinning, io::RootFileSaver const &saver){
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
   
   saver.save(*hist_purity,"hist_purity");
   saver.save(*hist_N_GenRec,"hist_N_GenRec");
   saver.save(*hist_N_Rec,"hist_N_Rec");
   
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
   
   saver.save(*hist_stability,"hist_stability");
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
   else if(syst.type() == Systematic::bFrag || syst.type() == Systematic::bSemilep || syst.type() == Systematic::meFacScale || syst.type() == Systematic::meRenScale){
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

void fillMCHists(TUnfoldBinning* const &detectorDistribution, TUnfoldBinning* const &signalBinning, io::RootFileSaver const &saver, Systematic::Systematic const &syst){
   TString minTreePath_nominal = TString::Format("%s/100.0/%s/",cfg.minTreePath.Data(),"Nominal");
   TString minTreePath_syst = TString::Format("%s/100.0/%s/",cfg.minTreePath.Data(),syst.name().Data());

   TH2 *histMCGenRec=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning,detectorDistribution,"histMCGenRec");
   TH2 *histMCGenRec_purity=TUnfoldBinning::CreateHistogramOfMigrations(signalBinning,signalBinning,"histMCGenRec_purity");
   TH1 *histMCRec_fakes=detectorDistribution->CreateHistogram("histMCRec_fakes");
   TH1 *histMCRec_fakes_coarse=signalBinning->CreateHistogram("histMCRec_fakes_coarse");
   TH1 *histMCRec_bkg=detectorDistribution->CreateHistogram("histMCRec_bkg");
   TH1 *histMCRec_bkg_coarse=signalBinning->CreateHistogram("histMCRec_bkg_coarse");
   
   std::vector<std::string> samples = cfg.tunfold_bkgSamples_ttbar;
   samples.insert(samples.begin(),cfg.tunfold_ResponseSample.Data());
   samples.insert(samples.end(),cfg.tunfold_bkgSamples_other.begin(),cfg.tunfold_bkgSamples_other.end());
   
   for (TString sample : samples){     // loop over all MC samples
      bool isSignal = (sample == "TTbar_diLepton");
      bool isBKGother = (std::find(cfg.tunfold_bkgSamples_other.begin(),cfg.tunfold_bkgSamples_other.end(), sample) != cfg.tunfold_bkgSamples_other.end());
      
      auto [minTreePath_current, currentSample, sfUnc] = getPath_SampleName_SF(sample,minTreePath_nominal,minTreePath_syst,isSignal,isBKGother,syst);
      
      for (int i=1; i<=cfg.getTotalFileNR(currentSample.Data()); i++){     // loop over each minTree of sample
   
         TFile *signalFile=new TFile(minTreePath_current+currentSample+"_"+std::to_string(i)+".root","read");
         TTree *signalTree=(TTree *) signalFile->Get("ttbar_res100.0/"+currentSample);
            
         if(!signalTree) {
            cout<<"could not read 'signal' tree\n";
         }
         
         Float_t phiRec,metRec,phiGen,metGen,mcWeight,recoWeight,genMET,reweight;
         UInt_t genDecayMode;
        
         signalTree->ResetBranchAddresses();
         signalTree->SetBranchAddress("Phi_recPuppi",&phiRec);
         signalTree->SetBranchAddress("PuppiMET",&metRec);
         if(cfg.tunfold_withDNN) {
            signalTree->SetBranchAddress("DNN_MET_dPhi_nextLep",&phiRec);
            signalTree->SetBranchAddress("DNN_MET_pT",&metRec);
         }
         signalTree->SetBranchAddress("Phi_NuNu",&phiGen);
         signalTree->SetBranchAddress("PtNuNu",&metGen);
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
            
            // change reco weight (used for lumi and xsec unc.)
            recoWeight*=sfUnc;
            
            // ~mcWeight=mcWeight*(137191.0/35867.05998);   //scale to full Run2 Lumi
            
            //only bin to bin migration
            // ~if(metRec<0 || genDecayMode>3 || metGen<0) continue;
            
            //ignore acceptance
            // ~if(metRec<0) continue;
            
            //ignore fakes
            // ~if(genDecayMode>3 || metGen<0) continue;

            // bin number on reconstructed level
            // bin number 0 corresponds to non-reconstructed events
            Int_t recBin=0;
            Int_t recBin_purity=0;
            recBin=detectorDistribution->GetGlobalBinNumber(metRec,phiRec);
            recBin_purity=signalBinning->GetGlobalBinNumber(metRec,phiRec);
            
            if(isBKGother){
               histMCRec_bkg->Fill(recBin,mcWeight*recoWeight);
               histMCRec_bkg_coarse->Fill(recBin_purity,mcWeight*recoWeight);
               continue;
            }
            
            if (metGen<0 || genDecayMode>3 || (genDecayMode!=3 && metGen<40) || !isSignal) {    // store fakes and tt_other backgrounds for signal fraction
               histMCRec_fakes->Fill(recBin,mcWeight*recoWeight);
               histMCRec_fakes_coarse->Fill(recBin_purity,mcWeight*recoWeight);
               continue;
            }
                        
            // bin number on generator level for signal
            Int_t genBin=signalBinning->GetGlobalBinNumber(metGen,phiGen);

            histMCGenRec->Fill(genBin,recBin,mcWeight*recoWeight);
            histMCGenRec_purity->Fill(genBin,recBin_purity,mcWeight*recoWeight);
            
            // count fraction of events which have a different weight on truth and reco (in reco underflow bin)
            // this is required for TUnfold to function properly
            histMCGenRec->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
            histMCGenRec_purity->Fill(genBin,0.,mcWeight-(mcWeight*recoWeight));
         }
         io::log<<"";
      }
   }
  
   saver.save(*histMCGenRec,"histMCGenRec");
   saver.save(*histMCGenRec_purity,"histMCGenRec_sameBins");
   saver.save(*histMCRec_fakes,"histMCRec_fakes");
   saver.save(*histMCRec_bkg,"histMCRec_bkg");
   saver.save(*histMCRec_bkg_coarse,"histMCRec_bkg_coarse");
   saver.save(*(histMCGenRec->ProjectionX()),"histMCGenRec_projX");
   saver.save(*(histMCGenRec->ProjectionY()),"histMCGenRec_projY");
   
   //Calculate Signal fration
   TH1* hist_SignalFraction = deriveSignalFraction(histMCGenRec,histMCRec_fakes,detectorDistribution,saver,"hist_SignalFraction");
   saver.save(*hist_SignalFraction,"hist_SignalFraction");
   
   //Calculate Signal fration for coarse binning
   TH1* hist_SignalFraction_coarse = deriveSignalFraction(histMCGenRec_purity,histMCRec_fakes_coarse,signalBinning,saver,"hist_SignalFraction_coarse");
   saver.save(*hist_SignalFraction_coarse,"hist_SignalFraction_coarse");
   
   //Save normalized reco distribution for fake and non fake events
   TH1D histGen_fakes=*(histMCGenRec_purity->ProjectionY("Fakes",1,1));
   TH1D histGen_noFakes=*(histMCGenRec_purity->ProjectionY("Fakes",2,-1));
   saver.save(histGen_fakes,"histGen_fakes");
   saver.save(histGen_noFakes,"histGen_noFakes");
   
   //efficiency
   TH1* hist_efficiency = deriveEfficiency(histMCGenRec,signalBinning);
   saver.save(*hist_efficiency,"hist_efficiency");
   
   // calculate bin purities and stabilities
   derivePurityStability(histMCGenRec_purity,signalBinning,saver);
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
   io::RootFileSaver saver(TString::Format(!cfg.tunfold_withScaleFactor ? "TUnfold/TUnfold_%s_%.1f.root" : "TUnfold/TUnfold_SF91_%s_%.1f.root",currentSystematic.name().Data(),cfg.processFraction*100),save_path);
   
   //=======================================================
   // Step 2: define and save binnings
   std::vector<double> metBinsCoarse_vector = {0,40,80,120,160,230};
   std::vector<double> phiBinsCoarse_vector = {0,0.7,1.4,3.141};
   
   std::vector<double> metBinsFine_vector={0,20,40,60,80,100,120,140,160,195,230};
   std::vector<double> phiBinsFine_vector={0,0.35,0.7,1.05,1.4,2.27,3.141};
   
   // get TUnfoldBinnings for detector and generator level
   auto [detectorDistribution, signalBinning] = defineBinnings2D(metBinsCoarse_vector,phiBinsCoarse_vector,metBinsFine_vector,phiBinsFine_vector);
   
   saver.save(*detectorDistribution,"detector_binning");
   saver.save(*signalBinning,"generator_binning");

   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //=======================================================
   // Step 3: book and fill data histograms (only for nominal case)
   if (isNominal) {
      fillDataHists(detectorDistribution,signalBinning,saver);
   }
   
   //=======================================================
   // Step 4: book and fill histogram of migrations
   //         it receives events from both signal MC and background MC (depends on choosen systematic)
   fillMCHists(detectorDistribution,signalBinning,saver,currentSystematic);
}
