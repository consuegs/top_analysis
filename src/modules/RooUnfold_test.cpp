//Script to test different unfolding methods from RooUnfold
#include <iostream>
#include <cmath>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include "TUnfoldDensity.h"
#include <TRandom3.h>
#include <TProfile.h>


#include "RooUnfold/RooUnfoldResponse.h"
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
   // unfolded sample
   // ~TString sample="MadGraph";
   TString sample="dilepton";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="dilepton";
   
   // switch on histogram errors
   TH1::SetDefaultSumw2();

   //==============================================
   // step 1 : open output file
   io::RootFileSaver saver(TString::Format("RooUnfold%.1f.root",cfg.processFraction*100),"RooUnfold_results_"+sample+"_"+sample_response);

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format("TUnfold%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;


   // read histograms
   TH1 *histDataReco=histReader.read<TH1>(input_loc+"/histDataReco");
   TH1 *histDataTruth=histReader.read<TH1>(input_loc+"/histDataTruth");
   TH2 *histMCGenRec=histReader.read<TH2>(input_loc+"/histMCGenRec");

   if((!histDataReco)||(!histDataTruth)||(!histMCGenRec)) {
      cout<<"problem to read input histograms\n";
   }

   //========================
   // Step 3: unfolding
   
   RooUnfoldResponse response(2, 3, 4);
   

   //===========================
   // Step 4: retreive and plot unfolding results

}
