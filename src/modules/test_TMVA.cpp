// Script for first test of TMVA in framework
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>
#include <TROOT.h>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"

Config const &cfg=Config::get();

extern "C"
void run()
{
   

   TMVA::Tools::Instance();

   auto inputFile = TFile::Open("https://raw.githubusercontent.com/iml-wg/tmvatutorials/master/inputdata.root");
   auto outputFile = TFile::Open("../output/TMVAOutputCV.root", "RECREATE");

   TMVA::Factory factory("TMVAClassification", outputFile,
                      "!V:ROC:!Correlations:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );
   
   TMVA::DataLoader loader("dataset");
   loader.AddVariable("var1");
   loader.AddVariable("var2");
   loader.AddVariable("var3");
   

   TTree *tsignal, *tbackground;
   inputFile->GetObject("Sig", tsignal);
   inputFile->GetObject("Bkg", tbackground);

   TCut mycuts, mycutb;

   loader.AddSignalTree    (tsignal,     1.0);   //signal weight  = 1
   loader.AddBackgroundTree(tbackground, 1.0);   //background weight = 1 
   loader.PrepareTrainingAndTestTree(mycuts, mycutb,
                                   "nTrain_Signal=1000:nTrain_Background=1000:SplitMode=Random:NormMode=NumEvents:!V" );
                                   

   

   //Boosted Decision Trees
   factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT",
                      "!V:NTrees=200:MinNodeSize=2.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

   //Multi-Layer Perceptron (Neural Network)
   factory.BookMethod(&loader, TMVA::Types::kMLP, "MLP",
                      "!H:!V:NeuronType=tanh:VarTransform=N:NCycles=100:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

   

   factory.TrainAllMethods();
   
   factory.TestAllMethods();
   factory.EvaluateAllMethods();
   
   auto c1 = factory.GetROCCurve(&loader);
   c1->Draw();
   c1->SaveAs("test.pdf");
}


