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

   auto inputFile = TFile::Open("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0_new.root");
   auto outputFile = TFile::Open("../output/test_TMVA.root", "RECREATE");

   TMVA::Factory factory("TMVAClassification", outputFile, "!V:ROC:!Correlations:!Silent:Color:!DrawProgressBar:AnalysisType=Classification" );
   
   TMVA::DataLoader loader("dataset");
   for (TString var :{"n_Interactions","nJets", "HT","dPhiMETbJet_Puppi","dPhiMETleadJet_Puppi","dPhiMETlead2Jet_Puppi","dPhiLep1bJet","Jet1_pt","Jet2_pt",
                        "ratioMET_sqrtMETunc_Puppi","dPhiLep1Lep2","dPhiLep1Jet1","ratio_pTj1_vecsum_pT_l1_l2_bjet"}){
     loader.AddVariable(var);
   }   
   
   TTree *t_all;

   inputFile->GetObject("ttbar_res100.0/ttbar_res_dilepton_CP5", t_all);
   TCut mycuts, mycutb;
   TCut met6_sig = "PuppiMET>230 && absmetres_PUPPI<40";
   TCut met6_bg = "PuppiMET>230 && absmetres_PUPPI>40";
   loader.SetInputTrees(t_all,met6_sig, met6_bg);
   
   loader.PrepareTrainingAndTestTree(mycuts, mycutb,
                                   // ~"nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V" );
                                   "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V" );
                                    
 
 
   // improved neural network implementation
   TString tanh_v1 ("Layout=TANH|(N+100)*2,LINEAR"); 
   TString tanh_v2 ("Layout=TANH|50,TANH|20,LINEAR");
   TString tanh_v3 ("Layout=TANH|(N+100)*2,TANH|N+100,LINEAR");
   TString tanh_v4 ("Layout=TANH|100,TANH|50,TANH|10,LINEAR");
   TString tanh_v5 ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
   TString tanh_v6 ("Layout=TANH|300,TANH|100,TANH|30,TANH|10,LINEAR");

   TString softsign_v1 ("Layout=SOFTSIGN|(N+100)*2,LINEAR"); 
   TString softsign_v2 ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
   TString softsign_v3 ("Layout=SOFTSIGN|(N+100)*2,SOFTSIGN|N+100,LINEAR");
   TString softsign_v4 ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|10,LINEAR");
   TString softsign_v5 ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
   TString softsign_v6 ("Layout=SOFTSIGN|300,SOFTSIGN|100,SOFTSIGN|30,SOFTSIGN|10,LINEAR");

   TString relu_v1 ("Layout=RELU|(N+100)*2,LINEAR"); 
   TString relu_v2 ("Layout=RELU|50,RELU|20,LINEAR");
   TString relu_v3 ("Layout=RELU|(N+100)*2,RELU|N+100,LINEAR");
   TString relu_v4 ("Layout=RELU|100,RELU|50,RELU|10,LINEAR");
   TString relu_v5 ("Layout=RELU|50,RELU|30,RELU|20,RELU|10,LINEAR");//does not work, too few nodes
   TString relu_v6 ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");

   TString sigmoid_v1 ("Layout=SIGMOID|(N+100)*2,LINEAR"); //best ROC-Integ .666
   TString sigmoid_v2 ("Layout=SIGMOID|50,SIGMOID|20,LINEAR");
   TString sigmoid_v3 ("Layout=SIGMOID|(N+100)*2,SIGMOID|N+100,LINEAR");
   TString sigmoid_v4 ("Layout=SIGMOID|100,SIGMOID|50,SIGMOID|10,LINEAR");
   TString sigmoid_v5 ("Layout=SIGMOID|50,SIGMOID|30,SIGMOID|20,SIGMOID|10,LINEAR");
   TString sigmoid_v6 ("Layout=SIGMOID|300,SIGMOID|100,SIGMOID|30,SIGMOID|10,LINEAR");

   TString symmrelu_v6 ("Layout=SYMMRELU|300,SYMMRELU|100,SYMMRELU|30,SYMMRELU|10,LINEAR");
   TString sigrelu_4l ("Layout=SIGMOID|300,RELU|100,RELU|30,SIGMOID|10,LINEAR");

   //~ TString training0 ("LearningRate=1e-1,Momentum=0.0,Repetitions=1,ConvergenceSteps=300,BatchSize=20,TestRepetitions=15,WeightDecay=0.001,Regularization=NONE,DropConfig=0.0+0.5+0.5+0.5,DropRepetitions=1,Multithreading=True");
   //~ TString training1 ("LearningRate=1e-2,Momentum=0.5,Repetitions=1,ConvergenceSteps=300,BatchSize=30,TestRepetitions=7,WeightDecay=0.001,Regularization=L2,Multithreading=True,DropConfig=0.0+0.1+0.1+0.1,DropRepetitions=1");
   //~ TString training2 ("LearningRate=1e-2,Momentum=0.3,Repetitions=1,ConvergenceSteps=300,BatchSize=40,TestRepetitions=7,WeightDecay=0.0001,Regularization=L2,Multithreading=True");
   //~ TString training3 ("LearningRate=1e-3,Momentum=0.1,Repetitions=1,ConvergenceSteps=200,BatchSize=70,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE,Multithreading=True");

	//modified training 
   TString training0 ("LearningRate=1e-2,Momentum=0.0,Repetitions=1,ConvergenceSteps=300,BatchSize=20,TestRepetitions=15,WeightDecay=0.001,Regularization=NONE,DropConfig=0.0+0.5+0.5+0.5,DropRepetitions=1,Multithreading=True");
   TString training1 ("LearningRate=1e-3,Momentum=0.5,Repetitions=1,ConvergenceSteps=200,BatchSize=30,TestRepetitions=7,WeightDecay=0.0001,Regularization=L2,Multithreading=True,DropConfig=0.0+0.1+0.1+0.1,DropRepetitions=1");
   //~ TString training0 ("LearningRate=1e-3,Momentum=0.0,Repetitions=1,ConvergenceSteps=300,BatchSize=10,TestRepetitions=15,WeightDecay=0.001,Regularization=L2,DropConfig=0.0+0.5+0.5+0.5,DropRepetitions=1,Multithreading=True");
   //~ TString training1 ("LearningRate=1e-4,Momentum=0.5,Repetitions=1,ConvergenceSteps=200,BatchSize=20,TestRepetitions=7,WeightDecay=0.0001,Regularization=L2,Multithreading=True,DropRepetitions=1");


	// ~for (TString var:{tanh_v1, tanh_v2, tanh_v3, tanh_v4, tanh_v5, tanh_v6}){
	//~ for (TString var:{sigmoid_v1, sigmoid_v2, sigmoid_v3, sigmoid_v4, sigmoid_v5, sigmoid_v6}){
   //~ for (TString var:{softsign_v1, softsign_v2, softsign_v3, softsign_v4, softsign_v5, softsign_v6}){
   //~ for (TString var:{relu_v1, relu_v2, relu_v3, relu_v4, relu_v5, relu_v6}){
   for (TString var:{sigmoid_v1}){
      TString trainingStrategyString ("TrainingStrategy=");
   //~ trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;
      trainingStrategyString += training0 + "|" + training1 ;

      //TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CROSSENTROPY");
      TString nnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
      //TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
      nnOptions.Append (":");
    
      nnOptions.Append (var);

      nnOptions.Append (":");
      nnOptions.Append (trainingStrategyString);
      //~ factory.BookMethod(&loader, TMVA::Types::kDNN, "DNN_RELU_v5_t2m2_PUPPI_met6", nnOptions ); // NN
      factory.BookMethod(&loader, TMVA::Types::kDNN,"DNN", nnOptions ); // NN
   }
   factory.TrainAllMethods();
   
   factory.TestAllMethods();
   factory.EvaluateAllMethods();
   
   auto c1 = factory.GetROCCurve(&loader);
   c1->Draw();
   c1->SaveAs("../output/test_dnn_all6tanh_t2_PUPPI.pdf");
}
