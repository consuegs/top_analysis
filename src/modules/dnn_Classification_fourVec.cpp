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
#include <TSystem.h>

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

Config const &cfg=Config::get();

extern "C"
void run()
{
    TMVA::Tools::Instance();

    TFile *input = TFile::Open("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0_new.root");
    TMVA::DataLoader *loader=new TMVA::DataLoader("dataset");
    std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

    // --- Register the training and test trees

    TTree *signal     = (TTree*)input->Get("ttbar_res100.0/ttbar_res_dilepton_CP5");
    TTree *background = (TTree*)input->Get("ttbar_res100.0/ttbar_res_dilepton_CP5");

    // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
    TString outfileName("../output/TMVA_4Vec.root");
    TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

//     TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
//             "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
    TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
            "!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );

    loader->AddVariable("Lep1_pt",'F');
    loader->AddVariable("Lep1_phi",'F');
    loader->AddVariable("Lep1_eta",'F');
    loader->AddVariable("Lep1_E",'F');
    loader->AddVariable("Lep1_flavor",'I',1,2);
    loader->AddVariable("Lep2_pt",'F');
    loader->AddVariable("Lep2_phi",'F');
    loader->AddVariable("Lep2_eta",'F');
    loader->AddVariable("Lep2_E",'F');
    loader->AddVariable("Lep2_flavor",'I',1,2);
    loader->AddVariable("Jet1_pt",'F');
    loader->AddVariable("Jet1_phi",'F');
    loader->AddVariable("Jet1_eta",'F');
    loader->AddVariable("Jet1_E",'F');
    loader->AddVariable("Jet1_bTagScore",'F');
    loader->AddVariable("Jet2_pt",'F');
    loader->AddVariable("Jet2_phi",'F');
    loader->AddVariable("Jet2_eta",'F');
    loader->AddVariable("Jet2_E",'F');
    loader->AddVariable("Jet2_bTagScore",'F');
    loader->AddVariable("nJets",'I');

    // global event weights per tree (see below for setting event-wise weights)
    Double_t signalWeight     = 1.0;
    Double_t backgroundWeight = 1.0;

    // You can add an arbitrary number of signal or background trees
    loader->AddSignalTree    ( signal,     signalWeight     );
    loader->AddBackgroundTree( background, backgroundWeight );

//     loader->SetBackgroundWeightExpression( "weight" );

    // Apply additional cuts on the signal and background samples (can be different)
    TCut met6_sig = "PuppiMET>230 && absmetres_PUPPI<40";
    TCut met6_bg = "PuppiMET>230 && absmetres_PUPPI>40";


    loader->PrepareTrainingAndTestTree( met6_sig, met6_bg,
                                         "nTrain_Signal=10000:nTrain_Background=10000:SplitMode=Random:NormMode=NumEvents:!V" );


    // improved neural network implementation
//       TString layoutString ("Layout=TANH|(N+100)*2,LINEAR");
//       TString layoutString ("Layout=SOFTSIGN|100,SOFTSIGN|50,SOFTSIGN|20,LINEAR");
//       TString layoutString ("Layout=RELU|300,RELU|100,RELU|30,RELU|10,LINEAR");
//       TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|30,SOFTSIGN|20,SOFTSIGN|10,LINEAR");
//       TString layoutString ("Layout=TANH|50,TANH|30,TANH|20,TANH|10,LINEAR");
//       TString layoutString ("Layout=SOFTSIGN|50,SOFTSIGN|20,LINEAR");
    TString layoutString ("Layout=TANH|100,TANH|50,TANH|10,LINEAR");

    TString training0 ("LearningRate=1e-1,Momentum=0.0,Repetitions=1,ConvergenceSteps=300,BatchSize=20,TestRepetitions=15,WeightDecay=0.001,Regularization=NONE,DropConfig=0.0+0.5+0.5+0.5,DropRepetitions=1,Multithreading=True");
    TString training1 ("LearningRate=1e-2,Momentum=0.5,Repetitions=1,ConvergenceSteps=300,BatchSize=30,TestRepetitions=7,WeightDecay=0.001,Regularization=L2,Multithreading=True,DropConfig=0.0+0.1+0.1+0.1,DropRepetitions=1");
    TString training2 ("LearningRate=1e-2,Momentum=0.3,Repetitions=1,ConvergenceSteps=300,BatchSize=40,TestRepetitions=7,WeightDecay=0.0001,Regularization=L2,Multithreading=True");
    TString training3 ("LearningRate=1e-3,Momentum=0.1,Repetitions=1,ConvergenceSteps=200,BatchSize=70,TestRepetitions=7,WeightDecay=0.0001,Regularization=NONE,Multithreading=True");

    TString trainingStrategyString ("TrainingStrategy=");
    trainingStrategyString += training0 + "|" + training1 + "|" + training2 + "|" + training3;


//       TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CROSSENTROPY");
    TString nnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=G:WeightInitialization=XAVIERUNIFORM");
//       TString nnOptions ("!H:V:VarTransform=Normalize:ErrorStrategy=CHECKGRADIENTS");
    nnOptions.Append (":");
    nnOptions.Append (layoutString);
    nnOptions.Append (":");
    nnOptions.Append (trainingStrategyString);

    factory->BookMethod(loader, TMVA::Types::kDNN, "DNN", nnOptions ); // NN


    // Train MVAs using the set of training events
    factory->TrainAllMethods();

    // ---- Evaluate all MVAs using the set of test events
    factory->TestAllMethods();

    // ----- Evaluate and compare performance of all configured MVAs
    factory->EvaluateAllMethods();

    // --------------------------------------------------------------

    // Save the output
    outputFile->Close();

    std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
    std::cout << "==> TMVAClassification is done!" << std::endl;

    delete factory;
    delete loader;
    // Launch the GUI for the root macros
    if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
}
