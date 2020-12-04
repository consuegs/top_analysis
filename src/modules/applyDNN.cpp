//Small script to save DNN output to existing Tree
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
#include <TTree.h>

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"

Config const &cfg=Config::get();

extern "C"
void run()
{   
    // ~TString sample="dilepton";
    TString sample="T2tt_650_350";
    
    TFile *dataFile=new TFile(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/output/DNNapplied/ttbar_res%.1f_new.root",cfg.processFraction*100),"update");
    TTree *Tree=(TTree *) dataFile->Get(TString::Format("ttbar_res%.1f",cfg.processFraction*100)+"/ttbar_res_"+sample);
    float PuppiMET,PFMET,METunc_Puppi,HT,nJets,Lep1_pt,Lep1_phi,Lep1_eta,Lep1_E,Lep1_flavor,Lep2_pt,Lep2_phi,Lep2_eta,Lep2_E,Lep2_flavor,Jet1_pt,Jet1_phi,Jet1_eta,Jet1_E,Jet2_pt,Jet2_phi,Jet2_eta,Jet2_E,n_Interactions, genMET;
    Float_t DNN_regression;
    
    TBranch* newBranch=Tree->Branch("DNN_regression",&DNN_regression,"DNN_regression/f");
    Tree->SetBranchAddress("genMET",&genMET);
    
    Tree->SetBranchAddress("MET",&PFMET);
    Tree->SetBranchAddress("PuppiMET",&PuppiMET);
    Tree->SetBranchAddress("METunc_Puppi",&METunc_Puppi);
    Tree->SetBranchAddress("HT",&HT);
    Tree->SetBranchAddress("nJets",&nJets);
    Tree->SetBranchAddress("n_Interactions",&n_Interactions);
    Tree->SetBranchAddress("Lep1_pt",&Lep1_pt);
    Tree->SetBranchAddress("Lep1_phi",&Lep1_phi);
    Tree->SetBranchAddress("Lep1_eta",&Lep1_eta);
    Tree->SetBranchAddress("Lep1_E",&Lep1_E);
    Tree->SetBranchAddress("Lep1_flavor",&Lep1_flavor);
    Tree->SetBranchAddress("Lep2_pt",&Lep2_pt);
    Tree->SetBranchAddress("Lep2_phi",&Lep2_phi);
    Tree->SetBranchAddress("Lep2_eta",&Lep2_eta);
    Tree->SetBranchAddress("Lep2_E",&Lep2_E);
    Tree->SetBranchAddress("Lep2_flavor",&Lep2_flavor);
    Tree->SetBranchAddress("Jet1_pt",&Jet1_pt);
    Tree->SetBranchAddress("Jet1_phi",&Jet1_phi);
    Tree->SetBranchAddress("Jet1_eta",&Jet1_eta);
    Tree->SetBranchAddress("Jet1_E",&Jet1_E);
    Tree->SetBranchAddress("Jet2_pt",&Jet2_pt);
    Tree->SetBranchAddress("Jet2_phi",&Jet2_phi);
    Tree->SetBranchAddress("Jet2_eta",&Jet2_eta);
    Tree->SetBranchAddress("Jet2_E",&Jet2_E);
    
    TMVA::PyMethodBase::PyInitialize();
    TMVA::Reader* reader_TMVA_Bin1=new TMVA::Reader("Color:!Silent");
    TMVA::Reader* reader_TMVA_Bin2=new TMVA::Reader("Color:!Silent");
    TMVA::Reader* reader_TMVA_Bin3=new TMVA::Reader("Color:!Silent");
    TMVA::Reader* reader_TMVA_Bin4=new TMVA::Reader("Color:!Silent");
    TMVA::Reader* reader_TMVA_Bin5=new TMVA::Reader("Color:!Silent");
    TMVA::Reader* reader_TMVA_Bin6=new TMVA::Reader("Color:!Silent");
    for(TMVA::Reader* tempreader:{reader_TMVA_Bin1, reader_TMVA_Bin2, reader_TMVA_Bin3, reader_TMVA_Bin4, reader_TMVA_Bin5, reader_TMVA_Bin6}){
      tempreader->AddVariable("PuppiMET", &PuppiMET);
      tempreader->AddVariable("METunc_Puppi", &METunc_Puppi);
      tempreader->AddVariable("MET", &PFMET);
      tempreader->AddVariable("HT", &HT);
      tempreader->AddVariable("nJets", &nJets);
      tempreader->AddVariable("n_Interactions", &n_Interactions);
      tempreader->AddVariable("Lep1_flavor", &Lep1_flavor);
      tempreader->AddVariable("Lep2_flavor", &Lep2_flavor);
      tempreader->AddVariable("Lep1_pt", &Lep1_pt);
      tempreader->AddVariable("Lep1_phi", &Lep1_phi);
      tempreader->AddVariable("Lep1_eta", &Lep1_eta);
      tempreader->AddVariable("Lep1_E", &Lep1_E);
      tempreader->AddVariable("Lep2_pt", &Lep2_pt);
      tempreader->AddVariable("Lep2_phi", &Lep2_phi);
      tempreader->AddVariable("Lep2_eta", &Lep2_eta);
      tempreader->AddVariable("Lep2_E", &Lep2_E);
      tempreader->AddVariable("Jet1_pt", &Jet1_pt);
      tempreader->AddVariable("Jet1_phi", &Jet1_phi);
      tempreader->AddVariable("Jet1_eta", &Jet1_eta);
      tempreader->AddVariable("Jet1_E", &Jet1_E);
      tempreader->AddVariable("Jet2_pt", &Jet2_pt);
      tempreader->AddVariable("Jet2_phi", &Jet2_phi);
      tempreader->AddVariable("Jet2_eta", &Jet2_eta);
      tempreader->AddVariable("Jet2_E", &Jet2_E);
      tempreader->AddSpectator("PuppiMET", &PuppiMET);   //Placeholder
      tempreader->AddSpectator("genMET", &Jet2_E);    //Placeholder
    }
    reader_TMVA_Bin1->BookMVA("PyKerasBin1", "dataset/weights/TMVARegression_PyKerasBin1.weights.xml");
    reader_TMVA_Bin2->BookMVA("PyKerasBin2", "dataset/weights/TMVARegression_PyKerasBin2.weights.xml");
    reader_TMVA_Bin3->BookMVA("PyKerasBin3", "dataset/weights/TMVARegression_PyKerasBin3.weights.xml");
    reader_TMVA_Bin4->BookMVA("PyKerasBin4", "dataset/weights/TMVARegression_PyKerasBin4.weights.xml");
    reader_TMVA_Bin5->BookMVA("PyKerasBin5", "dataset/weights/TMVARegression_PyKerasBin5.weights.xml");
    reader_TMVA_Bin6->BookMVA("PyKerasBin6", "dataset/weights/TMVARegression_PyKerasBin6.weights.xml");
    
    int totalEntries=Tree->GetEntriesFast();
    int iEv=0;
    for(Int_t ievent=0;ievent<totalEntries;ievent++) {
        if(Tree->GetEntry(ievent)<=0) break;
        iEv++;
        if (iEv%(std::max(totalEntries/10,1))==0){
            io::log*".";
            io::log.flush();
        }
        if(PuppiMET<0) DNN_regression=-1;
        else if(PuppiMET<40) DNN_regression=reader_TMVA_Bin1->EvaluateRegression("PyKerasBin1")[0];
        else if(PuppiMET<80) DNN_regression=reader_TMVA_Bin2->EvaluateRegression("PyKerasBin2")[0];
        else if(PuppiMET<120) DNN_regression=reader_TMVA_Bin3->EvaluateRegression("PyKerasBin3")[0];
        else if(PuppiMET<160) DNN_regression=reader_TMVA_Bin4->EvaluateRegression("PyKerasBin4")[0];
        else if(PuppiMET<230) DNN_regression=reader_TMVA_Bin5->EvaluateRegression("PyKerasBin5")[0];
        else DNN_regression=reader_TMVA_Bin6->EvaluateRegression("PyKerasBin6")[0];
        
        newBranch->Fill();
    }
        
    dataFile->cd(TString::Format("ttbar_res%.1f",cfg.processFraction*100));
    Tree->Write("ttbar_res_"+sample,TObject::kOverwrite);
}
