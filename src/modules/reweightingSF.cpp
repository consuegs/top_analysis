//Script to derive weights for reweighting
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

void saveProfileFrom2D(TH2F* hist, TString name,io::RootFileSaver* saver){
   TProfile profile;
   profile=*(hist->ProfileX("Profile"));
   saver->save(*hist,name);
   profile.SetTitle("Profile;pt_NuNu;mean weight");
   saver->save(profile,name+"_profile");
}

Config const &cfg=Config::get();

extern "C"
void run()
{   
    //Use constant weight per met-bin
    // ~bool binWeighting=true;
    bool binWeighting=false;
    
    //Define reweighting variables
    // ~float A=1.0;
    // ~float B=0.55;
    // ~float C=5.5;
    // ~float D=0.;
    // ~float E=0.;
    
    // ~float slope=0.001;
    float slope=0.002;
    // ~float slope=0.004;
    // ~float slope=-0.001;
    
    //Define binning
    int NBIN_MET_COARSE=4;
    int NBIN_PHI_COARSE=3;
    Double_t metBinsCoarse[NBIN_MET_COARSE+1]={0,40,80,120,230};
    Double_t phiBinsCoarse[NBIN_PHI_COARSE+1]={0,0.7,1.4,3.141};
    TUnfoldBinning *generatorBinning=new TUnfoldBinning("generator");
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
    TH1 *histGen=signalBinning->CreateHistogram("histGen");
    TH1 *histGen_reweight=signalBinning->CreateHistogram("histGen_reweight");
    
    
    // ~TFile *dataFile=new TFile("/net/data_cms1b/user/dmeuser/top_analysis/output/"+TString::Format("ttbar_res%.1f.root",cfg.processFraction*100));
    TFile *dataFile=new TFile(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/ttbar_res%.1f.root",cfg.processFraction*100),"update");
    TTree *PowhegTree=(TTree *) dataFile->Get(TString::Format("ttbar_res%.1f",cfg.processFraction*100)+"/ttbar_res_dilepton");
    Float_t metGen,phiGen,mcWeight,metRec,genMET;
    UInt_t genDecayMode;
    ULong64_t evtNo;
    float weight;
    
    // ~PowhegTree->ResetBranchAddresses();
    TBranch* bpt=PowhegTree->Branch("reweight_PTnunu",&weight,"weight/f");
    PowhegTree->SetBranchAddress("N",&mcWeight);
    PowhegTree->SetBranchAddress("PtNuNu",&metGen);
    PowhegTree->SetBranchAddress("Phi_NuNu",&phiGen);
    PowhegTree->SetBranchAddress("genDecayMode",&genDecayMode);
    PowhegTree->SetBranchAddress("MET",&metRec);
    PowhegTree->SetBranchAddress("evtNo",&evtNo);
    PowhegTree->SetBranchAddress("genMET",&genMET);
    
    for(Int_t ievent=0;ievent<PowhegTree->GetEntriesFast();ievent++) {
        if(PowhegTree->GetEntry(ievent)<=0) break;
        // ~if(metRec<0 || genDecayMode>3 || metGen<0) continue;
        
        Int_t genBin=signalBinning->GetGlobalBinNumber(metGen,phiGen);
        histGen->Fill(genBin,mcWeight);
        
        float x=genMET/2000.;
        // ~weight=(A*std::pow(x,B)*pow(1-x,C)*(1+D*x+E*x*x));
        weight=std::max((1.+(metGen-100.)*slope)*(1.+(metGen-100.)*slope),0.1);
        // ~weight=std::max((1.+(metGen-100.)*slope),0.1);
        if(metGen<0) weight=0;
        
        if(binWeighting) weight=((genBin-1)%4)+1;
        
        histGen_reweight->Fill(genBin,mcWeight*weight);
    }
    
    std::cout<<"okay"<<std::endl;
    
    float norm=1.0/histGen_reweight->Integral()*histGen->Integral();
    histGen_reweight->Scale(norm);
    
    for(int i=1; i<=12; i++){
        std::cout<<histGen_reweight->GetBinContent(i)<<"   "<<histGen->GetBinContent(i)<<std::endl;
    }
    
    TH1* hist_mean=(TH1*)histGen->Clone();
    hist_mean->Add(histGen_reweight);
    
    histGen_reweight->Add(histGen,-1.);
    histGen_reweight->Divide(hist_mean);
    
    for(int i=1; i<=12; i++){
        std::cout<<histGen_reweight->GetBinContent(i)<<std::endl;
    }
    
    TH2F hist2D_weights("","",100,0,800,1000,0,10);
    TH2F hist2D_weights_binning("","",13,-0.5,12.5,1000,-10,10);
    
    for(Int_t ievent=0;ievent<PowhegTree->GetEntriesFast();ievent++) {
        if(PowhegTree->GetEntry(ievent)<=0) break;
        Int_t genBin=signalBinning->GetGlobalBinNumber(metGen,phiGen);
        float x=genMET/2000.;
        // ~weight=(A*std::pow(x,B)*pow(1-x,C)*(1+D*x+E*x*x))*norm;
        weight=std::max((1.+(metGen-100.)*slope)*(1.+(metGen-100.)*slope),0.1)*norm;
        // ~weight=std::max((1.+(metGen-100.)*slope),0.1)*norm;
        if(metGen<0) weight=0;
                
        if(binWeighting) weight=(((genBin-1)%4)+1)*norm;
        
        hist2D_weights.Fill(metGen,weight);
        
        hist2D_weights_binning.Fill(genBin,weight);
        
        bpt->Fill();
    }
        
    dataFile->cd(TString::Format("ttbar_res%.1f",cfg.processFraction*100));
    PowhegTree->Write("ttbar_res_dilepton_PTreweight",TObject::kOverwrite);
    
    io::RootFileSaver saver(TString::Format("TUnfold%.1f.root",cfg.processFraction*100),"reweightingSF");
    saveProfileFrom2D(&hist2D_weights,"weights/weights_ptNuNu_2D",&saver);
    saveProfileFrom2D(&hist2D_weights_binning,"weights/weights_binning_2D",&saver);
}
