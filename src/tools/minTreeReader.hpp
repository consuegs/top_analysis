#ifndef MINITREEREADER_H
#define MINITREEREADER_H

#include <vector>
#include <TROOT.h>
#include <TTree.h>
#include <TBranch.h>
#include "tools/systematics.hpp"


using namespace std;

const UInt_t length = 20;
const UInt_t lengthLong = 150;

class minTreeReader {

    public:
    minTreeReader(TTree& tree, const Systematic::Systematic& syst);
        
    UInt_t ee;
    UInt_t mumu;
    UInt_t emu;
    float MET;
    float PtNuNu;
    float Phi_rec;
    float Phi_gen;
    float Phi_NuNu;
    float dPhiMETnearJet;
    float dPhiMETfarJet;
    float dPhiMETleadJet;
    float dPhiMETlead2Jet;
    float dPhiMETbJet;
    float dPhiLep1Lep2;
    float dPhiJet1Jet2;
    float METsig;
    float N;
    float SF;
    float totalWeight;
    UInt_t runNo;
    UInt_t lumNo;
    ULong64_t evtNo;
    UInt_t runEra;
    UInt_t genDecayMode;
    float genMet;
    float PuppiMET;
    float XYcorrMET;
    float HT;
    float HT_phi;
    float MHT;
    float MT;
    float genMT;
    float MT_nextLep;
    float genMT_nextLep;
    float n_Interactions;
    UInt_t n_Interactions_gen;
    float dPhiPtnunuMet;
    float leadTop_pT;
    float dPhiNuNu;
    float Phi_recPuppi;
    float Phi_recXYcorr;
    UInt_t looseLeptonVeto;
    float nJets;
    float dPhiMETnearJet_Puppi;
    float dPhiMETfarJet_Puppi;
    float dPhiMETleadJet_Puppi;
    float dPhiMETlead2Jet_Puppi;
    float dPhiMETbJet_Puppi;
    float dPhiLep1bJet;
    float dPhiLep1Jet1;
    float ratioMET_sqrtMETunc_Puppi;
    float ratio_pTj1_vecsum_pT_l1_l2_bjet;
    float METunc_Puppi;
    float METunc_PF;
    float absmetres_PUPPI;
    float Lep1_pt;
    float Lep1_phi;
    float Lep1_eta;
    float Lep1_E;
    float Lep1_flavor;
    float Lep2_pt;
    float Lep2_phi;
    float Lep2_eta;
    float Lep2_E;
    float Lep2_flavor;
    float Jet1_pt;
    float Jet1_phi;
    float Jet1_eta;
    float Jet1_E;
    float Jet1_bTagScore;
    float Jet1_unc;
    float Jet2_pt;
    float Jet2_phi;
    float Jet2_eta;
    float Jet2_E;
    float Jet2_bTagScore;
    float Jet2_unc;
    float mLL;
    float PFMET_phi;
    float PuppiMET_phi;
    float CaloMET;
    float CaloMET_phi;
    float genMET_phi;
    float PtNuNu_phi;
    float DNN_MET_pT;
    float DNN_MET_phi;
    float DNN_MET_dPhi_nextLep;
    float MT2;
    float vecsum_pT_allJet;
    float vecsum_pT_l1l2_allJet;
    float mass_l1l2_allJet;
    float ratio_vecsumpTlep_vecsumpTjet;
    float mjj;
    
    std::map<Systematic::Type,std::pair<float,float>> systTotalWeights;
    
};
#endif
