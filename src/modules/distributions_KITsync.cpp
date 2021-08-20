#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"

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

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"

#include "cppflow/ops.h"
#include "cppflow/model.h"

Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   TString dssName_multi="";
   //Read TriggerSF hists
   io::RootFileReader triggerSF_ee(cfg.trigger_SF_ee,"");
   TH2F* triggerSF_ee_hist = (TH2F*)(triggerSF_ee.read<TH2F>("eff_histo"));
   io::RootFileReader triggerSF_mumu(cfg.trigger_SF_mumu,"");
   TH2F* triggerSF_mumu_hist = (TH2F*)(triggerSF_mumu.read<TH2F>("eff_histo"));
   io::RootFileReader triggerSF_emu(cfg.trigger_SF_emu,"");
   TH2F* triggerSF_emu_hist = (TH2F*)(triggerSF_emu.read<TH2F>("eff_histo"));
      
   for (TString ds_name: cfg.datasets.getDatasetNames()){
      auto ds=cfg.datasets.getDataset(ds_name);
      
      io::RootFileSaver ttbar_res_saver(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees_KITsync/%.1f/%s.root",cfg.year.Data(),cfg.treeVersion.Data(),cfg.processFraction*100,TString(ds_name).Data()),TString::Format("ttbar_res%.1f",cfg.processFraction*100),true,false);
      TTree ttbar_res("ttbar_res","ttbar_res");
      
      int runEra=0;     //int to store run Era in minimal Trees
      
      for (auto dss: cfg.datasets.getDatasubsets({ds.name})){   
         TFile file(dss.getPath(),"read");
         if (file.IsZombie()) {
            return;
         }
         io::log * ("Processing '"+dss.name+"' ");

         bool const isData=dss.isData;
         bool const isSignal=dss.isSignal;
         int year_int=cfg.year_int;
         
         runEra++;
         if(!isData) runEra=0;
         
         //Lumi weight for current sample
         float lumi_weight=dss.xsec/float(dss.Ngen)*cfg.lumi;
         
         //Check if current sample is TTbar powheg dilepton
         bool ttBar_dilepton=false;
         if (dss.datasetName=="TTbar_diLepton") ttBar_dilepton=true;
         
         //Check if current sample is TTbar powheg dilepton tau
         bool ttBar_dilepton_tau=false;
         if (dss.datasetName=="TTbar_diLepton_tau") ttBar_dilepton_tau=true;
         
         //Check if current sample is TTbar powheg semileptonic
         bool ttBar_singleLepton=false;
         if (dss.datasetName=="TTbar_singleLepton") ttBar_singleLepton=true;
         
         //Check if current sample is TTbar powheg hadronic
         bool ttBar_hadronic=false;
         if (dss.datasetName=="TTbar_hadronic") ttBar_hadronic=true;
         
         //Check if current sample is TTbar amc@NLO
         bool ttBar_amc=false;
         if (dss.datasetName=="TTbar_amcatnlo") ttBar_amc=true;
         
         //Check if current sample is selectedSusy scenario 
         bool SUSY_T2tt_650_350=false;
         if (dss.datasetName=="T2tt_650_350") SUSY_T2tt_650_350=true;
         
         //Check if current sample is selectedDM scenario 
         bool DM_scalar_1_200=false;
         if (dss.datasetName=="DM_scalar_1_200") DM_scalar_1_200=true;
         
         //Check if current sample is SUSY scenario
         bool SUSY_scen=false;
         if (dss.name.find("SMS")!=std::string::npos) SUSY_scen=true;
         
         //Check which PD is currently used
         bool SingleElectron=false;
         bool SingleMuon=false;
         bool DoubleEG=false;
         bool DoubleMuon=false;
         bool MuonEG=false;
         bool EGamma=false;
         if (dss.datasetName.find("SingleElectron")!=std::string::npos) SingleElectron=true;
         else if (dss.datasetName.find("SingleMuon")!=std::string::npos) SingleMuon=true;
         else if (dss.datasetName.find("DoubleEG")!=std::string::npos) DoubleEG=true;
         else if (dss.datasetName.find("DoubleMuon")!=std::string::npos) DoubleMuon=true;
         else if (dss.datasetName.find("MuonEG")!=std::string::npos) MuonEG=true;
         else if (dss.datasetName.find("EGamma")!=std::string::npos) EGamma=true;
         
         //Check if current sample is Run2016H
         bool Run2016H=false;
         if (dss.name.find("Run2016H")!=std::string::npos) Run2016H=true;
         
         //Check if current sample is Run2017AB
         bool Run2017AB=false;
         if (dss.name.find("Run2017A")!=std::string::npos) Run2017AB=true;
         else if (dss.name.find("Run2017B")!=std::string::npos) Run2017AB=true;
         
         //Set boolean for savin minimalTree
         bool minimalTree=true;
         
         //Ntuple and file to save minimal ttbar tree used for binning studies
         float minTree_MET, minTree_PtNuNu, minTree_PhiRec, minTree_PhiGen, minTree_PhiNuNu, minTree_PhiMetNearJet, minTree_PhiMetFarJet, minTree_PhiMetLeadJet, minTree_PhiMetLead2Jet,
         minTree_PhiMetbJet, minTree_dPhiLep1Lep2, minTree_dPhiJet1Jet2, minTree_METsig, minTree_N, minTree_SF, minTree_genMet, minTree_PuppiMet, minTree_XYcorrMet, 
         minTree_HT, minTree_HT_phi, minTree_MHT, minTree_MT, minTree_genMT, minTree_MT_nextLep, minTree_genMT_nextLep,
         minTree_PhiPtnunuMet, minTree_leadTop, minTree_dPhiNuNu, minTree_PhiRecPuppi, minTree_PhiRecXYcorr, minTree_PhiMetNearJet_Puppi, minTree_PhiMetFarJet_Puppi,
         minTree_PhiMetLeadJet_Puppi, minTree_PhiMetLead2Jet_Puppi, minTree_PhiMetbJet_Puppi, minTree_dPhiLep1bJet, minTree_dPhiLep1Jet1, minTree_ratioMET_sqrtMETunc_Puppi,
         minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet, minTree_METunc_Puppi, minTree_METunc_PF, minTree_absmetres_PUPPI,
         minTree_Lep1_pt, minTree_Lep1_phi, minTree_Lep1_eta, minTree_Lep1_E, minTree_Lep1_flavor,
         minTree_Lep2_pt, minTree_Lep2_phi, minTree_Lep2_eta, minTree_Lep2_E, minTree_Lep2_flavor,
         minTree_Jet1_pt, minTree_Jet1_phi, minTree_Jet1_eta, minTree_Jet1_E, minTree_Jet1_bTagScore, minTree_Jet1_unc,
         minTree_Jet2_pt, minTree_Jet2_phi, minTree_Jet2_eta, minTree_Jet2_E, minTree_Jet2_bTagScore, minTree_Jet2_unc,
         minTree_PFMET_phi, minTree_PuppiMET_phi, minTree_CaloMET, minTree_CaloMET_phi, minTree_genMET_phi, minTree_PtNuNu_phi, minTree_nJets, minTree_n_Interactions, minTree_DNN_MET_pT, minTree_DNN_MET_phi,minTree_DNN_MET_dPhi_nextLep,
         minTree_mLL, minTree_MT2, minTree_vecsum_pT_allJet, minTree_vecsum_pT_l1l2_allJet, minTree_mass_l1l2_allJet, minTree_ratio_vecsumpTlep_vecsumpTjet, minTree_mjj;
         UInt_t minTree_runNo, minTree_lumNo, minTree_genDecayMode, minTree_n_Interactions_gen, minTree_looseLeptonVeto, minTree_NpromptNeutrinos, minTree_NnonpromptNeutrinos, minTree_ee, minTree_mumu, minTree_emu;
         ULong64_t minTree_evtNo;
         if(minimalTree){
            ttbar_res.Branch("ee",&minTree_ee,"ee/i");
            ttbar_res.Branch("mumu",&minTree_mumu,"mumu/i");
            ttbar_res.Branch("emu",&minTree_emu,"emu/i");
            ttbar_res.Branch("MET",&minTree_MET,"MET/f");
            ttbar_res.Branch("PtNuNu",&minTree_PtNuNu,"PtNuNu/f");
            ttbar_res.Branch("Phi_rec",&minTree_PhiRec,"Phi_rec/f");
            ttbar_res.Branch("Phi_gen",&minTree_PhiGen,"Phi_gen/f");
            ttbar_res.Branch("Phi_NuNu",&minTree_PhiNuNu,"Phi_NuNu/f");
            ttbar_res.Branch("dPhiMETnearJet",&minTree_PhiMetNearJet,"dPhiMETnearJet/f");
            ttbar_res.Branch("dPhiMETfarJet",&minTree_PhiMetFarJet,"dPhiMETfarJet/f");
            ttbar_res.Branch("dPhiMETleadJet",&minTree_PhiMetLeadJet,"dPhiMETleadJet/f");
            ttbar_res.Branch("dPhiMETlead2Jet",&minTree_PhiMetLead2Jet,"dPhiMETlead2Jet/f");
            ttbar_res.Branch("dPhiMETbJet",&minTree_PhiMetbJet,"dPhiMETbJet/f");
            ttbar_res.Branch("dPhiLep1Lep2",&minTree_dPhiLep1Lep2,"dPhiLep1Lep2/f");
            ttbar_res.Branch("dPhiJet1Jet2",&minTree_dPhiJet1Jet2,"dPhiJet1Jet2/f");
            ttbar_res.Branch("METsig",&minTree_METsig,"METsig/f");
            ttbar_res.Branch("N",&minTree_N,"N/f");
            ttbar_res.Branch("SF",&minTree_SF,"SF/f");
            ttbar_res.Branch("runNo",&minTree_runNo,"runNo/i");
            ttbar_res.Branch("lumNo",&minTree_lumNo,"lumNo/i");
            ttbar_res.Branch("evtNo",&minTree_evtNo,"evtNo/l");
            ttbar_res.Branch("runEra",&runEra,"runEra/i");
            ttbar_res.Branch("genDecayMode",&minTree_genDecayMode,"genDecayMode/i");
            ttbar_res.Branch("genMET",&minTree_genMet,"genMET/f");
            ttbar_res.Branch("PuppiMET",&minTree_PuppiMet,"PuppiMET/f");
            ttbar_res.Branch("XYcorrMET",&minTree_XYcorrMet,"XYcorrMET/f");
            ttbar_res.Branch("HT",&minTree_HT,"HT/f");
            ttbar_res.Branch("HT_phi",&minTree_HT_phi,"HT_phi/f");
            ttbar_res.Branch("MHT",&minTree_MHT,"MHT/f");
            ttbar_res.Branch("MT",&minTree_MT,"MT/f");
            ttbar_res.Branch("genMT",&minTree_genMT,"genMT/f");
            ttbar_res.Branch("MT_nextLep",&minTree_MT_nextLep,"MT_nextLep/f");
            ttbar_res.Branch("genMT_nextLep",&minTree_genMT_nextLep,"genMT_nextLep/f");
            ttbar_res.Branch("n_Interactions",&minTree_n_Interactions,"n_Interactions/f");
            ttbar_res.Branch("n_Interactions_gen",&minTree_n_Interactions_gen,"n_Interactions_gen/i");
            ttbar_res.Branch("dPhiPtnunuMet",&minTree_PhiPtnunuMet,"dPhiPtnunuMet/f");
            ttbar_res.Branch("leadTop_pT",&minTree_leadTop,"leadTop_pT/f");
            ttbar_res.Branch("dPhiNuNu",&minTree_dPhiNuNu,"dPhiNuNu/f");
            ttbar_res.Branch("Phi_recPuppi",&minTree_PhiRecPuppi,"Phi_recPuppi/f");
            ttbar_res.Branch("Phi_recXYcorr",&minTree_PhiRecXYcorr,"Phi_recXYcorr/f");
            ttbar_res.Branch("looseLeptonVeto",&minTree_looseLeptonVeto,"looseLeptonVeto/i");
            ttbar_res.Branch("nJets",&minTree_nJets,"nJets/f");
            ttbar_res.Branch("dPhiMETnearJet_Puppi",&minTree_PhiMetNearJet_Puppi,"dPhiMETnearJet_Puppi/f");
            ttbar_res.Branch("dPhiMETfarJet_Puppi",&minTree_PhiMetFarJet_Puppi,"dPhiMETfarJet_Puppi/f");
            ttbar_res.Branch("dPhiMETleadJet_Puppi",&minTree_PhiMetLeadJet_Puppi,"dPhiMETleadJet_Puppi/f");
            ttbar_res.Branch("dPhiMETlead2Jet_Puppi",&minTree_PhiMetLead2Jet_Puppi,"dPhiMETlead2Jet_Puppi/f");
            ttbar_res.Branch("dPhiMETbJet_Puppi",&minTree_PhiMetbJet_Puppi,"dPhiMETbJet_Puppi/f");
            ttbar_res.Branch("dPhiLep1bJet",&minTree_dPhiLep1bJet,"dPhiLep1bJet/f");
            ttbar_res.Branch("dPhiLep1Jet1",&minTree_dPhiLep1Jet1,"dPhiLep1Jet1/f");
            ttbar_res.Branch("ratioMET_sqrtMETunc_Puppi",&minTree_ratioMET_sqrtMETunc_Puppi,"ratioMET_sqrtMETunc_Puppi/f");
            ttbar_res.Branch("ratio_pTj1_vecsum_pT_l1_l2_bjet",&minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet,"ratio_pTj1_vecsum_pT_l1_l2_bjet/f");
            ttbar_res.Branch("METunc_Puppi",&minTree_METunc_Puppi,"METunc_Puppi/f");
            ttbar_res.Branch("METunc_PF",&minTree_METunc_PF,"METunc_PF/f");
            ttbar_res.Branch("absmetres_PUPPI",&minTree_absmetres_PUPPI,"absmetres_PUPPI/f");
            ttbar_res.Branch("Lep1_pt",&minTree_Lep1_pt,"Lep1_pt/f");
            ttbar_res.Branch("Lep1_phi",&minTree_Lep1_phi,"Lep1_phi/f");
            ttbar_res.Branch("Lep1_eta",&minTree_Lep1_eta,"Lep1_eta/f");
            ttbar_res.Branch("Lep1_E",&minTree_Lep1_E,"Lep1_E/f");
            ttbar_res.Branch("Lep1_flavor",&minTree_Lep1_flavor,"Lep1_flavor/f");
            ttbar_res.Branch("Lep2_pt",&minTree_Lep2_pt,"Lep2_pt/f");
            ttbar_res.Branch("Lep2_phi",&minTree_Lep2_phi,"Lep2_phi/f");
            ttbar_res.Branch("Lep2_eta",&minTree_Lep2_eta,"Lep2_eta/f");
            ttbar_res.Branch("Lep2_E",&minTree_Lep2_E,"Lep2_E/f");
            ttbar_res.Branch("Lep2_flavor",&minTree_Lep2_flavor,"Lep2_flavor/f");
            ttbar_res.Branch("Jet1_pt",&minTree_Jet1_pt,"Jet1_pt/f");
            ttbar_res.Branch("Jet1_phi",&minTree_Jet1_phi,"Jet1_phi/f");
            ttbar_res.Branch("Jet1_eta",&minTree_Jet1_eta,"Jet1_eta/f");
            ttbar_res.Branch("Jet1_E",&minTree_Jet1_E,"Jet1_E/f");
            ttbar_res.Branch("Jet1_bTagScore",&minTree_Jet1_bTagScore,"Jet1_bTagScore/f");
            ttbar_res.Branch("Jet1_unc",&minTree_Jet1_unc,"Jet1_unc/f");
            ttbar_res.Branch("Jet2_pt",&minTree_Jet2_pt,"Jet2_pt/f");
            ttbar_res.Branch("Jet2_phi",&minTree_Jet2_phi,"Jet2_phi/f");
            ttbar_res.Branch("Jet2_eta",&minTree_Jet2_eta,"Jet2_eta/f");
            ttbar_res.Branch("Jet2_E",&minTree_Jet2_E,"Jet2_E/f");
            ttbar_res.Branch("Jet2_bTagScore",&minTree_Jet2_bTagScore,"Jet2_bTagScore/f");
            ttbar_res.Branch("Jet2_unc",&minTree_Jet2_unc,"Jet2_unc/f");
            ttbar_res.Branch("mLL",&minTree_mLL,"mLL/f");
            ttbar_res.Branch("PFMET_phi",&minTree_PFMET_phi,"PFMET_phi/f");
            ttbar_res.Branch("PuppiMET_phi",&minTree_PuppiMET_phi,"PuppiMET_phi/f");
            ttbar_res.Branch("CaloMET",&minTree_CaloMET,"CaloMET/f");
            ttbar_res.Branch("CaloMET_phi",&minTree_CaloMET_phi,"CaloMET_phi/f");
            ttbar_res.Branch("genMET_phi",&minTree_genMET_phi,"genMET_phi/f");
            ttbar_res.Branch("PtNuNu_phi",&minTree_PtNuNu_phi,"PtNuNu_phi/f");
            ttbar_res.Branch("NpromptNeutrinos",&minTree_NpromptNeutrinos,"NpromptNeutrinos/i");
            ttbar_res.Branch("NnonpromptNeutrinos",&minTree_NnonpromptNeutrinos,"NnonpromptNeutrinos/i");
            ttbar_res.Branch("DNN_MET_pT",&minTree_DNN_MET_pT,"DNN_MET_pT/f");
            ttbar_res.Branch("DNN_MET_phi",&minTree_DNN_MET_phi,"DNN_MET_phi/f");
            ttbar_res.Branch("DNN_MET_dPhi_nextLep",&minTree_DNN_MET_dPhi_nextLep,"DNN_MET_dPhi_nextLep/f");
            ttbar_res.Branch("MT2",&minTree_MT2,"MT2/f");
            ttbar_res.Branch("vecsum_pT_allJet",&minTree_vecsum_pT_allJet,"vecsum_pT_allJet/f");
            ttbar_res.Branch("vecsum_pT_l1l2_allJet",&minTree_vecsum_pT_l1l2_allJet,"vecsum_pT_l1l2_allJet/f");
            ttbar_res.Branch("mass_l1l2_allJet",&minTree_mass_l1l2_allJet,"mass_l1l2_allJet/f");
            ttbar_res.Branch("ratio_vecsumpTlep_vecsumpTjet",&minTree_ratio_vecsumpTlep_vecsumpTjet,"ratio_vecsumpTlep_vecsumpTjet/f");
            ttbar_res.Branch("mjj",&minTree_mjj,"mjj/f");
         }
         
         //Set Tree Input variables
         TTreeReader reader(cfg.treeName, &file);
         TTreeReaderValue<float> w_pu(reader, "pu_weight");
         TTreeReaderValue<UInt_t> runNo(reader, "runNo");
         TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
         TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
         TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
         TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
         TTreeReaderValue<float> w_topPT(reader, "topPTweight");
         TTreeReaderValue<float> w_bTag(reader, (year_int==1)? "bTagWeight_DeepCSV" : "bTagWeight");     //Use DeepCSV for 2016 at the moment
         TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
         TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
         // ~TTreeReaderValue<std::vector<tree::Electron>> electrons_add(reader, "electrons_add");
         // ~TTreeReaderValue<std::vector<tree::Muon>>     muons_add    (reader, "muons_add");
         TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
         TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
         TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
         TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
         TTreeReaderValue<tree::MET> MET(reader, "met");
         TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
         TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
         TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
         TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
         TTreeReaderValue<tree::MET> MET_XYcorr(reader, "metXYcorr");
         TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
         TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
         TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
         TTreeReaderValue<int> n_Interactions(reader, "nPV");
         TTreeReaderValue<int> n_Interactions_gen(reader, "true_nPV");
         TTreeReaderValue<float> HTgen(reader, "genHt");
         TTreeReaderValue<bool> is_ee   (reader, "ee");
         TTreeReaderValue<bool> is_emu   (reader, "emu");
         TTreeReaderValue<bool> is_mumu   (reader, "mumu");
         TTreeReaderValue<bool> addLepton   (reader, "addLepton");
         TTreeReaderValue<float> mll   (reader, "mll");
         TTreeReaderValue<float> mt2   (reader, "MT2");
         TTreeReaderValue<float> genMT2   (reader, "genMT2");
         TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
         TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
         TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
         TTreeReaderValue<bool> muonTrigg1(reader, cfg.muonTrigg1);
         TTreeReaderValue<bool> muonTrigg2(reader, cfg.muonTrigg2);
         TTreeReaderValue<bool> muonTrigg3(reader, cfg.muonTrigg3);
         TTreeReaderValue<bool> muonTrigg4(reader, cfg.muonTrigg4);
         TTreeReaderValue<bool> singleMuonTrigg1(reader, cfg.singleMuonTrigg1);
         TTreeReaderValue<bool> singleMuonTrigg2(reader, cfg.singleMuonTrigg2);
         TTreeReaderValue<bool> eleTrigg1(reader, cfg.eleTrigg1);
         TTreeReaderValue<bool> eleTrigg2(reader, cfg.eleTrigg2);
         TTreeReaderValue<bool> eleMuTrigg1(reader, cfg.eleMuTrigg1);
         TTreeReaderValue<bool> eleMuTrigg2(reader, cfg.eleMuTrigg2);
         TTreeReaderValue<bool> eleMuTrigg3(reader, cfg.eleMuTrigg3);
         TTreeReaderValue<bool> eleMuTrigg4(reader, cfg.eleMuTrigg4);
         TTreeReaderValue<bool> singleEleTrigg(reader, cfg.singleEleTrigg);
         TTreeReaderValue<TLorentzVector> genTop(reader, "pseudoTop");
         TTreeReaderValue<TLorentzVector> genAntiTop(reader, "pseudoAntiTop");
         TTreeReaderValue<TLorentzVector> genLepton(reader, "pseudoLepton");
         TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "pseudoAntiLepton");
         TTreeReaderValue<TLorentzVector> genTau(reader, "pseudoTau");
         TTreeReaderValue<TLorentzVector> genAntiTau(reader, "pseudoAntiTau");
         TTreeReaderValue<int> genLeptonPdgId(reader, "pseudoLeptonPdgId");
         TTreeReaderValue<int> genAntiLeptonPdgId(reader, "pseudoAntiLeptonPdgId");
         TTreeReaderValue<TLorentzVector> genB(reader, "pseudoBJet");
         TTreeReaderValue<TLorentzVector> genAntiB(reader, "pseudoAntiBJet");
         TTreeReaderValue<TLorentzVector> genNeutrino(reader, "pseudoNeutrino");
         TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "pseudoAntiNeutrino");
         TTreeReaderValue<TLorentzVector> genWMinus(reader, "pseudoWMinus");
         TTreeReaderValue<TLorentzVector> genWPlus(reader, "pseudoWPlus");
         TTreeReaderValue<int> genDecayMode_pseudo(reader, "ttbarPseudoDecayMode");
         // ~TTreeReaderValue<TLorentzVector> genTop(reader, "genTop");    //Alternative genParticles, which are not based on RIVET
         // ~TTreeReaderValue<TLorentzVector> genAntiTop(reader, "genAntiTop");
         // ~TTreeReaderValue<TLorentzVector> genLepton(reader, "genLepton");
         // ~TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "genAntiLepton");
         // ~TTreeReaderValue<TLorentzVector> genTau(reader, "genTau");
         // ~TTreeReaderValue<TLorentzVector> genAntiTau(reader, "genAntiTau");
         // ~TTreeReaderValue<int> genLeptonPdgId(reader, "genLeptonPdgId");
         // ~TTreeReaderValue<int> genAntiLeptonPdgId(reader, "genAntiLeptonPdgId");
         // ~TTreeReaderValue<TLorentzVector> genB(reader, "genB");
         // ~TTreeReaderValue<TLorentzVector> genAntiB(reader, "genAntiB");
         // ~TTreeReaderValue<TLorentzVector> genNeutrino(reader, "genNeutrino");
         // ~TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "genAntiNeutrino");
         // ~TTreeReaderValue<TLorentzVector> genWMinus(reader, "genWMinus");
         // ~TTreeReaderValue<TLorentzVector> genWPlus(reader, "genWPlus");
         TTreeReaderValue<int> genDecayMode(reader, "ttbarDecayMode");
      
      
         int iEv=0;
         int processEvents=cfg.processFraction*dss.entries;
         while (reader.Next()){
            iEv++;
            if (iEv>processEvents) break;
            if (iEv%(std::max(processEvents/10,1))==0){
               io::log*".";
               io::log.flush();
            }
            
            float met=MET->p.Pt();
            float const met_puppi=MET_Puppi->p.Pt();
            float const genMet=GENMET->p.Pt();
            
            //Booleans for reco and pseudo selection
            // ~bool rec_selection=false;
            bool rec_selection=false;
            
            //Do not use tau events in signal sample
            if (ttBar_dilepton && *genDecayMode>3) continue;
            
            //Do only use ee,emu,mumu in in amc ttbar
            if (ttBar_amc && (*genDecayMode>3 || *genDecayMode==0)) continue;
            
            //Trigger selection
            std::vector<bool> diElectronTriggers={*eleTrigg1,*eleTrigg2,*singleEleTrigg};
            std::vector<bool> diMuonTriggers={*muonTrigg1,*muonTrigg2,*muonTrigg3,*muonTrigg4,*singleMuonTrigg1,*singleMuonTrigg2};
            std::vector<bool> electronMuonTriggers={*eleMuTrigg1,*eleMuTrigg2,*eleMuTrigg3,*eleMuTrigg4,*singleMuonTrigg1,*singleMuonTrigg2,*singleEleTrigg};
            std::vector<bool> channel={*is_ee,*is_mumu,*is_emu};
            std::vector<bool> PD={SingleElectron,DoubleEG,SingleMuon,DoubleMuon,MuonEG,EGamma};
            bool triggerMC=true;
            
            if (!isData){
               triggerMC=selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,false,year_int);
            }
            else{
               if(!selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,true,year_int,PD,Run2016H,Run2017AB)) continue;
            }
            
            //Baseline selection (separation into ee, emu, mumu already done at TreeWriter)
            TLorentzVector p_l1;
            TLorentzVector p_l2;
            int flavor_l1=0;  //1 for electron and 2 for muon
            int flavor_l2=0;
            bool muonLead=true; //Boolean for emu channel
            TString cat="";
            if(isSignal) *w_topPT=1.;  //Set top weight for signals to 1 
            
            rec_selection=selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,cat,muonLead);
            
            if (triggerMC==false) rec_selection=false;
            
            
            std::vector<tree::Jet> cjets;
            std::vector<tree::Jet> BJets;
            std::vector<bool> kitSyncSelection=selection::kitSyncSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
            
            if(!std::all_of(kitSyncSelection.begin(), kitSyncSelection.end(), [](bool v) { return v; })) rec_selection=false;
                                          
            // end reco baseline selection
            
            // Get pT of Neutrino Pair, which is further changed in case of BSM scenarios!!
            TLorentzVector neutrinoPair(0,0,0,0);
            neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
            
            
            if(rec_selection==false) continue;  // only proceed with events selected
            
            //Set Event weights
            float fEventWeight = 1.;
            float SFWeight = 1.;
            float triggerSF = 1.;
            float p_l1_trigg = std::min(199.,p_l1.Pt());
            float p_l2_trigg = std::min(199.,p_l2.Pt());
            if(!isData){
               if (*is_ee) triggerSF = triggerSF_ee_hist->GetBinContent(triggerSF_ee_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_ee_hist->GetYaxis()->FindBin(p_l2_trigg));
               else if (*is_mumu) triggerSF = triggerSF_mumu_hist->GetBinContent(triggerSF_mumu_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_mumu_hist->GetYaxis()->FindBin(p_l2_trigg));
               else if (*is_emu){
                  if(muonLead)triggerSF = triggerSF_emu_hist->GetBinContent(triggerSF_emu_hist->GetXaxis()->FindBin(p_l1_trigg), triggerSF_emu_hist->GetYaxis()->FindBin(p_l2_trigg));
                  else triggerSF = triggerSF_emu_hist->GetBinContent(triggerSF_emu_hist->GetXaxis()->FindBin(p_l2_trigg), triggerSF_emu_hist->GetYaxis()->FindBin(p_l1_trigg));
               }
               fEventWeight=*w_pu * *w_mc;     //Set event weight 
               SFWeight=*sf_lep1 * *sf_lep2 * *w_topPT * *w_bTag * triggerSF;     //Set combined SF weight
            }
            
            //Calculate MET before parton shower and hadronization for DM scenario and SUSY scenarios
            //And get number of (non)prompt neutrinos in event
            int NpromptNeutrinos=0;
            int NnonpromptNeutrinos=0;
            for (auto const &genParticle : *genParticles){
               if((!SUSY_scen && abs(genParticle.pdgId)>100000) || (SUSY_scen && abs(genParticle.pdgId)==1000022)){
                  neutrinoPair+=genParticle.p;
               }
               if(abs(genParticle.pdgId)==12 || abs(genParticle.pdgId)==14 || abs(genParticle.pdgId)==16){
                  if(genParticle.isPrompt) NpromptNeutrinos++;
                  else {
                     NnonpromptNeutrinos++;
                  }
               }
            }
            
            //Get DeltaPhi between MET (or genMet or neutrino pT) and nearest Lepton
            float dPhiMETnearLep=4;
            float dPhiMETnearLepPuppi=4;
            float dPhiMETnearLepXYcorr=4;
            float mt_MetNextLep=0; 
            float mt_NuNuNextLep=0; 
            for (TLorentzVector const lep : {p_l1,p_l2}){
               const float dPhi=MET->p.DeltaPhi(lep);
               const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(lep);
               const float dPhi_XYcorr=MET_XYcorr->p.DeltaPhi(lep);
               if (std::abs(dPhi) < std::abs(dPhiMETnearLep)) {
                  dPhiMETnearLep=dPhi;
                  mt_MetNextLep=phys::M_T(MET->p,lep);
               }
               if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearLepPuppi)) {
                  dPhiMETnearLepPuppi=dPhi_Puppi;
               }
               if (std::abs(dPhi_XYcorr) < std::abs(dPhiMETnearLepXYcorr)) {
                  dPhiMETnearLepXYcorr=dPhi_XYcorr;
               }
            }
            float dPhigenMETnearLep=4;
            for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
               const float dPhi=GENMET->p.DeltaPhi(lep);
               if (std::abs(dPhi) < std::abs(dPhigenMETnearLep)) dPhigenMETnearLep=dPhi;
            }
            float dPhiPtNunearLep=4;
            for (TLorentzVector const lep : {*genLepton,*genAntiLepton}) {
               const float dPhi=neutrinoPair.DeltaPhi(lep);
               if (std::abs(dPhi) < std::abs(dPhiPtNunearLep)) {
                  dPhiPtNunearLep=dPhi;
                  mt_NuNuNextLep=phys::M_T(neutrinoPair,lep);
               }
            }
            
            //Get DeltaPhi between MET nearest/leading/farest (b) Jet and HT
            float dPhiMETnearJet=4;
            float dPhiMETfarJet=0;
            float dPhiMETleadJet=4;
            float dPhiMETlead2Jet=4;
            float dPhiMetBJet=4;
            float dPhiMETnearJet_Puppi=4;
            float dPhiMETfarJet_Puppi=0;
            float dPhiMETleadJet_Puppi=4;
            float dPhiMETlead2Jet_Puppi=4;
            float dPhiMetBJet_Puppi=4;
            float dPhiLep1Lep2=4;
            float dPhiJet1Jet2=4;
            float dPhiLep1bJet=4;
            float dPhiLep1Jet1=4;
            TLorentzVector MHT(0.,0.,0.,0.);
            float HT=0;
            float ratio_pTj1_vecsum_pT_l1_l2_bjet=0;
            if(cjets.size()>1){
               if (rec_selection){ 
                  for (tree::Jet const &jet : cjets) {
                     const float dPhi=MET->p.DeltaPhi(jet.p);
                     const float dPhi_Puppi=MET_Puppi->p.DeltaPhi(jet.p);
                     if (std::abs(dPhi) < std::abs(dPhiMETnearJet)) dPhiMETnearJet=dPhi;
                     if (std::abs(dPhi) > std::abs(dPhiMETfarJet)) dPhiMETfarJet=dPhi;
                     if (std::abs(dPhi_Puppi) < std::abs(dPhiMETnearJet_Puppi)) dPhiMETnearJet_Puppi=dPhi_Puppi;
                     if (std::abs(dPhi_Puppi) > std::abs(dPhiMETfarJet_Puppi)) dPhiMETfarJet_Puppi=dPhi_Puppi;
                     HT+=jet.p.Pt();
                     MHT+=jet.p;
                  }
                  dPhiMETleadJet=MET->p.DeltaPhi(cjets[0].p);
                  dPhiMETlead2Jet=MET->p.DeltaPhi(cjets[1].p);
                  dPhiMETleadJet_Puppi=MET_Puppi->p.DeltaPhi(cjets[0].p);
                  dPhiMETlead2Jet_Puppi=MET_Puppi->p.DeltaPhi(cjets[1].p);
                  dPhiLep1Lep2=p_l1.DeltaPhi(p_l2);
                  dPhiJet1Jet2=cjets[0].p.DeltaPhi(cjets[1].p);
                  dPhiLep1Jet1=p_l1.DeltaPhi(cjets[0].p);
                  
                  if(BJets.size()>0){
                     dPhiMetBJet=MET->p.DeltaPhi(BJets[0].p);
                     dPhiMetBJet_Puppi=MET_Puppi->p.DeltaPhi(BJets[0].p);
                     ratio_pTj1_vecsum_pT_l1_l2_bjet=cjets[0].p.Pt()/(p_l1+p_l2+BJets[0].p).Pt();
                     dPhiLep1bJet=p_l1.DeltaPhi(BJets[0].p);
                  }
               }
            }
            
            //Further variables
            float mt_MetLep1=phys::M_T(MET->p,p_l1);
            float mt_MetLep2=phys::M_T(MET->p,p_l2);
            float conMt_Lep1Lep2=phys::conM_T(p_l1,p_l2);
            float dPhiNuNu=genNeutrino->DeltaPhi(*genAntiNeutrino);
            
            TLorentzVector leadGenLepton;
            if (genLepton->Pt()>genAntiLepton->Pt()) leadGenLepton=*genLepton;
            else leadGenLepton=*genAntiLepton;
            
            float mt_genMetLep1=phys::M_T(GENMET->p,leadGenLepton);
            float mt_genNeutrinosLep1=phys::M_T(neutrinoPair,leadGenLepton);
            
            //Sort genTop pT
            std::vector<TLorentzVector> gen_tops;
            gen_tops.push_back(*genTop);
            gen_tops.push_back(*genAntiTop);
            sort(gen_tops.begin(), gen_tops.end(), tree::PtGreaterLorentz);
            float pT_top1=0;
            float pT_top2=0;
            pT_top1=gen_tops[0].Pt();
            pT_top2=gen_tops[1].Pt();
            
            //Fill minimal tree for TTbar resolution used in binning/unfolding studies
            if (minimalTree){
               minTree_ee=*is_ee;
               minTree_mumu=*is_mumu;
               minTree_emu=*is_emu;
               minTree_MET=met;
               minTree_PtNuNu=neutrinoPair.Pt();
               minTree_PhiRec=abs(dPhiMETnearLep);
               minTree_PhiGen=abs(dPhigenMETnearLep);
               minTree_PhiNuNu=abs(dPhiPtNunearLep);
               minTree_PhiMetNearJet=abs(dPhiMETnearJet);
               minTree_PhiMetFarJet=abs(dPhiMETfarJet);
               minTree_PhiMetLeadJet=abs(dPhiMETleadJet);
               minTree_PhiMetLead2Jet=abs(dPhiMETlead2Jet);
               minTree_PhiMetbJet=abs(dPhiMetBJet);
               minTree_dPhiLep1Lep2=abs(dPhiLep1Lep2);
               minTree_dPhiJet1Jet2=abs(dPhiJet1Jet2);
               minTree_METsig=MET->sig;
               minTree_N=lumi_weight*fEventWeight;
               minTree_SF=SFWeight;
               minTree_runNo=*runNo;
               minTree_lumNo=*lumNo;
               minTree_evtNo=*evtNo;
               minTree_genDecayMode=*genDecayMode_pseudo;
               minTree_genMet=genMet;
               minTree_PuppiMet=met_puppi;
               minTree_XYcorrMet=MET_XYcorr->p.Pt();
               minTree_HT=HT;
               minTree_HT_phi=MHT.Phi();
               minTree_MHT=MHT.M();
               minTree_MT=mt_MetLep1;
               minTree_genMT=mt_genNeutrinosLep1;
               minTree_n_Interactions=*n_Interactions;
               minTree_n_Interactions_gen=*n_Interactions_gen;
               minTree_MT_nextLep=mt_MetNextLep;
               minTree_genMT_nextLep=mt_NuNuNextLep;
               minTree_leadTop=pT_top1;
               minTree_dPhiNuNu=abs(dPhiNuNu);
               minTree_PhiRecPuppi=abs(dPhiMETnearLepPuppi);
               minTree_PhiRecXYcorr=abs(dPhiMETnearLepXYcorr);
               minTree_looseLeptonVeto= *addLepton? 1: 0;
               minTree_nJets=cjets.size();
               minTree_PhiMetNearJet_Puppi=abs(dPhiMETnearJet_Puppi);
               minTree_PhiMetFarJet_Puppi=abs(dPhiMETfarJet_Puppi);
               minTree_PhiMetLeadJet_Puppi=abs(dPhiMETleadJet_Puppi);
               minTree_PhiMetLead2Jet_Puppi=abs(dPhiMETlead2Jet_Puppi);
               minTree_PhiMetbJet_Puppi=abs(dPhiMetBJet_Puppi);
               minTree_dPhiLep1bJet=abs(dPhiLep1bJet);
               minTree_dPhiLep1Jet1=abs(dPhiLep1Jet1);
               minTree_ratioMET_sqrtMETunc_Puppi=MET_Puppi->uncertainty/sqrt(met_puppi);
               minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet=ratio_pTj1_vecsum_pT_l1_l2_bjet;
               minTree_METunc_Puppi=MET_Puppi->uncertainty;
               minTree_METunc_PF=MET->uncertainty;
               minTree_absmetres_PUPPI=abs(met_puppi-genMet);
               minTree_PFMET_phi=MET->p.Phi();
               minTree_PuppiMET_phi=MET_Puppi->p.Phi();
               minTree_CaloMET=MET_Calo->p.Pt();
               minTree_CaloMET_phi=MET_Calo->p.Phi();
               minTree_genMET_phi=GENMET->p.Phi();
               minTree_PtNuNu_phi=neutrinoPair.Phi();
               minTree_NpromptNeutrinos=NpromptNeutrinos;
               minTree_NnonpromptNeutrinos=NnonpromptNeutrinos;
               minTree_MT2=*mt2;
               minTree_vecsum_pT_allJet=MHT.Pt();
               minTree_vecsum_pT_l1l2_allJet=(MHT+p_l1+p_l2).Pt();
               minTree_mass_l1l2_allJet=(MHT+p_l1+p_l2).M();
               minTree_ratio_vecsumpTlep_vecsumpTjet=(p_l1+p_l2).Pt()/MHT.Pt();
               
               minTree_DNN_MET_pT=-1;
               minTree_DNN_MET_phi=-4;
               minTree_DNN_MET_dPhi_nextLep=-1;
               
               if (rec_selection){
                  minTree_Lep1_pt=p_l1.Pt();
                  minTree_Lep1_phi=p_l1.Phi();
                  minTree_Lep1_eta=p_l1.Eta();
                  minTree_Lep1_E=p_l1.E();
                  minTree_Lep1_flavor=flavor_l1;
                  minTree_Lep2_pt=p_l2.Pt();
                  minTree_Lep2_phi=p_l2.Phi();
                  minTree_Lep2_eta=p_l2.Eta();
                  minTree_Lep2_E=p_l2.E();
                  minTree_Lep2_flavor=flavor_l2;
                  minTree_mLL=(p_l1+p_l2).M();
                  if(cjets.size()>1){
                     minTree_Jet1_pt=cjets[0].p.Pt();
                     minTree_Jet1_phi=cjets[0].p.Phi();
                     minTree_Jet1_eta=cjets[0].p.Eta();
                     minTree_Jet1_E=cjets[0].p.E();
                     minTree_Jet1_bTagScore=cjets[0].bTagDeepCSV;
                     minTree_Jet2_pt=cjets[1].p.Pt();
                     minTree_Jet2_phi=cjets[1].p.Phi();
                     minTree_Jet2_eta=cjets[1].p.Eta();
                     minTree_Jet2_E=cjets[1].p.E();
                     minTree_Jet2_bTagScore=cjets[1].bTagDeepCSV;
                     minTree_mjj=(cjets[0].p+cjets[1].p).M();
                  }
               }
               else{
                  minTree_MET=-1.;
                  minTree_PhiRec=-1.;
                  minTree_PhiMetNearJet=-1.;
                  minTree_PhiMetFarJet=-1.;
                  minTree_PhiMetLeadJet=-1.;
                  minTree_PhiMetLead2Jet=-1.;
                  minTree_PhiMetbJet=-1.;
                  minTree_PhiMetbJet=-1.;
                  minTree_METsig=-1.;
                  minTree_PuppiMet=-1.;
                  minTree_XYcorrMet=-1.;
                  minTree_MHT=-1.;
                  minTree_HT=-1.;
                  minTree_HT_phi=-4.;
                  minTree_MT=-1.;
                  minTree_MT_nextLep=-1.;
                  minTree_SF=0.;
                  minTree_PhiPtnunuMet=-1;
                  minTree_PhiRecPuppi=-1;
                  minTree_PhiRecXYcorr=-1;
                  minTree_n_Interactions=0;
                  minTree_looseLeptonVeto=5;
                  minTree_nJets=-1;
                  minTree_ratioMET_sqrtMETunc_Puppi=-20;
                  minTree_ratio_pTj1_vecsum_pT_l1_l2_bjet=-1;
                  minTree_METunc_PF=-1;
                  minTree_METunc_Puppi=-1;
                  minTree_absmetres_PUPPI=-1;
                  minTree_Lep1_pt=-1;
                  minTree_Lep1_phi=-5;
                  minTree_Lep1_eta=-3;
                  minTree_Lep1_E=-1;
                  minTree_Lep1_flavor=-1;
                  minTree_Lep2_pt=-1;
                  minTree_Lep2_phi=-5;
                  minTree_Lep2_eta=-3;
                  minTree_Lep2_E=-1;
                  minTree_Lep2_flavor=-1;
                  minTree_Jet1_pt=-1;
                  minTree_Jet1_phi=-5;
                  minTree_Jet1_eta=-3;
                  minTree_Jet1_E=-1;
                  minTree_Jet1_bTagScore=-5;
                  minTree_Jet1_unc=-1;
                  minTree_Jet2_pt=-1;
                  minTree_Jet2_phi=-5;
                  minTree_Jet2_eta=-3;
                  minTree_Jet2_E=-1;
                  minTree_Jet2_unc=-1;
                  minTree_CaloMET=-1;
                  minTree_CaloMET_phi=-5;
                  minTree_PuppiMET_phi=-5;
                  minTree_PFMET_phi=-5;
                  minTree_MT2=-1;
                  minTree_vecsum_pT_allJet=-1;
                  minTree_vecsum_pT_l1l2_allJet=-1;
                  minTree_mass_l1l2_allJet=-1;
                  minTree_ratio_vecsumpTlep_vecsumpTjet=-1;
                  minTree_mjj=-1;

               }
            }
            
            if (minimalTree) ttbar_res.Fill();
                              
         }// evt loop
         io::log<<"";
         file.Close();
         
         //Save ntuple for TTbar resolution used in binning studies (only if last file of sample e.g. for SingleTop)
         if (minimalTree && cfg.datasets.getDatasubsets({dss.datasetName}).back().name==dss.name) {
            ttbar_res_saver.save(ttbar_res,dss.datasetName);
            ttbar_res.Reset();
         }
         
         //For multi save dss name
         dssName_multi=TString(dss.datasetName);
      
      }  // datasubset loop
            
   } // dataset loop
   
}
