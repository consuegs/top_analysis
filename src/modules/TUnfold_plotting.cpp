//Script to plot the result of TUnfold_unfolding
#include <iostream>
#include <iomanip>
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
#include <TVectorD.h>
#include <TParameter.h>

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/tunfoldPlottingHelper.hpp"

using namespace std;
using namespace Systematic;
using namespace tunfoldplotting;

Config const &cfg=Config::get();

extern "C"
void run()
{
   TString sample = cfg.tunfold_InputSamples[0];
   TString sample_response = cfg.tunfold_ResponseSample;
   
   // Use pT reweighted
   bool withPTreweight = cfg.tunfold_withPTreweight;
   TString scale = cfg.tunfold_scalePTreweight;
   
   // Use phi reweighted
   bool withPHIreweight = cfg.tunfold_withPHIreweight;
   TString scale_phi = cfg.tunfold_scalePHIreweight;
   
   // Use DNN instead of pfMET
   bool withDNN = cfg.tunfold_withDNN;
   
   // Use pf instead of PuppiMET
   bool withPF = cfg.tunfold_withPF;
   
   // Use puppi instead of pfMET
   bool withPuppi = !withDNN && !withPF;
   
   // Use same bin numbers for gen/true
   bool withSameBins = cfg.tunfold_withSameBins;
   
   // include signal to pseudo data
   bool withBSM = cfg.tunfold_withBSM;
   TString scale_BSM = cfg.tunfold_scaleBSM;
   
   //Use scale factor
   bool withScaleFactor = cfg.tunfold_withScaleFactor;
   
   //Plot comparison
   bool plotComparison = cfg.tunfold_plotComparison;
   
   //Plot toy studies
   bool plotToyStudies = cfg.tunfold_plotToyStudies;
   
   // ~//Use alternative pseudo data (amcAtNLO)
   // ~bool useAltReco = true;
   bool useAltReco = false;
   
   //Plot only theory prediction (also MC)
   // ~bool onlyTheo = true;
   bool onlyTheo = false;
   
   //Plot fixed order theory prediction
   // ~bool plotTheo = true;
   bool plotTheo = false;
   
   //Use deltaR gen level cut
   // ~bool deltaRgen = true;
   bool deltaRgen = false;
   
   //Use deltaR gen level cut
   // ~bool jetVetoMaps = true;
   bool jetVetoMaps = false;
   
   //Use real data
   // ~bool useRealData = false;
   bool useRealData = true;
   
   //Use Single Top DS
   // ~bool useSingleTopDS = false;
   bool useSingleTopDS = true;
   
   //////////////////////////
   // Define Distributions //
   //////////////////////////
   
   std::vector<distrUnfold> vecDistr;
   // ~vecDistr.push_back({"2D_dPhi_pTnunu",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new40",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new40_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"pTnunu",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"dPhi",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.1f",false});
   // ~vecDistr.push_back({"pTnunu_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"dPhi_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.1f",false});
   // ~vecDistr.push_back({"pTll",0,400,";p_{T}^{ll} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   // ~vecDistr.push_back({"pTnunu_new",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"pTnunu_new_singleLast_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"pTnunu_new",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   
   vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.2f",false});
   vecDistr.push_back({"inclusive",0,1,";Signal Bin;d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});

   
   vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",false,true});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];1/#sigma d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)]","%.2f",false,true});
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",true,true});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_mergedBins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur2PhiBins14_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur2PhiBins16_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_diffDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true,true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_diffDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true,true});
   
   //////////////////////////////////
   // Set Systematic Uncertainties //
   //////////////////////////////////
   
   // Nominal (no envelopes, required for combination)
   std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_DCTR_UP","MATCH_DCTR_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","TWDS","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Nominal (with envelopes, required for individual results)
   // ~std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Nominal (with envelopes, required for individual results) no PDF
   // ~std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Split JES
   // ~std::vector<TString> systVec = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN"};
   // ~std::vector<TString> systVec = {"JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"CR1","CR2","ERDON"};
   // ~std::vector<TString> systVec = {"MATCH_UP","MATCH_DOWN"};
   // ~std::vector<TString> systVec = {"MATCH_DCTR_UP","MATCH_DCTR_DOWN"};
   // ~std::vector<TString> systVec = {"MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN",};
   // ~std::vector<TString> systVec = {"XSEC_ST_UP","XSEC_ST_DOWN","MESCALE_UP","MESCALE_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP"};
   // ~std::vector<TString> systVec = {"MESCALE_UP"};
   // ~std::vector<TString> systVec = {"MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP"};
   // ~std::vector<TString> systVec = {"TWDS"};
   // ~std::vector<TString> systVec = {"MTOP169p5","MTOP175p5"};
   // ~std::vector<TString> systVec = {"TWDS","XSEC_ST_UP","XSEC_ST_DOWN"};
   // ~std::vector<TString> systVec = {"JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systVec = {"BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN"};
   // ~std::vector<TString> systVec = {"BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN"};
   // ~std::vector<TString> systVec = {"JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"JETPILEUPID_DOWN"};
   // ~std::vector<TString> systVec = {"JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN"};
   // ~std::vector<TString> systVec = {"JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN"};
   // ~std::vector<TString> systVec = {"MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN"};
   // ~std::vector<TString> systVec = {"ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN"};
   // ~std::vector<TString> systVec = {"JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"};
   // ~std::vector<TString> systVec = {"MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN"};
   // ~std::vector<TString> systVec = {"JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"};
   // ~std::vector<TString> systVec = {"CR1","CR2","ERDON"};
   // ~std::vector<TString> systVec = {};
/*
   //Remove HEM unc. for all year except 2018
   auto itr =std::find(systVec.begin(), systVec.end(), "JESUserDefinedHEM1516_DOWN");
   if (itr != systVec.end() && cfg.year_int != 3){
      systVec.erase(itr);
   }
*/
   ////////////////////////////
   // Set Merged Systematics //
   ////////////////////////////
   // ~std::map<TString,std::vector<TString>> systCombinations = {};
   std::map<TString,std::vector<TString>> systCombinations = {
      {"JES",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JESUserDefinedHEM1516"}},
      {"JER",{"JEREta0","JEREta1"}},
      {"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      {"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      {"PS",{"PSISRSCALE","PSFSRSCALE"}},
      {"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
   };
   
   for (distrUnfold &dist : vecDistr){
   
      //==============================================
      // step 1 : open output file
      TString extraStoreFolder = "";
      if (deltaRgen) extraStoreFolder = "applyGenLevel_DeltaRcut";
      else if (jetVetoMaps) extraStoreFolder = "applyJetVetoMaps";
      else if (useRealData) extraStoreFolder = "realData";
      if (useSingleTopDS) extraStoreFolder += "_SingleTopDS";
      if (withBSM) extraStoreFolder += "_BSM_"+scale_BSM;
      io::RootFileSaver saver(TString::Format("TUnfold/plots%.1f.root",cfg.processFraction*100),TString::Format(!withScaleFactor ? "TUnfold_plotting%.1f/%s/%s" : "TUnfold_plotting_SF_%.1f/%s/%s",cfg.processFraction*100,extraStoreFolder.Data(),dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));
      
      // Saver for result used for combination of years
      io::RootFileSaver resultSaver(TString::Format("TUnfold/results%.1f.root",cfg.processFraction*100),TString::Format("%s/%s",extraStoreFolder.Data(),dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));

      //==============================================
      // step 2 : read binning schemes and input histograms
      TString sampleResult = sample;
      if (useRealData) sampleResult = "realData";
      TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
      TString input_loc_result=(useAltReco)? "TUnfold_results_TTbar_amcatnlo_"+sample_response : "TUnfold_results_"+sampleResult+"_"+sample_response;
      TString saveName=(useAltReco)? "TTbar_amcatnlo_"+sample_response : sampleResult+"_"+sample_response;
      TString saveName2D=(useAltReco)? "correlations/TTbar_amcatnlo_"+sample_response : "correlations/"+sampleResult+"_"+sample_response;
      if (useSingleTopDS) {
         input_loc_result+="_SingleTopDS";
         saveName+="_SingleTopDS";
         saveName2D+="_SingleTopDS";
      }
      if (withBSM) {
         TString addString = "_BSM";
         if (std::stof(scale_BSM.Data()) != 1.) addString = "_BSM"+scale_BSM;
         input_loc+=addString;
         input_loc_result+=addString;
         saveName+=addString;
         saveName2D+=addString;
      }
      if (withPF) {
         input_loc+="_PF";
         input_loc_result+="_PF";
         saveName+="_PF";
         saveName2D+="_PF";
      }
      if (withPuppi) {
         input_loc+="_Puppi";
         input_loc_result+="_Puppi";
         saveName+="_Puppi";
         saveName2D+="_Puppi";
      }
      if (withDNN) {
         input_loc+="_DNN";
         input_loc_result+="_DNN";
         saveName+="_DNN";
         saveName2D+="_DNN";
      }
      if (withSameBins) {
         input_loc+="_SameBins";
         input_loc_result+="_SameBins";
         saveName+="_SameBins";
         saveName2D+="_SameBins";
      }
      if (withPTreweight) {
         input_loc+="_PTreweight"+scale;
         input_loc_result+="_PTreweight"+scale;
         saveName+="_PTreweight"+scale;
         saveName2D+="_PTreweight"+scale;
      }
      else if (withPHIreweight) {
         input_loc+="_PHIreweight"+scale_phi;
         input_loc_result+="_PHIreweight"+scale_phi;
         saveName+="_PHIreweight"+scale_phi;
         saveName2D+="_PHIreweight"+scale_phi;
      }
      if (plotComparison) saveName="compareMethods/"+saveName;
      if (onlyTheo) saveName="theoryComparison/"+saveName;
      
      TString inputName = "Nominal";
      if (deltaRgen) inputName = "applyGenLevel_DeltaRcut";
      else if (jetVetoMaps) inputName = "applyJetVetoMaps";
      io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",input_loc.Data(),inputName.Data(),cfg.processFraction*100));
      
      TString input_loc_old = input_loc;
      input_loc += "/"+dist.varName;
      input_loc_result += "/"+dist.varName;
      
      // Print total cross section
      if (dist.varName == "inclusive") {
         TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
         TH1F *realDis=histReader.read<TH1F>((useAltReco)? input_loc+"/histDataTruthAlt" : input_loc+"/histDataTruth");
         TH1F *realDisAlt=histReader.read<TH1F>(input_loc+"/histDataTruthAlt");
         TH1F *realDisHerwig=histReader.read<TH1F>(input_loc+"/histDataTruthHerwig");
         TH2D *covMatrix_bbb=histReader.read<TH2D>(input_loc_result+"/BBB/cov_output");
         
         unfolded_bbb->Scale(1./cfg.lumi);
         realDis->Scale(1./cfg.lumi);
         realDisAlt->Scale(1./cfg.lumi);
         realDisHerwig->Scale(1./cfg.lumi);
         
         std::map<TString,TH1F> indShifts_bbb;
         std::map<TString,TH2F> indResponse_bbb;
         std::map<TString,TH2F> indResponseAlt_bbb;
         std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb,indResponse_bbb,indResponseAlt_bbb);         
         // ~std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,false,withScaleFactor,indShifts_bbb,indResponse_bbb);         
         
         std::cout<<"----------------------------------------------------------"<<std::endl;
         std::cout<<"Measured total cross section:"<<unfolded_bbb->GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section:"<<realDis->GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section (amc@NLO):"<<realDisAlt->GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section (Herwig):"<<realDisHerwig->GetBinContent(1)<<std::endl;
         plot_systBreakdown(indShifts_bbb,&saver,"","BBB",dist.varName,useRealData,systCombinations);
         std::cout<<"----------------------------------------------------------"<<std::endl;
         
         // Scale back to get event count for combination
         unfolded_bbb->Scale(cfg.lumi);
         realDis->Scale(cfg.lumi);
         realDisAlt->Scale(cfg.lumi);
         realDisHerwig->Scale(cfg.lumi);
         
         if (!(withPTreweight || withPHIreweight)){
            resultSaver.save(*realDis,"Truth");
            resultSaver.save(*realDisAlt,"TruthAlt");
            resultSaver.save(*realDisHerwig,"TruthHerwig");
            resultSaver.save(*unfolded_bbb,"BBB/unfolded");
            resultSaver.save(*covMatrix_bbb,"BBB/cov");
            
            for (const auto &[key, value]: indShifts_bbb){
               resultSaver.save(indShifts_bbb[key],"BBB/unfolded_"+key);
               
               // Store shifted MC predictions
               if(std::find(Systematic::ttbarTypes.begin(), Systematic::ttbarTypes.end(), Systematic::convertType(key)) != Systematic::ttbarTypes.end()){
                  if (Systematic::convertType(key) != Systematic::pdf_envelope){
                  TH1D temp_realDis_syst;
                  TH1D temp_realDisAlt_syst;
                  TH1D temp_realDisHerwig_syst;
                  temp_realDis_syst = *indResponse_bbb[key].ProjectionX();
                  temp_realDisAlt_syst = *indResponseAlt_bbb[key].ProjectionX(); 
                  temp_realDisHerwig_syst = *indResponseAlt_bbb[key].ProjectionX(); 
                     if(dist.norm) {   // Strange scaling, but this matches the one used for the nominal true hist and is fixed in combineYears
                        temp_realDis_syst.Scale(1./temp_realDis_syst.Integral());
                        temp_realDisAlt_syst.Scale(1./temp_realDis_syst.Integral());
                        temp_realDisHerwig_syst.Scale(1./temp_realDis_syst.Integral());
                        temp_realDis_syst.Scale(cfg.lumi);
                        temp_realDisAlt_syst.Scale(cfg.lumi);
                        temp_realDisHerwig_syst.Scale(cfg.lumi);
                     }
                  resultSaver.save(temp_realDis_syst,"Truth_"+key);
                  resultSaver.save(temp_realDisAlt_syst,"TruthAlt_"+key);
                  resultSaver.save(temp_realDisHerwig_syst,"TruthHerwig_"+key);
                  }
                  else {   //shift for pdf envelope has to be derived separately since it does not follow logic as other uncertainties (enevelope done before year combination)
                     TH1F temp_realDis_syst = getPDFenvelope(input_loc+"/histMCGenRec_projX",input_loc_old,realDis,Systematic::convertVariation(key) == Systematic::up,false,false);
                     resultSaver.save(temp_realDis_syst,"Truth_"+key);
                  }
               }
            }            

         }
         continue;
      }

      TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
      TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

      if((!detectorBinning)||(!generatorBinning)) {
         cout<<"problem to read binning schemes\n";
      }

      // read histograms
      TH1F *realDis=histReader.read<TH1F>(input_loc+"/histDataTruth");
      TH1F *realDisAlt=histReader.read<TH1F>((withBSM)? ((useAltReco)? input_loc+"/histDataTruthAltBSM" : input_loc+"/histDataTruthBSM") : input_loc+"/histDataTruthAlt");
      TH1F *realDisHerwig=histReader.read<TH1F>((withBSM)? ((useAltReco)? input_loc+"/histDataTruthHerwigBSM" : input_loc+"/histDataTruthBSM") : input_loc+"/histDataTruthHerwig");
      TH1F *realDis_response=histReader.read<TH1F>(input_loc+"/histMCGenRec_projX");
      TH1F *unfolded=histReader.read<TH1F>(input_loc_result+"/hist_unfoldedResult");
      TH1F *unfolded_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_unfoldedResult");
      TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
      TH2F *corrMatrix=histReader.read<TH2F>(input_loc_result+"/corr_matrix");
      TH2F *corrMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/corr_matrix");
      TH2D *covMatrix=histReader.read<TH2D>(input_loc_result+"/cov_output");
      TH2D *covMatrix_reg=histReader.read<TH2D>(input_loc_result+"/reg/cov_output");
      TH2D *covMatrix_bbb=histReader.read<TH2D>(input_loc_result+"/BBB/cov_output");
      TH2F *responseMatrix=histReader.read<TH2F>(input_loc+"/histMCGenRec_sameBins");
      TH2F *responseMatrix_fine=histReader.read<TH2F>(input_loc+"/histMCGenRec");
      TH1F *purity=histReader.read<TH1F>(input_loc+"/hist_purity");
      TH1F *stability=histReader.read<TH1F>(input_loc+"/hist_stability");
      TH1F *efficiency=histReader.read<TH1F>(input_loc+"/hist_efficiency");
      TParameter<float> *tau_par = histReader.read<TParameter<float>>(input_loc_result+"/reg/tau");
      
      if((!realDis)||(!realDis_response)||(!unfolded)) {
         cout<<"problem to read input histograms\n";
      }
      
      // divide by lumi to get xsec or by integral to get normalized xsec
      if(dist.norm){
         realDis->Scale(1./realDis->Integral());
         realDisAlt->Scale(1./realDisAlt->Integral());
         realDisHerwig->Scale(1./realDisHerwig->Integral());
         realDis_response->Scale(1./realDis_response->Integral());
         unfolded->Scale(1./unfolded->Integral());
         unfolded_reg->Scale(1./unfolded_reg->Integral());
         unfolded_bbb->Scale(1./unfolded_bbb->Integral());
      }
      else{
         realDis->Scale(1./cfg.lumi);
         realDisAlt->Scale(1./cfg.lumi);
         realDisHerwig->Scale(1./cfg.lumi);
         realDis_response->Scale(1./cfg.lumi);
         unfolded->Scale(1./cfg.lumi);
         unfolded_reg->Scale(1./cfg.lumi);
         unfolded_bbb->Scale(1./cfg.lumi);
      }
      
      // Get syst unc
      std::map<TString,TH1F> indShifts;
      std::map<TString,TH1F> indShifts_reg;
      std::map<TString,TH1F> indShifts_bbb;
      std::map<TString,TH2F> indResponse;
      std::map<TString,TH2F> indResponse_reg;
      std::map<TString,TH2F> indResponse_bbb;
      std::map<TString,TH2F> indResponseAlt;
      std::map<TString,TH2F> indResponseAlt_reg;
      std::map<TString,TH2F> indResponseAlt_bbb;
      std::pair<TH1F*,TH1F*> unfolded_total = getSystUnc(unfolded,input_loc_result+"/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts,indResponse,indResponseAlt,dist.norm);
      std::pair<TH1F*,TH1F*> unfolded_reg_total = getSystUnc(unfolded_reg,input_loc_result+"/reg/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_reg,indResponse_reg,indResponseAlt_reg,dist.norm);
      std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb,indResponse_bbb,indResponseAlt_bbb,dist.norm);
      std::pair<TH1F*,TH1F*> unfolded_dummy_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb,indResponse_bbb,indResponseAlt_bbb,dist.norm); // mock unfolded_bbb_total used for MC uncertainty since not derived for single year
      
      TH1F* unfoldedClone = (TH1F*)unfolded->Clone();    // needed for storing results for combination
      TH1F* unfoldedClone_reg = (TH1F*)unfolded_reg->Clone();    // needed for storing results for combination
      TH1F* unfoldedClone_bbb = (TH1F*)unfolded_bbb->Clone();    // needed for storing results for combination
      TH1F* realDisClone = (TH1F*)realDis->Clone();    // needed for storing results for combination
      TH1F* realDisCloneAlt = (TH1F*)realDisAlt->Clone();    // needed for storing results for combination
      TH1F* realDisCloneHerwig = (TH1F*)realDisHerwig->Clone();    // needed for storing results for combination
      TH2D* covMatrixClone = (TH2D*)covMatrix->Clone();    // needed for storing results for combination
      TH2D* covMatrixClone_reg = (TH2D*)covMatrix_reg->Clone();    // needed for storing results for combination
      TH2D* covMatrixClone_bbb = (TH2D*)covMatrix_bbb->Clone();    // needed for storing results for combination
      
      int num_bins;
      std::ofstream chi2_file;
      chi2_file.open ("chi2_values_data.txt");
      
      // Plot results
      std::vector<double> xbins_vec = plot_UnfoldedResult(generatorBinning,unfolded,unfolded_reg,unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,unfolded_dummy_total,unfolded_dummy_total,unfolded_dummy_total,tau_par->GetVal(),realDis,(withPTreweight || withPHIreweight)? realDis_response : realDisAlt,realDisHerwig,covMatrix,dist,plotComparison,saveName,&saver,num_bins,(withPTreweight || withPHIreweight),onlyTheo,plotTheo,chi2_file);
      
      chi2_file.close();
      
      Double_t* xbins = &xbins_vec[0];    //Convert binning vector to array used for setting binnings of final results
      
      
      //Plot response matrix
      plot_response(responseMatrix,saveName,&saver,dist.is2D);
      plot_correlation(corrMatrix,saveName2D,&saver);
      plot_correlation(corrMatrix_reg,saveName2D+"_reg",&saver);
      
      // Plot purity, stability and efficiency
      if (!dist.is2D){
         purity->SetBins(num_bins,xbins);
         stability->SetBins(num_bins,xbins);
         efficiency->SetBins(num_bins,xbins);
      }
      plot_pur_stab_eff(purity,stability,efficiency,saveName,unfolded->GetXaxis()->GetTitle(),dist,generatorBinning,&saver);
      
      //Plot syst. unc. breakdown
      if (!dist.is2D){ // set correct binning for plotting syst breakdown (currently only 1D)
         for(std::pair<TString,TH1F> const &shift : indShifts){
            indShifts[shift.first].SetBins(num_bins,xbins);
            indShifts_reg[shift.first].SetBins(num_bins,xbins);
            indShifts_bbb[shift.first].SetBins(num_bins,xbins);
         }
      }
      unfolded->SetTitle(dist.title);
      plot_systBreakdown(indShifts,&saver,saveName,"Nominal",unfolded->GetXaxis()->GetTitle(),useRealData,systCombinations);
      plot_systBreakdown(indShifts_reg,&saver,saveName,"Reg",unfolded->GetXaxis()->GetTitle(),useRealData,systCombinations);
      plot_systBreakdown(indShifts_bbb,&saver,saveName,"BBB",unfolded->GetXaxis()->GetTitle(),useRealData,systCombinations);
      
      //Plot diff in response for syst.
      // ~plot_response_diff(indResponse,get_response(responseMatrix),&saver,saveName);
      // ~plot_response_diff(indResponse,responseMatrix_fine,&saver,saveName);
      
      //Plot toy studies
      // ~if (plotComparison) saveName=sample+"_"+sample_response;
      if (plotToyStudies) {
         TProfile *profPull=histReader.read<TProfile>(input_loc_result+"/prof_pull");
         TProfile *profPull_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_pull");
         TProfile *profPull_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_pull");
         TProfile *profRes=histReader.read<TProfile>(input_loc_result+"/prof_res");
         TProfile *profRes_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_res");
         TProfile *profRes_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_res");
         TH1F *histPull=histReader.read<TH1F>(input_loc_result+"/hist_pull");
         TH1F *histPull_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_pull");
         TH1F *histPull_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_pull");
         TH1F *histRes=histReader.read<TH1F>(input_loc_result+"/hist_res");
         TH1F *histRes_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_res");
         TH1F *histRes_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_res");
         TH1F *histChi=histReader.read<TH1F>(input_loc_result+"/hist_chi");
         TH1F *histChi_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_chi");
         TH1F *histChi_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_chi");
         for (TString type : {"Pull","Residual"}) {   //Plot bin-wise pull and residual
            TCanvas canToy;
            canToy.cd();
            
            //Set temp profiles/hists
            TProfile* temp=0;
            TProfile* temp_reg=0;
            TProfile* temp_bbb=0;
            if (type=="Pull") {
               temp=profPull;
               temp_reg=profPull_reg;
               temp_bbb=profPull_bbb;
            }
            else {
               temp=profRes;
               temp_reg=profRes_reg;
               temp_bbb=profRes_bbb;
            }
            
            //Get RMS from profile
            TH1D tempRMS=*(temp->ProjectionX());
            TH1D tempRMS_reg=*(TH1D*)tempRMS.Clone();
            TH1D tempRMS_bbb=*(TH1D*)tempRMS.Clone();
            for (int i=1; i<=profPull->GetNbinsX(); i++){
               tempRMS.SetBinContent(i,temp->GetBinError(i)*sqrt(temp->GetBinEntries(i)));
               tempRMS_reg.SetBinContent(i,temp_reg->GetBinError(i)*sqrt(temp_reg->GetBinEntries(i)));
               tempRMS_bbb.SetBinContent(i,temp_bbb->GetBinError(i)*sqrt(temp_bbb->GetBinEntries(i)));
            }
            
            temp->SetStats(0);
            temp->GetXaxis()->SetTitle("Bin");
            if (type=="Pull") temp->SetMaximum(1.2);
            else temp->SetMaximum(1.5);
            
            temp->SetMarkerColor(kBlack);
            temp_reg->SetMarkerColor(kGreen+2);
            temp_bbb->SetMarkerColor(kViolet);
            tempRMS.SetMarkerColor(kBlack);
            tempRMS_reg.SetMarkerColor(kGreen+2);
            tempRMS_bbb.SetMarkerColor(kViolet);
            
            temp->SetLineWidth(0);
            temp_reg->SetLineWidth(0);
            temp_bbb->SetLineWidth(0);
            tempRMS.SetLineWidth(0);
            tempRMS_reg.SetLineWidth(0);
            tempRMS_bbb.SetLineWidth(0);
            
            temp->SetMarkerSize(1.3);
            temp_reg->SetMarkerSize(1.3);
            temp_bbb->SetMarkerSize(1.3);
            tempRMS.SetMarkerSize(1.3);
            tempRMS_reg.SetMarkerSize(1.3);
            tempRMS_bbb.SetMarkerSize(1.3);
            
            tempRMS.SetMarkerStyle(4);
            tempRMS_reg.SetMarkerStyle(4);
            tempRMS_bbb.SetMarkerStyle(4);
            
            temp->Draw("");
            temp_reg->Draw("same");
            temp_bbb->Draw("same");
            tempRMS.Draw("same");
            tempRMS_reg.Draw("same");
            tempRMS_bbb.Draw("same");
            
            // ~aline->DrawLine(0.5,0,num_bins+0.5,0);
            // ~if (type=="Pull") aline->DrawLine(0.5,1,num_bins+0.5,1);
            
            gfx::LegendEntries legE_pull;
            legE_pull.append(*temp,"Nominal "+type+" Mean","p");
            legE_pull.append(tempRMS,"Nominal "+type+" RMS","p");
            legE_pull.append(*temp_reg,"Regularized "+type+" Mean","p");
            legE_pull.append(tempRMS_reg,"Regularized "+type+" RMS","p");
            legE_pull.append(*temp_bbb,"BBB "+type+" Mean","p");
            legE_pull.append(tempRMS_bbb,"BBB "+type+" RMS","p");
            TLegend leg_pull=legE_pull.buildLegend(.2,.45,0.8,.65,2);
            leg_pull.SetTextSize(0.04);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/prof"+type,true,true,true);
         }
         
         for (TString type : {"Pull","Residual","Chi"}) {
            TCanvas canToy;
            canToy.cd();
            
            //Set temp profiles/hists
            TH1* temp=0;
            TH1* temp_reg=0;
            TH1* temp_bbb=0;
            if (type=="Pull") {
               temp=histPull;
               temp_reg=histPull_reg;
               temp_bbb=histPull_bbb;
            }
            else if (type=="Residual") {
               temp=histRes;
               temp_reg=histRes_reg;
               temp_bbb=histRes_bbb;
            }
            else {
               temp=histChi;
               temp_reg=histChi_reg;
               temp_bbb=histChi_bbb;
            }
            
            //Rebin
            // ~temp->Rebin(2);
            // ~temp_reg->Rebin(2);
            // ~temp_bbb->Rebin(2);
            
            //Fit expected distribution
            TF1* fitFunction;
            TFitResultPtr fitResult;
            if (type=="Pull") {
               fitResult = temp->Fit("gaus");
               fitFunction = temp->GetFunction("gaus");
            }
            else if(type=="Chi") {
               fitFunction = new TF1("chi2_func", "[0]*ROOT::Math::chisquared_pdf(x,[1])",0,30);
               fitFunction->SetParameter(0, temp->Integral());
               fitFunction->SetParameter(1, temp->GetMean());
               fitResult = temp->Fit(fitFunction);
               // ~fitFunction = temp->GetFunction("chi2_func");
            }
            
            temp->SetLineColor(kBlack);
            temp_reg->SetLineColor(kGreen+2);
            temp_bbb->SetLineColor(kViolet);
            temp->SetLineWidth(2);
            temp_reg->SetLineWidth(2);
            temp_bbb->SetLineWidth(2);
            temp_bbb->SetMaximum(1.3*temp_bbb->GetMaximum());
            
            temp_bbb->SetStats(0);
            temp_bbb->Draw("hist");
            temp->Draw("hist same");
            temp_reg->Draw("hist same");

            gfx::LegendEntries legE_pull;
            if(type=="Chi") {
               legE_pull.append(*temp,TString::Format("Nominal (#mu=%.1f)",temp->GetMean()),"l");
               legE_pull.append(*temp_reg,TString::Format("Regularized (#mu=%.1f)",temp_reg->GetMean()),"l");
               legE_pull.append(*temp_bbb,TString::Format("BBB (#mu=%.1f)",temp_bbb->GetMean()),"l");
            }
            else{
               legE_pull.append(*temp,TString::Format("Nominal (#mu=%.1f #sigma=%.1f)",temp->GetMean(),temp->GetRMS()),"l");
               legE_pull.append(*temp_reg,TString::Format("Regularized (#mu=%.1f #sigma=%.1f)",temp_reg->GetMean(),temp_reg->GetRMS()),"l");
               legE_pull.append(*temp_bbb,TString::Format("BBB (#mu=%.1f #sigma=%.1f)",temp_bbb->GetMean(),temp_bbb->GetRMS()),"l");
            }
            
            if(fitFunction) {
               fitFunction->SetLineColor(kBlue);
               fitFunction->Draw("same");
               if(type=="Chi") {
                  legE_pull.append(*fitFunction,TString::Format("Fit (#mu=%.1f)",fitFunction->GetParameter(1)),"l");
               }
               else{
                  legE_pull.append(*fitFunction,TString::Format("Fit (#mu=%.1f #sigma=%.1f)",fitFunction->GetParameter(1),fitFunction->GetParameter(2)),"l");
               }
            }
            
            TLegend leg_pull=legE_pull.buildLegend(.59,.75,0.93,0.94,1);
            leg_pull.SetTextSize(0.03);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/hist"+type,true,true,true);
         }
            
      }
      
      // Store results for combination of years
      
      if(!(withPTreweight || withPHIreweight)){
         //Scale with lumi to get event count
         realDisClone->Scale(cfg.lumi);
         realDisCloneAlt->Scale(cfg.lumi);
         realDisCloneHerwig->Scale(cfg.lumi);
         
         if (!dist.is2D){  // for 2D Signal bins should be used
            unfoldedClone->SetBins(num_bins,xbins);
            unfoldedClone_reg->SetBins(num_bins,xbins);
            unfoldedClone_bbb->SetBins(num_bins,xbins);
            // ~realDisClone->SetBins(num_bins,xbins);
            // ~realDisCloneAlt->SetBins(num_bins,xbins);
         }
         
         resultSaver.save(*realDisClone,"Truth");
         resultSaver.save(*realDisCloneAlt,"TruthAlt");
         resultSaver.save(*realDisCloneHerwig,"TruthHerwig");
         resultSaver.save(*responseMatrix,"ResponseMatrix");
         resultSaver.save(*responseMatrix_fine,"ResponseMatrix_fine");
         
         //Scale with lumi to get event count
         unfoldedClone->Scale(cfg.lumi);
         unfoldedClone_reg->Scale(cfg.lumi);
         unfoldedClone_bbb->Scale(cfg.lumi);
               
         resultSaver.save(*unfoldedClone,"NoReg/unfolded");
         resultSaver.save(*unfoldedClone_reg,"Reg/unfolded");
         resultSaver.save(*unfoldedClone_bbb,"BBB/unfolded");
         
         resultSaver.save(*covMatrixClone,"NoReg/cov");
         resultSaver.save(*covMatrixClone_reg,"Reg/cov");
         resultSaver.save(*covMatrixClone_reg,"BBB/cov");
         
         for (const auto &[key, value]: indShifts){
            if (!dist.is2D){
               indShifts_reg[key].SetBins(num_bins,xbins);
               indShifts_bbb[key].SetBins(num_bins,xbins);
            }
            
            // Store shifted unfolded results
            resultSaver.save(indShifts[key],"NoReg/unfolded_"+key);
            resultSaver.save(indShifts_reg[key],"Reg/unfolded_"+key);
            resultSaver.save(indShifts_bbb[key],"BBB/unfolded_"+key);
            
            // Store shifted MC predictions
            if(std::find(Systematic::ttbarTypes.begin(), Systematic::ttbarTypes.end(), Systematic::convertType(key)) != Systematic::ttbarTypes.end()){
               if (Systematic::convertType(key) != Systematic::pdf_envelope){
                  TH1D temp_realDis_syst = *indResponse[key].ProjectionX();
                  TH1D temp_realDisAlt_syst = *indResponseAlt[key].ProjectionX();
                  TH1D temp_realDisHerwig_syst = *indResponseAlt[key].ProjectionX();
                  if(dist.norm) {   // Strange scaling, but this matches the one used for the nominal true hist and is fixed in combineYears
                     temp_realDis_syst.Scale(1./temp_realDis_syst.Integral());
                     temp_realDisAlt_syst.Scale(1./temp_realDis_syst.Integral());
                     temp_realDisHerwig_syst.Scale(1./temp_realDis_syst.Integral());
                     temp_realDis_syst.Scale(cfg.lumi);
                     temp_realDisAlt_syst.Scale(cfg.lumi);
                     temp_realDisHerwig_syst.Scale(cfg.lumi);
                  }
                  resultSaver.save(temp_realDis_syst,"Truth_"+key);
                  resultSaver.save(temp_realDisAlt_syst,"TruthAlt_"+key);
                  resultSaver.save(temp_realDisHerwig_syst,"TruthHerwig_"+key);
               }
               else {   //shift for pdf envelope has to be derived separately since it does not follow logic as other uncertainties (enevelope done before year combination)
                  TH1F temp_realDis_syst = getPDFenvelope(input_loc+"/histMCGenRec_projX",input_loc_old,realDisClone,Systematic::convertVariation(key) == Systematic::up,dist.norm,false);
                  resultSaver.save(temp_realDis_syst,"Truth_"+key);
               }
            }
         }
         
         resultSaver.save(*generatorBinning,"generatorBinning");
      }
   }
}
