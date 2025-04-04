//Script to plot and combine and plot the unfolded results for all four periods
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
#include "tools/util.hpp"

using namespace std;
using namespace Systematic;
using namespace tunfoldplotting;

Config &cfg=Config::get();

extern "C"
void run()
{
   cfg.setLumi(137650.);   // only used for plotting reasons
   
   bool jesComparison = false;   // produces plots comparing split and regrouped JES
   // ~bool jesComparison = true;   // produces plots comparing split and regrouped JES
   
   //Plot fixed order theory prediction
   bool plotTheo = true;
   // ~bool plotTheo = false;
   
   // include signal to pseudo data
   bool withBSM = cfg.tunfold_withBSM;
   TString scale_BSM = cfg.tunfold_scaleBSM;
   
   //Use real data
   // ~bool useRealData = false;
   bool useRealData = true;  // set this flag to false for Figure 6 of paper, true for Figures 4/7/8

   //Use Single Top DS
   // ~bool useSingleTopDS = false;
   bool useSingleTopDS = true;
   
   //Set CR uncertainty in last 2D bin based on result from merged 2D distribution (MC stats to low)
   bool fixCRunc = true;
   // ~bool fixCRunc = false;
   float CR_up_lastBin = 0.;
   float CR_down_lastBin = -0.0382661;
      
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
   
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.2f",false});
   // ~vecDistr.push_back({"inclusive",0,1,";Signal Bin;d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",true,true});
   vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",false,true});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];1/#sigma d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)]","%.2f",false,true});
   
   // ~vecDistr.push_back({"pTll",0,400,";p_{T}^{ll} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_diffDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true,true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur9Bins_diffDet_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true,true});
   
   //////////////////////////////////
   // Set Systematic Uncertainties //
   //////////////////////////////////
   
   // Nominal
   std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_DCTR_UP","MATCH_DCTR_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","TWDS","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   //Experimental unc
   // ~std::vector<TString> systVec = {"BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","TRIG_UP","TRIG_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN"};
   
   //Theory unc
    // ~std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","MATCH_DCTR_UP","MATCH_DCTR_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","TOP_PT","TWDS","UETUNE_UP","UETUNE_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Split JES
   // ~std::vector<TString> systVec = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   //"five most relevant"
    // ~std::vector<TString> systVec = {"BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systVec = {"TRIG_UP","TRIG_DOWN"};
   // ~std::vector<TString> systVec = {"TWDS"};
   // ~std::vector<TString> systVec = {"JEREta1_UP","JEREta1_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","CR1","CR2","ERDON"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"MATCH_UP","MATCH_DOWN"};
   // ~std::vector<TString> systVec = {"MATCH_DCTR_UP","MATCH_DCTR_DOWN"};
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN","XSEC_ST_UP","XSEC_ST_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systVec = {"PSFSRSCALE_UP","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systVec = {"ELECTRON_ID_UP","ELECTRON_ID_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsoluteYear_UP","JESAbsoluteYear_DOWN"};
   // ~std::vector<TString> systVec = {"JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsoluteYear_UP"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"};
   // ~std::vector<TString> systVec = {"JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   // ~std::vector<TString> systVec = {"ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN"};
   // ~std::vector<TString> systVec = {}; // for Figure 6 of paper systVec should be kept empty
   
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"JES",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JESUserDefinedHEM1516"}},
      // ~{"JER",{"JEREta0","JEREta1"}},
      // ~{"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      // ~{"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      // ~{"PS",{"PSISRSCALE","PSFSRSCALE"}},
      // ~{"PDF",{"PDF_ALPHAS","PDF_ENVELOPE"}},
      // ~{"BKG XSEC",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
   // ~};
   std::map<TString,std::vector<TString>> systCombinations = {
      {"JET",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JEREta0","JEREta1"}},
      {"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      {"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      {"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
      {"xOTHER EXP",{"L1PREFIRING","PU","TRIG","UNCLUSTERED","JETPILEUPID"}},
      {"ME_PS",{"MATCH_DCTR","MESCALE_ENVELOPE","PSISRSCALE","PSFSRSCALE","TOP_PT"}},
      {"xOTHER THEO",{"BSEMILEP","MTOP","PDF_ALPHAS","PDF_ENVELOPE","TOP_PT"}},
   };
   
   // add ingredients for CR and ME envelope if required
   if (std::find(systVec.begin(), systVec.end(), "CR_ENVELOPE_UP") != systVec.end() || std::find(systVec.begin(), systVec.end(), "CR_ENVELOPE_DOWN") != systVec.end()){
      systVec = util::addVectors(systVec,{"CR1","CR2","ERDON"});
   }
   if (std::find(systVec.begin(), systVec.end(), "MESCALE_ENVELOPE_UP") != systVec.end() || std::find(systVec.begin(), systVec.end(), "MESCALE_ENVELOPE_DOWN") != systVec.end()){
      systVec = util::addVectors(systVec,{"MESCALE_UP","MESCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN"});
   }
   if (std::find(systVec.begin(), systVec.end(), "MTOP_UP") != systVec.end() || std::find(systVec.begin(), systVec.end(), "MTOP_DOWN") != systVec.end()){
      systVec = util::addVectors(systVec,{"MTOP175p5","MTOP169p5"});
   }
   
   // Changes for jesComparison
   if(jesComparison){
      
      systVec = {"JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN",
      
      "JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN"
      };
      
      systCombinations = {
         {"Absolute",{"JESAbsoluteMPFBias","JESAbsoluteScale","JESFragmentation","JESPileUpDataMC","JESPileUpPtRef","JESRelativeFSR","JESSinglePionECAL","JESSinglePionHCAL"}},
         {"AbsoluteYear",{"JESAbsoluteStat","JESRelativeStatFSR","JESTimePtEta"}},
         {"BBEC1",{"JESPileUpPtBB","JESPileUpPtEC1","JESRelativePtBB"}},
         {"BBEC1Year",{"JESRelativeJEREC1","JESRelativePtEC1","JESRelativeStatEC"}},
         {"JES Split",{"Absolute","AbsoluteYear","BBEC1","BBEC1Year","JESFlavorRealistic","JESRelativeBal","JESRelativeSample"}},
   
         {"JES Regrouped",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1"}}
      };
   }
   
   std::vector<TString> systVec_realDis;
   
   std::vector<Config> cfg_years = {Config("2016_preVFP"),Config("2016_postVFP"),Config("2017"),Config("2018")};
   
   TString resultFolder = "";    //Define input and starge folder
   if (useRealData) resultFolder = "realData";
   if (useSingleTopDS) resultFolder += "_SingleTopDS";
   if (withBSM) resultFolder += "_BSM_"+scale_BSM;
   
   //derive version string for all years
   TString versionString = "";
   for (int i=0; i<4; i++)versionString += cfg_years[i].treeVersion.Data();
   
   std::ofstream chi2_file; //file to store chi2 values
   chi2_file.open(TString::Format("%s/../../../Combined/%s/TUnfold/chi2_values_data%s.txt",cfg.outputDirectory.Data(),versionString.Data(),(useSingleTopDS)?"_SingleTopDS" : ""));
   
   for (distrUnfold &dist : vecDistr){
      std::vector<TH1F> vec_Unfolded(4);
      std::vector<TH1F> vec_Unfolded_reg(4);
      std::vector<TH1F> vec_Unfolded_bbb(4);
      
      std::vector<TH1F> vec_realDis(4);
      std::vector<TH1F> vec_realDisAlt(4);
      std::vector<TH1F> vec_realDisHerwig(4);
      std::vector<std::map<TString,TH1F>> vec_realDis_systShifts(4);
      std::vector<std::map<TString,TH1F>> vec_realDisAlt_systShifts(4);
      std::vector<std::map<TString,TH1F>> vec_realDisHerwig_systShifts(4);
      
      std::vector<TH2F> vec_response(4);
      std::vector<TH2F> vec_response_fine(4);
      
      std::vector<TH2D> vec_cov_stat(4);
      
      std::vector<std::map<TString,TH1F>> vec_systShifts(4);
      std::vector<std::map<TString,TH1F>> vec_systShifts_reg(4);
      std::vector<std::map<TString,TH1F>> vec_systShifts_bbb(4);
      
      
      TUnfoldBinning generatorBinning;
      
      // Inclusive cross section
      if (dist.varName == "inclusive") {
         for (int i=0; i<4; i++){
            io::RootFileReader histReader(TString::Format("%s/TUnfold/results%.1f.root",cfg_years[i].outputDirectory.Data(),cfg.processFraction*100),resultFolder.Data(),false);
            vec_Unfolded_bbb[i] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded",dist.varName.Data()));
            
            vec_realDis[i] = *histReader.read<TH1F>(TString::Format("%s/Truth",dist.varName.Data()));
            vec_realDisAlt[i] = *histReader.read<TH1F>(TString::Format("%s/TruthAlt",dist.varName.Data()));
            vec_realDisHerwig[i] = *histReader.read<TH1F>(TString::Format("%s/TruthHerwig",dist.varName.Data()));
            
            for (const TString& syst : systVec){
               if (i!=3 && syst == "JESUserDefinedHEM1516_DOWN") continue;    // HEM only used for 2018
               if (Systematic::convertType(syst) == Systematic::CR_envelope || Systematic::convertType(syst) == Systematic::meScale_envelope || Systematic::convertType(syst) == Systematic::mtop) {  // Envelope type systematic are derived within this script and not read from result.root
                  systVec_realDis.push_back(syst);
                  continue;
               }
               vec_systShifts_bbb[i][syst] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded_%s",dist.varName.Data(),syst.Data()));
               vec_systShifts_bbb[i][syst].Multiply(&vec_Unfolded_bbb[i]);   // only relative shifts stored in input
               if(std::find(Systematic::ttbarTypes.begin(), Systematic::ttbarTypes.end(), Systematic::convertType(syst)) != Systematic::ttbarTypes.end()){  //read realDis shifts
                  vec_realDis_systShifts[i][syst] = phys::getSystShift(vec_realDis[i],*histReader.read<TH1F>(TString::Format("%s/Truth_%s",dist.varName.Data(),syst.Data())));
                  vec_realDisAlt_systShifts[i][syst] = phys::getSystShift(vec_realDisAlt[i],*histReader.read<TH1F>(TString::Format("%s/TruthAlt_%s",dist.varName.Data(),syst.Data())));
                  vec_realDisHerwig_systShifts[i][syst] = phys::getSystShift(vec_realDisHerwig[i],*histReader.read<TH1F>(TString::Format("%s/TruthHerwig_%s",dist.varName.Data(),syst.Data())));
                  systVec_realDis.push_back(syst);
               }
            }
         }
         
         TH1F unfolded_bbb = vec_Unfolded_bbb[0];
         TH1F realDis = vec_realDis[0];
         TH1F realDisAlt = vec_realDisAlt[0];
         TH1F realDisHerwig = vec_realDisHerwig[0];
         TH2D cov_syst_bbb = vec_cov_stat[0];
         TH2D cov_syst_realDis = vec_cov_stat[0];
         cov_syst_bbb.Reset();
         cov_syst_realDis.Reset();
         for (int i=1; i<4; i++){
            std::cout<<unfolded_bbb.GetBinContent(1)<<std::endl;
            std::cout<<vec_Unfolded_bbb[i].GetBinContent(1)<<std::endl;
            unfolded_bbb.Add(&vec_Unfolded_bbb[i]);
            realDis.Add(&vec_realDis[i]);
            realDisAlt.Add(&vec_realDisAlt[i]);
            realDisHerwig.Add(&vec_realDisHerwig[i]);
         }
         
         std::map<TString,TH1F> map_combinedShifts_bbb = getCombinedUnc(vec_systShifts_bbb,systVec,unfolded_bbb,vec_Unfolded_bbb,cov_syst_bbb,dist.norm);
         std::map<TString,TH1F> map_combinedShifts_realDis_bbb = getCombinedUnc(vec_realDis_systShifts,systVec_realDis,realDis,vec_realDis,cov_syst_realDis,dist.norm);
         
         for (const auto &[key, value]: map_combinedShifts_bbb){
            map_combinedShifts_bbb[key].Divide(&unfolded_bbb);
         }
         for (const auto &[key, value]: map_combinedShifts_realDis_bbb){
            map_combinedShifts_realDis_bbb[key].Divide(&realDis);
         }
         
         io::RootFileSaver saver(TString::Format("../../../Combined/%s/TUnfold/plots%.1f.root",versionString.Data(),cfg.processFraction*100),TString::Format("%s/%s",dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data(),resultFolder.Data()));
         
         std::cout<<"----------------------------------------------------------"<<std::endl;
         std::cout<<"Measured total cross section:"<<unfolded_bbb.GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section:"<<realDis.GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section (amc@NLO):"<<realDisAlt.GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section (Herwig):"<<realDisHerwig.GetBinContent(1)<<std::endl;
         std::cout<<"-----------------------------Uncertainty on unfolded result------------------------"<<std::endl;
         plot_systBreakdown(map_combinedShifts_bbb,&saver,"SystBreakdown","BBB",dist.varName,useRealData,systCombinations);
         std::cout<<"-----------------------------Uncertainty on powheg prediction------------------------"<<std::endl;
         plot_systBreakdown(map_combinedShifts_realDis_bbb,&saver,"SystBreakdown_realDis","BBB",dist.varName,useRealData,systCombinations);
         std::cout<<"----------------------------------------------------------"<<std::endl;
         
         continue;
      }
      
      systVec_realDis = {};   // reset to work if inclusive is not plotted
            
      
      // Read histograms from input files
      for (int i=0; i<4; i++){
         io::RootFileReader histReader(TString::Format("%s/TUnfold/results%.1f.root",cfg_years[i].outputDirectory.Data(),cfg.processFraction*100),resultFolder.Data(),false);
                  
         vec_Unfolded[i] = *histReader.read<TH1F>(TString::Format("%s/NoReg/unfolded",dist.varName.Data()));
         vec_Unfolded_reg[i] = *histReader.read<TH1F>(TString::Format("%s/Reg/unfolded",dist.varName.Data()));
         vec_Unfolded_bbb[i] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded",dist.varName.Data()));
         
         vec_realDis[i] = *histReader.read<TH1F>(TString::Format("%s/Truth",dist.varName.Data()));
         vec_realDisAlt[i] = *histReader.read<TH1F>(TString::Format("%s/TruthAlt",dist.varName.Data()));
         vec_realDisHerwig[i] = *histReader.read<TH1F>(TString::Format("%s/TruthHerwig",dist.varName.Data()));
         
         vec_response[i] = *histReader.read<TH2F>(TString::Format("%s/ResponseMatrix",dist.varName.Data()));
         vec_response_fine[i] = *histReader.read<TH2F>(TString::Format("%s/ResponseMatrix_fine",dist.varName.Data()));
         
         vec_cov_stat[i] = *histReader.read<TH2D>(TString::Format("%s/NoReg/cov",dist.varName.Data()));
         
         for (const TString& syst : systVec){
            if (i!=3 && syst == "JESUserDefinedHEM1516_DOWN") continue;    // HEM only used for 2018
            // ~if (Systematic::convertType(syst) == Systematic::CR_envelope || Systematic::convertType(syst) == Systematic::meScale_envelope || Systematic::convertType(syst) == Systematic::mtop) {  // Envelope type systematic are derived within this script and not read from result.root
            if (Systematic::convertType(syst) == Systematic::CR_envelope || Systematic::convertType(syst) == Systematic::meScale_envelope || Systematic::convertType(syst) == Systematic::mtop) {  // Envelope type systematic are derived within this script and not read from result.root
               if (Systematic::convertType(syst) == Systematic::meScale_envelope) systVec_realDis.push_back(syst);
               continue;
            }
            vec_systShifts[i][syst] = *histReader.read<TH1F>(TString::Format("%s/NoReg/unfolded_%s",dist.varName.Data(),syst.Data()));
            vec_systShifts_reg[i][syst] = *histReader.read<TH1F>(TString::Format("%s/Reg/unfolded_%s",dist.varName.Data(),syst.Data()));
            vec_systShifts_bbb[i][syst] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded_%s",dist.varName.Data(),syst.Data()));
            
            // Fix CR uncertainty to result from combined 2D binning to avoid small MC stats
            // Since envelope will be derived at the end allocation of UP and Down to the three different sources can be chosen randomly
            if (fixCRunc){
               if (dist.varName == "2D_dPhi_pTnunu_new_30StabPur12Bins_DNN") {
                  if (syst == "CR1") {
                     vec_systShifts[i][syst].SetBinContent(11,CR_up_lastBin);
                     vec_systShifts[i][syst].SetBinContent(12,CR_up_lastBin);
                  }
                  else if (syst == "CR2") {
                     vec_systShifts[i][syst].SetBinContent(11,CR_down_lastBin);
                     vec_systShifts[i][syst].SetBinContent(12,CR_down_lastBin);
                  }
                  else if (syst == "ERDON"){
                     vec_systShifts[i][syst].SetBinContent(11,0.);
                     vec_systShifts[i][syst].SetBinContent(12,0.);
                  }
               }
            }
            
            vec_systShifts[i][syst].Multiply(&vec_Unfolded[i]);   // only relative shifts stored in input
            vec_systShifts_reg[i][syst].Multiply(&vec_Unfolded_reg[i]);   // only relative shifts stored in input
            vec_systShifts_bbb[i][syst].Multiply(&vec_Unfolded_bbb[i]);   // only relative shifts stored in input
            if(std::find(Systematic::ttbarTypes.begin(), Systematic::ttbarTypes.end(), Systematic::convertType(syst)) != Systematic::ttbarTypes.end()){  //read realDis shifts
               vec_realDis_systShifts[i][syst] = phys::getSystShift(vec_realDis[i],*histReader.read<TH1F>(TString::Format("%s/Truth_%s",dist.varName.Data(),syst.Data())));
               vec_realDisAlt_systShifts[i][syst] = phys::getSystShift(vec_realDisAlt[i],*histReader.read<TH1F>(TString::Format("%s/TruthAlt_%s",dist.varName.Data(),syst.Data())));
               vec_realDisHerwig_systShifts[i][syst] = phys::getSystShift(vec_realDisHerwig[i],*histReader.read<TH1F>(TString::Format("%s/TruthHerwig_%s",dist.varName.Data(),syst.Data())));
               systVec_realDis.push_back(syst);
            }
            systVec_realDis.push_back("XSEC_TTSIGNAL");
         }
                  
         if (i==0) generatorBinning = *histReader.read<TUnfoldBinning>(TString::Format("%s/generatorBinning",dist.varName.Data()));
      }
      
      // Combine nominal results, truth level and response and scale to
      TH1F unfolded = vec_Unfolded[0];
      TH1F unfolded_reg = vec_Unfolded_reg[0];
      TH1F unfolded_bbb = vec_Unfolded_bbb[0];
      TH1F realDis = vec_realDis[0];
      TH1F realDisAlt = vec_realDisAlt[0];
      TH1F realDisHerwig = vec_realDisHerwig[0];
      TH2F response = vec_response[0];
      TH2F response_fine = vec_response_fine[0];
      TH2D cov_stat = vec_cov_stat[0];
      TH2D cov_syst = vec_cov_stat[0];
      TH2D cov_syst_reg = vec_cov_stat[0];
      TH2D cov_syst_bbb = vec_cov_stat[0];
      TH2D cov_total = vec_cov_stat[0];
      TH2D cov_syst_realDis = vec_cov_stat[0];
      TH2D cov_syst_realDisAlt = vec_cov_stat[0];
      TH2D cov_syst_realDisHerwig = vec_cov_stat[0];
      cov_syst.Reset();
      cov_syst_reg.Reset();
      cov_syst_bbb.Reset();
      cov_total.Reset();
      cov_syst_realDis.Reset();
      cov_syst_realDisAlt.Reset();
      cov_syst_realDisHerwig.Reset();
      
      for (int i=1; i<4; i++){
         unfolded.Add(&vec_Unfolded[i]);
         unfolded_reg.Add(&vec_Unfolded_reg[i]);
         unfolded_bbb.Add(&vec_Unfolded_bbb[i]);
         realDis.Add(&vec_realDis[i]);
         realDisAlt.Add(&vec_realDisAlt[i]);
         realDisHerwig.Add(&vec_realDisHerwig[i]);
         response.Add(&vec_response[i]);
         response_fine.Add(&vec_response_fine[i]);
         cov_stat.Add(&vec_cov_stat[i]);
      }
      
      // Combined systematic uncertainties and add stat/total uncertainty
      std::map<TString,TH1F> map_combinedShifts = getCombinedUnc(vec_systShifts,systVec,unfolded,vec_Unfolded,cov_syst,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_reg = getCombinedUnc(vec_systShifts_reg,systVec,unfolded_reg,vec_Unfolded_reg,cov_syst_reg,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_bbb = getCombinedUnc(vec_systShifts_bbb,systVec,unfolded_bbb,vec_Unfolded_bbb,cov_syst_bbb,dist.norm);
      
      // Combined systematic uncertainties on MC prediction
      std::map<TString,TH1F> map_combinedShifts_realDis = getCombinedUnc(vec_realDis_systShifts,systVec_realDis,realDis,vec_realDis,cov_syst_realDis,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_realDisAlt = getCombinedUnc(vec_realDisAlt_systShifts,systVec_realDis,realDisAlt,vec_realDisAlt,cov_syst_realDisAlt,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_realDisHerwig = getCombinedUnc(vec_realDisHerwig_systShifts,systVec_realDis,realDisHerwig,vec_realDisHerwig,cov_syst_realDisHerwig,dist.norm);
            
      //Clone hists before plotting and scaling for syst breakdown (to avoid changes)
      TH1F* unfoldedClone = (TH1F*)unfolded.Clone();
      TH1F* unfoldedClone_reg = (TH1F*)unfolded_reg.Clone();
      TH1F* unfoldedClone_bbb = (TH1F*)unfolded_bbb.Clone();
      TH1F* realDisClone = (TH1F*)realDis.Clone();
      TH1F* realDisAltClone = (TH1F*)realDisAlt.Clone();
      TH1F* realDisHerwigClone = (TH1F*)realDisHerwig.Clone();
      
      // Scale to lumi to get cross section
      if(dist.norm){
         // Get cov for normalized results
         cov_stat = getNormCov(cov_stat,unfolded);         
      
         //Get normalized histograms with uncertainties taking bin-to-bin correlations into account
         realDis = hist::getNormalizedHist(realDis);
         realDisAlt = hist::getNormalizedHist(realDisAlt);
         realDisHerwig = hist::getNormalizedHist(realDisHerwig);
         unfolded = hist::getNormalizedHist(unfolded);
         unfolded_reg = hist::getNormalizedHist(unfolded_reg);
         unfolded_bbb = hist::getNormalizedHist(unfolded_bbb);
      }
      else{
         realDis.Scale(1./cfg.lumi);
         realDisAlt.Scale(1./cfg.lumi);
         realDisHerwig.Scale(1./cfg.lumi);
         unfolded.Scale(1./cfg.lumi);
         unfolded_reg.Scale(1./cfg.lumi);
         unfolded_bbb.Scale(1./cfg.lumi);
         cov_stat.Scale(1./(cfg.lumi*cfg.lumi));
         cov_syst.Scale(1./(cfg.lumi*cfg.lumi));
         cov_syst_reg.Scale(1./(cfg.lumi*cfg.lumi));
         cov_syst_bbb.Scale(1./(cfg.lumi*cfg.lumi));
      }
      
      // Get Combined covariance matrix
      cov_total.Add(&cov_stat);
      cov_total.Add(&cov_syst);
      
      // Plot combined result
      int num_bins;
      io::RootFileSaver saver(TString::Format("../../../Combined/%s/TUnfold/plots%.1f.root",versionString.Data(),cfg.processFraction*100),TString::Format("%s/%s",dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data(),resultFolder.Data()));
      std::pair<TH1F*,TH1F*> unfolded_total = getTotalShifts(map_combinedShifts,unfolded,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> unfolded_reg_total = getTotalShifts(map_combinedShifts_reg,unfolded_reg,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> unfolded_bbb_total = getTotalShifts(map_combinedShifts_bbb,unfolded_bbb,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> realDis_syst_total = getTotalShifts(map_combinedShifts_realDis,realDis,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> realDisAlt_syst_total = getTotalShifts(map_combinedShifts_realDisAlt,realDisAlt,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> realDisHerwig_syst_total = getTotalShifts(map_combinedShifts_realDisHerwig,realDisHerwig,dist.norm,cfg.lumi);
      
      for (int i=1; i<=realDis.GetNbinsX(); i++){
         std::cout<<realDis.GetBinContent(i)<<"   "<<realDis_syst_total.first->GetBinContent(i)<<"   "<<realDis_syst_total.second->GetBinContent(i)<<std::endl;
         std::cout<<unfolded.GetBinContent(i)<<"   "<<unfolded_total.first->GetBinContent(i)<<"   "<<unfolded_total.second->GetBinContent(i)<<std::endl;
      }
            
      // ~plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,-1.,&realDis,&realDisAlt,&cov_stat,dist,true,"CombinedResults_Compare",&saver,num_bins,false,false,plotTheo);
      // ~plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,-1.,&realDis,&realDisAlt,&cov_stat,dist,false,"CombinedResults_noreg",&saver,num_bins,false,false,plotTheo,false,false);
      plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,realDis_syst_total,realDisAlt_syst_total,realDisHerwig_syst_total,-1.,&realDis,&realDisAlt,&realDisHerwig,&cov_total,dist,true,(withBSM)? "CombinedResults_Compare_BSM_"+scale_BSM : "CombinedResults_Compare",&saver,num_bins,false,false,plotTheo,chi2_file);
      plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,realDis_syst_total,realDisAlt_syst_total,realDisHerwig_syst_total,-1.,&realDis,&realDisAlt,&realDisHerwig,&cov_total,dist,false,"CombinedResults_noreg",&saver,num_bins,false,false,plotTheo,chi2_file,false,false);

      if(dist.norm){
         unfoldedClone->Scale(1./unfoldedClone->Integral());
         unfoldedClone_reg->Scale(1./unfoldedClone_reg->Integral());
         unfoldedClone_bbb->Scale(1./unfoldedClone_bbb->Integral());
         realDisClone->Scale(1./realDisClone->Integral());
         realDisAltClone->Scale(1./realDisAltClone->Integral());
         realDisHerwigClone->Scale(1./realDisHerwigClone->Integral());
      }
      
      // Divide by nominal for plot of relative syst uncertainties
      for (const auto &[key, value]: map_combinedShifts){
         map_combinedShifts[key].Divide(unfoldedClone);
         map_combinedShifts_reg[key].Divide(unfoldedClone_reg);
         map_combinedShifts_bbb[key].Divide(unfoldedClone_bbb);
      }
      for (const auto &[key, value]: map_combinedShifts_realDis){
         map_combinedShifts_realDis[key].Divide(realDisClone);
         map_combinedShifts_realDisAlt[key].Divide(realDisClone);
         map_combinedShifts_realDisHerwig[key].Divide(realDisClone);
      }
      
      // Plot systematic uncertainty breakdown
      unfolded.SetTitle(dist.title);
      
      plot_systBreakdown(map_combinedShifts,&saver,"SystBreakdown","NoReg",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_reg,&saver,"SystBreakdown","Reg",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_bbb,&saver,"SystBreakdown","BBB",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_realDis,&saver,"SystBreakdown","MC_pred",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_realDisAlt,&saver,"SystBreakdown","MC_predAlt",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_realDisHerwig,&saver,"SystBreakdown","MC_predHerwig",unfolded.GetXaxis()->GetTitle(),useRealData,systCombinations,jesComparison);
      
      // Print CR Uncertainty for merged 2D distribution (used for nominal result since CR unc is limited by mc stats)
      if (dist.varName == "2D_dPhi_pTnunu_new_30StabPur9Bins_sameDet_DNN" && !dist.norm){
         std::cout<<"-------------------CR uncertainty for merged 2D distribution (absolute)--------------------------"<<std::endl;
         std::cout<<"CR_ENVELOPE_UP: "<<map_combinedShifts["CR_ENVELOPE_UP"].GetBinContent(9)<<std::endl;
         std::cout<<"CR_ENVELOPE_DOWN: "<<map_combinedShifts["CR_ENVELOPE_DOWN"].GetBinContent(9)<<std::endl;
         std::cout<<"--------------------------------------------------------------------------------------"<<std::endl;
      }
      
   }
   chi2_file.close();      

}
