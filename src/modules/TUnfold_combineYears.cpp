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
   // ~bool plotTheo = true;
   bool plotTheo = false;
      
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
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.2f",false});
   // ~vecDistr.push_back({"inclusive",0,1,";Signal Bin;d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true,true});
   vecDistr.push_back({"2D_dPhi_pTnunu_new_30StabPur12Bins_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",true,true});
   // ~vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false,true});
   vecDistr.push_back({"dPhi_new_DNN",0,3.2,";min[#Delta#phi(p_{T}^{#nu#nu},l)];1/#sigma d#sigma/dmin[#Delta#phi(p_{T}^{#nu#nu},l)] (pb)","%.2f",false,true});
   
   //////////////////////////////////
   // Set Systematic Uncertainties //
   //////////////////////////////////
   
   // Nominal
   std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   //Experimental unc
   // ~std::vector<TString> systVec = {"BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","TRIG_UP","TRIG_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN"};
   
   //Theory unc
    // ~std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","TOP_PT","UETUNE_UP","UETUNE_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Split JES
   // ~std::vector<TString> systVec = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR1","CR2","ERDON","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","MATCH_UP","MATCH_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN","MTOP169p5","MTOP175p5","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systVec = {"TRIG_UP","TRIG_DOWN"};
   // ~std::vector<TString> systVec = {"TWDS"};
   // ~std::vector<TString> systVec = {"JEREta1_UP","JEREta1_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","CR1","CR2","ERDON"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN","XSEC_ST_UP","XSEC_ST_DOWN"};
   // ~std::vector<TString> systVec = {"JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systVec = {"PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systVec = {"PSFSRSCALE_UP","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systVec = {};
   
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"JES",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JESUserDefinedHEM1516"}},
      // ~{"JER",{"JEREta0","JEREta1"}},
      // ~{"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      // ~{"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      // ~{"PS",{"PSISRSCALE","PSFSRSCALE"}},
      // ~{"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
   // ~};
   std::map<TString,std::vector<TString>> systCombinations = {
      {"JET",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JEREta0","JEREta1","JETPILEUPID"}},
      {"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      {"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      {"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
      {"OTHER EXP",{"L1PREFIRING","LUMI","PU","TRIG","UNCLUSTERED"}},
      {"OTHER THEO",{"BSEMILEP","CR_ENVELOPE","MATCH","MTOP","PDF_ALPHAS","PDF_ENVELOPE","PSFSRSCALE","PSISRSCALE","TOP_PT"}},
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
   
   std::vector<Config> cfg_years = {Config("2016_preVFP"),Config("2016_postVFP"),Config("2017"),Config("2018")};
   
   
   for (distrUnfold &dist : vecDistr){
      std::vector<TH1F> vec_Unfolded(4);
      std::vector<TH1F> vec_Unfolded_reg(4);
      std::vector<TH1F> vec_Unfolded_bbb(4);
      
      std::vector<TH1F> vec_realDis(4);
      std::vector<TH1F> vec_realDisAlt(4);
      
      std::vector<TH2F> vec_response(4);
      std::vector<TH2F> vec_response_fine(4);
      
      std::vector<std::map<TString,TH1F>> vec_systShifts(4);
      std::vector<std::map<TString,TH1F>> vec_systShifts_reg(4);
      std::vector<std::map<TString,TH1F>> vec_systShifts_bbb(4);
      
      TString versionString = "";
      
      TUnfoldBinning generatorBinning;
      
      // Inclusive cross section
      if (dist.varName == "inclusive") {
         for (int i=0; i<4; i++){
            io::RootFileReader histReader(TString::Format("%s/TUnfold/results%.1f.root",cfg_years[i].outputDirectory.Data(),cfg.processFraction*100),"",false);
            vec_Unfolded_bbb[i] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded",dist.varName.Data()));
            
            vec_realDis[i] = *histReader.read<TH1F>(TString::Format("%s/Truth",dist.varName.Data()));
            vec_realDisAlt[i] = *histReader.read<TH1F>(TString::Format("%s/TruthAlt",dist.varName.Data()));
            
            for (const TString& syst : systVec){
               if (i!=3 && syst == "JESUserDefinedHEM1516_DOWN") continue;    // HEM only used for 2018
               vec_systShifts_bbb[i][syst] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded_%s",dist.varName.Data(),syst.Data()));
               vec_systShifts_bbb[i][syst].Multiply(&vec_Unfolded_bbb[i]);   // only relative shifts stored in input
            }
            versionString += cfg_years[i].treeVersion.Data();
         }
         
         TH1F unfolded_bbb = vec_Unfolded_bbb[0];
         TH1F realDis = vec_realDis[0];
         TH1F realDisAlt = vec_realDisAlt[0];
         for (int i=1; i<4; i++){
            unfolded_bbb.Add(&vec_Unfolded_bbb[i]);
            realDis.Add(&vec_realDis[i]);
            realDisAlt.Add(&vec_realDisAlt[i]);
         }
         
         std::map<TString,TH1F> map_combinedShifts_bbb = getCombinedUnc(vec_systShifts_bbb,systVec,unfolded_bbb,vec_Unfolded_bbb,dist.norm);
         
         for (const auto &[key, value]: map_combinedShifts_bbb){
            map_combinedShifts_bbb[key].Divide(&unfolded_bbb);
         }
         
         io::RootFileSaver saver(TString::Format("../../../Combined/%s/TUnfold/plots%.1f.root",versionString.Data(),cfg.processFraction*100),TString::Format("%s",dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));
         
         std::cout<<"----------------------------------------------------------"<<std::endl;
         std::cout<<"Measured total cross section:"<<unfolded_bbb.GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section:"<<realDis.GetBinContent(1)<<std::endl;
         plot_systBreakdown(map_combinedShifts_bbb,&saver,"SystBreakdown","BBB",dist.varName,systCombinations);
         std::cout<<"----------------------------------------------------------"<<std::endl;
         
         continue;
      }
            
            
      
      // Read histograms from input files
      for (int i=0; i<4; i++){
         io::RootFileReader histReader(TString::Format("%s/TUnfold/results%.1f.root",cfg_years[i].outputDirectory.Data(),cfg.processFraction*100),"",false);
         
         vec_Unfolded[i] = *histReader.read<TH1F>(TString::Format("%s/NoReg/unfolded",dist.varName.Data()));
         vec_Unfolded_reg[i] = *histReader.read<TH1F>(TString::Format("%s/Reg/unfolded",dist.varName.Data()));
         vec_Unfolded_bbb[i] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded",dist.varName.Data()));
         
         vec_realDis[i] = *histReader.read<TH1F>(TString::Format("%s/Truth",dist.varName.Data()));
         vec_realDisAlt[i] = *histReader.read<TH1F>(TString::Format("%s/TruthAlt",dist.varName.Data()));
         
         vec_response[i] = *histReader.read<TH2F>(TString::Format("%s/ResponseMatrix",dist.varName.Data()));
         vec_response_fine[i] = *histReader.read<TH2F>(TString::Format("%s/ResponseMatrix_fine",dist.varName.Data()));
         
         for (const TString& syst : systVec){
            if (i!=3 && syst == "JESUserDefinedHEM1516_DOWN") continue;    // HEM only used for 2018
            vec_systShifts[i][syst] = *histReader.read<TH1F>(TString::Format("%s/NoReg/unfolded_%s",dist.varName.Data(),syst.Data()));
            vec_systShifts_reg[i][syst] = *histReader.read<TH1F>(TString::Format("%s/Reg/unfolded_%s",dist.varName.Data(),syst.Data()));
            vec_systShifts_bbb[i][syst] = *histReader.read<TH1F>(TString::Format("%s/BBB/unfolded_%s",dist.varName.Data(),syst.Data()));
            vec_systShifts[i][syst].Multiply(&vec_Unfolded[i]);   // only relative shifts stored in input
            vec_systShifts_reg[i][syst].Multiply(&vec_Unfolded_reg[i]);   // only relative shifts stored in input
            vec_systShifts_bbb[i][syst].Multiply(&vec_Unfolded_bbb[i]);   // only relative shifts stored in input
         }
         
         versionString += cfg_years[i].treeVersion.Data();
         
         if (i==0) generatorBinning = *histReader.read<TUnfoldBinning>(TString::Format("%s/generatorBinning",dist.varName.Data()));
      }
      
      // Combine nominal results, truth level and response and scale to
      TH1F unfolded = vec_Unfolded[0];
      TH1F unfolded_reg = vec_Unfolded_reg[0];
      TH1F unfolded_bbb = vec_Unfolded_bbb[0];
      TH1F realDis = vec_realDis[0];
      TH1F realDisAlt = vec_realDisAlt[0];
      TH2F response = vec_response[0];
      TH2F response_fine = vec_response_fine[0];
      
      for (int i=1; i<4; i++){
         unfolded.Add(&vec_Unfolded[i]);
         unfolded_reg.Add(&vec_Unfolded_reg[i]);
         unfolded_bbb.Add(&vec_Unfolded_bbb[i]);
         realDis.Add(&vec_realDis[i]);
         realDisAlt.Add(&vec_realDisAlt[i]);
         response.Add(&vec_response[i]);
         response_fine.Add(&vec_response_fine[i]);
      }
      
      // Combined systematic uncertainties and add stat/total uncertainty
      std::map<TString,TH1F> map_combinedShifts = getCombinedUnc(vec_systShifts,systVec,unfolded,vec_Unfolded,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_reg = getCombinedUnc(vec_systShifts_reg,systVec,unfolded_reg,vec_Unfolded_reg,dist.norm);
      std::map<TString,TH1F> map_combinedShifts_bbb = getCombinedUnc(vec_systShifts_bbb,systVec,unfolded_bbb,vec_Unfolded_bbb,dist.norm);
      
      //Clone hists before plotting and scaling for syst breakdown (to avoid changes)
      TH1F* unfoldedClone = (TH1F*)unfolded.Clone();
      TH1F* unfoldedClone_reg = (TH1F*)unfolded_reg.Clone();
      TH1F* unfoldedClone_bbb = (TH1F*)unfolded_bbb.Clone();
      
      // Scale to lumi to get cross section
      if(dist.norm){
         realDis.Scale(1./realDis.Integral());
         realDisAlt.Scale(1./realDisAlt.Integral());
         unfolded.Scale(1./unfolded.Integral());
         unfolded_reg.Scale(1./unfolded_reg.Integral());
         unfolded_bbb.Scale(1./unfolded_bbb.Integral());
      }
      else{
         realDis.Scale(1./cfg.lumi);
         realDisAlt.Scale(1./cfg.lumi);
         unfolded.Scale(1./cfg.lumi);
         unfolded_reg.Scale(1./cfg.lumi);
         unfolded_bbb.Scale(1./cfg.lumi);
      }
      
      // Plot combined result
      int num_bins;
      io::RootFileSaver saver(TString::Format("../../../Combined/%s/TUnfold/plots%.1f.root",versionString.Data(),cfg.processFraction*100),TString::Format("%s",dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));
      std::pair<TH1F*,TH1F*> unfolded_total = getTotalShifts(map_combinedShifts,unfolded,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> unfolded_reg_total = getTotalShifts(map_combinedShifts_reg,unfolded_reg,dist.norm,cfg.lumi);
      std::pair<TH1F*,TH1F*> unfolded_bbb_total = getTotalShifts(map_combinedShifts_bbb,unfolded_bbb,dist.norm,cfg.lumi);
      plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,-1.,&realDis,&realDisAlt,dist,true,"CombinedResults_Compare",&saver,num_bins,false,false,plotTheo);
      plot_UnfoldedResult(&generatorBinning,&unfolded,&unfolded_reg,&unfolded_bbb,unfolded_total,unfolded_reg_total,unfolded_bbb_total,-1.,&realDis,&realDisAlt,dist,false,"CombinedResults_noreg",&saver,num_bins,false,false,plotTheo,false,false);
      
      

      if(dist.norm){
         unfoldedClone->Scale(1./unfoldedClone->Integral());
         unfoldedClone_reg->Scale(1./unfoldedClone_reg->Integral());
         unfoldedClone_bbb->Scale(1./unfoldedClone_bbb->Integral());
      }
      
      // Divide by nominal for plot of relative syst uncertainties
      for (const auto &[key, value]: map_combinedShifts){
         map_combinedShifts[key].Divide(unfoldedClone);
         map_combinedShifts_reg[key].Divide(unfoldedClone_reg);
         map_combinedShifts_bbb[key].Divide(unfoldedClone_bbb);
      }
      
      // Plot systematic uncertainty breakdown
      unfolded.SetTitle(dist.title);
      
      plot_systBreakdown(map_combinedShifts,&saver,"SystBreakdown","NoReg",unfolded.GetXaxis()->GetTitle(),systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_reg,&saver,"SystBreakdown","Reg",unfolded.GetXaxis()->GetTitle(),systCombinations,jesComparison);
      plot_systBreakdown(map_combinedShifts_bbb,&saver,"SystBreakdown","BBB",unfolded.GetXaxis()->GetTitle(),systCombinations,jesComparison);
      
   }
   
}
