//Script to combine and plot distributions from distributions.cpp for all periods (taking correct correlations into account)

#include "tools/distributionsPlottingHelper.hpp"

using namespace std::chrono;
using namespace distributionsplotting;

Config &cfg=Config::get();
   
extern "C"
void run()
{
   // Define systematics
   
   //Nominal
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESUserDefinedHEM1516_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_IND_UP","MTOP_IND_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP","PU_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_ST_UP","XSEC_ST_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_ST_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","ELECTRON_ID_UP","ELECTRON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","ELECTRON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","LUMI_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN"};
   std::vector<TString> systToPlot = {"Nominal","MTOP_IND_UP","MTOP_IND_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","ERDON"};
   // ~std::vector<TString> systToPlot = {"Nominal"};
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr.push_back({selection+channel,"Lep1_pt",0.,360.,30});
         // ~vecDistr.push_back({selection+channel,"Lep2_pt",0.,300.,25});
         vecDistr.push_back({selection+channel,"mLL",0,200,25});
         // ~vecDistr.push_back({selection+channel,"pTbJet",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"Jet1_pt",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"Jet2_pt",0.,400.,20});
         // ~vecDistr.push_back({selection+channel,"nJets",1.5,9.5,8});
         // ~vecDistr.push_back({selection+channel,"nBjets",-0.5,4.5,5});
         // ~vecDistr.push_back({selection+channel,"Lep1_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Lep2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet1_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet2_eta",-2.5,2.5,25});
      }
   }
   
   cfg.setLumi(137650.);   // only used for plotting reasons
   
   std::vector<Config> cfg_years = {Config("2016_preVFP"),Config("2016_postVFP"),Config("2017"),Config("2018")};
   
   TString versionString = "";   // combined string for all versions
   
   std::vector<TString> mcSamples_merged={};    // placeholder for samples to plot
   
   std::vector<std::vector<systHists*>> systHists_vec_all ={{},{},{},{}};
   std::vector<hist::Histograms<TH1F>*> hs_vec(4);
   
   // add ingredients for CR and ME envelope if required
   if (std::find(systToPlot.begin(), systToPlot.end(), "CR_ENVELOPE_UP") != systToPlot.end() || std::find(systToPlot.begin(), systToPlot.end(), "CR_ENVELOPE_DOWN") != systToPlot.end()){
      systToPlot = util::addVectors(systToPlot,{"CR1","CR2","ERDON"});
   }
   if (std::find(systToPlot.begin(), systToPlot.end(), "MESCALE_ENVELOPE_UP") != systToPlot.end() || std::find(systToPlot.begin(), systToPlot.end(), "MESCALE_ENVELOPE_DOWN") != systToPlot.end()){
      systToPlot = util::addVectors(systToPlot,{"MESCALE_UP","MESCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN"});
   }
   
   // loop over periods to get input hists
   for (int i=0; i<4; i++){
      versionString += cfg_years[i].treeVersion.Data();
      
      // Define Samples per period to use
      std::vector<TString> mcSamples={};
      std::vector<TString> dataSamples={};
      std::vector<TString> ttbarSamples={};
      std::vector<TString> signalSamples={};
      std::vector<TString> stSamples={};
      std::vector<TString> bsmSamples={};
      getSampleVectors(cfg_years[i].year_int,mcSamples,dataSamples,ttbarSamples,signalSamples,stSamples,bsmSamples);
      std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
      
      // Setup systematics
      for (TString syst : systToPlot){
                           
         if (syst == "JESUserDefinedHEM1516_DOWN" && i!=3) continue;   //use HEM1516 only for 2018
         
         if (syst.BeginsWith("CR_ENVELOPE") || syst.BeginsWith("MESCALE_ENVELOPE")) continue;   //Envelope uncertainties handled in getTotalSystCombined (only ingredients added here)
         
         systHists* temp = new systHists(syst,TString::Format("%s/multiHists/%s/histograms_merged_%s.root",cfg_years[i].outputDirectory.Data(),syst.Data(),cfg_years[i].treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal" || syst=="met40Cut")? samplesToPlot : mcSamples,(syst=="CR_ENVELOPE_UP" || syst=="CR_ENVELOPE_DOWN" || syst=="MTOP_DOWN" || syst=="MTOP_UP")? signalSamples : ttbarSamples, signalSamples,stSamples);
         systHists_vec_all[i].push_back(temp);
      }
      
      for (auto const & distr_ : vecDistr){
         std::cout<<"Read input for "<<distr_.path+distr_.name<<" in period "<<i<<std::endl;
         importHists(systHists_vec_all[i],samplesToPlot,mcSamples,{distr_},{},false);
      }
      
      // Combine lepton channels
      for (auto &current : systHists_vec_all[i]){
         current->combineChannel();
      }
      
      // Define hist collection for nominal
      hs_vec[i] = &(systHists_vec_all[i][0]->hists_);
      
      // combine different samples to improve readability
      combineAllSamples(cfg_years[i].year_int,hs_vec[i],mcSamples_merged);
      
   }
   
   
   // combine years (nominals)
   std::vector<TString> samples_merged = util::addVectors(mcSamples_merged,{"data","MC"});
   hist::Histograms<TH1F> hs_combined(samples_merged);
   getCombinedDistributions(hs_combined,hs_vec,samples_merged);
   
   // Define color map
   std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   io::RootFileSaver saver(TString::Format("../../../Combined/%s/Distributions/plots%.1f.root",versionString.Data(),cfg.processFraction*100),"plot_distributions");
   
   for (auto const & distr_ : vecDistr){
      
      // get combined syst uncertainty
      getTotalSystCombined(systHists_vec_all,distr_.path+distr_.name);
      
      plotHistograms(distr_.path,distr_.name,&hs_combined,mcSamples_merged,colormap,systHists_vec_all,saver,false,false);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // for cutflow plot
            plotHistograms(combinedPath,"all",&hs_combined,mcSamples_merged,colormap,systHists_vec_all,saver,false,false);
         }
         else {
            plotHistograms(combinedPath,distr_.name,&hs_combined,mcSamples_merged,colormap,systHists_vec_all,saver,false,false);
         }
      }
   }
   
   // Print total yields
   printTotalYields(&hs_combined,systHists_vec_all,mcSamples_merged);
   
   /*
   // Print uncertainties
   printUncBreakDown(hs,systHists_vec,mcSamples);
      
   // Print shift per sample
   printShiftBySample(hs,systHists_vec,mcSamples);
   */
}
