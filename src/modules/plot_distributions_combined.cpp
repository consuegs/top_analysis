//Script to combine and plot distributions from distributions.cpp for all periods (taking correct correlations into account) Since reading of inputs histograms is slowing down with progress, one should only plot 1-2 distributions at a time

#include "tools/distributionsPlottingHelper.hpp"

using namespace std::chrono;
using namespace distributionsplotting;

Config &cfg=Config::get();
   
extern "C"
void run()
{
   //Plot BSM
   // ~bool plotBSM = true;
   bool plotBSM = false;
   
   //Select which tW sample is used (DS or DR, DS is default)
   bool useDR = false;
   // ~bool useDR = true;
   
   // Define systematics
   
   //Nominal
   std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_DCTR_UP","MATCH_DCTR_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   // ~//Nominal (without PDF)
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_DCTR_UP","MATCH_DCTR_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
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
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP175p5","MTOP169p5"};
   // ~std::vector<TString> systToPlot = {"Nominal","ERDON"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_DCTR_UP","MATCH_DCTR_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","TOP_PT"};
   // ~std::vector<TString> systToPlot = {"Nominal","TWDS"};
   // ~std::vector<TString> systToPlot = {"TOP_PT"};
   // ~std::vector<TString> systToPlot = {"Nominal"};
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr.push_back({selection+channel,"Lep1_pt",0.,360.,20});
         // ~vecDistr.push_back({selection+channel,"Lep2_pt",0.,306.,17});
         // ~vecDistr.push_back({selection+channel,"mLL",0,204,17});
         // ~vecDistr.push_back({selection+channel,"pTbJet",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"Jet1_pt",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"Jet2_pt",0.,400.,20});
         // ~vecDistr.push_back({selection+channel,"nJets",1.5,8.5,7});
         // ~vecDistr.push_back({selection+channel,"nBjets",0.5,4.5,4});
         // ~vecDistr.push_back({selection+channel,"Lep1_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Lep2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet1_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet2_eta",-2.5,2.5,25});
         
         // ~vecDistr.push_back({selection+channel,"MET_xy",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET_xy",0.,500.,25});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500,50});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0.,3.2,64});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,18,{0,20,40,54,68,84,100,120,140,168,196,228,260,296,332,371,410,455,500}});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0,3.2,12,{0.,0.2,0.4,0.64,0.88,1.12,1.36,1.6,1.84,2.1,2.36,2.74,3.2}});
         
         // DNN inputs
         // ~vecDistr.push_back({selection+channel,"PuppiMET_xy*cos(PuppiMET_xy_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET_xy*sin(PuppiMET_xy_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"MET_xy*cos(MET_xy_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"MET_xy*sin(MET_xy_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_allJet*cos(HT_phi)",-400.,400.,40});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_allJet*sin(HT_phi)",-400.,400.,40});
         // ~vecDistr.push_back({selection+channel,"mass_l1l2_allJet",0.,2000.,100});
         // ~vecDistr.push_back({selection+channel,"Jet1_pt*cos(Jet1_phi)",-400.,400.,40});
         // ~vecDistr.push_back({selection+channel,"Jet1_pt*sin(Jet1_phi)",-400.,400.,40});
         // ~vecDistr.push_back({selection+channel,"Lep1_pt*cos(Lep1_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"Lep1_pt*sin(Lep1_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"MHT",0.,2000.,100});  //M_allJet
         // ~vecDistr.push_back({selection+channel,"CaloMET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"MT2",0.,200.,50});
         // ~vecDistr.push_back({selection+channel,"mjj",0.,2000.,100});
         // ~vecDistr.push_back({selection+channel,"nJets",1.5,8.5,7});
         // ~vecDistr.push_back({selection+channel,"Jet1_E",0,500,50});
         // ~vecDistr.push_back({selection+channel,"HT",0,1000,40});
         // ~vecDistr.push_back({selection+channel,"Jet2_pt*cos(Jet2_phi)",-250.,250.,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_pt*sin(Jet2_phi)",-250.,250.,50});
         
      }
   }
   
   // 2D plots
   std::vector<distr2D> vecDistr2D;
   for(TString selection:{"baseline"}){
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,8,0.,3.2,6,{0,40,70,97.5,125,162.5,200,300,400},{0,0.28,0.56,0.82,1.08,2.14,3.2}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,4,0.,3.2,3,{0,70,125,200,400},{0,0.56,1.08,3.2}});
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
   if (std::find(systToPlot.begin(), systToPlot.end(), "MTOP_UP") != systToPlot.end() || std::find(systToPlot.begin(), systToPlot.end(), "MTOP_DOWN") != systToPlot.end()){
      systToPlot = util::addVectors(systToPlot,{"MTOP175p5","MTOP169p5"});
   }
   
   // loop over periods to get input hists
   for (int i=0; i<4; i++){
      
      std::cout<<"Start importing period "<<cfg_years[i].year<<std::endl;
      
      versionString += cfg_years[i].treeVersion.Data();
      
      // Define Samples per period to use
      std::vector<TString> mcSamples={};
      std::vector<TString> dataSamples={};
      std::vector<TString> ttbarSamples={};
      std::vector<TString> signalSamples={};
      std::vector<TString> stSamples={};
      std::vector<TString> bsmSamples={};
      getSampleVectors(cfg_years[i].year_int,mcSamples,dataSamples,ttbarSamples,signalSamples,stSamples,bsmSamples,useDR);
      std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
      samplesToPlot = util::addVectors(samplesToPlot,bsmSamples);
      
      // Setup systematics
      for (TString syst : systToPlot){
                           
         if (syst == "JESUserDefinedHEM1516_DOWN" && i!=3) continue;   //use HEM1516 only for 2018
         
         if (syst.BeginsWith("CR_ENVELOPE") || syst.BeginsWith("MESCALE_ENVELOPE") || syst.BeginsWith("MTOP_")) continue;   //Envelope and mTop uncertainties handled in getTotalSystCombined (only ingredients added here)
         
         systHists* temp = new systHists(syst,TString::Format("%s/multiHists/%s/histograms_merged_%s.root",cfg_years[i].outputDirectory.Data(),syst.Data(),cfg_years[i].treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal" || syst=="met40Cut")? samplesToPlot : mcSamples, ttbarSamples, signalSamples,stSamples);
         systHists_vec_all[i].push_back(temp);
      }
      
      // Import histograms
      importHists(systHists_vec_all[i],samplesToPlot,mcSamples,vecDistr,vecDistr2D,false);
      
      // Combine lepton channels
      for (auto &current : systHists_vec_all[i]){
         current->combineChannel();
      }
      
      // Define hist collection for nominal
      hs_vec[i] = &(systHists_vec_all[i][0]->hists_);
      
      // combine different samples to improve readability
      combineAllSamples(cfg_years[i].year_int,hs_vec[i],mcSamples_merged,useDR);
      
   }
   
   // combine years (nominals)
   std::vector<TString> samples_merged = util::addVectors(mcSamples_merged,{"data","MC","T2tt_525_350","T2tt_525_438"});
   hist::Histograms<TH1F> hs_combined(samples_merged);
   getCombinedDistributions(hs_combined,hs_vec,samples_merged);
   
   // Define color map and printName map
   std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   std::map<const TString,TString> printNameMap = {{"TTbar_diLepton","t#bar{t} (ll)"},{"tt other","t#bar{t} other"},{"SingleTop","Single top DR"},{"SingleTop_DS","Single top DS"},{"DrellYan_comb","DY+jets"},{"ttW/Z","t#bar{t}W/Z"},{"WJetsToLNu","W+jets"}};
   
   io::RootFileSaver saver(TString::Format("../../../Combined/%s/Distributions/plots%.1f.root",versionString.Data(),cfg.processFraction*100),TString::Format("plot_distributions%s",(useDR)? "_DR" : ""));
   
   // plot 1D distributions
   for (auto const & distr_ : vecDistr){
      
      // get combined syst uncertainty
      // ~auto start = high_resolution_clock::now();      
      plotHistograms(distr_.path,distr_.name,&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,false);
      
      // ~auto stop = high_resolution_clock::now();
      // ~auto duration = duration_cast<microseconds>(stop - start);
      // ~std::cout<<duration.count()<<std::endl;
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // for cutflow plot
            plotHistograms(combinedPath,"all",&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,false);
         }
         else {
            plotHistograms(combinedPath,distr_.name,&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,false);
         }
      }
   }
   
   // plot 2D distributions
   for (auto const & distr_ : vecDistr2D){
      
      plotHistograms(distr_.path,distr_.name,&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,true,distr_.binEdgesY);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // for cutflow plot
            plotHistograms(combinedPath,"all",&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,true,distr_.binEdgesY);
         }
         else {
            plotHistograms(combinedPath,distr_.name,&hs_combined,mcSamples_merged,colormap,printNameMap,systHists_vec_all,saver,plotBSM,true,distr_.binEdgesY);
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
