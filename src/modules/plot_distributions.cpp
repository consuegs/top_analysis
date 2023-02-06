//Script to (re-)plot distributions from distributions.cpp 

#include "tools/distributionsPlottingHelper.hpp"

using namespace std::chrono;
using namespace distributionsplotting;

Config const &cfg=Config::get();
   
extern "C"
void run()
{
   
   std::vector<TString> mcSamples={};
   std::vector<TString> dataSamples={};
   std::vector<TString> ttbarSamples={};
   std::vector<TString> signalSamples={};
   std::vector<TString> stSamples={};
   std::vector<TString> bsmSamples={};
   getSampleVectors(cfg.year_int,mcSamples,dataSamples,ttbarSamples,signalSamples,stSamples,bsmSamples);
   std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
   samplesToPlot = util::addVectors(samplesToPlot,bsmSamples);
   
   //Nominal
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   //Nominal (without PDF)
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JETPILEUPID_UP","JETPILEUPID_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JER_UP","JER_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","ELECTRON_ID_UP","ELECTRON_ID_DOWN","MUON_ID_UP","MUON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ISO_UP","MUON_ISO_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ID_UP","MUON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ISO_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGBC_UP","BTAGBC_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGL_UP","BTAGL_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP","PU_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","UNCLUSTERED_UP","UNCLUSTERED_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","UNCLUSTERED_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MEFACSCALE_UP","MEFACSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MERENSCALE_UP","MERENSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MESCALE_UP","MESCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MEFACSCALE_UP","MEFACSCALE_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MESCALE_UP","MESCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BFRAG_UP","BFRAG_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BFRAG_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR1","CR2","ERDON"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_UP","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","UETUNE_UP","UETUNE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","UETUNE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PDF_50_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_ST_UP","XSEC_ST_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_DY_UP","XSEC_DY_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PDF_ALPHAS_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESRelativeSample_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESBBEC1Year_UP","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESBBEC1Year_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","TWDS"};
   std::vector<TString> systToPlot = {"Nominal"};
   // ~std::vector<TString> systToPlot = {"jetPileupIDapplied"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps_subleading"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps_leading"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps_cleanedJets"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps_loose"};
   // ~std::vector<TString> systToPlot = {"applyJetVetoMaps_HEM1516"};
   // ~std::vector<TString> systToPlot = {"useDNNnoMETcut"};
   // ~std::vector<TString> systToPlot = {"useDNNmumu"};
   // ~std::vector<TString> systToPlot = {"useDNNnoMetCutDY"};
   // ~std::vector<TString> systToPlot = {"removeMetCut"};
   
   // ~for (int i=1;i<=50;i++){
       // ~systToPlot.push_back(TString::Format("PDF_%i_UP",i));
       // ~systToPlot.push_back(TString::Format("PDF_%i_DOWN",i));
   // ~}
   
   // Remove HEM for all years except 2018
   auto itr =std::find(systToPlot.begin(), systToPlot.end(), "JESUserDefinedHEM1516_DOWN");
   if (itr != systToPlot.end() && cfg.year_int != 3){
      systToPlot.erase(itr);
   }
   
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
   
   // Define saver
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("plot_distributions%s",(systToPlot[0]=="Nominal")? "" : ("_"+systToPlot[0]).Data()));
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr.push_back({selection+channel,"Lep_e_pt",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep_mu_pt",0.,600.,50});
         vecDistr.push_back({selection+channel,"Lep1_pt",0.,360.,30});
         vecDistr.push_back({selection+channel,"Lep2_pt",0.,300.,25});
         // ~vecDistr.push_back({selection+channel,"MET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET",0.,500.,7,{0,40,70,110,170,260,370,500}});
         // ~vecDistr.push_back({selection+channel,"PuppiMET",0.,500.,7,{0,20,40,55,70,90,110,140,170,215,260,315,370,435,500}});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0,3.2,40});
         // ~vecDistr.push_back({selection+channel,"met1000",0.,1000.,50});
         vecDistr.push_back({selection+channel,"mLL",0,200,25});
         // ~vecDistr.push_back({selection+channel,"pTsumlep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"sumpTlep",0.,600.,30});
         vecDistr.push_back({selection+channel,"pTbJet",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"bJet_eta",-2.5,2.5,25});
         vecDistr.push_back({selection+channel,"Jet1_pt",0.,600.,60});
         vecDistr.push_back({selection+channel,"Jet2_pt",0.,400.,40});
         // ~vecDistr.push_back({selection+channel,"dPhiMETnearJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETleadJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETlead2Jet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi",0.,3.2,8});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi",0.,3.2,16});
         // ~vecDistr.push_back({selection+channel,"COSdphi_metNearLep",-1.,1,50});
         // ~vecDistr.push_back({selection+channel,"SINdphi_metNearLep",0.,1,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETbJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1bJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_bJetLep2",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_bJetnearLep",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_b1b2",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dR_b1b2",0.,5.,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLep1",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLep2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLepsum",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1Lep2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dR_Lep1Lep2",0.,5.,100});
         vecDistr.push_back({selection+channel,"nJets",1.5,9.5,8});
         vecDistr.push_back({selection+channel,"nBjets",-0.5,4.5,5});
         // ~vecDistr.push_back({selection+channel,"MT2",0.,200.,50});
         // ~vecDistr.push_back({selection+channel,"MT",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"mt_MetLep2",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"mt_MetNextLep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"conMt_Lep1Lep2",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"ST",0.,1500.,50});
         // ~vecDistr.push_back({selection+channel,"HT",0.,2500.,50});
         // ~vecDistr.push_back({selection+channel,"sum_STHT",0.,4000.,50});
         // ~vecDistr.push_back({selection+channel,"sum_mlb",0.,3000.,50});
         // ~vecDistr.push_back({selection+channel,"METunc_Puppi",0.,50.,50});
         vecDistr.push_back({selection+channel,"n_Interactions",0.,100.,50});
         // ~vecDistr.push_back({selection+channel,"Lep1_flavor",0.5,2.5,2});
         // ~vecDistr.push_back({selection+channel,"Lep2_flavor",0.5,2.5,2});
         // ~vecDistr.push_back({selection+channel,"Lep1_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"Lep2_phi",-3.2,3.2,50});
         vecDistr.push_back({selection+channel,"Lep1_eta",-2.5,2.5,25});
         vecDistr.push_back({selection+channel,"Lep2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Lep1_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep2_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Jet1_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_phi",-3.2,3.2,50});
         vecDistr.push_back({selection+channel,"Jet1_eta",-2.5,2.5,25});
         vecDistr.push_back({selection+channel,"Jet2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet1_E",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_E",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETfarJet",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiJet1Jet2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"METsig",0.,200.,50});
         // ~vecDistr.push_back({selection+channel,"MHT",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"looseLeptonVeto",-0.5,1.5,2});
         // ~vecDistr.push_back({selection+channel,"dPhiMETnearJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETfarJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETleadJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETlead2Jet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETbJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1Jet1",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"PFMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"CaloMET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"CaloMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_allJet",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_l1l2_allJet",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"mass_l1l2_allJet",0.,3000.,50});
         // ~vecDistr.push_back({selection+channel,"ratio_vecsumpTlep_vecsumpTjet",0.,20.,50});
         // ~vecDistr.push_back({selection+channel,"mjj",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"Lep1_pt*cos(Lep1_phi)",-250,250,25});
         // ~vecDistr.push_back({"baseline_GOF2D"+channel,"PuppiMET_xy*sin(PuppiMET_xy_phi)_VS_MET_xy*sin(MET_xy_phi)",0.5,36.5,36});
         
         vecDistr.push_back({selection+channel,"MET_xy",0.,500.,50});
         vecDistr.push_back({selection+channel,"PuppiMET_xy",0.,500.,50});
         vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,18,{0,20,40,54,68,84,100,120,140,168,196,228,260,296,332,371,410,455,500}});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,9,{0,40,68,100,140,196,260,332,410,500}});
         vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0,3.2,12,{0.,0.2,0.4,0.64,0.88,1.12,1.36,1.6,1.84,2.1,2.36,2.74,3.2}});
         vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi_xy",0,3.2,12,{0.,0.2,0.4,0.64,0.88,1.12,1.36,1.6,1.84,2.1,2.36,2.74,3.2}});
         vecDistr.push_back({selection+channel,"dphi_metNearLep_xy",0,3.2,12,{0.,0.2,0.4,0.64,0.88,1.12,1.36,1.6,1.84,2.1,2.36,2.74,3.2}});
         
         // ~vecDistr.push_back({selection+channel,"mLL",0,200,25});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,18,{0,20,40,54,68,84,100,120,140,168,196,228,260,296,332,371,410,455,500}});
      }
   }
   
   // 2D plots
   std::vector<distr2D> vecDistr2D;
   for(TString selection:{"baseline"}){
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_Puppi",0.,400.,6,0.,3.14,3,{0,40,80,120,160,230,400},{0,0.7,1.4,3.14}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_Puppi",0.,400.,6,0.,3.14,3,{0,20,40,60,80,100,120,140,160,195,230,400},{0,0.35,0.7,1.05,1.4,2.27,3.14}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,14,0.,3.2,6,{0,20,40,52.5,65,80,95,110,125,142.5,160,180,200,300,400},{0,0.32,0.64,0.96,1.28,2.24,3.2}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,7,0.,3.2,3,{0,40,65,95,125,160,200,400},{0,0.64,1.28,3.2}});
         
         
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,7,0.,3.2,2,{0,65,95,135,175,250,295,400},{0,0.64,3.2}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,7,0.,3.2,2,{0,40,65,80,95,115,135,155,175,212.5,250,272.5,295,347.5,400},{0,0.32,0.64,1.92,3.2}});
         
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,4,0.,3.2,3,{0,70,125,200,400},{0,0.56,1.08,3.2}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,8,0.,3.2,6,{0,40,70,97.5,125,162.5,200,300,400},{0,0.28,0.56,0.82,1.08,2.14,3.2}});
      }
   }
   
   // Setup systematics
   std::vector<systHists*> systHists_vec;
   for (TString syst : systToPlot){
      
      bool isNominal = (std::find(Systematic::nominalTypes.begin(), Systematic::nominalTypes.end(), Systematic::convertType(syst)) != Systematic::nominalTypes.end());
      
      if (syst.BeginsWith("CR_ENVELOPE") || syst.BeginsWith("MESCALE_ENVELOPE") || syst.BeginsWith("MTOP_")) continue;   //Envelope and mTop uncertainties handled in getTotalSystCombined (only ingredients added here)
      
      systHists* temp = new systHists(syst,TString::Format("multiHists/%s/histograms_merged_%s.root",syst.Data(),cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(isNominal)? samplesToPlot : mcSamples, ttbarSamples, signalSamples, stSamples);
      systHists_vec.push_back(temp);
   }
   
   std::vector<TString> mcSamples_merged={};
   
   // Import hists 1D
   importHists(systHists_vec,samplesToPlot,mcSamples,vecDistr,vecDistr2D);
   
   // Combine lepton channels
   for (auto &current : systHists_vec){
      current->combineChannel();
   }
   
   // Define color map
   std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   // Define hist collection for nominal
   hist::Histograms<TH1F>* hs;
   hs = &(systHists_vec[0]->hists_);
   
   // Combine different samples to improve readability
   combineAllSamples(cfg.year_int,hs,mcSamples_merged);
   
   // Loop over distributions to plot
   for (auto const & distr_ : vecDistr){
      plotHistograms(distr_.path,distr_.name,hs,mcSamples_merged,colormap,{systHists_vec},saver,false,false);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // fot cutflow plot
            plotHistograms(combinedPath,"all",hs,mcSamples_merged,colormap,{systHists_vec},saver,false,false);
         }
         else {
            plotHistograms(combinedPath,distr_.name,hs,mcSamples_merged,colormap,{systHists_vec},saver,false,false);
         }
      }
   }
   
   // Loop over 2Ddistributions to plot
   for (auto const & distr_ : vecDistr2D){
      plotHistograms(distr_.path,distr_.name,hs,mcSamples_merged,colormap,{systHists_vec},saver,false,true,distr_.binEdgesY);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // fot cutflow plot
            plotHistograms(combinedPath,"all",hs,mcSamples_merged,colormap,{systHists_vec},saver,false,true,distr_.binEdgesY);
         }
         else {
            plotHistograms(combinedPath,distr_.name,hs,mcSamples_merged,colormap,{systHists_vec},saver,false,true,distr_.binEdgesY);
         }
      }
   }
   
   // ~hist::Histograms<TH1F>* hs;
   // ~hs = &(systHists_vec[0]->hists_);
   // Print total yields
   printTotalYields(hs,{systHists_vec},mcSamples_merged);
   
   // Print uncertainties
   printUncBreakDown(hs,systHists_vec,mcSamples);
      
   // Print shift per sample
   printShiftBySample(hs,systHists_vec,mcSamples);
}
