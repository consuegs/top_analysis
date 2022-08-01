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
   
   //Use scale factor
   bool withScaleFactor = cfg.tunfold_withScaleFactor;
   
   //Plot comparison
   bool plotComparison = cfg.tunfold_plotComparison;
   
   //Plot toy studies
   bool plotToyStudies = cfg.tunfold_plotToyStudies;
   
   //Use alternative pseudo data (amcAtNLO)
   // ~bool useAltReco = true;
   bool useAltReco = false;
   
   //////////////////////////
   // Define Distributions //
   //////////////////////////
   
   std::vector<distrUnfold> vecDistr;
   vecDistr.push_back({"2D_dPhi_pTnunu",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new40",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"2D_dPhi_pTnunu_new40_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",true});
   // ~vecDistr.push_back({"pTnunu",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"dPhi",0,3.2,";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);d#sigma/d|#Delta#phi|(p_{T}^{#nu#nu},nearest l) (pb)","%.1f",false});
   // ~vecDistr.push_back({"pTnunu_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.0f",false});
   // ~vecDistr.push_back({"dPhi_DNN",0,3.2,";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);d#sigma/d|#Delta#phi|(p_{T}^{#nu#nu},nearest l) (pb)","%.1f",false});
   // ~vecDistr.push_back({"pTll",0,400,";p_{T}^{ll} (GeV);d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   // ~vecDistr.push_back({"inclusive",0,1,";Signal Bin;d#sigma/dp_{T}^{#nu#nu} (pb GeV^{-1})","%.1f",false});
   
   /*
   vecDistr.push_back({"pTnunu_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",false,true});
   // ~vecDistr.push_back({"pTnunu_new_DNN",0,500.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",false,true});
   vecDistr.push_back({"dPhi_DNN",0,3.2,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/d|#Delta#phi|(p_{T}^{#nu#nu},nearest l)","%.0f",false,true});
   vecDistr.push_back({"2D_dPhi_pTnunu_new_DNN",0,400.,";p_{T}^{#nu#nu} (GeV);1/#sigma d#sigma/dp_{T}^{#nu#nu} (GeV^{-1})","%.0f",true,true});
   */
   
   //////////////////////////////////
   // Set Systematic Uncertainties //
   //////////////////////////////////
   
   // Nominal
   std::vector<TString> systVec = {"BSEMILEP_UP","BSEMILEP_DOWN","BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESBBEC1_UP","JESBBEC1_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESUserDefinedHEM1516_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","LUMI_UP","LUMI_DOWN","MATCH_UP","MATCH_DOWN","MESCALE_ENVELOPE_UP","MESCALE_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN","MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PU_UP","PU_DOWN","TOP_PT","TRIG_UP","TRIG_DOWN","UETUNE_UP","UETUNE_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN"};
   // ~std::vector<TString> systVec = {"JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN",};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP","XSEC_DY_DOWN"};
   // ~std::vector<TString> systVec = {"XSEC_DY_UP"};
   // ~std::vector<TString> systVec = {"JESUserDefinedHEM1516_DOWN"};
   // ~std::vector<TString> systVec = {"BTAGBC_CORR_UP","BTAGBC_CORR_DOWN","BTAGBC_UNCORR_UP","BTAGBC_UNCORR_DOWN","BTAGL_CORR_UP","BTAGL_CORR_DOWN","BTAGL_UNCORR_UP","BTAGL_UNCORR_DOWN"};
   // ~std::vector<TString> systVec = {"BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN"};
   // ~std::vector<TString> systVec = {"JER_UP","JER_DOWN"};
   // ~std::vector<TString> systVec = {"JEREta0_UP","JEREta0_DOWN","JEREta1_UP","JEREta1_DOWN"};
   // ~std::vector<TString> systVec = {"MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN"};
   // ~std::vector<TString> systVec = {"MUON_ISO_STAT_UP","MUON_ISO_STAT_DOWN","MUON_ISO_SYST_UP","MUON_ISO_SYST_DOWN","MUON_ID_STAT_UP","MUON_ID_STAT_DOWN","MUON_ID_SYST_UP","MUON_ID_SYST_DOWN"};
   // ~std::vector<TString> systVec = {};
   
   //Remove HEM unc. for all year except 2018
   auto itr =std::find(systVec.begin(), systVec.end(), "JESUserDefinedHEM1516_DOWN");
   if (itr != systVec.end() && cfg.year_int != 3){
      systVec.erase(itr);
   }
   
   ////////////////////////////
   // Set Merged Systematics //
   ////////////////////////////
   std::map<TString,std::vector<TString>> systCombinations = {};
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"JES/JER_regr",{"JESRelativeBalreg","JESFlavorRealistic","JESRelativeSampleYear","JESAbsoluteYear","JESAbsolute","JESBBEC1Year","JESBBEC1","JEREta0","JEREta1"}},
      // ~{"BTAG",{"BTAGBC_CORR","BTAGL_CORR","BTAGBC_UNCORR","BTAGL_UNCORR"}},
      // ~{"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID_STAT","MUON_ID_SYST","MUON_ISO_STAT","MUON_ISO_SYST","MUON_SCALE"}},
      // ~{"PS",{"PSISRSCALE","PSFSRSCALE"}},
      // ~{"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
   // ~};
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"JES/JER",{"JESTotal","JER"}},
      // ~{"BTAG",{"BTAGBC","BTAGL"}},
      // ~{"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID","MUON_ISO","MUON_SCALE"}},
      // ~{"BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
      // ~{"THEORY",{"PSISRSCALE","PSFSRSCALE","MESCALE_ENVELOPE","UETUNE","MATCH","BFRAG","BSEMILEP","PDF_ALPHAS","TOP_PT","CR_ENVELOPE","MTOP","PDF_ENVELOPE",}},
   // ~};
   
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"Absolute",{"JESAbsoluteMPFBias","JESAbsoluteScale","JESFragmentation","JESPileUpDataMC","JESPileUpPtRef","JESRelativeFSR","JESSinglePionECAL","JESSinglePionHCAL"}},
      // ~{"AbsoluteYear",{"JESAbsoluteStat","JESRelativeStatFSR","JESTimePtEta"}},
      // ~{"BBEC1",{"JESPileUpPtBB","JESPileUpPtEC1","JESRelativePtBB"}},
      // ~{"BBEC1Year",{"JESRelativeJEREC1","JESRelativePtEC1","JESRelativeStatEC"}}
   // ~};
   
   for (distrUnfold &dist : vecDistr){
   
      //==============================================
      // step 1 : open output file
      io::RootFileSaver saver(TString::Format("TUnfold/plots%.1f.root",cfg.processFraction*100),TString::Format(!withScaleFactor ? "TUnfold_plotting%.1f/%s" : "TUnfold_plotting_SF_%.1f/%s",cfg.processFraction*100,dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));

      //==============================================
      // step 2 : read binning schemes and input histograms
      TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
      TString input_loc_result=(useAltReco)? "TUnfold_results_TTbar_amcatnlo_"+sample_response : "TUnfold_results_"+sample+"_"+sample_response;
      if (withBSM) {
         input_loc+="_BSM";
         input_loc_result+="_BSM";
      }
      if (withPF) {
         input_loc+="_PF";
         input_loc_result+="_PF";
      }
      if (withPuppi) {
         input_loc+="_Puppi";
         input_loc_result+="_Puppi";
      }
      if (withDNN) {
         input_loc+="_DNN";
         input_loc_result+="_DNN";
      }
      if (withSameBins) {
         input_loc+="_SameBins";
         input_loc_result+="_SameBins";
      }
      if (withPTreweight) {
         input_loc+="_PTreweight"+scale;
         input_loc_result+="_PTreweight"+scale;
      }
      io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold/%s/TUnfold_Nominal_%.1f.root" : "TUnfold/%s/TUnfold_SF91_Nominal_%.1f.root",input_loc.Data(),cfg.processFraction*100));
      
      TString input_loc_old = input_loc;
      input_loc += "/"+dist.varName;
      input_loc_result += "/"+dist.varName;
      
      // Print total cross section
      if (dist.varName == "inclusive") {
         TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
         TH1F *realDis=histReader.read<TH1F>((useAltReco)? input_loc+"/histDataTruthAlt" : input_loc+"/histDataTruth");
         
         unfolded_bbb->Scale(1./cfg.lumi);
         realDis->Scale(1./cfg.lumi);
         
         std::map<TString,TH1F> indShifts_bbb;
         std::map<TString,TH2F> indResponse_bbb;
         // ~std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb,indResponse_bbb);         
         std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,false,withScaleFactor,indShifts_bbb,indResponse_bbb);         
         
         std::cout<<"----------------------------------------------------------"<<std::endl;
         std::cout<<"Measured total cross section:"<<unfolded_bbb->GetBinContent(1)<<std::endl;
         std::cout<<"Expected total cross section:"<<realDis->GetBinContent(1)<<std::endl;
         plot_systBreakdown(indShifts_bbb,&saver,"","BBB",dist.varName,systCombinations);
         std::cout<<"----------------------------------------------------------"<<std::endl;
         
         continue;
      } 

      TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
      TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

      if((!detectorBinning)||(!generatorBinning)) {
         cout<<"problem to read binning schemes\n";
      }

      // read histograms
      TH1F *realDis=histReader.read<TH1F>(input_loc+"/histDataTruth");
      TH1F *realDisAlt=histReader.read<TH1F>(input_loc+"/histDataTruthAlt");
      TH1F *realDis_response=histReader.read<TH1F>(input_loc+"/histMCGenRec_projX");
      TH1F *unfolded=histReader.read<TH1F>(input_loc_result+"/hist_unfoldedResult");
      TH1F *unfolded_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_unfoldedResult");
      TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
      TH2F *corrMatrix=histReader.read<TH2F>(input_loc_result+"/corr_matrix");
      TH2F *corrMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/corr_matrix");
      TH2F *covMatrix=histReader.read<TH2F>(input_loc_result+"/cov_output");
      TH2F *covMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/cov_output");
      TH2F *covMatrix_bbb=histReader.read<TH2F>(input_loc_result+"/BBB/cov_output");
      TH2F *responseMatrix=histReader.read<TH2F>(input_loc+"/histMCGenRec_sameBins");
      TH2F *responseMatrix_fine=histReader.read<TH2F>(input_loc+"/histMCGenRec");
      TParameter<float> *tau_par = histReader.read<TParameter<float>>(input_loc_result+"/reg/tau");
      
      if((!realDis)||(!realDis_response)||(!unfolded)) {
         cout<<"problem to read input histograms\n";
      }
      
      // divide by lumi to get xsec or by integral to get normalized xsec
      if(dist.norm){
         realDis->Scale(1./realDis->Integral());
         realDisAlt->Scale(1./realDisAlt->Integral());
         realDis_response->Scale(1./realDis_response->Integral());
         unfolded->Scale(1./unfolded->Integral());
         unfolded_reg->Scale(1./unfolded_reg->Integral());
         unfolded_bbb->Scale(1./unfolded_bbb->Integral());
      }
      else{
         realDis->Scale(1./cfg.lumi);
         realDisAlt->Scale(1./cfg.lumi);
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
      std::pair<TH1F*,TH1F*> unfolded_total = getSystUnc(unfolded,input_loc_result+"/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts,indResponse,dist.norm);
      std::pair<TH1F*,TH1F*> unfolded_reg_total = getSystUnc(unfolded_reg,input_loc_result+"/reg/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_reg,indResponse_reg,dist.norm);
      std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb,indResponse_bbb,dist.norm);
      
      //========================
      // Step 3: plotting
      
      gfx::SplitCan can;
      can.can_.SetWindowSize(1800,600);
      can.pU_.SetLogy();
      can.pU_.cd();  //Use upper part of canvas
            
      //Initialize proper binning for plotting
      TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
      
      int num_bins = binning_met.GetNoElements();
      int num_bins_met = binning_met.GetNoElements();
      
      if (dist.is2D){
         TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning(1)));
         num_bins *= (binning_phi.GetNoElements()-1);
      }
      
      binning_met.ResizeTo(binning_met.GetNoElements()+1);
      binning_met[binning_met.GetNoElements()-1] = dist.xMax;  //Plotting end for overflow bin
      
      Double_t xbins[num_bins+1];
      Double_t xbins_minus[num_bins+1];
      Double_t xbins_plus[num_bins+1];
      xbins[0] = 0;   //Plotting start
      xbins_minus[0] = -1.*(binning_met[1]/4.);
      xbins_plus[0] = binning_met[1]/4.;
      
      
      int phi_bin = 0;
      for (int i=0; i<(num_bins); i++)   {
         xbins[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax;
         xbins_minus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-binning_met[1]/4.;
         xbins_plus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+binning_met[1]/4.;
         if (i%num_bins_met==num_bins_met-1) phi_bin++;
      }
      
      unfolded->SetBins(num_bins,xbins);
      unfolded_total.first->SetBins(num_bins,xbins);
      unfolded_total.second->SetBins(num_bins,xbins);
      if (plotComparison){
         unfolded_reg->SetBins(num_bins,xbins_minus);
         unfolded_reg_total.first->SetBins(num_bins,xbins_minus);
         unfolded_reg_total.second->SetBins(num_bins,xbins_minus);
         unfolded_bbb->SetBins(num_bins,xbins_plus);
         unfolded_bbb_total.first->SetBins(num_bins,xbins_plus);
         unfolded_bbb_total.second->SetBins(num_bins,xbins_plus);
      }
      else{
         unfolded_reg->SetBins(num_bins,xbins);
         unfolded_reg_total.first->SetBins(num_bins,xbins);
         unfolded_reg_total.second->SetBins(num_bins,xbins);
      }
      
      if (!dist.is2D){ // set correct binning for plotting syst breakdown (currently only 1D)
         for(std::pair<TString,TH1F> const &shift : indShifts){
            indShifts[shift.first].SetBins(num_bins,xbins);
            indShifts_reg[shift.first].SetBins(num_bins,xbins);
            indShifts_bbb[shift.first].SetBins(num_bins,xbins);
         }
      }
      
      realDis->SetBins(num_bins,xbins);
      realDisAlt->SetBins(num_bins,xbins);
      realDis_response->SetBins(num_bins,xbins);
      for (int i=1; i<=num_bins; i++) {  //Set proper label for x axis
         int bin_label_no = (i-1)%num_bins_met+1;
         TString label;
         if (bin_label_no == num_bins_met) label = TString::Format(">"+dist.labelFormat,binning_met[bin_label_no-1]);
         else label = TString::Format(dist.labelFormat+"-"+dist.labelFormat,binning_met[bin_label_no-1],binning_met[bin_label_no]);
         unfolded->GetXaxis()->SetBinLabel(i,label);
         realDis->GetXaxis()->SetBinLabel(i,label);
      }
      
      // Divide by bin width
      hist::divideByBinWidth(*realDis);
      hist::divideByBinWidth(*realDisAlt);
      hist::divideByBinWidth(*realDis_response);
      hist::divideByBinWidth(*unfolded_reg);
      hist::divideByBinWidth(*unfolded_reg_total.first);
      hist::divideByBinWidth(*unfolded_reg_total.second);
      hist::divideByBinWidth(*unfolded_bbb);
      hist::divideByBinWidth(*unfolded_bbb_total.first);
      hist::divideByBinWidth(*unfolded_bbb_total.second);
      hist::divideByBinWidth(*unfolded);
      hist::divideByBinWidth(*unfolded_total.first);
      hist::divideByBinWidth(*unfolded_total.second);
      
      // Plotting options
      unfolded->GetXaxis()->SetTickLength(0.);
      unfolded->GetYaxis()->SetTickLength(0.008);
      unfolded->GetXaxis()->SetTitleOffset(1.5);
      unfolded->GetYaxis()->SetTitleOffset(1.5);
      unfolded->GetYaxis()->SetTitleSize(0.04);
      unfolded->GetXaxis()->CenterLabels(false);
         
      unfolded->LabelsOption("v");
      realDis->LabelsOption("v");
      unfolded->SetMaximum(4*unfolded->GetMaximum());
      unfolded->SetMinimum(0.5*unfolded->GetMinimum());
      unfolded->SetLineColor(kBlack);
      unfolded->SetTitle(dist.title);
      unfolded->SetStats(false);
      realDis->SetTitle(dist.title);
      realDis->SetStats(false);
      
      
      unfolded_reg->SetLineColor(kGreen+2);
      unfolded_bbb->SetLineColor(kViolet);
      unfolded_reg->SetMarkerColor(kGreen+2);
      unfolded_bbb->SetMarkerColor(kViolet);
      
      // Setup unc. plotting
      TGraphAsymmErrors unfolded_totalGraph = hist::getErrorGraph(unfolded_total.first,unfolded_total.second,unfolded,true,true);
      TGraphAsymmErrors unfolded_reg_totalGraph = hist::getErrorGraph(unfolded_reg_total.first,unfolded_reg_total.second,unfolded_reg,true,true);
      TGraphAsymmErrors unfolded_bbb_totalGraph = hist::getErrorGraph(unfolded_bbb_total.first,unfolded_bbb_total.second,unfolded_bbb,true,true);
      unfolded_totalGraph.SetLineColor(kBlack);
      unfolded_reg_totalGraph.SetLineColor(kGreen+2);
      unfolded_bbb_totalGraph.SetLineColor(kViolet);
      TGraphErrors unfolded_graph(unfolded);
      TGraphErrors unfolded_reg_graph(unfolded_reg);
      TGraphErrors unfolded_bbb_graph(unfolded_bbb);
      for (int i=0; i<unfolded->GetNbinsX(); i++){    //change x error for plotting
         unfolded_graph.SetPointError(i,binning_met[1]/10.,unfolded->GetBinError(i+1));
         unfolded_reg_graph.SetPointError(i,binning_met[1]/10.,unfolded_reg->GetBinError(i+1));
         unfolded_bbb_graph.SetPointError(i,binning_met[1]/10.,unfolded_bbb->GetBinError(i+1));
      }
      unfolded_graph.SetFillStyle(1001);
      unfolded_reg_graph.SetFillStyle(1001);
      unfolded_bbb_graph.SetFillStyle(1001);
      unfolded_graph.SetLineWidth(0);
      unfolded_reg_graph.SetLineWidth(0);
      unfolded_bbb_graph.SetLineWidth(0);
      unfolded_graph.SetFillColor(kGray);
      unfolded_reg_graph.SetFillColor(kGreen-9);
      unfolded_bbb_graph.SetFillColor(kMagenta-9);
      unfolded_graph.SetMarkerSize(0.4);
      unfolded_reg_graph.SetMarkerSize(0.4);
      unfolded_bbb_graph.SetMarkerSize(0.4);
      
      unfolded_totalGraph.SetLineWidth(1.);
      unfolded_reg_totalGraph.SetLineWidth(1.);
      unfolded_bbb_totalGraph.SetLineWidth(1.);
      unfolded_totalGraph.SetMarkerSize(0.4);
      unfolded_reg_totalGraph.SetMarkerSize(0.4);
      unfolded_bbb_totalGraph.SetMarkerSize(0.4);
      
      realDis->SetLineColor(kRed-6);
      realDis->SetFillColor(kRed-6);
      realDisAlt->SetLineColor(kBlue-6);
      realDisAlt->SetFillColor(kBlue-6);
      
      if (plotComparison) {
         
         // first draw used for axis and stuff
         // ~unfolded->Draw("px0");
         unfolded->Draw("axis");
         
         realDis->Draw("hist same");   //Draw into current canvas (axis are not drawn again due to "same")
         realDisAlt->Draw("hist same");
         
         //draw stat. unc.
         unfolded_graph.Draw("pe2 same");
         unfolded_reg_graph.Draw("pe2 same");
         unfolded_bbb_graph.Draw("pe2 same");
         
         // draw total unc.
         unfolded_totalGraph.Draw("p same");
         unfolded_reg_totalGraph.Draw("p same");
         unfolded_bbb_totalGraph.Draw("p same");
      }
      else {
         unfolded->Draw("axis");  //Draw into current canvas
         
         realDis->Draw("hist same");
         realDisAlt->Draw("hist same");
         
         unfolded_reg_graph.Draw("pe2 same");
         unfolded_reg_totalGraph.Draw("p same");
      }
      
      realDis_response->SetLineColor(kBlue);
      realDis_response->SetFillColor(kBlue);
      // ~realDis_response->Draw("hist same");
      
      TLatex * atext = new TLatex();
      TLine * aline = new TLine();
      if (dist.is2D){
         //Draw vertical lines and binning ranges for deltaPhi
         atext->SetTextSize(0.03);
         aline->SetLineWidth(2);
         aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
         aline->DrawLine(400,unfolded->GetMinimum(),400,unfolded->GetMaximum());
         aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
         atext->DrawLatex(75,0.5*unfolded->GetMaximum(),"0<|#Delta#phi|<0.64");
         atext->DrawLatex(475,0.5*unfolded->GetMaximum(),"0.64<|#Delta#phi|<1.2");
         atext->DrawLatex(875,0.5*unfolded->GetMaximum(),"1.2<|#Delta#phi|<3.14");
      }
      
      //Get Chi2 and NDF
      auto Chi2Pair = getChi2NDF(unfolded,realDis);
      auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
      
      //Draw legend
      gfx::LegendEntries legE;
      if (plotComparison) {
         auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
         auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
         auto Chi2Pair_bbb = getChi2NDF(unfolded_bbb,realDis);
         auto Chi2Pair_corr_bbb = getChi2NDF_withCorr(unfolded_bbb,realDis,covMatrix_bbb);
         // ~legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair.first,Chi2Pair.second,Chi2Pair_corr.first,Chi2Pair_corr.second),"pe");
         // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
         // ~legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_bbb.first,Chi2Pair_bbb.second,Chi2Pair_corr_bbb.first,Chi2Pair_corr_bbb.second),"pe");
         legE.append(unfolded_totalGraph,"NoReg","pe");
         legE.append(unfolded_reg_totalGraph,"Reg","pe");
         legE.append(unfolded_bbb_totalGraph,"BBB","pe");
      }
      else {
         // ~legE.append(*unfolded,"Unfolded","pe");
         // ~atext->DrawLatex(30,10,TString::Format("#chi^{2}/NDF=%.1f/%i",Chi2Pair.first,Chi2Pair.second));
         // ~atext->DrawLatex(30,3,TString::Format("#chi^{2}/NDF(corr.)=%.1f/%i",Chi2Pair_corr.first,Chi2Pair_corr.second));
         
         auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
         auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
         // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
         legE.append(unfolded_reg_totalGraph,TString::Format("Unfolded (#tau=%.5f)",tau_par->GetVal()),"pe");
      }
      legE.append(*realDis,"Powheg","l");
      legE.append(*realDisAlt,"MadGraph","l");
      // ~legE.append(*realDis_response,"MC true ttbar (response)","l");
      TLegend leg=legE.buildLegend(.16,.4,0.28,.57,1);
      leg.SetTextSize(0.03);
      leg.Draw();
      
      unfolded->Draw("axis same");
      
      if (withScaleFactor) {
         atext->DrawLatex(75,10,"With ScaleFactor");
      }
      
      //Change to lower part of canvas
      can.pL_.cd();
      can.pL_.SetBottomMargin(0.57);
      can.pL_.SetTickx(0);
      TH1F ratio;
      TH1F ratio_alt;
      TH1F ratio_response;
      TH1F ratio_unfolded;
      TH1F ratio_unfolded_reg;
      TH1F ratio_unfolded_bbb;
      
      // derive ratios
      ratio=hist::getRatio(*realDis,*realDis,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
      ratio_alt=hist::getRatio(*realDisAlt,*realDis,"ratio",hist::NOERR);
      ratio_response=hist::getRatio(*realDis_response,*realDis,"ratio",hist::NOERR);
      ratio_unfolded=hist::getRatio(*unfolded,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_reg=hist::getRatio(*unfolded_reg,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_bbb=hist::getRatio(*unfolded_bbb,*realDis,"ratio",hist::ONLY1);
      
      // derive syst. ratios
      TGraphAsymmErrors ratio_totalGraph = hist::getRatioAsymmGraph(*unfolded_total.first,*unfolded_total.second,*unfolded,*realDis);
      TGraphAsymmErrors ratio_reg_totalGraph = hist::getRatioAsymmGraph(*unfolded_reg_total.first,*unfolded_reg_total.second,*unfolded_reg,*realDis);
      TGraphAsymmErrors ratio_bbb_totalGraph = hist::getRatioAsymmGraph(*unfolded_bbb_total.first,*unfolded_bbb_total.second,*unfolded_bbb,*realDis);
      
      // setup stat. unc. plotting
      TGraphErrors ratio_graph(&ratio_unfolded);
      TGraphErrors ratio_reg_graph(&ratio_unfolded_reg);
      TGraphErrors ratio_bbb_graph(&ratio_unfolded_bbb);
      for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){    //change x error for plotting
         ratio_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded.GetBinError(i+1));
         ratio_reg_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_reg.GetBinError(i+1));
         ratio_bbb_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_bbb.GetBinError(i+1));
      }
      
      ratio_graph.SetFillStyle(1001);
      ratio_reg_graph.SetFillStyle(1001);
      ratio_bbb_graph.SetFillStyle(1001);
      ratio_graph.SetLineWidth(0);
      ratio_reg_graph.SetLineWidth(0);
      ratio_bbb_graph.SetLineWidth(0);
      ratio_graph.SetFillColor(kGray);
      ratio_reg_graph.SetFillColor(kGreen-9);
      ratio_bbb_graph.SetFillColor(kMagenta-9);
      ratio_graph.SetMarkerSize(0.4);
      ratio_reg_graph.SetMarkerSize(0.4);
      ratio_bbb_graph.SetMarkerSize(0.4);
      
      // setup axis
      ratio.SetMaximum(1.19);
      ratio.SetMinimum(0.81);
      ratio.SetLineColor(kRed-6);
      ratio.SetMarkerColor(kRed-6);
      ratio.GetYaxis()->SetTitleOffset(0.4);
      ratio.GetXaxis()->SetTitleOffset(1.7);
      ratio.GetXaxis()->SetTickLength(0.);
      // ~if (plotComparison) ratio.Draw("hist");
      // ~else ratio.Draw();
      ratio.Draw("hist");
      
      ratio_response.SetLineColor(kBlue);
      ratio_response.SetMarkerColor(kBlue);
      // ~ratio_response.Draw("same");
      
      ratio_alt.SetLineColor(kBlue-6);
      ratio_alt.SetMarkerColor(kBlue-6);
      ratio_alt.Draw("same");
      
      ratio_totalGraph.SetLineWidth(1.);
      ratio_reg_totalGraph.SetLineWidth(1.);
      ratio_bbb_totalGraph.SetLineWidth(1.);
      ratio_totalGraph.SetMarkerSize(0.4);
      ratio_reg_totalGraph.SetMarkerSize(0.4);
      ratio_bbb_totalGraph.SetMarkerSize(0.4);
      
      if (plotComparison) {
         
         // first draw used for axis and stuff
         // ~ratio_unfolded.Draw("pex0 same");
         
         //draw stat. unc
         ratio_graph.Draw("pe2 same");
         ratio_reg_graph.Draw("pe2 same");
         ratio_bbb_graph.Draw("pe2 same");
         
         //draw tot. unc.
         ratio_totalGraph.Draw("p same");
         ratio_reg_totalGraph.Draw("p same");
         ratio_bbb_totalGraph.Draw("p same");
         
      }
      else {
         ratio_reg_graph.Draw("pe2 same");
         ratio_reg_totalGraph.Draw("p same");
         ratio.Draw("axis same");
      }
      
      if (dist.is2D){
         aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
         aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
         aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
      }
      
      //Print rel. uncertainties:
      std::cout<<"---------------------Uncertainties for "<<dist.varName<<"-------------------------------"<<std::endl;
      float meanRelError = 0;
      for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){
         meanRelError = roundf((ratio_totalGraph.GetErrorYlow(i)+ratio_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
         std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
      }
      std::cout<<"-----------------------"<<std::endl;
      for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
         meanRelError = roundf((ratio_reg_totalGraph.GetErrorYlow(i)+ratio_reg_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
         std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
      }
      std::cout<<"-----------------------"<<std::endl;
      for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
         meanRelError = roundf((ratio_bbb_totalGraph.GetErrorYlow(i)+ratio_bbb_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
         std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
      }
      
      //===========================
      // Step 4: save plot
      TString saveName=(useAltReco)? "TTbar_amcatnlo_"+sample_response : sample+"_"+sample_response;
      TString saveName2D=(useAltReco)? "correlations/TTbar_amcatnlo_"+sample_response : "correlations/"+sample+"_"+sample_response;
      if (withBSM) {
         saveName+="_BSM";
         saveName2D+="_BSM";
      }
      if (withPF)  {
         saveName+="_PF";
         saveName2D+="_PF";
      }
      if (withPuppi)  {
         saveName+="_Puppi";
         saveName2D+="_Puppi";
      }
      if (withDNN)  {
         saveName+="_DNN";
         saveName2D+="_DNN";
      }
      if (withSameBins) {
         saveName+="_SameBins";
         saveName2D+="_SameBins";
      }
      if (withPTreweight) {
         saveName+="_PTreweight"+scale;
         saveName2D+="_PTreweight"+scale;
      }
      if (plotComparison) saveName="compareMethods/"+saveName;
      saver.save(can,saveName,true,true);
      
      //Plot response matrix
      plot_response(responseMatrix,saveName,&saver,dist.is2D);
      plot_correlation(corrMatrix,saveName2D,&saver);
      plot_correlation(corrMatrix_reg,saveName2D+"_reg",&saver);
      
      //Plot syst. unc. breakdown
      plot_systBreakdown(indShifts,&saver,saveName,"Nominal",ratio.GetXaxis()->GetTitle(),systCombinations);
      plot_systBreakdown(indShifts_reg,&saver,saveName,"Reg",ratio.GetXaxis()->GetTitle(),systCombinations);
      plot_systBreakdown(indShifts_bbb,&saver,saveName,"BBB",ratio.GetXaxis()->GetTitle(),systCombinations);
      
      //Plot diff in response for syst.
      // ~plot_response_diff(indResponse,get_response(responseMatrix),&saver,saveName);
      // ~plot_response_diff(indResponse,responseMatrix_fine,&saver,saveName);
      
      //Plot toy studies
      if (plotComparison) saveName=sample+"_"+sample_response;
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
            
            temp->SetMarkerColor(kBlue);
            temp_reg->SetMarkerColor(kRed);
            temp_bbb->SetMarkerColor(kBlack);
            tempRMS.SetMarkerColor(kBlue);
            tempRMS_reg.SetMarkerColor(kRed);
            tempRMS_bbb.SetMarkerColor(kBlack);
            
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
            legE_pull.append(*temp,"NoReg "+type+" Mean","p");
            legE_pull.append(tempRMS,"NoReg "+type+" RMS","p");
            legE_pull.append(*temp_reg,"Reg "+type+" Mean","p");
            legE_pull.append(tempRMS_reg,"Reg "+type+" RMS","p");
            legE_pull.append(*temp_bbb,"BBB "+type+" Mean","p");
            legE_pull.append(tempRMS_bbb,"BBB "+type+" RMS","p");
            TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.8,.65,2);
            leg_pull.SetTextSize(0.04);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/prof"+type);
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
            
            temp->SetLineColor(kBlue);
            temp_reg->SetLineColor(kRed);
            temp_bbb->SetLineColor(kBlack);
            temp->SetLineWidth(2);
            temp_reg->SetLineWidth(2);
            temp_bbb->SetLineWidth(2);
            
            temp_bbb->SetStats(0);
            temp_bbb->Draw("hist");
            temp->Draw("hist same");
            temp_reg->Draw("hist same");
            gfx::LegendEntries legE_pull;
            legE_pull.append(*temp,TString::Format("NoReg [#mu=%.3f #sigma=%.3f]",temp->GetMean(),temp->GetRMS()),"l");
            legE_pull.append(*temp_reg,TString::Format("Reg [#mu=%.3f #sigma=%.3f]",temp_reg->GetMean(),temp_reg->GetRMS()),"l");
            legE_pull.append(*temp_bbb,TString::Format("BBB [#mu=%.3f #sigma=%.3f]",temp_bbb->GetMean(),temp_bbb->GetRMS()),"l");
            TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.4,.65,1);
            leg_pull.SetTextSize(0.04);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/hist"+type);
         }
            
      }
      
      // Store results for combination of years
      io::RootFileSaver resultSaver(TString::Format("TUnfold/results%.1f.root",cfg.processFraction*100),TString::Format("%s",dist.norm? (dist.varName+"_norm").Data() : dist.varName.Data()));
      
      resultSaver.save(*realDis,"Truth");
      resultSaver.save(*realDisAlt,"TruthAlt");
      resultSaver.save(*responseMatrix,"ResponseMatrix");
      resultSaver.save(*responseMatrix_fine,"ResponseMatrix_fine");
      
      resultSaver.save(*unfolded,"NoReg/unfolded");
      resultSaver.save(*unfolded_reg,"Reg/unfolded");
      resultSaver.save(*unfolded_bbb,"BBB/unfolded");
      
      for (const auto &[key, value]: indShifts){
         resultSaver.save(indShifts[key],"NoReg/unfolded_"+key);
         resultSaver.save(indShifts_reg[key],"Reg/unfolded_"+key);
         resultSaver.save(indShifts_bbb[key],"BBB/unfolded_"+key);
      }
   }
}
