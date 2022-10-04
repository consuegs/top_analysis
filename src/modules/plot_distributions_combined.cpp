//Script to combine and plot distributions from distributions.cpp for all periods (taking correct correlations into account)

#include "tools/distributionsPlottingHelper.hpp"

using namespace std::chrono;
using namespace distributionsplotting;

Config &cfg=Config::get();
   
extern "C"
void run()
{
   // Define systematics
   std::vector<TString> systToPlot = {"Nominal"};
   
   // Remove HEM for all years except 2018
   auto itr =std::find(systToPlot.begin(), systToPlot.end(), "JESUserDefinedHEM1516_DOWN");
   if (itr != systToPlot.end() && cfg.year_int != 3){
      systToPlot.erase(itr);
   }
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         vecDistr.push_back({selection+channel,"Lep2_pt",0.,300.,25});
      }
   }
   
   cfg.setLumi(137650.);   // only used for plotting reasons
   
   std::vector<Config> cfg_years = {Config("2016_preVFP"),Config("2016_postVFP"),Config("2017"),Config("2018")};
   
   TString versionString = "";   // combined string for all versions
   
   std::vector<TString> mcSamples_merged={};    // placeholder for samples to plot
   
   std::vector<std::vector<systHists*>> systHists_vec_all ={{},{},{},{}};
   std::vector<hist::Histograms<TH1F>*> hs_vec(4);
   
   // loop over periods to get input hists
   for (int i=0; i<4; i++){
      versionString += cfg_years[i].treeVersion.Data();
      
      // Define Samples per period to use
      std::vector<TString> mcSamples={};
      std::vector<TString> dataSamples={};
      std::vector<TString> ttbarSamples={};
      std::vector<TString> signalSamples={};
      getSampleVectors(cfg_years[i].year_int,mcSamples,dataSamples,ttbarSamples,signalSamples);
      std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
      
      // Setup systematics
      for (TString syst : systToPlot){
         systHists* temp = new systHists(syst,TString::Format("%s/multiHists/%s/histograms_merged_%s.root",cfg_years[i].outputDirectory.Data(),syst.Data(),cfg_years[i].treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal" || syst=="met40Cut")? samplesToPlot : mcSamples,(syst=="CR_ENVELOPE_UP" || syst=="CR_ENVELOPE_DOWN" || syst=="MTOP_DOWN" || syst=="MTOP_UP")? signalSamples : ttbarSamples, signalSamples);
         systHists_vec_all[i].push_back(temp);
      }
      
      for (auto const & distr_ : vecDistr){
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
   
   // combine years
   std::vector<TString> samples_merged = util::addVectors(mcSamples_merged,{"data","MC"});
   hist::Histograms<TH1F> hs_combined(samples_merged);
   getCombinedDistributions(hs_combined,hs_vec,samples_merged);
   
   // Define color map
   std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   io::RootFileSaver saver(TString::Format("../../../Combined/%s/Distributions/plots%.1f.root",versionString.Data(),cfg.processFraction*100),"plot_distributions");
   
   for (auto const & distr_ : vecDistr){
      
      plotHistograms(distr_.path,distr_.name,&hs_combined,mcSamples_merged,colormap,systHists_vec_all[3],saver,false,false);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // for cutflow plot
            plotHistograms(combinedPath,"all",&hs_combined,mcSamples_merged,colormap,systHists_vec_all[3],saver,false,false);
         }
         else {
            plotHistograms(combinedPath,distr_.name,&hs_combined,mcSamples_merged,colormap,systHists_vec_all[3],saver,false,false);
         }
      }
   }
   
   // Print total yields
   printTotalYields(&hs_combined,systHists_vec_all[3],mcSamples_merged);
   
   /*
   // Print uncertainties
   printUncBreakDown(hs,systHists_vec,mcSamples);
      
   // Print shift per sample
   printShiftBySample(hs,systHists_vec,mcSamples);
   */
}
