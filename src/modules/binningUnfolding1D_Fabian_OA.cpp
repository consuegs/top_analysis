//~ //Raw script for 1D Binning Studies (Bachelor Thesis Fabian)
//~ //Implementation of the optimization algorithm for the binning (1D)

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TTreeReader.h>
#include <TVectorF.h>

#include <iomanip> 
#include <iostream>
#include <fstream>
using namespace std;

TH1D getRMS(const TH2D* hist2D){
   TProfile::SetDefaultSumw2();
   TProfile TTbar_profile_RMS=*(hist2D->ProfileX("ProfileRMS",1, -1,"s"));
   TH1D RMS = *(hist2D->ProjectionX("RMS"));
   for (int i=1; i<=TTbar_profile_RMS.GetNbinsX(); i++){
      RMS.SetBinContent(i,TTbar_profile_RMS.GetBinError(i));
      RMS.SetBinError(i,0.0001);
   }
   return RMS;
}

//Calculates the next upper edge of a bin
//hist2D 		Histogram of GenMET and PuppiMET (both Selection)
//resolution	Histogram of PuppiMET-GenMET against PuppiMET
//edges_fine	Vector of the binning of PuppiMET-GenMET
//bothSelection	Histogram of pTnunu with selection on generator and reconstruction level
//genSelection	Histogram of pTnunu with selection on generator level
//bin_start		Index of the bin at which to start the routine
//bin_end		Index of the bin at which the previous routine ended. Currently not used
//*_min			Threshold values for the binning
//logger		Logging obeject to write into the log file
//is_logging	True, when logging is enabled
//edge_40GeV	True, when it should be tried to set a bin edge at 40GeV
//merge_OverflowTrue, when overflow bins should be considered
int get_next_end(const TH2D* hist2D, const TH2D resolution, std::vector<float> edges_fine,  const TH1D bothSelection, const TH1D genSelection, int bin_start, int bin_end, float p_min, float s_min, float N_min, float eff_min, ofstream &logger, bool is_logging, bool edge_40GeV, bool merge_Overflow){
	int current_bin_start = bin_start;
	int current_bin_end = bin_end;
	
	int max_bin;
	if (not merge_Overflow) max_bin = hist2D->GetNbinsX();
	else max_bin = -1;
	TH1D reco = *(hist2D->ProjectionX("xProj", 0, max_bin));
	TH1D gen = *(hist2D->ProjectionY("yProj", 0, max_bin));
	
	//Iterates over bins until a requirements are fulfilled or the last bin is reached
	for(int bin_i = current_bin_start; bin_i <= hist2D->GetNbinsX(); bin_i++) {
		//Try to set a bin edge at 40GeV
		if (edge_40GeV) {
			int Nbins = hist2D->GetNbinsX();
			if (bin_i > 20*Nbins/500 and bin_i < 40*Nbins/500) {
				bin_i = int(40*Nbins/500);
				if(is_logging) logger << "Set bin_i to " + std::to_string(bin_i) + " to try enforcing a bin edge at 40GeV. The current value of bin_i corresponds to a bin edge at " + std::to_string(bin_i*500/Nbins) +"GeV.\n";
			}
		}
		//Efficiency
		std::vector<double> region;
		region.push_back(reco.GetBinLowEdge(current_bin_start));
		region.push_back(reco.GetBinLowEdge(bin_i) + reco.GetBinWidth(bin_i));
		std::vector<float> region_f;
		region_f.push_back(region.at(0));
		region_f.push_back(region.at(1));
		bool mergeOverflow = false;
		if (region.at(1) == 500 and merge_Overflow) mergeOverflow = true;
		TH1D bothSelection_eff = hist::rebinned_double(bothSelection, region, mergeOverflow, false);
		TH1D genSelection_eff = hist::rebinned_double(genSelection, region, mergeOverflow, false);
		TH2D res_rebinned = hist::rebinned_double(resolution, region_f, edges_fine, false, false);
		
		//N_rec, N_gen and N_recgen
		float N_rec = reco.Integral(current_bin_start, bin_i);
		float N_gen = gen.Integral(current_bin_start, bin_i);
		float N_recgen = 0;
		float RMS = 0;
		N_recgen = hist2D->Integral(current_bin_start, bin_i, current_bin_start, bin_i);
		
		
		float p = 0;
		float s = 0;
		float eff = 0;
		if (N_rec > 0) p = N_recgen/N_rec;
		if (N_gen > 0) s = N_recgen/N_gen;
		if (genSelection_eff.GetBinContent(1) > 0) {
			bothSelection_eff.Divide(&genSelection_eff);
			eff = bothSelection_eff.GetBinContent(1);
		}
		RMS = getRMS(&res_rebinned).GetBinContent(1);
		
		bool res_leq_binWidth = false;
		if (RMS <= hist::getWidths(region_f).at(0)) res_leq_binWidth = true;
		if(is_logging) logger << "Current values in get_next_end: current_bin_start=" + std::to_string(current_bin_start)   + ", bin_i=" + std::to_string(bin_i) +  ", N_rec=" + std::to_string(N_rec)+ ", N_gen=" + std::to_string(N_gen)+ ", N_recgen=" + std::to_string(N_recgen) + "\n";
		if(is_logging) logger <<"\t purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", RMS=" + std::to_string(RMS) + "\n";
		//Check if requirements are met
		if (p >= p_min and s >= s_min and N_rec >= N_min and eff >= eff_min and res_leq_binWidth) {
			current_bin_end = bin_i;
			if(is_logging) logger <<"Break: Requirements for bin are satisfied.\n";
			break;
		}
		current_bin_end = bin_i;
	}
	return current_bin_end;
}

//calculates the next lower edge of a bin
//hist2D 		Histogram of GenMET and PuppiMET (both Selection)
//resolution	Histogram of PuppiMET-GenMET against PuppiMET
//edges_fine	Vector of the binning of PuppiMET-GenMET
//bothSelection	Histogram of pTnunu with selection on generator and reconstruction level
//genSelection	Histogram of pTnunu with selection on generator level
//bin_start		Index of the bin at which to start the routine
//bin_end		Index of the bin at which the previous routine ended. Currently not used
//*_min			Threshold values for the binning
//logger		Logging obeject to write into the log file
//is_logging	True, when logging is enabled
//edge_40GeV	True, when it should be tried to set a bin edge at 40GeV
//merge_OverflowTrue, when overflow bins should be considered
int get_next_start(const TH2D* hist2D, const TH2D resolution, std::vector<float> edges_fine,  const TH1D bothSelection, const TH1D genSelection, int bin_start, int bin_end, float p_min, float s_min, float N_min, float eff_min, ofstream &logger, bool is_logging, bool edge_40GeV, bool merge_Overflow){
	int current_bin_start = bin_start;
	int current_bin_end = bin_end;
	
	int max_bin;
	if (not merge_Overflow) max_bin = hist2D->GetNbinsX();
	else max_bin = -1;
	TH1D reco = *(hist2D->ProjectionX("xProj", 0, max_bin));
	TH1D gen = *(hist2D->ProjectionY("yProj", 0, max_bin));
	
	//Iterates over bins until a requirements are fulfilled or the first bin is reached
	for(int bin_i = current_bin_end; bin_i >= 1; bin_i--) {
		//Try to set a bin edge at 40GeV
		if (edge_40GeV) {
			int Nbins = hist2D->GetNbinsX();
			if (bin_i < 60*Nbins/500 and bin_i > 40*Nbins/500) {
				bin_i = int(40*Nbins/500) + 1;
				if(is_logging) logger << "Set bin_i to " + std::to_string(bin_i) + " to try enforcing a bin edge at 40GeV. The current value of bin_i corresponds to a bin edge at " + std::to_string((bin_i-1)*500/Nbins) +"GeV.\n";
			}
		}
		//Efficiency
		std::vector<double> region;
		region.push_back(reco.GetBinLowEdge(bin_i));
		region.push_back(reco.GetBinLowEdge(current_bin_end) + reco.GetBinWidth(current_bin_end));
		std::vector<float> region_f;
		region_f.push_back(region.at(0));
		region_f.push_back(region.at(1));
		bool mergeOverflow = false;
		if (region.at(1) == 500 and merge_Overflow) mergeOverflow = true;
		TH1D bothSelection_eff = hist::rebinned_double(bothSelection, region, mergeOverflow, false);
		TH1D genSelection_eff = hist::rebinned_double(genSelection, region, mergeOverflow, false);
		TH2D res_rebinned = hist::rebinned_double(resolution, region_f, edges_fine, false, false);
		
		//N_rec, N_gen and N_recgen
		float N_rec = reco.Integral(bin_i, current_bin_end);
		float N_gen = gen.Integral(bin_i, current_bin_end);
		float N_recgen = 0;
		float RMS = 0;
		N_recgen = hist2D->Integral(bin_i, current_bin_end, bin_i, current_bin_end);
		
		float p = 0;
		float s = 0;
		float eff = 0;
		if (N_rec > 0) p = N_recgen/N_rec;
		if (N_gen > 0) s = N_recgen/N_gen;
		if (genSelection_eff.GetBinContent(1) > 0) {
			bothSelection_eff.Divide(&genSelection_eff);
			eff = bothSelection_eff.GetBinContent(1);
		}
		RMS = getRMS(&res_rebinned).GetBinContent(1);
		
		//Check if requirements are met
		bool res_leq_binWidth = false;
		if (RMS <= hist::getWidths(region_f).at(0)) res_leq_binWidth = true;
		if(is_logging) logger << "Current values in get_next_end: current_bin_end=" + std::to_string(current_bin_end)   + ", bin_i=" + std::to_string(bin_i) +  ", N_rec=" + std::to_string(N_rec)+ ", N_gen=" + std::to_string(N_gen)+ ", N_recgen=" + std::to_string(N_recgen) + "\n";
		if(is_logging) logger <<"\t purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", RMS=" + std::to_string(RMS) + "\n";
		if (p >= p_min and s >= s_min and N_rec >= N_min and eff >= eff_min and res_leq_binWidth) {
			current_bin_start = bin_i;
			if(is_logging) logger <<"Break: Requirements for bin are satisfied.\n";
			break;
		}
		current_bin_start = bin_i;
	}
	return current_bin_start;
}

Config const &cfg=Config::get();

extern "C"



// Run
void run()
{
	bool is_logging = true;
	if(is_logging) std::cout<<"Logging is performed! See Log_OA.txt for the log."<<std::endl;
	else std::cout<<"No Logging is performed."<<std::endl;
	
   TString treeName="TTbar_diLepton";     //input sample name
   TFile file(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/100.0/Nominal/"+treeName+"_merged.root",cfg.year.Data(),cfg.treeVersion.Data()),"read");    //open input file
   TTreeReader reader("ttbar_res100.0/"+treeName, &file);      //open tree in root file
   
   
   //define variables taken from tree
   TTreeReaderValue<float> PuppiMET   (reader, "PuppiMET_xy");
   TTreeReaderValue<float> DNN_MET   (reader, "DNN_MET_pT");
   TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
   TTreeReaderValue<float> MC_weight   (reader, "N");
   TTreeReaderValue<float> SF_weight   (reader, "SF");
   
   //define additional variables
   // ~float step = 5;			//initial bin width in GeV
   float step = 2;			//initial bin width in GeV
   //Threshold values for purity, stability, efficiency and N_rec
   float p_min = 0.5;
   float s_min = 0.5;
   float N_min = 50;
   float eff_min = 0.2;
   //Add. options
   // ~bool merge_Overflow = false;
   bool merge_Overflow = true;
   // ~bool forward = false;
   bool forward = true;
   bool consider_cut = true;
   bool DNN = true;
   // ~bool DNN = false;
      
   //Calculate Edges of the initial binning and the RMS plot
   std::vector<float> edges;
   for (float i = 0; i<=500; i+=step) {
	   edges.push_back(i);
   }
   std::vector<float> edges_fine;
   for (int i = -500; i<=500; i+=5) edges_fine.push_back(i);
      
   //define empty histograms
   
   //2D histogram
   std::cout<<treeName<<std::endl;
   hist::Histograms<TH2D> hs2D({treeName});
   hs2D.addHist("Binning_Optimization/InitialBinning"   ,";%MET;%pTnunu;EventsBIN"           , edges, hist::getWidths(edges), edges, hist::getWidths(edges));
   hs2D.addHist("Binning_Optimization/Resolution"   ,";%pTnunu;%MET - p_{#scale[.8]{T}}(#nu #bar{#nu});EventsBIN"           , edges, hist::getWidths(edges), edges_fine, hist::getWidths(edges_fine));
      
   //Temporary: 1D histogram
   hist::Histograms<TH1D> hs({treeName});
   hs.addHist("Temp/bothSelection_RecMET", ";%MET;;EventsBIN", edges, hist::getWidths(edges));
   hs.addHist("Temp/bothSelection_GenMET", ";%pTnunu;;EventsBIN", edges, hist::getWidths(edges));
      
   //Histogram without weights for efficiency
   hist::Histograms<TH1D> hs_NoWeight({treeName});
   hs_NoWeight.addHist("NoWeight/bothSelection_pTnunu", ";%pTnunu;;EventsBIN", edges, hist::getWidths(edges));
   hs_NoWeight.addHist("NoWeight/genSelection_pTnunu", ";%pTnunu;;EventsBIN", edges, hist::getWidths(edges));
     
   //set current sample for hist selection
   hs2D.setCurrentSample(treeName); 
   hs.setCurrentSample(treeName); 
   hs_NoWeight.setCurrentSample(treeName); 
      
   //define booleans for selection
   bool rec_selection;
   bool gen_selection;
   
   //loop over events in tree
   int processEvents=cfg.processFraction*reader.GetEntries(true);
   int iEv=0;
   
   while (reader.Next()){
      iEv++;
      if (iEv>processEvents) break;
      if (iEv%(std::max(processEvents/10,1))==0){      //logging stuff
         io::log*".";
         io::log.flush();
      }
      
      //set weight for current event
      hs2D.setFillWeight((*SF_weight)*(*MC_weight));
      hs.setFillWeight((*SF_weight)*(*MC_weight));
     
      // set booleans for selection
      gen_selection=false;
      rec_selection=false;
      
      if (DNN) {
		  if(*PtNuNu>0 and (merge_Overflow or *PtNuNu <= 500)) gen_selection=true;
		  if(*DNN_MET>0 and (merge_Overflow or *DNN_MET <= 500)) rec_selection=true;
		  
		  //Fill Histogramms
		  if(gen_selection && rec_selection){
			  hs2D.fill("Binning_Optimization/InitialBinning", *DNN_MET, *PtNuNu);
			  hs2D.fill("Binning_Optimization/Resolution", *DNN_MET, *DNN_MET - *PtNuNu);
			  hs_NoWeight.fill("NoWeight/bothSelection_pTnunu", *PtNuNu); 
		  } 
		  if(gen_selection){
			  hs_NoWeight.fill("NoWeight/genSelection_pTnunu", *PtNuNu);
		  }
		  
		  if(*DNN_MET>0 and *PtNuNu>0){
			  hs.fill("Temp/bothSelection_RecMET", *DNN_MET);
			  hs.fill("Temp/bothSelection_GenMET", *PtNuNu);
		  }
	  } else {
		  if(*PtNuNu>0 and (merge_Overflow or *PtNuNu <= 500)) gen_selection=true;
		  if(*PuppiMET>0 and (merge_Overflow or *PuppiMET <= 500)) rec_selection=true;
		  
		  //Fill Histogramms
		  if(gen_selection && rec_selection){
			  hs2D.fill("Binning_Optimization/InitialBinning", *PuppiMET, *PtNuNu);
			  hs2D.fill("Binning_Optimization/Resolution", *PuppiMET, *PuppiMET - *PtNuNu);
			  hs_NoWeight.fill("NoWeight/bothSelection_pTnunu", *PtNuNu); 
		  } 
		  if(gen_selection){
			  hs_NoWeight.fill("NoWeight/genSelection_pTnunu", *PtNuNu);
		  }
		  
		  if(*PuppiMET>0 and *PtNuNu>0){
			  hs.fill("Temp/bothSelection_RecMET", *PuppiMET);
			  hs.fill("Temp/bothSelection_GenMET", *PtNuNu);
		  }
	  }
   }
   
   io::log<<"";
   
   //Merge overflow bins
   if (merge_Overflow) {
	    hs2D.mergeOverflow();
	    hs_NoWeight.mergeOverflow(); 
   }
   file.Close();
   
   TString dir = cfg.outputDirectory+"/binningOpt/1D/";
   io::ensurePathForFile(dir);
   //Plot histograms
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"binningUnfolding1D_Fabian_OA");
   TCanvas can;
   //~ can.SetLogz();
   for(TString var : hs2D.getVariableNames()){
      TH2D* tempHist=hs2D.getHistogram(var,{treeName});
      tempHist->SetStats(false);
      tempHist->GetZaxis()->SetTitle("Events/BIN");
      gPad->SetRightMargin(0.2);
      tempHist->Draw("colz");
      saver.save(can,var);
      if (var == "Binning_Optimization/InitialBinning" and not DNN) can.SaveAs(dir + std::to_string(int(step)) + "GeV_" + "Initial2D.pdf");
      if (var == "Binning_Optimization/InitialBinning" and DNN) can.SaveAs(dir + std::to_string(int(step)) + "GeV_DNN_" + "Initial2D.pdf");
   }
   
   
   //Build file appendix
   TString fileAppendix= "varBins_";
   if (forward) fileAppendix += "fw_";
   else fileAppendix += "bw_";
   fileAppendix += std::to_string(int(step)) + "GeV_";
   if (merge_Overflow) fileAppendix += "mo_";
   else fileAppendix += "nmo_";
   if (consider_cut) fileAppendix += "cc_";
   else fileAppendix += "ncc_";
   if (DNN) fileAppendix +=  "DNN_";
   if (not (N_min == 200)) fileAppendix += "N" + std::to_string(int(N_min)) + "_";
   
   
   //Algorithm to optimize the binning
   ofstream log_file;
   if(is_logging) {
	   log_file.open(dir+fileAppendix + "Log.txt", ios::trunc);
	   log_file << "START\n";
	   if (DNN) log_file << "DNN is used\n";
	   else log_file << "PuppiMET ist used\n";
	   log_file << "---------------- Parameters ----------------\n";
	   if (forward) log_file << "Operating Mode: Forward\n" ;
	   else log_file << "Operating Mode: Backward\n" ;
	   log_file << "Initial Binning: " + std::to_string(step) + "GeV steps \n";
	   log_file << "Threshold Values: \n";
	   log_file << "\tPurity >= " + std::to_string(p_min) + " \n";
	   log_file << "\tStability >= " + std::to_string(s_min) + " \n";
	   log_file << "\tEfficiency >= " + std::to_string(eff_min) + " \n";
	   log_file << "\tN_rec >= " + std::to_string(N_min) + " \n";
	   if (merge_Overflow) log_file << "The overflow bin is merged.\n" ;
	   else log_file << "The overflow bin is not merged.\n" ;
	   log_file << "--------------------------------------------\n";
	}
   
   int current_bin_start = 1;
   int current_bin_end = 1;
   
   TH2D* hist=hs2D.getHistogram("Binning_Optimization/InitialBinning",{treeName});
   float n_bins = hist->GetNbinsX();
   
   //Container for the optimized bin edges
   std::vector<float> edges_optimized;
   
   //Calculate bin edges
   if (forward) {
	   edges_optimized.push_back(0);
	   while (current_bin_end < n_bins) {
		   if(is_logging) log_file << "Start merging routine at bin " + std::to_string(current_bin_start) + "\n";
		   current_bin_end = get_next_end(hist,*(hs2D.getHistogram("Binning_Optimization/Resolution",{treeName})), edges_fine,  *(hs_NoWeight.getHistogram("NoWeight/bothSelection_pTnunu", {treeName})), *(hs_NoWeight.getHistogram("NoWeight/genSelection_pTnunu", {treeName})), current_bin_start, current_bin_end, p_min, s_min, N_min, eff_min, log_file, is_logging, consider_cut, merge_Overflow);
		   if(is_logging) log_file << "End merging routine at bin " + std::to_string(current_bin_end) + "\n";
		   if(is_logging) log_file << "\n";
		   
		   edges_optimized.push_back(hist->GetXaxis()->GetBinLowEdge(current_bin_end) + hist->GetXaxis()->GetBinWidth(current_bin_end));
		   current_bin_start = current_bin_end + 1;
	   }
    } else {
	   edges_optimized.push_back(500);
	   current_bin_start = n_bins;
	   current_bin_end = n_bins;
	   while (current_bin_start > 1) {
		   if(is_logging) log_file << "Start merging routine at bin " + std::to_string(current_bin_end) + "\n";
		   current_bin_start = get_next_start(hist,*(hs2D.getHistogram("Binning_Optimization/Resolution",{treeName})), edges_fine,  *(hs_NoWeight.getHistogram("NoWeight/bothSelection_pTnunu", {treeName})), *(hs_NoWeight.getHistogram("NoWeight/genSelection_pTnunu", {treeName})), current_bin_start, current_bin_end, p_min, s_min, N_min, eff_min, log_file, is_logging, consider_cut, merge_Overflow);
		   if(is_logging) log_file << "End merging routine at bin " + std::to_string(current_bin_start) + "\n";
		   if(is_logging) log_file << "\n";
		   
		   edges_optimized.insert(edges_optimized.begin(), hist->GetXaxis()->GetBinLowEdge(current_bin_start));
		   current_bin_end = current_bin_start - 1;
	   }
	}

	
	//Check last bin 
	TH2D rebinned = hist::rebinned_double( *(hs2D.getHistogram("Binning_Optimization/InitialBinning", {treeName})), edges_optimized, edges_optimized, merge_Overflow, false);
	TH2D rebinned_res = hist::rebinned_double( *(hs2D.getHistogram("Binning_Optimization/Resolution", {treeName})), edges_optimized, edges_fine, merge_Overflow, false);
	TH1D rebinned_bothSelection_eff = hist::rebinned_double( *(hs_NoWeight.getHistogram("NoWeight/bothSelection_pTnunu", {treeName})), edges_optimized, hist::getWidths(edges_optimized), merge_Overflow, false);
	TH1D rebinned_genSelection_eff = hist::rebinned_double( *(hs_NoWeight.getHistogram("NoWeight/genSelection_pTnunu", {treeName})), edges_optimized, hist::getWidths(edges_optimized), merge_Overflow, false);
	
	
	float p = 0;
	float s = 0;
	float eff = 0;
	float res = 0;
	float N_rec = 0;
	float N_gen = 0;
	float N_recgen = 0;
	bool res_is_ok = false;
	
	int max_bin;
	if (not merge_Overflow) max_bin = rebinned.GetNbinsX();
	else max_bin = -1;
	
	int bin_to_check;
	if (forward) bin_to_check = rebinned.GetNbinsX();
	else bin_to_check = 1;
	
	//Calculate values of the last/first bin
	N_rec = rebinned.ProjectionX("N_rec", 0, max_bin)->GetBinContent(bin_to_check);
	N_gen = rebinned.ProjectionY("N_gen", 0, max_bin)->GetBinContent(bin_to_check);
	N_recgen = rebinned.GetBinContent(bin_to_check, bin_to_check);
	if (N_rec > 0) p = N_recgen/N_rec;
	if (N_gen > 0) s = N_recgen/N_gen;
	rebinned_bothSelection_eff.Divide(&rebinned_genSelection_eff);	
	eff = rebinned_bothSelection_eff.GetBinContent(bin_to_check);	
	res = getRMS(&rebinned_res).GetBinContent(bin_to_check);
	if (forward) {
		if (res <= hist::getWidths(edges_optimized).at(edges_optimized.size()-2)) res_is_ok = true;
	} else {
		if (res <= hist::getWidths(edges_optimized).at(0)) res_is_ok = true;
	}

	//Perform merging steps until requirements in the last/first bin are met
	while(not(p >= p_min and s >= s_min and N_rec >= N_min and eff >= eff_min and res_is_ok)){
		if (forward) {
			if(is_logging) log_file << "Last bin does not fulfill restrictions: purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", N_rec=" + std::to_string(N_rec) + ", RMS=" + std::to_string(res) + "\n";
			if(is_logging) log_file << "Merging it with bin " + std::to_string(int(round(edges_optimized.at(edges_optimized.size()-3)))) + "..." + std::to_string(int(round(edges_optimized.at(edges_optimized.size()-2)))) + "GeV\n";
			edges_optimized.erase(edges_optimized.end() - 2);
			bin_to_check -= 1;
		} else {
			if(is_logging) log_file << "First bin does not fulfill restrictions: purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", N_rec=" + std::to_string(N_rec) + ", RMS=" + std::to_string(res) + "\n";
			if(is_logging) log_file << "Merging it with bin " + std::to_string(int(round(edges_optimized.at(1)))) + "..." + std::to_string(int(round(edges_optimized.at(2)))) + "GeV\n";
			edges_optimized.erase(edges_optimized.begin() + 1);
			for (int i = 0; i < edges_optimized.size(); i++)  log_file << std::to_string(int(round(edges_optimized.at(i)))) + " ";
		}
		rebinned = hist::rebinned_double( *(hs2D.getHistogram("Binning_Optimization/InitialBinning", {treeName})), edges_optimized, edges_optimized, merge_Overflow, false);
		rebinned_res = hist::rebinned_double( *(hs2D.getHistogram("Binning_Optimization/Resolution", {treeName})), edges_optimized, edges_fine, merge_Overflow, false);
		rebinned_bothSelection_eff = hist::rebinned_double( *(hs_NoWeight.getHistogram("NoWeight/bothSelection_pTnunu", {treeName})), edges_optimized, hist::getWidths(edges_optimized), merge_Overflow, false);
		rebinned_genSelection_eff = hist::rebinned_double( *(hs_NoWeight.getHistogram("NoWeight/genSelection_pTnunu", {treeName})), edges_optimized, hist::getWidths(edges_optimized), merge_Overflow, false);
		
		if (not merge_Overflow) max_bin = rebinned.GetNbinsX();
		else max_bin = -1;
		N_rec = rebinned.ProjectionX("N_rec", 0, max_bin)->GetBinContent(bin_to_check);
		N_gen = rebinned.ProjectionY("N_gen", 0, max_bin)->GetBinContent(bin_to_check);
		N_recgen = rebinned.GetBinContent(bin_to_check, bin_to_check);
		
		if (N_rec > 0) p = N_recgen/N_rec;
		if (N_gen > 0) s = N_recgen/N_gen;
		rebinned_bothSelection_eff.Divide(&rebinned_genSelection_eff);	
		eff = rebinned_bothSelection_eff.GetBinContent(bin_to_check);	
		res = getRMS(&rebinned_res).GetBinContent(bin_to_check);
		if (forward) {
			if (res <= hist::getWidths(edges_optimized).at(edges_optimized.size()-2)) res_is_ok = true;
		} else {
			if (res <= hist::getWidths(edges_optimized).at(0)) res_is_ok = true;
		}
	}
	
	
	if(is_logging) {
		if (forward) log_file << "Last bin fulfills restrictions: purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", N_rec=" + std::to_string(N_rec) + ", RMS=" + std::to_string(res) + "\n";
		else log_file << "First bin fulfills restrictions: purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + ", N_rec=" + std::to_string(N_rec) + ", RMS=" + std::to_string(res) + "\n";
		//~ log_file <<"\t N_rec from bothSelecton: " +  std::to_string(N_rec_F.GetBinContent(1)) + "\n";
		//~ log_file <<"\t Sum of Events(float): " +  std::to_string(N_rec_F.Integral(1, N_rec_F.GetNbinsX())) + "\n";
		//~ log_file <<"\t Sum of Events(double): " +  std::to_string(rebinned.ProjectionX("N_rec", 0, max_bin)->Integral(1, rebinned.ProjectionX("N_rec", 0, max_bin)->GetNbinsX())) + "\n";
		log_file << "Final bin edges in GeV: ";
		for (int i = 0; i < edges_optimized.size(); i++)  log_file << std::to_string(int(edges_optimized.at(i))) + " ";
		log_file << "\n";
		log_file << "\n";
		log_file << "Input for the plotting script:\n";
		log_file << "\t TString PDFfileAppendix = " ;
		log_file << '"' + fileAppendix + '"' + ";\n";
		log_file << "\t std::vector<float> edges = {";
		for (int i = 0; i < edges_optimized.size() - 1; i++) log_file << std::to_string(int(edges_optimized.at(i))) + ", ";
		log_file << std::to_string(int(edges_optimized.at(edges_optimized.size() - 1))) + "}; \n";
		log_file << "END";
		log_file.close();
	}
	
	//Draw final 2D histogram of normalized rec/Gen Events
	//Nominalization of the bins with respect to N_rec
	TH1D* N_rec_plot = rebinned.ProjectionX("N_rec_proj", 1, max_bin);
	for (int i=1; i <= N_rec_plot->GetNbinsX(); i++) {
		for (int j=1; j <= N_rec_plot->GetNbinsX(); j++) {
			rebinned.SetBinContent(i, j, rebinned.GetBinContent(i, j)/N_rec_plot->GetBinContent(i));
			rebinned.SetBinError(i, j, rebinned.GetBinError(i, j)/N_rec_plot->GetBinContent(i));
		}
	}
	std::cout<<"Overflow Events (N_rec): " + std::to_string(hs.getHistogram("Temp/bothSelection_RecMET", {treeName})->GetBinContent(hs.getHistogram("Temp/bothSelection_RecMET", {treeName})->GetNbinsX() + 1))<<std::endl;
	std::cout<<"Overflow Events (N_gen): " + std::to_string(hs.getHistogram("Temp/bothSelection_GenMET", {treeName})->GetBinContent(hs.getHistogram("Temp/bothSelection_GenMET", {treeName})->GetNbinsX() + 1))<<std::endl;
		
	rebinned.SetStats(false);
	rebinned.GetZaxis()->SetTitle("Events(norm.)/BIN");
	//~ rebinned.Draw("colz text");
	rebinned.Draw("colz");
	saver.save(can,"Binning_Optimization/Optmized");
	can.SaveAs(dir + fileAppendix + "Normalized2D.pdf");
	
	//Draw resolution plot of final binning
	can.cd();
	gPad->SetRightMargin(0.1);
	TProfile* profile_res = rebinned_res.ProfileX("_pfx",1,-1,"s");
	TH1D* hist_res = rebinned_res.ProjectionX();
	hist_res->Reset();
	TH1D* binWidths = (TH1D*)hist_res->Clone();
	
	for (int i=1; i<=profile_res->GetNbinsX(); i++){
		hist_res->SetBinContent(i,profile_res->GetBinError(i));
		binWidths->SetBinContent(i,edges_optimized[i]-edges_optimized[i-1]);
	}
	
	hist_res->GetYaxis()->SetTitle("RMS or Width (GeV)");
	hist_res->SetMaximum(1.5*binWidths->GetMaximum());
	hist_res->Draw("hist");
	binWidths->SetLineColor(kRed);
	binWidths->Draw("hist same");
	
	gfx::LegendEntries legE;
	legE.append(*hist_res,"RMS","l");
	legE.append(*binWidths,"Bin width","l");
	TLegend leg=legE.buildLegend(.75,.8,0.9,.95,1);
	leg.SetTextSize(0.035);
	leg.Draw();
	
	saver.save(can,"Resolution_final",true,true,true);
}
