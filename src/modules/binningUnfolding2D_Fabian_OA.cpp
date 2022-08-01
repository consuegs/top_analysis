//~ //Raw script for 2D Binning Studies (Bachelor Thesis Fabian)
//~ //Implementation of the optimization algorithm for the binning
//This scribt optmizes the binnings for 2D and produces also the plots for stability/purity/effciency, RMS, N_rec and 
//stacked SM contributions + BSM events, S/sqrt(B) with BSM signal process and SM background; and proportion of SM events in the total events with the optimized binning.
//

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
#include <TLegendEntry.h>
#include <TFrame.h>

#include <iomanip> 
#include <iostream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <math.h>
using namespace std;

void print_container(const std::vector<float>& c) 
{
    for (auto &i : c) {
        std::cout << i << " ";
    }
    std::cout << '\n';
}

TH1D HistTrafo_2D(TH2D* const &hist2D, std::vector<float> binedges_x, std::vector<float> binedges_y){   
	*hist2D=hist::rebinned_double(*hist2D, binedges_x, binedges_y);

	int numBins_x=binedges_x.size()-1;
	int numBins_y=binedges_y.size()-1;
	int numBins=numBins_x*numBins_y;
	double binedges_1d[numBins+1];
	binedges_1d[0]=0;
	int phi_bin = 0;
	for (int i=0; i<(numBins); i++)   {
		binedges_1d[i+1] = binedges_x[i%numBins_x+1]+phi_bin*binedges_x[numBins_x];
		if (i%numBins_x==numBins_x-1) phi_bin++;
	}

	TH1D* tempHist= new TH1D("", "", numBins_x*numBins_y, binedges_1d);
	int binNew = 1;
	for (int j=1; j<=numBins_y; j++){
		for (int i=1; i<=numBins_x; i++){
			tempHist->SetBinContent(binNew, hist2D->GetBinContent(i,j));
			tempHist->SetBinError(binNew, hist2D->GetBinError(i,j));
			binNew++;
		}
	}

	return *tempHist;
}

TString getPtEdges(std::vector<float> edges_pT) {
	TString s = "pT Edges: ";
	for(int i=0; i<edges_pT.size(); i++) {
			s+= " " + std::to_string(int(edges_pT.at(i))) + " ";
	}
	s += "\n";
	return s;
} 

TString getPhiEdges(std::vector<float> edges_dPhi) {
	TString s = "dPhi Edges: ";
	for(int i=0; i<edges_dPhi.size(); i++) {
			s+= " " + TString::Format("%.2f", edges_dPhi.at(i)) + " ";
	}
	s += "\n";
	return s;
}	
	
TH1D getSqrt(const TH1D hist){
	TH1D::SetDefaultSumw2();
	TH1D Sqrt(hist);
	for (int i=1; i<=hist.GetNbinsX(); i++){
		Sqrt.SetBinContent(i, sqrt(hist.GetBinContent(i)));
		Sqrt.SetBinError(i, 0.5/sqrt(hist.GetBinContent(i)));
	}
	return Sqrt;
}

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

//Calculates Index for Rec/Gen 2D plot
int getIndex(float phi, float pT, std::vector<float> edgesPhi, std::vector<float> edgesPT, bool forward){  
	int index = 1;
	for (int i=0; i<edgesPhi.size(); i++) {
		if ((phi >= edgesPhi.at(i)) and (0 != edgesPhi.at(i))) index+=edgesPT.size()-1;
	} 
	for (int i=0; i<edgesPT.size(); i++) {
		if (pT >= edgesPT.at(i) and (0 != edgesPT.at(i))) index+=1;
	} 
	//If the direction is backwards, the indices have to be counted the other way round
	if (not forward) index = (edgesPT.size()-1)*(edgesPhi.size()-1) + 1 - index;
	return index;
} 

TH1D calculateStability(TH2D* const &hist2D) {
	int Nbins = hist2D->GetNbinsX() * hist2D->GetNbinsY();
	std::string name(hist2D->GetName());
	std::string nameG ="N_gen_proj_stab" + name;
	std::string nameRG ="N_recGen_proj_stab" + name;
	std::string nameStab ="Stab" + name;
	TH1D* N_Gen = hist2D->ProjectionX(nameG.c_str(), 0, Nbins);
	TH1D* N_RecGen = hist2D->ProjectionY(nameRG.c_str(), 0, Nbins);

	for (int i=1; i<= N_RecGen->GetNbinsX(); i++){
	   N_RecGen->SetBinContent(i, hist2D->GetBinContent(i,i));
	   N_RecGen->SetBinError(i, hist2D->GetBinError(i,i));
	}
	
	TH1D* stability = hist2D->ProjectionX(nameStab.c_str(), 0, Nbins);
	stability->Divide(N_RecGen, N_Gen);
	
	return *stability;
}

TH1D calculatePurity(TH2D* const &hist2D) {
	int Nbins = hist2D->GetNbinsX() * hist2D->GetNbinsY();
	std::string name(hist2D->GetName());
	std::string nameR ="N_rec_proj_pur" + name;
	std::string nameRG ="N_recGen_proj_pur" + name;
	std::string namePur ="Pur" + name;
	TH1D* N_Rec = hist2D->ProjectionY(nameR.c_str(), 0, Nbins);
	TH1D* N_RecGen = hist2D->ProjectionY(nameRG.c_str(), 0, Nbins);
	
	for (int i=1; i<= N_RecGen->GetNbinsX(); i++){
	   N_RecGen->SetBinContent(i, hist2D->GetBinContent(i,i));
	   N_RecGen->SetBinError(i, hist2D->GetBinError(i,i));
	}
	TH1D* purity= hist2D->ProjectionY(namePur.c_str(), 0, Nbins);
	purity->Divide(N_RecGen, N_Rec);
	
	return *purity;
}

TH1D calculateEfficiency(TH1D* const &bothSel, TH1D* const &genSel) {
	TH1D eff = *bothSel;
	eff.Divide(genSel);
	std::string name(bothSel->GetName());
	name += "_Eff";
	eff.SetName(name.c_str());
	
	return eff;
}

//Merges a phi bin (2D plot)
//i			Number of pT Bins in RecoGen
//m			Bin index of the dPhi bin which has to be merged (lower one of the two) 
TH2D mergePhiRecoGen(TH2D RecoGen, int i, int m) {
	//Add the Bin contents of the lower dPhi Bin to the higher dPhi Bin
	//Merge Gen Bins
	for (int j = (m-1)*i+1; j <= m*i; j++) {
		for (int k = 1; k <= RecoGen.GetNbinsY(); k++) {
			RecoGen.SetBinContent(j+i, k, RecoGen.GetBinContent(j+i, k) + RecoGen.GetBinContent(j, k));
			RecoGen.SetBinError(j+i, k, sqrt(pow(RecoGen.GetBinError(j+i, k),2) + pow(RecoGen.GetBinError(j, k), 2)));
		}
	}
	//Merge Reco Bins
	for (int j = 1; j <= RecoGen.GetNbinsX(); j++) {
		for (int k = (m-1)*i+1; k <= m*i; k++) {
			RecoGen.SetBinContent(j, k+i, RecoGen.GetBinContent(j, k+i) + RecoGen.GetBinContent(j, k));
			RecoGen.SetBinError(j, k+i, sqrt(pow(RecoGen.GetBinError(j, k+i),2) + pow(RecoGen.GetBinError(j, k), 2)));
		}
	}
	
	//Create a new 2D RecoGen Plot with the merged Phi Bin. This leads to i less bins on each axis
	//For the 2D Plot, each bin content will be copied but the dPhi bin m will be skipped.
	std::string name(RecoGen.GetTitle());
	name += "_RGdPhi";
	TH2D newRecoGen(name.c_str(), "", RecoGen.GetNbinsX()-i, 0, RecoGen.GetNbinsX()-i, RecoGen.GetNbinsY()-i, 0, RecoGen.GetNbinsY()-i);
	int jnew=1;
	int knew=1;
	for (int j = 1; j <= RecoGen.GetNbinsX(); j++) {
		for (int k = 1; k <= RecoGen.GetNbinsY(); k++) {
			if (not ( (j >= (m-1)*i+1 	and		 j <= m*i)	or		(k >= (m-1)*i+1 	and		 k <= m*i)) ) {
				newRecoGen.SetBinContent(jnew, knew, RecoGen.GetBinContent(j, k));
				newRecoGen.SetBinError(jnew, knew, RecoGen.GetBinError(j, k));
				if (knew == newRecoGen.GetNbinsY()) {
					knew = 1;
					jnew += 1;
				} else knew +=1;
			}
		}
	}
	newRecoGen.SetStats(false);
	return newRecoGen;
}

//merges a pT bin (2D plot)
//i			Number of pT Bins in RecoGen
//n			Bin number of the pT bin which has to be merged (lower one of the two) 
TH2D mergePtRecoGen(TH2D RecoGen, int i, int n) {
	//Construct the bin edges for the RecoGen hist. 
	std::vector<float> binEdges;
	for (int j=0; j <= RecoGen.GetNbinsX(); j++) {
		if ( not ((j%i) == n)) binEdges.push_back(j); 
	}
	//Rebinning the histogram with merged bins. 
	TH2D tempHist = hist::rebinned_double(RecoGen, binEdges, binEdges, false, false);
	//New histogram since the bin edges are not equidistant
	std::string name(RecoGen.GetTitle());
	name += "_RGpT";
	TH2D newRecoGen(name.c_str(), "", tempHist.GetNbinsX(), 0, tempHist.GetNbinsX(), tempHist.GetNbinsY(), 0, tempHist.GetNbinsY());
	for (int j = 1; j <= tempHist.GetNbinsX(); j++) {
		for (int k = 1; k <= tempHist.GetNbinsY(); k++) {
			newRecoGen.SetBinContent(j, k, tempHist.GetBinContent(j, k));
			newRecoGen.SetBinError(j, k, tempHist.GetBinError(j, k));
		}
	}
	newRecoGen.SetStats(false);
	return newRecoGen;
}

//merges dPhi bin (1D plot)
//i			Number of pT Bins in RecoGen
//m			Bin number of the dPhi bin which has to be merged (lower one of the two) 
TH1D mergePhi1D(TH1D hist, int i, int m) {
	//Add the Bin contents of the lower dPhi Bin to the higher dPhi Bin
	for (int j = (m-1)*i+1; j <= m*i; j++) {
		hist.SetBinContent(j+i, hist.GetBinContent(j+i) + hist.GetBinContent(j));
		hist.SetBinError(j+i, sqrt(pow(hist.GetBinError(j+i),2) + pow(hist.GetBinError(j), 2)));
	}
	
	//Create new 1D Plot with the merged Phi Bin. This leads to i less bins
	//For the 1D Plot, each bin content will be copied but the dPhi bin m will be skipped.
	std::string name(hist.GetTitle());
	name += "_newPhi";
	TH1D newHist(name.c_str(), "", hist.GetNbinsX()-i, 0, hist.GetNbinsX()-i);
	int jnew=1;
	for (int j = 1; j <= hist.GetNbinsX(); j++) {
			if (not ( j >= (m-1)*i+1 	and		 j <= m*i )) {
				newHist.SetBinContent(jnew, hist.GetBinContent(j));
				newHist.SetBinError(jnew, hist.GetBinError(j));
				jnew+=1;
		}
	}
	newHist.SetStats(false);
	return newHist;
}

//merges pT bin (1D plot)
//i			Number of pT Bins in RecoGen
//n			Bin number of the pT bin which has to be merged (lower one of the two) 
TH1D mergePt1D(TH1D hist, int i, int n) {
	//Construct the bin edges for the RecoGen hist. 
	std::vector<double> binEdges;
	for (int j=0; j <= hist.GetNbinsX(); j++) {
		if ( not ((j%i) == n)) binEdges.push_back(j); 
	}
	//Rebinning the histogram with merged bins. 
	TH1D tempHist = hist::rebinned_double(hist, binEdges, false, false);
	//New histogram since the new bin edges are not equidistant with width 1
	std::string name(hist.GetTitle());
	name += "_newPt";
	TH1D newHist(name.c_str(), "", tempHist.GetNbinsX(), 0, tempHist.GetNbinsX());
	for (int j = 1; j <= tempHist.GetNbinsX(); j++) {
		newHist.SetBinContent(j, tempHist.GetBinContent(j));
		newHist.SetBinError(j, tempHist.GetBinError(j));
	}
	newHist.SetStats(false);
	return newHist;
}

//merges dPhi bin (2D RMS plot)
//i			Number of pT Bins in RecoGen
//m			Bin number of the dPhi bin which has to be merged (lower one of the two) 
TH2D mergePhiRes(TH2D Res, int i, int m) {
	//Add the Bin contents of the lower dPhi Bin to the higher dPhi Bin
	for (int j = (m-1)*i+1; j <= m*i; j++) {
		for (int k = 1; k <= Res.GetNbinsY(); k++){
			Res.SetBinContent(j+i, k, Res.GetBinContent(j+i, k) + Res.GetBinContent(j, k));
			Res.SetBinError(j+i, k, sqrt(pow(Res.GetBinError(j+i, k),2) + pow(Res.GetBinError(j, k), 2)));
		}
	}
	
	//Create new 2D Plot with the merged Phi Bin. This leads to i less bins on the x axis
	//For the 2D Plot, each bin content will be copied but the dPhi bin m will be skipped.
	std::string name(Res.GetTitle());
	name += "_newPhiRes";
	TH2D newRes(name.c_str(), "", Res.GetNbinsX()-i, 0, Res.GetNbinsX()-i, Res.GetNbinsY(), Res.GetYaxis()->GetBinLowEdge(1), Res.GetYaxis()->GetBinLowEdge(Res.GetNbinsY()) + Res.GetYaxis()->GetBinWidth(Res.GetNbinsY()));
	int jnew=1;
	for (int j = 1; j <= Res.GetNbinsX(); j++) {
		if (not ( j >= (m-1)*i+1 	and		 j <= m*i )) {
			for (int k = 1; k <= Res.GetNbinsY(); k++) {
				newRes.SetBinContent(jnew, k, Res.GetBinContent(j,k));
				newRes.SetBinError(jnew, k, Res.GetBinError(j,k));
			}
			jnew+=1;
		}
	}
	newRes.SetStats(false);
	return newRes;
}

//merges pT bin (2D RMS plot)
//i			Number of pT Bins in RecoGen
//n			Bin number of the pT bin which has to be merged (lower one of the two) 
TH2D mergePtRes(TH2D Res, int i, int n) {
	//Construct the bin edges for the RecoGen Res. 
	std::vector<float> binEdges_x;
	std::vector<float> binEdges_y;
	for (int j=0; j <= Res.GetNbinsX(); j++) {
		if ( not ((j%i) == n)) binEdges_x.push_back(j); 
	}
	for (int j=1; j <= Res.GetNbinsY(); j++) {
		binEdges_y.push_back(Res.GetYaxis()->GetBinLowEdge(j)); 
	}
	binEdges_y.push_back(Res.GetYaxis()->GetBinLowEdge(Res.GetNbinsY()) + Res.GetYaxis()->GetBinWidth(Res.GetNbinsY()));
	//Rebinning the histogram with merged bins. 
	TH2D tempRes = hist::rebinned_double(Res, binEdges_x, binEdges_y, false, false);
	//New Histogram since the new bin edges are not equidistant with width 1
	std::string name(Res.GetTitle());
	name += "_newPtRes";
	TH2D newRes(name.c_str(), "", tempRes.GetNbinsX(), 0, tempRes.GetNbinsX(), tempRes.GetNbinsY(), binEdges_y.at(0), binEdges_y.at(binEdges_y.size()-1));
	for (int j = 1; j <= tempRes.GetNbinsX(); j++) {
		for (int k = 1; k <= Res.GetNbinsY(); k++) {
			newRes.SetBinContent(j, k, tempRes.GetBinContent(j, k));
			newRes.SetBinError(j, k, tempRes.GetBinError(j, k));
		}
	}
	newRes.SetStats(false);
	return newRes;
}

void PlotStatistics(TH2D* const &hist2D, std::vector<float> binedges_x, std::vector<float> binedges_y, bool Puppi, bool savePDF, TString dir, TString filename, io::RootFileSaver* saver){
	
	TCanvas canvas;
	canvas.SetLogy();
	gPad->SetLeftMargin(0.17);

	TH1D stat = HistTrafo_2D(hist2D, binedges_x, binedges_y);
	stat.SetMaximum(stat.GetMaximum()*1.5);
	stat.SetStats(false);
	stat.GetXaxis()->SetTitle("%MET");
	stat.GetYaxis()->SetTitle("Events/Bin");
	stat.SetLineColor(kBlack);
	stat.SetMarkerSize(0);
	stat.Draw("hist e");
	canvas.Update();

	int phiBins = binedges_y.size() - 1;
	stat.GetXaxis()->SetNdivisions(8*phiBins, false);
	stat.GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
	for (int i=0; i<=phiBins-1; i++){
		stat.GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
		stat.GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
		stat.GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
		stat.GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
		stat.GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
		stat.GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
		stat.GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
		stat.GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
	}
	TLine * aline = new TLine();
	TLatex * atext = new TLatex();
	atext->SetTextSize(0.03);
	aline->SetLineWidth(2);
	aline->SetLineColor(kBlack);
	
	//GetY returns the values in the form of 10^x where x is the return
	double_t upperEndS=canvas.GetFrame()->GetY2();
	double_t lowerEndS=canvas.GetFrame()->GetY1();
	double labelPos = pow(10, lowerEndS + 0.96*(upperEndS-lowerEndS));
	upperEndS = pow(10, upperEndS);
	lowerEndS = pow(10, lowerEndS);

	for (int i=1; i<=phiBins-1; i++){
		aline->DrawLine(i*400,lowerEndS,i*400,upperEndS);
	}
	for (int i=0;i < phiBins; i++){
		if (i==0) atext->DrawLatex(100,labelPos,TString::Format("|#Delta#phi|<%.1f",binedges_y.at(1))); 
		else {
			if (i==phiBins-1) atext->DrawLatex(400*i+100, labelPos,TString::Format("%.1f<|#Delta#phi|",binedges_y.at(i)));
			else atext->DrawLatex(400*i+100, labelPos,TString::Format("%.1f<|#Delta#phi|<%.1f",binedges_y.at(i), binedges_y.at(i+1)));
		}
	}
	TLatex labelStatPuppi = gfx::cornerLabel("PUPPI", 3);
	TLatex labelStatDNN = gfx::cornerLabel("DNN", 3);

	if (Puppi) labelStatPuppi.Draw();
	else labelStatDNN.Draw();

    saver->save(canvas, filename);
    if (savePDF) canvas.SaveAs(dir + filename.ReplaceAll("/", "_") + ".pdf");
    canvas.Close();
}

void PlotStabPurEff(TH2D* const &pT_dPhi, TH2D* const &RecoGen, TH1D* const &bothSel, TH1D* const &genSel, std::vector<float> binedges_x, std::vector<float> binedges_y, bool Puppi, bool savePDF, TString dir, TString filename, io::RootFileSaver* saver){
	TCanvas canvas;
	int Nbins = (binedges_x.size()-1)*(binedges_y.size()-1);
	
	TH1D stab_final = HistTrafo_2D(pT_dPhi, binedges_x, binedges_y);
	TH1D stab = calculateStability(RecoGen);
	for (int i=1; i <=Nbins; i++) {
		stab_final.SetBinContent(i, stab.GetBinContent(i));
		stab_final.SetBinError(i, stab.GetBinError(i));
	}
	
	gPad->SetLeftMargin(0.17);
	float lowerEnd=0.18;
	float upperEnd=0.78;
	
	stab_final.SetName("stab_final");
	stab_final.SetLineColor(kRed);
	stab_final.SetFillColor(kRed);
	stab_final.SetMarkerStyle(20);
	stab_final.SetMarkerSize(0.8);
	stab_final.SetMarkerColor(kRed);
	stab_final.GetYaxis()->SetRangeUser(lowerEnd, upperEnd);
	stab_final.SetStats(false);
	stab_final.GetXaxis()->SetTitle("%pTnunu");
	stab_final.GetYaxis()->SetTitle("");
	stab_final.Draw("pe1");

	int phiBins = binedges_y.size() - 1;
	stab_final.GetXaxis()->SetNdivisions(8*phiBins, false);
	stab_final.GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
	for (int i=0; i<=phiBins-1; i++){
		stab_final.GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
		stab_final.GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
		stab_final.GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
		stab_final.GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
		stab_final.GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
		stab_final.GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
		stab_final.GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
		stab_final.GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
	}
	
	TH1D pur_final = HistTrafo_2D(pT_dPhi, binedges_x, binedges_y);
	TH1D pur = calculatePurity(RecoGen);
	for (int i=1; i <=Nbins; i++) {
		pur_final.SetBinContent(i, pur.GetBinContent(i));
		pur_final.SetBinError(i, pur.GetBinError(i));
	}
	
	pur_final.SetName("pur_final");
	pur_final.SetLineColor(kBlue);
	pur_final.SetFillColor(kBlue);
	pur_final.SetMarkerStyle(21);
	pur_final.SetMarkerSize(0.8);
	pur_final.SetMarkerColor(kBlue);
	pur_final.Draw("same pe1");
	
	TH1D eff_final = HistTrafo_2D(pT_dPhi, binedges_x, binedges_y);
	TH1D eff = calculateEfficiency(bothSel, genSel);
	for (int i=1; i <=Nbins; i++) {
		eff_final.SetBinContent(i, eff.GetBinContent(i));
		eff_final.SetBinError(i, eff.GetBinError(i));
	}
	eff_final.SetName("eff_final");
	eff_final.SetLineColor(kGreen+2);
	eff_final.SetFillColor(kGreen+2);
	eff_final.SetMarkerStyle(22);
	eff_final.SetMarkerSize(0.8);
	eff_final.SetMarkerColor(kGreen+2);
	eff_final.Draw("same pe1");
	
	TLine * aline = new TLine();
	TLatex * atext = new TLatex();
	atext->SetTextSize(0.03);
	aline->SetLineWidth(2);
	aline->SetLineColor(kBlack);

	auto legend = new TLegend(gPad->GetLeftMargin(),0.91,1-gPad->GetRightMargin(),1-gPad->GetTopMargin());
	legend->AddEntry("stab_final", "Stability", "lep");
	legend->AddEntry("pur_final", "Purity", "lep");
	legend->AddEntry("eff_final", "Efficiency", "lep");
	legend->AddEntry((TObject*)0, " ", "");
	TLegendEntry * entry;
	if (Puppi) entry = legend->AddEntry((TObject*)0, "Puppi", "");
    else entry = legend->AddEntry((TObject*)0, "DNN", "");
    entry->SetTextFont(62);
	legend->SetShadowColor(0);
	legend->SetBorderSize(1);
	legend->SetFillColor(kWhite);
	legend->SetNColumns(5);
	
	for (int i=1; i<=phiBins-1; i++){
		aline->DrawLine(i*400,lowerEnd,i*400, upperEnd);
	}
	//Set Label at 90% of y-axis
	for (int i=0;i < phiBins; i++){
		if (i==0) atext->DrawLatex(100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("|#Delta#phi|<%.1f",binedges_y.at(1))); 
		else {
			if (i==phiBins-1) atext->DrawLatex(400*i+100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|",binedges_y.at(i)));
			else atext->DrawLatex(400*i+100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|<%.1f",binedges_y.at(i), binedges_y.at(i+1)));
		}
	}
	legend->Draw();
    TLatex labelSPEPuppi = gfx::cornerLabel("PuppiMET", 3);
	TLatex labelSPEDNN = gfx::cornerLabel("DNN MET", 3);
	
    //~ if (Puppi) labelSPEPuppi.Draw();
    //~ else labelSPEDNN.Draw();
    
    saver->save(canvas, filename);
    if (savePDF) canvas.SaveAs(dir + filename.ReplaceAll("/", "_") + "1D.pdf");
    canvas.Close();
}

void PlotRMS(TH2D* const &pT_dPhi, TH2D* const &res, std::vector<float> binedges_x, std::vector<float> binedges_y, bool Puppi, bool RMSpT, bool savePDF, TString dir, TString filename, io::RootFileSaver* saver){
	TCanvas canvas;
	int Nbins = (binedges_x.size()-1)*(binedges_y.size()-1);
	int NbinsPt = binedges_x.size()-1;
	int NbinsDphi = binedges_y.size()-1;
	canvas.SetLogy(false);

	TH1D RMS = getRMS(res);
	TH1D RMS_final = HistTrafo_2D(pT_dPhi, binedges_x, binedges_y);
	for (int i=1; i <=Nbins; i++) {
		RMS_final.SetBinContent(i, RMS.GetBinContent(i));
		RMS_final.SetBinError(i, RMS.GetBinError(i));
	}
	RMS_final.SetStats(false);
	RMS_final.GetXaxis()->SetTitle("%pTnunu");
	gPad->SetLeftMargin(0.17);
	if (RMSpT) RMS_final.GetYaxis()->SetTitle("RMS(%MET-p_{#scale[.8]{T}}^{#scale[.8]{#nu #bar{#nu}}})");
	else RMS_final.GetYaxis()->SetTitle("RMS(#Delta#phi_{rec}-#Delta#phi_{gen})");
	RMS_final.SetLineColor(kBlack);
	auto legend_RMS = new TLegend(gPad->GetLeftMargin(),0.9,1-gPad->GetRightMargin(),1-gPad->GetTopMargin());
	TH1D* res_copy = (TH1D*) RMS_final.Clone("copy");
	res_copy->SetLineColor(kBlue+2);
	res_copy->SetLineStyle(2);
	//~ res_copy->SetLineWidth(RMS_final.GetLineWidth()*1.5);
	res_copy->SetMarkerSize(0);
	TLegendEntry * RMSentry;
	if (RMSpT) {
		for (int i=1; i<= res_copy->GetNbinsX(); i++){
			res_copy->SetBinContent(i, binedges_x.at((i-1)%NbinsPt +1) - binedges_x.at((i-1)%NbinsPt));
		}
		RMS_final.Draw("hist");
		res_copy->Draw("same hist");
		RMS_final.GetYaxis()->SetRangeUser(15, 220);
		   
		legend_RMS->AddEntry(&RMS_final, "RMS p_{#scale[.8]{T}}","l");
		legend_RMS->AddEntry(res_copy, "Bin Width p_{#scale[.8]{T}}","l");
	} else {
		int phiBin = 0;
		for (int i=1; i<= res_copy->GetNbinsX(); i++){
			res_copy->SetBinContent(i, binedges_y.at(phiBin +1) - binedges_y.at(phiBin));
			if (i%NbinsPt == 0) phiBin++;
		}
		RMS_final.Draw("hist");
		res_copy->Draw("same hist");
		RMS_final.GetYaxis()->SetRangeUser(0.05, 2.2);
		   
		legend_RMS->AddEntry(&RMS_final, "RMS #Delta#phi","l");
		legend_RMS->AddEntry(res_copy, "Bin Width #Delta#phi","l");
	}
	legend_RMS->AddEntry((TObject*)0, " ", "");
	if (Puppi) RMSentry = legend_RMS->AddEntry((TObject*)0, "PUPPI", "");
	else RMSentry = legend_RMS->AddEntry((TObject*)0, "DNN", "");
	RMSentry->SetTextFont(62);
	legend_RMS->SetShadowColor(0);
	legend_RMS->SetBorderSize(1);
	legend_RMS->SetFillColor(kWhite);
	legend_RMS->SetNColumns(4);
	legend_RMS->SetTextSize(0.03);


	int phiBins = binedges_y.size() - 1;
	RMS_final.GetXaxis()->SetNdivisions(8*phiBins, false);
	RMS_final.GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
	for (int i=0; i<=phiBins-1; i++){
		RMS_final.GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
		RMS_final.GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
		RMS_final.GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
		RMS_final.GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
		RMS_final.GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
		RMS_final.GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
		RMS_final.GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
		RMS_final.GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
	}
	
	TLine * aline = new TLine();
	TLatex * atext = new TLatex();
	atext->SetTextSize(0.03);
	aline->SetLineWidth(2);
	aline->SetLineColor(kBlack);
	float upperEnd=1.05*RMS_final.GetMaximum();
	float lowerEnd=0.9*RMS_final.GetMinimum();
	RMS_final.GetYaxis()->SetRangeUser(lowerEnd, upperEnd);

	for (int i=1; i<=phiBins-1; i++){
		aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
	}
	
	legend_RMS->Draw();
	for (int i=0;i < phiBins; i++){
		if (i==0) atext->DrawLatex(100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("|#Delta#phi|<%.1f",binedges_y.at(1))); 
		else {
			if (i==phiBins-1) atext->DrawLatex(400*i+100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|",binedges_y.at(i)));
			else atext->DrawLatex(400*i+100,0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|<%.1f",binedges_y.at(i), binedges_y.at(i+1)));
		}
	}
	//~ TLatex labelRMSPuppi = gfx::cornerLabel("PUPPI", 4);
	//~ TLatex labelRMSDNN = gfx::cornerLabel("DNN", 4);

	//~ if (Puppi) labelRMSPuppi.Draw();
	//~ else labelRMSDNN.Draw();

    saver->save(canvas, filename);
    if (savePDF) canvas.SaveAs(dir + filename.ReplaceAll("/", "_") + ".pdf");
    canvas.Close();
}


//minVal		Vector of threshold values
//newVal		Vector of values after the merging process
//oldVal		Vector of values before merging 
//widthFactor	Ratio of the bin width to the range of the corresponding axis
//MoIOption		States, which definition of measurement to use
float getMeasureOfImprovement(std::vector<float> minVal, std::vector<float> newVal, std::vector<float> oldVal, float widthFactor, int MoIOption) {
	float measure;
	widthFactor*=1000;
	widthFactor = round(widthFactor);
	widthFactor = widthFactor/1000;
	
	switch (MoIOption) {
		case 0: //Relative Improvement
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure *= 1 + (newVal.at(i) - oldVal.at(i)) / minVal.at(i);
			}
			return measure;
		case 1: //Relative Improvement and missing Value
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure *= 1 + (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / (pow(minVal.at(i), 2));
			}
			return measure;
		case 2: //Relative Improvement and missing Value and WidthFactor
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure *= 1 + (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / (pow(minVal.at(i), 2)*widthFactor);
			}
			return measure;
		case 3: //Relative Improvement and missing Value and WidthFactor at the end
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure *= 1 + (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / pow(minVal.at(i), 2);
			}
			return measure/widthFactor;
		case 4: //Relative Improvement and missing Value and WidthFactor at the end and only increase Allowed
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				if (newVal.at(i) < oldVal.at(i)) return 0;
				measure *= 1 + (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / pow(minVal.at(i), 2);
			}
			return measure/widthFactor;
		case 5: //Relative Improvement(max.contribution missing Val) and missing Value and WidthFactor at the end
			measure = 1;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure *= 1 + std::min(newVal.at(i) - oldVal.at(i), minVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / pow(minVal.at(i), 2);
			}
			return measure/widthFactor;
		//Additive Calculations
		case 10: //F-Score Approach (does not take inf and negative terms into account)
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure += (minVal.at(i) - oldVal.at(i))/(newVal.at(i) - oldVal.at(i));
			}
			return 1/(measure*widthFactor);
		case 11: //Same as 3 but as inverse sum (most probably not correct approach)
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure += (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / (pow(minVal.at(i), 2));
			}
			return 1/(measure*widthFactor);
		case 12: //Same as 3 but as sum
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure += (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / (pow(minVal.at(i), 2));	
			}
			return measure/widthFactor;
		case 13: //Same as 12 but only increase Allowed
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				if (newVal.at(i) < oldVal.at(i)) return 0;
				measure += (newVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / (pow(minVal.at(i), 2));	
			}
			return measure/widthFactor;
		case 14: //Same as 5 but as sum
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				measure += std::min(newVal.at(i) - oldVal.at(i), minVal.at(i) - oldVal.at(i))*(minVal.at(i) - oldVal.at(i)) / pow(minVal.at(i), 2);
			}
			return measure/widthFactor;
		case 15: //F-Score Approach (does take inf and negative terms into account)
			//If the new value is worse or equal than the old value, do not add anything. To compensate for this, the result ist weighted by the number of terms.
			int k=0;
			measure = 0;
			for (int i=0; i<minVal.size(); i++) {
				if (oldVal.at(i) >= minVal.at(i)) continue;
				if (newVal.at(i) <= oldVal.at(i)) continue;
				measure += (minVal.at(i) - oldVal.at(i))/(newVal.at(i) - oldVal.at(i));
				k++;
			}
			//In case measure=0 because all newVal are allready fulfilled or worse than the old value, return a very small MoI which still takes the widthFactor into account
			if (measure==0) return 0.0001/widthFactor;
			else return k/(measure*widthFactor);
	}
	return 0;
}

//Determines the next merging direction. See flowchart in bachelor thesis

//Returns 			-1 when minVals are already fulfilled and RMS pT/dPhi is ok
//					0 when merging in pT and some minVals are not fulfilled or the RMS pT is too great
//					1 when merging in pT and minVals are fulfilled and RMS pT/dPhi is ok
//					10 when merging in dPhi and some minVals are not fulfilled or the RMS dPhi is too great
//					11 when merging in dPhi and minVals are fulfilled and RMS pT/dPhi is ok

//RecoGen 			Hist of RecoGen
//ResPt				Histogram of Pt Resolution
//ResPhi			Histogram of dPhi Resolution
//bothSelection		Histogram of pTnunu with selection on generator and reconstruction level
//genSelection		Histogram of pTnunu with selection on generator level
//currentBin		Current Bin in the algorithm 
//i					Number of pT bins
//*_min				Threshold values for the binning
//binWidthPt/Phi 	Bin Width of the current bin
//binWidthPt_PT		Bin Width of the current pT bin when merged with the next pT Bin
//binWidthPhi_PHI	Bin Width of the current dPhi bin when merged with the next dPhi Bin
//rangePt/Phi		Full range of the pT/dPhi axis
//MoI				Denotes which measure of improvement to use
//logger			Logging object to write into the log file
//is_logging		True, when logging is enabled
//edge_40GeV		True, when it should be tried to set a bin edge at 40GeV	(currently not used)
//merge_Overflow	True, when overflow bins should be considered				(currently not used)
//acMoI				true if always the MoI should be considered
//lastBinMergeRule	true, if the current bin is the last bin on one axis, only merges are performed on the other axis (Including Check if N_rec can be fulfilled with subsequent bins)
//counter			Used determine the merging when N=0 in a bin and its surronding
int getNextDirection(TH2D* const RecoGenHist, TH2D* const ResPtHist , TH2D* const ResPhiHist,  TH1D* const bothSelectionHist, TH1D* const genSelectionHist, int currentBin, int i, float p_min, float s_min, float N_min, float eff_min, float binWidthPt, float binWidthPhi, float binWidthPt_PT, float binWidthPhi_PHI ,float rangePt, float rangePhi, int MoI, ofstream &logger, bool is_logging, bool edge_40GeV, bool merge_Overflow, bool acMoI, bool lastBinMergeRule, int counter){
	int n = currentBin;
	int m = 1;
	while (n > i) {
		n -= i;
		m++;
	}
	std::vector<float> minVal = {p_min, s_min, N_min, eff_min};
	
	//Current Binning
	TH2D RecoGen = *RecoGenHist;
	TH2D ResPt = *ResPtHist;
	TH2D ResPhi = *ResPhiHist;
	TH1D bothSelection = *bothSelectionHist;
	TH1D genSelection = *genSelectionHist;
	
	TH1D purity = calculatePurity(&RecoGen);
	TH1D stability = calculateStability(&RecoGen);
	TH1D efficiency = calculateEfficiency(&bothSelection, &genSelection);

	float p = purity.GetBinContent(currentBin);
	float s = stability.GetBinContent(currentBin);
	float eff = efficiency.GetBinContent(currentBin);
	float N_rec = RecoGen.ProjectionY(("yProj" + std::to_string(counter)).c_str(), 1, RecoGen.GetNbinsY())->GetBinContent(currentBin);
	float RMS_pT = getRMS(&ResPt).GetBinContent(currentBin);
	float RMS_Phi = getRMS(&ResPhi).GetBinContent(currentBin);
	std::vector<float> oldVal = {p, s, N_rec, eff};
	
	if(is_logging) logger << "Current values in getNextDirection: currentBin=" + std::to_string(currentBin)   +  ", N_rec=" + std::to_string(N_rec)+ ", purity=" + std::to_string(p)  + ", stability=" + std::to_string(s)  + ", efficiency=" + std::to_string(eff) + "\n";
	if(is_logging) logger << "\tCheck RMS values: RMS pT=" + std::to_string(RMS_pT) + ", bin width pT=" + std::to_string(binWidthPt) + "\t, RMS dPhi=" + std::to_string(RMS_Phi) + ", bin width dPhi=" + std::to_string(binWidthPhi) + "\n";
	
	//Check the binWidth (counter dependent which merge is checked first)
	if(counter%2 == 0) {
		if (RMS_Phi >= binWidthPhi) {
			if(is_logging) logger << "->dPhi Merge. The RMS dPhi is too great.\n";
			return 10;
		}
		if (RMS_pT >= binWidthPt) {
			if(is_logging) logger << "->pT Merge. The RMS pT is too great.\n";
			return 0;
		}
	} else {
		if (RMS_pT >= binWidthPt) {
			if(is_logging) logger << "->pT Merge. The RMS pT is too great.\n";
			return 0;
		}
		if (RMS_Phi >= binWidthPhi) {
			if(is_logging) logger << "->dPhi Merge. The RMS dPhi is too great.\n";
			return 10;
		}
	}
	
	if(is_logging) logger << "\tRMS values of current bin are ok.\n";
	if(s>= s_min and p >= p_min and N_rec >= N_min and eff >= eff_min) {
		if(is_logging) logger << "->No Merge, current Bin is ok. The minVals are already fulfilled.\n";
		return -1;
	}
	
	//Calculate the indices for the merging steps
	//ATTENTION: m and n are no longer corresponding to the currentBin
	//If it is the last pT bin, the merging can only be executed with the next to last bin.
	int currentBin_PT = currentBin - m + 1;
	bool lastPtBin = false;
	if (n == i) {
		if(is_logging) logger << "\tLast pT bin: n=" + std::to_string(n) + " . Merging it with the next to last pT bin.\n";
		n -= 1;
		currentBin_PT -= 1;
		lastPtBin = true;
	}
	
	//If it is the last pT bin, the merging can only be executed with the next to last bin.
	int currentBin_PHI = currentBin;
	bool lastPhiBin = false;
	if (currentBin > RecoGen.GetNbinsX()-i) {
		if(is_logging) logger << "\tLast dPhi bin: m=" + std::to_string(m) + " . Merging it with the next to last dPhi bin.\n";
		m -= 1;
		currentBin_PHI -= i;
		lastPhiBin = true;
	}
	//Merging in both directions
	//Merge in pT direction
	TH2D RecoGen_PT = mergePtRecoGen(*RecoGenHist, i, n);
	TH2D ResPt_PT = mergePtRes(*ResPtHist, i, n);
	TH2D ResPhi_PT = mergePtRes(*ResPhiHist, i, n);
	TH1D bothSelection_PT = mergePt1D(*bothSelectionHist, i, n);
	TH1D genSelection_PT = mergePt1D(*genSelectionHist, i, n);
	RecoGen_PT.SetName(("RG_PT"+ std::to_string(counter)).c_str());
	ResPt_PT.SetName(("ResPt_PT"+ std::to_string(counter)).c_str());
	ResPhi_PT.SetName(("ResPhi_PT"+ std::to_string(counter)).c_str());
	bothSelection_PT.SetName(("bS_PT"+ std::to_string(counter)).c_str());
	genSelection_PT.SetName(("gS_PT"+ std::to_string(counter)).c_str());
	
	TH1D purity_PT = calculatePurity(&RecoGen_PT);
	TH1D stability_PT = calculateStability(&RecoGen_PT);
	TH1D efficiency_PT = calculateEfficiency(&bothSelection_PT, &genSelection_PT);
	float p_PT = purity_PT.GetBinContent(currentBin_PT);
	float s_PT = stability_PT.GetBinContent(currentBin_PT);
	float eff_PT = efficiency_PT.GetBinContent(currentBin_PT);
	float N_rec_PT = RecoGen_PT.ProjectionY(("yProjPT"+ std::to_string(counter)).c_str(), 1, RecoGen_PT.GetNbinsY())->GetBinContent(currentBin_PT);
	float RMS_pT_PT = getRMS(&ResPt_PT).GetBinContent(currentBin_PT);
	float RMS_Phi_PT = getRMS(&ResPhi_PT).GetBinContent(currentBin_PT);
	std::vector<float> newVal_PT = {p_PT, s_PT, N_rec_PT, eff_PT};
	
	//Merge in dPhi direction
	TH2D RecoGen_PHI = mergePhiRecoGen(*RecoGenHist, i, m);
	TH2D ResPt_PHI = mergePhiRes(*ResPtHist, i, m);
	TH2D ResPhi_PHI = mergePhiRes(*ResPhiHist, i, m);
	TH1D bothSelection_PHI = mergePhi1D(*bothSelectionHist, i, m);
	TH1D genSelection_PHI = mergePhi1D(*genSelectionHist, i, m);

	TH1D purity_PHI = calculatePurity(&RecoGen_PHI);
	TH1D stability_PHI = calculateStability(&RecoGen_PHI);
	TH1D efficiency_PHI = calculateEfficiency(&bothSelection_PHI, &genSelection_PHI);

	float p_PHI = purity_PHI.GetBinContent(currentBin_PHI);
	float s_PHI = stability_PHI.GetBinContent(currentBin_PHI);
	float eff_PHI = efficiency_PHI.GetBinContent(currentBin_PHI);
	float N_rec_PHI = RecoGen_PHI.ProjectionY(("yProjPhi"+ std::to_string(counter)).c_str(), 1, RecoGen_PHI.GetNbinsY())->GetBinContent(currentBin_PHI);
	float RMS_pT_PHI = getRMS(&ResPt_PHI).GetBinContent(currentBin_PHI);
	float RMS_Phi_PHI = getRMS(&ResPhi_PHI).GetBinContent(currentBin_PHI);
	std::vector<float> newVal_PHI = {p_PHI, s_PHI, N_rec_PHI, eff_PHI};	
	
	if(is_logging) logger << "\tMerging the nth pT Bin:\t\t n=" + std::to_string(n)   +  ", N_rec=" + std::to_string(N_rec_PT)+ ", purity=" + std::to_string(p_PT)  + ", stability=" + std::to_string(s_PT)  + ", efficiency=" + std::to_string(eff_PT) + ", RMS pT=" + std::to_string(RMS_pT_PT) + ", RMS dPhi=" + std::to_string(RMS_Phi_PT) + ", new bin width pT=" + std::to_string(binWidthPt_PT) + "\n";
	
	if(is_logging) logger << "\tMerging the mth dPhi Bin:\t m=" + std::to_string(m)   +  ", N_rec=" + std::to_string(N_rec_PHI)+ ", purity=" + std::to_string(p_PHI)  + ", stability=" + std::to_string(s_PHI)  + ", efficiency=" + std::to_string(eff_PHI) + ", RMS pT=" + std::to_string(RMS_pT_PHI) + ", RMS dPhi=" + std::to_string(RMS_Phi_PHI) + ", new bin width dPhi=" + std::to_string(binWidthPhi_PHI) + "\n";
	
	//Check if params are ok
	bool ok_PT = s_PT >= s_min and p_PT >= p_min and N_rec_PT >= N_min and eff_PT >= eff_min and RMS_pT_PT<binWidthPt_PT and RMS_Phi_PT<binWidthPhi;
	bool ok_PHI = s_PHI >= s_min and p_PHI >= p_min and N_rec_PHI >= N_min and eff_PHI >= eff_min and RMS_Phi_PHI<binWidthPhi_PHI and RMS_pT_PHI<binWidthPt;
	
	//If last Bin merge rule is etablished, additional checks have to be performed
	if (lastBinMergeRule) {
		//Check if N_rec requirement can be reached, when currently at the last bin of pT OR dPhi
		//Case: CurrentBin is last pT bin thus it has to be checked, if phi merges are able to establish N_rec > N_min
		if(lastPtBin and (not lastPhiBin)){
			//Calculate N_rec for the current bin added to the sum of N_Rec in the higher dPhi bins with fixed last pT bin
			TH1D* Stat = RecoGen.ProjectionY("RecoProj", 1, RecoGen.GetNbinsY());
			float N_recMax = 0;
			for (int term = currentBin; term <= Stat->GetNbinsX(); term += i) {
				N_recMax += Stat->GetBinContent(term);
				//~ std::cout<<term<<std::endl;
			}
			if (N_recMax < N_min) {
				if (ok_PT) {
					if(is_logging) logger << "->pT Merge because current bin is the last pT bin but N_min can not be fulfilled with the following dPhi bins. The minVals are fulfilled and RMS is ok. N_recMax=" + std::to_string(N_recMax) + "\n";
					return 1;
				} else {
					if(is_logging) logger << "->pT Merge because current bin is the last pT bin but N_min can not be fulfilled with the following dPhi bins. The minVals are not fulfilled. N_recMax=" + std::to_string(N_recMax) + "\n";
					return 0;
				}
			}
		}
		//Case: CurrentBin is last dPhi bin thus it has to be checked, if pT merges are able to establish N_rec > N_min
		if(lastPhiBin and (not lastPtBin)){
			//Calculate N_rec for the current bin added to the sum of N_Rec in the higher pT bins with fixed last pT bin
			TH1D* Stat = RecoGen.ProjectionY("RecoProj", 1, RecoGen.GetNbinsY());
			float N_recMax = 0;
			//~ std::cout<<currentBin<<std::endl;
			//~ std::cout<<i*m<<std::endl;
			for (int term = currentBin; term <= i*(m+1); term++) {
				N_recMax += Stat->GetBinContent(term);
				//~ std::cout<<term<<std::endl;
			}
			if (N_recMax < N_min) {
				if (ok_PHI) {
					if(is_logging) logger << "->dPhi Merge because current bin is the last dPhi bin but N_min can not be fulfilled with the following pT bins. The minVals are fulfilled and RMS is ok. N_recMax=" + std::to_string(N_recMax) + "\n";
					return 11;
				} else {
					if(is_logging) logger << "->dPhi Merge because current bin is the last dPhi bin but N_min can not be fulfilled with the following pT bins. The minVals are not fulfilled. N_recMax=" + std::to_string(N_recMax) + "\n";
					return 10;
				}
			}
		}
		
		
		//If the algorithm is at the last pT bin, merge only in dPhi unless it is also the last dPhi bin and the same in the oppposite case
		//Case: Last pT bin but not last dPhi bin -> Merge in dPhi
		if(lastPtBin and (not lastPhiBin)) {
			if (ok_PHI) {
				if(is_logging) logger << "->dPhi Merge because current bin is the last pT bin. The minVals are fulfilled and RMS is ok.\n";
				return 11;
			} else {
				if(is_logging) logger << "->dPhi Merge because current bin is the last pT bin. The minVals are not fulfilled.\n";
				return 10;
			}
		}
		//Case: Last dPhi bin but not last pT bin -> Merge in pT
		if(lastPhiBin and (not lastPtBin)) {
			if (ok_PT) {
				if(is_logging) logger << "->pT Merge because current bin is the last dPhi bin. The minVals are fulfilled and RMS is ok.\n";
				return 1;
			} else {
				if(is_logging) logger << "->pT Merge because current bin is the last dPhi bin. The minVals are not fulfilled.\n";
				return 0;
			}
		}
	}
	
	//Case: N_rec = 0 in the surronding bins, choose accordingly to the counter value		
	if(N_rec_PT==0 and N_rec_PHI==0) {
		if (counter%2 == 0) {
			if(is_logging) logger << "->dPhi Merge because N_rec=0 in the surronding bins and counter=" + std::to_string(counter) + ".\n";
			return 10;
		} else {
			if(is_logging) logger << "->pT Merge because N_rec=0 in the surronding bins and counter=" + std::to_string(counter) + ".\n";
			return 0;
		}
	}
	
	//Other cases are handled by comparing the measure of improvement
	
	//Calculate Measure of Improvement		Currently decrease in values is allowed
	float MoI_PT = getMeasureOfImprovement(minVal, newVal_PT, oldVal, (binWidthPt_PT-binWidthPt)/rangePt, MoI);
	float MoI_PHI= getMeasureOfImprovement(minVal, newVal_PHI, oldVal, (binWidthPhi_PHI-binWidthPhi)/rangePhi, MoI);
	if(is_logging) logger << "\tMeasure of Improvement: For Merging in pT: " + std::to_string(MoI_PT) + "\t For Merging in dPhi: " + std::to_string(MoI_PHI) + "\n";
	
	
	if (acMoI) {
		//Case: MoI equal
		if(MoI_PT == MoI_PHI) {
			if (counter%2 == 0) {
				if(is_logging) logger << "->dPhi Merge because both measures are equal and counter=" + std::to_string(counter) + ".\n";
				return 10;
			} else {
				if(is_logging) logger << "->pT Merge because both measures are equal and counter=" + std::to_string(counter) + ".\n";
				return 0;
			}
		}
		if (MoI_PT < MoI_PHI) {
			if (ok_PHI) {
				if(is_logging) logger << "->dPhi Merge. The minVals are fulfilled and RMS is ok.\n";
				return 11;
			} else {
				if(is_logging) logger << "->dPhi Merge. The minVals are not fulfilled.\n";
				return 10;
			}
		} else {
			if (ok_PT) {
				if(is_logging) logger << "->pT Merge. The minVals are fulfilled and RMS is ok.\n";
				return 1;
			} else {
				if(is_logging) logger << "->pT Merge. The minVals are not fulfilled.\n";
				return 0;
			}
		}
	} else {
		//Case: At least one merging direction fulfilles the requirements
		if(ok_PT or ok_PHI) {
			if(ok_PT and ok_PHI) {
				//Case: Both merges yield satisfying results -> use the one with greater MoI
				if (MoI_PT <= MoI_PHI) {
					if(is_logging) logger << "->dPhi Merge. The minVals are fulfilled and RMS is ok. Also a merge in pT would have been sufficient.\n";
					return 11;
				} else {
					if(is_logging) logger << "->pT Merge. The minVals are fulfilled and RMS is ok. Also a merge in dPhi would have been sufficient.\n";
					return 1;
				}
			} else {
				//Case: only one direction fulfilles all requirements
				if (ok_PT) {
					if(is_logging) logger << "->pT Merge. The minVals are fulfilled and RMS is ok.\n";
					return 1;
				}
				if (ok_PHI) {
					if(is_logging) logger << "->dPhi Merge. The minVals are fulfilled and RMS is ok.\n";
					return 11;
				}
			}
		}
		
		//Case: MoI equal
		if(MoI_PT == MoI_PHI) {
			if (counter%2 == 0) {
				if(is_logging) logger << "->dPhi Merge because both measures are equal and counter=" + std::to_string(counter) + ".\n";
				return 10;
			} else {
				if(is_logging) logger << "->pT Merge because both measures are equal and counter=" + std::to_string(counter) + ".\n";
				return 0;
			}
		}
		if (MoI_PT < MoI_PHI) {
			if(is_logging) logger << "->dPhi Merge. The minVals are not fulfilled.\n";
			return 10;
		} else {
			if(is_logging) logger << "->pT Merge. The minVals are not fulfilled.\n";
			return 0;
		}
	}
}

Config const &cfg=Config::get();

extern "C"

// Run
int run()
{	// Loop over varioues options to test them subsequently without running the scribt for each set of options separatly
	//The varied options are 'forward/backward', 'nextBinPt/nextBinPhi' and '(n)acMoI' 
	std:: vector<vector<bool>> options
    {
        {false, false, false},
        {false, false, true},
        {false, true, false},
        {false, true, true},
        {true, false, false},
        {true, false, true},
        {true, true, false},
        {true, true, true},
    };
    std::vector<int> MoI_IDs = {3, 4, 5, 12, 13, 14, 15};
    
    //Used optimized binnings in the thesis
    //N=25
	//~ std:: vector<vector<bool>> options
    //~ {
        //~ {true, true, false},
    //~ };
    //~ //N=10
	// ~std:: vector<vector<bool>> options{{true, false, true},};
	
	//Used MoI version
    // ~std::vector<int> MoI_IDs = {5};
	bool is_logging = true;
	
	
	int initialBins = 40; //number of equidistant initial bins in each dimension
	// ~int initialBins = 80; //number of equidistant initial bins in each dimension
	float p_min = 0.2;
	float s_min = 0.2;
	float N_min = 10;
	float eff_min = 0.2;
	// ~bool merge_Overflow = false; //Currently not used
	bool merge_Overflow = true; //Currently not used
	bool consider_cut = false;	 //Currently not used
	// ~bool DNN = true;		//If false, PUPPI is used
	bool DNN = false;		//If false, PUPPI is used
	//Loop over sets of options
	for (int r = 0; r < MoI_IDs.size(); r++){
		int MoI = MoI_IDs.at(r);		//See getMeasurementOfImprovement for an explanation of the values
		
		//Two Logs are produced, one comprises a detailed description of the algorithm steps for a certain set of options. (See below)
		//The other Log contains the optimized binnings for each set of options which is processed in the script (See here)
		TString dir = cfg.outputDirectory+"/binningOpt/2D/";
		io::ensurePathForFile(dir);
		ofstream resultLogger;
		TString fileAppendixResults= "2D0A_";
		fileAppendixResults += "iB" + std::to_string(initialBins) + "_";
		fileAppendixResults += "MoI" + std::to_string(MoI) + "_";
		if (merge_Overflow) fileAppendixResults += "mo_";
		if (DNN) fileAppendixResults +=  "DNN_";
		else fileAppendixResults +=  "Puppi_";
		if (not (N_min == 200)) fileAppendixResults += "N" + std::to_string(int(N_min)) + "_";
		if(is_logging) {
			resultLogger.open(dir + fileAppendixResults + "Log.txt", ios::trunc);
		}
		
		//Do the calculations for each option defined in 'options'
		for(int k = 0; k < options.size(); k++) {
			
			bool forward = options.at(k).at(0);
			bool nextBinPt = options.at(k).at(1); 		//true if, after a bin fulfilles the requirements, the next bin is the next pT Bin. If false, the next bin ist the next dPhi bin.
			bool alwaysConsiderMoI = options.at(k).at(2);
			bool lastBinMergeRule = not alwaysConsiderMoI;
			
			TString treeName="TTbar_diLepton";     //input sample name
			TFile file(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/100.0/Nominal/"+treeName+"_merged.root",cfg.year.Data(),cfg.treeVersion.Data()),"read");    //open input file
			TTreeReader reader("ttbar_res100.0/"+treeName, &file);      //open tree in root file

			//define histograms
			hist::Histograms<TH2D> hs2D({treeName});
			hs2D.addHist("PuppiMET/Pt_DeltaPhi"   ,"PuppiPD;%MET;|#Delta#phi|(p_{#scale[.8]{T}}^{#scale[.8]{miss}}, nearest l);EventsBIN"           ,initialBins, 0, 400, initialBins, 0, 3.2);
			hs2D.addHist("PuppiMET/RecoGen"   ,"PuppiRG;Gen BinNumber;Reco BinNumber;EventsBIN"           ,pow(initialBins, 2), 0, pow(initialBins, 2), pow(initialBins, 2), 0, pow(initialBins, 2));
			hs2D.addHist("PuppiMET/Resolution_pT"   ,"PuppiResPt;Reco Bin Number;%MET - p_{#scale[.8]{T}}(#nu #bar{#nu});EventsBIN"           , pow(initialBins, 2),0,pow(initialBins, 2),160,-400,400);
			hs2D.addHist("PuppiMET/Resolution_dPhi"   ,"PuppiResdPhi;Reco Bin Number;#Delta#phi_{rec}-#Delta#phi_{gen};EventsBIN"           , pow(initialBins, 2),0,pow(initialBins, 2),160,-M_PI,M_PI);

			hs2D.addHist("DNN_MET/RecoGen"   ,"DNNRG;Gen BinNumber;Reco BinNumber;EventsBIN"           ,pow(initialBins, 2), 0, pow(initialBins, 2),pow(initialBins, 2), 0, pow(initialBins, 2));
			hs2D.addHist("DNN_MET/Pt_DeltaPhi"   ,"DNNPD;%MET;|#Delta#phi|(p_{#scale[.8]{T}}^{#scale[.8]{miss}}, nearest l);EventsBIN"           ,initialBins, 0, 400, initialBins, 0, 3.2);
			hs2D.addHist("DNN_MET/Resolution_pT"   ,"DNNRespT;Reco Bin Number;%MET - p_{#scale[.8]{T}}(#nu #bar{#nu});EventsBIN"           , pow(initialBins, 2),0,pow(initialBins, 2),160,-400,400);
			hs2D.addHist("DNN_MET/Resolution_dPhi"   ,"DNNResPhi;Reco Bin Number;#Delta#phi_{rec}-#Delta#phi_{gen};EventsBIN"           , pow(initialBins, 2),0,pow(initialBins, 2),160,-M_PI,M_PI);

			hist::Histograms<TH1D> hs_nW({treeName});
			hs_nW.addHist("PtNuNu/bothSel_nW", "NuNubS;%pTnunu;Events/Bin", pow(initialBins, 2), 0, pow(initialBins, 2));
			hs_nW.addHist("PtNuNu/genSel_nW", "NuNugS;%pTnunu;Events/Bin", pow(initialBins, 2), 0, pow(initialBins, 2));
		   
			std::vector<float> edges_pT;
			edges_pT.push_back(0);
			std::vector<float> edges_dPhi;
			edges_dPhi.push_back(0);
			for (int i = 1; i<=initialBins; i++) {
				edges_pT.push_back(hs2D.getHistogram("PuppiMET/Pt_DeltaPhi", treeName)->GetXaxis()->GetBinUpEdge(i));
				edges_dPhi.push_back(hs2D.getHistogram("PuppiMET/Pt_DeltaPhi", treeName)->GetYaxis()->GetBinUpEdge(i));
			}
			//Reverse the vektor if the direction is backwards
			if (not forward) {
				std::reverse(edges_pT.begin(),edges_pT.end()); 
				std::reverse(edges_dPhi.begin(),edges_dPhi.end()); 
			}
			
			//Logger for detailed algorithm steps. The abbrevitions denote the following options
			//2DOA						Identifier for 2D optimized binnings
			//fw/bw						Forward/Backward direction
			//PtBin/PhiBin				Next bin after the current bin is completed
			//iB						Number of equidistant initial bins in each dimension
			//MoI						Used MoI version (See getMeasureOfImprovement(...) used MoI IDs)
			//nmo/mo					Not used (merging overflow)
			//cc/ncc					Not used (consider cut)
			//acMoI/nacMoI				(not) always considering MoI
			//nlbmr/lbmr				last bin merging rule, is currently true if nacMoI is used. If the current bin is the last in a dimension, merging steps in the other are preferred with lbmr
			//DNN/PUPPI					DNN or PUPPI is used
			//Nx						minimum required number of reconstructed events. If not present in a file name, N_min = 200 holds.
			
			//Build file appendix
			TString fileAppendix= "2D0A_";
			if (forward) fileAppendix += "fw_";
			else fileAppendix += "bw_";
			if (nextBinPt) fileAppendix += "PtBin_";
			else fileAppendix += "PhiBin_";
			fileAppendix += "iB" + std::to_string(initialBins) + "_";
			fileAppendix += "MoI" + std::to_string(MoI) + "_";
			if (merge_Overflow) fileAppendix += "mo_";
			else fileAppendix += "nmo_";
			if (consider_cut) fileAppendix += "cc_";
			else fileAppendix += "ncc_";
			if (alwaysConsiderMoI) fileAppendix += "acMoI_";
			else fileAppendix += "nacMoI_";
			if (lastBinMergeRule) fileAppendix += "lbmr_";
			else fileAppendix += "nlbmr_";
			if (DNN) fileAppendix +=  "DNN_";
			else fileAppendix +=  "Puppi_";
			if (not (N_min == 200)) fileAppendix += "N" + std::to_string(int(N_min)) + "_";
			
			if(is_logging) std::cout<<"Logging is performed! See " + dir + fileAppendix + ".txt for the log."<<std::endl;
			else std::cout<<"No Logging is performed."<<std::endl;
			
			ofstream logger;
			if(is_logging) {
				logger.open(dir + fileAppendix + "Log.txt", ios::trunc);
				logger << "START\n";
				logger << "---------------- Parameters ----------------\n";
				if (DNN) logger << "DNN MET is used\n";
				else logger << "PuppiMET ist used\n";
				if (forward) logger << "Operating Mode: Forward\n" ;
				else logger << "Operating Mode: Backward\n" ;
				if (nextBinPt) logger << "Next Bin: pT\n" ;
				else logger << "NextBin: dPhi\n" ;
				logger << "Measure of improvement: " + std::to_string(MoI) + "\n";
				logger << "Initial Bins (each axis): " + std::to_string(initialBins) + "\n";
				logger << "Initial pT bin width: " + std::to_string(400/initialBins) + " GeV\n";
				logger << "Initial dPhi bin width: " + std::to_string(M_PI/initialBins) + "\n";
				logger << "Threshold Values: \n";
				logger << "\tPurity >= " + std::to_string(p_min) + " \n";
				logger << "\tStability >= " + std::to_string(s_min) + " \n";
				logger << "\tEfficiency >= " + std::to_string(eff_min) + " \n";
				logger << "\tN_rec >= " + std::to_string(N_min) + " \n";
				if (merge_Overflow) logger << "The overflow bin is merged.\n" ;
				else logger << "The overflow bin is not merged.\n" ;
				logger << "--------------------------------------------\n";
			}
		 
			
			//define variables taken from tree
			TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
			TTreeReaderValue<float> PuppiMET   (reader, "PuppiMET_xy");
			TTreeReaderValue<float> DNN_MET   (reader, "DNN_MET_pT");
			TTreeReaderValue<float> Phi_Puppi  (reader, "Phi_recPuppi_xy");
			TTreeReaderValue<float> Phi_NuNu   (reader, "Phi_NuNu");
			TTreeReaderValue<float> Phi_DNN   (reader, "DNN_MET_dPhi_nextLep");
			TTreeReaderValue<float> MC_weight   (reader, "N");
			TTreeReaderValue<float> SF_weight   (reader, "SF");
			
				//define booleans for selection
			bool rec_selectionPuppi;
			bool rec_selectionDNN;
			bool gen_selection;
			
			//set current sample for hist selection
			hs2D.setCurrentSample(treeName);
			hs_nW.setCurrentSample(treeName);
			
			//loop over events in tree
			int processEvents=cfg.processFraction*reader.GetEntries(true);
			int iEv=0;

			//~ int cc = 0;
			while (reader.Next()){
				iEv++;
				if (iEv>processEvents) break;
				if (iEv%(std::max(processEvents/10,1))==0){      //logging stuff
				 io::log*".";
				 io::log.flush();
				}
			  
				//set weight for current event
				hs2D.setFillWeight((*MC_weight)*(*SF_weight));

				// set booleans for selection
				gen_selection=false;
				rec_selectionPuppi=false;
				rec_selectionDNN=false;
				if(*PtNuNu>0 and (merge_Overflow or *PtNuNu <= 400)) gen_selection=true;
				if(*DNN_MET>0 and (DNN and (merge_Overflow or *DNN_MET <= 400))) rec_selectionDNN=true;
				if(*PuppiMET>0 and (not DNN and (merge_Overflow or *PuppiMET <= 400))) rec_selectionPuppi=true;
				
				int column = 1;
				int rowPuppi = 1;
				int rowDNN = 1;
				
				if (gen_selection) {
					column = getIndex(*Phi_NuNu, *PtNuNu, edges_dPhi, edges_pT, forward);
					hs_nW.fill("PtNuNu/genSel_nW", column-0.5);
				}
				
				if (gen_selection and rec_selectionPuppi) {
					//GenMET
					hs_nW.fill("PtNuNu/bothSel_nW", column-0.5);
					
					//PuppiMET
					hs2D.fill("PuppiMET/Pt_DeltaPhi", *PuppiMET, *Phi_Puppi);
					rowPuppi = getIndex(*Phi_Puppi, *PuppiMET, edges_dPhi, edges_pT, forward);
					hs2D.fill("PuppiMET/RecoGen", column-0.5, rowPuppi-0.5);
					hs2D.fill("PuppiMET/Resolution_pT", rowPuppi-0.5, *PuppiMET-*PtNuNu);
					hs2D.fill("PuppiMET/Resolution_dPhi", rowPuppi-0.5, *Phi_Puppi-*Phi_NuNu);
					
				}
				
				if (gen_selection and rec_selectionDNN) {
					//GenMET
					hs_nW.fill("PtNuNu/bothSel_nW", column-0.5);
					
					//DNN MET
					hs2D.fill("DNN_MET/Pt_DeltaPhi", *DNN_MET, *Phi_DNN);
					rowDNN = getIndex(*Phi_DNN, *DNN_MET, edges_dPhi, edges_pT, forward);
					hs2D.fill("DNN_MET/RecoGen", column-0.5, rowDNN-0.5);
					hs2D.fill("DNN_MET/Resolution_pT", rowDNN-0.5, *DNN_MET-*PtNuNu);
					hs2D.fill("DNN_MET/Resolution_dPhi", rowDNN-0.5, *Phi_DNN-*Phi_NuNu);
				}	
			}
			io::log<<"";
			file.Close();

			//Merge overflow bins
			if (merge_Overflow) {
			   hs2D.mergeOverflow();
			   hs_nW.mergeOverflow();
			}

			

			io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"binningUnfolding2D_Fabian_OA");
			TCanvas can;
			can.SetLogz();
			for(TString var : hs2D.getVariableNames()){
				TH2D* tempHist=hs2D.getHistogram(var,{treeName});
				tempHist->SetStats(false);
				tempHist->GetZaxis()->SetTitle("Events/BIN");
				//~ gPad->SetRightMargin(0.2);
				if (var.Contains("RecoGen") or var.Contains("Resolution")) can.SetLogz(false);
				else can.SetLogz(true);
				tempHist->Draw("colz");
				if (var.Contains("RecoGen")) {
					TLine* line = new TLine();
					line->SetLineWidth(2);
					line->SetLineColor(kGreen);
					line->SetLineStyle(2);
					for (int i=1; i < initialBins; i++) {
						line->DrawLine(i*(edges_pT.size()-1), 0, i*(edges_pT.size()-1),pow(initialBins, 2));
						line->DrawLine(0, i*( edges_pT.size()-1), pow(initialBins, 2), i*(edges_pT.size()-1));
					}
				}
				saver.save(can,var);
			}
			
			//currentBin is the bin index in the hists of the current bin
			//i is the current number of pT Bins
			//m is the current dPhi bin
			//n is the current pT bin
			//counter is an int which ensures that the rebinned hists name's are unique and is used for arbitrarily merging decisions
			int currentBin = 1;
			int i = initialBins;
			int m = 1;
			int n = 1;
			int counter = 1;
			
			TH2D* RecoGen;
			TH2D* ResPt;
			TH2D* ResPhi;
			TH1D* bothSel;
			TH1D* genSel; 
			
			if (DNN) {
				RecoGen = hs2D.getHistogram("DNN_MET/RecoGen", treeName);
				ResPt = hs2D.getHistogram("DNN_MET/Resolution_pT", treeName);
				ResPhi = hs2D.getHistogram("DNN_MET/Resolution_dPhi", treeName);
				bothSel = hs_nW.getHistogram("PtNuNu/bothSel_nW", treeName);
				genSel = hs_nW.getHistogram("PtNuNu/genSel_nW", treeName);
			} else {
				RecoGen = hs2D.getHistogram("PuppiMET/RecoGen", treeName);
				ResPt = hs2D.getHistogram("PuppiMET/Resolution_pT", treeName);
				ResPhi = hs2D.getHistogram("PuppiMET/Resolution_dPhi", treeName);
				bothSel = hs_nW.getHistogram("PtNuNu/bothSel_nW", treeName);
				genSel = hs_nW.getHistogram("PtNuNu/genSel_nW", treeName);
			}	
			
			//used range in both dimensions
			float rangePt = std::abs(edges_pT.at(edges_pT.size()-1) - edges_pT.at(0));
			float rangePhi = std::abs(edges_dPhi.at(edges_dPhi.size()-1) - edges_dPhi.at(0));
			
			//Breaks out of the loop when algorithm ends (I'm aware that this is not the best way to it)
			// See the flowchart for a clear overview
			while (true) {
				logger.flush();
				//Check if the binning consists of only one bin in pT or dPhi. If so, stop the algorithm
				if (edges_pT.size()==2 or edges_dPhi.size()==2) {
					if(is_logging) logger << "Break, only one bin left in one quantity.";
					break;
				}
				std::cout<<"-----------------------------------------------------------------" <<std::endl;
				std::cout<<"Begin while, cB=" + std::to_string(currentBin) + ", i=" + std::to_string(i)+ ", m=" + std::to_string(m)+ ", n=" + std::to_string(n) + ", #Bins=" + std::to_string(RecoGen->GetNbinsX())<<std::endl;
				if(is_logging){ 
					logger << "Start cycle with properties currentBin=" + std::to_string(currentBin) + ", i=" + std::to_string(i) + ", m=" + std::to_string(m) + ", n=" + std::to_string(n) + ", Number of dPhi Bins=" + std::to_string(RecoGen->GetNbinsX()/i) + "\n";
					logger << getPtEdges(edges_pT);
					logger << getPhiEdges(edges_dPhi);
					logger << "-> Run getNextDirection( ... )\t";
				}
				
				//Make sure each Hist has a unique title in order to prevent replacing hists due to missing hist names
				RecoGen->SetTitle(("RG" + std::to_string(counter)).c_str());
				ResPt->SetTitle(("ResPt" + std::to_string(counter)).c_str());
				ResPhi->SetTitle(("ResPhi" + std::to_string(counter)).c_str());
				bothSel->SetTitle(("bS" + std::to_string(counter)).c_str());
				genSel->SetTitle(("gS" + std::to_string(counter)).c_str());
				//Calculate the binWidths for current bins and candidates for a merge (always check if it's the last bin)
				float binWidthPt = std::abs(edges_pT.at(n) - edges_pT.at(n-1));
				float binWidthPhi = std::abs(edges_dPhi.at(m) - edges_dPhi.at(m-1));
				float binWidthPt_PT=0;
				float binWidthPhi_PHI=0;
				if (not (n==i)) binWidthPt_PT = std::abs(edges_pT.at(n+1) - edges_pT.at(n-1));
				else binWidthPt_PT = std::abs(edges_pT.at(n) - edges_pT.at(n-2));
				if (m != RecoGen->GetNbinsX()/i) binWidthPhi_PHI = std::abs(edges_dPhi.at(m+1) - edges_dPhi.at(m-1));
				else binWidthPhi_PHI = std::abs(edges_dPhi.at(m) - edges_dPhi.at(m-2));
				
				int nextDir = getNextDirection(RecoGen, ResPt, ResPhi, bothSel, genSel, currentBin, i, p_min, s_min, N_min, eff_min, binWidthPt, binWidthPhi, binWidthPt_PT, binWidthPhi_PHI, rangePt, rangePhi, MoI, logger, is_logging, consider_cut, merge_Overflow, alwaysConsiderMoI, lastBinMergeRule, counter);
				
				std::cout<<"nextDir=" + std::to_string(nextDir)<<std::endl;
				std::cout<<"-----------------------------------------------------------------" <<std::endl;
				//Case: Current bin fulfilles all requirements
				if (nextDir == -1) {
					//Case: It's the very last bin
					if (currentBin == RecoGen->GetNbinsX()) {
						if(is_logging) logger << "Break, all bins are processed.";
						break;
					}
					// Decide if the next bin is the next pT or dPhi bin
					if (nextBinPt) {
						//Case: The checked bin was the last Pt bin; else: It was not the last bin
						if (i == n) {
							//Set the pT index to one and increment the dPhi index
							n=1;
							m++;
						} else {
							//Increment pT index, dPhi index remains unchanged
							n++;
						}
					} else {
						//Case: The checked bin was the last dPhi bin; else: It was not the last bin
						if (m == RecoGen->GetNbinsX()/i) {
							//Set the dPhi index to one and increment the pT index
							n++;
							m=1;
						} else {
							//Increment dPhi index, pT index remains unchanged
							m++;
						}
					}
					//The current bin can be calculated with m, n, i independent of nextBinPt
					currentBin = (m-1)*i + n;
					counter++;
					if(is_logging) logger << "End cycle with properties currentBin=" + std::to_string(currentBin) + ", i=" + std::to_string(i) + ", m=" + std::to_string(m) + ", n=" + std::to_string(n) + ", nextDir= " + std::to_string(nextDir) + "\n\n";
					continue;
				}
				
				//Case: Merging in pT
				if (nextDir == 0 or nextDir == 1) {
					//Merge and calculate the indices that way that the additional step for nextDir=1 is not yet taken into account
					//Check if it's the last pT bin
					if (n != i) {
						RecoGen = (TH2D*)mergePtRecoGen(*RecoGen, i, n).Clone(("RG" + std::to_string(counter)).c_str());
						ResPt = (TH2D*)mergePtRes(*ResPt, i, n).Clone(("ResPt" + std::to_string(counter)).c_str());
						ResPhi = (TH2D*)mergePtRes(*ResPhi, i, n).Clone(("ResPhi" + std::to_string(counter)).c_str());
						bothSel = (TH1D*)mergePt1D(*bothSel, i, n).Clone(("bS" + std::to_string(counter)).c_str());
						genSel = (TH1D*)mergePt1D(*genSel, i, n).Clone(("gS" + std::to_string(counter)).c_str());

						edges_pT.erase(edges_pT.begin() + n);
						//Correct the currentBin by number of merged pT bins in the previous dPhi bins
						currentBin = currentBin - m + 1;
						
					} else {
						RecoGen = (TH2D*)mergePtRecoGen(*RecoGen, i, n-1).Clone(("RG" + std::to_string(counter)).c_str());
						ResPt = (TH2D*)mergePtRes(*ResPt, i, n-1).Clone(("ResPt" + std::to_string(counter)).c_str());
						ResPhi = (TH2D*)mergePtRes(*ResPhi, i, n-1).Clone(("ResPhi" + std::to_string(counter)).c_str());
						bothSel = (TH1D*)mergePt1D(*bothSel, i, n-1).Clone(("bS" + std::to_string(counter)).c_str());
						genSel = (TH1D*)mergePt1D(*genSel, i, n-1).Clone(("gS" + std::to_string(counter)).c_str());
						
						edges_pT.erase(edges_pT.end()-2);
						//Correct the currentBin by number of merged pT bins in the previous and current dPhi bin and decrement the pT bin index
						currentBin = currentBin - m;
						n--;
					}
					//Decrement number of pT bins
					i--;
					//Here, calculate the indices for the case that the current bin is ok. If the current bin is the last one, nothing has to be done here, 
					if (nextDir == 1) {
						//In case that the current bin is now the last bin, the algoritm ends because getNextDirection(...) already confirmed that this bin is ok
						if (currentBin == RecoGen->GetNbinsX()) {
							if(is_logging) logger << "Break, all bins are processed.";
							break;
						} 
						
						// Decide if the next bin is the next pT or dPhi bin
						if (nextBinPt) {
							//Case: The checked bin was the last Pt bin; else: It was not the last bin
							if (i == n) {
								//Set the pT index to one and increment the dPhi index
								n=1;
								m++;
							} else {
								//Increment pT index, dPhi index remains unchanged
								n++;
							}
						} else {
							//Case: The checked bin was the last dPhi bin; else: It was not the last bin
							if (m == RecoGen->GetNbinsX()/i) {
								//Set the dPhi index to one and increment the pT index
								n++;
								m=1;
							} else {
								//Increment dPhi index, pT index remains unchanged
								m++;
							}
						}
						//The current bin can be calculated with m, n, i independent of nextBinPt
						currentBin = (m-1)*i + n;
					}
					counter++;
					if(is_logging) logger << "End cycle with properties currentBin=" + std::to_string(currentBin) + ", i=" + std::to_string(i) + ", m=" + std::to_string(m) + ", n=" + std::to_string(n) + ", nextDir= " + std::to_string(nextDir) + "\n\n";
					continue;
				}
				
				//Case: Merging in dPhi
				if (nextDir == 10 or nextDir == 11) {
					//Merge and calculate the indices that way that the additional step for nextDir=10 is not yet taken into account
					//Check if it's the last dPhi bin
					if (m != RecoGen->GetNbinsX()/i) {
						RecoGen = (TH2D*)mergePhiRecoGen(*RecoGen, i, m).Clone(("RG" + std::to_string(counter)).c_str());
						ResPt = (TH2D*)mergePhiRes(*ResPt, i, m).Clone(("ResPt" + std::to_string(counter)).c_str());
						ResPhi = (TH2D*)mergePhiRes(*ResPhi, i, m).Clone(("ResPhi" + std::to_string(counter)).c_str());
						bothSel = (TH1D*)mergePhi1D(*bothSel, i, m).Clone(("bS" + std::to_string(counter)).c_str());
						genSel = (TH1D*)mergePhi1D(*genSel, i, m).Clone(("gS" + std::to_string(counter)).c_str());

						edges_dPhi.erase(edges_dPhi.begin() + m);
						//No changes in the indices
					} else {
						RecoGen = (TH2D*)mergePhiRecoGen(*RecoGen, i, m-1).Clone(("RG" + std::to_string(counter)).c_str());
						ResPt = (TH2D*)mergePhiRes(*ResPt, i, m-1).Clone(("ResPt" + std::to_string(counter)).c_str());
						ResPhi = (TH2D*)mergePhiRes(*ResPhi, i, m-1).Clone(("ResPhi" + std::to_string(counter)).c_str());
						bothSel = (TH1D*)mergePhi1D(*bothSel, i, m-1).Clone(("bS" + std::to_string(counter)).c_str());
						genSel = (TH1D*)mergePhi1D(*genSel, i, m-1).Clone(("gS" + std::to_string(counter)).c_str());
						
						edges_dPhi.erase(edges_dPhi.end()-2);
						//currentBin is shifted left by the number of pT bins and the current dPhi bin index is decremented.
						currentBin -= i;
						m--;
					}
					//Here, calculate the indices for the case that the current bin is ok. If the current bin is the last one, nothing has to be done here, 
					if (nextDir == 1) {
						//In case that the current bin is now the last bin, the algoritm ends because getNextDirection(...) already confirmed that this bin is ok
						if (currentBin == RecoGen->GetNbinsX()) {
							if(is_logging) logger << "Break, all bins are processed.";
							break;
						} 
						// Decide if the next bin is the next pT or dPhi bin
						if (nextBinPt) {
							//Case: The checked bin was the last Pt bin; else: It was not the last bin
							if (i == n) {
								//Set the pT index to one and increment the dPhi index
								n=1;
								m++;
							} else {
								//Increment pT index, dPhi index remains unchanged
								n++;
							}
						} else {
							//Case: The checked bin was the last dPhi bin; else: It was not the last bin
							if (m == RecoGen->GetNbinsX()/i) {
								//Set the dPhi index to one and increment the pT index
								n++;
								m=1;
							} else {
								//Increment dPhi index, pT index remains unchanged
								m++;
							}
						}
						//The current bin can be calculated with m, n, i independent of nextBinPt
						currentBin = (m-1)*i + n;
					}
					counter++;
					if(is_logging) logger << "End cycle with properties currentBin=" + std::to_string(currentBin) + ", i=" + std::to_string(i) + ", m=" + std::to_string(m) + ", n=" + std::to_string(n) + ", nextDir= " + std::to_string(nextDir) + "\n\n";
					continue;
				}
				if(is_logging) logger << "Error: Loop should not be reaching this point. Properties currentBin=" + std::to_string(currentBin) + ", i=" + std::to_string(i) + ", m=" + std::to_string(m) + ", n=" + std::to_string(n) + ", nextDir= " + std::to_string(nextDir) + "\n\n";
			}
			print_container(edges_pT);
			print_container(edges_dPhi);			
			
			if (not forward) {
				//The bin indices and plots have to be reversed
				std::reverse(edges_pT.begin(),edges_pT.end()); 
				std::reverse(edges_dPhi.begin(),edges_dPhi.end()); 
				
				TH2D* revRecoGen = (TH2D*) RecoGen->Clone("RevRecogen");
				TH2D* revResPt = (TH2D*) ResPt->Clone("RevResPt");
				TH2D* revResPhi = (TH2D*) ResPhi->Clone("RevResPhi");
				TH1D* revGenSel = (TH1D*) genSel->Clone("RevGenSel");
				TH1D* revBothSel = (TH1D*) bothSel->Clone("RevBothSel");
				
				int Nbins = RecoGen->GetNbinsX();
				for (int j=1; j<=Nbins; j++) {
					for (int k=1; k<=Nbins; k++) {
						revRecoGen->SetBinContent(Nbins - j+1 , Nbins - k+1, RecoGen->GetBinContent(j, k));
						revRecoGen->SetBinError(Nbins - j+1 , Nbins - k+1, RecoGen->GetBinError(j, k));
					}
					for (int k=1; k<=ResPt->GetNbinsY(); k++) {
						revResPt->SetBinContent(Nbins - j+1 , k, ResPt->GetBinContent(j, k));
						revResPhi->SetBinContent(Nbins - j+1 , k, ResPhi->GetBinContent(j, k));
						revResPt->SetBinError(Nbins - j+1 , k, ResPt->GetBinContent(j, k));
						revResPhi->SetBinError(Nbins - j+1 , k, ResPhi->GetBinContent(j, k));
					}
					revGenSel->SetBinContent(Nbins - j+1, genSel->GetBinContent(j));
					revBothSel->SetBinContent(Nbins - j+1, bothSel->GetBinContent(j));
					revGenSel->SetBinError(Nbins - j+1, genSel->GetBinError(j));
					revBothSel->SetBinError(Nbins - j+1, bothSel->GetBinError(j));
				}
				RecoGen = revRecoGen;
				ResPt = revResPt;
				ResPhi = revResPhi;
				genSel = revGenSel;
				bothSel = revBothSel;
			}
			
			if(is_logging) {
				logger << "\n\nFinal binnings:\n";
				logger << getPtEdges(edges_pT);
				logger << getPhiEdges(edges_dPhi);
				logger << "\nChecking all bins again.\n";
				
				resultLogger << fileAppendix + "\n";
				resultLogger << getPtEdges(edges_pT);
				resultLogger << getPhiEdges(edges_dPhi);
				resultLogger << "Total Bins =" + std::to_string((edges_pT.size()-1) * (edges_dPhi.size()-1)) + "\n";
			}

			TH2D N_recRecLevel;
			if (DNN) N_recRecLevel = hist::rebinned_double( *(hs2D.getHistogram("DNN_MET/Pt_DeltaPhi", treeName)), edges_pT, edges_dPhi, merge_Overflow, false);
			else N_recRecLevel = hist::rebinned_double( *(hs2D.getHistogram("PuppiMET/Pt_DeltaPhi", treeName)), edges_pT, edges_dPhi, merge_Overflow, false);
			bool checksOK = true;
				
			
			PlotStatistics(&N_recRecLevel, edges_pT, edges_dPhi, false, true, dir.ReplaceAll("Logs/", "Plots/"), fileAppendix + "stat", &saver);
			PlotRMS(&N_recRecLevel, ResPt, edges_pT, edges_dPhi, false, true ,true, dir.ReplaceAll("Logs/", "Plots/"), fileAppendix + "RMSpT", &saver);
			PlotRMS(&N_recRecLevel, ResPhi, edges_pT, edges_dPhi, false, false ,true, dir.ReplaceAll("Logs/", "Plots/"), fileAppendix + "RMSphi", &saver);
			PlotStabPurEff(&N_recRecLevel, RecoGen, bothSel, genSel, edges_pT, edges_dPhi, false, true, dir.ReplaceAll("Logs/", "Plots/"), fileAppendix + "SPE", &saver);
			
			//Check if all bins fulfil the requirements
			float stabMin = calculateStability(RecoGen).GetBinContent(calculateStability(RecoGen).GetMinimumBin());
			if (stabMin >= s_min) {
				if(is_logging) logger << "Stability overall fulfilled. Minimal Stability=" + std::to_string(stabMin) + "\n";
			} else {
				if(is_logging) logger << "Stability NOT fulfilled. Minimal Stability=" + std::to_string(stabMin) + "\n";
				checksOK = false;
			}
			
			float purMin = calculatePurity(RecoGen).GetBinContent(calculatePurity(RecoGen).GetMinimumBin());
			if (purMin >= p_min) {
				if(is_logging) logger << "Purity overall fulfilled. Minimal Purity=" + std::to_string(purMin) + "\n";
			} else {
				if(is_logging) logger << "Purity NOT fulfilled. Minimal Purity=" + std::to_string(purMin) + "\n";
				checksOK = false;
			}
			
			float effMin = calculateEfficiency(bothSel, genSel).GetBinContent(calculateEfficiency(bothSel, genSel).GetMinimumBin());
			if (effMin >= eff_min) {
				if(is_logging) logger << "Efficiency overall fulfilled. Minimal Efficiency=" + std::to_string(effMin) + "\n";
			} else {
				if(is_logging) logger << "Efficiency NOT fulfilled. Minimal Efficiency=" + std::to_string(effMin) + "\n";
				checksOK = false;
			}
			
			float Nmin = RecoGen->ProjectionY("yProjCheck", 1, RecoGen->GetNbinsY())->GetBinContent(RecoGen->ProjectionY("yProjC", 1, RecoGen->GetNbinsY())->GetMinimumBin());
			if (Nmin >= N_min) {
				if(is_logging) logger << "N_rec overall fulfilled. Minimal N_rec=" + std::to_string(Nmin) + "\n";
			} else {
				if(is_logging) logger << "N_rec NOT fulfilled. Minimal N_rec=" + std::to_string(Nmin) + "\n";
				checksOK = false;
			}
			
			bool RMSveto = false;
			for (int j=1; j <= ResPt->GetNbinsX(); j++) {
				//Use new variables for the pT and dPhi bin
				int pTbin = ((j-1) % i) + 1;
				int phibin = floor( (j-1)/i ) + 1;
				float widthPt = edges_pT.at(pTbin) - edges_pT.at(pTbin -1);
				float widthPhi = edges_dPhi.at(phibin) - edges_dPhi.at(phibin -1);
				if (widthPt < ResPt->GetBinContent(j)) {
					if(is_logging) logger << "RMS pT too great at bin j=" + std::to_string(j) + "\n";
					RMSveto = true;
				}
				if (widthPhi < ResPhi->GetBinContent(j)) {
					if(is_logging) logger << "RMS dPhi too great at bin j=" + std::to_string(j) + "\n";
					RMSveto = true;
				}
			}
			if (RMSveto == false) {
				if(is_logging) logger << "RMS ok in all bins.\n";
			} else {
				if(is_logging) logger << "RMS NOT ok in all bins.\n";
				checksOK = false;
			}
			
			if (checksOK) {
				if(is_logging) resultLogger << "Final Checks: OK\n";
			} else {
				if(is_logging) resultLogger << "Final Checks: NOT OK\n";
			}
			
			
			//Ploting BSM, SBR, Ratio
			int Nbins = RecoGen->GetNbinsX();
			
			std::map<TString, TString> le { {"DrellYan_NLO", "DY"}, {"SingleTop", "Single t"}, {"TTbar_diLepton", "t#bar{t} ll"}, {"TTbar_diLepton_tau", "t#bar{t} ll#tau" }, {"TTbar_hadronic", "t#bar{t} hadronic"}, {"TTbar_singleLepton", "t#bar{t} single l"}, {"ttG", "ttG"}, {"ttW", "ttW"}, {"ttZ", "ttZ"}, {"WJetsToLNu", "W+jets"}, {"WW", "WW"}, {"WZ", "WZ"}, {"ZZ", "ZZ"}, {"ttX", "ttX"}, {"Diboson", "Diboson"}, {"tt other", "t#bar{t} other"}, {"T1tttt_1200_800", "T1tttt_1200_800"}, {"T2tt_650_350", "T2tt_650_350"}, {"DM_scalar_1_200", "DM_scalar_1_200"}};
			std::vector<TString> SMtrees = {"TTbar_hadronic","ZZ", "WZ", "WW", "WJetsToLNu", "ttW", "ttZ", "ttG", "TTbar_singleLepton", "SingleTop", "DrellYan_NLO", "TTbar_diLepton_tau", "TTbar_diLepton"};
			std::vector<TString> BSMtrees = {"T1tttt_1200_800", "T2tt_650_350", "DM_scalar_1_200"};
			std::vector<TString> BSMtrees2 = {"T1tttt_1500_100"};
			std::vector<TString> alltrees;
			alltrees.reserve( SMtrees.size() + BSMtrees.size() + BSMtrees2.size()); 
			alltrees.insert( alltrees.end(), SMtrees.begin(), SMtrees.end() );
			alltrees.insert( alltrees.end(), BSMtrees.begin(), BSMtrees.end() );
			alltrees.insert( alltrees.end(), BSMtrees2.begin(), BSMtrees2.end() );
			std::map<const TString,Color_t> colormap { {"DrellYan_NLO", kOrange}, {"SingleTop", kMagenta+3}, {"TTbar_diLepton", kRed-6}, {"TTbar_diLepton_tau", kRed-4 }, {"TTbar_hadronic", kRed}, {"TTbar_singleLepton", kRed-2}, {"ttG", kGreen}, {"ttW", kGreen+3}, {"ttZ", kGreen-9}, {"WJetsToLNu", kGray}, {"WW", kCyan}, {"WZ", kCyan+3}, {"ZZ", kCyan+4}, {"ttX", kGreen}, {"Diboson", kCyan}, {"tt other", kRed}, {"T1tttt_1200_800", kBlue+2}, {"T2tt_650_350", kPink-3}, {"DM_scalar_1_200", kYellow+2}};
			
			hist::Histograms<TH1D> hs(alltrees);
			
			hs.addHist("newBinning/N_rec", ";;Events/Bin", Nbins, 0, Nbins);
			
			//Fill 1D hists
			for (TString sSample : alltrees){
				hs.setCurrentSample(sSample);
				TFile file2(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/100.0/"+sSample+".root",cfg.year.Data(),cfg.treeVersion.Data()),"read");    //open input file
				TTreeReader reader2("ttbar_res100.0/"+sSample, &file2);      //open tree in root file
			   
				//define variables taken from tree
				TTreeReaderValue<float> PtNuNu2   (reader2, "PtNuNu");
				TTreeReaderValue<float> PuppiMET2   (reader2, "PuppiMET");
				TTreeReaderValue<float> DNN_MET2   (reader2, "DNN_MET_pT");
				TTreeReaderValue<float> Phi_Puppi2  (reader2, "Phi_recPuppi");
				TTreeReaderValue<float> Phi_NuNu2   (reader2, "Phi_NuNu");
				TTreeReaderValue<float> Phi_DNN2   (reader2, "DNN_MET_dPhi_nextLep");
				TTreeReaderValue<float> MC_weight2   (reader2, "N");
				TTreeReaderValue<float> SF_weight2   (reader2, "SF");
				//~ std::cin.ignore();
				//loop over events in tree
				processEvents=cfg.processFraction*reader2.GetEntries(true);
				iEv=0;
				io::log*"Current Tree: ";
				io::log* sSample;
				int row = 0;
				bool rec_selection=false;
				while (reader2.Next()){
					iEv++;
					if (iEv>processEvents) break;
					if (iEv%(std::max(processEvents/10,1))==0){      //logging stuff
						io::log*".";
						io::log.flush();
					}
				  
					//set weight for current event
					hs.setFillWeight((*MC_weight2)*(*SF_weight2));
					 
					// set booleans for selection
					rec_selection=false;
					gen_selection=false;
					if(*PtNuNu2>0 and (merge_Overflow or *PtNuNu2 <= 400)) gen_selection=true;
					if(*PuppiMET2>0 and (not DNN and (merge_Overflow or *PuppiMET2 <= 400))) rec_selection=true;
					if(*DNN_MET2>0 and (DNN and (merge_Overflow or *DNN_MET2 <= 400))) rec_selection=true;
					if (rec_selection) {
						if (DNN) {
							row = getIndex(*Phi_DNN2, *DNN_MET2, edges_dPhi, edges_pT, true);
							hs.fill("newBinning/N_rec", row-0.5);
						} else {
							row = getIndex(*Phi_Puppi2, *PuppiMET2, edges_dPhi, edges_pT, true);
							hs.fill("newBinning/N_rec", row-0.5);
						}
					}	
				}
				if (merge_Overflow) hs.mergeOverflow();
				io::log<<"";	
				file2.Close();
			}
			
			
			std::vector<TString> cs = {"WJetsToLNu", "Diboson", "ttX", "DrellYan_NLO", "SingleTop", "tt other", "TTbar_diLepton"};
			std::vector<TString> BGtreesSM = {"tt other", "SingleTop", "DrellYan_NLO", "ttX", "Diboson", "WJetsToLNu"};
			hs.combineSamples("Diboson",{"WW","WZ","ZZ"});
			hs.combineSamples("ttX",{"ttW","ttZ", "ttG"});
			hs.combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
			hs.combineSamples("SMtrees", SMtrees);
			
			for(TString sSample : BGtreesSM) {
				TH1D* tempHist=hs.getHistogram("newBinning/N_rec", sSample);
				tempHist->SetStats(false);
				tempHist->SetLineColor(colormap.find(sSample)->second);
				tempHist->SetFillColor(colormap.find(sSample)->second);
				tempHist->Draw("hist e");
				saver.save(can, "N_rec/" + sSample);
			}
			
			//Signal over sqrt(Background), Signal BSM variables, Background SM
			can.SetLogy(false);
			TH1D* S_T1 = hs.getHistogram("newBinning/N_rec", "T1tttt_1200_800");
			TH1D* S2_T1 = hs.getHistogram("newBinning/N_rec", "T1tttt_1500_100");
			TH1D* S_T2 = hs.getHistogram("newBinning/N_rec", "T2tt_650_350");
			TH1D* S_DM = hs.getHistogram("newBinning/N_rec", "DM_scalar_1_200");
			TH1D  Background = getSqrt(*(hs.getHistogram("newBinning/N_rec", "SMtrees")));
			S_T1->Divide(&Background);
			S2_T1->Divide(&Background);
			S_T2->Divide(&Background);
			S_DM->Divide(&Background);
			TH1D SBR_T1 = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
			TH1D SBR_T2 = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
			TH1D SBR_DM = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
			for (int j = 1; j <= Nbins; j++) {
				SBR_T1.SetBinContent(j, S_T1->GetBinContent(j));
				SBR_T1.SetBinError(j, S_T1->GetBinError(j));
				SBR_T2.SetBinContent(j, S_T2->GetBinContent(j));
				SBR_T2.SetBinError(j, S_T2->GetBinError(j));
				SBR_DM.SetBinContent(j, S_DM->GetBinContent(j));
				SBR_DM.SetBinError(j, S_DM->GetBinError(j));
			}
			resultLogger << "Max. BSM SBR: T1tttt_1200_800: " +  TString::Format("%.2f", SBR_T1.GetMaximum()) + "\t T2tt_650_350: " +  TString::Format("%.2f", SBR_T2.GetMaximum()) + "\t DM_scalar_1_200: " +  TString::Format("%.2f", SBR_DM.GetMaximum()) + "\t T1tttt_1500_100: " +  TString::Format("%.2f", S2_T1->GetMaximum())+ "\n\n";
			SBR_T1.SetName("SBR_T1");
			SBR_T2.SetName("SBR_T2");
			SBR_DM.SetName("SBR_DM");
			SBR_T1.SetStats(false);
			SBR_T1.GetXaxis()->SetTitle("%MET");
			SBR_T1.GetYaxis()->SetTitle("S/#sqrt{B}");
			SBR_T1.SetLineColor(colormap.find("T1tttt_1200_800")->second);
			SBR_T1.SetFillColor(colormap.find("T1tttt_1200_800")->second);
			SBR_T1.SetMarkerColor(colormap.find("T1tttt_1200_800")->second);
			SBR_T2.SetLineColor(colormap.find("T2tt_650_350")->second);
			SBR_T2.SetFillColor(colormap.find("T2tt_650_350")->second);
			SBR_T2.SetMarkerColor(colormap.find("T2tt_650_350")->second);
			SBR_DM.SetLineColor(colormap.find("DM_scalar_1_200")->second);
			SBR_DM.SetFillColor(colormap.find("DM_scalar_1_200")->second);
			SBR_DM.SetMarkerColor(colormap.find("DM_scalar_1_200")->second);
			SBR_T1.SetMarkerStyle(20);
			SBR_T2.SetMarkerStyle(21);
			SBR_DM.SetMarkerStyle(22);
			SBR_T1.SetMarkerSize(0.8);
			SBR_T2.SetMarkerSize(0.8);
			SBR_DM.SetMarkerSize(0.8);
			float upperEnd=0.81;
			float lowerEnd=0.002;
			SBR_T1.GetYaxis()->SetRangeUser(lowerEnd, upperEnd);
			SBR_T1.Draw("pe1");
			SBR_T2.Draw("same pe1");
			SBR_DM.Draw("same pe1");
			
			int phiBins = edges_dPhi.size() - 1;
			SBR_T1.GetXaxis()->SetNdivisions(8*phiBins, false);
			SBR_T1.GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
			for (int i=0; i<=phiBins-1; i++){
				SBR_T1.GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
				SBR_T1.GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
			}
			TLine * aline = new TLine();
			TLatex * atext = new TLatex();
			atext->SetTextSize(0.03);
			aline->SetLineWidth(2);
			aline->SetLineColor(kBlack);

			for (int i=1; i<=phiBins-1; i++){
				aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
			}
			for (int i=0;i < phiBins; i++){
				if (i==0) atext->DrawLatex(100, 0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("|#Delta#phi|<%.1f",edges_dPhi.at(1))); 
				else {
					if (i==phiBins-1) atext->DrawLatex(400*i+100, 0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|",edges_dPhi.at(i)));
					else atext->DrawLatex(400*i+100, 0.9*(upperEnd-lowerEnd)+lowerEnd,TString::Format("%.1f<|#Delta#phi|<%.1f",edges_dPhi.at(i), edges_dPhi.at(i+1)));
				}
			}
			
			auto legend_SBR = new TLegend(gPad->GetLeftMargin(),0.91,1-gPad->GetRightMargin(),1-gPad->GetTopMargin());
			legend_SBR->AddEntry("SBR_T1", "T1tttt_1200_800","lep");
			legend_SBR->AddEntry("SBR_T2", "T2tt_650_350","lep");
			legend_SBR->AddEntry("SBR_DM", "DM_scalar_1_200","lep");
			legend_SBR->AddEntry((TObject*)0, " ","");
			TLegendEntry *leSBR;
			if (DNN) leSBR = legend_SBR->AddEntry((TObject*)0, "DNN", "");
			else leSBR = legend_SBR->AddEntry((TObject*)0, "PUPPI", "");
			leSBR->SetTextFont(62);
			legend_SBR->SetShadowColor(0);
			legend_SBR->SetBorderSize(1);
			legend_SBR->SetFillColor(kWhite);
			legend_SBR->SetNColumns(5);
			legend_SBR->Draw();

			saver.save(can, fileAppendix + "SBR");
			can.SaveAs(dir.ReplaceAll("Logs/", "Plots/") + fileAppendix + "SBR.pdf");
			
			
			//Proportion of Events in the total of Events (SM)
			can.SetLogy(true);
			std::vector<int> marker {20, 21, 22, 23, 29, 33, 47};
			auto proplegend = new TLegend(gPad->GetLeftMargin(),0.88,1-gPad->GetRightMargin(),1-gPad->GetTopMargin());
			proplegend->SetShadowColor(0);
			proplegend->SetFillColor(kWhite);
			proplegend->SetBorderSize(1);
			proplegend->SetNColumns(8);
			
			TH1D* SignalTemp = hs.getHistogram("newBinning/N_rec", "TTbar_diLepton");
			SignalTemp->Divide(hs.getHistogram("newBinning/N_rec", "SMtrees"));
			TH1D Signal = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
			for (int j = 1; j <= Nbins; j++) {
				Signal.SetBinContent(j, SignalTemp->GetBinContent(j));
				Signal.SetBinError(j, SignalTemp->GetBinError(j));
			}
			Signal.Scale(100);
			Signal.SetName("Signal");
			Signal.SetStats(false);
			Signal.GetXaxis()->SetTitle("%MET");
			Signal.GetYaxis()->SetTitle("%");
			Signal.GetYaxis()->SetRangeUser(0.03, 310);
			Signal.SetLineColor(colormap.find("TTbar_diLepton")->second);
			Signal.SetFillColor(colormap.find("TTbar_diLepton")->second);
			Signal.SetMarkerStyle(marker.at(0));
			marker.erase(marker.begin());
			Signal.SetMarkerSize(0.8);
			Signal.SetMarkerColor(colormap.find("TTbar_diLepton")->second);
			Signal.Draw("pe1");
			TString legendText = le.find("TTbar_diLepton")->second + "";
			proplegend->AddEntry("Signal", legendText, "lep");
			for (TString sSample : BGtreesSM){
				TH1D* tempHistBG = hs.getHistogram("newBinning/N_rec", sSample);
				tempHistBG->Divide(hs.getHistogram("newBinning/N_rec", "SMtrees"));
				TH1D BG = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
				for (int j = 1; j <= Nbins; j++) {
					BG.SetBinContent(j, tempHistBG->GetBinContent(j));
					BG.SetBinError(j, tempHistBG->GetBinError(j));
				}
				BG.Scale(100);
				std::string name(sSample);
				name += "Ratio";
				BG.SetName(name.c_str());
				BG.SetStats(false);
				BG.SetLineColor(colormap.find(sSample)->second);
				BG.SetFillColor(colormap.find(sSample)->second);
				BG.SetMarkerStyle(marker.at(0));
				marker.erase(marker.begin());
				BG.SetMarkerSize(0.8);
				BG.SetMarkerColor(colormap.find(sSample)->second);
				BG.DrawCopy("same pe1");
				name+="_copy";
				legendText = le.find(sSample)->second + " ";
				proplegend->AddEntry(name.c_str(), legendText, "lep");
			}
			
			
			Signal.GetXaxis()->SetNdivisions(8*phiBins, false);
			Signal.GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
			for (int i=0; i<=phiBins-1; i++){
				Signal.GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
				Signal.GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
				Signal.GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
				Signal.GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
				Signal.GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
				Signal.GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
				Signal.GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
				Signal.GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
			}
			TLine * alineR = new TLine();
			TLatex * atextR = new TLatex();
			atextR->SetTextSize(0.03);
			alineR->SetLineWidth(2);
			alineR->SetLineColor(kBlack);
			float upperEndR=180;
			float lowerEndR=0.03;

			for (int i=1; i<=phiBins-1; i++){
				alineR->DrawLine(i*400,lowerEndR,i*400,upperEndR);
			}
			for (int i=0;i < phiBins; i++){
				if (i==0) atextR->DrawLatex(100, 0.56*upperEndR,TString::Format("|#Delta#phi|<%.1f",edges_dPhi.at(1))); 
				else {
					if (i==phiBins-1) atextR->DrawLatex(400*i+100, 0.56*upperEndR,TString::Format("%.1f<|#Delta#phi|",edges_dPhi.at(i)));
					else atextR->DrawLatex(400*i+100, 0.56*upperEndR,TString::Format("%.1f<|#Delta#phi|<%.1f",edges_dPhi.at(i), edges_dPhi.at(i+1)));
				}
			}
			
			TLegendEntry *leProp;
			if (DNN) leProp = proplegend->AddEntry((TObject*)0, "DNN", "");
			else leProp = proplegend->AddEntry((TObject*)0, "MET", "");
			leProp->SetTextFont(62);
			proplegend->Draw();
			
			saver.save(can, fileAppendix + "Ratio");
			can.SaveAs(dir.ReplaceAll("Logs/", "Plots/") + fileAppendix + "Ratio.pdf");
			
			//Plot Stacked Plots 
			TString fS("WJetsToLNu");
			THStack *Stack = new THStack("stack", "");
			auto StackLegend = new TLegend(gPad->GetLeftMargin(), 0.86,1-gPad->GetRightMargin(),1-gPad->GetTopMargin());
			StackLegend->SetShadowColor(0);
			StackLegend->SetFillColor(kWhite);
			StackLegend->SetBorderSize(1);
			StackLegend->SetNColumns(6);
			TH1D* fSTemp = (TH1D*) hs.getHistogram("newBinning/N_rec", fS)->Clone();
			TH1D firstStack = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
			for (int j = 1; j <= Nbins; j++) {
				firstStack.SetBinContent(j, fSTemp->GetBinContent(j));
				firstStack.SetBinError(j, fSTemp->GetBinError(j));
			}
			firstStack.SetFillColor(colormap.find(fS)->second);
			firstStack.SetFillStyle(1001);
			firstStack.SetLineColor(kBlack);
			firstStack.SetLineWidth(1);
			Stack->Add(&firstStack, "hist");
			//~ StackLegend->AddEntry(&firstStack, le.find(fS)->second, "f");
			int leCounter = 1;
			
			std::vector<TH1D> histVector;
			
			for (TString sSample : cs){
				if (sSample == fS) continue;
				TH1D* tempHistStack = hs.getHistogram("newBinning/N_rec", sSample);
				TH1D st = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
				for (int j = 1; j <= Nbins; j++) {
					st.SetBinContent(j, tempHistStack->GetBinContent(j));
					st.SetBinError(j, tempHistStack->GetBinError(j));
				}
				std::string name(sSample);
				name += "Stack";
				st.SetName(name.c_str());
				st.SetStats(false);
				st.SetFillColor(colormap.find(sSample)->second);
				st.SetFillStyle(1001);
				st.SetLineColor(kBlack);
				st.SetLineWidth(1);
				histVector.push_back(st);
			}
			int p = 0;
			for (TString sSample : cs){
				if (sSample == fS) continue;
				vector <TH1D>::iterator it = histVector.begin() + p;
				TH1D * HistToAdd = &(*it);
				p++;
				Stack->Add(HistToAdd, "hist");
				leCounter++;
			}
			
			Stack->SetMinimum(1);
			Stack->SetMaximum(Stack->GetMaximum()*3.6);
			Stack->Draw();
			can.Update();
			Stack->GetYaxis()->SetTitle("Events/Bin");
			Stack->GetXaxis()->SetTitle("%MET");
			
			for (TString sSample : BSMtrees){
				TH1D* tempHistBSM=hs.getHistogram("newBinning/N_rec", sSample);
				TH1D BSM = HistTrafo_2D(&N_recRecLevel, edges_pT, edges_dPhi);
				for (int j = 1; j <= Nbins; j++) {
					BSM.SetBinContent(j, tempHistBSM->GetBinContent(j));
					BSM.SetBinError(j, tempHistBSM->GetBinError(j));
				}
				std::string name(sSample);
				name += "BSMstack";
				BSM.SetName(name.c_str());
				BSM.SetLineColor(colormap.find(sSample)->second);
				BSM.SetMarkerSize(0);
				BSM.SetStats(false);
				BSM.DrawCopy("same hist");
				name += "_copy";
				//~ StackLegend->AddEntry(name.c_str(), le.find(sSample)->second, "l");
			}
			//####
			TLegendEntry * leStack;
			std::vector<TString> csNewOrder = {"TTbar_diLepton", "SingleTop", "ttX", "WJetsToLNu", "T2tt_650_350", "tt other", "DrellYan_NLO", "Diboson", "T1tttt_1200_800", "DM_scalar_1_200"};
			int r=0;
			for (TString sSample : csNewOrder){
				if (r == 5) {
					leStack = StackLegend->AddEntry((TObject*)0, "DNN  ", "");		
					leStack->SetTextFont(62);
				}
				TH1D* tempHist=hs.getHistogram("newBinning/N_rec", sSample);
				tempHist->SetLineColor(colormap.find(sSample)->second);
				if (not sSample.Contains("0")){
					tempHist->SetFillColor(colormap.find(sSample)->second);
					tempHist->SetFillStyle(1001);
					tempHist->SetLineColor(kBlack);
					tempHist->SetLineWidth(1);
					StackLegend->AddEntry(tempHist, le.find(sSample)->second, "f");
				} else {
					StackLegend->AddEntry(tempHist, le.find(sSample)->second, "l");
				}
				r++;
			}
			//####
			
			Stack->GetXaxis()->SetNdivisions(8*phiBins, false);
			Stack->GetXaxis()->ChangeLabel(8*phiBins+1,-1,-1,-1,-1,-1," ");
			for (int i=0; i<=phiBins-1; i++){
				Stack->GetXaxis()->ChangeLabel(i*8+1,-1,-1,-1,-1,-1," ");
				Stack->GetXaxis()->ChangeLabel(i*8+2,-1,-1,-1,-1,-1," ");
				Stack->GetXaxis()->ChangeLabel(i*8+3,-1,-1,-1,-1,-1,"100");
				Stack->GetXaxis()->ChangeLabel(i*8+4,-1,-1,-1,-1,-1," ");
				Stack->GetXaxis()->ChangeLabel(i*8+5,-1,-1,-1,-1,-1," ");
				Stack->GetXaxis()->ChangeLabel(i*8+6,-1,-1,-1,-1,-1," ");
				Stack->GetXaxis()->ChangeLabel(i*8+7,-1,-1,-1,-1,-1,"300");
				Stack->GetXaxis()->ChangeLabel(i*8+8,-1,-1,-1,-1,-1," ");
			}
			TLine * alineS = new TLine();
			TLatex * atextS = new TLatex();
			atextS->SetTextSize(0.03);
			alineS->SetLineWidth(2);
			alineS->SetLineColor(kBlack);
			//GetY returns the values in the form of 10^x where x is the return
			double_t upperEndS=can.GetFrame()->GetY2();
			double_t lowerEndS=can.GetFrame()->GetY1();
			double labelPos = pow(10, lowerEndS + 0.84*(upperEndS-lowerEndS));
			upperEndS = pow(10, upperEndS);
			lowerEndS = pow(10, lowerEndS);

			for (int i=1; i<=phiBins-1; i++){
				alineS->DrawLine(i*400,lowerEndS,i*400,upperEndS);
			}
			for (int i=0;i < phiBins; i++){
				if (i==0) atextS->DrawLatex(100, labelPos,TString::Format("|#Delta#phi|<%.1f",edges_dPhi.at(1))); 
				else {
					if (i==phiBins-1) atextS->DrawLatex(400*i+100, labelPos,TString::Format("%.1f<|#Delta#phi|",edges_dPhi.at(i)));
					else atextR->DrawLatex(400*i+100, labelPos,TString::Format("%.1f<|#Delta#phi|<%.1f",edges_dPhi.at(i), edges_dPhi.at(i+1)));
				}
			}
			
			StackLegend->Draw();
			saver.save(can, fileAppendix + "Stack");
			can.SaveAs(dir.ReplaceAll("Logs/", "Plots/") + fileAppendix + "Stack.pdf");
			
			//End
			dir.ReplaceAll("Plots/", "Logs/");
			if(is_logging) logger << "END";
			logger.close();
		}
		
		resultLogger.close();
	}
	return 0;
}
