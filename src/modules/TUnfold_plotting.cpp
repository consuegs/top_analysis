//Script to plot the result of TUnfold_unfolding
#include <iostream>
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


#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;

Config const &cfg=Config::get();

extern "C"
void run()
{
   // unfolded sample
   TString sample="MadGraph";
   // ~TString sample="dilepton";
   
   // response sample
   // ~TString sample_response="MadGraph";
   TString sample_response="dilepton";
   
   //==============================================
   // step 1 : open output file
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("TUnfold_plotting%.1f",cfg.processFraction*100));

   //==============================================
   // step 2 : read binning schemes and input histograms
   io::RootFileReader histReader(TString::Format("TUnfold%.1f.root",cfg.processFraction*100));
   TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
   TString input_loc_result="TUnfold_results_"+sample+"_"+sample_response;

   TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
   TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

   if((!detectorBinning)||(!generatorBinning)) {
      cout<<"problem to read binning schemes\n";
   }

   // read histograms
   TH1F *realDis=histReader.read<TH1F>(input_loc+"/histDataTruth");
   TH1F *unfolded=histReader.read<TH1F>(input_loc_result+"/hist_unfoldedResult");

   if((!realDis)||(!unfolded)) {
      cout<<"problem to read input histograms\n";
   }

   //========================
   // Step 3: plotting
   
   // ~TCanvas can("c1","c1", 1200, 600);
   
   gfx::SplitCan can;
   can.can_.SetWindowSize(1800,600);
   can.pU_.SetLogy();
   can.pU_.cd();
   
   // ~can.SetLogy();

   
   Double_t xbins[14] = {-30,0,60,120,230,400,460,520,630,800,860,920,1030,1200};
   unfolded->SetBins(13,xbins);
   realDis->SetBins(13,xbins);
   // ~realDis_signal->SetBins(12,xbins);
   unfolded->GetXaxis()->SetTickLength(0.);
   unfolded->GetYaxis()->SetTickLength(0.008);
   unfolded->GetXaxis()->SetTitleOffset(1.5);
   unfolded->GetXaxis()->CenterLabels(false);
   unfolded->GetXaxis()->SetBinLabel(1,"fakes");
   unfolded->GetXaxis()->SetBinLabel(2,"0-60");
   unfolded->GetXaxis()->SetBinLabel(3,"60-120");
   unfolded->GetXaxis()->SetBinLabel(4,"120-230");
   unfolded->GetXaxis()->SetBinLabel(5,">230");
   unfolded->GetXaxis()->SetBinLabel(6,"0-60");
   unfolded->GetXaxis()->SetBinLabel(7,"60-120");
   unfolded->GetXaxis()->SetBinLabel(8,"120-230");
   unfolded->GetXaxis()->SetBinLabel(9,">230");
   unfolded->GetXaxis()->SetBinLabel(10,"0-60");
   unfolded->GetXaxis()->SetBinLabel(11,"60-120");
   unfolded->GetXaxis()->SetBinLabel(12,"120-230");
   unfolded->GetXaxis()->SetBinLabel(13,">230");
   unfolded->LabelsOption("v");
   unfolded->SetMaximum(4*unfolded->GetMaximum());
   unfolded->SetMinimum(1);
   unfolded->SetLineColor(kBlack);
   unfolded->SetTitle(";p_{T}^{#nu#nu} (GeV);arbitrary unit");
   unfolded->SetStats(false);
   
   // ~unfolded->Scale(1.0/(10e5));
   // ~realDis->Scale(1.0/(10e5));
   // ~realDis_signal->Scale(1.0/(10e5));
   
   unfolded->Draw("pe1");
   realDis->SetLineColor(kRed-6);
   realDis->SetFillColor(kRed-6);
   realDis->Draw("hist same");
   // ~realDis_signal->SetLineColor(kBlue);
   // ~realDis_signal->Draw("hist same");
   
   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.03);
   //~ atext->SetTextFont(42);
   aline->SetLineWidth(2);
   aline->DrawLine(0,1,0,unfolded->GetMaximum());
   aline->DrawLine(800,1,800,unfolded->GetMaximum());
   aline->DrawLine(400,1,400,unfolded->GetMaximum());
   aline->DrawLine(800,1,800,unfolded->GetMaximum());
   aline->SetLineStyle(2);
   atext->DrawLatex(50,0.5*unfolded->GetMaximum(),"0<|#Delta#phi|(genMET,nearest l)<0.7");
   atext->DrawLatex(450,0.5*unfolded->GetMaximum(),"0.7<|#Delta#phi|(genMET,nearest l)<1.4");
   atext->DrawLatex(850,0.5*unfolded->GetMaximum(),"1.4<|#Delta#phi|(genMET,nearest l)<3.14");
   
   
   gfx::LegendEntries legE;
   legE.append(*unfolded,"Unfolded","pe");
   legE.append(*realDis,"MC true ttbar","l");
   // ~legE.append(*realDis_signal,"MC true signal","l");
   TLegend leg=legE.buildLegend(.2,.25,0.35,.4,1);
   leg.SetTextSize(0.035);
   leg.Draw();
   
   // ~can.pU_.RedrawAxis();
   unfolded->Draw("axis same");
   
   can.pL_.cd();
   TH1F ratio=hist::getRatio(*unfolded,*realDis,"ratio",hist::ONLY1);
   ratio.Draw();
   
   //===========================
   // Step 4: save plot
   saver.save(can,sample+"_"+sample_response);

}
