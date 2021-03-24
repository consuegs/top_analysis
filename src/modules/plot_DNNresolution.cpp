//Script to plot DNN resolution studies

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_DNNresolution");
   
   TCanvas can;
   can.cd();
   // ~io::RootFileReader histReader(TString::Format("binningUnfolding_diLepton%.1f.root",cfg.processFraction*100));
   io::RootFileReader histReader(TString::Format("binningUnfolding_T2tt_650_350%.1f.root",cfg.processFraction*100));
   
   std::vector<float> axisLimits_low={0,0,0,50,50,50};
   std::vector<float> axisLimits_high={150,150,200,250,300,500};
   
   std::vector<TString> metRegions={"p_{T}^{miss}<40GeV","40GeV<p_{T}^{miss}<80GeV","80GeV<p_{T}^{miss}<120GeV","120GeV<p_{T}^{miss}<160GeV","160GeV<p_{T}^{miss}<230GeV","230GeV<p_{T}^{miss}"};
   
   int i=0;
   for (TString bin :{"bin1","bin2","bin3","bin4","bin5","bin6"}){
      
      //Distribution
      TH1F* genMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/GenMET_"+bin);
      TH1F* PuppiMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMET_"+bin);
      TH1F* PuppiMETcorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMETcorr_"+bin);
      TH1F* PFMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PFMET_"+bin);
      
      genMET_hist->Rebin(2);
      PuppiMET_hist->Rebin(2);
      PuppiMETcorr_hist->Rebin(2);
      PFMET_hist->Rebin(2);
      // ~genMET_hist->Rebin(20);
      // ~PuppiMET_hist->Rebin(20);
      // ~PuppiMETcorr_hist->Rebin(20);
      // ~PFMET_hist->Rebin(20);
      
      genMET_hist->SetLineColor(kGreen);
      PuppiMET_hist->SetLineColor(kRed);
      PuppiMETcorr_hist->SetLineColor(kBlue);
      PFMET_hist->SetLineColor(kBlack);
      
      genMET_hist->SetStats(0);
      genMET_hist->GetXaxis()->SetTitle("MET (GeV)");
      genMET_hist->GetYaxis()->SetTitle("Events/Bin");
      genMET_hist->GetYaxis()->SetRangeUser(0,1.5*PuppiMET_hist->GetMaximum());
      genMET_hist->GetXaxis()->SetRangeUser(axisLimits_low[i],axisLimits_high[i]);
      
      genMET_hist->Draw("hist");
      PuppiMET_hist->Draw("hist same");
      PuppiMETcorr_hist->Draw("hist same");
      PFMET_hist->Draw("hist same");
      
      gfx::LegendEntries le;
      le.append(*PuppiMET_hist,"PuppiMET","l");
      le.append(*PFMET_hist,"PFMET","l");
      le.append(*PuppiMETcorr_hist,"DNN","l");
      le.append(*genMET_hist,"genMET","l");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel(metRegions[i],1);
      label.Draw();
      saver.save(can,"genMETtarget/dist_"+bin,true,true);
      
      //Resolution
      TH1F* diff_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diff_"+bin);
      TH1F* diffPF_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffPF_"+bin);
      TH1F* diffCorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffcorr_"+bin);
      TH1F* diffScaled_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffscaled_"+bin);
      
      // ~diff_hist->Rebin(2);
      // ~diffPF_hist->Rebin(2);
      // ~diffCorr_hist->Rebin(2);
      // ~diffScaled_hist->Rebin(2);
      diff_hist->Rebin(20);
      diffPF_hist->Rebin(20);
      diffCorr_hist->Rebin(20);
      diffScaled_hist->Rebin(20);
      
      diff_hist->SetLineColor(kRed);
      diffPF_hist->SetLineColor(kBlack);
      diffCorr_hist->SetLineColor(kBlue);
      diffScaled_hist->SetLineColor(kGray);
      
      diff_hist->SetStats(0);
      diff_hist->GetYaxis()->SetRangeUser(0,1.4*diffCorr_hist->GetMaximum());
      diff_hist->GetXaxis()->SetTitle("genMET-recoMET (GeV)");
      diff_hist->GetYaxis()->SetTitle("Events/Bin");
      
      diff_hist->Draw("hist");
      diffPF_hist->Draw("hist same");
      diffCorr_hist->Draw("hist same");
      // ~diffScaled_hist->Draw("hist same");
      
      gfx::LegendEntries le2;
      le2.append(*diff_hist,TString::Format("PuppiMET(#mu=%.1f #sigma=%.1f)",diff_hist->GetMean(),diff_hist->GetRMS()),"l");
      le2.append(*diffPF_hist,TString::Format("PFMET(#mu=%.1f #sigma=%.1f)",diffPF_hist->GetMean(),diffPF_hist->GetRMS()),"l");
      le2.append(*diffCorr_hist,TString::Format("DNN(#mu=%.1f #sigma=%.1f)",diffCorr_hist->GetMean(),diffCorr_hist->GetRMS()),"l");
      // ~le2.append(*diffScaled_hist,TString::Format("genMET-scaled PuppiMET(#mu=%.1f #sigma=%.1f)",diffScaled_hist->GetMean(),diffScaled_hist->GetRMS()),"l");
      TLegend leg2=le2.buildLegend(.45,.75,1-1.5*gPad->GetRightMargin(),-1,1);
      leg2.Draw();
      
      TLatex label2=gfx::cornerLabel(metRegions[i],1);
      label2.Draw();
      saver.save(can,"genMETtarget/diff_"+bin,true,true);
      
      i++;
   }
}
