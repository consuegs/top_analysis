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
   TCanvas canCombined;
   gfx::LegendEntries leCombined;
   can.cd();
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_addVariables%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_normDistr%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_normDistr_recoEvents%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_InclusiveDNN_Patched_LR_EP%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_genMetTarget_normDistr%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_test_cppflow_new%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_normDistr_emu_recoEvents%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_normDistr_100EP_%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_pTnunuTarget_normDistr%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_genMETweighted%.1f.root",cfg.processFraction*100);
   
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_test_cppflow_new_2D%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_xyTarget_2D_%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_xyTarget_2D_bothCorr_%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_xyTarget_30EP_2D_bothCorr_%.1f.root",cfg.processFraction*100);
   TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_xyTarget_30EP_2D_bothCorr_noOverflow_%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_xyTarget_JetLepXY_50EP_2D_bothCorr_noOverflow_%.1f.root",cfg.processFraction*100);
   
   io::RootFileReader histReader(inputFile);
   
   std::vector<float> axisLimits_low={0,0,0,50,50,50};
   std::vector<float> axisLimits_high={150,150,200,250,300,500};
   
   std::vector<TString> metRegions={"p_{T}^{miss}<40GeV","40GeV<p_{T}^{miss}<80GeV","80GeV<p_{T}^{miss}<120GeV","120GeV<p_{T}^{miss}<160GeV","160GeV<p_{T}^{miss}<230GeV","230GeV<p_{T}^{miss}",""};
   
   int i=0;
   // ~for(TString target :{"","_ptnunu"}){
   // ~for(TString target :{""}){
   for(TString target :{"","_phi"}){
      
      double axis_limit=500;
      double axis_limit_low=0;
      double axis_limit_diff=100;
      int numBins_diff=100;
      if(target=="_phi") {
         axis_limit=3.2;
         axis_limit_low=-3.2;
         axis_limit_diff=1.6;
         numBins_diff=400;
      }
      
      TH1F genMET_histCombined("","",250,axis_limit_low,axis_limit);
      TH1F PuppiMET_histCombined("","",250,axis_limit_low,axis_limit);
      TH1F PuppiMETcorr_histCombined("","",250,axis_limit_low,axis_limit);
      TH1F PFMET_histCombined("","",250,axis_limit_low,axis_limit);
      
      TH1F PuppiMET_diffCombined("","",numBins_diff,-1.*axis_limit_diff,axis_limit_diff);
      TH1F PuppiMETcorr_diffCombined("","",numBins_diff,-1.*axis_limit_diff,axis_limit_diff);
      TH1F PFMET_diffCombined("","",numBins_diff,-1.*axis_limit_diff,axis_limit_diff);
      
      for (TString bin :{"bin1","bin2","bin3","bin4","bin5","bin6","combined"}){
         //Distribution
         can.cd();
         TH1F* genMET_hist;
         TH1F* PuppiMET_hist;
         TH1F* PuppiMETcorr_hist;
         TH1F* PFMET_hist;
         
         if(bin!="combined") {
            if(target=="" || target=="_phi") genMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/GenMET"+target+"_"+bin);
            else genMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/Ptnunu_"+bin);
            PuppiMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMET"+target+"_"+bin);
            PuppiMETcorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMETcorr"+target+"_"+bin);
            PFMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PFMET"+target+"_"+bin);
         }
         
         genMET_hist->Rebin(2);
         PuppiMET_hist->Rebin(2);
         PuppiMETcorr_hist->Rebin(2);
         PFMET_hist->Rebin(2);
         // ~genMET_hist->Rebin(20);
         // ~PuppiMET_hist->Rebin(20);
         // ~PuppiMETcorr_hist->Rebin(20);
         // ~PFMET_hist->Rebin(20);
         
         hist::mergeOverflow(*genMET_hist,true);
         hist::mergeOverflow(*PuppiMET_hist,true);
         hist::mergeOverflow(*PuppiMETcorr_hist,true);
         hist::mergeOverflow(*PFMET_hist,true);
         
         if(bin!="combined") {
            genMET_histCombined.Add(genMET_hist);     //add for combined plot
            PuppiMET_histCombined.Add(PuppiMET_hist);
            PuppiMETcorr_histCombined.Add(PuppiMETcorr_hist);
            PFMET_histCombined.Add(PFMET_hist);
         }
         else{
            genMET_hist=&genMET_histCombined;
            PuppiMET_hist=&PuppiMET_histCombined;
            PuppiMETcorr_hist=&PuppiMETcorr_histCombined;
            PFMET_hist=&PFMET_histCombined;
         }
         
         genMET_hist->SetLineColor(kGreen);
         PuppiMET_hist->SetLineColor(kRed);
         PuppiMETcorr_hist->SetLineColor(kBlue);
         PFMET_hist->SetLineColor(kBlack);
         
         genMET_hist->SetStats(0);
         if(target=="")genMET_hist->GetXaxis()->SetTitle("MET (GeV)");
         else genMET_hist->GetXaxis()->SetTitle("#phi(MET)");
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
         if(target=="" || target=="_phi") le.append(*genMET_hist,"genMET","l");
         else le.append(*genMET_hist,"pTnunu","l");
         TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TLatex label=gfx::cornerLabel(metRegions[i],1);
         label.Draw();
         saver.save(can,inputFile+"/dist"+target+"_"+bin,true,true);
         
         //Resolution
         TH1F* diff_hist;
         TH1F* diffPF_hist;
         TH1F* diffCorr_hist;
         
         if(bin!="combined") {
            diff_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diff"+target+"_"+bin);
            diffPF_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffPF"+target+"_"+bin);
            diffCorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffcorr"+target+"_"+bin);
         }
         
         if(target==""){
            diff_hist->Rebin(2);
            diffPF_hist->Rebin(2);
            diffCorr_hist->Rebin(2);
            // ~diff_hist->Rebin(20);
            // ~diffPF_hist->Rebin(20);
            // ~diffCorr_hist->Rebin(20);
         }
         
         if(bin!="combined") {
            PuppiMET_diffCombined.Add(diff_hist);     //add for combined plot
            PuppiMETcorr_diffCombined.Add(diffPF_hist);
            PFMET_diffCombined.Add(diffCorr_hist);
         }
         else{
            diff_hist=&PuppiMET_diffCombined;
            diffPF_hist=&PuppiMETcorr_diffCombined;
            diffCorr_hist=&PFMET_diffCombined;
         }
         
         diff_hist->SetLineColor(kRed);
         diffPF_hist->SetLineColor(kBlack);
         diffCorr_hist->SetLineColor(kBlue);
         
         diff_hist->SetStats(0);
         diff_hist->GetYaxis()->SetRangeUser(0,1.4*diffCorr_hist->GetMaximum());
         if(target=="")diff_hist->GetXaxis()->SetTitle("genMET-recoMET (GeV)");
         else if (target=="_phi") diff_hist->GetXaxis()->SetTitle("#Delta#phi(genMET,recoMET)");
         diff_hist->GetYaxis()->SetTitle("Events/Bin");
         
         diff_hist->Draw("hist");
         diffPF_hist->Draw("hist same");
         diffCorr_hist->Draw("hist same");
         
         gfx::LegendEntries le2;
         le2.append(*diff_hist,TString::Format("PuppiMET(#mu=%.2f #sigma=%.2f)",diff_hist->GetMean(),diff_hist->GetRMS()),"l");
         le2.append(*diffPF_hist,TString::Format("PFMET(#mu=%.2f #sigma=%.2f)",diffPF_hist->GetMean(),diffPF_hist->GetRMS()),"l");
         le2.append(*diffCorr_hist,TString::Format("DNN(#mu=%.2f #sigma=%.2f)",diffCorr_hist->GetMean(),diffCorr_hist->GetRMS()),"l");
         TLegend leg2=le2.buildLegend(.45,.75,1-1.5*gPad->GetRightMargin(),-1,1);
         leg2.Draw();
         
         TLatex label2=gfx::cornerLabel(metRegions[i],1);
         label2.Draw();
         saver.save(can,inputFile+"/diff"+target+"_"+bin,true,true);
         
         //DNN Output
         if(target==""){
            if(bin!="combined") {
               canCombined.cd();
               gPad->SetLogy();
               TH1F* DNN_output = histReader.read<TH1F>("binningUnfolding_DNN/DNN/DNNoutput_"+bin);
               
               hist::mergeOverflow(*DNN_output,true);
               
               DNN_output->Scale(1./DNN_output->Integral());
               
               DNN_output->SetLineColor(Color::next());
               DNN_output->SetStats(0);
               DNN_output->GetYaxis()->SetRangeUser(0.001,0.3);
               DNN_output->GetXaxis()->SetTitle("DNN Output");
               DNN_output->GetYaxis()->SetTitle("normalized distributions");
               
               if (i==0) DNN_output->Draw("hist");
               else DNN_output->Draw("same hist");
               
               leCombined.append(*DNN_output,bin,"l");
            }
         }
         i++;
      }
      i=0;
   }
   
   canCombined.cd();
   TLegend leg=leCombined.buildLegend(.45,.75,1-1.5*gPad->GetRightMargin(),-1,2);
   leg.Draw();
   saver.save(canCombined,inputFile+"/DNN_Output",true,true);
   
}
