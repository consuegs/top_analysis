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
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_InclusiveDNN_Patched_LR_EP%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_genMetTarget_normDistr%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_test_cppflow_new%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_test_cppflow_new2%.1f.root",cfg.processFraction*100);
   // ~TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_diffTarget_normDistr_100EP_%.1f.root",cfg.processFraction*100);
   TString inputFile=TString::Format("binningUnfolding_diLepton_amcatnlo_pTnunuTarget_normDistr%.1f.root",cfg.processFraction*100);
   
   io::RootFileReader histReader(inputFile);
   
   std::vector<float> axisLimits_low={0,0,0,50,50,50};
   std::vector<float> axisLimits_high={150,150,200,250,300,500};
   
   std::vector<TString> metRegions={"p_{T}^{miss}<40GeV","40GeV<p_{T}^{miss}<80GeV","80GeV<p_{T}^{miss}<120GeV","120GeV<p_{T}^{miss}<160GeV","160GeV<p_{T}^{miss}<230GeV","230GeV<p_{T}^{miss}",""};
   
   TH1F genMET_histCombined("","",250,0,500);
   TH1F PuppiMET_histCombined("","",250,0,500);
   TH1F PuppiMETcorr_histCombined("","",250,0,500);
   TH1F PFMET_histCombined("","",250,0,500);
   
   TH1F PuppiMET_diffCombined("","",100,-100,100);
   TH1F PuppiMETcorr_diffCombined("","",100,-100,100);
   TH1F PFMET_diffCombined("","",100,-100,100);
   
   int i=0;
   for(TString target :{"","_ptnunu"}){
   // ~for(TString target :{""}){
      for (TString bin :{"bin1","bin2","bin3","bin4","bin5","bin6","combined"}){
         //Distribution
         can.cd();
         TH1F* genMET_hist;
         TH1F* PuppiMET_hist;
         TH1F* PuppiMETcorr_hist;
         TH1F* PFMET_hist;
         
         if(bin!="combined") {
            if(target=="") genMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/GenMET_"+bin);
            else genMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/Ptnunu_"+bin);
            PuppiMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMET_"+bin);
            PuppiMETcorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PuppiMETcorr_"+bin);
            PFMET_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/PFMET_"+bin);
         }
         
         genMET_hist->Rebin(2);
         PuppiMET_hist->Rebin(2);
         PuppiMETcorr_hist->Rebin(2);
         PFMET_hist->Rebin(2);
         // ~genMET_hist->Rebin(20);
         // ~PuppiMET_hist->Rebin(20);
         // ~PuppiMETcorr_hist->Rebin(20);
         // ~PFMET_hist->Rebin(20);
         
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
         if(target=="") le.append(*genMET_hist,"genMET","l");
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
         TH1F* diffScaled_hist;
         
         if(bin!="combined") {
            diff_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diff"+target+"_"+bin);
            diffPF_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffPF"+target+"_"+bin);
            diffCorr_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffcorr"+target+"_"+bin);
            diffScaled_hist = histReader.read<TH1F>("binningUnfolding_DNN/DNN/diffscaled"+target+"_"+bin);
         }
         
         diff_hist->Rebin(2);
         diffPF_hist->Rebin(2);
         diffCorr_hist->Rebin(2);
         diffScaled_hist->Rebin(2);
         // ~diff_hist->Rebin(20);
         // ~diffPF_hist->Rebin(20);
         // ~diffCorr_hist->Rebin(20);
         // ~diffScaled_hist->Rebin(20);
         
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
         diffScaled_hist->SetLineColor(kGray);
         
         diff_hist->SetStats(0);
         diff_hist->GetYaxis()->SetRangeUser(0,1.4*diffCorr_hist->GetMaximum());
         if(target=="")diff_hist->GetXaxis()->SetTitle("genMET-recoMET (GeV)");
         else diff_hist->GetXaxis()->SetTitle("pTnunu-recoMET (GeV)");
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
         saver.save(can,inputFile+"/diff"+target+"_"+bin,true,true);
         
         //DNN Output
         if(target==""){
            if(bin!="combined") {
               canCombined.cd();
               gPad->SetLogy();
               TH1F* DNN_ouput = histReader.read<TH1F>("binningUnfolding_DNN/DNN/DNNoutput_"+bin);
               
               DNN_ouput->Scale(1./DNN_ouput->Integral());
               
               DNN_ouput->SetLineColor(Color::next());
               DNN_ouput->SetStats(0);
               DNN_ouput->GetYaxis()->SetRangeUser(0.001,0.3);
               DNN_ouput->GetXaxis()->SetTitle("DNN Output");
               DNN_ouput->GetYaxis()->SetTitle("normalized distributions");
               
               if (i==0) DNN_ouput->Draw("hist");
               else DNN_ouput->Draw("same hist");
               
               leCombined.append(*DNN_ouput,bin,"l");
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
