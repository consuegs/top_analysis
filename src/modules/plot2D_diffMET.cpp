//Script to plot 2D distributions for different METS and others variables

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_diffMET");
   io::RootFileReader histReader(TString::Format("binningUnfolding_diLepton%.1f.root",cfg.processFraction*100),"binningUnfolding");
   TCanvas can;
   
   for (TString sVar :{"METvsGenMET", "GenMETvsPTnunu"}) {
      
      TH2F *hist;
      hist = (TH2F*) histReader.read<TH2F>(sVar);
      hist->Scale(1.0/(hist->Integral()));
            
      // Plotting stuff in the following
      TGraph diag;
      diag.SetPoint(1,0,0);
      diag.SetPoint(2,600,600);
      can.cd();
      can.SetLogz();
      
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.11);
      
      hist->GetYaxis()->SetTitleOffset(1.3);
      hist->GetXaxis()->SetTitleOffset(0.9);
      hist->GetZaxis()->SetTitleOffset(1.3);
      hist->GetYaxis()->SetTitleSize(0.05);
      hist->GetXaxis()->SetTitleSize(0.05);
      hist->GetZaxis()->SetTitleSize(0.05);
      hist->GetYaxis()->SetLabelSize(0.04);
      hist->GetXaxis()->SetLabelSize(0.04);
      hist->GetZaxis()->SetLabelSize(0.04);
      
      hist->GetZaxis()->SetTitle("normalized distribution");
      
      hist->SetMinimum(0.000001);
      hist->SetMaximum(0.1);

      hist->SetStats(false);
      hist->Draw("colz");
      diag.Draw("same");
      
      TString cat_label="all";
      TLatex label=gfx::cornerLabel(cat_label,1);
      label.Draw();
      
      can.RedrawAxis();
      TString plotLoc=sVar;
      saver.save(can,plotLoc,true,true);
      can.Clear();
   }
   
   // This part was used before, when distributions produced the resolution plots
   // ~io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot2D_diffMET");
   // ~io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   // ~TCanvas can;
   
   // ~for (TString sSelection : {"genParticles","genParticles_Met200","baseline","baseline_Met200"}){
      // ~for (TString sSample :{"TTbar","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ",
      // ~"T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){//,"TTbar_diLepton"}){
         
         // ~std::vector<TString> Plots_2d={"2d_nunuVSgenMet","MetVSgenMet","2d_nunuVSMet","2d_nunuVSDMgenMet","MetVSDMgenMet","2d_MT2VSgenMT2","2d_MT2VSgenMT2neutrino",
            // ~"2d_MT_l1METVSgenMT","2d_MT_l1METVSgenMTneutrino","2d_dPhiMetNearLep_genResponse","2d_dPhiMetNearLep_Response"};
         // ~if (sSelection == "baseline" or sSelection == "baseline_Met200") Plots_2d={"2d_MT2VSdPhiMetNearLep","2d_MetVSdPhiMetNearLep"};
         
         // ~for (TString sVar : Plots_2d) {
            
            // ~TH2F hist_add;
            // ~TH2F *hist;
            
            // ~for (TString cat:{"ee","emu","mumu","all"}){
               
               // ~//Add histograms for category all, otherwise use only hist (not hist_add)
               // ~if (cat!="all"){
                  // ~hist = (TH2F*) histReader.read<TH2F>(sSelection+"/"+cat+"/"+sVar+"/"+sSample);
                  // ~if (cat=="ee") hist_add=(TH2F) *(histReader.read<TH2F>(sSelection+"/"+cat+"/"+sVar+"/"+sSample));
                  // ~else hist_add.Add(hist);
               // ~}
               // ~else {
                  // ~hist = &hist_add;
               // ~}
                  
               
               // ~hist->RebinX(2);  //Use wider binning to increase statistics for BSM samples
               // ~hist->RebinY(2);
               
               // ~for (TString norm:{"","LineNorm"}){    //Plot nominal response and response normalized in each row
                  
                  // ~if (norm==""){
                     // ~hist->Scale(1.0/(hist->Integral()));
                  // ~}
                  // ~else{
                     // ~//Normalize each individual row of diagram
                     // ~float sum;
                     // ~for (int y=1; y<=hist->GetNbinsY(); y++){
                        // ~sum=hist->Integral(1,hist->GetNbinsX(),y,y);
                        // ~if (sum==0) continue;
                        // ~for (int x=1; x<=hist->GetNbinsX(); x++){
                           // ~hist->SetBinContent(x,y,hist->GetBinContent(x,y)/sum);
                        // ~}
                     // ~}
                  // ~}
                  
                  // ~// Plotting stuff in the following
                  // ~TGraph diag;
                  // ~diag.SetPoint(1,0,0);
                  // ~diag.SetPoint(2,600,600);
                  // ~can.cd();
                  // ~can.SetLogz();
                  
                  // ~gPad->SetRightMargin(0.2);
                  // ~gPad->SetLeftMargin(0.13);
                  // ~gPad->SetBottomMargin(0.11);
                  
                  // ~hist->GetYaxis()->SetTitleOffset(1.3);
                  // ~hist->GetXaxis()->SetTitleOffset(0.9);
                  // ~hist->GetZaxis()->SetTitleOffset(1.3);
                  // ~hist->GetYaxis()->SetTitleSize(0.05);
                  // ~hist->GetXaxis()->SetTitleSize(0.05);
                  // ~hist->GetZaxis()->SetTitleSize(0.05);
                  // ~hist->GetYaxis()->SetLabelSize(0.04);
                  // ~hist->GetXaxis()->SetLabelSize(0.04);
                  // ~hist->GetZaxis()->SetLabelSize(0.04);
                  
                  // ~if (norm==""){
                     // ~hist->GetZaxis()->SetTitle("normalized distribution");
                     
                     // ~hist->SetMinimum(0.000001);
                     // ~hist->SetMaximum(0.1);
                  // ~}
                  // ~else{
                     // ~hist->GetZaxis()->SetTitle("line normalized distribution");
                     
                     // ~hist->SetMinimum(0.00001);
                     // ~hist->SetMaximum(1.);
                  // ~}

                  // ~hist->SetStats(false);
                  // ~hist->Draw("colz");
                  // ~diag.Draw("same");
                  
                  // ~TString cat_label=cat;
                  // ~if (cat=="emu") cat_label="e#mu";
                  // ~else if (cat=="mumu") cat_label="#mu#mu";
                  // ~else if (cat=="all") cat_label="all";
                  // ~TString labelString=cat_label+"    "+sSample;
                  // ~if (sSelection=="genParticles_Met200" or sSelection=="baseline_Met200") labelString=cat_label+"    "+sSample+"    p_{T}^{miss}>200 GeV";
                  // ~TLatex label=gfx::cornerLabel(labelString,1);
                  // ~label.Draw();
                  
                  // ~can.RedrawAxis();
                  // ~TString plotLoc=sSample+"/"+sVar+norm+"/"+cat;
                  // ~if (sSelection=="genParticles_Met200" or sSelection=="baseline_Met200") plotLoc="Met200/"+sSample+"/"+sVar+norm+"/"+cat;
                  // ~saver.save(can,plotLoc,true,true);
                  // ~can.Clear();
               // ~}
            // ~}
            
         // ~}
      // ~}
   // ~}
}
