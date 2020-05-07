//Script to plot studies from diff_MET script

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_metStudies");
   
   TCanvas can;
   can.cd();
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()));
   
   for (TString var :{"GenMetDiffMETRel_dPhiMETLep","MetSig_dPhiMETLep"}){
      
      // ~TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"/TTbar_diLepton");
      // ~TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"/T2tt_650_350");
      // ~TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"/TTbar_diLepton");
      // ~TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"/T2tt_650_350");
      // ~TH2F* TTbar_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"_Puppi/TTbar_diLepton");
      // ~TH2F* OtherSample_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met120_230/"+var+"_Puppi/T2tt_650_350");
      TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"/TTbar_diLepton");
      TH2F* OtherSample_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"/T2tt_650_350");
      TH2F* TTbar_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Puppi/TTbar_diLepton");
      TH2F* OtherSample_2D_Puppi = histReader.read<TH2F>("diff_MET100.0/baseline_met230/"+var+"_Puppi/T2tt_650_350");
      
      TProfile TTbar_profile;
      TProfile OtherSample_profile;
      TProfile TTbar_profile_Puppi;
      TProfile OtherSample_profile_Puppi;
      TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
      OtherSample_profile=*(OtherSample_2D->ProfileX("Profile"));
      TTbar_profile_Puppi=*(TTbar_2D_Puppi->ProfileX("Profile"));
      OtherSample_profile_Puppi=*(OtherSample_2D_Puppi->ProfileX("Profile"));
      
      TTbar_profile.SetLineColor(kRed-6);
      OtherSample_profile.SetLineColor(kBlue);
      // ~TTbar_profile.SetMarkerSize(0);
      // ~OtherSample_profile.SetMarkerSize(0);
      TTbar_profile_Puppi.SetLineColor(kRed-6);
      OtherSample_profile_Puppi.SetLineColor(kBlue);
      TTbar_profile_Puppi.SetMarkerStyle(24);
      OtherSample_profile_Puppi.SetMarkerStyle(24);
      
      TTbar_profile.SetStats(0);
      if (var=="GenMetDiffMETRel_dPhiMETLep") {
            TTbar_profile.SetMaximum(0.1);
            TTbar_profile.SetMinimum(-0.5);
            TTbar_profile.GetYaxis()->SetTitle("mean[genMET-p_{T}^{miss}]/genMET]");
      }
      else if (var=="MetSig_dPhiMETLep") {
            // ~TTbar_profile.SetMaximum(200);
            TTbar_profile.SetMaximum(400);
            TTbar_profile.SetMinimum(20);
            TTbar_profile.GetYaxis()->SetTitle("mean[metSig]");
      }
      
      gPad->SetLeftMargin(0.13);
      TTbar_profile.GetYaxis()->SetTitleOffset(1.0);
      TTbar_profile.Draw("e1");
      OtherSample_profile.Draw("same e1");
      TTbar_profile_Puppi.Draw("same e1");
      OtherSample_profile_Puppi.Draw("same e1");
      
      gfx::LegendEntries le;
      le.append(TTbar_profile,"TTbar","pl");
      le.append(OtherSample_profile,"T2tt_650_350","pl");
      le.append(TTbar_profile_Puppi,"TTbar(Puppi)","pl");
      le.append(OtherSample_profile_Puppi,"T2tt_650_350(Puppi)","pl");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      // ~TLatex label=gfx::cornerLabel("p_{T}^{miss}>120GeV",1);
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>230GeV",1);
      // ~TLatex label=gfx::cornerLabel("120GeV<p_{T}^{miss}<230GeV",1);
      label.Draw();
      saver.save(can,var,true,true);
   }
}
