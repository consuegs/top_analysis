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
      
      TH2F* TTbar_2D = histReader.read<TH2F>("diff_MET100.0/baseline_met120/"+var+"/TTbar_diLepton");
      
      TProfile TTbar_profile;
      TTbar_profile=*(TTbar_2D->ProfileX("Profile"));
      
      TTbar_profile.SetLineColor(kRed-6);
      OtherSample_profile.SetLineColor(kBlue);
      
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
      
      gfx::LegendEntries le;
      le.append(TTbar_profile,"TTbar","pl");
      TLegend leg=le.buildLegend(.6,.7,1-1.5*gPad->GetRightMargin(),-1,1);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel("p_{T}^{miss}>120GeV",1);
      // ~TLatex label=gfx::cornerLabel("p_{T}^{miss}>230GeV",1);
      // ~TLatex label=gfx::cornerLabel("120GeV<p_{T}^{miss}<230GeV",1);
      label.Draw();
      saver.save(can,var,true,true);
   }
}
