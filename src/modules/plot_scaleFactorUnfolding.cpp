//Script to plot purity and stability as a function of an additional scaleFactor applied to ptmiss

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_scaleFactorUnfolding");
   TCanvas can;
   TCanvas can_multi;
   gfx::LegendEntries le_multi;
   int color=2;
   
   //Plot all graphs for each sample
   for (TString sSample :{"","dilepton","MadGraph"}){
      
      can.cd();
      io::RootFileReader histReader((sSample=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sSample+"%.1f.root",cfg.processFraction*100));
      
      TGraph* purity_last = histReader.read<TGraph>("binningUnfolding/scaleFactor_Purity_lastbin");
      TGraph* purity_mean = histReader.read<TGraph>("binningUnfolding/scaleFactor_Purity_mean");
      TGraph* stability_last = histReader.read<TGraph>("binningUnfolding/scaleFactor_Stability_lastbin");
      TGraph* stability_mean = histReader.read<TGraph>("binningUnfolding/scaleFactor_Stability_mean");
      
      purity_last->SetLineWidth(2);
      purity_mean->SetLineWidth(2);
      stability_last->SetLineWidth(2);
      stability_mean->SetLineWidth(2);
      
      purity_last->SetLineColor(kBlue);
      purity_last->SetLineStyle(kDashed);
      purity_mean->SetLineColor(kBlue);
      stability_last->SetLineColor(kRed);
      stability_last->SetLineStyle(kDashed);
      stability_mean->SetLineColor(kRed);
      
      stability_mean->SetMaximum(1.5*stability_last->GetMaximum());
      stability_mean->GetXaxis()->SetTitle("Scale Factor");
      stability_mean->GetXaxis()->SetRangeUser(0.5,1.04);
      stability_mean->GetYaxis()->SetTitle("Purity/Stability");
      
      stability_mean->Draw("AL");
      stability_last->Draw("same L");
      purity_last->Draw("same L");
      purity_mean->Draw("same L");
      
      gfx::LegendEntries le;
      le.append(*purity_last,"Purity last bin","l");
      le.append(*purity_mean,"Purity mean","l");
      le.append(*stability_last,"Stability last bin","l");
      le.append(*stability_mean,"Stability mean","l");
      TLegend leg=le.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
      leg.Draw();
      
      TLatex label=gfx::cornerLabel(sSample,1);
      label.Draw();
      saver.save(can,"ScaleFactor_"+sSample,true,true);
      can.Clear();
      
      can_multi.cd();
      stability_mean->SetLineColor(color);
      purity_mean->SetLineColor(color);
      
      if (sSample=="") le_multi.append(*stability_mean,"Powheg incl.","l");
      else if (sSample=="dilepton") le_multi.append(*stability_mean,"Powheg dilep.","l");
      else le_multi.append(*stability_mean,"MadGraph dilep.","l");
      
      stability_mean->GetXaxis()->SetRangeUser(0.8,1.04);
      if (sSample=="") stability_mean->Draw("AL");
      else stability_mean->Draw("same L");
      purity_mean->Draw("same L");
      
      color++;
   }
   
   can_multi.cd();
   TLine aline(0.91,0.,0.91,1.);
   aline.Draw();
   TLegend leg=le_multi.buildLegend(.4,.7,1-gPad->GetRightMargin(),-1,2);
   leg.Draw();
   saver.save(can_multi,"ScaleFactor_compareSample",true,true);
   
}
