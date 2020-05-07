//Script to plot distributions for fakes and nofake events

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
#include <iomanip>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_diff_fakes");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("diff_fakes%.1f",cfg.processFraction*100));
   TCanvas can;
   
   std::vector<TString> msPresel_vVars={"met","met1000","pTlep1","pTlep2","etalep1","etalep2","pTbJet","etabJet","pTjet1","pTjet2","etajet1","etajet2","deepCSV"};
   
   for (TString sVar:msPresel_vVars){
         
      TH1F* hist_nofakes;
      TH1F* hist_fakes;
      TH1F* hist_tau;
      
      hist_nofakes=(TH1F*) histReader.read<TH1F>("NoFakes/"+sVar+"/TTbar_diLepton");
      hist_fakes=(TH1F*) histReader.read<TH1F>("Fakes/"+sVar+"/TTbar_diLepton");
      hist_tau=(TH1F*) histReader.read<TH1F>("Taus/"+sVar+"/TTbar_diLepton");
      // ~hist_nofakes=(TH1F*) histReader.read<TH1F>("NoFakes/"+sVar+"/TTbar");
      // ~hist_fakes=(TH1F*) histReader.read<TH1F>("Fakes/"+sVar+"/TTbar");
      
      if (sVar=="deepCSV") {
            
            std::cout<<hist_nofakes->Integral()<<"   "<<std::setprecision(7)<<hist_nofakes->GetEntries()<<std::endl;
            std::cout<<hist_fakes->Integral()<<"   "<<hist_fakes->GetEntries()<<std::endl;
            std::cout<<hist_tau->Integral()<<"   "<<hist_tau->GetEntries()<<std::endl;
      }
         
      hist_nofakes->Scale(1.0/(hist_nofakes->Integral()));   //Normalize the hist to the integral
      hist_fakes->Scale(1.0/(hist_fakes->Integral()));
      hist_tau->Scale(1.0/(hist_tau->Integral()));
      
      can.cd();
      can.SetLogz();
      
      gPad->SetLeftMargin(0.15);
      hist_fakes->GetYaxis()->SetTitle("normalized distribution");
      hist_nofakes->SetLineColor(kBlue);
      hist_fakes->SetLineColor(kRed);
      hist_tau->SetLineColor(kGreen+2);
      
      hist_fakes->SetStats(0);
      hist_nofakes->SetStats(0);
      hist_tau->SetStats(0);
      
      hist_fakes->Draw("hist");
      hist_nofakes->Draw("same hist");
      hist_tau->Draw("same hist");
      
      hist_fakes->SetMaximum(1.3*std::max(hist_fakes->GetMaximum(),hist_nofakes->GetMaximum()));
      
      gfx::LegendEntries legE;
      legE.append(*hist_nofakes,"No Fakes","l");
      legE.append(*hist_fakes,"Fakes","l");
      legE.append(*hist_tau,"Taus","l");
      TLegend leg=legE.buildLegend(.5,.65,0.75,.9,1);
      leg.SetTextSize(0.035);
      leg.Draw();
      
      can.RedrawAxis();
      TString plotLoc=sVar;
      // ~TString plotLoc="InclusivePowHeg/"+sVar;
      saver.save(can,plotLoc,true,true);
      can.Clear();
   }
}
