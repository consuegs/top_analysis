//Script to perform 2D unfolding
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TUnfoldDensity.h>
#include <TUnfold.h>
Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileReader histReader(TString::Format("binningUnfolding_dilepton%.1f.root",cfg.processFraction*100),"binningUnfolding");
   io::RootFileReader histReaderMadGraph(TString::Format("binningUnfolding_MadGraph%.1f.root",cfg.processFraction*100),"binningUnfolding");
   io::RootFileReader histReaderSignal(TString::Format("binningUnfolding_T2tt_650_350%.1f.root",cfg.processFraction*100),"binningUnfolding");
   io::RootFileSaver saver_hist(TString::Format("unfolding%.1f.root",cfg.processFraction*100),TString::Format("unfolding%.1f",cfg.processFraction*100),false);
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("unfolding%.1f",cfg.processFraction*100));
   
   TCanvas can2D;
   TH2* Rho2D;
   TH1* unfolded;
   
   TH2F* migrMatrix=(TH2F*)histReader.read<TH2F>("Binning_Met1_Phi1/Unfolding/response");
   TH1F* fakeData=(TH1F*)histReaderMadGraph.read<TH1F>("Binning_Met1_Phi1/Unfolding/recoDistribution");
   fakeData->Add((TH1F*)histReaderSignal.read<TH1F>("Binning_Met1_Phi1/Unfolding/recoDistribution"));
   TH1F* realDis=(TH1F*)histReaderMadGraph.read<TH1F>("Binning_Met1_Phi1/Unfolding/trueDistribution");
   TH1F* realDis_signal=(TH1F*)histReaderSignal.read<TH1F>("Binning_Met1_Phi1/Unfolding/trueDistribution");
   TUnfoldDensity unfold(migrMatrix,TUnfold::kHistMapOutputVert);
   if(unfold.SetInput(fakeData)>=10000) {
      std::cout<<"Unfolding result may be wrong\n";
   }
   
   // ~// scan L curve and find best point for regularisation
   // ~Int_t nScan=30;
   // ~// use automatic L-curve scan: start with taumin=taumax=0.0
   // ~Double_t tauMin=0.0;
   // ~Double_t tauMax=0.0;
   // ~Int_t iBest;
   // ~TSpline *logTauX,*logTauY;
   // ~TGraph *lCurve;
   // ~iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
   
   unfold.DoUnfold(0.);

   unfolded=unfold.GetOutput("Unfolded");
   Rho2D=unfold.GetRhoIJtotal("Rho2D");
   double tau=unfold.GetTau();
   
   std::cout<<tau<<std::endl;
   
   saver_hist.save(*fakeData,"fakeData");
   saver_hist.save(*unfolded,"unfolded");
   saver_hist.save(*realDis,"true");
   saver_hist.save(*realDis_signal,"true_signal");
   saver_hist.save(*Rho2D,"Rho2D");
   
   
   
   //Plot Unfolded, real distribution and pseudo unfolded
   TCanvas can("c1","c1", 1200, 600);
   // ~gfx::SplitCan can;
   // ~can.can_.SetWindowSize(1200,600);
   can.SetLogy();
   // ~can.pU_.SetLogy();
   
   Double_t xbins[13] = {0,60,120,230,400,460,520,630,800,860,920,1030,1200};
   unfolded->SetBins(12,xbins);
   realDis->SetBins(12,xbins);
   realDis_signal->SetBins(12,xbins);
   unfolded->GetXaxis()->SetTickLength(0.);
   unfolded->GetXaxis()->SetTitleOffset(1.5);
   unfolded->GetXaxis()->CenterLabels(false);
   unfolded->GetXaxis()->SetBinLabel(1,"0-60");
   unfolded->GetXaxis()->SetBinLabel(2,"60-120");
   unfolded->GetXaxis()->SetBinLabel(3,"120-230");
   unfolded->GetXaxis()->SetBinLabel(4,">230");
   unfolded->GetXaxis()->SetBinLabel(5,"0-60");
   unfolded->GetXaxis()->SetBinLabel(6,"60-120");
   unfolded->GetXaxis()->SetBinLabel(7,"120-230");
   unfolded->GetXaxis()->SetBinLabel(8,">230");
   unfolded->GetXaxis()->SetBinLabel(9,"0-60");
   unfolded->GetXaxis()->SetBinLabel(10,"60-120");
   unfolded->GetXaxis()->SetBinLabel(11,"120-230");
   unfolded->GetXaxis()->SetBinLabel(12,">230");
   unfolded->LabelsOption("v");
   unfolded->SetMaximum(4*unfolded->GetMaximum());
   unfolded->SetMinimum(10);
   unfolded->SetLineColor(kBlack);
   unfolded->SetTitle(";genMET (GeV);arbitrary unit");
   unfolded->SetStats(false);
   
   unfolded->Scale(1.0/(10e5));
   realDis->Scale(1.0/(10e5));
   realDis_signal->Scale(1.0/(10e5));
   
   unfolded->Draw("pe1");
   realDis->SetLineColor(kRed-6);
   realDis->SetFillColor(kRed-6);
   realDis->Draw("hist same");
   realDis_signal->SetLineColor(kBlue);
   realDis_signal->Draw("hist same");
   
   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.03);
   //~ atext->SetTextFont(42);
   aline->SetLineWidth(2);
   aline->DrawLine(400,10,400,unfolded->GetMaximum());
   aline->DrawLine(800,10,800,unfolded->GetMaximum());
   aline->SetLineStyle(2);
   atext->DrawLatex(50,0.5*unfolded->GetMaximum(),"0<|#Delta#phi|(genMET,nearest l)<0.7");
   atext->DrawLatex(450,0.5*unfolded->GetMaximum(),"0.7<|#Delta#phi|(genMET,nearest l)<1.4");
   atext->DrawLatex(850,0.5*unfolded->GetMaximum(),"1.4<|#Delta#phi|(genMET,nearest l)<3.14");
   
   
   gfx::LegendEntries legE;
   legE.append(*unfolded,"Unfolded","pe");
   legE.append(*realDis,"MC true ttbar","l");
   legE.append(*realDis_signal,"MC true signal","l");
   TLegend leg=legE.buildLegend(.2,.25,0.35,.4,1);
   leg.SetTextSize(0.035);
   leg.Draw();
   
   // ~can.pU_.RedrawAxis();
   unfolded->Draw("axis same");
      
   saver.save(can,"/MadGraph",true);
   
   /*
   //Plot Rho2D
   can2D.cd();
   gPad->SetRightMargin(0.2);
   gPad->SetLeftMargin(0.13);
   gPad->SetBottomMargin(0.11);
   
   Rho2D->GetYaxis()->SetTitleOffset(1.3);
   Rho2D->GetXaxis()->SetTitleOffset(0.9);
   Rho2D->GetZaxis()->SetTitleOffset(1.3);
   Rho2D->GetYaxis()->SetTitleSize(0.05);
   Rho2D->GetXaxis()->SetTitleSize(0.05);
   Rho2D->GetZaxis()->SetTitleSize(0.05);
   Rho2D->GetYaxis()->SetLabelSize(0.04);
   Rho2D->GetXaxis()->SetLabelSize(0.04);
   Rho2D->GetZaxis()->SetLabelSize(0.04);
   Rho2D->GetZaxis()->SetLabelOffset(0.01);

   Rho2D->SetStats(false);
   Rho2D->SetTitle(";genMET (GeV);genMET (GeV)");
   Rho2D->Draw("colz");
   
   TLatex label2=gfx::cornerLabel(cat_label,1);
   label2.Draw();
   
   saver.save(can2D,sSample+"/Rho2D/"+cat,true,true);
   */
}
