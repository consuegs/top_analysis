//First script to try some unfolding stuff and get used to the unfolding framework
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
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileSaver saver_hist(TString::Format("unfolding%.1f.root",cfg.processFraction*100),TString::Format("unfolding%.1f",cfg.processFraction*100),false);
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("unfolding%.1f",cfg.processFraction*100));
   
   TCanvas can;
   TCanvas can2D;
   TH2* Rho2D;
   TH1* unfolded;
      
   for (TString cat: {"ee","emu","mumu"}){
      TString sSelection="genParticles/"+cat+"/MetVSgenMet";
      TString sSelection2="baseline/"+cat+"/met";
      TString sSelection3="genParticles/"+cat+"/genMet";
      for (TString sSample:{"T1tttt_1200_800","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
         if (sSample.Contains("DM")) sSelection3="genParticles/"+cat+"/DMgenMet";
         TH2F* migrMatrix=(TH2F*)histReader.read<TH2F>(sSelection+"/TTbar");
         TH1F* fakeData=(TH1F*)histReader.read<TH1F>(sSelection2+"/"+sSample);
         TH1F* realDis=(TH1F*)histReader.read<TH1F>(sSelection3+"/"+sSample);
         migrMatrix->RebinY(10);
         migrMatrix->RebinX(2);
         fakeData->Rebin(2);
         realDis->Rebin(10);
         TUnfoldDensity unfold(migrMatrix,TUnfold::kHistMapOutputVert);
         if(unfold.SetInput(fakeData)>=10000) {
            std::cout<<"Unfolding result may be wrong\n";
         }
         
         // scan L curve and find best point for regularisation
         //~ Int_t nScan=30;
         //~ // use automatic L-curve scan: start with taumin=taumax=0.0
         //~ Double_t tauMin=0.0;
         //~ Double_t tauMax=0.0;
         //~ Int_t iBest;
         //~ TSpline *logTauX,*logTauY;
         //~ TGraph *lCurve;
         //~ iBest=unfold.ScanLcurve(nScan,tauMin,tauMax,&lCurve,&logTauX,&logTauY);
         
         unfold.DoUnfold(0.);
  
         unfolded=unfold.GetOutput("Unfolded");
         Rho2D=unfold.GetRhoIJtotal("Rho2D");
         double tau=unfold.GetTau();
         
         std::cout<<tau<<std::endl;
         
         saver_hist.save(*unfolded,sSample+"/unfolded/"+cat);
         saver_hist.save(*realDis,sSample+"/genMet/"+cat);
         saver_hist.save(*Rho2D,sSample+"/Rho2D/"+cat);
         
         //PseudoUnfolding with bin-by-bin scale factors
         TH1F* pseudo2;
         TH1F* pseudo5;
         for (int binFactor :{2,10}){
            TH1F* metTTbar=(TH1F*)histReader.read<TH1F>(sSelection2+"/TTbar");
            TH1F* genMetTTbar=(TH1F*)histReader.read<TH1F>(sSelection3+"/TTbar");
            metTTbar->Rebin(binFactor);
            genMetTTbar->Rebin(binFactor);
            
            TH1F* scaleFactors=genMetTTbar;
            scaleFactors->Divide(metTTbar);
            
            TH1F* metBSM=(TH1F*)histReader.read<TH1F>(sSelection2+"/"+sSample);
            metBSM->Rebin(binFactor);
            metBSM->Multiply(scaleFactors);
            if (binFactor==2){
               metBSM->Rebin(5);
               pseudo2=metBSM;
            }
            else pseudo5=metBSM;
            saver_hist.save(*metBSM,sSample+"/PseudoUnfolding_"+std::to_string(binFactor)+"x"+std::to_string(binFactor)+"/"+cat);
         }
         
         
         //Plot Unfolded, real distribution and pseudo unfolded
         can.cd();
         unfolded->SetLineColor(kBlack);
         unfolded->SetTitle(";genMET (GeV);Events/Bin");
         unfolded->SetStats(false);
         unfolded->Draw("pe1");
         realDis->SetLineColor(kBlue);
         realDis->Draw("hist same");
         pseudo2->SetLineColor(kRed);
         pseudo5->SetLineColor(kGreen);
         // ~pseudo2->Draw("pe1 same");
         // ~pseudo5->Draw("pe1 same");
         gfx::LegendEntries legE;
         legE.append(*unfolded,"Unfolded pseudo Data","pe");
         // ~legE.append(*pseudo2,"Scaled pseudo Data (50x50)","pe");
         // ~legE.append(*pseudo5,"Scaled pseudo Data (10x10)","pe");
         legE.append(*realDis,"MC True","l");
         TLegend leg=legE.buildLegend(.5,.65,0.75,.9,1);
         leg.SetTextSize(0.035);
         leg.Draw();
         
         TString cat_label=cat;
         if (cat=="emu") cat_label="e#mu";
         else if (cat=="mumu") cat_label="#mu#mu";
         else if (cat=="all") cat_label="all";
         cat_label+="   "+sSample;
         TLatex label=gfx::cornerLabel(cat_label,1);
         label.Draw();
            
         saver.save(can,sSample+"/unfolded/"+cat,true,true);
         
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
         
      }
   }
}
