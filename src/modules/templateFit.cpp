// Script to study if a possible signal can still be found after unfolding (currently only based on genMET and not yet pT_nunu+BSM!!!)
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
#include <TFractionFitter.h>

#include <RooGlobalFunc.h>
#include <RooRealVar.h>
#include <RooDataSet.h>
#include <RooDataHist.h>
#include <RooGaussian.h>
#include <RooConstVar.h>
#include <RooChebychev.h>
#include <RooAddPdf.h>
#include <RooWorkspace.h>
#include <RooPlot.h>
#include <RooHistPdf.h>
#include <RooFitResult.h>
Config const &cfg=Config::get();

extern "C"
void run()
{
   
   int scale=10;     //Scale to be used for BSM samples (otherwise the signal yield would just be too small)
   
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),TString::Format("templateFit%.1f",cfg.processFraction*100));
   io::RootFileSaver saver_templateFit(TString::Format("templateFit%.1f.root",cfg.processFraction*100),TString::Format("templateFit%.1f",cfg.processFraction*100));
   
   std::ofstream out;
   out.open("../output/txt/templateFit.txt");   //File to save fitting results
   out.precision(4);
   
   TCanvas can;
   can.SetLogy();
      
   for (TString cat: {"ee","emu","mumu"}){
      TString sSelection="genParticles/"+cat+"/MetVSgenMet";
      TString sSelectionMET="baseline/"+cat+"/met";
      TString sSelectionGenMET="genParticles/"+cat+"/genMet";
      
      TH2F* migrMatrix=(TH2F*)histReader.read<TH2F>(sSelection+"/TTbar");  //Response matrix used for the unfolding
      migrMatrix->RebinY(10);
      migrMatrix->RebinX(2);
      for (TString sSample:{"T1tttt_1200_800","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200"}){
         if (sSample.Contains("DM")) sSelectionGenMET="genParticles/"+cat+"/DMgenMet";    //Use DMgenMET in case of DM samples
         TH1F* genMET_ttbar=(TH1F*)histReader.read<TH1F>(sSelectionGenMET+"/TTbar");
         TH1F* genMET_BSM=(TH1F*)histReader.read<TH1F>(sSelectionGenMET+"/"+sSample);
         TH1F* MET_ttbar=(TH1F*)histReader.read<TH1F>(sSelectionMET+"/TTbar");
         TH1F* MET_BSM=(TH1F*)histReader.read<TH1F>(sSelectionMET+"/"+sSample);
         genMET_ttbar->Rebin(10);
         genMET_BSM->Rebin(10);
         MET_ttbar->Rebin(2);
         MET_BSM->Rebin(2);
         
         //Scale BSM sample
         genMET_BSM->Scale(scale);
         MET_BSM->Scale(scale);
         
         //Summed distributions
         TH1F pseudoData_MET=*(MET_ttbar);
         pseudoData_MET.Add(MET_BSM);
         TH1F pseudoData_genMET=*(genMET_ttbar);
         pseudoData_genMET.Add(genMET_BSM);
         
         //Unfolding
         TUnfoldDensity unfold(migrMatrix,TUnfold::kHistMapOutputVert);
         if(unfold.SetInput(&pseudoData_MET)>=10000) {
            std::cout<<"Unfolding result may be wrong\n";
         }
         unfold.DoUnfold(0.);
         TH1* unfolded=unfold.GetOutput("Unfolded");
         
         //Plot original genMET distribution
         can.cd();
         genMET_ttbar->SetMarkerSize(0);
         genMET_BSM->SetMarkerSize(0);
         unfolded->SetLineColor(kBlack);
         unfolded->Draw("pe1");
         pseudoData_genMET.SetLineColor(kGray);
         pseudoData_genMET.Draw("hist same");
         genMET_ttbar->Draw("hist e same");
         genMET_BSM->Draw("hist e same");
         saver_templateFit.save(can,"genMET_org/"+cat+"/"+sSample);
         
         //////////////////////////////////////////////////////////
         ///////Template fit of gen level distributions///////
         //////////////////////////////////////////////////////////
         RooRealVar x_gen("MET","MET",0,600);
         RooRealVar nsig_gen("nsig","fitted number of signal events", 0, 5000*scale) ;
         RooRealVar nbkg_gen("nbkg","fitted number of bkg events",0, 1000000);
         
         RooDataHist hist_sig_gen("hist_sig","hist_sig",x_gen,genMET_BSM);
         RooDataHist hist_bkg_gen("hist_bkg","hist_bkg",x_gen,genMET_ttbar);
         RooDataHist hist_data_gen("hist_data","hist_data",x_gen,unfolded);
         
         RooHistPdf sig_pdf_gen("sig_pdf", "" , x_gen, hist_sig_gen);
         RooHistPdf bkg_pdf_gen("bkg_pdf", "" , x_gen, hist_bkg_gen);
         
         RooAddPdf model_gen("model","model",RooArgList(sig_pdf_gen,bkg_pdf_gen),RooArgList(nsig_gen,nbkg_gen));
         
         RooFitResult* fitres_gen=model_gen.fitTo(hist_data_gen,RooFit::SumW2Error(kFALSE),RooFit::Save(), RooFit::PrintLevel(-1));    //Template fit
         
         //RooFit plotting stuff
         RooPlot* xframe_gen = x_gen.frame();
         hist_data_gen.plotOn(xframe_gen,RooFit::Name("data"));
         model_gen.plotOn(xframe_gen,RooFit::Name("fit"));
         model_gen.plotOn(xframe_gen, RooFit::Components(bkg_pdf_gen), RooFit::LineStyle(ELineStyle::kDashed),RooFit::Name("bkg"));
         model_gen.plotOn(xframe_gen, RooFit::Components(sig_pdf_gen), RooFit::LineStyle(ELineStyle::kDashed), RooFit::LineColor(kRed),RooFit::Name("sig"));
         
         TString cat_label=cat;
         if (cat=="emu") cat_label="e#mu";
         else if (cat=="mumu") cat_label="#mu#mu";
         TString labelString=cat_label+"    "+sSample;
         model_gen.paramOn(xframe_gen,RooFit::Label(labelString),RooFit::Layout(0.6,0.95,0.75));
         
         xframe_gen->SetTitle("");
         xframe_gen->SetXTitle("genMET (GeV)");
         xframe_gen->SetYTitle("Events / Bin");
         xframe_gen->SetMinimum(0.1);
         xframe_gen->Draw();
         
         TLegend leg1_gen(0.65,0.75,0.86,0.9);
         leg1_gen.AddEntry("bkg","Background","L");
         leg1_gen.AddEntry("sig","Signal","L");
         leg1_gen.AddEntry("data","Data","P");
         leg1_gen.AddEntry("fit","Fit","L");
         leg1_gen.SetBorderSize(0);
         leg1_gen.Draw();
         
         saver_templateFit.save(can,"genMET_unf/"+cat+"/"+sSample);
         
         fitres_gen->Print();
         can.Clear();
         
         //////////////////////////////////////////////////////////
         ///////Template fit of detector level distributions///////
         //////////////////////////////////////////////////////////
         RooRealVar x("MET","MET",0,600);
         RooRealVar nsig("nsig","fitted number of signal events", 0, 5000*scale) ;
         RooRealVar nbkg("nbkg","fitted number of bkg events",0, 1000000);
         
         //~ MET_ttbar->Rebin(5);
         //~ MET_BSM->Rebin(5);
         //~ pseudoData_MET.Rebin(5);
         
         RooDataHist hist_sig("hist_sig","hist_sig",x,MET_BSM);
         RooDataHist hist_bkg("hist_bkg","hist_bkg",x,MET_ttbar);
         RooDataHist hist_data("hist_data","hist_data",x,&pseudoData_MET);
         
         RooHistPdf sig_pdf("sig_pdf", "" , x, hist_sig);
         RooHistPdf bkg_pdf("bkg_pdf", "" , x, hist_bkg);
         
         RooAddPdf model("model","model",RooArgList(sig_pdf,bkg_pdf),RooArgList(nsig,nbkg));
         
         RooFitResult* fitres=model.fitTo(hist_data,RooFit::SumW2Error(kFALSE),RooFit::Save(), RooFit::PrintLevel(-1));
         
         RooPlot* xframe = x.frame();
         hist_data.plotOn(xframe,RooFit::Name("data"));
         model.plotOn(xframe,RooFit::Name("fit"));
         std::cout<<xframe->chiSquare()<<std::endl;
         model.plotOn(xframe, RooFit::Components(bkg_pdf), RooFit::LineStyle(ELineStyle::kDashed),RooFit::Name("bkg"));
         model.plotOn(xframe, RooFit::Components(sig_pdf), RooFit::LineStyle(ELineStyle::kDashed), RooFit::LineColor(kRed),RooFit::Name("sig"));
         
         model.paramOn(xframe,RooFit::Label(labelString),RooFit::Layout(0.6,0.95,0.75));
         
         xframe->SetTitle("");
         xframe->SetXTitle("p_{T}^{miss} (GeV)");
         xframe->SetYTitle("Events / Bin");
         xframe->SetMinimum(0.1);
         xframe->Draw();
         
         TLegend leg1(0.65,0.75,0.86,0.9);
         leg1.AddEntry("bkg","Background","L");
         leg1.AddEntry("sig","Signal","L");
         leg1.AddEntry("data","Data","P");
         leg1.AddEntry("fit","Fit","L");
         leg1.SetBorderSize(0);
         leg1.Draw();
         
         saver_templateFit.save(can,"MET/"+cat+"/"+sSample);
         
         fitres->Print();
         can.Clear();
         
         float diff=(nsig.getValV()-nsig_gen.getValV())/sqrt(nsig.getError()*nsig.getError()+nsig_gen.getError()*nsig_gen.getError());
         
         //Output of fitting results
         out<<sSample<<" & "<<cat_label<<" & $"<<nsig.getValV()<<"\\pm"<<nsig.getError()<<"$ & $"<<nsig_gen.getValV()<<"\\pm"<<nsig_gen.getError()<<"$ & $"<<diff<<"$ \\\\"<<std::endl;
         
      }
   }
   
   out.close();
}
