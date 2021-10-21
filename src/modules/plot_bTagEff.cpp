//Script to plot trigger efficiency studies

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/systematics.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TEfficiency.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN"};
   
   for (TString systName : systToPlot){
      std::cout<<"---------------"<<systName<<"----------------"<<std::endl;
      Systematic::Systematic currentSystematic(systName);
      TString inputLoc = TString::Format("bTagEff/%s/%s_merged_%s.root",currentSystematic.name().Data(),"TTbar_diLepton",cfg.treeVersion.Data());
      io::RootFileReader histReader(inputLoc,TString::Format("bTagEff%.1f",cfg.processFraction*100));
      io::RootFileSaver saver(TString::Format("bTagEff/plots_bTagEff%.1f.root",cfg.processFraction*100),TString::Format("%s/",currentSystematic.name().Data()));
      io::RootFileSaver histSaver(TString::Format("data/bTagEff/%s.root",currentSystematic.name().Data()),"");
      
      //Plot 2D bTag Eff.
      TCanvas can;
      can.SetLogx();
      for(TString channel:{"ee","mumu","emu"}){
         for(TString flavor:{"B","C","Light"}){
            for(TString tagger:{"DeepJet_loose","DeepCSV_loose"}){
               TH2F* all=histReader.read<TH2F>("baseline/"+channel+"/"+flavor+"_all/TTbar_diLepton");
               TH2F* tagged=histReader.read<TH2F>("baseline/"+channel+"/"+flavor+"_"+tagger+"/TTbar_diLepton");
               
               std::vector<float> Edges_eta={0,0.5,1.0,2.4};
               std::vector<float> Edges_pT={30,50,70,100,150,200,300,1000};
               *all=hist::rebinned(*all,Edges_pT,Edges_eta);
               *tagged=hist::rebinned(*tagged,Edges_pT,Edges_eta);

               TEfficiency eff(*tagged,*all);   // Will complain about weights due to rebinning using only float precision
               eff.SetUseWeightedEvents(kFALSE);
               TH2* eff_hist = eff.CreateHistogram();
               
               for(int i=0; i<=eff_hist->GetNbinsX();i++){
                  for(int j=0; j<=eff_hist->GetNbinsY();j++){
                     eff_hist->SetBinError(i,j,eff.GetEfficiencyErrorUp(eff.GetGlobalBin(i,j)));
                  }
               }
               TString cat;
               if (channel.Contains("ee")) cat="ee";
               else if (channel.Contains("emu")) cat="e#mu";
               else if (channel.Contains("mumu")) cat="#mu#mu";
               TLatex label=gfx::cornerLabel(cat,1);
               
               gPad->SetRightMargin(0.2);
               gPad->SetLeftMargin(0.13);
               
               eff_hist->GetYaxis()->SetTitleOffset(1.0);
               eff_hist->GetZaxis()->SetLabelOffset(0.02);
               if(flavor=="B"){
                  eff_hist->SetMaximum(1.0);
                  eff_hist->SetMinimum(0.6);
               }
               else if(flavor=="C"){
                  eff_hist->SetMaximum(0.7);
                  eff_hist->SetMinimum(0.1);
               }
               else {
                  eff_hist->SetMaximum(0.5);
                  eff_hist->SetMinimum(0.0);
               }
               eff_hist->SetMarkerSize(1.2);
               eff_hist->Draw("colz text e");
               label.Draw();
               saver.save(can,"baseline/"+channel+"/"+flavor+"_"+tagger,true,true,true);
               histSaver.save(*eff_hist,channel+"/"+flavor+"_"+tagger,false,false);
            }
         }
      }
   }
}
