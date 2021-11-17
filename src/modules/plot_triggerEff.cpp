//Script to plot trigger efficiency studies

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
#include <TEfficiency.h>

Config const &cfg=Config::get();

// Function to derive and plot eff. for data and MC in 1D
void eff1D(io::RootFileReader const &histReader, io::RootFileSaver const &saver, TString const &dataName, TString const &mcName){
   TCanvas can;
   for(TString selection:{"baselineTrigger"}){
      for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         for(TString channel:{"ee","mumu","emu"}){
            for(TString var:{"pTl1","pTl2","etal1","etal2"}){
               // Read needed histograms
               TH1F* base_data=histReader.read<TH1F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+dataName);
               TH1F* data=histReader.read<TH1F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+dataName);
               TH1F* base_MC=histReader.read<TH1F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+mcName);
               TH1F* MC=histReader.read<TH1F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+mcName);
               
               // Derive integrated efficiency with correct unc.
               float efficiency_data=data->Integral()/float(base_data->Integral());
               float e_u_data= TEfficiency::ClopperPearson(base_data->Integral(),data->Integral(),0.682689492137,true) -efficiency_data;
               float e_d_data=-TEfficiency::ClopperPearson(base_data->Integral(),data->Integral(),0.682689492137,false)+efficiency_data;
               float efficiency_MC=MC->Integral()/float(base_MC->Integral());
               float e_u_MC= TEfficiency::ClopperPearson(base_MC->Integral(),MC->Integral(),0.682689492137,true) -efficiency_MC;
               float e_d_MC=-TEfficiency::ClopperPearson(base_MC->Integral(),MC->Integral(),0.682689492137,false)+efficiency_MC;
               
               // Get binned efficiency
               TH1F axis=hist::getRatio(*MC,*base_MC,"Efficiency",hist::ONLY1);
               TEfficiency eff_data_hist(*data,*base_data);
               TEfficiency eff_MC_hist(*MC,*base_MC);
               
               //Set drawing options
               eff_MC_hist.SetLineColor(kRed);
               eff_MC_hist.SetMarkerColor(kRed);
               axis.SetMaximum(1.1);
               axis.SetMinimum(0.75);
               axis.SetStats(0);
               axis.GetYaxis()->SetTitleOffset(0.8);
               
               // Draw eff. for MC and data
               axis.Draw("axis");
               eff_data_hist.Draw("pe1 same");
               eff_MC_hist.Draw("same");
               
               TLatex label=gfx::cornerLabel(channel+", "+trigg,3);
               label.Draw();
               
               gfx::LegendEntries legE;
               legE.append(eff_data_hist,TString::Format(dataName+" (eff.=%.2f^{+%.2f}_{-%.2f}%%)",efficiency_data*100,e_u_data*100,e_d_data*100),"lep");
               legE.append(eff_MC_hist,TString::Format(mcName+" (eff.=%.2f^{+%.2f}_{-%.2f}%%)",efficiency_MC*100,e_u_MC*100,e_d_MC*100),"lep");
               TLegend leg=legE.buildLegend(.52,.7,1-(gPad->GetRightMargin()+0.02),-1);
               leg.Draw();
               
               saver.save(can,selection+"/"+trigg+"/"+channel+"/"+var,true,false);
            }
         }
      }
   }
}

// Function to derive and plot eff. and SF for data and MC in 2D
void SF2D(io::RootFileReader const &histReader, io::RootFileSaver const &saver, io::RootFileSaver const &saverHist, TString const &dataName, TString const &mcName){
   TCanvas can;
   for(TString selection:{"baselineTrigger"}){
      for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         for(TString channel:{"ee","mumu","emu"}){
            TString var=(channel!="emu")? "pTl1_pTl2":"pTlmu_pTle";
            TH2F* base_data=histReader.read<TH2F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+dataName);
            TH2F* data=histReader.read<TH2F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+dataName);
            TH2F* base_MC=histReader.read<TH2F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+mcName);
            TH2F* MC=histReader.read<TH2F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+mcName);
            
            std::vector<float> Edges={20,40,60,80,100,150,200};
            *base_data=hist::rebinned(*base_data,Edges,Edges);
            *data=hist::rebinned(*data,Edges,Edges);
            *base_MC=hist::rebinned(*base_MC,Edges,Edges);
            *MC=hist::rebinned(*MC,Edges,Edges);

            TEfficiency eff_data(*data,*base_data);
            TH2* eff_data_hist_eUp = eff_data.CreateHistogram();
            TH2* eff_data_hist_eDown = eff_data.CreateHistogram();
            TEfficiency eff_MC(*MC,*base_MC);
            TH2* eff_MC_hist_eUp = eff_MC.CreateHistogram();
            TH2* eff_MC_hist_eDown = eff_MC.CreateHistogram();
            //Set correct error calculating (need since due to hist::rebinned TEfficiency thinks weights are used due to float precision)
            eff_data.SetUseWeightedEvents(false);
            eff_MC.SetUseWeightedEvents(false);
            eff_data.SetStatisticOption(TEfficiency::kFCP);
            eff_MC.SetStatisticOption(TEfficiency::kFCP);
            
            for(int i=0; i<=eff_MC_hist_eUp->GetNbinsX();i++){
               for(int j=0; j<=eff_MC_hist_eUp->GetNbinsY();j++){
                  eff_MC_hist_eUp->SetBinError(i,j,eff_MC.GetEfficiencyErrorUp(eff_MC.GetGlobalBin(i,j)));
                  eff_data_hist_eUp->SetBinError(i,j,eff_data.GetEfficiencyErrorUp(eff_data.GetGlobalBin(i,j)));
                  eff_MC_hist_eDown->SetBinError(i,j,eff_MC.GetEfficiencyErrorLow(eff_MC.GetGlobalBin(i,j)));
                  eff_data_hist_eDown->SetBinError(i,j,eff_data.GetEfficiencyErrorLow(eff_data.GetGlobalBin(i,j)));
               }
            }
            
            TString cat;
            if (channel.Contains("ee")) cat="ee";
            else if (channel.Contains("emu")) cat="e#mu";
            else if (channel.Contains("mumu")) cat="#mu#mu";
            TLatex label=gfx::cornerLabel(cat,1);
            
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            
            // Plot eff. with up and down stat error
            std::vector<TH2*> histVec = {eff_MC_hist_eUp,eff_data_hist_eUp,eff_MC_hist_eDown,eff_data_hist_eDown};
            std::vector<TString> nameVec = {"_MC_statUp","_data_statUp","_MC_statDown","_data_statDown"};
            int i = 0;
            for (auto tempHist : histVec){
               can.Clear();
               tempHist->GetYaxis()->SetTitleOffset(1.0);
               tempHist->GetZaxis()->SetLabelOffset(0.02);
               tempHist->SetMaximum(1.0);
               tempHist->SetMinimum(0.85);
               tempHist->SetMarkerSize(1.2);
               tempHist->Draw("colz text e");
               label.Draw();
               saver.save(can,selection+"/"+trigg+"/"+channel+"/"+var+nameVec[i],true,false);
               i++;
            }
            
            // Derive SF
            eff_data_hist_eUp->Divide(eff_MC_hist_eUp);
            eff_data_hist_eDown->Divide(eff_MC_hist_eDown);
            
            histVec = {eff_data_hist_eUp,eff_data_hist_eDown};
            nameVec = {"_SF_statUp","_SF_statDown"};
            
            // Plot SF with up and down stat error
            i = 0;
            for (auto tempHist : histVec){
               can.Clear();
               saverHist.save(*tempHist,selection+"/"+trigg+"/"+channel+"/"+var+nameVec[i]);      //Save hist for further unc. calculation
               tempHist->SetMaximum(1.05);
               tempHist->SetMinimum(0.9);
               tempHist->GetYaxis()->SetTitleOffset(1.0);
               tempHist->GetZaxis()->SetLabelOffset(0.02);
               tempHist->SetMarkerSize(1.2);
               tempHist->Draw("colz text e");
               label.Draw();
               saver.save(can,selection+"/"+trigg+"/"+channel+"/"+var+nameVec[i],true,false);
               i++;
            }
            
         }
      }
   }
}

// Function to derive correlation between dilepton and baseline triggers
void alpha(io::RootFileReader const &histReader, io::RootFileSaver const &saver, io::RootFileSaver const &saverHist, TString const &mcName){
   TCanvas can;
   for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
      for(TString channel:{"ee","mumu","emu"}){
         TString var=(channel!="emu")? "pTl1_pTl2":"pTlmu_pTle";
         TH2F* base_MC=histReader.read<TH2F>("noTrigger/baselineTrigger/"+channel+"/"+var+"/"+mcName);
         TH2F* all_MC=histReader.read<TH2F>("noTrigger/all/"+channel+"/"+var+"/"+mcName);
         TH2F* diLepton_MC=histReader.read<TH2F>("noTrigger/"+trigg+"/"+channel+"/"+var+"/"+mcName);
         TH2F* both_MC=histReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"/"+mcName);
         
         std::vector<float> Edges={20,40,60,80,100,150,200};
         *all_MC=hist::rebinned(*all_MC,Edges,Edges);
         *base_MC=hist::rebinned(*base_MC,Edges,Edges);
         *diLepton_MC=hist::rebinned(*diLepton_MC,Edges,Edges);
         *both_MC=hist::rebinned(*both_MC,Edges,Edges);

         TEfficiency eff_base(*base_MC,*all_MC);
         TEfficiency eff_diLepton(*diLepton_MC,*all_MC);
         TEfficiency eff_both(*both_MC,*all_MC);
         TH2* eff_base_hist = eff_base.CreateHistogram();
         TH2* eff_diLepton_hist = eff_diLepton.CreateHistogram();
         TH2* eff_both_hist = eff_both.CreateHistogram();
         TH2* alpha = eff_diLepton.CreateHistogram();
         
         TString cat;
         if (channel.Contains("ee")) cat="ee";
         else if (channel.Contains("emu")) cat="e#mu";
         else if (channel.Contains("mumu")) cat="#mu#mu";
         
         
         gPad->SetRightMargin(0.2);
         gPad->SetLeftMargin(0.13);
         
         // Plot eff. used for correlation estimation
         std::vector<TH2*> histVec = {eff_base_hist,eff_diLepton_hist,eff_both_hist};
         std::vector<TString> nameVec = {"_eff_baseline","_eff_trigger","_eff_both"};
         int i = 0;
         for (auto tempHist : histVec){
            can.Clear();
            TLatex label=gfx::cornerLabel(cat+", "+nameVec[i],1);
            tempHist->GetYaxis()->SetTitleOffset(1.0);
            tempHist->GetZaxis()->SetLabelOffset(0.02);
            // ~tempHist->SetMaximum(1.0);
            // ~tempHist->SetMinimum(0.85);
            tempHist->SetMarkerSize(1.2);
            tempHist->Draw("colz text e");
            label.Draw();
            saver.save(can,"alpha/"+trigg+"/"+channel+"/"+var+nameVec[i],true,false);
            i++;
         }
         
         TLatex label=gfx::cornerLabel(cat,1);
         alpha->Multiply(eff_base_hist);
         alpha->Divide(eff_both_hist);
         
         saverHist.save(*alpha,"alpha/"+trigg+"/"+channel+"/"+var+"_alpha");      //Save hist for further unc. calculation
         
         can.Clear();
         alpha->SetMaximum(1.05);
         alpha->SetMinimum(0.9);
         alpha->GetYaxis()->SetTitleOffset(1.0);
         alpha->GetZaxis()->SetLabelOffset(0.02);
         alpha->GetZaxis()->SetTitle("alpha");
         alpha->SetMarkerSize(1.2);
         alpha->Draw("colz text e");
         label.Draw();
         saver.save(can,"alpha/"+trigg+"/"+channel+"/"+var+"_alpha",true,false);
      }
   }
}

void lumiWeightedSF(io::RootFileReader const &histReader, io::RootFileSaver const &saver, io::RootFileSaver const &saverHist, TString const &mcName, TString const &datasetName){
   std::map<int,std::map<TString,float>> lumiFractionsPerEra = {
                                                            {2,{{"B",0.116},{"C",0.233},{"D",0.102},{"E",0.223},{"F",0.326}}},   //2017B,2017C,2017D,2017E,2017F
                                                            {3,{{"A",0.233},{"B",0.118},{"C",0.116},{"D",0.532}}}};        //2018A,2018B,2018C,2018D
   std::map<int,std::vector<TString>> eraNames = {
                                                            {2,{"B","C","D","E","F"}},   //2017B,2017C,2017D,2017E,2017F
                                                            {3,{"A","B","C","D"}}};        //2018A,2018B,2018C,2018D
   
   TCanvas can;
   std::vector<float> Edges={20,40,60,80,100,150,200};
   for(TString selection:{"baselineTrigger"}){
      for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         for(TString channel:{"ee","mumu","emu"}){
            TString var=(channel!="emu")? "pTl1_pTl2":"pTlmu_pTle";
            int i = 0;
            TH2F lumiWeightSF("","",Edges.size(),&Edges[0],Edges.size(),&Edges[0]);
            for (auto dsName: cfg.datasets.getDatasubsetNames({datasetName})){
               TH2F* base_data=histReader.read<TH2F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+dsName);
               TH2F* data=histReader.read<TH2F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+dsName);
               TH2F* base_MC=histReader.read<TH2F>(selection+"/baselineTrigger/"+channel+"/"+var+"/"+mcName);
               TH2F* MC=histReader.read<TH2F>(selection+"/"+trigg+"/"+channel+"/"+var+"/"+mcName);
               
               *base_data=hist::rebinned(*base_data,Edges,Edges);
               *data=hist::rebinned(*data,Edges,Edges);
               *base_MC=hist::rebinned(*base_MC,Edges,Edges);
               *MC=hist::rebinned(*MC,Edges,Edges);

               TEfficiency eff_data(*data,*base_data);
               TH2* eff_data_hist = eff_data.CreateHistogram();
               TEfficiency eff_MC(*MC,*base_MC);
               TH2* eff_MC_hist = eff_MC.CreateHistogram();
               //Set correct error calculating (need since due to hist::rebinned TEfficiency thinks weights are used due to float precision)
               eff_data.SetUseWeightedEvents(false);
               eff_MC.SetUseWeightedEvents(false);
               eff_data.SetStatisticOption(TEfficiency::kFCP);
               eff_MC.SetStatisticOption(TEfficiency::kFCP);
               
               for(int i=0; i<=eff_MC_hist->GetNbinsX();i++){
                  for(int j=0; j<=eff_MC_hist->GetNbinsY();j++){
                     eff_MC_hist->SetBinError(i,j,eff_MC.GetEfficiencyErrorUp(eff_MC.GetGlobalBin(i,j)));
                     eff_data_hist->SetBinError(i,j,eff_data.GetEfficiencyErrorUp(eff_data.GetGlobalBin(i,j)));
                  }
               }
               
               TString cat;
               if (channel.Contains("ee")) cat="ee";
               else if (channel.Contains("emu")) cat="e#mu";
               else if (channel.Contains("mumu")) cat="#mu#mu";
               TLatex label=gfx::cornerLabel(cat+", Era"+eraNames[cfg.year_int][i],1);
               
               gPad->SetRightMargin(0.2);
               gPad->SetLeftMargin(0.13);
               
               //Plot eff
               can.Clear();
               eff_data_hist->GetYaxis()->SetTitleOffset(1.0);
               eff_data_hist->GetZaxis()->SetLabelOffset(0.02);
               eff_data_hist->SetMaximum(1.0);
               eff_data_hist->SetMinimum(0.85);
               eff_data_hist->SetMarkerSize(1.2);
               eff_data_hist->Draw("colz text e");
               label.Draw();
               saver.save(can,selection+"/"+trigg+"/"+channel+"/Era"+eraNames[cfg.year_int][i]+"/"+var+"_data",true,false);

               
               // Derive SF
               eff_data_hist->Divide(eff_MC_hist);
               
               can.Clear();
               eff_data_hist->SetMaximum(1.05);
               eff_data_hist->SetMinimum(0.9);
               eff_data_hist->GetYaxis()->SetTitleOffset(1.0);
               eff_data_hist->GetZaxis()->SetLabelOffset(0.02);
               eff_data_hist->SetMarkerSize(1.2);
               eff_data_hist->Draw("colz text e");
               label.Draw();
               saver.save(can,selection+"/"+trigg+"/"+channel+"/Era"+eraNames[cfg.year_int][i]+"/"+var+"_SF",true,false);
               
               // Store for lumi weighted sum
               eff_data_hist->Scale(lumiFractionsPerEra[cfg.year_int][eraNames[cfg.year_int][i]]);
               lumiWeightSF.Add(eff_data_hist);
               if((i+1)==eraNames[cfg.year_int].size()){  //store if last era is reached
                  can.Clear();
                  saverHist.save(lumiWeightSF,selection+"/"+trigg+"/"+channel+"/"+var+"_SF_lumiWeighted");
                  lumiWeightSF.SetMaximum(1.05);
                  lumiWeightSF.SetMinimum(0.9);
                  lumiWeightSF.GetYaxis()->SetTitleOffset(1.0);
                  lumiWeightSF.GetZaxis()->SetLabelOffset(0.02);
                  lumiWeightSF.SetMarkerSize(1.2);
                  lumiWeightSF.Draw("colz text e");
                  TLatex labelLumi=gfx::cornerLabel(cat+", lumiWeighted",1);
                  labelLumi.Draw();
                  saver.save(can,selection+"/"+trigg+"/"+channel+"/"+var+"_SF_lumiWeighted",true,false);
               }
               i++;
            }
         }
      }
   }
}

void finalSFunc(io::RootFileReader const &sFReader, io::RootFileSaver const &saver, io::RootFileSaver const &saverHist, io::RootFileSaver const &saverSF){
   for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
      for(TString channel:{"ee","mumu","emu"}){
         TString var=(channel!="emu")? "pTl1_pTl2":"pTlmu_pTle";
         TH2F* alpha=sFReader.read<TH2F>("alpha/"+trigg+"/"+channel+"/"+var+"_alpha");
         TH2F* nominalSF_alpha=sFReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_statUp");
         TH2F* nominalSF_lumiWeight=sFReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_statUp");
         TH2F* lumiWeighted=sFReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_lumiWeighted");
         TH2F* totalUp=sFReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_statUp");
         TH2F* totalDown=sFReader.read<TH2F>("baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_statDown");
         
         //Derive unc. based on alpha
         for(int i=0; i<=alpha->GetNbinsX();i++){
            for(int j=0; j<=alpha->GetNbinsY();j++){
               nominalSF_alpha->SetBinError(i,j,abs(alpha->GetBinContent(i,j)-1.));
            }
         }
         saverHist.save(*nominalSF_alpha,"baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_systAlpha");
         if(trigg=="analysisTrigg") saverSF.save(*nominalSF_alpha,channel+"_SF_systAlpha");
         
         //Derive unc. based on lumiWeightedSF
         for(int i=0; i<=lumiWeighted->GetNbinsX();i++){
            for(int j=0; j<=lumiWeighted->GetNbinsY();j++){
               nominalSF_lumiWeight->SetBinError(i,j,abs(lumiWeighted->GetBinContent(i,j)-nominalSF_lumiWeight->GetBinContent(i,j)));
            }
         }
         saverHist.save(*nominalSF_lumiWeight,"baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_systLumiWeight");
         if(trigg=="analysisTrigg") saverSF.save(*nominalSF_lumiWeight,channel+"_SF_systLumiWeight");
         
         //Derive total syst unc.
         for(int i=0; i<=alpha->GetNbinsX();i++){
            for(int j=0; j<=alpha->GetNbinsY();j++){
               nominalSF_alpha->SetBinError(i,j,sqrt(pow(nominalSF_alpha->GetBinError(i,j),2)+pow(nominalSF_lumiWeight->GetBinError(i,j),2)));
            }
         }
         saverHist.save(*nominalSF_alpha,"baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_systTotal");
         if(trigg=="analysisTrigg") saverSF.save(*nominalSF_alpha,channel+"_SF_systTotal");
         
         //Derive total unc.
         for(int i=0; i<=alpha->GetNbinsX();i++){
            for(int j=0; j<=alpha->GetNbinsY();j++){
               totalUp->SetBinError(i,j,sqrt(pow(nominalSF_alpha->GetBinError(i,j),2)+pow(totalUp->GetBinError(i,j),2)));
               totalDown->SetBinError(i,j,sqrt(pow(nominalSF_alpha->GetBinError(i,j),2)+pow(totalDown->GetBinError(i,j),2)));
            }
         }
         
         saverHist.save(*totalUp,"baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_totalUp");
         saverHist.save(*totalDown,"baselineTrigger/"+trigg+"/"+channel+"/"+var+"_SF_totalDown");
         if(trigg=="analysisTrigg"){
             saverSF.save(*totalUp,channel+"_SF_totalUp");
             saverSF.save(*totalDown,channel+"_SF_totalDown");
         }
         
      }
   }
}
         
         
         

extern "C"
void run()
{   
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("triggerEff%.1f",cfg.processFraction*100));
   io::RootFileSaver saver(TString::Format("triggerEff/plots_triggerEff%.1f.root",cfg.processFraction*100),"plot_triggerEff");
   io::RootFileSaver saverHist(TString::Format("triggerEff/hists_triggerEff%.1f.root",cfg.processFraction*100),"triggerEff");
   io::RootFileReader sfReader(TString::Format("triggerEff/hists_triggerEff%.1f.root",cfg.processFraction*100),"triggerEff");
   io::RootFileSaver saverSF(TString::Format("data/TriggerSF_%s.root",cfg.year.Data()),"");

   
   // ~//Plot 1D efficiencies
   // ~eff1D(histReader,saver,"TTbar_diLepton_CP5up","DoubleMuon");
   
   //Plot 2D scale factors (Data/MC)
   SF2D(histReader,saver,saverHist,"TTbar_diLepton_CP5up","DoubleMuon");
   
   // Derive correlation between dilepton and reference trigger
   alpha(histReader,saver,saverHist,"TTbar_diLepton_CP5up");
   
   // Derive SF per era and compare to nominal
   lumiWeightedSF(histReader,saver,saverHist,"TTbar_diLepton_CP5up","DoubleMuon");
   
   // Derive final uncertainties
   finalSFunc(sfReader,saver,saverHist,saverSF);

}
