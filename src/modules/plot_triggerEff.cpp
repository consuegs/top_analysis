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

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_triggerEff");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("triggerEff%.1f",cfg.processFraction*100));
   
   TCanvas can;
   for(TString selection:{"baseline","Zpeak","Zpeak_noJetRequ"}){
      for(TString trigg:{"analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         for(TString channel:{"ee","mumu","emu"}){
            for(TString var:{"pTl1","pTl2","etal1","etal2"}){
               TH1F* base_data=histReader.read<TH1F>(selection+"/referenceTrigg/"+channel+"/"+var+"/MET_split");
               TH1F* data=histReader.read<TH1F>(selection+"/"+trigg+"/"+channel+"/"+var+"/MET_split");
               // ~TH1F* base_MC=histReader.read<TH1F>(selection+"/referenceTrigg/"+channel+"/"+var+"/DrellYan_NLO_M50");
               // ~TH1F* MC=histReader.read<TH1F>(selection+"/analysisTrigg/"+channel+"/"+var+"/DrellYan_NLO_M50");
               TH1F* base_MC=histReader.read<TH1F>(selection+"/referenceTrigg/"+channel+"/"+var+"/TTbar_diLepton");
               TH1F* MC=histReader.read<TH1F>(selection+"/analysisTrigg/"+channel+"/"+var+"/TTbar_diLepton");
               
               float efficiency_data=data->Integral()/float(base_data->Integral());
               float e_u_data= TEfficiency::ClopperPearson(base_data->Integral(),data->Integral(),0.682689492137,true) -efficiency_data;
               float e_d_data=-TEfficiency::ClopperPearson(base_data->Integral(),data->Integral(),0.682689492137,false)+efficiency_data;
               float efficiency_MC=MC->Integral()/float(base_MC->Integral());
               float e_u_MC= TEfficiency::ClopperPearson(base_MC->Integral(),MC->Integral(),0.682689492137,true) -efficiency_MC;
               float e_d_MC=-TEfficiency::ClopperPearson(base_MC->Integral(),MC->Integral(),0.682689492137,false)+efficiency_MC;

               TH1F axis=hist::getRatio(*MC,*base_MC,"Efficiency",hist::ONLY1);
               TEfficiency eff_data_hist(*data,*base_data);
               TEfficiency eff_MC_hist(*MC,*base_MC);
               
               eff_MC_hist.SetLineColor(kRed);
               eff_MC_hist.SetMarkerColor(kRed);
               axis.SetMaximum(1.1);
               axis.SetMinimum(0.75);
               axis.SetStats(0);
               axis.GetYaxis()->SetTitleOffset(0.8);
               
               axis.Draw("axis");
               eff_data_hist.Draw("pe1 same");
               eff_MC_hist.Draw("same");
               
               TLatex label=gfx::cornerLabel(channel+", "+selection,3);
               label.Draw();
               
               gfx::LegendEntries legE;
               legE.append(eff_data_hist,TString::Format("data (eff.=%.2f^{+%.2f}_{-%.2f}%%)",efficiency_data*100,e_u_data*100,e_d_data*100),"lep");
               legE.append(eff_MC_hist,TString::Format("MC (eff.=%.2f^{+%.2f}_{-%.2f}%%)",efficiency_MC*100,e_u_MC*100,e_d_MC*100),"lep");
               TLegend leg=legE.buildLegend(.62,.7,1-(gPad->GetRightMargin()+0.02),-1);
               leg.Draw();
               
               saver.save(can,selection+"/"+trigg+"/"+channel+"/"+var,true,false);
            }
         }
      }
   }
   
   
   // ~gfx::SplitCan spcan;                    
   
   // ~TH1F MetBaseline_MuPt_ratio = TH1F(MetBaseline_MuPt);
   // ~TH1F Signal_MuPt_ratio = TH1F(Signal_MuPt);
   // ~MetBaseline_MuPt.Scale(1,"width");
   // ~Signal_MuPt.Scale(1,"width");
   // ~spcan.cdUp();
   // ~can.SetRightMargin (.14);
   // ~can.SetBottomMargin(.14);
   // ~can.SetLogy();
   // ~MetBaseline_MuPt.SetMarkerSize(0);
   // ~MetBaseline_MuPt.Draw("hist");
   // ~MetBaseline_MuPt.SetStats(0);
   // ~Signal_MuPt.SetMarkerSize(0);
   // ~Signal_MuPt.SetLineWidth(0);
   // ~Signal_MuPt.SetLineColor(kWhite);
   // ~Signal_MuPt.SetFillColor(kRed);
   // ~Signal_MuPt.SetFillStyle(3002);
   // ~Signal_MuPt.Draw("hist same b");
   // ~Signal_MuPt.SetStats(0);;
   // ~gPad->SetLogy(true);
   // ~gfx::LegendEntries legE;
   // ~legE.append(MetBaseline_MuPt,"HT-baseline","l");
   // ~legE.append(Signal_MuPt,"HT-baseline\n & signal","f");
   // ~TLegend leg=legE.buildLegend(.6,.75);
   // ~leg.SetTextSize(0.028);
   // ~leg.Draw();


   // ~spcan.cdLow();
   // ~TH1F histRatioFrame=hist::getRatio(MetBaseline_MuPt,MetBaseline_MuPt,"Effizienz",hist::ONLY1);
   // ~histRatioFrame.Draw("axis");
   // ~histRatioFrame.Draw("same axig");
   // ~histRatioFrame.SetMinimum(0);
   // ~histRatioFrame.SetMaximum(1.1);

   // ~TEfficiency eff(Signal_MuPt_ratio,MetBaseline_MuPt_ratio);

   // ~TH1 *htot =eff.GetCopyTotalHisto();
   // ~TH1 *hpass=eff.GetCopyPassedHisto();
   // ~int lowerBin=0;
   // ~if(name.Contains("Mu50")){lowerBin=htot->FindBin(56);}
   // ~if(name.Contains("isoMu22")){lowerBin=htot->FindBin(26);}
   // ~if(name.Contains("Mu45eta2p1")){lowerBin=htot->FindBin(51);}
   // ~int upperBin=htot->FindBin(490);
   // ~int tot=htot ->Integral(lowerBin,upperBin);
   // ~int pas=hpass->Integral(lowerBin,upperBin);
   // ~float efficiency=pas/float(tot);
   // ~float e_u= TEfficiency::ClopperPearson(tot,pas,0.682689492137,true) -efficiency;
   // ~float e_d=-TEfficiency::ClopperPearson(tot,pas,0.682689492137,false)+efficiency;
   // ~float lowEdge=htot->GetBinLowEdge(lowerBin);
   // ~float uppEdge=htot->GetBinLowEdge(upperBin+1);
   // ~TLatex effLabel(.5*(lowEdge+uppEdge),.05,TString::Format("%.1f^{#plus%.1f}_{#minus%.1f}%% ",efficiency*100,e_u*100,e_d*100));
   // ~effLabel.SetNDC(false);
   // ~effLabel.SetTextFont(42);
   // ~effLabel.SetTextSize(0.04);
   // ~effLabel.SetTextAlign(21);
   // ~effLabel.DrawClone();
   // ~TBox box;
   // ~box.SetFillColor(kRed);
   // ~box.SetFillStyle(1001);
   // ~box.DrawBox(lowEdge,efficiency-e_d,uppEdge,efficiency+e_u);
   // ~TLine line;
   // ~line.SetLineWidth(2);
   // ~line.SetLineStyle(2);
   // ~line.SetLineColor(kRed+2);
   // ~line.DrawLine(lowEdge,efficiency,uppEdge,efficiency);

   // ~eff.Draw("same E1");
   // ~eff.SetLineColor(kBlack);
   // ~eff.SetMarkerSize(0.5);

   // ~saver.save(spcan,"MuPt_Signal_MetBaseline_Eff",false);
}
