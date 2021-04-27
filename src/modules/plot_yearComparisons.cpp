//Script to plot comparisons between different years

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

Config &cfg=Config::get();

TH1F getRMS(const TH2F* hist2D){
   TProfile::SetDefaultSumw2();
   // ~TProfile TTbar_profile_RMS=*(hist2D->ProfileX("ProfileRMS",1,-1,"s"));
   TProfile TTbar_profile_RMS=*(hist2D->ProfileX("ProfileRMS",2700,3300,"s"));
   std::cout<<"RMS only taken for diff<50!!!!!!!!!!!!!!"<<std::endl;
   // ~TH1F RMS("","",100,0,100);
   TH1F RMS("","",hist2D->GetNbinsX(),0,hist2D->GetXaxis()->GetXmax());
   for (int i=1; i<=TTbar_profile_RMS.GetNbinsX(); i++){
      RMS.SetBinContent(i,TTbar_profile_RMS.GetBinError(i));
      RMS.SetBinError(i,0.0001);
   }
   return RMS;
}

TH1F getAbsoluteValues(const TProfile* hist){
   TH1F hist_abs("","",hist->GetNbinsX(),0,hist->GetXaxis()->GetXmax());
   hist_abs.GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
   for (int i=1; i<=hist_abs.GetNbinsX(); i++){
      hist_abs.SetBinContent(i,abs(hist->GetBinContent(i)));
   }
   return hist_abs;
}

extern "C"
void run()
{
   cfg.setOutput("Comparison");
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_metStudies",false);
   
   TCanvas can;
   TCanvas can_comp;
   can.cd();
   io::RootFileReader histReader18("../2018/histograms_v01.root");
   io::RootFileReader histReader18_oldTune("../2018/histograms_vXX.root");
   io::RootFileReader histReader17("../2017/histograms_v01.root");
   io::RootFileReader histReader17_oldTune("../2017/histograms_vXX.root");
   io::RootFileReader histReader16("../2016/histograms_v24.root");
   
   
   std::vector<TString> comb_names;
   std::vector<TH1F> hist_comb;
   
   //Plotting MET resolution as a function of nInteractions and genMET
   for (std::string plot : {"nVertex_vs_MetRes","Gen_nVertex_vs_MetRes","nVertex_vs_MetResPF","Gen_nVertex_vs_MetResPF","genMET_vs_MetRes","genMET_vs_MetResPF"}){
      bool PF=false;
      bool genMET=false;
      if(plot.find("PF")!=std::string::npos) PF=true;
      if(plot.find("genMET")!=std::string::npos) genMET=true;
      
      can.cd();
      can.Clear();
      TH2F* TTbar_2D_18 = histReader18.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+plot+"/TTbar_diLepton");
      TH2F* TTbar_2D_18_oldTune = histReader18_oldTune.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+plot+"/TTbar_diLepton");
      TH2F* TTbar_2D_17 = histReader17.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+plot+"/TTbar_diLepton");
      TH2F* TTbar_2D_17_oldTune = histReader17_oldTune.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+plot+"/TTbar_diLepton");
      TH2F* TTbar_2D_16 = histReader16.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+plot+"/TTbar_diLepton");
      
      if(genMET){
         TTbar_2D_18->Rebin2D(4,1);
         TTbar_2D_18_oldTune->Rebin2D(4,1);
         TTbar_2D_17->Rebin2D(4,1);
         TTbar_2D_17_oldTune->Rebin2D(4,1);
         TTbar_2D_16->Rebin2D(4,1);
      }
      
      std::cout<<"MEAN only taken for diff<50!!!!!!!!!!!!!!"<<std::endl;
      TH1F Mean_18=getAbsoluteValues((TTbar_2D_18->ProfileX("Profile",2700,3300)));
      TH1F Res_18=getRMS(TTbar_2D_18);
      
      TH1F Mean_17=getAbsoluteValues(TTbar_2D_17->ProfileX("Profile",2700,3300));
      TH1F Res_17=getRMS(TTbar_2D_17);
      
      TH1F Mean_18_oldTune=getAbsoluteValues(TTbar_2D_18_oldTune->ProfileX("Profile",2700,3300));
      TH1F Res_18_oldTune=getRMS(TTbar_2D_18_oldTune);
      
      TH1F Mean_17_oldTune=getAbsoluteValues(TTbar_2D_17_oldTune->ProfileX("Profile",2700,3300));
      TH1F Res_17_oldTune=getRMS(TTbar_2D_17_oldTune);
      
      TH1F Mean_16=getAbsoluteValues(TTbar_2D_16->ProfileX("Profile",2700,3300));
      TH1F Res_16=getRMS(TTbar_2D_16);

      gPad->SetLeftMargin(0.13);
      Mean_18.SetStats(0);
      Mean_18_oldTune.SetStats(0);
      Mean_18.GetYaxis()->SetTitleOffset(1.1);
      Mean_18_oldTune.GetYaxis()->SetTitleOffset(1.1);
      Mean_18.GetYaxis()->SetTitle("p_{T}^{miss}-genMET (GeV)");
      Mean_18_oldTune.GetYaxis()->SetTitle("p_{T}^{miss}-genMET (GeV)");
      if (!genMET) Mean_18.GetXaxis()->SetRangeUser(10,60);
      if (!genMET) Mean_18.GetYaxis()->SetRangeUser(0,50);
      else{
         Mean_18.GetYaxis()->SetRangeUser(0,80);
         Mean_18.GetXaxis()->SetTitle("GenMET (GeV)");
         Mean_18_oldTune.GetYaxis()->SetRangeUser(0,80);
         Mean_18_oldTune.GetXaxis()->SetTitle("GenMET (GeV)");
      }
      Mean_18.SetMarkerSize(0);
      Mean_17.SetMarkerSize(0);
      Mean_18_oldTune.SetMarkerSize(0);
      Mean_17_oldTune.SetMarkerSize(0);
      Mean_16.SetMarkerSize(0);
      Mean_18.SetLineColor(kRed-6);
      Mean_17.SetLineColor(kBlue+2);
      Mean_18_oldTune.SetLineColor(kCyan+2);
      Mean_17_oldTune.SetLineColor(kMagenta+2);
      Mean_16.SetLineColor(kGreen+2);
      Res_18.SetLineColor(kRed-6);
      Res_17.SetLineColor(kBlue+2);
      Res_18_oldTune.SetLineColor(kCyan+2);
      Res_17_oldTune.SetLineColor(kMagenta+2);
      Res_16.SetLineColor(kGreen+2);
      Res_18.SetLineStyle(2);
      Res_17.SetLineStyle(2);
      Res_18_oldTune.SetLineStyle(2);
      Res_17_oldTune.SetLineStyle(2);
      Res_16.SetLineStyle(2);
      Mean_18.Draw("hist");
      Res_18.Draw("same hist");
      Mean_17.Draw("same hist");
      Res_17.Draw("same hist");
      if(!PF){
         Mean_18_oldTune.Draw("same hist");
         Mean_17_oldTune.Draw("same hist");
         Res_18_oldTune.Draw("same hist");
         Res_17_oldTune.Draw("same hist");
      }
      Mean_16.Draw("same hist");
      Res_16.Draw("same hist");
      
      gfx::LegendEntries le;
      le.append(Mean_18,"|Mean|18UL","lep");
      le.append(Res_18,"RMS 18","l");
      le.append(Mean_17,"|Mean|17UL","lep");
      le.append(Res_17,"RMS 17","l");
      le.append(Mean_16,"|Mean|16(v11)","lep");
      le.append(Res_16,"RMS 16","l");
      if(!PF){
         le.append(Mean_18_oldTune,"|Mean|18UL(v11)","lep");
         le.append(Res_18_oldTune,"RMS 18","l");
         le.append(Mean_17_oldTune,"|Mean|17UL(v11)","lep");
         le.append(Res_17_oldTune,"RMS 17","l");
      }
      TLegend leg=le.buildLegend(.4,.7,1-1.5*gPad->GetRightMargin(),-1,2);
      leg.Draw();
      
      TString label_string="PUPPI";
      if(PF) label_string="PF";
      
      TLatex label=gfx::cornerLabel(label_string,1);
      label.Draw();
               
      saver.save(can,plot,true,true);
      
      //Compare PF and Puppi
      if (genMET){
         
         if(PF){
             Mean_18.SetLineColor(kGreen+2);
             Res_18.SetLineColor(kGreen+2);
             Mean_17.SetLineColor(kOrange+2);
             Res_17.SetLineColor(kOrange+2);
         }
         
         hist_comb.push_back(Mean_18);
         hist_comb.push_back(Res_18);
         hist_comb.push_back(Mean_17);
         hist_comb.push_back(Res_17);
         if(!PF){
            hist_comb.push_back(Mean_18_oldTune);
            hist_comb.push_back(Res_18_oldTune);
            hist_comb.push_back(Mean_17_oldTune);
            hist_comb.push_back(Res_17_oldTune);
            comb_names.push_back("Puppi |Mean| 18UL");
            comb_names.push_back("Puppi RMS 18UL");
            comb_names.push_back("Puppi |Mean| 17UL");
            comb_names.push_back("Puppi RMS 17UL");
            comb_names.push_back("PUPPI |Mean| 18UL(v11)");
            comb_names.push_back("PUPPI RMS 18UL(v11)");
            comb_names.push_back("PUPPI |Mean| 17UL(v11)");
            comb_names.push_back("PUPPI RMS 17UL(v11)");
         }
         else{
            comb_names.push_back("PF |Mean| 18UL");
            comb_names.push_back("PF RMS 18UL");
            comb_names.push_back("PF |Mean| 17UL");
            comb_names.push_back("PF RMS 17UL");
         }
      }
      
   }
   
   //Plot PF Puppi Comparison (v15)
   can_comp.cd();
   gfx::LegendEntries le_comb;
   hist_comb[0].Draw("hist");
   le_comb.append(hist_comb[0],(TString)comb_names[0],"le");
   for(int i=1; i<hist_comb.size();i++){
      if(i>3 && i<8) continue;
      hist_comb[i].Draw("same hist");
      le_comb.append(hist_comb[i],(TString)comb_names[i],"le");
   }
   TLegend leg_comb=le_comb.buildLegend(.4,.6,1-1.5*gPad->GetRightMargin(),0.9,2);
   leg_comb.Draw();
   TString label_string="|recoMET-genMET|<50 GeV";
   TLatex label=gfx::cornerLabel(label_string,1);
   label.Draw();
   saver.save(can_comp,"Compare_PF_Puppi_vsGenMET",true,true);
   
   //Plot PF Puppi Comparison (v11)
   can_comp.Clear();
   gfx::LegendEntries le_comb2;
   hist_comb[4].Draw("hist");
   le_comb2.append(hist_comb[4],(TString)comb_names[4],"le");
   for(int i=5; i<hist_comb.size();i++){
      hist_comb[i].Draw("same hist");
      le_comb2.append(hist_comb[i],(TString)comb_names[i],"le");
   }
   TLegend leg_comb2=le_comb2.buildLegend(.4,.6,1-1.5*gPad->GetRightMargin(),0.9,2);
   leg_comb2.Draw();
   label.Draw();
   saver.save(can_comp,"Compare_PF_PuppiV11_vsGenMET",true,true);
   
   
   //Plotting MET resolution
   for (std::string channel : {"ee/","emu/","mumu/",""}){
      for (std::string plot : {"nVertex_vs_MetRes","nVertex_vs_MetResPF"}){
         bool PF=false;
         if(plot.find("PF")!=std::string::npos) PF=true;
         
         can.Clear();
         TH2F* TTbar_2D_18 = histReader18.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+channel+plot+"/TTbar_diLepton");
         TH2F* TTbar_2D_18_oldTune = histReader18_oldTune.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+channel+plot+"/TTbar_diLepton");
         TH2F* TTbar_2D_17 = histReader17.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+channel+plot+"/TTbar_diLepton");
         TH2F* TTbar_2D_17_oldTune = histReader17_oldTune.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+channel+plot+"/TTbar_diLepton");
         TH2F* TTbar_2D_16 = histReader16.read<TH2F>(TString::Format("diff_MET%.1f/",cfg.processFraction*100)+"baseline/"+channel+plot+"/TTbar_diLepton");
         
         TH1D Hist_18=*(TTbar_2D_18->ProjectionY("Projection"));
         TH1D Hist_17=*(TTbar_2D_17->ProjectionY("Projection"));
         TH1D Hist_18_oldTune=*(TTbar_2D_18_oldTune->ProjectionY("Projection"));
         TH1D Hist_17_oldTune=*(TTbar_2D_17_oldTune->ProjectionY("Projection"));
         TH1D Hist_16=*(TTbar_2D_16->ProjectionY("Projection"));
         
         Hist_18.Scale(1./Hist_18.Integral());
         Hist_17.Scale(1./Hist_17.Integral());
         Hist_18_oldTune.Scale(1./Hist_18_oldTune.Integral());
         Hist_17_oldTune.Scale(1./Hist_17_oldTune.Integral());
         Hist_16.Scale(1./Hist_16.Integral());

         gPad->SetLeftMargin(0.15);
         Hist_18.SetStats(0);
         Hist_18.GetXaxis()->SetTitle("p_{T}^{miss}-genMET (GeV)");
         Hist_18.GetYaxis()->SetTitle("Normlized Distribution");
         Hist_18.GetXaxis()->SetRangeUser(-100,100);
         Hist_18.GetYaxis()->SetRangeUser(0,0.005);
         Hist_18.GetXaxis()->SetTitleOffset(0.9);
         Hist_18.SetMarkerSize(0);
         Hist_17.SetMarkerSize(0);
         Hist_18_oldTune.SetMarkerSize(0);
         Hist_17_oldTune.SetMarkerSize(0);
         Hist_16.SetMarkerSize(0);
         Hist_17.SetLineColor(kBlue+2);
         Hist_18_oldTune.SetLineColor(kCyan+2);
         Hist_17_oldTune.SetLineColor(kMagenta+2);
         Hist_16.SetLineColor(kGreen+2);
         Hist_18.Draw("e0");
         Hist_17.Draw("same e0");
         if(!PF){
            Hist_18_oldTune.Draw("same e0");
            Hist_17_oldTune.Draw("same e0");
         }
         Hist_16.Draw("same e0");
         
         gfx::LegendEntries le;
         le.append(Hist_18,TString::Format("18UL(v15) (#mu=%.1f #sigma=%.1f)",Hist_18.GetMean(),Hist_18.GetRMS()),"lep");
         le.append(Hist_17,TString::Format("17UL(v15) (#mu=%.1f #sigma=%.1f)",Hist_17.GetMean(),Hist_17.GetRMS()),"lep");
         le.append(Hist_16,TString::Format("16ReReco(v11) (#mu=%.1f #sigma=%.1f)",Hist_16.GetMean(),Hist_16.GetRMS()),"lep");
         if(!PF){
            le.append(Hist_18_oldTune,TString::Format("18UL(v11) (#mu=%.1f #sigma=%.1f)",Hist_18_oldTune.GetMean(),Hist_18_oldTune.GetRMS()),"lep");
            le.append(Hist_17_oldTune,TString::Format("17UL(v11) (#mu=%.1f #sigma=%.1f)",Hist_17_oldTune.GetMean(),Hist_17_oldTune.GetRMS()),"lep");
         }
         TLegend leg=le.buildLegend(.4,.7,1-1.5*gPad->GetRightMargin(),-1,1);
         leg.Draw();
         
         TString label_string=channel+"PUPPI";
         if(PF) label_string=channel+"PF";
         
         TLatex label=gfx::cornerLabel(label_string,1);
         label.Draw();
         
         TString name="PF_res";
         if(!PF) name="PUPPI_res";
         saver.save(can,channel+name,true,true);
      }
   }
}
