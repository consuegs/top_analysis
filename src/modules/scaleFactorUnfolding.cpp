//Script to study the effect of a scale factors used for the measured ptmiss regarding stability and purity

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
#include <TTreeReader.h>
#include <TFitResult.h>
#include <math.h>

#include <iomanip> 


Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<float> met_bins={0,60,120,230,400};
   std::vector<float> phi_bins={0,0.7,1.4,3.14};
   
   TH2F N_gen;
   TH2F N_rec;
   TH2F N_genrec;
   
   TH2F Evt_gen;  //Histograms with correct weights and not just MC N events
   TH2F Evt_rec;
   
   TH2F eff_gen;  //Histograms to calculate correct efficiency
   TH2F eff_gen_recSom;
   
   TH2F migration;   //Histogram to study the migration
   
   TGraph scaleFactor_Purity_lastbin;
   TGraph scaleFactor_Purity_mean;
   TGraph scaleFactor_Stability_lastbin;
   TGraph scaleFactor_Stability_mean;
   
   int bin_gen;
   int bin_rec;
   
   int nBinsMet=met_bins.size()-1;
   int nBinsPhi=phi_bins.size()-1;
   int nBins=nBinsMet*nBinsPhi;
   
   // ~TString sampleName="";
   // ~TString sampleName="dilepton";
   TString sampleName="MadGraph";
   TFile file("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0.root","read");
   
   for (float scaleFactor=0.5; scaleFactor<=1.04; scaleFactor+=0.01) {
      
      std::cout<<"ScaleFactor="<<scaleFactor<<std::endl;
   
      N_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      N_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      N_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      
      Evt_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      Evt_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      
      eff_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      eff_gen_recSom=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      
      migration=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));         
      
      TTreeReader reader((sampleName=="") ? "ttbar_res100.0/ttbar_res" : "ttbar_res100.0/ttbar_res_"+sampleName, &file);
      
      
      TTreeReaderValue<float> MET   (reader, "MET");
      TTreeReaderValue<float> PtNuNu_true   (reader, "PtNuNu");
      TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
      TTreeReaderValue<float> Phi_rec   (reader, "Phi_rec");
      TTreeReaderValue<float> Phi_gen   (reader, "Phi_NuNu");
      TTreeReaderValue<float> dPhiMETnearJet   (reader, "dPhiMETnearJet");
      TTreeReaderValue<float> dPhiMETfarJet   (reader, "dPhiMETfarJet");
      TTreeReaderValue<float> dPhiMETleadJet   (reader, "dPhiMETleadJet");
      TTreeReaderValue<float> dPhiMETlead2Jet   (reader, "dPhiMETlead2Jet");
      TTreeReaderValue<float> dPhiMETbJet   (reader, "dPhiMETbJet");
      TTreeReaderValue<float> dPhiLep1Lep2   (reader, "dPhiLep1Lep2");
      TTreeReaderValue<float> METsig   (reader, "METsig");
      TTreeReaderValue<float> N   (reader, "N");
      TTreeReaderValue<UInt_t> runNo   (reader, "runNo");
      TTreeReaderValue<UInt_t> lumNo   (reader, "lumNo");
      TTreeReaderValue<ULong64_t> evtNo   (reader, "evtNo");
      TTreeReaderValue<UInt_t> genDecayMode   (reader, "genDecayMode");
      TTreeReaderValue<float> genMET   (reader, "genMET");
      TTreeReaderValue<float> PuppiMET   (reader, "PuppiMET");
      TTreeReaderValue<float> HT_tree   (reader, "HT");
      TTreeReaderValue<UInt_t> n_Interactions(reader, "n_Interactions");
         
      while (reader.Next()){
         
         if(*MET>=met_bins.back()) *MET=met_bins.back()-0.01;     //Handel overflow correctly
         if(*PtNuNu>=met_bins.back()) *PtNuNu=met_bins.back()-0.01;
         if(*Phi_rec>=phi_bins.back()) *Phi_rec=phi_bins.back()-0.01;
         if(*Phi_gen>=phi_bins.back()) *Phi_gen=phi_bins.back()-0.01;
         
         if (*PtNuNu>-1 && *Phi_gen>-1){
            eff_gen.Fill(*PtNuNu,*Phi_gen);
            if (*MET>-1 && *Phi_rec>-1) eff_gen_recSom.Fill(*PtNuNu,*Phi_gen);
         }
         
         if (*MET<met_bins[0] || *PtNuNu<met_bins[0] || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
         
         if(*genDecayMode>3) continue;    //Remove tau events
         
         ///////////////////////////
         /////Ptmiss Scale Factor///
         ///////////////////////////
         *MET=scaleFactor * *MET;
         
         bin_gen=N_gen.Fill(*PtNuNu,*Phi_gen);
         bin_rec=N_rec.Fill(*MET,*Phi_rec);
         
         if(bin_gen==bin_rec) N_genrec.Fill(*PtNuNu,*Phi_gen);
         
         Evt_gen.Fill(*PtNuNu,*Phi_gen,*N);
         Evt_rec.Fill(*MET,*Phi_rec,*N);
         
      }
            
      TH2F stability=N_genrec;
      TH2F purity=N_genrec;
      TH2F efficiency=eff_gen_recSom;
      TH2F efficiency_RG=N_rec;
      TH2F hist_res_phi=N_rec;
      TH2F hist_res_met=N_rec;
      
      stability.Divide(&N_gen);
      purity.Divide(&N_rec);
      efficiency.Divide(&eff_gen);
      efficiency_RG.Divide(&eff_gen);
      
      scaleFactor_Purity_lastbin.SetPoint(scaleFactor_Purity_lastbin.GetN(),scaleFactor,purity.GetBinContent(4,3));
      scaleFactor_Purity_mean.SetPoint(scaleFactor_Purity_lastbin.GetN(),scaleFactor,purity.Integral()/(purity.GetNbinsX()*purity.GetNbinsY()));
      scaleFactor_Stability_lastbin.SetPoint(scaleFactor_Purity_lastbin.GetN(),scaleFactor,stability.GetBinContent(4,3));
      scaleFactor_Stability_mean.SetPoint(scaleFactor_Purity_lastbin.GetN(),scaleFactor,stability.Integral()/(stability.GetNbinsX()*stability.GetNbinsY()));
   }
   file.Close();
   
   io::RootFileSaver saver((sampleName=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sampleName+"%.1f.root",cfg.processFraction*100),"binningUnfolding");
   
   saver.save(scaleFactor_Purity_lastbin,"scaleFactor_Purity_lastbin",true,true);
   saver.save(scaleFactor_Purity_mean,"scaleFactor_Purity_mean",true,true);
   saver.save(scaleFactor_Stability_lastbin,"scaleFactor_Stability_lastbin",true,true);
   saver.save(scaleFactor_Stability_mean,"scaleFactor_Stability_mean",true,true);
   
}
