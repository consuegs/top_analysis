//Script to test different binnings regarding stability and purity for unfolding based on minimal ttbarMC tree "~/top_analysis/framework/output/ttbar_res100.0.root"

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


Config const &cfg=Config::get();

extern "C"
void run()
{
   
   std::vector<float> met_bins1={0,100,200,300,400};
   std::vector<float> phi_bins1={0,0.8,1.6,2.4,3.2};
   
   std::vector<float> met_bins2={0,50,100,200,400};
   std::vector<float> phi_bins2={0,0.4,0.8,1.6,3.2};
   
   std::vector<float> met_bins3={0,70,140,250,400};
   std::vector<float> phi_bins3={0,0.4,0.8,1.2,3.14};
   
   std::vector<float> met_bins4={0,80,160,260,400};
   std::vector<float> phi_bins4={0,0.7,1.4,3.14};
   
   TH2F N_gen;
   TH2F N_rec;
   TH2F N_genrec;
   
   TH2F Evt_gen;  //Histograms with correct weights and not just MC N events
   TH2F Evt_rec;
   
   TH2F eff_gen;  //Histograms to calculate correct efficiency
   TH2F eff_gen_recSom;
   
   int bin_gen;
   int bin_rec;
   
   int numberBinningSchemeMet=1;
   int numberBinningSchemePhi=1;
   
   for(std::vector<float> met_bins : {met_bins1,met_bins2,met_bins3,met_bins4}){    //Test every possible combination of the binning defined above
      for(std::vector<float> phi_bins : {phi_bins1,phi_bins2,phi_bins3,phi_bins4}){
   
         std::vector<TH1F> met_res;    //Resolution for each bin
         std::vector<TH1F> phi_res;
         int nBinsMet=met_bins.size()-1;
         int nBinsPhi=phi_bins.size()-1;
         int nBins=nBinsMet*nBinsPhi;
         
         for(int i=0;i<nBins;i++){
            TH1F temp_met("","",200,-100,100);
            TH1F temp_phi("","",200,-3.2,3.2);
            met_res.push_back(temp_met);
            phi_res.push_back(temp_phi);
         }
         
         
         N_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         N_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         N_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         Evt_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         Evt_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         eff_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         eff_gen_recSom=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         TFile file("~/top_analysis/framework/output/ttbar_res100.0.root","read");
         TTreeReader reader("ttbar_res100.0/ttbar_res", &file);
         TTreeReaderValue<float> MET   (reader, "MET");
         TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
         TTreeReaderValue<float> Phi_rec   (reader, "Phi_rec");
         TTreeReaderValue<float> Phi_gen   (reader, "Phi_gen");
         TTreeReaderValue<float> N   (reader, "N");
         
         while (reader.Next()){
            
            if(*MET<60) *MET=-1;    //Additional cuts, which might solve problem of poor dPhi resolution
            if(*PtNuNu<60) *PtNuNu=-1;
            
            if(*MET>=met_bins.back()) *MET=met_bins.back()-0.01;     //Handel overflow correctly
            if(*PtNuNu>=met_bins.back()) *PtNuNu=met_bins.back()-0.01;
            if(*Phi_rec>=phi_bins.back()) *Phi_rec=phi_bins.back()-0.01;
            if(*Phi_gen>=phi_bins.back()) *Phi_gen=phi_bins.back()-0.01;
            
            if (*PtNuNu>-1 && *Phi_gen>-1){
               eff_gen.Fill(*PtNuNu,*Phi_gen);
               if (*MET>-1 && *Phi_rec>-1) eff_gen_recSom.Fill(*PtNuNu,*Phi_gen);
            }
            
            if (*MET<0 || *PtNuNu<0 || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
            
            bin_gen=N_gen.Fill(*PtNuNu,*Phi_gen);
            bin_rec=N_rec.Fill(*MET,*Phi_rec);
            
            if(bin_gen==bin_rec) N_genrec.Fill(*PtNuNu,*Phi_gen);
            
            Evt_gen.Fill(*PtNuNu,*Phi_gen,*N);
            Evt_rec.Fill(*MET,*Phi_rec,*N);
            
            int realBin=bin_rec-((bin_rec/(nBinsMet+2)-1)*2+nBinsMet+3);   //Take into account over and underflow bin counting
            
            met_res[realBin].Fill(*PtNuNu-*MET);
            phi_res[realBin].Fill(*Phi_gen-*Phi_rec);
            
         }
         file.Close();
         
         io::RootFileSaver saver(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100),"binningUnfolding");
         
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
         
         std::vector<TFitResultPtr> r_met;
         std::vector<TFitResultPtr> r_phi;
         
         TString Binning="Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi);
         for(int i=0;i<nBins;i++){
            r_met.push_back(met_res[i].Fit("gaus","SQ"));
            r_phi.push_back(phi_res[i].Fit("gaus","SQ"));
            saver.save(met_res[i],Binning+"/Resolution/met_bin"+std::to_string(i+1));
            saver.save(phi_res[i],Binning+"/Resolution/phi_bin"+std::to_string(i+1));
         }
         
         int i=0;
         for(int y=1;y<=nBinsPhi;y++){
            for(int x=1;x<=nBinsMet;x++){
               std::cout<<"-------------Bin Number "<<i+1<<"----Width_MET="<<hist::getWidths(met_bins)[i%nBinsMet]<<"----Width_PHI="<<hist::getWidths(phi_bins)[i%nBinsPhi]<<"-------------"<<std::endl;
               if(r_met[i]!=-1) std::cout<<"Res_MET="<<r_met[i]->Parameter(2)<<"+-"<<r_met[i]->ParError(2)<<std::endl;
               if(r_phi[i]!=-1) std::cout<<"Res_PHI="<<r_phi[i]->Parameter(2)<<"+-"<<r_phi[i]->ParError(2)<<std::endl;
               std::cout<<"Purity="<<purity.GetBinContent(x,y)<<std::endl;
               std::cout<<"Stability="<<stability.GetBinContent(x,y)<<std::endl;
               std::cout<<"N_Events="<<met_res[i].GetEntries()<<std::endl;
               std::cout<<"N_DataEvents="<<Evt_gen.GetBinContent(x,y)<<std::endl;
               
               if(r_met[i]!=-1) hist_res_phi.SetBinContent(x,y,r_phi[i]->Parameter(2));
               if(r_phi[i]!=-1) hist_res_met.SetBinContent(x,y,r_met[i]->Parameter(2));
               
               i++;
            }
         }
         
         std::vector<TString> z_axis={"stability","purity","efficiency","efficiency_RG","N_gen","N_rec","Evt_gen","Evt_rec","res_phi","res_met"};
         i=0;
         
         for(TH2F current_hist: {stability,purity,efficiency,efficiency_RG,N_gen,N_rec,Evt_gen,Evt_rec,hist_res_phi,hist_res_met}){
            TCanvas can;
            can.cd();
            // ~can.SetLogz();
            
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.11);
            
            current_hist.GetYaxis()->SetTitleOffset(1.3);
            current_hist.GetXaxis()->SetTitleOffset(0.9);
            // ~current_hist.GetZaxis()->SetTitleOffset(1.3);
            current_hist.GetYaxis()->SetTitleSize(0.05);
            current_hist.GetXaxis()->SetTitleSize(0.05);
            current_hist.GetZaxis()->SetTitleSize(0.05);
            current_hist.GetYaxis()->SetLabelSize(0.04);
            current_hist.GetXaxis()->SetLabelSize(0.04);
            current_hist.GetZaxis()->SetLabelSize(0.04);
            current_hist.GetZaxis()->SetTitle(z_axis[i]);
            
            current_hist.SetStats(false);
            current_hist.SetMarkerColor(kRed);
            current_hist.SetMarkerSize(1.5);
            current_hist.Draw("colz text");
            
            TString plotLoc=Binning+"/"+z_axis[i];
            saver.save(can,plotLoc,true,true);
            can.Clear();
            i++;
         }
         
         numberBinningSchemePhi++;
      }
      
      numberBinningSchemeMet++;
      numberBinningSchemePhi=1;
   }
   
}
