//Script to test different binnings regarding stability and purity for unfolding based on minimal ttbarMC tree "~/top_analysis/framework/output/ttbar_res100.0.root"

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
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
   double scaleFactor = 1.0;
   
   std::vector<float> met_bins1={0,100,200,300,400};
   std::vector<float> phi_bins1={0,0.8,1.6,2.4,3.2};
   
   std::vector<float> met_bins2={0,50,100,200,400};
   std::vector<float> phi_bins2={0,0.4,0.8,1.6,3.2};
   
   std::vector<float> met_bins3={0,70,140,250,400};
   std::vector<float> phi_bins3={0,0.4,0.8,1.2,3.14};
   
   std::vector<float> met_bins4={0,60,120,230,400};
   std::vector<float> phi_bins4={0,0.7,1.4,3.14};
   // ~std::vector<float> met_bins4={40,60,90,120,230,400};
   // ~std::vector<float> phi_bins4={0,0.7,1.4,3.14};
   
   TH2F N_gen;
   TH2F N_rec;
   TH2F N_genrec;
   
   TH2F Evt_gen;  //Histograms with correct weights and not just MC N events
   TH2F Evt_rec;
   
   TH2F eff_gen;  //Histograms to calculate correct efficiency
   TH2F eff_gen_recSom;
   
   TH2F migration;   //Histogram to study the migration
   
   int bin_gen;
   int bin_rec;
   
   int numberBinningSchemeMet=1;
   int numberBinningSchemePhi=1;
   
   float diffMET=0.;
   float diffPHI=0.;
   
   float diffMET_normed=0.;
   float diffPHI_normed=0.;
   
   TH1F dPhiMETnearJet_badReso("dphi_metnearJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETfarJet_badReso("dphi_metfarJet"   ,";|#Delta#phi|(p_{T}^{miss},farest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETleadJet_badReso("dphi_metleadJet"   ,";|#Delta#phi|(p_{T}^{miss},lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETlead2Jet_badReso("dphi_metlead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETbJet_badReso("dphi_metbJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiLep1Lep2_badReso("dphi_Lep1Lep2"   ,";|#Delta#phi|(l_1,l_2);EventsBIN"           ,100,0,3.2);
   TH1F METsig_badReso("METsig"   ,";METsig;EventsBIN"           ,100,0,200);
   TH1F diffPuppi_badReso("diffPuppi"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,0.5);
   
   TH1F dPhiMETnearJet_goodReso("dphi_metnearJet_good"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETfarJet_goodReso("dphi_metfarJet_good"   ,";|#Delta#phi|(p_{T}^{miss},farest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETleadJet_goodReso("dphi_metleadJet_good"   ,";|#Delta#phi|(p_{T}^{miss},lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETlead2Jet_goodReso("dphi_metlead2Jet_good"   ,";|#Delta#phi|(p_{T}^{miss},2nd lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETbJet_goodReso("dphi_metbJet_good"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiLep1Lep2_goodReso("dphi_Lep1Lep2_good"   ,";|#Delta#phi|(l_1,l_2);EventsBIN"           ,100,0,3.2);
   TH1F METsig_goodReso("METsig_good"   ,";METsig;EventsBIN"           ,100,0,200);
   TH1F diffPuppi_goodReso("diffPuppi_good"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,0.5);
   
   //Histograms for 2D unfolding (current benchmark)
   std::vector<float>met_bins_unfold={0,30,60,90,120,175,230,315,400};
   std::vector<float>phi_bins_unfold={0,0.35,0.7,1.05,1.4,2.27,3.14};
   TH2F response("response",";binNumber;EventsBIN",48,-0.5,47.5,12,-0.5,11.5);
   TH2F response_sameBins("response_sameBins",";binNumber;EventsBIN",12,-0.5,11.5,12,-0.5,11.5);
   TH1F trueDistributions("trueDistributions",";binNumber;EventsBIN",12,-0.5,11.5);
   TH1F recoDistributions("recoDistributions",";binNumber;EventsBIN",48,-0.5,47.5);
   TH2F reco2D=hist::fromWidths_2d("reco2D",";p_{T}^{miss}(GeV);|#Delta#phi|(p_{T}^{miss},nearest l);",met_bins_unfold,hist::getWidths(met_bins_unfold),phi_bins_unfold,hist::getWidths(phi_bins_unfold));
   
   //Hist for METdiff vs. measured MET
   TH2F METdiff("METdiff_MET",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",1000,0,400,600,-1,3);
   
   
   // ~for(std::vector<float> met_bins : {met_bins1,met_bins2,met_bins3,met_bins4}){    //Test every possible combination of the binning defined above
      // ~for(std::vector<float> phi_bins : {phi_bins1,phi_bins2,phi_bins3,phi_bins4}){
   for(std::vector<float> met_bins : {met_bins4}){    //Test every possible combination of the binning defined above
      for(std::vector<float> phi_bins : {phi_bins4}){
   
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
         
         migration=hist::fromWidths_2d("",";p_{T}^{#nu#nu(+BSM)}(GeV);|#Delta#phi|(p_{T}^{#nu#nu(+BSM)},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         // ~N_gen=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         // ~N_rec=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         // ~N_genrec=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         // ~Evt_gen=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         // ~Evt_rec=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         // ~eff_gen=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         // ~eff_gen_recSom=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         // ~migration=hist::fromWidths_2d("",";genMET(GeV);|#Delta#phi|(genMET,nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         // ~TString sampleName="";
         TString sampleName="dilepton";
         // ~TString sampleName="MadGraph";
         // ~TString sampleName="T2tt_650_350";
         TFile file("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0.root","read");
         TTreeReader reader((sampleName=="") ? "ttbar_res100.0/ttbar_res" : "ttbar_res100.0/ttbar_res_"+sampleName, &file);
         
         
         TTreeReaderValue<float> MET   (reader, "MET");
         // ~TTreeReaderValue<float> MET   (reader, "PuppiMET");
         TTreeReaderValue<float> PtNuNu_true   (reader, "PtNuNu");
         TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
         // ~TTreeReaderValue<float> PtNuNu   (reader, "genMET");
         TTreeReaderValue<float> Phi_rec   (reader, "Phi_rec");
         TTreeReaderValue<float> Phi_gen   (reader, "Phi_NuNu");
         // ~TTreeReaderValue<float> Phi_gen   (reader, "Phi_gen");
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
         // ~TTreeReaderValue<float> PuppiMET   (reader, "MET");
         TTreeReaderValue<float> HT_tree   (reader, "HT");
         TTreeReaderValue<UInt_t> n_Interactions(reader, "n_Interactions");
         
          int migrated=0;
         
         
         
         while (reader.Next()){
            
            diffMET=*PtNuNu-*MET;//Save difference befor overflow handling
            diffPHI=*Phi_gen-*Phi_rec;
            diffMET_normed=diffMET/(*MET);
            diffPHI_normed=diffPHI/(*Phi_rec);
            
            if(*MET>=met_bins.back()) *MET=met_bins.back()-0.01;     //Handel overflow correctly
            if(*PtNuNu>=met_bins.back()) *PtNuNu=met_bins.back()-0.01;
            if(*Phi_rec>=phi_bins.back()) *Phi_rec=phi_bins.back()-0.01;
            if(*Phi_gen>=phi_bins.back()) *Phi_gen=phi_bins.back()-0.01;
            
            if (*PtNuNu>-1 && *Phi_gen>-1){
               eff_gen.Fill(*PtNuNu,*Phi_gen);
               if (*MET>-1 && *Phi_rec>-1) eff_gen_recSom.Fill(*PtNuNu,*Phi_gen);
            }
            
            if(*genDecayMode>3) continue;    //Remove tau events
            
            if (*MET<met_bins[0] || *PtNuNu<met_bins[0] || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
            
            // ~if(*MET<60) *MET=-1;    //Additional cuts, which might solve problem of poor dPhi resolution
            // ~if(*PtNuNu<60) *PtNuNu=-1;
            // ~if(*dPhiMETleadJet<1.0) continue;
            // ~if(*dPhiMETnearJet>0.5) continue;
            // ~if(*METsig<80.) continue;
            
            // ~if(*dPhiMETfarJet>2.64 || *dPhiMETnearJet<0.5) continue;  //Remove Events with met back to back/close with jet
            // ~if(*dPhiMETnearJet<0.5) continue;
            // ~if(*dPhiMETleadJet>2.5 || *dPhiMETlead2Jet>2.5) continue;
            
            // ~if(*dPhiMETlead2Jet<0.8) continue;
            // ~if(*dPhiMETnearJet<0.8) continue;
            // ~if(*dPhiLep1Lep2<2.1) continue;
            // ~if(*dPhiMETleadJet<1.6) continue;
            // ~if(*dPhiMETbJet>2.2) continue;
            
            float HT=*HT_tree;
            // ~if((abs(*MET-*PuppiMET)/HT)>0.05) continue;
            // ~if(abs(*MET-*PuppiMET)>30) continue;
            // ~if(*n_Interactions>35) continue;
            
            ///////////////////////////
            /////Ptmiss Scale Factor///
            ///////////////////////////
            *MET=scaleFactor * *MET;
            
            bin_gen=N_gen.Fill(*PtNuNu,*Phi_gen);
            bin_rec=N_rec.Fill(*MET,*Phi_rec);
            
            if(bin_gen==bin_rec) N_genrec.Fill(*PtNuNu,*Phi_gen);
            
            Evt_gen.Fill(*PtNuNu,*Phi_gen,*N);
            Evt_rec.Fill(*MET,*Phi_rec,*N);
            
            int realBin=bin_rec-((bin_rec/(nBinsMet+2)-1)*2+nBinsMet+3);   //Take into account over and underflow bin counting
            int realBin_gen=bin_gen-((bin_gen/(nBinsMet+2)-1)*2+nBinsMet+3);
            
            met_res[realBin].Fill(diffMET);
            phi_res[realBin].Fill(diffPHI);
            
            METdiff.Fill(diffMET/diffMET_normed,1.+diffMET_normed);  //Same as Ptnunu/MET
            
            // ~if (realBin==11 && realBin_gen!=11 && abs(*MET-*PtNuNu_true)>100) {
            if (realBin==11 && realBin_gen!=11) {
               migration.Fill(*PtNuNu,*Phi_gen);
               // ~std::cout<<"-------------------------------"<<std::endl;
               // ~std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
               // ~std::cout<<"pT_NuNu="<<*PtNuNu<<"   "<<"MET="<<*MET<<"   "<<"genMET="<<*genMET<<"   "<<"puppiMET="<<*PuppiMET<<std::endl;
               // ~std::cout<<"phi_rec="<<*Phi_rec<<"   "<<"phi_gen="<<*Phi_gen<<std::endl;
               // ~std::cout<<*genDecayMode<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETleadJet)<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETlead2Jet)<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETbJet)<<std::endl;
               // ~std::cout<<abs(*dPhiMETnearJet)<<std::endl;
               
               dPhiMETnearJet_badReso.Fill(*dPhiMETnearJet);
               dPhiMETfarJet_badReso.Fill(*dPhiMETfarJet);
               dPhiMETleadJet_badReso.Fill(*dPhiMETleadJet);
               dPhiMETlead2Jet_badReso.Fill(*dPhiMETlead2Jet);
               dPhiMETbJet_badReso.Fill(*dPhiMETbJet);
               dPhiLep1Lep2_badReso.Fill(*dPhiLep1Lep2);
               METsig_badReso.Fill(*METsig);
               diffPuppi_badReso.Fill(abs(*MET-*PuppiMET)/HT);
               migrated++;
            }
            
            else if (realBin==11 && realBin_gen==11) {
               dPhiMETnearJet_goodReso.Fill(*dPhiMETnearJet);
               dPhiMETfarJet_goodReso.Fill(*dPhiMETfarJet);
               dPhiMETleadJet_goodReso.Fill(*dPhiMETleadJet);
               dPhiMETlead2Jet_goodReso.Fill(*dPhiMETlead2Jet);
               dPhiMETbJet_goodReso.Fill(*dPhiMETbJet);
               dPhiLep1Lep2_goodReso.Fill(*dPhiLep1Lep2);
               METsig_goodReso.Fill(*METsig);
               diffPuppi_goodReso.Fill(abs(*MET-*PuppiMET)/HT);
            }
            
            //Fill histograms for unfolding
            int bin_rec_unfold=reco2D.Fill(*MET,*Phi_rec);
            int realBin_rec_unfold=bin_rec_unfold-((bin_rec_unfold/(2*nBinsMet+2)-1)*2+2*nBinsMet+3);
            trueDistributions.Fill(realBin_gen,*N);
            recoDistributions.Fill(realBin_rec_unfold,*N);
            response.Fill(realBin_rec_unfold,realBin_gen,*N);
            response_sameBins.Fill(realBin,realBin_gen);
            
            // ~std::cout<<realBin_gen<<"   "<<*PtNuNu<<"   "<<*Phi_gen<<std::endl;
            // ~std::cout<<realBin_rec_unfold<<"   "<<*MET<<"   "<<*Phi_rec<<std::endl;
            // ~std::cout<<"--------------------------------------"<<std::endl;
            
            
         }
         file.Close();
         std::cout<<migrated<<std::endl;
         
         // ~io::RootFileSaver saver(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100),"binningUnfolding");
         // ~io::RootFileSaver saver(TString::Format("binningUnfolding_dilepton%.1f.root",cfg.processFraction*100),"binningUnfolding");
         io::RootFileSaver saver((sampleName=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sampleName+"%.1f.root",cfg.processFraction*100),"binningUnfolding");
         
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
         
         
         TString Binning= (scaleFactor==1.0) ? "Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)
                                                : "Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+"_SF"+std::to_string(scaleFactor);
         for(int i=0;i<nBins;i++){
            r_met.push_back(met_res[i].Fit("gaus","SQ"));
            r_phi.push_back(phi_res[i].Fit("gaus","SQ"));
            saver.save(met_res[i],Binning+"/Resolution/met_bin"+std::to_string(i+1));
            saver.save(phi_res[i],Binning+"/Resolution/phi_bin"+std::to_string(i+1));
         }
         
         if (!r_met.empty() && !r_phi.empty()){
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
         }
         
         std::vector<TString> z_axis={"stability","purity","efficiency","efficiency_RG","N_gen","N_rec","N_gen_rec","Evt_gen","Evt_rec","res_phi","res_met","migration"};
         int i=0;
         
         for(TH2F current_hist: {stability,purity,efficiency,efficiency_RG,N_gen,N_rec,N_genrec,Evt_gen,Evt_rec,hist_res_phi,hist_res_met,migration}){
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
            current_hist.GetZaxis()->SetLabelOffset(0.02);
            
            if (z_axis[i]=="Evt_rec") current_hist.GetZaxis()->SetTitle("expected events");
            else if (z_axis[i]=="res_phi") current_hist.GetZaxis()->SetTitle("#Delta#phi resolution");
            else if (z_axis[i]=="res_met") current_hist.GetZaxis()->SetTitle("p_{T}^{miss} resolution (GeV)");
            else current_hist.GetZaxis()->SetTitle(z_axis[i]);
            
            current_hist.SetStats(false);
            current_hist.SetMarkerColor(kRed);
            current_hist.SetMarkerSize(2.5);
            current_hist.Draw("colz text");
            
            TString plotLoc=Binning+"/"+z_axis[i];
            saver.save(can,plotLoc,true,true);
            can.Clear();
            i++;
         }
         
         std::vector<TH1F> goodHists={dPhiMETnearJet_goodReso,dPhiMETfarJet_goodReso,dPhiMETleadJet_goodReso,dPhiMETlead2Jet_goodReso,dPhiMETbJet_goodReso,dPhiLep1Lep2_goodReso,METsig_goodReso,diffPuppi_goodReso};
         std::vector<TH1F> badHists={dPhiMETnearJet_badReso,dPhiMETfarJet_badReso,dPhiMETleadJet_badReso,dPhiMETlead2Jet_badReso,dPhiMETbJet_badReso,dPhiLep1Lep2_badReso,METsig_goodReso,diffPuppi_badReso};
         std::vector<TString> names={"dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","METsig","diffPuppi"};
         for (int i=0;i<goodHists.size();i++){
            TCanvas can;
            can.cd();
            
            // ~badHists[i].Rebin(4);
            // ~goodHists[i].Rebin(4);
            badHists[i].Rebin(5);
            goodHists[i].Rebin(5);
               
            goodHists[i].Scale(1.0/(goodHists[i].Integral()));   //Normalize the hist to the integral
            badHists[i].Scale(1.0/(badHists[i].Integral()));
            
            can.cd();
            can.SetLogz();
            
            badHists[i].GetYaxis()->SetTitle("normalized distribution");
            goodHists[i].SetLineColor(kBlue);
            badHists[i].SetLineColor(kRed);
            
            badHists[i].SetStats(0);
            goodHists[i].SetStats(0);
            
            badHists[i].Draw("hist");
            goodHists[i].Draw("same hist");
            
            badHists[i].SetMaximum(1.3*std::max(badHists[i].GetMaximum(),goodHists[i].GetMaximum()));
            
            gfx::LegendEntries legE;
            // ~legE.append(*hist_good,"good resolution (<0.3)","l");
            // ~legE.append(*hist_bad,"bad resolution (>0.3)","l");
            legE.append(goodHists[i],"true","l");
            legE.append(badHists[i],"migrated","l");
            TLegend leg=legE.buildLegend(.5,.65,0.75,.9,1);
            leg.SetTextSize(0.035);
            leg.Draw();
            
            TLatex label=gfx::cornerLabel("Last bin",1);
            // ~TLatex label=gfx::cornerLabel("MET>230 GeV, dPhi>1.4",1);
            // ~TLatex label2=gfx::cornerLabel("|#Delta#phi(met,near jet)|>0.5",2);
            label.Draw();
            // ~label2.Draw();
            
            can.RedrawAxis();
            // ~TString plotLoc=sPresel+"/"+sVar;
            TString plotLoc=Binning+"/METresolution/"+names[i];
            saver.save(can,plotLoc,true,true);
            can.Clear();
         }
         
         //Save Graph for METdiff vs met
         TProfile METdiff_profile;
         METdiff_profile=*(METdiff.ProfileX("s"));
         saver.save(METdiff,"METdiff_vs_MET");
         saver.save(METdiff_profile,"METdiff_vs_MET_profile");
         
         //Save histograms for unfolding
         for (int i=0; i<trueDistributions.GetNbinsX(); i++){
            trueDistributions.SetBinError(i,sqrtf(trueDistributions.GetBinContent(i)));
         }
         saver.save(trueDistributions,Binning+"/Unfolding/trueDistribution");
         saver.save(recoDistributions,Binning+"/Unfolding/recoDistribution");
         saver.save(response,Binning+"/Unfolding/response");
         saver.save(response_sameBins,Binning+"/Unfolding/response_sameBins");
         saver.save(reco2D,Binning+"/Unfolding/reco2D");
         
         numberBinningSchemePhi++;
      }
      
      numberBinningSchemeMet++;
      numberBinningSchemePhi=1;
   }
   
}
