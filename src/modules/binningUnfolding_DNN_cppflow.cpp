//Script to test  stability and purity with DNN application using cppflow

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TTreeReader.h>
#include <TFitResult.h>
#include <math.h>
#include "cppflow/ops.h"
#include "cppflow/model.h"

#include <iomanip> 


Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<float> met_bins={0,40,80,120,160,230,400};
   std::vector<float> phi_bins={0,0.7,1.4,3.14};
   
   TH2F N_gen;
   TH2F N_rec;
   TH2F N_genrec;
   
   TH2F Evt_gen;  //Histograms with correct weights and not just MC N events
   TH2F Evt_rec;
   TH2F Evt_genrec;
   
   TH2F eff_gen;  //Histograms to calculate correct efficiency
   TH2F eff_gen_recSom;
   
   TH2F migration;   //Histogram to study the migration
      
   int bin_gen;
   int bin_rec;
   
   float PuppiMet_org;
   float PuppiMetcorr_org;
   float PuppiMetscaled_org;
   float diffMET;
   float diffPHI;
   
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
   
   std::vector<TH1F> genmet_hists;    //Hist to evaluate DNN performance
   std::vector<TH1F> puppimet_hists;
   std::vector<TH1F> pfmet_hists;
   std::vector<TH1F> ptnunu_hists;
   std::vector<TH1F> puppimetCorr_hists;
   std::vector<TH1F> puppimetScaled_hists;
   std::vector<TH1F> DNNoutput_hists;
   std::vector<TH1F> diffPuppi_hists;
   std::vector<TH1F> diffPF_hists;
   std::vector<TH1F> diffPuppiCorr_hists;
   std::vector<TH1F> diffPuppiScaled_hists;
   std::vector<TH1F> diffPuppi_ptnunu_hists;
   std::vector<TH1F> diffPF_ptnunu_hists;
   std::vector<TH1F> diffPuppiCorr_ptnunu_hists;
   std::vector<TH1F> diffPuppiScaled_ptnunu_hists;
   for(int i=0;i<nBins;i++){
      TH1F temp_met("","",500,0,500);
      // ~TH1F temp_dnn("","",100,0,5);
      TH1F temp_dnn("","",500,-100,200);
      TH1F temp_res("","",200,-150,150);
      genmet_hists.push_back(temp_met);
      ptnunu_hists.push_back(temp_met);
      puppimet_hists.push_back(temp_met);
      pfmet_hists.push_back(temp_met);
      puppimetCorr_hists.push_back(temp_met);
      puppimetScaled_hists.push_back(temp_met);
      DNNoutput_hists.push_back(temp_dnn);
      diffPuppi_hists.push_back(temp_res);
      diffPF_hists.push_back(temp_res);
      diffPuppiCorr_hists.push_back(temp_res);
      diffPuppiScaled_hists.push_back(temp_res);
      diffPuppi_ptnunu_hists.push_back(temp_res);
      diffPF_ptnunu_hists.push_back(temp_res);
      diffPuppiCorr_ptnunu_hists.push_back(temp_res);
      diffPuppiScaled_ptnunu_hists.push_back(temp_res);
   }
   
   N_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   N_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   N_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   
   Evt_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   Evt_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   Evt_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   
   eff_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   eff_gen_recSom=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
   
   migration=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
      
   // ~TString sampleName="";
   // ~TString sampleName="diLepton";
   // ~TString sampleName="dilepton_CP5";
   // ~TString sampleName="MadGraph";
   TString sampleName="T2tt_650_350";
   TString treeName="TTbar_"+sampleName;
   if (sampleName=="T2tt_650_350") treeName=sampleName;
   TFile file(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/100.0/"+treeName+".root",cfg.year.Data(),cfg.treeVersion.Data()),"read");
   // ~TFile file(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/1.0/"+treeName+".root",cfg.year.Data(),cfg.treeVersion.Data()),"read");
   TTreeReader reader((sampleName=="") ? "ttbar_res100.0/ttbar_res" : "ttbar_res100.0/"+treeName, &file);
   // ~TTreeReader reader((sampleName=="") ? "ttbar_res1.0/ttbar_res" : "ttbar_res1.0/"+treeName, &file);
   
   
   // ~TTreeReaderValue<float> MET   (reader, "MET");
   TTreeReaderValue<float> MET   (reader, "PuppiMET");
   // ~TTreeReaderValue<float> MET   (reader, "XYcorrMET");
   TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
   // ~TTreeReaderValue<float> PtNuNu   (reader, "genMET");
   // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_rec");
   TTreeReaderValue<float> Phi_rec   (reader, "Phi_recPuppi");
   // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_recXYcorr");
   TTreeReaderValue<float> Phi_gen   (reader, "Phi_NuNu");
   // ~TTreeReaderValue<float> Phi_gen   (reader, "Phi_gen");
   // ~TTreeReaderValue<float> dPhiMETnearJet   (reader, "dPhiMETnearJet");
   // ~TTreeReaderValue<float> dPhiMETfarJet   (reader, "dPhiMETfarJet");
   // ~TTreeReaderValue<float> dPhiMETleadJet   (reader, "dPhiMETleadJet");
   // ~TTreeReaderValue<float> dPhiMETlead2Jet   (reader, "dPhiMETlead2Jet");
   // ~TTreeReaderValue<float> dPhiMETbJet   (reader, "dPhiMETbJet");
   // ~TTreeReaderValue<float> dPhiLep1Lep2   (reader, "dPhiLep1Lep2");
   // ~TTreeReaderValue<float> METsig   (reader, "METsig");
   TTreeReaderValue<float> N   (reader, "N");
   TTreeReaderValue<UInt_t> runNo   (reader, "runNo");
   TTreeReaderValue<UInt_t> lumNo   (reader, "lumNo");
   TTreeReaderValue<ULong64_t> evtNo   (reader, "evtNo");
   TTreeReaderValue<UInt_t> genDecayMode   (reader, "genDecayMode");
   TTreeReaderValue<float> genMET   (reader, "genMET");
   TTreeReaderValue<float> PuppiMET   (reader, "PuppiMET");
   TTreeReaderValue<float> PFMET   (reader, "MET");
   TTreeReaderValue<float> Phi_recPuppi   (reader, "Phi_recPuppi");
   // ~TTreeReaderValue<float> PuppiMET   (reader, "MET");
   TTreeReaderValue<float> HT_tree   (reader, "HT");
   TTreeReaderValue<float> n_Interactions(reader, "n_Interactions");
   TTreeReaderValue<float> dPhiPtnunuMet(reader, "dPhiPtnunuMet");
   TTreeReaderValue<float> leadTop_pT(reader, "leadTop_pT");
   TTreeReaderValue<float> dPhiNuNu(reader, "dPhiNuNu");
   // ~TTreeReaderValue<UInt_t> looseLeptonVeto(reader, "looseLeptonVeto");
   // ~TTreeReaderValue<float> DNNregression(reader, "DNN_regression");
   TTreeReaderValue<UInt_t> emu(reader, "emu");
   
   TTreeReaderValue<float> nJets   (reader, "nJets");
   TTreeReaderValue<float> METunc_Puppi   (reader, "METunc_Puppi");
   TTreeReaderValue<float> METunc_PF   (reader, "METunc_PF");
   // ~TTreeReaderValue<float> PuppiMET_phi   (reader, "PuppiMET_phi");
   // ~TTreeReaderValue<float> PFMET_phi   (reader, "PFMET_phi");
   // ~TTreeReaderValue<float> CaloMET   (reader, "CaloMET");
   // ~TTreeReaderValue<float> CaloMET_phi   (reader, "CaloMET_phi");
   // ~TTreeReaderValue<float> MHT   (reader, "MHT");
   TTreeReaderValue<float> Lep1_pt   (reader, "Lep1_pt");
   TTreeReaderValue<float> Lep1_phi   (reader, "Lep1_phi");
   TTreeReaderValue<float> Lep1_eta   (reader, "Lep1_eta");
   TTreeReaderValue<float> Lep1_E   (reader, "Lep1_E");
   TTreeReaderValue<float> Lep1_flavor   (reader, "Lep1_flavor");
   TTreeReaderValue<float> Lep2_pt   (reader, "Lep2_pt");
   TTreeReaderValue<float> Lep2_phi   (reader, "Lep2_phi");
   TTreeReaderValue<float> Lep2_eta   (reader, "Lep2_eta");
   TTreeReaderValue<float> Lep2_E   (reader, "Lep2_E");
   TTreeReaderValue<float> Lep2_flavor   (reader, "Lep2_flavor");
   TTreeReaderValue<float> Jet1_pt   (reader, "Jet1_pt");
   TTreeReaderValue<float> Jet1_phi   (reader, "Jet1_phi");
   TTreeReaderValue<float> Jet1_eta   (reader, "Jet1_eta");
   TTreeReaderValue<float> Jet1_E   (reader, "Jet1_E");
   TTreeReaderValue<float> Jet1_unc   (reader, "Jet1_unc");
   TTreeReaderValue<float> Jet1_bTagScore   (reader, "Jet1_bTagScore");
   TTreeReaderValue<float> Jet2_pt   (reader, "Jet2_pt");
   TTreeReaderValue<float> Jet2_phi   (reader, "Jet2_phi");
   TTreeReaderValue<float> Jet2_eta   (reader, "Jet2_eta");
   TTreeReaderValue<float> Jet2_E   (reader, "Jet2_E");
   TTreeReaderValue<float> Jet2_unc   (reader, "Jet2_unc");
   TTreeReaderValue<float> Jet2_bTagScore   (reader, "Jet2_bTagScore");
   
   //additional variables
   TTreeReaderValue<float> dPhiMETnearJet   (reader,"dPhiMETnearJet");
   TTreeReaderValue<float> dPhiMETfarJet   (reader,"dPhiMETfarJet");
   TTreeReaderValue<float> dPhiMETleadJet   (reader,"dPhiMETleadJet");
   TTreeReaderValue<float> dPhiMETlead2Jet   (reader,"dPhiMETlead2Jet");
   TTreeReaderValue<float> dPhiMETbJet   (reader,"dPhiMETbJet");
   TTreeReaderValue<float> dPhiLep1Lep2   (reader,"dPhiLep1Lep2");
   TTreeReaderValue<float> dPhiJet1Jet2   (reader,"dPhiJet1Jet2");
   TTreeReaderValue<float> METsig   (reader,"METsig");
   TTreeReaderValue<float> MHT   (reader,"MHT");
   TTreeReaderValue<float> MT   (reader,"MT");
   TTreeReaderValue<UInt_t> looseLeptonVeto   (reader,"looseLeptonVeto");
   TTreeReaderValue<float> dPhiMETnearJet_Puppi   (reader,"dPhiMETnearJet_Puppi");
   TTreeReaderValue<float> dPhiMETfarJet_Puppi   (reader,"dPhiMETfarJet_Puppi");
   TTreeReaderValue<float> dPhiMETleadJet_Puppi   (reader,"dPhiMETleadJet_Puppi");
   TTreeReaderValue<float> dPhiMETlead2Jet_Puppi   (reader,"dPhiMETlead2Jet_Puppi");
   TTreeReaderValue<float> dPhiMETbJet_Puppi   (reader,"dPhiMETbJet_Puppi");
   TTreeReaderValue<float> dPhiLep1bJet   (reader,"dPhiLep1bJet");
   TTreeReaderValue<float> dPhiLep1Jet1   (reader,"dPhiLep1Jet1");
   TTreeReaderValue<float> mLL   (reader,"mLL");
   TTreeReaderValue<float> PFMET_phi   (reader,"PFMET_phi");
   TTreeReaderValue<float> PuppiMET_phi   (reader,"PuppiMET_phi");
   TTreeReaderValue<float> CaloMET   (reader,"CaloMET");
   TTreeReaderValue<float> CaloMET_phi   (reader,"CaloMET_phi");
   TTreeReaderValue<float> MT2   (reader,"MT2");
   TTreeReaderValue<float> vecsum_pT_allJet   (reader,"vecsum_pT_allJet");
   TTreeReaderValue<float> vecsum_pT_l1l2_allJet   (reader,"vecsum_pT_l1l2_allJet");
   TTreeReaderValue<float> mass_l1l2_allJet   (reader,"mass_l1l2_allJet");
   TTreeReaderValue<float> ratio_vecsumpTlep_vecsumpTjet   (reader,"ratio_vecsumpTlep_vecsumpTjet");
   TTreeReaderValue<float> mjj   (reader,"mjj");
   
   //DNN configuration
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo__PuppiMET-genMET_%s_20210428-1103normDistr",cfg.year.Data());     //2016
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo__PuppiMET-genMET_%s_20210427-1222normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo__genMET_%s_20210429-1119normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo_MSE__PuppiMET-genMET_%s_20210429-1529normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo_EP100__PuppiMET-genMET_%s_20210429-1718normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo__PtNuNu_%s_20210429-1832normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/Inlusive_amcatnlo_emu__PuppiMET-genMET_%s_20210503-1224normDistr",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/InlusivePatched_amcatnlo_30EP_sqrt_20Bins__PuppiMET-genMET_%s_20210506-1603genMETweighted",cfg.year.Data());
   // ~std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/InlusivePatched_amcatnlo_30EP_sqrt_305Bins__PuppiMET-genMET_%s_20210506-1546genMETweighted",cfg.year.Data());
   std::string modelName= (std::string) TString::Format("/home/home4/institut_1b/dmeuser/top_analysis/DNN_ttbar/trainedModel_Keras/2016/Inlusive_amcatnlo__PuppiMET-genMET_2016_20210609-1056normDistr",cfg.year.Data());
   // ~std::cout<<"------------------------------only using emu events---------------------------------"<<std::endl;
   cppflow::model model_Inclusive(modelName);
   // ~std::vector<float> input_vec(53);
   std::vector<double> input_vec(53);
   std::vector<int64_t> shape (2);
   shape[0]=1;
   shape[1]=53;
   
   
   int migrated=0;
   int totalEntries=reader.GetEntries(true);
   int iEv=0;
   while (reader.Next()){
      iEv++;
      if (iEv%(std::max(totalEntries/10,1))==0){
         io::log*".";
         io::log.flush();
      }
      
      // ~if (iEv>100000) break;
      // ~if(*emu!=1) continue;      //only selected emu events
      // ~if (*MET<120 || *MET>160) continue;
      
      //Set DNN input
      input_vec[0]=*PuppiMET;
      input_vec[1]=*METunc_Puppi;
      input_vec[2]=*PFMET;
      input_vec[3]=*HT_tree;
      input_vec[4]=*nJets;
      input_vec[5]=*n_Interactions;
      input_vec[6]=*Lep1_flavor;
      input_vec[7]=*Lep2_flavor;
      input_vec[8]=*Lep1_pt;
      input_vec[9]=*Lep1_phi;
      input_vec[10]=*Lep1_eta;
      input_vec[11]=*Lep1_E;
      input_vec[12]=*Lep2_pt;
      input_vec[13]=*Lep2_phi;
      input_vec[14]=*Lep2_eta;
      input_vec[15]=*Lep2_E;
      input_vec[16]=*Jet1_pt;
      input_vec[17]=*Jet1_phi;
      input_vec[18]=*Jet1_eta;
      input_vec[19]=*Jet1_E;
      input_vec[20]=*Jet2_pt;
      input_vec[21]=*Jet2_phi;
      input_vec[22]=*Jet2_eta;
      input_vec[23]=*Jet2_E;
      input_vec[24]=*dPhiMETnearJet;
      input_vec[25]=*dPhiMETfarJet;
      input_vec[26]=*dPhiMETleadJet;
      input_vec[27]=*dPhiMETlead2Jet;
      input_vec[28]=*dPhiMETbJet;
      input_vec[29]=*dPhiLep1Lep2;
      input_vec[30]=*dPhiJet1Jet2;
      input_vec[31]=*METsig;
      input_vec[32]=*MHT;
      input_vec[33]=*MT;
      input_vec[34]=*looseLeptonVeto;
      input_vec[35]=*dPhiMETnearJet_Puppi;
      input_vec[36]=*dPhiMETfarJet_Puppi;
      input_vec[37]=*dPhiMETleadJet_Puppi;
      input_vec[38]=*dPhiMETlead2Jet_Puppi;
      input_vec[39]=*dPhiMETbJet_Puppi;
      input_vec[40]=*dPhiLep1bJet;
      input_vec[41]=*dPhiLep1Jet1;
      input_vec[42]=*mLL;
      input_vec[43]=*PFMET_phi;
      input_vec[44]=*PuppiMET_phi;
      input_vec[45]=*CaloMET;
      input_vec[46]=*CaloMET_phi;
      input_vec[47]=*MT2;
      input_vec[48]=*vecsum_pT_allJet;
      input_vec[49]=*vecsum_pT_l1l2_allJet;
      input_vec[50]=*mass_l1l2_allJet;
      input_vec[51]=*ratio_vecsumpTlep_vecsumpTjet;
      input_vec[52]=*mjj;
      auto tensor = cppflow::tensor(input_vec, shape);
      
      PuppiMet_org=*PuppiMET;
      int metBin_org=1;
      float DNN_out=1;
      
      //Target ratio
      // ~if(*MET<0) *MET=0;
      // ~else if(*MET<40) {
         // ~PuppiMetscaled_org=*MET*1.28588;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=1;
      // ~}
      // ~else if(*MET<80){
         // ~PuppiMetscaled_org=*MET*0.94220;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=2;
      // ~}
      // ~else if(*MET<120){
         // ~PuppiMetscaled_org=*MET*0.88487;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=3;
      // ~}
      // ~else if(*MET<160){
         // ~PuppiMetscaled_org=*MET*0.87049;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=4;
      // ~}
      // ~else if(*MET<230){
         // ~PuppiMetscaled_org=*MET*0.88503;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=5;
      // ~}
      // ~else {
         // ~PuppiMetscaled_org=*MET*0.91246;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out*(*MET);
         // ~metBin_org=6;
      // ~}
      
      //target diff
      if(*MET<0) *MET=0;
      else if(*MET<40) {
         PuppiMetscaled_org=*MET*1.28588;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=1;
      }
      else if(*MET<80){
         PuppiMetscaled_org=*MET*0.94220;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=2;
      }
      else if(*MET<120){
         PuppiMetscaled_org=*MET*0.88487;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=3;
      }
      else if(*MET<160){
         PuppiMetscaled_org=*MET*0.87049;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=4;
      }
      else if(*MET<230){
         PuppiMetscaled_org=*MET*0.88503;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=5;
      }
      else {
         PuppiMetscaled_org=*MET*0.91246;
         auto output = model_Inclusive({{"serving_default_batch_normalization_input:0", tensor}},{"StatefulPartitionedCall:0"});
         DNN_out = output[0].get_data<float>()[0];
         *MET = (*MET)-DNN_out;
         metBin_org=6;
      }
      if(*MET<0) {
         std::cout<<"in"<<std::endl;
         *MET=0.0001;
      }
      
      //Target genMET
      // ~if(*MET<0) *MET=0;
      // ~else if(*MET<40) {
         // ~PuppiMetscaled_org=*MET*1.28588;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=1;
      // ~}
      // ~else if(*MET<80){
         // ~PuppiMetscaled_org=*MET*0.94220;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=2;
      // ~}
      // ~else if(*MET<120){
         // ~PuppiMetscaled_org=*MET*0.88487;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=3;
      // ~}
      // ~else if(*MET<160){
         // ~PuppiMetscaled_org=*MET*0.87049;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=4;
      // ~}
      // ~else if(*MET<230){
         // ~PuppiMetscaled_org=*MET*0.88503;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=5;
      // ~}
      // ~else {
         // ~PuppiMetscaled_org=*MET*0.91246;
         // ~auto output = model_Inclusive({{"serving_default_dense_input:0", tensor}},{"StatefulPartitionedCall:0"});
         // ~DNN_out = output[0].get_data<float>()[0];
         // ~*MET = DNN_out;
         // ~metBin_org=6;
      // ~}
      
      PuppiMetcorr_org=*MET;
      
      diffMET=*PtNuNu-*MET;//Save difference before overflow handling
      diffPHI=*Phi_gen-*Phi_rec;
      
      // ~if(*MET>=met_bins.back()) *MET=met_bins.back()-0.01;     //Handel overflow correctly
      // ~if(*PtNuNu>=met_bins.back()) *PtNuNu=met_bins.back()-0.01;
      // ~if(*Phi_rec>=phi_bins.back()) *Phi_rec=phi_bins.back()-0.01;
      // ~if(*Phi_gen>=phi_bins.back()) *Phi_gen=phi_bins.back()-0.01;
      if(*MET>=met_bins.back()) continue;     //Handel overflow correctly
      if(*PtNuNu>=met_bins.back()) continue;
      if(*Phi_rec>=phi_bins.back()) continue;
      if(*Phi_gen>=phi_bins.back()) continue;
      
      if (*PtNuNu>-1 && *Phi_gen>-1){
         eff_gen.Fill(*PtNuNu,*Phi_gen);
         if (*MET>-1 && *Phi_rec>-1) eff_gen_recSom.Fill(*PtNuNu,*Phi_gen);
      }
      
      //Fill hists for DNN performance evaluation with reco events
      if(PuppiMet_org>0){
         genmet_hists[metBin_org-1].Fill(*genMET);
         puppimet_hists[metBin_org-1].Fill(PuppiMet_org);
         pfmet_hists[metBin_org-1].Fill(*PFMET);
         puppimetCorr_hists[metBin_org-1].Fill(PuppiMetcorr_org);
         puppimetScaled_hists[metBin_org-1].Fill(PuppiMetscaled_org);
         DNNoutput_hists[metBin_org-1].Fill(DNN_out);
         diffPuppi_hists[metBin_org-1].Fill(*genMET-PuppiMet_org);
         diffPF_hists[metBin_org-1].Fill(*genMET-*PFMET);
         diffPuppiCorr_hists[metBin_org-1].Fill(*genMET-PuppiMetcorr_org);
         diffPuppiScaled_hists[metBin_org-1].Fill(*genMET-PuppiMetscaled_org);
      }
      
      // ~if(*looseLeptonVeto) continue;   //Remove Events with add. looser lepton
      if(*genDecayMode!=3 && *PtNuNu<40) continue;   //Remove SF events if ptNuNu is smaler than 40GeV
      if (*MET<met_bins[0] || *PtNuNu<met_bins[0] || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
      if(*genDecayMode>3) continue;    //Remove tau events
         
      bin_gen=N_gen.Fill(*PtNuNu,*Phi_gen);
      bin_rec=N_rec.Fill(*MET,*Phi_rec);
      
      if(bin_gen==bin_rec){
         N_genrec.Fill(*PtNuNu,*Phi_gen);
         Evt_genrec.Fill(*PtNuNu,*Phi_gen,*N);
      }
      
      Evt_gen.Fill(*PtNuNu,*Phi_gen,*N);
      Evt_rec.Fill(*MET,*Phi_rec,*N);
      
      int realBin=bin_rec-((bin_rec/(nBinsMet+2)-1)*2+nBinsMet+3);   //Take into account over and underflow bin counting
      int realBin_gen=bin_gen-((bin_gen/(nBinsMet+2)-1)*2+nBinsMet+3);
      
      met_res[realBin].Fill(diffMET);
      phi_res[realBin].Fill(diffPHI);
      
      //Fill hists for DNN performance evaluation with reco_gen events
      // ~genmet_hists[metBin_org-1].Fill(*genMET);
      // ~puppimet_hists[metBin_org-1].Fill(PuppiMet_org);
      // ~pfmet_hists[metBin_org-1].Fill(*PFMET);
      // ~puppimetCorr_hists[metBin_org-1].Fill(PuppiMetcorr_org);
      // ~puppimetScaled_hists[metBin_org-1].Fill(PuppiMetscaled_org);
      // ~DNNoutput_hists[metBin_org-1].Fill(DNN_out);
      // ~diffPuppi_hists[metBin_org-1].Fill(*genMET-PuppiMet_org);
      // ~diffPF_hists[metBin_org-1].Fill(*genMET-*PFMET);
      // ~diffPuppiCorr_hists[metBin_org-1].Fill(*genMET-PuppiMetcorr_org);
      // ~diffPuppiScaled_hists[metBin_org-1].Fill(*genMET-PuppiMetscaled_org);
      ptnunu_hists[metBin_org-1].Fill(*PtNuNu);
      diffPuppi_ptnunu_hists[metBin_org-1].Fill(*PtNuNu-PuppiMet_org);
      diffPF_ptnunu_hists[metBin_org-1].Fill(*PtNuNu-*PFMET);
      diffPuppiCorr_ptnunu_hists[metBin_org-1].Fill(*PtNuNu-PuppiMetcorr_org);
      diffPuppiScaled_ptnunu_hists[metBin_org-1].Fill(*PtNuNu-PuppiMetscaled_org);
   }
   file.Close();
   std::cout<<migrated<<std::endl;
   
   //get DNN name for proper saving
   std::string DNNpath = std::string(modelName);
   std::size_t botDirPos = DNNpath.find_last_of("/");
   std::string DNNname = DNNpath.substr(botDirPos+1, DNNpath.length());
   
   io::RootFileSaver saver(TString::Format("binningUnfolding/"+sampleName+"_"+DNNname+"_%.1f.root",cfg.processFraction*100),"binningUnfolding_DNN");
   
   TH2F stability=Evt_genrec;
   TH2F purity=Evt_genrec;
   TH2F efficiency=eff_gen_recSom;
   TH2F efficiency_RG=N_rec;
   TH2F hist_res_phi=N_rec;
   TH2F hist_res_met=N_rec;
   
   stability.Divide(&Evt_gen);
   purity.Divide(&Evt_rec);
   efficiency.Divide(&eff_gen);
   efficiency_RG.Divide(&eff_gen);
   
   std::vector<TFitResultPtr> r_met;
   std::vector<TFitResultPtr> r_phi;
   
   for(int i=0;i<nBins;i++){
      r_met.push_back(met_res[i].Fit("gaus","SQ"));
      r_phi.push_back(phi_res[i].Fit("gaus","SQ"));
      saver.save(met_res[i],"/Resolution/met_bin"+std::to_string(i+1));
      saver.save(phi_res[i],"/Resolution/phi_bin"+std::to_string(i+1));
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
   
   //Save hists for DNN performance evaluation
   for(int i=0;i<nBinsMet;i++){
      saver.save(genmet_hists[i],"/DNN/GenMET_bin"+std::to_string(i+1));
      saver.save(ptnunu_hists[i],"/DNN/Ptnunu_bin"+std::to_string(i+1));
      saver.save(puppimet_hists[i],"/DNN/PuppiMET_bin"+std::to_string(i+1));
      saver.save(pfmet_hists[i],"/DNN/PFMET_bin"+std::to_string(i+1));
      saver.save(puppimetCorr_hists[i],"/DNN/PuppiMETcorr_bin"+std::to_string(i+1));
      saver.save(puppimetScaled_hists[i],"/DNN/PuppiMETscaled_bin"+std::to_string(i+1));
      saver.save(DNNoutput_hists[i],"/DNN/DNNoutput_bin"+std::to_string(i+1));
      saver.save(diffPuppi_hists[i],"/DNN/diff_bin"+std::to_string(i+1));
      saver.save(diffPF_hists[i],"/DNN/diffPF_bin"+std::to_string(i+1));
      saver.save(diffPuppiCorr_hists[i],"/DNN/diffcorr_bin"+std::to_string(i+1));
      saver.save(diffPuppiScaled_hists[i],"/DNN/diffscaled_bin"+std::to_string(i+1));
      saver.save(diffPuppi_ptnunu_hists[i],"/DNN/diff_ptnunu_bin"+std::to_string(i+1));
      saver.save(diffPF_ptnunu_hists[i],"/DNN/diffPF_ptnunu_bin"+std::to_string(i+1));
      saver.save(diffPuppiCorr_ptnunu_hists[i],"/DNN/diffcorr_ptnunu_bin"+std::to_string(i+1));
      saver.save(diffPuppiScaled_ptnunu_hists[i],"/DNN/diffscaled_ptnunu_bin"+std::to_string(i+1));
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
      
      if (z_axis[i]=="purity" || z_axis[i]=="stability") current_hist.GetZaxis()->SetRangeUser(0.1,0.75);
      
      current_hist.SetStats(false);
      current_hist.SetMarkerColor(kRed);
      current_hist.SetMarkerSize(2.5);
      current_hist.Draw("colz text");
      
      TString plotLoc=z_axis[i];
      saver.save(can,plotLoc,true,true);
      can.Clear();
      i++;
   }
   
}
