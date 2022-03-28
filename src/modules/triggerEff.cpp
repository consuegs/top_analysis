#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"
#include "tools/systematics.hpp"
#include "tools/leptonCorrections.hpp"
#include "tools/jetCorrections.hpp"
#include "tools/bTagWeights.hpp"
#include "tools/leptonSF.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>

Config const &cfg=Config::get();

bool matchLepton(TLorentzVector recoLep, TLorentzVector genLep) {
   // ~return (abs(recoLep.DeltaR(genLep))<0.5) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.5); //probably wrong numbers
   return (abs(recoLep.DeltaR(genLep))<0.05) && ((abs(recoLep.Pt()-genLep.Pt())/recoLep.Pt())<0.1);
}

extern "C"
void run()
{
   std::vector<TString> baselineTriggerNames = {
         "HLT_PFHT300_PFMET110_v",
         "HLT_PFMET120_PFMHT120_IDTight_v",
         "HLT_PFMET170_HBHECleaned_v",
         "HLT_PFMET300_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
         "HLT_MET200_v",
         "HLT_MET200_v",   //Dummies to have same length as 2017_2018
         "HLT_MET200_v",
         "HLT_MET200_v",
         "HLT_MET200_v",
   };
   std::vector<TString> baselineTriggerNames2016 = {
         "HLT_PFHT300_PFMET110_v",
         "HLT_PFMET120_PFMHT120_IDTight_v",
         "HLT_PFMET170_HBHECleaned_v",
         "HLT_PFMET300_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
         "HLT_MET200_v",
         "HLT_MET200_v",   //Dummies to have same length as 2017_2018
         "HLT_MET200_v",
         "HLT_MET200_v",
         "HLT_MET200_v",
   };
   std::vector<TString> baselineTriggerNames2017_2018 = {
         "HLT_PFMET200_HBHECleaned_v",
         "HLT_PFMET200_HBHE_BeamHaloCleaned_v",
         "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v",
         "HLT_PFMET120_PFMHT120_IDTight_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
         "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
         "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
         "HLT_PFHT500_PFMET100_PFMHT100_IDTight_v",
         "HLT_PFHT700_PFMET85_PFMHT85_IDTight_v",
         "HLT_PFHT800_PFMET75_PFMHT75_IDTight_v",
   };
         
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   std::vector<std::string> saveDatasubsets = {};
   
   //Read systematic from command line
   Systematic::Systematic currentSystematic(cfg.systematic);
   
   //Define histograms in the following
   hist::Histograms<TH1F> hs(vsDatasubsets);
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   for(TString selection:{"baselineTrigger","noTrigger","nJets2","nJets3+","nPV30-","nPF30+","MET150-","MET150+"}){
      for(TString base:{"all","baselineTrigger","analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
         hs.addHist(selection+"/"+base+"/ee/pTl1"   ,";%pTl1;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/emu/pTle"   ,";%pTle;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/mumu/pTl1"   ,";%pTl1;EventsBIN"           ,10,0,200);
         
         hs.addHist(selection+"/"+base+"/ee/pTl2"   ,";%pTl2;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/emu/pTlmu"   ,";%pTlmu;EventsBIN"           ,10,0,200);
         hs.addHist(selection+"/"+base+"/mumu/pTl2"   ,";%pTl2;EventsBIN"           ,10,0,200);
         
         hs.addHist(selection+"/"+base+"/ee/etal1"   ,";#eta_{l1};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/emu/etale"   ,";#eta_{e};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/mumu/etal1"   ,";#eta_{l1};EventsBIN"           ,15,-2.4,2.4);
         
         hs.addHist(selection+"/"+base+"/ee/etal2"   ,";#eta_{l2};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/emu/etalmu"   ,";#eta_{mu};EventsBIN"           ,15,-2.4,2.4);
         hs.addHist(selection+"/"+base+"/mumu/etal2"   ,";#eta_{l2};EventsBIN"           ,15,-2.4,2.4);
         
         hs2d.addHist(selection+"/"+base+"/ee/pTl1_pTl2"   ,";%pTl1;%pTl2;EventsBIN"           ,40,0,200,40,0,200);
         hs2d.addHist(selection+"/"+base+"/emu/pTlmu_pTle"   ,";%pTlmu;%pTle;EventsBIN"           ,40,0,200,40,0,200);
         hs2d.addHist(selection+"/"+base+"/mumu/pTl1_pTl2"   ,";%pTl1;%pTl2;EventsBIN"           ,40,0,200,40,0,200);
      }
   }
   
   TString dssName_multi="";

   for (TString ds_name: cfg.datasets.getDatasetNames()){
      auto ds=cfg.datasets.getDataset(ds_name);
      if (ds.name.find("MET")==std::string::npos && ds.name.find("TTbar_diLepton")==std::string::npos) {
         std::cout<<ds.name<<" found in running schedule, but module should ony be running with ttbar MC or MET PD"<<std::endl;
         return;
      }
      for (auto dss: cfg.datasets.getDatasubsets({ds.name})){   
         TFile* file = TFile::Open(dss.getPath(),"read");
         if (file->IsZombie()) {
            return;
         }
         
         hs.setCurrentSample(dss.name);
         hs2d.setCurrentSample(dss.name);
         
         bool const isData=dss.isData;
         
         // store also individual eras if data
         if(isData) saveDatasubsets.push_back(dss.name);
         
         // Configure JER Corrections
         jerCorrections jerCorrector = jerCorrections(isData? cfg.jer_SF_data.Data() : cfg.jer_SF_mc.Data(),isData? cfg.jer_RES_data.Data() : cfg.jer_RES_mc.Data(),currentSystematic);
         
         // Configure lepton Correction
         leptonCorrections leptonCorretor = leptonCorrections(currentSystematic);
         
         // Configure leptonSF
         LeptonScaleFactors leptonSF = LeptonScaleFactors(cfg.electronID_file,cfg.electronID_hist,cfg.electronRECO_file,cfg.electronRECO_hist,
                                                         cfg.muonID_file,cfg.muonID_hist,cfg.muonISO_file,cfg.muonISO_hist,currentSystematic);
         
         // Configure bTag Weights
         BTagWeights bTagWeighter = BTagWeights(cfg.bTagSF_file.Data(),cfg.bTagEffPath.Data(),cfg.bTagger.Data(),BTagEntry::OperatingPoint(cfg.bTagWP),cfg.bTagWPcut,currentSystematic);
                  
         //Check if current sample is TTbar 2L sample (later used to veto tau events)
         bool ttBar_dilepton=dss.isTTbar2L;
               
         
         bool Run2016_preVFP=(cfg.year == "2016_preVFP");
         bool Run2016_postVFP=(cfg.year == "2016_postVFP");
         
         //Check if current sample is Run2016H
         bool Run2016H=false;
         if (dss.datasetName.find("Run2016H")!=std::string::npos) Run2016H=true;
         
         
         //Check if current sample is Run2017
         bool Run2017=(cfg.year=="2017");
         //Check if current sample is Run2017AB
         bool Run2017AB=false;
         if (dss.name.find("Run2017A")!=std::string::npos) Run2017AB=true;
         else if (dss.name.find("Run2017B")!=std::string::npos) Run2017AB=true;
         
         //Check if current sample is Run2018
         bool Run2018=(cfg.year=="2018");
         if (dss.name.find("Run2018")!=std::string::npos) Run2018=true;
         
         //Set correct baselineTriggers
         if (Run2016_preVFP || Run2016_postVFP) baselineTriggerNames=baselineTriggerNames2016;
         if (Run2018 || Run2017) baselineTriggerNames=baselineTriggerNames2017_2018;
         
         
         //Check if current sample is DY MC
         bool DY_MC=false;
         if (dss.datasetName.find("DY")!=std::string::npos) DY_MC=true;

         TTreeReader reader(cfg.treeName, file);
         TTreeReaderValue<float> w_pu(reader, "pu_weight");
         TTreeReaderValue<UInt_t> runNo(reader, "runNo");
         TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
         TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
         TTreeReaderValue<float> w_mc(reader, "mc_weight");
         TTreeReaderValue<float> fragCentralWeight(reader, "weightFragCentral");
         TTreeReaderValue<float> w_topPT(reader, "topPTweight");
         TTreeReaderValue<double> w_prefiring(reader, "prefiring_weight");
         TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
         TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
         TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
         TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
         TTreeReaderValue<tree::MET> MET(reader, "met");
         TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
         TTreeReaderValue<float> rho(reader, "rho");
         TTreeReaderValue<int> genDecayMode(reader, "ttbarDecayMode");
         TTreeReaderValue<int> n_Interactions(reader, "nPV");
         
         TTreeReaderValue<bool> muonTrigg1(reader, cfg.muonTrigg1);
         TTreeReaderValue<bool> muonTrigg2(reader, cfg.muonTrigg2);
         TTreeReaderValue<bool> muonTrigg3(reader, cfg.muonTrigg3);
         TTreeReaderValue<bool> muonTrigg4(reader, cfg.muonTrigg4);
         TTreeReaderValue<bool> singleMuonTrigg1(reader, cfg.singleMuonTrigg1);
         TTreeReaderValue<bool> singleMuonTrigg2(reader, cfg.singleMuonTrigg2);
         TTreeReaderValue<bool> eleTrigg1(reader, cfg.eleTrigg1);
         TTreeReaderValue<bool> eleTrigg2(reader, cfg.eleTrigg2);
         TTreeReaderValue<bool> eleMuTrigg1(reader, cfg.eleMuTrigg1);
         TTreeReaderValue<bool> eleMuTrigg2(reader, cfg.eleMuTrigg2);
         TTreeReaderValue<bool> eleMuTrigg3(reader, cfg.eleMuTrigg3);
         TTreeReaderValue<bool> eleMuTrigg4(reader, cfg.eleMuTrigg4);
         TTreeReaderValue<bool> singleEleTrigg(reader, cfg.singleEleTrigg);
         
         TTreeReaderValue<bool> baselineTrigg1(reader, baselineTriggerNames[0]);
         TTreeReaderValue<bool> baselineTrigg2(reader, baselineTriggerNames[1]);
         TTreeReaderValue<bool> baselineTrigg3(reader, baselineTriggerNames[2]);
         TTreeReaderValue<bool> baselineTrigg4(reader, baselineTriggerNames[3]);
         TTreeReaderValue<bool> baselineTrigg5(reader, baselineTriggerNames[4]);
         TTreeReaderValue<bool> baselineTrigg6(reader, baselineTriggerNames[5]);
         TTreeReaderValue<bool> baselineTrigg7(reader, baselineTriggerNames[6]);
         TTreeReaderValue<bool> baselineTrigg8(reader, baselineTriggerNames[7]);
         TTreeReaderValue<bool> baselineTrigg9(reader, baselineTriggerNames[8]);
         TTreeReaderValue<bool> baselineTrigg10(reader, baselineTriggerNames[9]);
         
         io::log * ("Processing '"+dss.name+"' ");
         int iEv=0;
         int processEvents=cfg.processFraction*dss.entries;
         // ~if(DY_MC) processEvents*=0.1;
         while (reader.Next()){
            iEv++;
            if (iEv>processEvents) break;
            if (iEv%(std::max(processEvents/10,1))==0){
               io::log*".";
               io::log.flush();
            }
            
            float met=MET->p.Pt();
            float const met_puppi=MET_Puppi->p.Pt();
            
            //Do not use tau events in signal sample
            if (ttBar_dilepton && *genDecayMode>3) continue;
            
            //Booleans for reco and zpeak selection
            bool rec_selection=false;
            
            // Construct vector of different METs for correction
            std::vector<tree::MET*> PFMETs = {&(*MET)};
            std::vector<tree::MET*> PuppiMETs = {&(*MET_Puppi)};
            
            // Correct and select leptons
            *muons = leptonCorretor.correctMuons(*muons,PFMETs,PuppiMETs);
            *electrons = leptonCorretor.correctElectrons(*electrons,PFMETs,PuppiMETs);
            
            //Apply JER smearing
            if(!isData){
               jerCorrector.smearCollection_Hybrid(*jets,*rho);
            }
            
            //Baseline selection (including separation into ee, emu, mumu)
            TLorentzVector p_l1;
            TLorentzVector p_l2;
            int flavor_l1=0;  //1 for electron and 2 for muon
            int flavor_l2=0;
            double etaSC_l1=0;  // etaSC needed for electron ID SF
            double etaSC_l2=0;
            bool muonLead=true; //Boolean for emu channel
            std::vector<bool> channel={false,false,false};     //ee, mumu, emu
            TString cat="";
            
            rec_selection=selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,etaSC_l1,etaSC_l2,cat,muonLead);
                  
            std::vector<tree::Jet> cjets;
            std::vector<tree::Jet> BJets;
            std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
            if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end(), [](bool v) { return v; })) rec_selection=false;
            
            //Set event weights
            float fEventWeight = 1.;
            float SFWeight = 1.;
            if(!isData){
               int channelID = 0;   // for bTag Weights
               if (channel[0]) channelID = 1;
               else if (channel[1]) channelID = 2;
               else if (channel[2]) channelID = 3;
               
               fEventWeight = *w_pu * *w_mc * *w_topPT;
               float leptonSFweight = leptonSF.getSFDilepton(p_l1,p_l2,flavor_l1,flavor_l2,etaSC_l1,etaSC_l2);
               float bTagWeight = bTagWeighter.getEventWeight(cjets,channelID);
               SFWeight = bTagWeight*leptonSFweight * *w_prefiring;
            }
            hs.setFillWeight(fEventWeight*SFWeight);
            hs2d.setFillWeight(fEventWeight*SFWeight);
            
            //Define trigger selections
            bool baselineTriggers=*baselineTrigg1 || *baselineTrigg2 || *baselineTrigg3 || *baselineTrigg4 || *baselineTrigg5 || *baselineTrigg6 || *baselineTrigg7 || *baselineTrigg8|| *baselineTrigg9 || *baselineTrigg10;
            bool diElectronTriggers=*eleTrigg1 || *eleTrigg2 || *singleEleTrigg;
            bool diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2;
            bool electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
            
            if (isData){
               if(Run2016_preVFP || Run2016_postVFP){
                  if(!Run2016H){ 
                     diMuonTriggers=*muonTrigg3 || *muonTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2; // no DZ
                     electronMuonTriggers=*eleMuTrigg3 || *eleMuTrigg4 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg; // no DZ
                  }else{ 
                     diMuonTriggers=*muonTrigg1 || *muonTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2; // with DZ
                     electronMuonTriggers=*eleMuTrigg1 || *eleMuTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg; // with DZ
                  }
               }
            }
            
            if(Run2017)diMuonTriggers=*muonTrigg2 || *singleMuonTrigg1 || *singleMuonTrigg2;
            if(Run2017AB){
               diMuonTriggers=*muonTrigg1|| *singleMuonTrigg1 || *singleMuonTrigg2;
            }
            
            //Fill hists
            TString path_cat="ee";
            bool trigger=diElectronTriggers;
            bool doubleTrigger=*eleTrigg1 || *eleTrigg2;
            bool doubleTrigger_DZ=*eleTrigg1 || *eleTrigg2;
            bool singleTrigger=*singleEleTrigg;
            if (channel[2]) {
               path_cat="emu";
               trigger=electronMuonTriggers;
               doubleTrigger=*eleMuTrigg3 || *eleMuTrigg4;
               doubleTrigger_DZ=*eleMuTrigg1 || *eleMuTrigg2;
               singleTrigger=*singleMuonTrigg1 || *singleMuonTrigg2 || *singleEleTrigg;
            }
            else if (channel[1]) {
               path_cat="mumu";
               trigger=diMuonTriggers;
               doubleTrigger=*muonTrigg3 || *muonTrigg4;
               doubleTrigger_DZ=*muonTrigg1 || *muonTrigg2;
               singleTrigger=*singleMuonTrigg1 || *singleMuonTrigg2;
            }
            
            std::vector<bool> triggerVec={true,baselineTriggers,trigger,doubleTrigger_DZ,doubleTrigger,singleTrigger};
            
            std::vector<bool> selectionBool;
            if(!(Run2016_preVFP || Run2016_postVFP)){
               selectionBool={baselineTriggers,true,
                                                (cjets.size()<3 && baselineTriggers),
                                                (cjets.size()>=3 && baselineTriggers),
                                                (*n_Interactions<30 && baselineTriggers),
                                                (*n_Interactions>=30 && baselineTriggers),
                                                (met_puppi<150 && baselineTriggers),
                                                (met_puppi>=150 && baselineTriggers)};
            }else{
               selectionBool={baselineTriggers,true,
                                                (cjets.size()<3 && baselineTriggers),
                                                (cjets.size()>=3 && baselineTriggers),
                                                (*n_Interactions<25 && baselineTriggers),
                                                (*n_Interactions>=25 && baselineTriggers),
                                                (met_puppi<150 && baselineTriggers),
                                                (met_puppi>=150 && baselineTriggers)};
            }
            
            int i_trigger=0;
            int i_selection=0;
            // ~for(TString selection:{"baselineTrigger/","noTrigger/"}){
            for(TString selection:{"baselineTrigger/","noTrigger/","nJets2/","nJets3+/","nPV30-/","nPF30+/","MET150-/","MET150+/"}){
               for(TString base:{"all","baselineTrigger","analysisTrigg","doubleTrigg_DZ","doubleTrigg","singleTrigg"}){
                  if(rec_selection==true && triggerVec[i_trigger]){
                     if(selectionBool[i_selection]){
                        if(!channel[2]){
                           hs.fill(selection+base+"/"+path_cat+"/pTl1",p_l1.Pt());
                           hs.fill(selection+base+"/"+path_cat+"/pTl2",p_l2.Pt());
                           hs.fill(selection+base+"/"+path_cat+"/etal1",p_l1.Eta());
                           hs.fill(selection+base+"/"+path_cat+"/etal2",p_l2.Eta());
                           hs2d.fill(selection+base+"/"+path_cat+"/pTl1_pTl2",p_l1.Pt(),p_l2.Pt());
                        }
                        else{
                           if(muonLead){
                              hs2d.fill(selection+base+"/"+path_cat+"/pTlmu_pTle",p_l1.Pt(),p_l2.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/pTle",p_l2.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/pTlmu",p_l1.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/etale",p_l2.Eta());
                              hs.fill(selection+base+"/"+path_cat+"/etalmu",p_l1.Eta());
                           }
                           else{
                              hs2d.fill(selection+base+"/"+path_cat+"/pTlmu_pTle",p_l2.Pt(),p_l1.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/pTle",p_l1.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/pTlmu",p_l2.Pt());
                              hs.fill(selection+base+"/"+path_cat+"/etale",p_l1.Eta());
                              hs.fill(selection+base+"/"+path_cat+"/etalmu",p_l2.Eta());
                           }
                        }
                     }
                  }
                  i_trigger++;
               }
               i_trigger=0;
               i_selection++;
            }
            
         }// evt loop
         io::log<<"";
         
         // ~hs.scaleLumi();
         // ~hs2d.scaleLumi();
         hs.mergeOverflow();
         hs2d.mergeOverflow();
         file->Close();
         
         //For multi save dss name
         dssName_multi=TString(dss.datasetName);
         
      }  // datasubset loop
   } // dataset loop
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs.combineFromSubsamples(samplesToCombine);
   hs2d.combineFromSubsamples(samplesToCombine);
   
   // Save 1d histograms
   TString loc=TString::Format("triggerEff/%s_%s_%s.root",currentSystematic.name().Data(),dssName_multi.Data(),cfg.treeVersion.Data());
   if(cfg.multi) loc=TString::Format("triggerEff/%s%s%s_%s.root",(currentSystematic.name()+"/").Data(),dssName_multi.Data(),(cfg.fileNR==0)?TString("").Data():TString("_"+std::to_string(cfg.fileNR)).Data(),cfg.treeVersion.Data());
   io::RootFileSaver saver_hist(loc,TString::Format("triggerEff%.1f",cfg.processFraction*100),false);
   hs.saveHistograms(saver_hist,samplesToCombine);
   hs2d.saveHistograms(saver_hist,samplesToCombine);
   hs.saveHistograms(saver_hist,saveDatasubsets);
   hs2d.saveHistograms(saver_hist,saveDatasubsets);
   
}
