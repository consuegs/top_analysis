//Script to derive bTag Efficiencies in MC
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"
#include "tools/jetCorrections.hpp"

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

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/PyMethodBase.h"

#include "tools/leptonSF.hpp"
#include "tools/leptonCorrections.hpp"
#include "tools/triggerSF.hpp"
#include "tools/mcWeights.hpp"

Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
   //Read systematic from command line
   Systematic::Systematic currentSystematic(cfg.systematic);
   bool isNominal = (currentSystematic.type()==Systematic::nominal);
   
   // Configure JES/JER Corrections
   jesCorrections jesCorrector = jesCorrections(cfg.getJESPath(0,false).Data(),currentSystematic);
   jesCorrections jesCorrector_puppi = jesCorrections(cfg.getJESPath(0,true).Data(),currentSystematic);
   jerCorrections jerCorrector = jerCorrections(cfg.jer_SF_mc.Data(),cfg.jer_RES_mc.Data(),currentSystematic);
   
   // Configure lepton Correction
   leptonCorrections leptonCorretor = leptonCorrections(currentSystematic);
   
   // Configure leptonSF
   LeptonScaleFactors leptonSF = LeptonScaleFactors(cfg.electronID_file,cfg.electronID_hist,cfg.electronRECO_file,cfg.electronRECO_hist,
                                                cfg.muonID_file,cfg.muonID_hist,cfg.muonISO_file,cfg.muonISO_hist,currentSystematic);
   leptonSF.setDYExtrapolationUncFactors(cfg.muonDYunc,cfg.electronDYunc);
   
   //Configure topPT reweighting
   bool applytopPTreweighting = checkTopPTreweighting(currentSystematic);
   
   //Configure mcWeights (should be also powheg samples in this script!!)
   mcWeights mcWeighter = mcWeights(currentSystematic,false);
   
   TString dssName_multi="";
   
   //Define 2D histograms in the following
   hist::Histograms<TH2F> hs2d(vsDatasubsets);
   
   for(TString channel:{"ee","emu","mumu"}){
      hs2d.addHist("baseline/"+channel+"/B_all", ";p_{T}^{b-jet} (GeV);|#eta^{b-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/B_DeepJet_loose", ";p_{T}^{b-jet} (GeV);|#eta^{b-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/B_DeepCSV_loose", ";p_{T}^{b-jet} (GeV);|#eta^{b-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      
      hs2d.addHist("baseline/"+channel+"/C_all", ";p_{T}^{c-jet} (GeV);|#eta^{c-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/C_DeepJet_loose", ";p_{T}^{c-jet} (GeV);|#eta^{c-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/C_DeepCSV_loose", ";p_{T}^{c-jet} (GeV);|#eta^{c-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      
      hs2d.addHist("baseline/"+channel+"/Light_all", ";p_{T}^{light-jet} (GeV);|#eta^{light-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/Light_DeepJet_loose", ";p_{T}^{light-jet} (GeV);|#eta^{light-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
      hs2d.addHist("baseline/"+channel+"/Light_DeepCSV_loose", ";p_{T}^{light-jet} (GeV);|#eta^{light-jet}|;Jets/Bin" ,100,30,1000,100,0,2.4);
   }
   
   for (TString ds_name: cfg.datasets.getDatasetNames()){
      auto ds=cfg.datasets.getDataset(ds_name);
      // Check if current systematic and sample match
      Systematic::checkAlternativeSample(currentSystematic,ds.systName,ds.name);
   }
   
   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile* file = TFile::Open(dss.getPath(),"read");
      if (file->IsZombie()) {
         return;
      }
      
      io::log * ("Processing '"+dss.name+"' ");

      int year_int=cfg.year_int;
      float DeepCSV_loose=cfg.bTagWPcut_alt;
      float DeepJet_loose=cfg.bTagWPcut;
      
      hs2d.setCurrentSample(dss.name);
      
      //Check if current sample is TTbar 2L sample (later used to veto tau events)
      bool ttBar_dilepton=dss.isTTbar2L;
      
      //Check if current sample is TTbar amc@NLO
      bool ttBar_amc=false;
      if (dss.datasetName=="TTbar_amcatnlo") ttBar_amc=true;
      
      //Set Tree Input variables
      TTreeReader reader(cfg.treeName, file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<UInt_t> runNo(reader, "runNo");
      TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
      TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
      TTreeReaderValue<float> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
      TTreeReaderValue<std::vector<float>> w_ps(reader, "ps_weights");
      TTreeReaderValue<float> w_topPT(reader, "topPTweight");
      TTreeReaderValue<float> fragUpWeight(reader, "weightFragUp");
      TTreeReaderValue<float> fragCentralWeight(reader, "weightFragCentral");
      TTreeReaderValue<float> fragDownWeight(reader, "weightFragDown");
      TTreeReaderValue<float> fragPetersonWeight(reader, "weightFragPeterson");
      TTreeReaderValue<float> semilepbrUpWeight(reader, "weightSemilepbrUp");
      TTreeReaderValue<float> semilepbrDownWeight(reader, "weightSemilepbrDown");
      TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<std::vector<tree::Jet>>      jets_puppi     (reader, "jetsPuppi");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<float> rho(reader, "rho");
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
      TTreeReaderValue<int> genDecayMode(reader, "ttbarDecayMode");
   
   
      int iEv=0;
      int processEvents=cfg.processFraction*dss.entries;
      while (reader.Next()){
         iEv++;
         if (iEv>processEvents) break;
         if (iEv%(std::max(processEvents/10,1))==0){
            io::log*".";
            io::log.flush();
         }
         
         //Do not use tau events in signal sample
         if (ttBar_dilepton && *genDecayMode>3) continue;
         
         //Do only use ee,emu,mumu in in amc ttbar
         if (ttBar_amc && (*genDecayMode>3 || *genDecayMode==0)) continue;
         
         // Correct and select leptons
         *muons = leptonCorretor.correctMuons(*muons,MET->p,MET_Puppi->p);
         *electrons = leptonCorretor.correctElectrons(*electrons,MET->p,MET_Puppi->p);
         
         //Baseline selection (including separation into ee, emu, mumu)
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         int flavor_l1=0;  //1 for electron and 2 for muon
         int flavor_l2=0;
         bool muonLead=true; //Boolean for emu channel
         std::vector<bool> channel={false,false,false};     //ee, mumu, emu
         TString cat="";
         
         if(!selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,cat,muonLead)) continue;
         
         //Trigger selection
         std::vector<bool> diElectronTriggers={*eleTrigg1,*eleTrigg2,*singleEleTrigg};
         std::vector<bool> diMuonTriggers={*muonTrigg1,*muonTrigg2,*muonTrigg3,*muonTrigg4,*singleMuonTrigg1,*singleMuonTrigg2};
         std::vector<bool> electronMuonTriggers={*eleMuTrigg1,*eleMuTrigg2,*eleMuTrigg3,*eleMuTrigg4,*singleMuonTrigg1,*singleMuonTrigg2,*singleEleTrigg};
         bool triggerMC=selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,false,year_int);

         if (triggerMC==false) continue;
         
         //Apply JES and JER systematics
         jesCorrector.applySystematics(*jets,MET->p);
         jesCorrector_puppi.applySystematics(*jets_puppi,MET_Puppi->p);    // Needed for correction of Puppi MET
         jerCorrector.smearCollection_Hybrid(*jets,*rho);
         
         float met=MET->p.Pt();
         float const met_puppi=MET_Puppi->p.Pt();
         float const genMet=GENMET->p.Pt();
         
         std::vector<tree::Jet> cjets;
         std::vector<tree::Jet> BJets;
         std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
         
         // Perform selection up to bTag requirement (last entry in bool vector)
         if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end()-1, [](bool v) { return v; })) continue;
         
         // Get weights
         if (!applytopPTreweighting) *w_topPT = 1.;   //For topPt syst set weight to 1
         int channelID = 0;
         if (channel[0]) channelID = 1;
         else if (channel[1]) channelID = 2;
         else if (channel[2]) channelID = 3;
         float leptonSFweight = leptonSF.getSFDilepton(p_l1,p_l2,flavor_l1,flavor_l2);
         // ~float triggerSF = triggerSFcalc.getTriggerSF(p_l1.Pt(),p_l2.Pt(),channelID,muonLead);    //Not used since bTagEff needed in triggerSF calc.
         float mcWeight = mcWeighter.getMCweight(*w_mc,*w_pdf,*w_ps,{*fragUpWeight,*fragCentralWeight,*fragDownWeight,*fragPetersonWeight,*semilepbrUpWeight,*semilepbrDownWeight});
         float fEventWeight=*w_pu * mcWeight;     //Set event weight 
         float SFWeight=leptonSFweight * *w_topPT;     //Set combined SF weight
         hs2d.setFillWeight(fEventWeight*SFWeight);
         
         TString path_cat="ee";
         if (channel[2]) path_cat="emu";
         else if (channel[1]) path_cat="mumu";
         
         for(auto jet : cjets){
            int flavor=jet.hadronFlavour;
            switch(flavor){
               case 5:
                  hs2d.fill("baseline/"+path_cat+"/B_all",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepJet>DeepJet_loose) hs2d.fill("baseline/"+path_cat+"/B_DeepJet_loose",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepCSV>DeepCSV_loose) hs2d.fill("baseline/"+path_cat+"/B_DeepCSV_loose",jet.p.Pt(),abs(jet.p.Eta()));
                  break;
               case 4:
                  hs2d.fill("baseline/"+path_cat+"/C_all",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepJet>DeepJet_loose) hs2d.fill("baseline/"+path_cat+"/C_DeepJet_loose",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepCSV>DeepCSV_loose) hs2d.fill("baseline/"+path_cat+"/C_DeepCSV_loose",jet.p.Pt(),abs(jet.p.Eta()));
                  break;
               default:
                  hs2d.fill("baseline/"+path_cat+"/Light_all",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepJet>DeepJet_loose) hs2d.fill("baseline/"+path_cat+"/Light_DeepJet_loose",jet.p.Pt(),abs(jet.p.Eta()));
                  if(jet.bTagDeepCSV>DeepCSV_loose) hs2d.fill("baseline/"+path_cat+"/Light_DeepCSV_loose",jet.p.Pt(),abs(jet.p.Eta()));
            }
         }
         
      
      }// evt loop
      io::log<<"";
      
      hs2d.mergeOverflow();
      file->Close();
      
      //For multi save dss name
      dssName_multi=TString(dss.datasetName);
      
   } // dataset loop
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs2d.combineFromSubsamples(samplesToCombine);
   
   // Save histograms
   TString loc=TString::Format("bTagEff/%s_%s_%s.root",currentSystematic.name().Data(),dssName_multi.Data(),cfg.treeVersion.Data());
   if(cfg.multi) loc=TString::Format("bTagEff/%s%s%s_%s.root",(currentSystematic.name()+"/").Data(),dssName_multi.Data(),(cfg.fileNR==0)?TString("").Data():TString("_"+std::to_string(cfg.fileNR)).Data(),cfg.treeVersion.Data());
   io::RootFileSaver saver_hist(loc,TString::Format("bTagEff%.1f",cfg.processFraction*100),false);
   hs2d.saveHistograms(saver_hist,samplesToCombine);
   
}
