//Script to derive bTag Efficiencies in MC
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"

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

Config const &cfg=Config::get();

extern "C"
void run()
{
   std::vector<std::string> vsDatasubsets(cfg.datasets.getDatasubsetNames());
   
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
   
   for (auto const &dss: cfg.datasets.getDatasubsets(true,true,true)){
      TFile file(dss.getPath(),"read");
      if (file.IsZombie()) {
         return;
      }
      io::log * ("Processing '"+dss.name+"' ");

      int year_int=cfg.year_int;
      float DeepCSV_loose=cfg.DeepCSV_loose;
      float DeepJet_loose=cfg.DeepJet_loose;
      
      hs2d.setCurrentSample(dss.name);
      
      //Check if current sample is TTbar powheg dilepton
      bool ttBar_dilepton=false;
      if (dss.datasetName=="TTbar_diLepton") ttBar_dilepton=true;
      
      //Check if current sample is TTbar amc@NLO
      bool ttBar_amc=false;
      if (dss.datasetName=="TTbar_amcatnlo") ttBar_amc=true;
      
      //Set Tree Input variables
      TTreeReader reader(cfg.treeName, &file);
      TTreeReaderValue<float> w_pu(reader, "pu_weight");
      TTreeReaderValue<UInt_t> runNo(reader, "runNo");
      TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
      TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
      TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
      TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
      TTreeReaderValue<float> w_topPT(reader, "topPTweight");
      TTreeReaderValue<float> w_bTag(reader, "bTagWeight");
      TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
      TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
      // ~TTreeReaderValue<std::vector<tree::Electron>> electrons_add(reader, "electrons_add");
      // ~TTreeReaderValue<std::vector<tree::Muon>>     muons_add    (reader, "muons_add");
      TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
      TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
      TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
      TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
      TTreeReaderValue<tree::MET> MET(reader, "met");
      TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
      TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
      TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
      TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
      TTreeReaderValue<tree::MET> MET_XYcorr(reader, "metXYcorr");
      TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
      TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
      TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
      TTreeReaderValue<int> n_Interactions(reader, "nPV");
      TTreeReaderValue<int> n_Interactions_gen(reader, "true_nPV");
      TTreeReaderValue<float> HTgen(reader, "genHt");
      TTreeReaderValue<bool> is_ee   (reader, "ee");
      TTreeReaderValue<bool> is_emu   (reader, "emu");
      TTreeReaderValue<bool> is_mumu   (reader, "mumu");
      TTreeReaderValue<bool> addLepton   (reader, "addLepton");
      TTreeReaderValue<float> mll   (reader, "mll");
      TTreeReaderValue<float> mt2   (reader, "MT2");
      TTreeReaderValue<float> genMT2   (reader, "genMT2");
      TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
      TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
      TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
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
      TTreeReaderValue<TLorentzVector> genTop(reader, "pseudoTop");
      TTreeReaderValue<TLorentzVector> genAntiTop(reader, "pseudoAntiTop");
      TTreeReaderValue<TLorentzVector> genLepton(reader, "pseudoLepton");
      TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "pseudoAntiLepton");
      TTreeReaderValue<TLorentzVector> genTau(reader, "pseudoTau");
      TTreeReaderValue<TLorentzVector> genAntiTau(reader, "pseudoAntiTau");
      TTreeReaderValue<int> genLeptonPdgId(reader, "pseudoLeptonPdgId");
      TTreeReaderValue<int> genAntiLeptonPdgId(reader, "pseudoAntiLeptonPdgId");
      TTreeReaderValue<TLorentzVector> genB(reader, "pseudoBJet");
      TTreeReaderValue<TLorentzVector> genAntiB(reader, "pseudoAntiBJet");
      TTreeReaderValue<TLorentzVector> genNeutrino(reader, "pseudoNeutrino");
      TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "pseudoAntiNeutrino");
      TTreeReaderValue<TLorentzVector> genWMinus(reader, "pseudoWMinus");
      TTreeReaderValue<TLorentzVector> genWPlus(reader, "pseudoWPlus");
      TTreeReaderValue<int> genDecayMode_pseudo(reader, "ttbarPseudoDecayMode");
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
         
         float met=MET->p.Pt();
         float const met_puppi=MET_Puppi->p.Pt();
         float const genMet=GENMET->p.Pt();
         
         //Do not use tau events in signal sample
         if (ttBar_dilepton && *genDecayMode>3) continue;
         
         //Do only use ee,emu,mumu in in amc ttbar
         if (ttBar_amc && (*genDecayMode>3 || *genDecayMode==0)) continue;
         
         //Trigger selection
         std::vector<bool> diElectronTriggers={*eleTrigg1,*eleTrigg2,*singleEleTrigg};
         std::vector<bool> diMuonTriggers={*muonTrigg1,*muonTrigg2,*muonTrigg3,*muonTrigg4,*singleMuonTrigg1,*singleMuonTrigg2};
         std::vector<bool> electronMuonTriggers={*eleMuTrigg1,*eleMuTrigg2,*eleMuTrigg3,*eleMuTrigg4,*singleMuonTrigg1,*singleMuonTrigg2,*singleEleTrigg};
         std::vector<bool> channel={*is_ee,*is_mumu,*is_emu};
         bool triggerMC=selection::triggerSelection(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,false,year_int);

         //Baseline selection (separation into ee, emu, mumu already done at TreeWriter)
         TLorentzVector p_l1;
         TLorentzVector p_l2;
         int flavor_l1=0;  //1 for electron and 2 for muon
         int flavor_l2=0;
         bool muonLead=true; //Boolean for emu channel
         TString cat="";
         
         if(!selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,cat,muonLead)) continue;
         
         if (triggerMC==false) continue;
                  
         std::vector<tree::Jet> cjets;
         std::vector<tree::Jet> BJets;
         std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
         
         // Perform selection up to bTag requirement (last entry in bool vector)
         if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end()-1, [](bool v) { return v; })) continue;
         
         TString path_cat="ee";
         if (*is_emu) path_cat="emu";
         else if (*is_mumu) path_cat="mumu";
         
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
      file.Close();
      
      //For multi save dss name
      dssName_multi=TString(dss.datasetName);
      
   } // dataset loop
   
   
   std::vector<TString> samplesToCombine=cfg.datasets.getDatasetNames();
   hs2d.combineFromSubsamples(samplesToCombine);
   
   // Save histograms
   TString loc=TString::Format("histograms_%s.root",cfg.treeVersion.Data());
   if(cfg.multi) loc=TString::Format("multiHists/histograms_%s_%s.root",dssName_multi.Data(),cfg.treeVersion.Data());
   io::RootFileSaver saver_hist(loc,TString::Format("bTagEff%.1f",cfg.processFraction*100),false);
   hs2d.saveHistograms(saver_hist,samplesToCombine);
   
}
