#include "selection.hpp"
#include "Config.hpp"
#include "dataTrigger2016.hpp"
#include "dataTrigger2017.hpp"
#include "dataTrigger2018.hpp"
#include <iostream>

#include <limits>

Config const &cfg=Config::get();

bool selection::triggerSelection(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                                 std::vector<bool> const &channel, bool const &isData, int const year, std::vector<bool> const &PD, bool const &is2016H, bool const &is2017AB)
{
   bool diElectronTriggers_combined=!std::none_of(diElectronTriggers.begin(), diElectronTriggers.end(), [](bool v) { return v; });
   bool diMuonTriggers_combined=!std::none_of(diMuonTriggers.begin(), diMuonTriggers.end(), [](bool v) { return v; });
   bool electronMuonTriggers_combined=!std::none_of(electronMuonTriggers.begin(), electronMuonTriggers.end(), [](bool v) { return v; });
   
   bool data;
   
   if (!isData){
      if(!(channel[0] && diElectronTriggers_combined) && !(channel[1] && diMuonTriggers_combined) && !(channel[2] && electronMuonTriggers_combined)) return false;
      else return true;
   }
   else{
      switch(year)
      {
         case 1:
            data=dataTrigger2016::DataTriggerSelection2016(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,PD,is2016H);
            break;
         case 2:
            data=dataTrigger2017::DataTriggerSelection2017(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,PD,is2017AB);
            break;
         case 3:
            data=dataTrigger2018::DataTriggerSelection2018(diElectronTriggers,diMuonTriggers,electronMuonTriggers,channel,PD);
            break;
      }
      return data;
   }
}

bool selection::diLeptonSelection(std::vector<tree::Electron> const &electrons, std::vector<tree::Muon> const &muons, std::vector<bool> &channel,
                        TLorentzVector &p_l1, TLorentzVector &p_l2, int &flavor_l1, int &flavor_l2, TString &cat, bool &muonLead)
{
   
   if ((electrons.size()+muons.size()) == 2){
      if (electrons.size() == 2){
         if (electrons[0].charge*electrons[1].charge!=-1) return false;
         else{
            channel[0] = true;
            p_l1 = electrons[0].p;
            p_l2 = electrons[1].p;
            flavor_l1 = 1;
            flavor_l2 = 1;
            cat = "ee";
         }
      }
      else if (muons.size() == 2){
         if (muons[0].charge*muons[1].charge != -1) return false;
         else{
            channel[1] = true;
            p_l1 = muons[0].p;
            p_l2 = muons[1].p;
            flavor_l1 = 2;
            flavor_l2 = 2;
            cat = "mumu";
         }
      }
      else {
         if (muons[0].charge*electrons[0].charge != -1) return false;
         else{
            if (muons[0].p.Pt() > electrons[0].p.Pt()){
               muonLead = true;
               p_l1 = muons[0].p;
               p_l2 = electrons[0].p;
               flavor_l1 = 2;
               flavor_l2 = 1;
            }
            else{
               muonLead = false;
               p_l1 = electrons[0].p;
               p_l2 = muons[0].p;
               flavor_l1 = 1;
               flavor_l2 = 2;
            }
            channel[2] = true;
            cat = "emu";
         }
      }
   }
   else return false;
   
   return (p_l1.Pt()>=25 && p_l2.Pt()>=20);
}

std::vector<bool> selection::ttbarSelection(TLorentzVector const &p_l1, TLorentzVector const &p_l2, float const &met, std::vector<bool> const &channel,
                                    std::vector<tree::Jet> const &jets, std::vector<tree::Jet> &cleanJets, std::vector<tree::Jet> &bJets)
{
   std::vector<bool> selection_vec={false,false,false,false};
      
   //mLL Cut
   float mll_corr=(p_l1+p_l2).M();
   if(mll_corr<20 || ((channel[0] || channel[1]) && mll_corr<106 && mll_corr>76)) return selection_vec;
   else selection_vec[0]=true;
   
   //Jet Cut
   cleanJets=phys::getCleanedJets(jets, p_l1, p_l2);
   if(cleanJets.size()<2) return selection_vec;
   else selection_vec[1]=true;
   
   //MET Cut
   if ((channel[0] || channel[1]) && met<40) return selection_vec;
   else selection_vec[2]=true;
   
   //bJet Cut
   bool bTag=false;
   if(cfg.year_int==1){
      for (tree::Jet const &jet : cleanJets) {
         if (jet.bTagDeepCSV>cfg.bTagWPcut_alt) {      //Loose working point for deepCSV
            bTag=true;
            bJets.push_back(jet);
         }
      }
   }
   else{
      for (tree::Jet const &jet : cleanJets) {
         if (jet.bTagDeepJet>cfg.bTagWPcut) {      //Loose working point for deepCSV
            bTag=true;
            bJets.push_back(jet);
         }
      }
   }
   if(!bTag) return selection_vec;
   else{
      selection_vec[3]=true;
      return selection_vec;
   }
}

std::vector<bool> selection::kitSyncSelection(TLorentzVector const &p_l1, TLorentzVector const &p_l2, float const &met, std::vector<bool> const &channel,
                                    std::vector<tree::Jet> const &jets, std::vector<tree::Jet> &cleanJets, std::vector<tree::Jet> &bJets)
{
   std::vector<bool> selection_vec={false};
   
   //Still check jets
   cleanJets=phys::getCleanedJets(jets,p_l1,p_l2);
   
   bool bTag=false;
   if(cfg.year_int==1){
      for (tree::Jet const &jet : cleanJets) {
         if (jet.bTagDeepCSV>cfg.bTagWPcut_alt) {      //Loose working point for deepCSV
            bTag=true;
            bJets.push_back(jet);
         }
      }
   }
   else{
      for (tree::Jet const &jet : cleanJets) {
         if (jet.bTagDeepJet>cfg.bTagWP) {      //Loose working point for deepCSV
            bTag=true;
            bJets.push_back(jet);
         }
      }
   }
   
   //mLL Cut
   float mll_corr=(p_l1+p_l2).M();
   if(mll_corr>106 || mll_corr<76) return selection_vec;
   else {
      selection_vec[0]=true;
      return selection_vec;
   }
}



