#include "selection.hpp"
#include "dataTrigger2016.hpp"
#include "dataTrigger2017.hpp"
#include "dataTrigger2018.hpp"
#include <iostream>

#include <limits>

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


/*
   Lepton Selection Function for 2016
*/

bool selection::diLeptonSelection(std::vector<tree::Electron> const &electrons, std::vector<tree::Muon> const &muons, std::vector<bool> const &channel,
                        TLorentzVector &p_l1, TLorentzVector &p_l2, int &flavor_l1, int &flavor_l2, TString &cat, bool &muonLead)
{
   int leadLepton=0;
   int subleadLepton=1;
   bool rec_selection=false;
   
   if (channel[0]){
      rec_selection=true;
      if(!electrons[0].isTight || !electrons[1].isTight) rec_selection=false; //currently double check since trees only have tight leptons!!
      if(abs(electrons[0].etaSC)>2.4 || abs(electrons[1].etaSC>2.4)) rec_selection=false; //To use same region as for muons, cut on supercluster eta
      if(electrons[0].p.Pt()*electrons[0].corr<electrons[1].p.Pt()*electrons[1].corr) leadLepton=1,subleadLepton=0;
      p_l1=electrons[leadLepton].p*electrons[leadLepton].corr;
      p_l2=electrons[subleadLepton].p*electrons[subleadLepton].corr;
      p_l1=electrons[leadLepton].p*electrons[leadLepton].corr;
      p_l2=electrons[subleadLepton].p*electrons[subleadLepton].corr;
      flavor_l1=1;
      flavor_l2=1;
      cat="ee";
   }
   else if (channel[1]){
      rec_selection=true;
      if(!muons[0].isTight || !muons[1].isTight) rec_selection=false;
      if(muons[0].rIso>0.15 || muons[1].rIso>0.15) rec_selection=false;
      if(abs(muons[0].p.Eta())>2.4 || abs(muons[1].p.Eta())>2.4) rec_selection=false;
      if(muons[0].p.Pt()*muons[0].rochesterCorrection<muons[1].p.Pt()*muons[1].rochesterCorrection) leadLepton=1,subleadLepton=0;
      p_l1=muons[leadLepton].p*muons[leadLepton].rochesterCorrection;
      p_l2=muons[subleadLepton].p*muons[subleadLepton].rochesterCorrection;
      flavor_l1=2;
      flavor_l2=2;
      cat="mumu";
   }
   else if (channel[2]){
      rec_selection=true;
      if(!muons[0].isTight || !electrons[0].isTight) rec_selection=false;
      if(muons[0].rIso>0.15 ) rec_selection=false;
      if(abs(muons[0].p.Eta())>2.4) rec_selection=false;
      if(abs(electrons[0].etaSC>2.4) ) rec_selection=false;
      if (muons[0].p.Pt()*muons[0].rochesterCorrection>electrons[0].p.Pt()*electrons[0].corr){
         p_l1=muons[0].p*muons[0].rochesterCorrection;
         p_l2=electrons[0].p*electrons[0].corr;
         flavor_l1=2;
         flavor_l2=1;
      }
      else {
         p_l1=electrons[0].p*electrons[0].corr;
         p_l2=muons[0].p*muons[0].rochesterCorrection;
         flavor_l1=1;
         flavor_l2=2;
         muonLead=false;
      }
      cat="emu";
   }
   
   if (p_l1.Pt()<25 || p_l2.Pt()<20) rec_selection=false;
   
   return rec_selection;
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
   cleanJets=phys::getCleanedJets(jets);
   if(cleanJets.size()<2) return selection_vec;
   else selection_vec[1]=true;
   
   //MET Cut
   if ((channel[0] || channel[1]) && met<40) return selection_vec;
   else selection_vec[2]=true;
   
   //bJet Cut
   bool bTag=false;
   for (tree::Jet const &jet : cleanJets) {
      if (jet.bTagDeepCSV>0.2217) {      //Loose working point for deepCSV
         bTag=true;
         bJets.push_back(jet);
      }
   }
   if(!bTag) return selection_vec;
   else{
      selection_vec[3]=true;
      return selection_vec;
   }
}

std::vector<bool> selection::ttbarSelection_looseJetID(TLorentzVector const &p_l1, TLorentzVector const &p_l2, float const &met, std::vector<bool> const &channel,
                                    std::vector<tree::Jet> const &jets, std::vector<tree::Jet> &cleanJets, std::vector<tree::Jet> &bJets)
{
   std::vector<bool> selection_vec={false,false,false,false};
   
   //mLL Cut
   float mll_corr=(p_l1+p_l2).M();
   if(mll_corr<20 || ((channel[0] || channel[1]) && mll_corr<106 && mll_corr>76)) return selection_vec;
   else selection_vec[0]=true;
   
   //Jet Cut
   cleanJets=phys::getCleanedJets_looseID(jets);
   if(cleanJets.size()<2) return selection_vec;
   else selection_vec[1]=true;
   
   //MET Cut
   if ((channel[0] || channel[1]) && met<40) return selection_vec;
   else selection_vec[2]=true;
   
   //bJet Cut
   bool bTag=false;
   for (tree::Jet const &jet : cleanJets) {
      if (jet.bTagDeepCSV>0.2217) {      //Loose working point for deepCSV
         bTag=true;
         bJets.push_back(jet);
      }
   }
   if(!bTag) return selection_vec;
   else{
      selection_vec[3]=true;
      return selection_vec;
   }
}



