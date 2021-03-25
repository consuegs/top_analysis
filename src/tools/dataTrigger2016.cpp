#include "dataTrigger2016.hpp"

/*
   Trigger Selection Function for 2016
   The vectors have the following structure
   diElectronTriggers: 
                        HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v     (placeholder, since other years have more triggers)
                        HLT_Ele27_WPTight_Gsf_v
   diMuonTriggers:
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
                        HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v
                        HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v
                        HLT_IsoMu24_v
                        HLT_IsoTkMu24_v
   electronMuonTriggers:
                        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
                        HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v
                        HLT_IsoMu24_v
                        HLT_IsoTkMu24_v
                        HLT_Ele27_WPTight_Gsf_v
   
   channel: ee,mumu,emu
   
   PD:
         SingleElectron
         DoubleEG
         SingleMuon
         DoubleMuon
         MuonEG
 */

bool dataTrigger2016::DataTriggerSelection2016(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                        std::vector<bool> const &channel, std::vector<bool> const &PD, bool const &is2016H)
{
   bool triggerData;
   bool triggerData_veto;
   
   //ee
   if(PD[0] && channel[0]){   //SingleElectron
      triggerData=diElectronTriggers[2];
      triggerData_veto=diElectronTriggers[0];
   }
   else if(PD[1] && channel[0]){    //DoubleEG
      triggerData=diElectronTriggers[0];
      triggerData_veto=false;
   }
   //mumu
   else if(PD[2] && channel[1]){    //SingleMuon
      triggerData=diMuonTriggers[4] || diMuonTriggers[5];
      if (is2016H) triggerData_veto=diMuonTriggers[0] || diMuonTriggers[1];
      else triggerData_veto=diMuonTriggers[0] || diMuonTriggers[1] || diMuonTriggers[2] || diMuonTriggers[3];
   }
   else if(PD[3] && channel[1]){    //DoubleMuon
      if (is2016H) triggerData=diMuonTriggers[0] || diMuonTriggers[1];
      else triggerData=diMuonTriggers[0] || diMuonTriggers[1] || diMuonTriggers[2] || diMuonTriggers[3];
      triggerData_veto=false;
   }
   //emu
   else if(PD[2] && channel[2]){    //SingleMuon
      triggerData=electronMuonTriggers[4] || electronMuonTriggers[5];
      if (is2016H) triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1] || diElectronTriggers[2];
      else triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3] || diElectronTriggers[2];
   }
   else if(PD[0] && channel[2]){    //SingleElectron
      triggerData=diElectronTriggers[2];
      if (is2016H) triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1];
      else triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3];
   }
   else if(PD[4] && channel[2]){    //MuonEG
      if (is2016H) triggerData=electronMuonTriggers[0] || electronMuonTriggers[1];
      else triggerData=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3];
      triggerData_veto=false;
   }
   else return false;
   
   if((triggerData && !triggerData_veto)) return true;
   else return false;
}
