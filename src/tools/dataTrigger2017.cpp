#include "dataTrigger2017.hpp"

/*
   Trigger Selection Function for 2017
   The vectors have the following structure
   diElectronTriggers: 
                        HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_v
                        HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
   diMuonTriggers:
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v   
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v     (placeholder, since other years have more triggers)
                        HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8_v     (placeholder, since other years have more triggers)
                        HLT_IsoMu27_v        
                        HLT_IsoMu27_v                                     (placeholder, since other years have more triggers)
   electronMuonTriggers:
                        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v
                        HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v
                        HLT_IsoMu27_v
                        HLT_IsoMu27_v                                     (placeholder, since other years have more triggers)
                        HLT_Ele32_WPTight_Gsf_L1DoubleEG_v
   
   channel: ee,mumu,emu
   
   PD:
         SingleElectron
         DoubleEG
         SingleMuon
         DoubleMuon
         MuonEG
         EGamma               (Placeholder for 2018)
 */

bool dataTrigger2017::DataTriggerSelection2017(std::vector<bool> const &diElectronTriggers, std::vector<bool> const &diMuonTriggers, std::vector<bool> const &electronMuonTriggers,
                        std::vector<bool> const &channel, std::vector<bool> const &PD, bool const &is2017AB)
{
   bool triggerData;
   bool triggerData_veto;
   
   //ee
   if(PD[0] && channel[0]){   //SingleElectron
      triggerData=diElectronTriggers[2];
      triggerData_veto=diElectronTriggers[0] || diElectronTriggers[1];
   }
   else if(PD[1] && channel[0]){    //DoubleEG
      triggerData=diElectronTriggers[0] || diElectronTriggers[1];
      triggerData_veto=false;
   }
   //mumu
   else if(PD[2] && channel[1]){    //SingleMuon
      triggerData=diMuonTriggers[4];
      if (is2017AB) triggerData_veto=diMuonTriggers[0];
      else triggerData_veto=diMuonTriggers[1];
   }
   else if(PD[3] && channel[1]){    //DoubleMuon
      if (is2017AB) triggerData=diMuonTriggers[0];
      else triggerData=diMuonTriggers[1];
      triggerData_veto=false;
   }
   //emu
   else if(PD[2] && channel[2]){    //SingleMuon
      triggerData=electronMuonTriggers[4];
      triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3] || diElectronTriggers[2];
   }
   else if(PD[0] && channel[2]){    //SingleElectron
      triggerData=diElectronTriggers[2];
      triggerData_veto=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3];
   }
   else if(PD[4] && channel[2]){    //MuonEG
      triggerData=electronMuonTriggers[0] || electronMuonTriggers[1] || electronMuonTriggers[2] || electronMuonTriggers[3];
      triggerData_veto=false;
   }
   else return false;
   
   if((triggerData && !triggerData_veto)) return true;
   else return false;
}
