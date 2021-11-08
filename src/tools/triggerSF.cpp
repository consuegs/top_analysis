#include "triggerSF.hpp"
#include "Config.hpp"
#include <iostream>

TriggerScaleFactors::TriggerScaleFactors(const std::string& fileName_ee,
                  const std::string& fileName_mumu,
                  const std::string& fileName_emu,
                  const Systematic::Systematic& systematic) :
h2_ee_(0),
h2_mumu_(0),
h2_emu_(0),
systematic_(systematic)
{
   std::cout<<"--- Beginning preparation of trigger scale factors\n";
   
   // Check if syst. variation should be used
   const Systematic::Type type = systematic_.type();
   if(systematic_.type()==Systematic::trig){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      std::cerr<<"Systematic not yet implemented"<<std::endl;
      exit(98);
   }
   else std::cout<<"Do not apply systematic variation\n";


   // Access scale factors
   h2_ee_ = this->prepareSF(fileName_ee, "eff_histo");
   h2_mumu_ = this->prepareSF(fileName_mumu, "eff_histo");
   h2_emu_ = this->prepareSF(fileName_emu, "eff_histo");

   std::cout<<"=== Finishing preparation of trigger scale factors\n\n";

}

TH2* TriggerScaleFactors::prepareSF(const std::string& sfInputFileName,
                               const std::string& histogramName)const
{
   
   io::RootFileReader triggerSF_ee(sfInputFileName.c_str(),"");

   // Access histogram containing scale factors
   TH2* h_scaleFactor(0);
   h_scaleFactor = (TH2*)(triggerSF_ee.read<TH2F>(histogramName));
   if(!h_scaleFactor){
      std::cerr<<"Error in TriggerScaleFactors::prepareSF()! TH2 for trigger factors not found in: "<<sfInputFileName
             <<"\n...break\n"<<std::endl;
      exit(39);
   }

   // Store histogram in memory and close file
   h_scaleFactor->SetDirectory(0);

   return h_scaleFactor;
}



const float TriggerScaleFactors::getTriggerSF(const float &pT_l1, const float &pT_l2,const int &channel, const bool &muonlead)
{
   switch(channel){
      case 1:
         return this->get2DSF(h2_ee_,pT_l1,pT_l2);
         break;
      case 2:
         return this->get2DSF(h2_mumu_,pT_l1,pT_l2);
         break;
      case 3:
         if(muonlead) return this->get2DSF(h2_emu_,pT_l1,pT_l2);
         else return this->get2DSF(h2_emu_,pT_l2,pT_l1);
         break;
      default:
         std::cerr<<"Error in TriggerScaleFactors::getTriggerSF()! Channel does not exist"<<std::endl;
         exit(49);
         return 1.;
   }
}

const float TriggerScaleFactors::get2DSF(const TH2* const histo, const float& x, const float& y)
{
   int xbin, ybin, dummy;
   histo->GetBinXYZ(histo->FindFixBin(x, y), xbin, ybin, dummy);
   //overflow to last bin
   xbin = std::min(xbin, histo->GetNbinsX());
   ybin = std::min(ybin, histo->GetNbinsY());
   float tempSF = histo->GetBinContent(xbin, ybin);
   return tempSF;
}

