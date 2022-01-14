#include "triggerSF.hpp"
#include "Config.hpp"
#include <iostream>

TriggerScaleFactors::TriggerScaleFactors(const std::string& fileName,
                  const Systematic::Systematic& systematic) :
h2_ee_(0),
h2_mumu_(0),
h2_emu_(0),
systematic_(systematic)
{
   std::cout<<"--- Beginning preparation of trigger scale factors\n";
   
   // Check if syst. variation should be used
   systematicInternal_ = nominal;
   const Systematic::Type type = systematic_.type();
   if(systematic_.type()==Systematic::trig){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up){ systematicInternal_ = vary_up; }
      else if(systematic.variation() == Systematic::down){ systematicInternal_ = vary_down; }
      else {
         std::cerr << "ERROR in constructor of triggerSF! Systematic variation is invalid: "
                  << Systematic::convertVariation(systematic.variation()) << "\n...break\n\n";
         exit(98);
      }
      if(systematicInternal_ == vary_up) {std::cout<<"Apply systematic variation: up\n";}
      else if(systematicInternal_ == vary_down) std::cout<<"Apply systematic variation: down\n";
   }
   else std::cout<<"Do not apply systematic variation\n";


   // Access scale factors
   if (systematicInternal_ == nominal || systematicInternal_ == vary_up){
      h2_ee_ = this->prepareSF(fileName, "ee_SF_totalUp");
      h2_mumu_ = this->prepareSF(fileName, "mumu_SF_totalUp");
      h2_emu_ = this->prepareSF(fileName, "emu_SF_totalUp");
   }
   else if (systematicInternal_ == vary_down){
      h2_ee_ = this->prepareSF(fileName, "ee_SF_totalDown");
      h2_mumu_ = this->prepareSF(fileName, "mumu_SF_totalDown");
      h2_emu_ = this->prepareSF(fileName, "emu_SF_totalDown");
   }

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
   float tempSF = 1.;
   if (systematicInternal_ == nominal) tempSF = histo->GetBinContent(xbin, ybin);
   else if (systematicInternal_ == vary_up) tempSF = histo->GetBinContent(xbin, ybin)+histo->GetBinError(xbin, ybin);
   else tempSF = histo->GetBinContent(xbin, ybin)-histo->GetBinError(xbin, ybin);
   return tempSF;
}

