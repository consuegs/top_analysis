//Partly copied from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/common/src/ScaleFactors.cc
#include "leptonSF.hpp"
#include "Config.hpp"
#include <iostream>

LeptonScaleFactors::LeptonScaleFactors(const std::string& electronFileNameID,
                  const std::string& electronHistNameID,
                  const std::string& electronFileNameRECO,
                  const std::string& electronHistNameRECO,
                  const std::string& muonFileNameID,
                  std::string muonHistNameID,
                  const std::string& muonFileNameISO,
                  std::string muonHistNameISO,
                  const Systematic::Systematic& systematic) :
muon_ISO_DY_unc_(0.),
electron_ID_DY_unc_(0.),
h2_electron_ID_histo_(0),
h2_electron_Reco_histo_(0),
h2_muon_ID_histo_(0),
h2_muon_ISO_histo_(0),
systematic_(systematic),
uncFactor_electron_RECO_(0.),
uncFactor_electron_ID_(0.),
uncFactor_muon_ISO_(0.),
uncFactor_muon_ID_(0.)
{
   std::cout<<"--- Beginning preparation of lepton scale factors\n";
   
   // Check if syst. variation should be used
   const Systematic::Type type = systematic_.type();
   if(std::find(Systematic::leptonsfTypes.begin(), Systematic::leptonsfTypes.end(), type) != Systematic::leptonsfTypes.end()){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic_.type()==Systematic::eleID) (systematic_.variation()==Systematic::up)? uncFactor_electron_ID_=1. : uncFactor_electron_ID_=-1.;
      else if(systematic_.type()==Systematic::eleReco) (systematic_.variation()==Systematic::up)? uncFactor_electron_RECO_=1. : uncFactor_electron_RECO_=-1.;
      else if(systematic_.type()==Systematic::muonID) (systematic_.variation()==Systematic::up)? uncFactor_muon_ID_=1. : uncFactor_muon_ID_=-1.;
      else if(systematic_.type()==Systematic::muonIso) (systematic_.variation()==Systematic::up)? uncFactor_muon_ISO_=1. : uncFactor_muon_ISO_=-1.;
      else if(systematic_.type()==Systematic::muonIDStat){
         (systematic_.variation()==Systematic::up)? uncFactor_muon_ID_=1. : uncFactor_muon_ID_=-1.;
         muonHistNameID += "_stat";
      }
      else if(systematic_.type()==Systematic::muonIDSyst){
         (systematic_.variation()==Systematic::up)? uncFactor_muon_ID_=1. : uncFactor_muon_ID_=-1.;
         muonHistNameID += "_syst";
      }
      else if(systematic_.type()==Systematic::muonIsoStat){
         (systematic_.variation()==Systematic::up)? uncFactor_muon_ISO_=1. : uncFactor_muon_ISO_=-1.;
         muonHistNameISO += "_stat";
      }
      else if(systematic_.type()==Systematic::muonIsoSyst){
         (systematic_.variation()==Systematic::up)? uncFactor_muon_ISO_=1. : uncFactor_muon_ISO_=-1.;
         muonHistNameISO += "_syst";
      }
      
   }
   else std::cout<<"Do not apply systematic variation\n";


   // Access electron scale factor
   std::cout<<"Electron:\n";
   h2_electron_ID_histo_ = this->prepareSF(electronFileNameID, electronHistNameID);
   if(h2_electron_ID_histo_) std::cout<<"ID scale factors found - will be used\n";
   h2_electron_Reco_histo_ = this->prepareSF(electronFileNameRECO, electronHistNameRECO);
   if(h2_electron_Reco_histo_) std::cout<<"Reconstruction scale factors found - will be used\n";
   
   // Access muon scale factor
   std::cout<<"Muon:\n";
   h2_muon_ID_histo_ = this->prepareSF(muonFileNameID, muonHistNameID);
   if(h2_muon_ID_histo_) std::cout<<"ID scale factors found - will be used\n";
   h2_muon_ISO_histo_ = this->prepareSF(muonFileNameISO, muonHistNameISO);
   if(h2_muon_ISO_histo_) std::cout<<"Iso scale factors found - will be used\n";

}

TH2* LeptonScaleFactors::prepareSF(const std::string& sfInputFileName,
                               const std::string& histogramName)const
{
   
   // Check if file exists
   if (!io::fileExists(sfInputFileName)){
      std::cerr<<"ERROR in LeptonScaleFactors::prepareSF! "
      <<sfInputFileName<<" does not exist\n...break\n"<<std::endl;
      exit(98);
   }

   TFile scaleFactorFile(sfInputFileName.c_str());

   // Access histogram containing scale factors
   TH2* h_scaleFactorPtEta(0);
   h_scaleFactorPtEta = dynamic_cast<TH2*>(scaleFactorFile.Get(histogramName.c_str()));
   if(!h_scaleFactorPtEta){
      std::cerr<<"Error in LeptonScaleFactors::prepareSF()! TH2 for lepton Id/Iso scale factors not found: "<<histogramName
             <<"\n...break\n"<<std::endl;
      exit(39);
   }

   // Store histogram in memory and close file
   h_scaleFactorPtEta->SetDirectory(0);
   scaleFactorFile.Close();

   return h_scaleFactorPtEta;
}



const float LeptonScaleFactors::getSFDilepton(const TLorentzVector &p_l1, const TLorentzVector &p_l2,
                                             const int &flavor_l1, const int &flavor_l2,
                                             const double &etaSC_l1, const double &etaSC_l2)
{
   float result(1.);

   result *= this->leptonSF(p_l1, flavor_l1, etaSC_l1);
   result *= this->leptonSF(p_l2, flavor_l2, etaSC_l2);
   return result;
}

const float LeptonScaleFactors::get2DSF(const TH2* const histo, const float& x, const float& y, const float uncFactor = 0., const float DYunc = 0.)
{
   int xbin, ybin, dummy;
   histo->GetBinXYZ(histo->FindFixBin(x, y), xbin, ybin, dummy);
   //overflow to last bin
   xbin = std::min(xbin, histo->GetNbinsX());
   ybin = std::min(ybin, histo->GetNbinsY());
   float tempSF = histo->GetBinContent(xbin, ybin);
   return (tempSF + uncFactor*histo->GetBinError(xbin, ybin) + uncFactor*DYunc*tempSF);
}

const float LeptonScaleFactors::leptonSF(const TLorentzVector p, const int flavor, const double etaSC)
{
   float result(1.);
   if(flavor==1) {
      result *= get2DSF(h2_electron_Reco_histo_, p.Eta(), p.Pt(), uncFactor_electron_RECO_);
      result *= get2DSF(h2_electron_ID_histo_, etaSC , p.Pt(), uncFactor_electron_ID_, electron_ID_DY_unc_);
   }
   else if (flavor==2){
      result *= get2DSF(h2_muon_ISO_histo_, fabs(p.Eta()), p.Pt(), uncFactor_muon_ISO_, muon_ISO_DY_unc_);
      result *= get2DSF(h2_muon_ID_histo_, fabs(p.Eta()), p.Pt(), uncFactor_muon_ID_);
   }

   return result;
}

void LeptonScaleFactors::setDYExtrapolationUncFactors(const float& muon_ISO_DY_unc, const float& electron_ID_DY_unc){
   
   if (systematic_.type()==Systematic::muonIso || systematic_.type()==Systematic::muonIsoSyst){
      muon_ISO_DY_unc_ = muon_ISO_DY_unc;
   }
   else if (systematic_.type()==Systematic::eleID){
      electron_ID_DY_unc_ = electron_ID_DY_unc;
   }
}


