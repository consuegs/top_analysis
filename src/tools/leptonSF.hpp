//Partly copied from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/common/include/ScaleFactors.h
#ifndef LEPTONSF_HPP__
#define LEPTONSF_HPP__

#include "tree/TreeParticles.hpp"
#include "tools/io.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

/// Class for application of lepton scale factors
class LeptonScaleFactors{

public:

   /// Constructor
   /// If input for specific lepton type is empty, set for this lepton type scale factors = 1.
   LeptonScaleFactors(const std::string& electronFileNameID,
                  const std::string& electronHistNameID,
                  const std::string& electronFileNameRECO,
                  const std::string& electronHistNameRECO,
                  const std::string& muonFileNameID,
                  std::string muonHistNameID,
                  const std::string& muonFileNameISO,
                  std::string muonHistNameISO,
                  const Systematic::Systematic& systematic);

   /// Destructor
   ~LeptonScaleFactors(){}



   /// Get lepton per-event scale factor for exactly two leptons
   const float getSFDilepton(const TLorentzVector &p_l1, const TLorentzVector &p_l2,
                   const int &flavor_l1, const int &flavor_l2);

   void setDYExtrapolationUncFactors(const float& muon_ISO_DY_unc, const float& electron_ID_DY_unc);


private:


   /// Enumeration for possible systematics
   enum SystematicInternal{nominal, vary_up, vary_down};
   
   const Systematic::Systematic systematic_;

   /// Return the scale factor histogram
   TH2* prepareSF(const std::string& inputFileName,
                   const std::string& histogramName)const;

   /// Get per-lepton scale factor
   const float leptonSF(const TLorentzVector p, const int flavor);
   
   /// Get SF from 2D hist
   const float get2DSF(const TH2* const histo, const float& x, const float& y, const float uncFactor, const float DYunc);
   
   /// Ingredients for DY extrapolation unc
   float muon_ISO_DY_unc_, electron_ID_DY_unc_;

   /// Histograms to store different SF
   TH2* h2_electron_ID_histo_;
   TH2* h2_electron_Reco_histo_;
   TH2* h2_muon_ID_histo_;
   TH2* h2_muon_ISO_histo_;
   
   /// Store shift direction
   float uncFactor_electron_RECO_;
   float uncFactor_electron_ID_;
   float uncFactor_muon_ISO_;
   float uncFactor_muon_ID_;
};




#endif /* LEPTONSF_HPP__ */
