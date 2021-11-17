#ifndef TRIGGERSF_HPP__
#define TRIGGERSF_HPP__

#include "tools/io.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

/// Class for application of trigger scale factors
class TriggerScaleFactors{

public:

   /// Constructor
   TriggerScaleFactors(const std::string& fileName,
                  const Systematic::Systematic& systematic);

   /// Destructor
   ~TriggerScaleFactors(){}



   /// Get trigger scale factor for given leptons pTs
   const float getTriggerSF(const float &pT_l1, const float &pT_l2,const int &channel, const bool &muonlead);


private:


   /// Enumeration for possible systematics
   enum SystematicInternal{nominal, vary_up, vary_down};
   SystematicInternal systematicInternal_;
   
   const Systematic::Systematic systematic_;

   /// Return the scale factor histogram
   TH2* prepareSF(const std::string& inputFileName,
                   const std::string& histogramName)const;

   
   /// Get SF from 2D hist
   const float get2DSF(const TH2* const histo, const float& x, const float& y);

   /// Histograms to store different SF
   TH2* h2_ee_;
   TH2* h2_mumu_;
   TH2* h2_emu_;
};




#endif /* TRIGGERSF_HPP__ */
