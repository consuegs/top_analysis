// Class to derive events weights for jet pileup ID (partly copied from bTagWeights)
#ifndef JETPILEUPIDSF_HPP__
#define JETPILEUPIDSF_HPP__

#include "tree/TreeParticles.hpp"
#include "tools/io.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

class JetPileupIDWeights
{
   public:
      
      JetPileupIDWeights(const std::string& sf_fileName, const std::string& sf_histName, const std::string& eff_histName, const Systematic::Systematic& systematic);
      
      /// Get event weight
      const float getEventWeight(const std::vector<tree::Jet>& jets);
            
      ~JetPileupIDWeights(){};
   
   private:
      /// Current systematic
      const Systematic::Systematic systematic_;
      
      /// Enumeration for possible systematics
      enum SystematicInternal{nominal, vary_up, vary_down};
      SystematicInternal systematicInternal_;
      
      /// Return the histograms
      void prepareHists(const std::string& sf_fileName, const std::string& sf_histName, const std::string& eff_histName,
                        TH2* &sf_hist_, TH2* &sfUnc_hist_, TH2* &eff_hist_);
      
      /// Helper to derive Eff from individual jet
      const float getJetEff(const float& pt, const float& eta);
      
      /// Helper to derive SF from individual jet
      const float getJetSF(const float& pt, const float& eta);
      
      /// Histograms with SF, Eff and unc
      TH2* sf_hist_;
      TH2* sfUnc_hist_;
      TH2* eff_hist_;
   
};

#endif /* JETPILEUPIDSF_HPP__ */
