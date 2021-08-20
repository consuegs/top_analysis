#ifndef JETCORRECTIONS_HPP__
#define JETCORRECTIONS_HPP__

#include "ext/JetCorrectorParameters.h"
#include "ext/JetCorrectionUncertainty.h"
#include "ext/JetResolution.h"
#include "tree/TreeParticles.hpp"
#include <vector>
#include <cmath>
#include <map>
#include <TRandom3.h>
#include "systematics.hpp"

class jesCorrections
{
   public:
      jesCorrections(const std::string&, const Systematic::Systematic&);
      
      /// Scale the jet collection
      void applySystematics(std::vector<tree::Jet>&, TLorentzVector&);
      
      /// Map to find correct JES unc names
      static const std::map<std::string, std::vector<std::string> > JESUncSourcesMap;
   
   private:
      /// Scale jet according to systematic variation
      void scaleJet(tree::Jet&);
      
      /// Enumeration for possible systematics
      enum SystematicInternal { nominal, vary_up, vary_down, undefined };
      
      /// Systematic variation
      SystematicInternal systematicInternal_;
      bool up=false;

      const Systematic::Systematic systematic_;
      
      /// Helpers to read scale from file
      JetCorrectionUncertainty* uncertainty_;
   
};

class jerCorrections
{
   public:
      jerCorrections(const std::string&, const std::string&, const Systematic::Systematic&);
      
      /// Smear the jet collection based on hybrid method
      void smearCollection_Hybrid(std::vector<tree::Jet>&, const float&);
   
   private:
      JME::JetResolution* m_resolution_;
      JME::JetResolutionScaleFactor* m_ScaleFactor_;
      
      float rho_;
      
      /// Smear jet based on hybrid method
      void smearJet_Hybrid(tree::Jet&);
      
      /// Variation the be applied
      Variation systematic_variation_;
      const Systematic::Systematic systematic_;
      
      /// Check if jet has to be taken into account for systematic
      bool checkApplySystematic(const TLorentzVector&);
      
      /// Random generator
      TRandom* randomGen() const;
      mutable TRandom* randomGen_;
   
      
   
};

#endif /* JETCORRECTIONS_HPP__ */
