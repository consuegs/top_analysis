#ifndef BTAGWEIGHTS_HPP__
#define BTAGWEIGHTS_HPP__

#include "../ext/BTagCalibrationStandalone.h"
#include "tree/TreeParticles.hpp"
#include "tools/io.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

//////////////////////////////////// btag Eff //////////////////////////////////////////////////////
class BTagEffMapFunctor
{
   public:
      BTagEffMapFunctor(const std::string& bTagEffPath,const std::string& tagger, const BTagEntry::OperatingPoint& WP, const Systematic::Systematic& systematic);
      const float getEff(const int& flavor, const double& pt, const double& eta, const int& channel);
      enum Channel { ee, mumu, emu};
  
   private:
      void readBTagEffHistos();
      
      const std::string tagger_;
      const BTagEntry::OperatingPoint WP_;
      const std::string bTagEffPath_;
      bool nominalEff_;
      const Systematic::Systematic systematic_;
      
      TFile effFile_;
      TH2F bEfficiencies_ee_;
      TH2F cEfficiencies_ee_;
      TH2F lightEfficiencies_ee_;

      TH2F bEfficiencies_mumu_;
      TH2F cEfficiencies_mumu_;
      TH2F lightEfficiencies_mumu_;

      TH2F bEfficiencies_emu_;
      TH2F cEfficiencies_emu_;
      TH2F lightEfficiencies_emu_;
};

//////////////////////////////////// btag Weights //////////////////////////////////////////////////////
class BTagWeights
{
   public:
      BTagWeights(const std::string& bTagSFFile, const std::string& bTagEffPath, const std::string& tagger, const BTagEntry::OperatingPoint& WP, const float& WPcut, const Systematic::Systematic& systematic);
      
      /// Get bTag event weight
      const float getEventWeight(const std::vector<tree::Jet>& jets, const int& channel);
      
      /// Map to find correct bTag unc names
      static const std::map<TString, std::vector<std::string> > BTagSourcesMap;
   
   private:
      /// Current systematic
      const Systematic::Systematic systematic_;
      
      /// Helpers to read SF from files
      BTagCalibration bTagCalibration_;
      BTagCalibrationReader calibrationReader_bJets_;
      BTagCalibrationReader calibrationReader_cJets_;
      BTagCalibrationReader calibrationReader_lightJets_;
      
      /// Cut for current WP
      const float WPcut_;
      
      /// Bool if current tagger is DeepJet
      bool isDeepJet_;
      
      /// Helper to read eff from file
      BTagEffMapFunctor bTagEff_;
      
      /// Helper to derive SF from individual jet
      const float getJetSF(const int& flavor, const float& eta, const float& pt);
   
};

#endif /* BTAGWEIGHTS_HPP__ */
