#ifndef MCWEIGHTS_HPP__
#define MCWEIGHTS_HPP__

#include "tree/TreeParticles.hpp"
#include "tools/io.hpp"
#include "Dataset.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

/// Class for application of mcWeights
class mcWeights{

public:

   /// Constructor
   mcWeights(const Systematic::Systematic& systematic, const bool &isMadGraph, const bool &isPythiaOnly=false);

   /// Destructor
   ~mcWeights(){}
   
   /// Prepare lumi weight by seeting dataset und lumi
   void prepareLumiWeight(const Datasubset &dss, const float &lumi);

   /// Get per-event mcWeight
   const float getMCweight(const float &nominalWeight, const std::vector<float> &w_pdf, const std::vector<float> &w_ps, const std::vector<float> &w_bFrag, const std::vector<float> &w_match);
   
   /// Get per-event mcWeight including lumi weight
   const float getMCweight_lumiWeighted(const float &nominalWeight, const std::vector<float> &w_pdf, const std::vector<float> &w_ps, const std::vector<float> &w_bFrag, const std::vector<float> &w_match);


private:
   
   const Systematic::Systematic systematic_;
      
   bool useNominal;
   bool upVariation;
   bool isMadGraph_;
   bool isPythiaOnly_;
   
   bool preparedLumiWeight;
   float lumiWeight = 1.;
   
   std::vector<int> meScaleBins;

};

#endif /* MCWEIGHTS_HPP__ */
