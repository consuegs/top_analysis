// Class to apply jet veto maps
#ifndef JETVETOMAPS_HPP__
#define JETVETOMAPS_HPP__

#include "tree/TreeParticles.hpp"
#include "tools/io.hpp"
#include <vector>
#include <cmath>
#include <map>
#include "systematics.hpp"
#include "TH2F.h"
#include "TFile.h"

class JetVetoMaps
{
   public:
      
      JetVetoMaps(const std::string& map_fileName, const std::string& map_histName, const std::string& map_histName_MC16, const Systematic::Systematic& systematic, const bool is2016=false);
      
      /// Check Veto application (false if veto should be applied)
      const bool checkVetoMap(const std::vector<tree::Jet>& jets);
            
      ~JetVetoMaps(){};
   
   private:
      /// Current systematic
      const Systematic::Systematic systematic_;
      
      /// Read the histogram
      void prepareHists(const std::string& map_fileName, const std::string& map_histName, const std::string& map_histName_MC16, TH2* &vetoMap_hist_, TH2* &vetoMap_hist_MC16_);
      
      /// Bool to store if 2016 is running (has two maps)
      bool is2016_;
      
      /// Histogram with veto map
      TH2* vetoMap_hist_;
      TH2* vetoMap_hist_MC16_;
   
};

#endif /* JETVETOMAPS_HPP__ */
