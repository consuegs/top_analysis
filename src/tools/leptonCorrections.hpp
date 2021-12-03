#ifndef LEPTONCORRECTIONS_HPP__
#define LEPTONCORRECTIONS_HPP__

#include "tree/TreeParticles.hpp"
#include <vector>
#include "systematics.hpp"

class leptonCorrections
{
   public:
      leptonCorrections(const Systematic::Systematic& systematic);
      
      /// Correct leptons based opn selected systematic
      std::vector<tree::Electron> correctElectrons(std::vector<tree::Electron>& electrons, TLorentzVector& met, 
                                                      TLorentzVector& metPuppi, const bool applySCcut = false);
      std::vector<tree::Muon> correctMuons(std::vector<tree::Muon>& muons, TLorentzVector& met, TLorentzVector& metPuppi);
   
   private:
      
      /// Correct individual lepton
      void correctEletron(tree::Electron& ele);
      void correctMuon(tree::Muon& mu);
      
      /// Sytematic the be applied
      const Systematic::Systematic systematic_;
      
      /// Shift direction to be applied
      enum SystematicInternal { nominal, vary_up, vary_down, undefined };
      SystematicInternal systematicInternal_;
      
      /// Enum for accessing correct shift
      enum ElectronCorrection { NominalEle, ScaleUp, ScaleDown, SigmaRhoUp, SigmaRhoDown, SigmaPhiUp, SigmaPhiDown};
      ElectronCorrection electronCorrection_;
      bool isTotalEle_ = false;
      
      enum MuonCorrection { NominalMu, Stat_RMS, Zpt, Ewk, DeltaM, Ewk2};
      MuonCorrection muonCorrection_;
      bool isTotalMu_ = false;
      
      /// Store shift direction for muons
      float uncFactor_muon_;
   
};

#endif /* LEPTONCORRECTIONS_HPP__ */
