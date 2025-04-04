#ifndef PHYSICS_HPP__
#define PHYSICS_HPP__

#include "tree/TreeParticles.hpp"
#include "tools/systematics.hpp"
#include "tools/MT2Functor.h"
#include "tools/jetCorrections.hpp"
#include <TH1F.h>

namespace phys
{
   std::vector<tree::Jet> getCleanedJets(std::vector<tree::Jet> &jets,TLorentzVector const &p_l1,TLorentzVector const &p_l2,jerCorrections &jerCorrector,const float& rho,std::vector<tree::MET*>& METs, const bool& applyPileupID=false);
   std::vector<tree::Jet> getCleanedJets_looseID(std::vector<tree::Jet> const &jets);
   float computeHT(std::vector<tree::Jet> const &jets);
   float METoverSqrtHT(float MET, float HT); // returns inf for HT=0.0

   // transverse mass (for massless daughters)
   float M_T (tree::Particle const &p1, tree::Particle const &p2);
   float M_T (TLorentzVector const &v1, TLorentzVector const &v2);
   float M_T_squared(TLorentzVector const &v1, TLorentzVector const &v2);
   
   // contransverse mass (for massless daugthers)
   float conM_T(TLorentzVector const &v1, TLorentzVector const &v2);

   // matching
   bool matchesGen(tree::Particle const &p,std::vector<tree::GenParticle> const &genP,int pdgId,float dR,float rel_dp);
   
   // invariant mass
   float invmass(tree::Particle const &p1, tree::Particle const &p2);
   float invmass(TVector3 const &v1, TVector3 const &v2);
   
   // sum mlb
   float sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<tree::Jet> &jets, const std::vector<tree::Jet> &bjets);
  
   // sum dPhi
   float dPhi(float &a, float &b);
   
   // MT2
   float MT2(TLorentzVector const &v1, TLorentzVector const &v2, TLorentzVector const &met);
   
   // systematics
   TH1F getSystShift(TH1F const &nominal, TH1F const &syst);
   
} // namespace phys

#endif /* PHYSICS_HPP__ */
