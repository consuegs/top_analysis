#ifndef PHYSICS_HPP__
#define PHYSICS_HPP__

#include "tree/TreeParticles.hpp"

namespace phys
{
   std::vector<tree::Jet> getCleanedJets(std::vector<tree::Jet> const &jets);
   std::vector<tree::Jet> getCleanedJets_looseID(std::vector<tree::Jet> const &jets);
   float computeHT(std::vector<tree::Jet> const &jets);
   float METoverSqrtHT(float MET, float HT); // returns inf for HT=0.0

   // transverse mass (for massless daughters)
   float M_T (tree::Particle const &p1, tree::Particle const &p2);
   float M_T (TLorentzVector const &v1, TLorentzVector const &v2);
   float M_T2(TLorentzVector const &v1, TLorentzVector const &v2);
   
   // contransverse mass (for massless daugthers)
   float conM_T(TLorentzVector const &v1, TLorentzVector const &v2);

   // matching
   bool matchesGen(tree::Particle const &p,std::vector<tree::GenParticle> const &genP,int pdgId,float dR,float rel_dp);
   
   // invariant mass
   float invmass(tree::Particle const &p1, tree::Particle const &p2);
   float invmass(TVector3 const &v1, TVector3 const &v2);
   
   // sum mlb
   float sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<tree::Jet> &jets, const std::vector<tree::Jet> &bjets);
} // namespace phys

#endif /* PHYSICS_HPP__ */
