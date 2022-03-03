#include "physics.hpp"

#include <limits>
#include <iostream>

std::vector<tree::Jet> phys::getCleanedJets(std::vector<tree::Jet> const &jets,TLorentzVector const &p_l1,TLorentzVector const &p_l2)
{
   std::vector<tree::Jet> cjets;
   for (const tree::Jet &j: jets){
      if (!j.TightIDlepVeto || j.p.Pt()<30 || fabs(j.p.Eta())>2.4) continue;
      
      // Check overlap with selected leptons
      if(j.p.DeltaR(p_l1)<0.4) continue;
      else if (j.p.DeltaR(p_l2)<0.4) continue;
      
      cjets.push_back(j);
   }
   sort(cjets.begin(), cjets.end(), tree::PtGreater);
   return cjets;
}

// ~std::vector<tree::Jet> phys::getCleanedJets_looseID(std::vector<tree::Jet> const &jets)
// ~{
   // ~std::vector<tree::Jet> cjets;
   // ~for (tree::Jet j: jets){
      // ~if (!j.isLoose || j.p.Pt()<30 || fabs(j.p.Eta())>2.4) continue;
      // ~if (j.hasElectronMatch || j.hasMuonMatch) continue;
      // ~cjets.push_back(j);
   // ~}
   // ~sort(cjets.begin(), cjets.end(), tree::PtGreater);
   // ~return cjets;
// ~}


/*
 * - compute HT using the given jet collection
 * - jets have to be cleaned of photons/electrons and
 *   quality criteria have to be applied before, if wanted
 * - any kinematic selection (pT, eta) has to be done before
 */
float phys::computeHT(std::vector<tree::Jet> const &jets)
{
   float HT=0.0;
   for (auto const &jet: jets) HT+=jet.p.Pt();
   return HT;
}

float phys::METoverSqrtHT(float MET, float HT)
{
   return (HT==0.0
           ? std::numeric_limits<float>::infinity()
           : MET/TMath::Sqrt(HT));
}

float phys::M_T_squared(TLorentzVector const &v1, TLorentzVector const &v2)
{
   return 2.0*v1.Pt()*v2.Pt()*(1-TMath::Cos(v1.DeltaPhi(v2)));
}

float phys::M_T(TLorentzVector const &v1, TLorentzVector const &v2)
{
   return TMath::Sqrt(M_T_squared(v1,v2));
}

float phys::M_T(tree::Particle const &p1, tree::Particle const &p2)
{
   return M_T(p1.p,p2.p);
}

float phys::conM_T(TLorentzVector const &v1, TLorentzVector const &v2)
{
   return TMath::Sqrt(2.0*v1.Pt()*v2.Pt()*(1+TMath::Cos(v1.DeltaPhi(v2))));
}

bool phys::matchesGen(tree::Particle const &p,std::vector<tree::GenParticle> const &genP,int pdgId,float dR,float rel_dp)
{
   for (tree::GenParticle const &gP: genP){
      float const rDP=fabs(p.p.Pt()-gP.p.Pt())/p.p.Pt();
      if (abs(gP.pdgId)==pdgId && gP.p.DeltaR(p.p)<dR && rDP<rel_dp) {
         return true;
      }
   }
   return false;
}

float phys::invmass(TVector3 const &v1, TVector3 const &v2)
{  
   float const p1 = v1.Mag();
   float const p2 = v2.Mag();
   TVector3 const pp = v1+v2;
   return TMath::Sqrt((p1+p2)*(p1+p2)-pp.Mag()*pp.Mag());;
}

float phys::invmass(tree::Particle const &p1, tree::Particle const &p2)
{
   return invmass(p1.p.Vect(),p2.p.Vect());
}

float phys::sumMlb(TLorentzVector &lepton1, TLorentzVector &lepton2, const std::vector<tree::Jet> &jets, const std::vector<tree::Jet> &bjets){
   float mlb_min = 1.E6;
   float mlb_max = 1.E6;

   float temp_mlb;
   float result_sum_mlb;

   TLorentzVector jet (0.,0.,0.,0.);
   TLorentzVector jet1 (0.,0.,0.,0.);
   TLorentzVector lepton (0.,0.,0.,0.);

   TLorentzVector leptList[2] = {lepton1,lepton2};

   int lmin = -1;

   std::vector< tree::Jet > jet1Coll;
   std::vector< tree::Jet > jet2Coll;

   // Calculate sum Mlb
   if (bjets.size() >= 1) jet1Coll = bjets;
   else jet1Coll = jets;

   if (bjets.size() >= 2) jet2Coll = bjets;
   else jet2Coll = jets;

   for (int il=0; il < 2; ++il){
      for (std::size_t ij=0; ij < jet1Coll.size(); ++ij){
         jet.SetPxPyPzE(jet1Coll.at(ij).p.Px(),jet1Coll.at(ij).p.Py(),jet1Coll.at(ij).p.Pz(),jet1Coll.at(ij).p.Energy());
         if (jet.Pt() > 30. && abs(jet.Eta()) < 2.4){
            lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
        
            temp_mlb = (jet+lepton).M();
            if (temp_mlb < mlb_min) {
               mlb_min = temp_mlb;
               lmin = il;
               jet1 = jet;
            }
         }
      }
   }
   for (int il=0; il < 2; ++il){
      if (il == lmin) continue;
      for (std::size_t ij=0; ij < jet2Coll.size(); ++ij){
         jet.SetPxPyPzE(jet2Coll.at(ij).p.Px(),jet2Coll.at(ij).p.Py(),jet2Coll.at(ij).p.Pz(),jet2Coll.at(ij).p.Energy());
         if (jet.Pt() > 30. && abs(jet.Eta()) < 2.4){
            if (bjets.size() == 1 && jet.DeltaR (jet1) < 0.1) continue;
            if ( (bjets.size() == 0 || bjets.size() >= 2) &&  jet == jet1) continue;  
            lepton.SetPxPyPzE(leptList[il].Px(),leptList[il].Py(),leptList[il].Pz(),leptList[il].Energy());
        
            temp_mlb = (jet+lepton).M();
            if (temp_mlb < mlb_max) {
               mlb_max = temp_mlb;
            }
         }
      }
   }
   
   if (mlb_min < 1.E6 and mlb_max < 1.E6) result_sum_mlb = mlb_min + mlb_max; 
   else result_sum_mlb = 1.E-6;
   
   return result_sum_mlb;
  
}

float phys::dPhi(float &a, float &b){
   float x=a-b;
   if (x >= M_PI) x -= 2*M_PI;
   else if (x < -M_PI) x += 2*M_PI;
   return x;
}

float phys::MT2(TLorentzVector const &v1, TLorentzVector const &v2, TLorentzVector const &met){
   MT2Functor fctMT2_;
   double pa[3] = {v1.M(),v1.Px(),v1.Py()};
   double pb[3] = {v2.M(),v2.Px(),v2.Py()};
   double pmiss[3] = {0.,met.Px(),met.Py()};
   
   fctMT2_.set_mn(0.);
   fctMT2_.set_momenta(pa,pb,pmiss);
   
   return fctMT2_.get_mt2();
}

TH1F phys::getSystShift(TH1F const &nominal,TH1F const &syst){
   
   TH1F shift = *(TH1F*)syst.Clone();   
   shift.Add(&nominal,-1.);
   
   return shift;
}

   
   

