#include "leptonCorrections.hpp"
#include <iostream>


leptonCorrections::leptonCorrections(const Systematic::Systematic& systematic) :
   systematic_(systematic)
{
   std::cout<<"--- Beginning preparation of Lepton corrections\n";
   
   // Set internal systematic
   systematicInternal_ = nominal;
   electronCorrection_ = NominalEle;
   muonCorrection_ = NominalMu;
   const Systematic::Type type = systematic.type();
   if(std::find(Systematic::leptonScaleResTypes.begin(), Systematic::leptonScaleResTypes.end(), type) != Systematic::leptonScaleResTypes.end()){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up){ systematicInternal_ = vary_up; }
      else if(systematic.variation() == Systematic::down){ systematicInternal_ = vary_down; }
      else {
         std::cerr << "ERROR in constructor of leptonCorrections! Systematic variation is invalid: "
                  << Systematic::convertVariation(systematic.variation()) << "\n...break\n\n";
         exit(98);
      }
      // Set which correction to be used
      if(Systematic::convertType(type).BeginsWith("ELECTRON")){
         switch(type){
            case Systematic::Type::eleScale:
               if (systematicInternal_ == vary_up) electronCorrection_ = ScaleUp;
               else electronCorrection_ = ScaleDown;
               break;
            case Systematic::Type::eleSmearingRho:
               if (systematicInternal_ == vary_up) electronCorrection_ = SigmaRhoUp;
               else electronCorrection_ = SigmaRhoDown;
               break;
            case Systematic::Type::eleSmearingPhi:
               if (systematicInternal_ == vary_up) electronCorrection_ = SigmaPhiUp;
               else electronCorrection_ = SigmaPhiDown;
               break;
            case Systematic::Type::eleScaleSmearing:
               isTotalEle_ = true;
               break;
         }
      }
      else if (Systematic::convertType(type).BeginsWith("MUON")){
         switch(type){
            case Systematic::Type::muonScaleStat:
               muonCorrection_ = Stat_RMS;
               break;
            case Systematic::Type::muonScaleZpt:
               muonCorrection_ = Zpt;
               break;
            case Systematic::Type::muonScaleEwk:
               muonCorrection_ = Ewk;
               break;
            case Systematic::Type::muonScaleDeltaM:
               muonCorrection_ = DeltaM;
               break;
            case Systematic::Type::muonScaleEwk2:
               muonCorrection_ = Ewk2;
               break;
            case Systematic::Type::muonScale:
               isTotalMu_ = true;
               break;
         }
         if (systematicInternal_ == vary_up) uncFactor_muon_ = 1.;      // factor for shift direction
         else uncFactor_muon_ = -1.;
      }
   }
   
   
   // Print which systematic is used
   if(systematicInternal_ == vary_up) {std::cout<<"Apply systematic variation: up\n";}
   else if(systematicInternal_ == vary_down) std::cout<<"Apply systematic variation: down\n";
   else std::cout<<"Do not apply systematic variation\n";

}

std::vector<tree::Electron> leptonCorrections::correctElectrons(std::vector<tree::Electron>& electrons, TLorentzVector& met, 
                                                      TLorentzVector& metPuppi, const bool applySCcut)
{
   for(size_t iEle=0; iEle<electrons.size(); ++iEle){   //Remove electrons, which will be corrected from MET
      TLorentzVector temp(0.,0.,0.,0.);
      temp.SetPtEtaPhiM(electrons.at(iEle).p.Pt(),0.,electrons.at(iEle).p.Phi(),0.);
      met = met + temp;
      metPuppi = metPuppi + temp;
   }
   
   std::vector<tree::Electron> cElectrons;
   for(size_t iEle=0; iEle<electrons.size(); ++iEle){
      this->correctEletron(electrons.at(iEle));
      if(electrons.at(iEle).p.Pt()>20 && (!applySCcut || electrons.at(iEle).etaSC<2.4)){     // Check only pT and maybe SC cut. Other cuts applied on tree level
         cElectrons.push_back(electrons.at(iEle));
      }
   }
   
   for(size_t iEle=0; iEle<electrons.size(); ++iEle){   //Add shifted electrons back to met
      TLorentzVector temp(0.,0.,0.,0.);
      temp.SetPtEtaPhiM(electrons.at(iEle).p.Pt(),0.,electrons.at(iEle).p.Phi(),0.);
      met = met - temp;
      metPuppi = metPuppi - temp;
   }
   sort(cElectrons.begin(), cElectrons.end(), tree::PtGreater);
   return cElectrons;
}

std::vector<tree::Muon> leptonCorrections::correctMuons(std::vector<tree::Muon>& muons, TLorentzVector& met, TLorentzVector& metPuppi)
{
   for(size_t iMu=0; iMu<muons.size(); ++iMu){   //Remove muons, which will be corrected from MET
      TLorentzVector temp(0.,0.,0.,0.);
      temp.SetPtEtaPhiM(muons.at(iMu).p.Pt(),0.,muons.at(iMu).p.Phi(),0.);
      met = met + temp;
      metPuppi = metPuppi + temp;
   }
   
   std::vector<tree::Muon> cMuons;
   for(size_t iMu=0; iMu<muons.size(); ++iMu){
      this->correctMuon(muons.at(iMu));
      if(muons.at(iMu).p.Pt()>20){     // Check only pT. Other cuts applied on tree level
         cMuons.push_back(muons.at(iMu));
      }
   }
   
   for(size_t iMu=0; iMu<muons.size(); ++iMu){   //Add shifted muons back to met
      TLorentzVector temp(0.,0.,0.,0.);
      temp.SetPtEtaPhiM(muons.at(iMu).p.Pt(),0.,muons.at(iMu).p.Phi(),0.);
      met = met - temp;
      metPuppi = metPuppi - temp;
   }
   sort(cMuons.begin(), cMuons.end(), tree::PtGreater);
   return cMuons;
}

void leptonCorrections::correctEletron(tree::Electron& ele)
{
   if(!isTotalEle_) ele.p *= ele.corrections[electronCorrection_];
   else{ // envelope for scale and smearing
      float nominalCorr = ele.corrections[NominalEle];
      float totalCorr = 0.;
      if (systematicInternal_ == vary_up) {
         totalCorr = nominalCorr+TMath::Sqrt(pow(nominalCorr-ele.corrections[ScaleUp],2)+
                                          pow(nominalCorr-ele.corrections[SigmaRhoUp],2)+
                                          pow(nominalCorr-ele.corrections[SigmaPhiUp],2));
      }
      else {
         totalCorr = nominalCorr-TMath::Sqrt(pow(nominalCorr-ele.corrections[ScaleDown],2)+
                                          // ~pow(nominalCorr-ele.corrections[SigmaPhiDown],2)+    //only one template for SigmaPhi (see https://twiki.cern.ch/twiki/bin/view/CMS/EgammaMiniAODV2#Applying_the_Energy_Scale_and_sm)
                                          pow(nominalCorr-ele.corrections[SigmaRhoDown],2));
      }
      ele.p *= totalCorr;
   }
}

void leptonCorrections::correctMuon(tree::Muon& mu)
{
   if(!isTotalMu_) mu.p *= mu.corrections[NominalMu]+uncFactor_muon_*abs(mu.corrections[NominalMu]-mu.corrections[muonCorrection_]);
   else{ // envelope for scale and smearing
      float nominalCorr = mu.corrections[NominalMu];
      float totalCorr = nominalCorr+uncFactor_muon_*TMath::Sqrt(pow(nominalCorr-mu.corrections[Stat_RMS],2)+
                                                               pow(nominalCorr-mu.corrections[Zpt],2)+
                                                               pow(nominalCorr-mu.corrections[Ewk],2)+
                                                               pow(nominalCorr-mu.corrections[DeltaM],2)+
                                                               pow(nominalCorr-mu.corrections[Ewk2],2));
      mu.p *= totalCorr;
   }
}
