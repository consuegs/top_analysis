//Partly copied from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/common/src/ScaleFactors.cc
#include "jetCorrections.hpp"
#include "Config.hpp"
#include <iostream>

//////////////////////////////////// Jet Energy Scales //////////////////////////////////////////////////////

jesCorrections::jesCorrections(const std::string& jesUncertaintySourceFile, const Systematic::Systematic& systematic) :
   systematicInternal_(undefined),
   systematic_(systematic)
{
   std::cout << "--- Beginning preparation of JES Corrections\n";
   // Set internal systematic
   systematicInternal_ = nominal;
   const Systematic::Type type = systematic.type();
   if(std::find(Systematic::jesTypes.begin(), Systematic::jesTypes.end(), type) != Systematic::jesTypes.end()){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up) systematicInternal_ = vary_up;
      else if(systematic.variation() == Systematic::down) systematicInternal_ = vary_down;
      else{
         std::cerr<<"ERROR in constructor of jesCorrections! Systematic variation is invalid: "
              <<Systematic::convertVariation(systematic.variation())<<"\n...break\n"<<std::endl;
         exit(98);
      }
      
      // Set restricted flavors
      restricttoflav_.clear();
      useRestrictFlav_ = false;
      if(std::find(Systematic::jesTypes_pureFlavor.begin(), Systematic::jesTypes_pureFlavor.end(), type) != Systematic::jesTypes_pureFlavor.end()){
         useRestrictFlav_ = true;
         switch(type){
            case Systematic::Type::jesFlavorPureGluon:
               restricttoflav_.push_back(21);
               restricttoflav_.push_back(0);
               std::cout<<"Restrict flavor to 21, 0"<<std::endl;
               break;
            case Systematic::Type::jesFlavorPureQuark:
               restricttoflav_.push_back(1);
               restricttoflav_.push_back(2);
               restricttoflav_.push_back(3);
               std::cout<<"Restrict flavor to 1, 2, 3"<<std::endl;
               break;
            case Systematic::Type::jesFlavorPureCharm:
               restricttoflav_.push_back(4);
               std::cout<<"Restrict flavor to 4"<<std::endl;
               break;
            case Systematic::Type::jesFlavorPureBottom:
               restricttoflav_.push_back(5);
               std::cout<<"Restrict flavor to 5"<<std::endl;
               break;
         }
      }
   }
   // Print which systematic is used and set boolean
   if(systematicInternal_ == vary_up) {std::cout<<"Apply systematic variation: up\n"; up=true;}
   else if(systematicInternal_ == vary_down) std::cout<<"Apply systematic variation: down\n";
   else std::cout<<"Do not apply systematic variation\n";
   
   // Check existence of uncertainty file
   if(jesUncertaintySourceFile.empty()){
   std::cerr<<"ERROR in constructor of jesCorrections! "
      <<"Systematic variation requested, but no uncertainty file specified\n...break\n"<<std::endl;
   exit(98);
   }
   
   // Configure helper
   useRealisticFlav_ = false;
   if(systematicInternal_!=nominal){
      if (type != Systematic::Type::jesFlavorRealistic){
         try {
            uncertainty_ = new JetCorrectionUncertainty(JetCorrectorParameters(jesUncertaintySourceFile,JESUncSourcesMap.at(systematic_.type_str())[0]));
         }
         catch(std::runtime_error &rte){
            std::cout << "JESBase::setFile: Uncertainty for source " << convertType(type) << " not found! Skipping\n";
            exit(98);
         }
      }
      else{    // for realistic flavor mixing, four different unc. sources are needed
         useRealisticFlav_ = true;
         std::cout<<"Use realistic flavor mixing"<<std::endl;
         try {
            uncertaintiesRelFlav_.push_back(new JetCorrectionUncertainty(JetCorrectorParameters(jesUncertaintySourceFile,JESUncSourcesMap.at(systematic_.type_str())[0])));
            uncertaintiesRelFlav_.push_back(new JetCorrectionUncertainty(JetCorrectorParameters(jesUncertaintySourceFile,JESUncSourcesMap.at(systematic_.type_str())[1])));
            uncertaintiesRelFlav_.push_back(new JetCorrectionUncertainty(JetCorrectorParameters(jesUncertaintySourceFile,JESUncSourcesMap.at(systematic_.type_str())[2])));
            uncertaintiesRelFlav_.push_back(new JetCorrectionUncertainty(JetCorrectorParameters(jesUncertaintySourceFile,JESUncSourcesMap.at(systematic_.type_str())[3])));
         }
         catch(std::runtime_error &rte){
            std::cout << "JESBase::setFile: Uncertainty for source " << convertType(type) << " not found! Skipping\n";
            exit(98);
         }
      }
   }
   
}
   

// ~void jesCorrections::applySystematics(std::vector<tree::Jet>& Jets, TLorentzVector& MET)
void jesCorrections::applySystematics(std::vector<tree::Jet>& Jets, std::vector<tree::MET*>& METs)
{
   if(systematicInternal_!=nominal){ // Apply shifts only if JES systematic is selected
      for(size_t iJet=0; iJet<Jets.size(); ++iJet){   //Remove jets, which will be corrected from MET
         if ((Jets.at(iJet).cef+Jets.at(iJet).nef)>0.9 || Jets.at(iJet).hasMuonMatch_loose) continue;
         TLorentzVector temp(0.,0.,0.,0.);
         temp.SetPtEtaPhiM(Jets.at(iJet).p.Pt(),0.,Jets.at(iJet).p.Phi(),0.);
         for (tree::MET* MET : METs){
            MET->p = MET->p + (temp-(temp*(Jets.at(iJet).uncorJecFactor_L1)));
         }
      }
      
      // This loop corrects the jet collection used for jet selections
      for(size_t iJet=0; iJet<Jets.size(); ++iJet){
         if (abs(Jets.at(iJet).p.Eta())>=5.4) continue;
         this->scaleJet(Jets.at(iJet));
      }
      
      for(size_t iJet=0; iJet<Jets.size(); ++iJet){   //Add shifted jets to MET
         if ((Jets.at(iJet).cef+Jets.at(iJet).nef)>0.9 || Jets.at(iJet).hasMuonMatch_loose) continue;
         TLorentzVector temp(0.,0.,0.,0.);
         temp.SetPtEtaPhiM(Jets.at(iJet).p.Pt(),0.,Jets.at(iJet).p.Phi(),0.);
         for (tree::MET* MET : METs){
            MET->p = MET->p - (temp-(temp*(Jets.at(iJet).uncorJecFactor_L1)));
         }
      }
   }

}

void jesCorrections::scaleJet(tree::Jet& jet)
{
   // Check jet matches restricted flavor
   if (useRestrictFlav_) {
      if(std::find(restricttoflav_.begin(),restricttoflav_.end(),abs(jet.hadronFlavour)) == restricttoflav_.end()){
         return;
      }
   }
   // Check if realistic flavor should be applied
   double uncert = 0.;
   if (!useRealisticFlav_){
      uncertainty_->setJetPt(jet.p.Pt());
      uncertainty_->setJetEta(jet.p.Eta());
      uncert = uncertainty_->getUncertainty(up);
   }
   else {
      switch(abs(jet.hadronFlavour)){
         case 0:  //Gluon
            uncertaintiesRelFlav_[0]->setJetPt(jet.p.Pt());
            uncertaintiesRelFlav_[0]->setJetEta(jet.p.Eta());
            uncert = uncertaintiesRelFlav_[0]->getUncertainty(up);
            break;
         case 4:  //Charm
            uncertaintiesRelFlav_[2]->setJetPt(jet.p.Pt());
            uncertaintiesRelFlav_[2]->setJetEta(jet.p.Eta());
            uncert = uncertaintiesRelFlav_[2]->getUncertainty(up);
            break;
         case 5:  //Bottom
            uncertaintiesRelFlav_[3]->setJetPt(jet.p.Pt());
            uncertaintiesRelFlav_[3]->setJetEta(jet.p.Eta());
            uncert = uncertaintiesRelFlav_[3]->getUncertainty(up);
            break;
         default: //Light Quark
            uncertaintiesRelFlav_[1]->setJetPt(jet.p.Pt());
            uncertaintiesRelFlav_[1]->setJetEta(jet.p.Eta());
            uncert = uncertaintiesRelFlav_[1]->getUncertainty(up);
            break;
      }
   }
         
   if(up){
      jet.p *= 1+uncert;
      jet.uncorJecFactor_L1 *= 1./(1+uncert);
   }
   else{
      jet.p *= 1-uncert;
      jet.uncorJecFactor_L1 *= 1./(1-uncert);
   }
}


const std::map<std::string, std::vector<std::string> > jesCorrections::JESUncSourcesMap = {

  {"JESAbsoluteStat",{"AbsoluteStat"}},
  {"JESAbsoluteScale",{"AbsoluteScale"}},
  {"JESAbsoluteFlavMap",{"AbsoluteFlavMap"}},
  {"JESAbsoluteMPFBias",{"AbsoluteMPFBias"}},
  {"JESFragmentation",{"Fragmentation"}},
  {"JESSinglePionECAL",{"SinglePionECAL"}},
  {"JESSinglePionHCAL",{"SinglePionHCAL"}},
  {"JESFlavorQCD",{"FlavorQCD"}},
  {"JESTimePtEta",{"TimePtEta"}},
  {"JESRelativeJEREC1",{"RelativeJEREC1"}},
  {"JESRelativeJEREC2",{"RelativeJEREC2"}},
  {"JESRelativeJERHF",{"RelativeJERHF"}},
  {"JESRelativePtBB",{"RelativePtBB"}},
  {"JESRelativePtEC1",{"RelativePtEC1"}},
  {"JESRelativePtEC2",{"RelativePtEC2"}},
  {"JESRelativePtHF",{"RelativePtHF"}},
  {"JESRelativeFSR",{"RelativeFSR"}},
  {"JESRelativeStatFSR",{"RelativeStatFSR"}},
  {"JESRelativeStatEC",{"RelativeStatEC"}},
  {"JESRelativeStatHF",{"RelativeStatHF"}},
  {"JESPileUpDataMC",{"PileUpDataMC"}},
  {"JESPileUpPtRef",{"PileUpPtRef"}},
  {"JESPileUpPtBB",{"PileUpPtBB"}},
  {"JESPileUpPtEC1",{"PileUpPtEC1"}},
  {"JESPileUpPtEC2",{"PileUpPtEC2"}},
  {"JESPileUpPtHF",{"PileUpPtHF"}},
  {"JESPileUpMuZero",{"PileUpMuZero"}},
  {"JESPileUpEnvelope",{"PileUpEnvelope"}},
  {"JESSubTotalPileUp",{"SubTotalPileUp"}},
  {"JESSubTotalRelative",{"SubTotalRelative"}},
  {"JESSubTotalPt",{"SubTotalPt"}},
  {"JESSubTotalScale",{"SubTotalScale"}},
  {"JESSubTotalMC",{"SubTotalMC"}},
  {"JESSubTotalAbsolute",{"SubTotalAbsolute"}},
  {"JESTotal",{"Total"}},
  {"JESTotalNoFlavor",{"TotalNoFlavor"}},
  {"JESFlavorZJet",{"FlavorZJet"}},
  {"JESFlavorPhotonJet",{"FlavorPhotonJet"}},
  {"JESFlavorPureGluon",{"FlavorPureGluon"}},
  {"JESFlavorPureQuark",{"FlavorPureQuark"}},
  {"JESFlavorPureCharm",{"FlavorPureCharm"}},
  {"JESFlavorPureBottom",{"FlavorPureBottom"}},
  {"JESFlavorRealistic",{"FlavorPureGluon","FlavorPureQuark","FlavorPureCharm","FlavorPureBottom"}},
  {"JESCorrelationGroupbJES",{"CorrelationGroupbJES"}},
  {"JESCorrelationGroupFlavor",{"CorrelationGroupFlavor"}},
  {"JESCorrelationGroupUncorrelated",{"CorrelationGroupUncorrelated"}},
  {"JESRelativeBal",{"RelativeBal"}},
  {"JESRelativeSample",{"RelativeSample"}},
  {"JESCorrelationGroupMPFInSitu",{"CorrelationGroupMPFInSitu"}},
  {"JESCorrelationGroupIntercalibration",{"CorrelationGroupIntercalibration"}},

  // JME combined sources
  {"JESAbsolute", {"Absolute"}},
  {"JESAbsoluteYear", {"AbsoluteYear"}},
  {"JESBBEC1", {"BBEC1"}},
  {"JESBBEC1Year", {"BBEC1Year"}},
  {"JESEC2", {"EC2"}},
  {"JESEC2Year", {"EC2Year"}},
  {"JESHF", {"HF"}},
  {"JESHFYear", {"HFYear"}},
};

//////////////////////////////////// Jet Energy Resolution //////////////////////////////////////////////////////



jerCorrections::jerCorrections(const std::string& jerSFSourceFile, const std::string& jerRESSourceFile, const Systematic::Systematic& systematic) :
   systematic_(systematic)
{
   std::cout<<"--- Beginning preparation of JER smearing\n";
   
   // Check existence of input files
   if(jerSFSourceFile.empty()){
   std::cerr<<"ERROR in constructor of jerCorrections! "
      <<"No ScaleFactor file specified\n...break\n"<<std::endl;
   exit(98);
   }
   if(jerRESSourceFile.empty()){
   std::cerr<<"ERROR in constructor of jerCorrections! "
      <<"No Resolution file specified\n...break\n"<<std::endl;
   exit(98);
   }
   
   m_resolution_ = new JME::JetResolution(jerRESSourceFile);
   m_ScaleFactor_ = new JME::JetResolutionScaleFactor(jerSFSourceFile);
   
   randomGen_ = new TRandom3();
   
   // Set internal systematic
   systematic_variation_ = Variation::NOMINAL;
   const Systematic::Type type = systematic.type();
   if(std::find(Systematic::jerTypes.begin(), Systematic::jerTypes.end(), type) != Systematic::jerTypes.end()){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up){ systematic_variation_ = Variation::UP; }
      else if(systematic.variation() == Systematic::down){ systematic_variation_ = Variation::DOWN; }
      else {
         std::cerr << "ERROR in constructor of jerCorrections! Systematic variation is invalid: "
                  << Systematic::convertVariation(systematic.variation()) << "\n...break\n\n";
         exit(98);
      }
   }
   
   // Print which systematic is used and set boolean
   if(systematic_variation_ == Variation::UP) {std::cout<<"Apply systematic variation: up\n";}
   else if(systematic_variation_ == Variation::DOWN) std::cout<<"Apply systematic variation: down\n";
   else std::cout<<"Do not apply systematic variation\n";

}

void jerCorrections::smearCollection_Hybrid(std::vector<tree::Jet>& Jets, const float& rho)
{
   // Set rho for smearing
   rho_ = rho;
   
   // This loop smeares the jet collection used for jet selections (for split JER check eta,pT bin)
   if(systematic_variation_==Variation::NOMINAL || systematic_.type()==Systematic::jer){
      for(size_t iJet=0; iJet<Jets.size(); ++iJet){
         this->smearJet_Hybrid(Jets.at(iJet));
      }
   }
   else{
      for(size_t iJet=0; iJet<Jets.size(); ++iJet){
         if(this->checkApplySystematic(Jets.at(iJet).p)) this->smearJet_Hybrid(Jets.at(iJet));
      }
   }
}

void jerCorrections::smearJet_Hybrid(tree::Jet& jet)
{
   double jer_sf= m_ScaleFactor_->getScaleFactor({{JME::Binning::JetEta, jet.p.Eta()}}, systematic_variation_);
   double jet_resolution= m_resolution_->getResolution({{JME::Binning::JetPt, jet.p.Pt()}, {JME::Binning::JetEta, jet.p.Eta()}, {JME::Binning::Rho, rho_}});
   float smearingFactor = 1.;
   
   const float radius = 0.4;
   const bool valid_genjet = (jet.matchedGenJet.Pt() > 0.);
   const bool match_DR = (jet.p.DeltaR(jet.matchedGenJet) < (0.5 * radius));
   const bool match_Dpt = (fabs(jet.p.Pt() - jet.matchedGenJet.Pt()) < (3. * jet_resolution * jet.p.Pt()));
   
   // ~const bool valid_genjet = true;
   // ~const bool match_DR = true;
   // ~const bool match_Dpt = true;
   
   // ~std::cout<<"Set genJET and seed manually for debugging"<<std::endl;
   // ~jet.seed = 1000;
   // ~jet.matchedGenJet.SetPtEtaPhiM(50.,1.2,2.5,100.);
   
   if(valid_genjet && match_DR && match_Dpt){
      smearingFactor = std::max(0., (1. + ((jer_sf - 1.) * ((jet.p.Pt() - jet.matchedGenJet.Pt()) / jet.p.Pt()))));
   }
   else {
      this->randomGen()->SetSeed(jet.seed);

      const float random_gauss = this->randomGen()->Gaus(0., jet_resolution);

      smearingFactor = std::max(0., (1. + (random_gauss * sqrt(std::max(((jer_sf*jer_sf) - 1.), 0.)))));
   }
   
   jet.p *= smearingFactor;
}

TRandom* jerCorrections::randomGen() const
{
  if(randomGen_ == nullptr)
  {
    throw std::runtime_error("JetEnergyResolutionScaleFactors::randomGen() -- null pointer to TRandom object");
  }

  return randomGen_;
}

bool jerCorrections::checkApplySystematic(const TLorentzVector& p)
{
   if(systematic_.type() == Systematic::jerEta0){
      return (std::abs(p.Eta()) < 1.93);
	}
	else if(systematic_.type() == Systematic::jerEta1){
	  return ((std::abs(p.Eta()) >= 1.93) and (std::abs(p.Eta()) < 2.5));
	}
	else if(systematic_.type() == Systematic::jerEta2Pt0){
	  return ((std::abs(p.Eta()) >= 2.5) and (std::abs(p.Eta()) < 3.0) and (p.Pt() < 50.));
	}
	else if(systematic_.type() == Systematic::jerEta2Pt1){
	  return ((std::abs(p.Eta()) >= 2.5) and (std::abs(p.Eta()) < 3.0) and (p.Pt() >= 50.));
	}
	else if(systematic_.type() == Systematic::jerEta3Pt0){
	  return ((std::abs(p.Eta()) >= 3.0) and (std::abs(p.Eta()) < 5.0) and (p.Pt() < 50.));
	}
	else if(systematic_.type() == Systematic::jerEta3Pt1){
	  return ((std::abs(p.Eta()) >= 3.0) and (std::abs(p.Eta()) < 5.0) and (p.Pt() >= 50.));
	}
   return false;
}
