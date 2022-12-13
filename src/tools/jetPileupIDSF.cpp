#include "jetPileupIDSF.hpp"
#include "Config.hpp"
#include <iostream>

JetPileupIDWeights::JetPileupIDWeights(const std::string& sf_fileName, const std::string& sf_histName, const std::string& eff_histName, const Systematic::Systematic& systematic) :
   systematic_(systematic),
   sf_hist_(0),
   sfUnc_hist_(0),
   eff_hist_(0)
{
   std::cout << "--- Beginning preparation of jet pileup ID weights\n";
   // Check existence of input files
   if(!io::fileExists(sf_fileName)){
   std::cerr<<"ERROR in constructor of JetPileupIDWeights! "
      <<"ScaleFactor file does not exist\n...break\n"<<std::endl;
   exit(200);
   }
   
   // Check selected systematic
   const Systematic::Type type = systematic.type();
   systematicInternal_ = nominal;
   if(type == Systematic::jetPileupID){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up){
         systematicInternal_ = vary_up; 
         std::cout<<"Apply systematic variation: up\n";
      }
      else if(systematic.variation() == Systematic::down){
         systematicInternal_ = vary_down;
         std::cout<<"Apply systematic variation: down\n";
      }
      else {
         std::cerr << "ERROR in constructor of JetPileupIDWeights! Systematic variation is invalid: "
                  << Systematic::convertVariation(systematic.variation()) << "\n...break\n\n";
         exit(201);
      }
   }
   else std::cout<<"Do not apply systematic variation\n";
      
   // Read required histograms
   this->prepareHists(sf_fileName,sf_histName,eff_histName,sf_hist_,sfUnc_hist_,eff_hist_);
   
}

void JetPileupIDWeights::prepareHists(const std::string& sf_fileName, const std::string& sf_histName, const std::string& eff_histName,
                        TH2* &sf_hist_, TH2* &sfUnc_hist_, TH2* &eff_hist_)
{
   
   io::RootFileReader sf_file(sf_fileName.c_str(),"",false);

   // Access histograms
   sf_hist_ = (TH2*)(sf_file.read<TH2F>(sf_histName));
   sfUnc_hist_ = (TH2*)(sf_file.read<TH2F>(sf_histName+"_Systuncty"));
   eff_hist_ = (TH2*)(sf_file.read<TH2F>(eff_histName));
   if(!sf_hist_ || !sfUnc_hist_ || !eff_hist_){
      std::cerr<<"Error in JetPileupIDWeights::prepareHists()! Not all TH2 found in: "<<sf_fileName
             <<"\n...break\n"<<std::endl;
      exit(202);
   }

   // Store histogram in memory and close file
   sf_hist_->SetDirectory(0);
   sfUnc_hist_->SetDirectory(0);    //stored as absolute shifts
   eff_hist_->SetDirectory(0);
}

const float JetPileupIDWeights::getEventWeight(const std::vector<tree::Jet>& jets)
{
   float P_MC = 1.;
   float P_Data = 1.;
   
   // Derive Scale Factor following 1a from bTagWeights (https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods)
   // For UL no mistagging SF available, thus loop runs only over matched Jets
   for(size_t iJet=0; iJet<jets.size(); ++iJet){
      if (jets.at(iJet).p.DeltaR(jets.at(iJet).matchedGenJet)>0.4 || jets.at(iJet).p.Pt()>50) continue;
      float SF = this->getJetSF(jets.at(iJet).p.Pt(), jets.at(iJet).p.Eta());
      float eff = this->getJetEff(jets.at(iJet).p.Pt(), jets.at(iJet).p.Eta());
      
      if(jets.at(iJet).PileupIDloose){
         P_MC = P_MC * eff;
         P_Data = P_Data * SF * eff;
      }
      else{
         P_MC = P_MC * (1 - eff);
         P_Data = P_Data  * (1 - SF * eff);
      }
   }
   if (P_MC > 0. && P_Data > 0.) return P_Data/P_MC;
   else return 1.;
}

const float JetPileupIDWeights::getJetSF(const float& pt, const float& eta)
{
   switch(systematicInternal_){
      case SystematicInternal::nominal:
         return sf_hist_->GetBinContent(sf_hist_->GetXaxis()->FindBin(pt),sf_hist_->GetYaxis()->FindBin(eta));
         break;
      case SystematicInternal::vary_up:
         return sf_hist_->GetBinContent(sf_hist_->GetXaxis()->FindBin(pt),sf_hist_->GetYaxis()->FindBin(eta))+sfUnc_hist_->GetBinContent(sfUnc_hist_->GetXaxis()->FindBin(pt),sfUnc_hist_->GetYaxis()->FindBin(eta));
         break;
      case SystematicInternal::vary_down:
         return sf_hist_->GetBinContent(sf_hist_->GetXaxis()->FindBin(pt),sf_hist_->GetYaxis()->FindBin(eta))-sfUnc_hist_->GetBinContent(sfUnc_hist_->GetXaxis()->FindBin(pt),sfUnc_hist_->GetYaxis()->FindBin(eta));
         break;
      default:
         return 1.;
   }
}

const float JetPileupIDWeights::getJetEff(const float& pt, const float& eta)
{
   return eff_hist_->GetBinContent(eff_hist_->GetXaxis()->FindBin(pt),eff_hist_->GetYaxis()->FindBin(eta));
   
}
