#include "bTagWeights.hpp"
#include "Config.hpp"
#include <iostream>

//Currently fixed to DeepJet loose WP!!!!!!!!!

//////////////////////////////////// btag Eff //////////////////////////////////////////////////////

BTagEffMapFunctor::BTagEffMapFunctor(const std::string& bTagEffPath,const std::string& tagger, const BTagEntry::OperatingPoint& WP, const Systematic::Systematic& systematic) :
   systematic_(systematic),
   bTagEffPath_(bTagEffPath),
   tagger_(tagger),
   WP_(WP)
{
   std::cout << "--- Beginning preparation of bTag Eff\n";
   
   //Check if alternative eff should be used
   const Systematic::Type type = systematic.type();
   if(std::find(Systematic::noBtagEffTypes.begin(), Systematic::noBtagEffTypes.end(), type) == Systematic::noBtagEffTypes.end()){
      nominalEff_ = false;
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<" for bTag Eff\n";
   }
   else{
      nominalEff_ = true;
      std::cout<<"Use nominal bTag Eff\n";
   }
   this->readBTagEffHistos();
}

void BTagEffMapFunctor::readBTagEffHistos()
{
   // Check if files exist and open file
   TString fileName = TString::Format("%s/%s.root",bTagEffPath_.c_str(),(nominalEff_)? "Nominal" : systematic_.name().Data());
   if (!io::fileExists(fileName.Data())){
      std::cerr<<"ERROR in readBTagEffHistos! "
      <<fileName<<" does not exist\n...break\n"<<std::endl;
      exit(98);
   }
   TFile file_(fileName);
   
   // Construct correct histogram name
   TString taggerWP = tagger_;
   switch(WP_){
      case BTagEntry::OP_LOOSE: taggerWP+="_loose";break;
      case BTagEntry::OP_MEDIUM: taggerWP+="_medium";break;
      case BTagEntry::OP_TIGHT: taggerWP+="_tight";break;
      default:
         std::cerr<<"ERROR in readBTagEffHistos! "<<" tagger not found in efficiency file\n...break\n"<<std::endl;
         exit(98);
         break;
   }
   
   // Read histograms from file
   bEfficiencies_ee_ = *(TH2F*)file_.Get("ee/B_"+taggerWP);
	cEfficiencies_ee_ = *(TH2F*)file_.Get("ee/C_"+taggerWP);
	lightEfficiencies_ee_ = *(TH2F*)file_.Get("ee/Light_"+taggerWP);
	
	bEfficiencies_mumu_ = *(TH2F*)file_.Get("mumu/B_"+taggerWP);
	cEfficiencies_mumu_ = *(TH2F*)file_.Get("mumu/C_"+taggerWP);
	lightEfficiencies_mumu_ = *(TH2F*)file_.Get("mumu/Light_"+taggerWP);
	
	bEfficiencies_emu_ = *(TH2F*)file_.Get("emu/B_"+taggerWP);
	cEfficiencies_emu_ = *(TH2F*)file_.Get("emu/C_"+taggerWP);
	lightEfficiencies_emu_ = *(TH2F*)file_.Get("emu/Light_"+taggerWP);
}

const float BTagEffMapFunctor::getEff(const int& flavor, const double& pt, const double& eta, const int& channel)
{
   float eta_abs = abs(eta);      // Efficiencies stored as function ab |eta|
   
   // Read efficiency from histogram
   if (flavor == 5){
      if(channel == Channel::ee)return bEfficiencies_ee_.GetBinContent(bEfficiencies_ee_.GetXaxis()->FindBin(pt),bEfficiencies_ee_.GetYaxis()->FindBin(eta_abs));
      else if(channel == Channel::mumu)return bEfficiencies_mumu_.GetBinContent(bEfficiencies_mumu_.GetXaxis()->FindBin(pt),bEfficiencies_mumu_.GetYaxis()->FindBin(eta_abs));
      else return bEfficiencies_emu_.GetBinContent(bEfficiencies_emu_.GetXaxis()->FindBin(pt),bEfficiencies_emu_.GetYaxis()->FindBin(eta_abs));
   }
   else if (flavor == 4){
      if(channel == Channel::ee)return cEfficiencies_ee_.GetBinContent(cEfficiencies_ee_.GetXaxis()->FindBin(pt),cEfficiencies_ee_.GetYaxis()->FindBin(eta_abs));
      else if(channel == Channel::mumu)return cEfficiencies_mumu_.GetBinContent(cEfficiencies_mumu_.GetXaxis()->FindBin(pt),cEfficiencies_mumu_.GetYaxis()->FindBin(eta_abs));
      else return cEfficiencies_emu_.GetBinContent(cEfficiencies_emu_.GetXaxis()->FindBin(pt),cEfficiencies_emu_.GetYaxis()->FindBin(eta_abs));
   }
   else {
      if(channel == Channel::ee)return lightEfficiencies_ee_.GetBinContent(lightEfficiencies_ee_.GetXaxis()->FindBin(pt),lightEfficiencies_ee_.GetYaxis()->FindBin(eta_abs));
      else if(channel == Channel::mumu)return lightEfficiencies_mumu_.GetBinContent(lightEfficiencies_mumu_.GetXaxis()->FindBin(pt),lightEfficiencies_mumu_.GetYaxis()->FindBin(eta_abs));
      else return lightEfficiencies_emu_.GetBinContent(lightEfficiencies_emu_.GetXaxis()->FindBin(pt),lightEfficiencies_emu_.GetYaxis()->FindBin(eta_abs));
   }
}


//////////////////////////////////// btag Weights //////////////////////////////////////////////////////
BTagWeights::BTagWeights(const std::string& bTagSFFile, const std::string& bTagEffPath, const std::string& tagger, const BTagEntry::OperatingPoint& WP, const float& WPcut, const Systematic::Systematic& systematic) :
   systematic_(systematic),
   bTagEff_(BTagEffMapFunctor(bTagEffPath, tagger, WP, systematic_)),
   WPcut_(WPcut)
{
   std::cout << "--- Beginning preparation of bTag weights\n";
   // Check existence of input files
   if(!io::fileExists(bTagSFFile)){
   std::cerr<<"ERROR in constructor of BTagWeights! "
      <<"ScaleFactor file does not exist\n...break\n"<<std::endl;
   exit(98);
   }
   
   // Check if curreent tagger is DeepJet
   isDeepJet_ = tagger=="DeepJet";
   
   // Check if syst. should be applied for bTag SF
   std::string bcSyst = "central";
   std::string lightSyst = "central";
   
   const Systematic::Type type = systematic.type();
   if(std::find(Systematic::btagTypes.begin(), Systematic::btagTypes.end(), type) != Systematic::btagTypes.end()){
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      bcSyst = BTagSourcesMap.at(systematic_.name())[0];
      lightSyst = BTagSourcesMap.at(systematic_.name())[1];
   }
   else std::cout<<"Do not apply systematic variation\n";
      
   
   // Configure SF readers
   bTagCalibration_ = BTagCalibration(tagger,bTagSFFile);
   
   calibrationReader_bJets_ = BTagCalibrationReader(WP,bcSyst);
   calibrationReader_cJets_ = BTagCalibrationReader(WP,bcSyst);
   calibrationReader_lightJets_ = BTagCalibrationReader(WP,lightSyst);
   
   calibrationReader_bJets_.load(bTagCalibration_,BTagEntry::FLAV_B,"mujets");
   calibrationReader_cJets_.load(bTagCalibration_,BTagEntry::FLAV_C,"mujets");
   calibrationReader_lightJets_.load(bTagCalibration_,BTagEntry::FLAV_UDSG,"incl");
   
}

const float BTagWeights::getJetSF(const int& flavor, const float& eta, const float& pt)
{
   
   if (flavor == 5) return calibrationReader_bJets_.eval(BTagEntry::FLAV_B, eta, pt); 
   else if (flavor == 4) return calibrationReader_cJets_.eval(BTagEntry::FLAV_C, eta, pt); 
   else return calibrationReader_lightJets_.eval(BTagEntry::FLAV_UDSG, eta, pt);
   
}
   

const float BTagWeights::getEventWeight(const std::vector<tree::Jet>& jets, const int& channel)
{
   float P_MC = 1.;
   float P_Data = 1.;
   
   //Derive Scale Factor following 1a (https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods)
   for(size_t iJet=0; iJet<jets.size(); ++iJet){
      float SF = getJetSF(jets.at(iJet).hadronFlavour, jets.at(iJet).p.Eta(), jets.at(iJet).p.Pt());
      float eff = bTagEff_.getEff(jets.at(iJet).hadronFlavour, jets.at(iJet).p.Pt(), jets.at(iJet).p.Eta(), channel);
      
      float bTagValue = (isDeepJet_)? jets.at(iJet).bTagDeepJet : jets.at(iJet).bTagDeepCSV;
      
      if(bTagValue > WPcut_){
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

const std::map<TString, std::vector<std::string> > BTagWeights::BTagSourcesMap = {

  {"BTAGBC_DOWN",{"down","central"}},
  {"BTAGBC_UP",{"up","central"}},
  {"BTAGL_DOWN",{"central","down"}},
  {"BTAGL_UP",{"central","up"}},
  {"BTAGBC_CORR_DOWN",{"down_correlated","central"}},
  {"BTAGBC_CORR_UP",{"up_correlated","central"}},
  {"BTAGL_CORR_DOWN",{"central","down_correlated"}},
  {"BTAGL_CORR_UP",{"central","up_correlated"}},
  {"BTAGBC_UNCORR_DOWN",{"down_uncorrelated","central"}},
  {"BTAGBC_UNCORR_UP",{"up_uncorrelated","central"}},
  {"BTAGL_UNCORR_DOWN",{"central","down_uncorrelated"}},
  {"BTAGL_UNCORR_UP",{"central","up_uncorrelated"}}
};
   
   
