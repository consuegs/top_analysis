#include "Config.hpp"

#include "tools/io.hpp"
#include "tools/util.hpp"
#include <stdio.h>
#include <stdlib.h>

#include <TROOT.h>

//static
Config& Config::get(std::string config_year){
   static Config instance(config_year);
   return instance;
}

void Config::setOutput(const std::string output){
   outputDirectory=output;
}

void Config::setLumi(const float newLumi){
   lumi=newLumi;
}

Config::Config(std::string config_year)
{
   // ~boost::property_tree::ptree pt;
   std::string cfgFile(std::string(getenv("PWD"))+"/../");
   if (config_year == ""){
      config_year=getenv("ANALYSIS_YEAR_CONFIG");
   }
   if (config_year!=NULL){
      std::cout<<"Running on year "+config_year+" (to change set ANALYSIS_YEAR_CONFIG variable)"<<std::endl;
      cfgFile+="config"+config_year+".ini";
   }
   else{
      std::cout<<"Analysis year is not defined (to change set ANALYSIS_YEAR_CONFIG variable)"<<std::endl;
      throw;
   }
   boost::property_tree::read_ini(cfgFile,pt);

   treeVersion=pt.get<std::string>("input.version");
   year=pt.get<std::string>("input.year");
   year_int=pt.get<int>("input.year_int");
   treeName=pt.get<std::string>("input.treeName");
   dataBasePath=pt.get<std::string>("input.dataBasePath")+treeVersion+"/nTuple/";
   minTreePath=pt.get<std::string>("input.dataBasePath")+treeVersion+"/minTrees/";
   gitHash=io::shellOutput("git log -1 --pretty=format:%h");
   if (TString(io::shellOutput("git status")).Contains("modified")) gitHash+="*";
   lumi=pt.get<float>("general.lumi");
   
   trigger_SF=pt.get<std::string>("triggerSF.trigger_SF");
   
   jer_SF_mc=pt.get<std::string>("jetCorrections.jer_SF_mc");
   jer_RES_mc=pt.get<std::string>("jetCorrections.jer_RES_mc");
   jer_SF_data=pt.get<std::string>("jetCorrections.jer_SF_data");
   jer_RES_data=pt.get<std::string>("jetCorrections.jer_RES_data");
   
   jes_Folder=pt.get<std::string>("jetCorrections.jes_Folder");
   jes_UNC_mc=pt.get<std::string>("jetCorrections.jes_UNC_mc");
   jes_UNC_mc_puppi=pt.get<std::string>("jetCorrections.jes_UNC_mc_puppi");
   jes_UNC_data=util::to_vector<std::string>(pt.get<std::string>("jetCorrections.jes_UNC_data"));
   jes_UNC_data_puppi=util::to_vector<std::string>(pt.get<std::string>("jetCorrections.jes_UNC_data_puppi"));
   jes_UNC_mc_regrouped=jes_Folder+pt.get<std::string>("jetCorrections.jes_UNC_mc_regrouped");
   jes_UNC_mc_puppi_regrouped=jes_Folder+pt.get<std::string>("jetCorrections.jes_UNC_mc_puppi_regrouped");
   
   muonTrigg1=pt.get<std::string>("trigger.muonTrigg1");
   muonTrigg2=pt.get<std::string>("trigger.muonTrigg2");
   muonTrigg3=pt.get<std::string>("trigger.muonTrigg3");
   muonTrigg4=pt.get<std::string>("trigger.muonTrigg4");
   singleMuonTrigg1=pt.get<std::string>("trigger.singleMuonTrigg1");
   singleMuonTrigg2=pt.get<std::string>("trigger.singleMuonTrigg2");
   eleTrigg1=pt.get<std::string>("trigger.eleTrigg1");
   eleTrigg2=pt.get<std::string>("trigger.eleTrigg2");
   eleMuTrigg1=pt.get<std::string>("trigger.eleMuTrigg1");
   eleMuTrigg2=pt.get<std::string>("trigger.eleMuTrigg2");
   eleMuTrigg3=pt.get<std::string>("trigger.eleMuTrigg3");
   eleMuTrigg4=pt.get<std::string>("trigger.eleMuTrigg4");
   singleEleTrigg=pt.get<std::string>("trigger.singleEleTrigg");
   
   bTagger=pt.get<std::string>("bTag.tagger");
   bTagWP=pt.get<int>("bTag.WP");
   bTagWPcut=pt.get<float>("bTag.WPcut");
   
   bTagger_alt=pt.get<std::string>("bTag_alternativ.tagger");
   bTagWP_alt=pt.get<int>("bTag_alternativ.WP");
   bTagWPcut_alt=pt.get<float>("bTag_alternativ.WPcut");
   
   electronID_file=pt.get<std::string>("leptonSF.electronFileNameID");
   electronID_hist=pt.get<std::string>("leptonSF.electronHistNameID");
   electronRECO_file=pt.get<std::string>("leptonSF.electronFileNameRECO");
   electronRECO_hist=pt.get<std::string>("leptonSF.electronHistNameRECO");
   muonID_file=pt.get<std::string>("leptonSF.muonFileNameID");
   muonID_hist=pt.get<std::string>("leptonSF.muonHistNameID");
   muonISO_file=pt.get<std::string>("leptonSF.muonFileNameISO");
   muonISO_hist=pt.get<std::string>("leptonSF.muonHistNameISO");
   muonDYunc=pt.get<float>("leptonSF.muonDYextrapolationunc");
   electronDYunc=pt.get<float>("leptonSF.electronDYextrapolationunc");
   
   jetPileupID_file=pt.get<std::string>("jetPileupID.jetPileupID_file");
   jetPileupID_effHist=pt.get<std::string>("jetPileupID.jetPileupID_effHist");
   jetPileupID_sfHist=pt.get<std::string>("jetPileupID.jetPileupID_sfHist");
   
   jetVetoMap_file=pt.get<std::string>("jetVetoMap.jetVetoMap_file");
   jetVetoMap_vetoMapHist=pt.get<std::string>("jetVetoMap.jetVetoMap_vetoMapHist");
   jetVetoMap_vetoMapHist_loose=pt.get<std::string>("jetVetoMap.jetVetoMap_vetoMapHist_loose");
   jetVetoMap_vetoMapHist_MC16=pt.get_optional<std::string>("jetVetoMap.jetVetoMap_vetoMapHist_MC16");
   
   applyDNN=pt.get<bool>("DNN.applyDNN");
   DNN_Path=pt.get<std::string>("DNN.DNN_Path");

   lumiText=TString::Format("%.1f fb^{-1}",lumi*1e-3);
   sqrtsText=pt.get<std::string>("general.sqrtsText");
   extraText=pt.get<std::string>("general.extraText");
   
   systUncFactor["LUMI"] = std::pair<float,std::vector<std::string>> (pt.get<float>("SystUnc.lumiUnc"),util::to_vector<std::string>(pt.get<std::string>("SystUnc.lumi_samples")));
   systUncFactor["XSEC_TTOTHER"] = std::pair<float,std::vector<std::string>> (pt.get<float>("SystUnc.xsec_ttother_unc"),util::to_vector<std::string>(pt.get<std::string>("SystUnc.xsec_ttother_samples")));
   systUncFactor["XSEC_DY"] = std::pair<float,std::vector<std::string>> (pt.get<float>("SystUnc.xsec_dy_unc"),util::to_vector<std::string>(pt.get<std::string>("SystUnc.xsec_dy_samples")));
   systUncFactor["XSEC_ST"] = std::pair<float,std::vector<std::string>> (pt.get<float>("SystUnc.xsec_st_unc"),util::to_vector<std::string>(pt.get<std::string>("SystUnc.xsec_st_samples")));
   systUncFactor["XSEC_OTHER"] = std::pair<float,std::vector<std::string>> (pt.get<float>("SystUnc.xsec_other_unc"),util::to_vector<std::string>(pt.get<std::string>("SystUnc.xsec_other_samples")));
   
   tunfold_InputSamples = util::to_vector<std::string>(pt.get<std::string>("TUnfold.inputSamples"));
   tunfold_ResponseSample = pt.get<std::string>("TUnfold.responseSample");
   tunfold_ResponseSampleAlt = pt.get<std::string>("TUnfold.responseSample_alt");
   tunfold_bkgSamples_ttbar = util::to_vector<std::string>(pt.get<std::string>("TUnfold.bkgSamples_ttbar"));
   tunfold_bkgSamples_other = util::to_vector<std::string>(pt.get<std::string>("TUnfold.bkgSamples_other"));
   tunfold_withPTreweight = pt.get<bool>("TUnfold.withPTreweight");
   tunfold_scalePTreweight = pt.get<std::string>("TUnfold.scalePTreweight");
   tunfold_withDNN = pt.get<bool>("TUnfold.withDNN");
   tunfold_withPF = pt.get<bool>("TUnfold.withPF");
   tunfold_withSameBins = pt.get<bool>("TUnfold.withSameBins");
   tunfold_withBSM = pt.get<bool>("TUnfold.withBSM");
   tunfold_withScaleFactor = pt.get<bool>("TUnfold.withScaleFactor");
   tunfold_plotComparison = pt.get<bool>("TUnfold.plotComparison");
   tunfold_plotToyStudies = pt.get<bool>("TUnfold.plotToyStudies");

   // Check if running on lx* with home directory excess, if not, output cannot be stored on /net/...
   TString hostname=getenv("HOSTNAME");
   if (hostname.Contains("lx")){
      outputDirectory=pt.get<std::string>("output.directory")+treeVersion+"/output_framework";
   }
   else{
      outputDirectory="../output_framework";
   }
   
   datasets=DatasetCollection(pt,dataBasePath);
   
   bTagSF_file=pt.get<std::string>("BTag_Weights.BTagSF_file");
   bTagEffPath=pt.get<std::string>("BTag_Weights.BTagEffPath");
}

TString Config::getJESPath(const int run_era, const bool isPuppi=false) const{
   if (run_era==0) {
      if(!isPuppi) return jes_Folder+jes_UNC_mc;
      else return jes_Folder+jes_UNC_mc_puppi;
   }
   else if (run_era<=jes_UNC_data.size()){
      if(!isPuppi) return jes_Folder+jes_UNC_data[run_era-1];
      else return jes_Folder+jes_UNC_data_puppi[run_era-1];
   }
   else {
      std::cerr<<"ERROR in selection for JES file, run_era does not match! "<<std::endl;
      exit(98);
   }
}

int Config::getTotalFileNR(std::string const sample) const{
   std::vector<std::string> tempVec = util::to_vector<std::string>(pt.get<std::string>(sample+".files"));
   return tempVec.size();
}

// static
Color& Color::instance_(){
   static Color instance;
   return instance;
}
// static
int Color::next(){
   return instance_().next_();
}
// static
void Color::reset(){
   instance_().pos_=-1;
}
int Color::next_(){
   pos_=(pos_+1)%size_;
   return cols_[pos_];
}

