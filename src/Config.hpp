#ifndef CONFIG_HPP__
#define CONFIG_HPP__

#include "Dataset.hpp"

#include <string>
#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#include <TString.h>

/* Singleton-like config class */
class Config
{
public:
   boost::property_tree::ptree pt;
   static Config& get(std::string = "");
   float lumi;
   TString trigger_SF;
   
   TString jer_SF_mc;
   TString jer_RES_mc;
   TString jer_SF_data;
   TString jer_RES_data;
   
   TString jes_Folder;
   TString jes_UNC_mc;
   TString jes_UNC_mc_puppi;
   std::vector<std::string> jes_UNC_data;
   std::vector<std::string> jes_UNC_data_puppi;
   std::string jes_UNC_mc_regrouped;
   std::string jes_UNC_mc_puppi_regrouped;
   TString getJESPath(const int,const bool) const;
   
   TString muonTrigg1;
   TString muonTrigg2;
   TString muonTrigg3;
   TString muonTrigg4;
   TString singleMuonTrigg1;
   TString singleMuonTrigg2;
   TString eleTrigg1;
   TString eleTrigg2;
   TString eleMuTrigg1;
   TString eleMuTrigg2;
   TString eleMuTrigg3;
   TString eleMuTrigg4;
   TString singleEleTrigg;
   
   TString bTagger;
   int bTagWP;
   float bTagWPcut;
   
   TString bTagger_alt;
   int bTagWP_alt;
   float bTagWPcut_alt;
   
   TString bTagSF_file;
   TString bTagEffPath;
   
   std::string electronID_file;
   std::string electronID_hist;
   std::string electronRECO_file;
   std::string electronRECO_hist;
   std::string muonID_file;
   std::string muonID_hist;
   std::string muonISO_file;
   std::string muonISO_hist;
   float muonDYunc;
   float electronDYunc;
   
   bool applyDNN;
   TString DNN_Path;

   TString treeVersion;
   TString year;
   int year_int;
   TString treeName;
   TString gitHash;
   std::string dataBasePath;
   std::string minTreePath;

   TString outputDirectory;

   float processFraction;
   bool  releaseMode;
   std::string datasetMC_single;
   std::string datasetDATA_single;
   std::string datasetSIGNAL_single;
   int fileNR;
   bool multi=false;
   std::string systematic;
   
   std::map<TString,std::pair<float,std::vector<std::string>>> systUncFactor;
   
   std::vector<std::string> tunfold_InputSamples;
   TString tunfold_ResponseSample;
   TString tunfold_ResponseSampleAlt;
   std::vector<std::string> tunfold_bkgSamples_ttbar;
   std::vector<std::string> tunfold_bkgSamples_other;
   bool tunfold_withPTreweight;
   TString tunfold_scalePTreweight;
   bool tunfold_withDNN;
   bool tunfold_withPF;
   bool tunfold_withSameBins;
   bool tunfold_withBSM;
   bool tunfold_withScaleFactor;
   bool tunfold_plotComparison;
   bool tunfold_plotToyStudies;

   std::vector<std::string> modules;

   DatasetCollection datasets;

   // canvas decoration
   TString lumiText ;
   TString sqrtsText;
   TString extraText;
   
   //manipulate
   void setOutput(const std::string);
   void setLumi(const float newLumi);
   
   int getTotalFileNR(std::string const) const;

   Config(std::string config_year="");
};

/* Singleton-like color class for automatic colors*/
class Color
{
public:
   static int next();
   static void reset();
private:
   std::vector<int> cols_={kGray,kGray+2,kAzure-8,kOrange-3,kCyan-3,kGreen-5,kRed-6,kGreen-6,kMagenta-6,kBlack};
   int size_=cols_.size();
   int pos_=-1;
   Color(){}
   static Color& instance_();
   int next_();
};
#endif /* CONFIG_HPP__ */
