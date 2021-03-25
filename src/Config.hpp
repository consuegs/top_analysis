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
   static Config& get();
   float lumi;
   TString trigger_SF_ee;
   TString trigger_SF_mumu;
   TString trigger_SF_emu;
   
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
   
   bool applyDNN;

   TString treeVersion;
   TString year;
   int year_int;
   TString treeName;
   TString gitHash;
   TString dataBasePath;

   TString outputDirectory;

   float processFraction;
   bool  releaseMode;
   std::string datasetMC_single;
   std::string datasetDATA_single;
   std::string datasetSIGNAL_single;
   bool multi=false;

   std::vector<std::string> modules;

   DatasetCollection datasets;

   // canvas decoration
   TString lumiText ;
   TString sqrtsText;
   TString extraText;
private:
   Config();
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
