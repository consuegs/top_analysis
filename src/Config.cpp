#include "Config.hpp"

#include "tools/io.hpp"
#include "tools/util.hpp"
#include <stdio.h>
#include <stdlib.h>

#include <TROOT.h>

//static
Config& Config::get(){
   static Config instance;
   return instance;
}

void Config::setOutput(const std::string output){
   outputDirectory=output;
}

Config::Config()
{
   // ~boost::property_tree::ptree pt;
   std::string cfgFile(CMAKE_SOURCE_DIR);
   std::string config_year=getenv("ANALYSIS_YEAR_CONFIG");
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
   dataBasePath=pt.get<std::string>("input.dataBasePath")+treeVersion+"/";
   gitHash=io::shellOutput("git log -1 --pretty=format:%h");
   if (TString(io::shellOutput("git status")).Contains("modified")) gitHash+="*";
   lumi=pt.get<float>("general.lumi");
   
   trigger_SF_ee=pt.get<std::string>("sf.trigger_SF_ee");
   trigger_SF_mumu=pt.get<std::string>("sf.trigger_SF_mumu");
   trigger_SF_emu=pt.get<std::string>("sf.trigger_SF_emu");
   
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
   
   DeepCSV_loose=pt.get<float>("bTag_WP.DeepCSV_loose");
   DeepJet_loose=pt.get<float>("bTag_WP.DeepJet_loose");
   
   applyDNN=pt.get<bool>("DNN.applyDNN");
   DNN_Path=pt.get<std::string>("DNN.DNN_Path");

   lumiText=TString::Format("%.1f fb^{-1}",lumi*1e-3);
   sqrtsText=pt.get<std::string>("general.sqrtsText");
   extraText=pt.get<std::string>("general.extraText");

   outputDirectory=pt.get<std::string>("output.directory");
   datasets=DatasetCollection(pt,dataBasePath);
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

