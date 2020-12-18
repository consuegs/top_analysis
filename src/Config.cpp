#include "Config.hpp"

#include "tools/io.hpp"
#include "tools/util.hpp"

#include <TROOT.h>

//static
Config& Config::get(){
   static Config instance;
   return instance;
}

Config::Config()
{
   boost::property_tree::ptree pt;
   std::string cfgFile(CMAKE_SOURCE_DIR);
   cfgFile+="config.ini";
   boost::property_tree::read_ini(cfgFile,pt);

   treeVersion=pt.get<std::string>("input.version");
   year=pt.get<std::string>("input.year");
   treeName=pt.get<std::string>("input.treeName");
   dataBasePath=pt.get<std::string>("input.dataBasePath")+treeVersion+"/";
   gitHash=io::shellOutput("git log -1 --pretty=format:%h");
   if (TString(io::shellOutput("git status")).Contains("modified")) gitHash+="*";
   lumi=pt.get<float>("general.lumi");
   
   trigger_SF_ee=pt.get<std::string>("sf.trigger_SF_ee");
   trigger_SF_mumu=pt.get<std::string>("sf.trigger_SF_mumu");
   trigger_SF_emu=pt.get<std::string>("sf.trigger_SF_emu");
   
   applyDNN=pt.get<bool>("others.applyDNN");

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

