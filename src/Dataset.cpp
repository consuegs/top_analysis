#include "Dataset.hpp"

#include "Config.hpp"
#include "tools/util.hpp"
#include "tools/io.hpp"

#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TTree.h>


/*
 * Arguments:
 * - `name`: simple name to use/display
 * - `files`: the files belonging to this dataset
 * - `xsec`: list of cross sections in pb ({0} for data files)
 */
Dataset::Dataset(std::string name,std::string label,std::string color,std::vector<std::string> files,std::vector<float> xsecs,float syst_unc,TString dataBasePath,bool isData,bool isSignal,bool isTTbar2L,bool isMadGraph,bool isPythiaOnly,std::string systName)
   : name(name)
   , label(label)
   , color(gROOT->ProcessLine((color+";").c_str()))
   , syst_unc(syst_unc)
   , isData(isData)
   , isSignal(isSignal)
   , isTTbar2L(isTTbar2L)
   , isMadGraph(isMadGraph)
   , isPythiaOnly(isPythiaOnly)
   , systName(systName)
{
   subsets.clear();
   for (uint i=0;i<files.size();i++){
      subsets.push_back(Datasubset(files[i],xsecs[i],dataBasePath,name,isData,isSignal,isTTbar2L,isMadGraph,isPythiaOnly));
   }
}

std::vector<TString> Dataset::getSubsetNames() const
{
   std::vector<TString> v;
   for (auto dss: subsets){
      v.push_back(dss.name);
   }
   return v;
}

/*
 * Arguments:
 * - `filename`: the filename of this datasubset
 * - `xsec`: cross section in pb (0 for data)
 */
Datasubset::Datasubset(std::string filename,float xsec,TString dataBasePath,std::string datasetName,bool isData,bool isSignal,bool isTTbar2L,bool isMadGraph,bool isPythiaOnly)
   : filename(filename)
   , xsec(xsec)
   , isData(isData)
   , isSignal(isSignal)
   , isTTbar2L(isTTbar2L)
   , isMadGraph(isMadGraph)
   , isPythiaOnly(isPythiaOnly)
   , datasetName(datasetName)
{
   std::vector<std::string> splitString;
   boost::split(splitString,filename,boost::is_any_of("."));
   assert(splitString.size()>=2);
   name=splitString[0];
   
   TFile* f = TFile::Open(dataBasePath+filename);
   if (!f->IsZombie()){
      TH1D*  h=(TH1D*)f->Get("TreeWriter/hCutFlow");
      TTree* t=(TTree*)f->Get("TreeWriter/eventTree");
      if (!h || !t) {
         debug_io<<TString::Format("%s is broken!",filename.c_str());
         Ngen=0;
         entries=0;
         Ngen_woWeight=0;
      } else {
         Ngen=h->GetBinContent(2);
         Ngen_woWeight=h->GetBinContent(1);
         entries=t->GetEntries();
      }
      f->Close();
   } else {
      debug_io<<TString::Format("%s is broken!",filename.c_str());
   }
}

/* return the full path to the file */
TString Datasubset::getPath() const{
   return Config::get().dataBasePath+filename;
}


// get sum of MC weight from Hist if Hist exists
double Datasubset::readNgenFromHist(const TString &histName, const int binNr) const
{
   TFile* f = TFile::Open(this->getPath());
   TH1D*  h=(TH1D*)f->Get("TreeWriter/"+histName);
   if (!h) {
      debug_io<<TString::Format("%s is broken! Can not read correct sum of mcWeights",histName.Data());
      exit(200);
   }
   else {
      return h->GetBinContent(binNr);
   }
}


// get correct sum of MC weights
double Datasubset::getNgen_syst(const Systematic::Systematic& systematic) const
{
   const Systematic::Type type = systematic.type();
   // if not mcWeight systematic return nominal Ngen
   if (std::find(Systematic::lumiRescalingTypes.begin(), Systematic::lumiRescalingTypes.end(), type) == Systematic::lumiRescalingTypes.end()){
      return Ngen;
   }
   else {
      // Check which variations is used
      bool upVariation = true;
      if(systematic.variation() == Systematic::down) upVariation = false;
      
      // BinID for facscaleUp, facScaleDown, renScaleUp, renScaleDown, scaleUp, scaleDown (differs between powheg and madgraph)
      std::vector<int> meScaleBins = {2,3,4,7,5,9};
      if(isMadGraph) meScaleBins = {4,7,2,3,5,9};
      
      switch(type){
         case Systematic::meFacScale:
            if(isPythiaOnly) return Ngen;    //no PDF weights available for pythiaOnly
            else return this->readNgenFromHist("hSystMCweight_PDF_",(upVariation)? meScaleBins[0] : meScaleBins[1]);
            break;
         case Systematic::meRenScale:
            if(isPythiaOnly) return Ngen;    //no PDF weights available for pythiaOnly
            else return this->readNgenFromHist("hSystMCweight_PDF_",(upVariation)? meScaleBins[2] : meScaleBins[3]);
            break;
         case Systematic::meScale:
            if(isPythiaOnly) return Ngen;    //no PDF weights available for pythiaOnly
            else return this->readNgenFromHist("hSystMCweight_PDF_",(upVariation)? meScaleBins[4] : meScaleBins[5]);
            break;
         case Systematic::alphasPdf:
            if(isPythiaOnly) return Ngen;    //no PDF weights available for pythiaOnly
            else return this->readNgenFromHist("hSystMCweight_PDF_",(upVariation)? 112 : 111);
            break;
         case Systematic::pdf:
            if (isPythiaOnly) return Ngen;  //no PDF weights available for pythiaOnly
            else return this->readNgenFromHist("hSystMCweight_PDF_",(upVariation)? 10+(2*systematic.variationNumber()) : 11+(2*systematic.variationNumber()));
            break;
         case Systematic::psISRScale:
            return this->readNgenFromHist("hSystMCweight_PS_",(upVariation)? 28 : 27);
            break;
         case Systematic::psFSRScale:
            return this->readNgenFromHist("hSystMCweight_PS_",(upVariation)? 6 : 5);
            break;
         case Systematic::bFrag:
            return this->readNgenFromHist("hSystMCweight_bFrag_",(upVariation)? 1 : 3);
            break;
         case Systematic::bSemilep:
            return this->readNgenFromHist("hSystMCweight_bFrag_",(upVariation)? 5 : 6);
            break;
         case Systematic::pu:
            return this->readNgenFromHist("hSystMCweight_PU_",(upVariation)? 2 : 3);
            break;
         default:
            std::cout<<"getNgen_syst: Correct normalization for "<<Systematic::convertType(type)<<" not found"<<std::endl;
            exit(200);
            break;
      }
   }
}


DatasetCollection::DatasetCollection(boost::property_tree::ptree const& pt,TString dataBasePath,bool single,std::string const &datasetMC_single,std::string const &datasetDATA_single,std::string const &datasetSIGNAL_single,int const fileNR)
{
   // MC
   std::vector<std::string> filenames;
   std::vector<float> xsecs;
   std::vector<float> kfacts;
   std::vector<float> feffs;
   std::string label;
   float syst_unc;
   bool isSplitSample=false;
   std::string systName;
   std::vector<std::string> mcDataset = util::to_vector<std::string>(pt.get<std::string>("input.mc_datasets"));
   if(single) mcDataset = util::to_vector<std::string>(datasetMC_single);  // Use only one dataset if single option is choosen
   for (std::string sDs: mcDataset){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      isSplitSample = pt.get_optional<bool>(sDs+".splitSample");
      if(isSplitSample && xsecs.size()==1) {   // Use split xsec for splitSamples
         xsecs[0] = xsecs[0]/(1.0*filenames.size());
         xsecs.resize(filenames.size(),xsecs[0]);
      }
      assert(filenames.size()==xsecs.size());
      kfacts = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      
      if (fileNR!=0 && single){ // use only selected file if single file option is choosen
         if(fileNR<=filenames.size()){
            filenames = {filenames[fileNR-1]};
            xsecs = {xsecs[fileNR-1]};
         }
         else{
            std::cerr<<"ERROR, fileNR too large for selected dataset"<<"\n...break\n"<<std::endl;
            exit(98);
         }
      }
      
      boost::optional<std::string> s_systName = pt.get_optional<std::string>(sDs+".systName");
      if (s_systName) systName = *s_systName;
      else systName="Nominal";
      
      boost::optional<bool> s_isTTbar2L = pt.get_optional<bool>(sDs+".isTTbar2L");     // Check if ttbar2l sample
      boost::optional<bool> s_isMadGraph = pt.get_optional<bool>(sDs+".isMadGraph");     // Check if madgraph sample
      boost::optional<bool> s_isPythiaOnly = pt.get_optional<bool>(sDs+".isPythiaOnly");     // Check if pythiaOnly sample
      
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      mc_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath,false,false,s_isTTbar2L,s_isMadGraph,s_isPythiaOnly,systName));
   }
   for (std::string sDs: util::to_vector<std::string>(pt.get<std::string>("input.mc_alternative_datasets"))){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      isSplitSample = pt.get_optional<bool>(sDs+".splitSample");
      if(isSplitSample && xsecs.size()==1) {   // Use split xsec for splitSamples
         xsecs[0] = xsecs[0]/(1.0*filenames.size());
         xsecs.resize(filenames.size(),xsecs[0]);
      }      assert(filenames.size()==xsecs.size());
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      
      boost::optional<std::string> s_systName = pt.get_optional<std::string>(sDs+".systName");
      if (s_systName) systName = *s_systName;
      else systName="Nominal";
      
      boost::optional<bool> s_isTTbar2L = pt.get_optional<bool>(sDs+".isTTbar2L");     // Check if ttbar2l sample
      boost::optional<bool> s_isMadGraph = pt.get_optional<bool>(sDs+".isMadGraph");     // Check if madgraph sample
      boost::optional<bool> s_isPythiaOnly = pt.get_optional<bool>(sDs+".isPythiaOnly");     // Check if pythiaOnly sample

      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      mc_alternative_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath,false,false,s_isTTbar2L,s_isMadGraph,s_isPythiaOnly,systName));
   }
   // Signals
   std::vector<std::string> signalDataset = util::to_vector<std::string>(pt.get<std::string>("input.signals"));
   if(single) signalDataset = util::to_vector<std::string>(datasetSIGNAL_single);  // Use only one dataset if single option is choosen
   for (std::string sDs: signalDataset){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs = util::to_vector<float>(pt.get<std::string>(sDs+".xsecs"));
      isSplitSample = pt.get_optional<bool>(sDs+".splitSample");
      if(isSplitSample && xsecs.size()==1) {   // Use split xsec for splitSamples
         xsecs[0] = xsecs[0]/(1.0*filenames.size());
         xsecs.resize(filenames.size(),xsecs[0]);
      }      assert(filenames.size()==xsecs.size());
      boost::optional<std::string> s_kfacts = pt.get_optional<std::string>(sDs+".kfact");
      if (s_kfacts) { // apply k factors
         kfacts=util::to_vector<float>(*s_kfacts);
         assert(kfacts.size()==xsecs.size());
         for (unsigned i=0; i<kfacts.size();i++) xsecs[i]*=kfacts[i];
      }
      boost::optional<std::string> s_feffs = pt.get_optional<std::string>(sDs+".filter_effs");
      if (s_feffs){ // apply filter efficiencies
         feffs = util::to_vector<float>(*s_feffs);
         assert(feffs.size()==xsecs.size());
         for (unsigned i=0; i<feffs.size();i++) xsecs[i]*=feffs[i];
      }
      
      if (fileNR!=0 && single){ // use only selected file if single file option is choosen
         if(fileNR<=filenames.size()){
            filenames = {filenames[fileNR-1]};
            xsecs = {xsecs[fileNR-1]};
         }
         else{
            std::cerr<<"ERROR, fileNR too large for selected dataset"<<"\n...break\n"<<std::endl;
            exit(98);
         }      
      }
      
      boost::optional<std::string> s_systName = pt.get_optional<std::string>(sDs+".systName");
      if (s_systName) systName = *s_systName;
      else systName="Nominal";
      
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      syst_unc=pt.get<float>(sDs+".syst_unc", 1e6); // default is something huge, so set it if you want to use it!
      signal_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,syst_unc,dataBasePath,false,true,false,false,false,systName));
   }
   // Data
   std::vector<std::string> dataDataset = util::to_vector<std::string>(pt.get<std::string>("input.data_streams"));
   if(single) dataDataset = util::to_vector<std::string>(datasetDATA_single);  // Use only one dataset if single option is choosen
   for (std::string sDs: dataDataset){
      filenames = util::to_vector<std::string>(pt.get<std::string>(sDs+".files"));
      xsecs=std::vector<float>(filenames.size(),-1);
      assert(filenames.size()==xsecs.size());
      
      if (fileNR!=0 && single){ // use only selected file if single file option is choosen
         if(fileNR<=filenames.size()){
            filenames = {filenames[fileNR-1]};
            xsecs = {xsecs[fileNR-1]};
         }
         else{
            std::cerr<<"ERROR, fileNR too large for selected dataset"<<"\n...break\n"<<std::endl;
            exit(98);
         }      
      }
      
      label=pt.get<std::string>(sDs+".label", sDs); // use Dataset name if no explicit label given
      data_datasets_.push_back(Dataset(sDs,label,pt.get<std::string>(sDs+".color"),filenames,xsecs,0,dataBasePath,true));
   }
   fillReferenceMaps();
}
DatasetCollection::DatasetCollection(std::vector<Dataset> mc_datasets,std::vector<Dataset> mc_alternative_datasets,std::vector<Dataset> data_datasets,std::vector<Dataset> signal_datasets)
   : mc_datasets_(mc_datasets)
   , mc_alternative_datasets_(mc_alternative_datasets)
   , data_datasets_(data_datasets)
   , signal_datasets_(signal_datasets)
{
   fillReferenceMaps();
}

void DatasetCollection::fillReferenceMaps()
{
   mDatasets_.clear();
   mDatasubsets_.clear();

   for (Dataset &ds: mc_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: mc_alternative_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: signal_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
   for (Dataset &ds: data_datasets_){
      mDatasets_.emplace(ds.name,std::ref(ds));
      for (Datasubset &dss:ds.subsets){
         mDatasubsets_.emplace(dss.name,std::ref(dss));
      }
   }
}

std::vector<Datasubset> DatasetCollection::getDatasubsets(bool mc, bool signal, bool data) const
{
   std::vector<Datasubset> vDss;
   if (mc)for (Dataset const &ds: mc_datasets_) // only defaults, no alternatives
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);
   if (signal) for (Dataset const &ds: signal_datasets_)
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);
   if (data)for (Dataset const &ds: data_datasets_)
      for (Datasubset const& dss: ds.subsets) vDss.push_back(dss);

   return vDss;
}

std::vector<Datasubset> DatasetCollection::getDatasubsets(std::vector<TString> dsNames) const
{
   std::vector<Datasubset> vDss;
   for (TString const &dsName: dsNames){
      for (Datasubset const &dss: getDataset(dsName).subsets){
         vDss.push_back(dss);
      }
   }
   return vDss;
}


Dataset DatasetCollection::getDataset(TString const &dsName) const
{
   return mDatasets_.at(dsName.Data());
}

Datasubset DatasetCollection::getDatasubset(TString const &dssName) const
{
   return mDatasubsets_.at(dssName.Data());
}

// ~std::vector<std::string> DatasetCollection::getDatasetNames() const
// ~{
   // ~std::vector<std::string> v;
   // ~for (auto ds: mDatasubsets_){
      // ~v.push_back(ds.first);
   // ~}
   // ~return v;
// ~}
std::vector<TString> DatasetCollection::getDatasetNames() const
{
   std::vector<TString> v;
   for (auto ds: mDatasets_){
      v.push_back(ds.first);
   }
   return v;
}
std::vector<std::string> DatasetCollection::getDatasubsetNames() const
{
   std::vector<std::string> v;
   for (auto const &dss: mDatasubsets_){
      v.push_back(dss.first);
   }
   return v;
}
std::vector<std::string> DatasetCollection::getDatasubsetNames(std::vector<TString> dsNames) const
{
   std::vector<std::string> v;
   for (auto const &dss: getDatasubsets(dsNames)){
      v.push_back(dss.name);
   }
   return v;
}

TString DatasetCollection::getLabel(TString const &dsName) const
{
   //~ Config const &cfg=Config::get();
   //~ if (dsName=="efake") return cfg.efake.label;
   return getDataset(dsName).label;
}

float DatasetCollection::getSystUncert(TString dsName) const
{
   //~ Config const &cfg=Config::get();
   //~ if (dsName=="efake") return cfg.efake.syst_unc;
   //~ if (dsName=="TTcomb") dsName="TTGJets";
   return getDataset(dsName).syst_unc;
}
