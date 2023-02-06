#include "jetVetoMaps.hpp"
#include "Config.hpp"
#include <iostream>

JetVetoMaps::JetVetoMaps(const std::string& map_fileName, const std::string& map_histName, const std::string& map_histName_MC16, const Systematic::Systematic& systematic, const int year_int, const bool is2016) :
   systematic_(systematic),
   vetoMap_hist_(0),
   vetoMap_hist_MC16_(0),
   is2016_(is2016)
{
   std::cout << "--- Beginning preparation of jet veto maps\n";
   // Check existence of input files
   if(!io::fileExists(map_fileName)){
   std::cerr<<"ERROR in constructor of JetVetoMaps! "
      <<"VetoMap file does not exist\n...break\n"<<std::endl;
   exit(220);
   }
   
   // ~HEM1516only = (systematic.type()==Systematic::applyJetVetoMaps_HEM1516);
   HEM1516only = (year_int==3);  //apply HEM1516 only if UL18 is used
      
   // Read required histogram
   this->prepareHists(map_fileName,map_histName,map_histName_MC16,vetoMap_hist_,vetoMap_hist_MC16_);
   
}

void JetVetoMaps::prepareHists(const std::string& map_fileName, const std::string& map_histName, const std::string& map_histName_MC16,TH2* &vetoMap_hist_, TH2* &vetoMap_hist_MC16_)
{
   
   io::RootFileReader vetoMap_file(map_fileName.c_str(),"",false);

   // Access histograms
   vetoMap_hist_ = (TH2*)(vetoMap_file.read<TH2F>(map_histName));
   if (is2016_) vetoMap_hist_MC16_ = (TH2*)(vetoMap_file.read<TH2F>(map_histName_MC16));
   if(!vetoMap_hist_ || (is2016_ && !vetoMap_hist_MC16_)){
      std::cerr<<"Error in JetVetoMaps::prepareHists()! TH2 not found in: "<<map_fileName
             <<"\n...break\n"<<std::endl;
      exit(221);
   }

   // Store histogram in memory
   vetoMap_hist_->SetDirectory(0);
   if (is2016_) vetoMap_hist_MC16_->SetDirectory(0);
}

const bool JetVetoMaps::checkVetoMap(const std::vector<tree::Jet>& jets)
{
   for(size_t iJet=0; iJet<jets.size(); ++iJet){
      if(HEM1516only){
         if (jets.at(iJet).p.Eta()>0 || jets.at(iJet).p.Phi()>0) continue;    // Dirty fix to take only the part of the vetomaps, where HEM1516 is shown
      }
      if (vetoMap_hist_->GetBinContent(vetoMap_hist_->GetXaxis()->FindBin(jets.at(iJet).p.Eta()),vetoMap_hist_->GetYaxis()->FindBin(jets.at(iJet).p.Phi()))>0) return false;
      if (is2016_){
         if (vetoMap_hist_MC16_->GetBinContent(vetoMap_hist_MC16_->GetXaxis()->FindBin(jets.at(iJet).p.Eta()),vetoMap_hist_MC16_->GetYaxis()->FindBin(jets.at(iJet).p.Phi()))>0) return false;
      }
   }
   return true;
}
