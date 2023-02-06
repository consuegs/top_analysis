#include "distributionsPlottingHelper.hpp"
#include "Config.hpp"
#include <iostream>

Config const &cfgDistr=Config::get();
using namespace std::chrono;

systHists::systHists(TString const &systematicName, TString filePath, TString const &histPath, std::vector<TString> const &datasets, std::vector<TString> const &datasets_ttBar, std::vector<TString> const &datasets_ttBar2L, std::vector<TString> const &datasets_st):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(hist::Histograms<TH1F>(datasets)),
      systematicName_(systematicName),
      datasets_(datasets),
      datasets_ttBar_(datasets_ttBar),
      datasets_ttBar2L_(datasets_ttBar2L),
      datasets_st_(datasets_st),
      filePath_(filePath),
      histPath_(histPath)
      {
         const std::vector<Systematic::Type> typeVec = Systematic::fileIndependentTypes;
         if (std::find(typeVec.begin(), typeVec.end(), systematic_.type()) == typeVec.end()){
            // ~histReader_ = new io::RootFileReader(filePath,histPath);
         }
         else {      //currently only used for lumi unc.
            hasRootFile_ = false;
            filePath_.ReplaceAll(systematic_.name(),"Nominal");
            // ~histReader_ = new io::RootFileReader(filePath,histPath);
            if(cfgDistr.systUncFactor.find(systematic_.type_str()) != cfgDistr.systUncFactor.end()){
               float unc = cfgDistr.systUncFactor.at(systematic_.type_str()).first;
               sf_ = (systematic_.variation() == Systematic::up)? 1.0+unc : 1.0-unc;
               datasets_SF = cfgDistr.systUncFactor.at(systematic_.type_str()).second;
            }
            else{
               std::cout<<"Error: Factor for "<<systematic_.type_str()<<" not found in config"<<std::endl;
            }
         }
         
         // Adapt datasets if systematic uses alternative samples
         const std::vector<Systematic::Type> typeVec_alt = Systematic::altSampleTypes;
         if (std::find(typeVec_alt.begin(), typeVec_alt.end(), systematic_.type()) != typeVec_alt.end()){
            altSampleType = true;
            if(systematic_.type() != Systematic::tw_ds) datasets_ = datasets_ttBar;
            else {
               datasets_ = datasets_st;
               onlyST = true;
            }
         }
         
         // Adapt datasets if syst. only changes ttbar (currently used for bFrag and bSemi)
         if (systematic_.type() == Systematic::bFrag || systematic_.type() == Systematic::bSemilep){
            onlyTTbar = true;
            datasets_ = datasets_ttBar;
         }
         
         // Adapt datasets for pdf unc. since only ttbar_dilepton is used
         if (systematic_.type() == Systematic::pdf || systematic_.type() == Systematic::pdf_envelope){
            datasets_ = datasets_ttBar2L_;
            datasets_ttBar_ = datasets_ttBar2L_;
            onlyTTbar2L = true;
         }
         
      }
      
void systHists::openFile(bool const &standardDict){
   histReader_ = new io::RootFileReader(filePath_,histPath_,standardDict,true);
}

void systHists::combineChannel(){
   hists_.combineChannel("/all",{"/ee","/mumu","/emu"});
}

// Return sample vectors to import
void distributionsplotting::getSampleVectors(int const &year_int, std::vector<TString> &mcSamples, std::vector<TString> &dataSamples, std::vector<TString> &ttbarSamples, std::vector<TString> &signalSamples, std::vector<TString> &stSamples, std::vector<TString> &bsmSamples){
   switch(year_int){
      case(3): //2018
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50_NLO","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      stSamples = {"SingleTop"};
      bsmSamples = {"T2tt_525_438","T2tt_525_350"};
      break;
      case(2): //2017
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50_NLO","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      stSamples = {"SingleTop"};
      break;
      case(1): //2016_postVFP 
      case(0): //2016_preVFP
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50_NLO","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      stSamples = {"SingleTop"};
      break;
   }
}

//Import hists for Nominal and systematics
void distributionsplotting::importHists(std::vector<systHists*> &systHists_vec, std::vector<TString> const &samplesToPlot, std::vector<TString> const &mcSamples,
                     std::vector<distr> const &vecDistr, std::vector<distr2D> const &vecDistr2D, bool const &standardDict) {
   for (systHists* &current : systHists_vec){
      std::vector<TString> inputSamples = current->datasets_;
      current->openFile(standardDict);
      
      std::cout<<"Import "<<current->systematicName_<<std::endl;
      
      for (TString sSample : inputSamples){
         current->hists_.setCurrentSample(sSample);
         TString sSample_alt = sSample;
         if (current->altSampleType) sSample_alt = sSample+"_"+current->systematicName_;    // get correct Sample name if alternative Sample Type
         //import 1D Hists
         for (auto const &distr_:vecDistr){
            TString loc;
            loc=distr_.path+distr_.name;
            TH1F* tempHist=current->histReader_->read<TH1F>(loc+"/"+sSample_alt);
            tempHist->SetDirectory(nullptr);
            
            TH1F rebinHist;   //rebin hist
            if (distr_.binEdges.empty()) {   // use same bin width
               rebinHist = hist::rebinned(*tempHist,distr_.xMin,distr_.xMax,distr_.nBins);
            }
            else{
               rebinHist = hist::rebinned(*tempHist,distr_.binEdges);
            }
                        
            if (!current->hasRootFile_){     //rescale only if dataset is affected
               if (std::find(current->datasets_SF.begin(),current->datasets_SF.end(),sSample) != current->datasets_SF.end()){
                  rebinHist.Scale(current->sf_);
               }
            }
            current->hists_.addFilledHist(loc,sSample,rebinHist);

         }
         //import 2D Hists
         for (auto const &distr_:vecDistr2D){
            TString loc;
            loc=distr_.path+distr_.name;
            TH2F* tempHist=current->histReader_->read<TH2F>(loc+"/"+sSample_alt);
            tempHist->SetDirectory(nullptr);
                        
            TH2F rebinHist;   //rebin hist
            if (distr_.binEdgesX.empty()) {   // use same bin width
               rebinHist = hist::rebinned(*tempHist,distr_.xMin,distr_.xMax,distr_.nBinsX,distr_.yMin,distr_.yMax,distr_.nBinsY);
            }
            else{
               rebinHist = hist::rebinned(*tempHist,distr_.binEdgesX,distr_.binEdgesY);
            }
            
            TH1F trafoHist = hist::histTrafo_2D(&rebinHist);
            
            if (!current->hasRootFile_){     //rescale only if dataset is affected
               if (std::find(current->datasets_SF.begin(),current->datasets_SF.end(),sSample) != current->datasets_SF.end()){
                  trafoHist.Scale(current->sf_);
               }
            }
            
            current->hists_.addFilledHist(loc,sSample,trafoHist);
            
         }
      }
     current->histReader_->closeFile();
   }
}

//Function to add the three different categories
void distributionsplotting::add_Categories(TString const path, io::RootFileReader const &reader_hist, TH1F &out_hist) {
   TH1F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH1F*) reader_hist.read<TH1F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH1F) *(reader_hist.read<TH1F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

//Function to combine samples to combined backgrounds
void distributionsplotting::combineAllSamples(int const &year_int, hist::Histograms<TH1F>* hs, std::vector<TString> &mcSamples_merged){
   hs->combineSamples("Diboson",{"WW","WZ","ZZ"});
   hs->combineSamples("DrellYan_comb",{"DrellYan_NLO","DrellYan_M10to50"});
   // ~hs->combineSamples("DrellYan_comb",{"DrellYan","DrellYan_M10to50"});
   // ~hs->combineSamples("DrellYan_comb",{"DrellYan_NLO","DrellYan_M10to50_NLO"});
   hs->combineSamples("ttZ",{"ttZ_2L","ttZ_QQ"});
   hs->combineSamples("ttW/Z",{"ttW","ttZ"});
   hs->combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
   switch(year_int){
      case(3): //2018
      mcSamples_merged = {"ttW/Z","WJetsToLNu","Diboson","DrellYan_comb","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","EGamma","MuonEG","SingleMuon"});
      hs->combineSamples("MC",mcSamples_merged);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_comb","ttZ","ttW"});
      break;
      case(2): //2017
      mcSamples_merged = {"ttW/Z","WJetsToLNu","Diboson","DrellYan_comb","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      hs->combineSamples("MC",mcSamples_merged);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_comb","ttZ","ttW"});
      break;
      case(0): //2016_preVFP
      case(1): //2016_postVFP
      mcSamples_merged = {"ttW/Z","WJetsToLNu","Diboson","DrellYan_comb","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      hs->combineSamples("MC",mcSamples_merged);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_comb","ttZ","ttW"});
      break;
   }
}

void distributionsplotting::addShifts(const TH1F &tempShift,TH1F* hist_shiftUP,TH1F* hist_shiftDOWN){
   for (int i=0; i<=tempShift.GetNbinsX(); i++){
      float content = tempShift.GetBinContent(i);
      if (content>0) hist_shiftUP->SetBinContent(i,hist_shiftUP->GetBinContent(i)+content*content);
      else hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
   }
}

// get up and down shift for set of systematics (for individual period)
std::pair<TH1F*,TH1F*> distributionsplotting::getTotalSyst(TH1F* const &nominal, std::vector<systHists*> const &systHists_vec, TString const loc, TString const sample, bool run2Combi, bool envelope){
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   TH1F* tempSys;
   TH1F tempShift;
   TH1F* nominal_ttbar;
   TH1F* nominal_ttbar2L;
   TH1F* nominal_st;
   
   std::vector<std::map<TString,TH1F>> vec_systShifts(1);   //needed for envelopes and top mass
   bool useCRenvelope = false;
   bool useMEenvelope = false;
   bool useMTop = false;
   
   for (auto &current : systHists_vec){
      if (std::find(Systematic::nominalTypes.begin(), Systematic::nominalTypes.end(), current->systematic_.type()) != Systematic::nominalTypes.end()){  //Store sum of ttBar for syst. with alt. samples
         if (sample == "") {
            nominal_ttbar = current->hists_.getSummedHist(loc,current->datasets_ttBar_);
            nominal_ttbar2L = current->hists_.getSummedHist(loc,current->datasets_ttBar2L_);   //for PDF unc.
            nominal_st = current->hists_.getSummedHist(loc,current->datasets_st_);   //for DS unc.
         }
         else {
            nominal_ttbar = current->hists_.getHistogram(loc,sample);
            nominal_ttbar2L = current->hists_.getHistogram(loc,sample);
            nominal_st = current->hists_.getSummedHist(loc,current->datasets_st_);   //for DS unc.
         }
         continue;
         
      }
            
      if (sample == ""){      //Option to get shift for individual sample
         tempSys = current->hists_.getSummedHist(loc,current->datasets_);
      }   
      else if (std::find(current->datasets_.begin(), current->datasets_.end(), sample) !=  current->datasets_.end()){
         tempSys=current->hists_.getHistogram(loc,sample);
      }
      else tempSys = nominal_ttbar;
            
      if ((current->altSampleType && current->onlyST==false) || current->onlyTTbar) tempShift= phys::getSystShift(*nominal_ttbar,*tempSys);  // Choose correct reference
      else if (current->onlyTTbar2L) tempShift= phys::getSystShift(*nominal_ttbar2L,*tempSys);  // for PDF unc.
      else if (current->onlyST) tempShift= phys::getSystShift(*nominal_st,*tempSys);  // for DS unc.
      else tempShift= phys::getSystShift(*nominal,*tempSys);
      
      // Part of deriving ME envelope, CR envelope and MTOP (only for single year plotting, otherwise done in getTotalSystCombined)
      if (!run2Combi){  // Store individual shifts to combine in the following
         if(std::find(Systematic::meTypes.begin(), Systematic::meTypes.end(), current->systematic_.type()) != Systematic::meTypes.end()){
            vec_systShifts[0][current->systematicName_] = tempShift;
            useMEenvelope = true;
            continue;
         }
         if(std::find(Systematic::crTypes.begin(), Systematic::crTypes.end(), current->systematic_.type()) != Systematic::crTypes.end()){
            vec_systShifts[0][current->systematicName_] = tempShift;
            useCRenvelope = true;
            continue;
         }
         if(std::find(Systematic::mTopTypes.begin(), Systematic::mTopTypes.end(), current->systematic_.type()) != Systematic::mTopTypes.end()){
            vec_systShifts[0][current->systematicName_] = tempShift;
            useMTop = true;
            continue;
         } 
      }
                
      
      if(envelope){
         for (int i=0; i<=tempShift.GetNbinsX(); i++){
            float content = tempShift.GetBinContent(i);
            if (content>0 && (content>hist_shiftUP->GetBinContent(i))) hist_shiftUP->SetBinContent(i,content);
            else if (content<0 && content<hist_shiftDOWN->GetBinContent(i)) hist_shiftDOWN->SetBinContent(i,abs(content));
         }
      }
      else{
         addShifts(tempShift,hist_shiftUP,hist_shiftDOWN);
      }
   }
   
   //Derive envelopes and mTOP uncertainty
   if(useMEenvelope){
      TH1F MEenvelope_down = tunfoldplotting::getMESCALEenvelopeCombined(vec_systShifts,false);
      addShifts(MEenvelope_down,hist_shiftUP,hist_shiftDOWN);
      TH1F MEenvelope_up = tunfoldplotting::getMESCALEenvelopeCombined(vec_systShifts,true);
      addShifts(MEenvelope_up,hist_shiftUP,hist_shiftDOWN);
   }
   if(useCRenvelope){
      TH1F CRenvelope_down = tunfoldplotting::getCRenvelopeCombined(vec_systShifts,false);
      addShifts(CRenvelope_down,hist_shiftUP,hist_shiftDOWN);
      TH1F CRenvelope_up = tunfoldplotting::getCRenvelopeCombined(vec_systShifts,true);
      addShifts(CRenvelope_up,hist_shiftUP,hist_shiftDOWN);
   }
   if(useMTop){
      TH1F mTOP_down = tunfoldplotting::getMTOPuncCombined(vec_systShifts,false);
      addShifts(mTOP_down,hist_shiftUP,hist_shiftDOWN);
      TH1F mTOP_up = tunfoldplotting::getMTOPuncCombined(vec_systShifts,true);
      addShifts(mTOP_up,hist_shiftUP,hist_shiftDOWN);
   }
   
   if(!envelope){
      hist::sqrtHist(*hist_shiftUP);
      hist::sqrtHist(*hist_shiftDOWN);
   }
   
   return std::make_pair(hist_shiftDOWN,hist_shiftUP);
}

// Method to print uncertainties
void distributionsplotting::printUnc(TString name, const float &down, const float &up, const float &nominal){
   name = name.ReplaceAll("ELECTRON_SCALESMEARING","EL_SCALESMEARING");
   TString out = TString::Format(" %s & $%.1f(%.1f)$ & $%.1f(%.1f)$\\\\\n",name.ReplaceAll("_","\\_").Data(),down,down/nominal*100,up,up/nominal*100);
   std::cout<<out;
}

// print total yields per contribution
void distributionsplotting::printTotalYields(hist::Histograms<TH1F>* hs, std::vector<std::vector<systHists*>> const &systHists_vec, const std::vector<TString> &mcSamples_merged){
   std::vector<TString> outputSamples(mcSamples_merged);
   std::reverse(outputSamples.begin(),outputSamples.end());
   outputSamples.push_back("MC");
   outputSamples.push_back("data");
   
   float mcYield_total;
   float mcYield_total_down;
   float mcYield_total_up;
   std::map<TString,float> totalMap;
   
   for (TString cat:{"ee","emu","mumu","all"}){    //Get the number of events per category
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      // ~TH1F* mc_total=hs->getHistogram("baseline/"+cat+"/nBjets","MC");
      
      std::pair<TH1F*,TH1F*> syst;
      if (systHists_vec.size() == 1){    // single period
         syst = getTotalSyst(mc_total,systHists_vec[0],"cutflow/"+cat);
      }
      else {   // combined
         syst = getTotalSystCombined(systHists_vec,"cutflow/"+cat);
      }
      
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      float mcYield = mc_total->GetBinContent(5);
      float mcYield_down = syst.first->GetBinContent(5);
      float mcYield_up = syst.second->GetBinContent(5);
      mcYield_total += mcYield;
      mcYield_total_down += mcYield_down;
      mcYield_total_up += mcYield_up;
      bool isMCtotal = false;
      for (TString sample:outputSamples){
         isMCtotal = (sample=="MC");
         TH1F* temp_hist=hs->getHistogram("cutflow/"+cat,sample);
         // ~TH1F* temp_hist=hs->getHistogram("baseline/"+cat+"/nBjets",sample);
         float sampleYield = temp_hist->GetBinContent(5);
         totalMap[sample.ReplaceAll("_","\\_")] += sampleYield;
         if(!isMCtotal){
            // ~std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<sampleYield<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
            std::cout<<std::fixed<<sample<<"&"<<sampleYield<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
         }
         else{
            TString yieldUnc = TString::Format("$%.1f^{+%.1f}_{-%.1f}$",sampleYield,mcYield_up,mcYield_down);
            std::cout<<"\\hline"<<std::endl;
            // ~std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<yieldUnc<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
            std::cout<<std::fixed<<sample<<"&"<<yieldUnc<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
         }
         // ~if (sample=="MC") std::cout<<mcYield_down<<"   "<<mcYield_up<<std::endl;
         // ~if (sample=="MC") std::cout<<mcYield_down/mcYield*100<<"   "<<mcYield_up/mcYield*100<<std::endl;
      }
   }
}

// print breakdown of syst uncertainties
void distributionsplotting::printUncBreakDown(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples) {
   std::cout<<std::endl<<"-----------------Uncertainties-------------------"<<std::endl;
   std::cout<<"!!!!!!!!!!!!!!!Uncertainty breakdown does not take ME and CR envelopes and mTop unc correctly into accout!!!!!!!!!!!!!!!!"<<std::endl;
   for (TString cat:{"ee","emu","mumu"}){
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      // ~TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"DrellYan_NLO");
      for (int i=1; i<(systHists_vec.size()); i++){
         if ((i+1)<systHists_vec.size()){    // check if still room for up and down shift
            if (systHists_vec[i]->systematic_.type_str()==systHists_vec[i+1]->systematic_.type_str()){     // check if up and down shift
               std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i],systHists_vec[i+1]};
               std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"",true);
               // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"DrellYan_NLO");
               printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(5),syst.second->GetBinContent(5),mc_total->GetBinContent(5));
               i++;
               continue;
            }
         }
         std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i]};   // print unc. without up and down shift
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"",true);
         // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"DrellYan_NLO");
         printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(5),syst.second->GetBinContent(5),mc_total->GetBinContent(5));
      }
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat);     // print total uncertainty
      // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat,"DrellYan_NLO");     // print total uncertainty
      std::cout<<"\\hline"<<std::endl;
      printUnc("TOTAL",syst.first->GetBinContent(5),syst.second->GetBinContent(5),mc_total->GetBinContent(5));
   }
}

// print uncertainty shift for individual samples
void distributionsplotting::printShiftBySample(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples){
   std::cout<<std::endl<<"-----------------Uncertainties per sample-------------------"<<std::endl;
   for (TString cat:{"ee","emu","mumu"}){
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      float mcYield = mc_total->GetBinContent(5);
      for (TString sample:mcSamples){
         TH1F* temp_hist=hs->getHistogram("cutflow/"+cat,sample);
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(temp_hist,systHists_vec,"cutflow/"+cat,sample);
         float sampleYield = temp_hist->GetBinContent(5);
         float sampleYield_down = syst.first->GetBinContent(5);
         float sampleYield_up = syst.second->GetBinContent(5);
         std::cout<<sample<<"   "<<sampleYield_down<<"   "<<sampleYield_up<<std::endl;
         std::cout<<sample<<"   "<<sampleYield_down/sampleYield*100<<"   "<<sampleYield_up/sampleYield*100<<std::endl;
      }
   }
}

void distributionsplotting::fixAxis2D(TAxis* axis, const int &nBinsY){     // only valid for x range in met between 0 and 400
   axis->SetNdivisions(4*nBinsY);   //for ticks per y bins
   axis->ChangeLabel(4*nBinsY+1,-1,-1,-1,-1,-1," ");
   for (int i=0; i<nBinsY; i++){
      axis->ChangeLabel(i*4+1,-1,-1,-1,-1,-1,"  ");
      axis->ChangeLabel(i*4+2,-1,-1,-1,-1,-1,"100");
      axis->ChangeLabel(i*4+3,-1,-1,-1,-1,-1," ");
      axis->ChangeLabel(i*4+4,-1,-1,-1,-1,-1,"300");
   }
}

void distributionsplotting::drawVertLines2D(std::vector<float> const &binEdgesY, float const &lowerEnd, float const &upperEnd,bool text){// only valid for x range in met between 0 and 400
   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.025);
   aline->SetLineWidth(1.8);
   for (int i=1; i<(binEdgesY.size()-1); i++){
      aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
      if (text){
         // ~atext->DrawLatex(75+((i-1)*400),1.7*lowerEnd,TString::Format("%.2f<|#Delta#phi|<%.2f",binEdgesY[i-1],binEdgesY[i]));
         atext->DrawLatex(75+((i-1)*400),1.7*lowerEnd,TString::Format("%.2f<|#Delta#phi|<%.2f",binEdgesY[i-1],binEdgesY[i]));
         if (i==(binEdgesY.size()-2)) atext->DrawLatex(75+(i*400),1.7*lowerEnd,TString::Format("%.2f<|#Delta#phi|<%.2f",binEdgesY[i],binEdgesY[i+1]));    // for last bin
      }
   }
}

void distributionsplotting::plotHistograms(TString const &sPresel, TString const &sVar, hist::Histograms<TH1F>* const &hs, std::vector<TString> const &mcSamples_merged, 
                     std::map<const TString,Color_t> const & colormap, std::vector<std::vector<systHists*>> const &systHists_vec, io::RootFileSaver const &saver, bool plotStatUncExpData, bool is2D, const std::vector<float> &binEdgesY){
   
   gfx::SplitCan sp_can;
   if(is2D) sp_can.can_.Size(1200,600);
   sp_can.pU_.SetLogy();
   sp_can.pU_.cd();
   TString loc;
   loc=sPresel+sVar;
   THStack st_mc=hs->getStack(loc,mcSamples_merged,colormap);
   gfx::LegendEntries le=hs->getLegendEntries();
   
   //systematics
   TH1F* hist_mc = hs->getHistogram(loc,{"MC"});
   std::pair<TH1F*,TH1F*> syst;
   if (systHists_vec.size() == 1){  // single period plot
      syst = getTotalSyst(hist_mc,systHists_vec[0],loc);
   }
   else{    // combined plot
      syst = getTotalSystCombined(systHists_vec,loc);
   }
   TGraphAsymmErrors systGraph = hist::getErrorGraph(syst.first,syst.second,hist_mc,true);
   
   TString cat;   // set channel label
   if (loc.Contains("ee")) cat="ee";
   else if (loc.Contains("emu")) cat="e#mu";
   else if (loc.Contains("mumu")) cat="#mu#mu";
   else if (loc.Contains("all")) cat="all";
   if (sPresel.Contains("Met200")) cat+="  p_{T}^{miss}>200 GeV";
   TLatex label=gfx::cornerLabel(cat,1);
   
   if (sVar.Contains("phi")){    // set plotting ranges
      st_mc.SetMinimum(1);
      st_mc.SetMaximum(1e6);
   }
   st_mc.SetMinimum(1);
   if (is2D) st_mc.SetMaximum(5e3*st_mc.GetMaximum());
   else st_mc.SetMaximum(1e3*st_mc.GetMaximum());
   st_mc.Draw();     // draw stack
   if(sPresel.Contains("cutflow")) st_mc.GetXaxis()->SetRangeUser(0.5,6.5);
   if(is2D) {     //needed because title is missing
      st_mc.SetTitle(";;Events/Bin;");
      st_mc.GetYaxis()->SetTitleOffset(0.75);
   }
   
   systGraph.SetFillStyle(3001);
   systGraph.SetLineWidth(0);
   systGraph.SetFillColor(kGray);
   systGraph.Draw("same 2");    // draw syst.
   
   // data plotting part
   auto hist_data = hs->getHistogram(loc,{"data"});
   hist_data->SetLineColor(kBlack);
   hist_data->SetMarkerSize((is2D)? 0.4 : 0.5);
   bool plotData = false;
   if(sPresel.Contains("GOF2D") || !(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep")|| sVar.Contains("C_em")|| sVar.Contains("Met"))) {
      hist_data->Draw("same");
      le.append(*hist_data,"data","lep");
      plotData = true;
   }
   else if (!is2D && !sPresel.Contains("GOF2D")) {    // for MET distributions only plot up to pT=140
      for (int i=1; i<=hist_data->GetNbinsX(); i++){
         if (hist_data->GetXaxis()->GetBinUpEdge(i)>140) hist_data->SetBinContent(i,0.);
      }
      hist_data->Draw("same");
      le.append(*hist_data,"data","lep");
      plotData = true;
   }
   else if (is2D){
      int width_index = hist_data->GetNbinsX()/(binEdgesY.size()-1);
      float width = hist_data->GetXaxis()->GetBinUpEdge(width_index);      //width of one dPhi bin in MET
      int dphi_bin = 1;
      for (int i=1; i<=hist_data->GetNbinsX(); i++){
         if (i>dphi_bin*width_index) dphi_bin++;
         if (hist_data->GetXaxis()->GetBinUpEdge(i)-(width*(dphi_bin-1))>140) hist_data->SetBinContent(i,0.);
      }
      hist_data->Draw("same");
      le.append(*hist_data,"data","lep");
      plotData = true;
   }
   
   // bsm plotting part
   auto hists_BSM = hs->getHistograms(loc,{"TTbar_diLepton"});     //placeholder
   auto BSM_legend = hs->getLegendEntries();    //placeholder
   TH1F* sqrtB_hist = (TH1F*)hist_mc->Clone();
   TH1F* S_sqrtB_hist;
   hist::sqrtHist(*sqrtB_hist);
   if (systHists_vec.size() == 1){  // single period plot
      if (cfgDistr.year_int==3){    //Plot BSM currently only available for 2018
         // ~hists_BSM=hs->getHistograms(loc,{"T2tt_525_438","T2tt_525_350"});
         hists_BSM=hs->getHistograms(loc,{"T2tt_525_350"});
         for (auto const &h: hists_BSM) {
            h->Draw("same hist");
            h->SetLineWidth(1.8);
            S_sqrtB_hist = (TH1F*)h->Clone();
            S_sqrtB_hist->Divide(sqrtB_hist);
         }
         le+=hs->getLegendEntries();
         BSM_legend=hs->getLegendEntries();
      }
   }
   
   if (is2D){
      fixAxis2D(st_mc.GetXaxis(),binEdgesY.size()-1);
      drawVertLines2D(binEdgesY,st_mc.GetHistogram()->GetMinimum(),st_mc.GetHistogram()->GetMaximum(),true);
   }
   
   // redraw axis, draw label and legend
   // ~auto hists_SM=hs->getHistograms(loc,{"TTbar_diLepton"});
   // ~TH1F axis=*(hists_SM[0]);
   // ~axis.SetStats(0);
   // ~axis.SetYTitle("Events/Bin");
   // ~axis.Draw("same axis");
   TLegend leg=le.buildLegend(.42,.72,1-(gPad->GetRightMargin()+0.02),-1,2);
   leg.SetFillStyle(1001);
   leg.Draw();
   label.Draw();
   
   // ratio part
   sp_can.pL_.cd();
   TH1F ratio = hist::getRatio(*hist_data,st_mc,"data/MC",hist::ONLY1);
   TH1F ratio_mc = hist::getRatio(st_mc,st_mc,"data/MC",hist::ONLY1);
   
   // syst unc.
   TH1F* hist_mc_sysDown = (TH1F*)hist_mc->Clone();
   TH1F* hist_mc_sysUp = (TH1F*)hist_mc->Clone();
   hist_mc_sysDown->Add(syst.second,-1.);
   hist_mc_sysUp->Add(syst.first);
            
   TH1F ratio_mc_systDown = hist::getRatio(*hist_mc_sysDown,*hist_mc,"data/MC",hist::ONLY1);
   TH1F ratio_mc_systUp = hist::getRatio(*hist_mc_sysUp,*hist_mc,"data/MC",hist::ONLY1);
   
   TGraphAsymmErrors systGraphRatio = hist::getErrorGraph(&ratio_mc_systUp,&ratio_mc_systDown,&ratio_mc,false);
   
   // total unc.
   TH1F* hist_mc_totalDown = (TH1F*)syst.second->Clone();
   TH1F* hist_mc_totalUp = (TH1F*)syst.first->Clone();
   
   TH1F* hist_mc_expDataStatDown = (TH1F*)hist_mc->Clone();
   TH1F* hist_mc_expDataStatUp = (TH1F*)hist_mc->Clone();
   
   for (int i=0; i<=hist_mc_totalUp->GetNbinsX(); i++){
      float stat2 = hist_mc->GetBinError(i)*hist_mc->GetBinError(i);
      float down2 = hist_mc_totalDown->GetBinContent(i)*hist_mc_totalDown->GetBinContent(i);
      float up2 = hist_mc_totalUp->GetBinContent(i)*hist_mc_totalUp->GetBinContent(i);
      hist_mc_totalDown->SetBinContent(i,sqrt(stat2+down2));
      hist_mc_totalUp->SetBinContent(i,sqrt(stat2+up2));
      
      if (plotStatUncExpData) {
         float binContent = hist_mc->GetBinContent(i);
         hist_mc_expDataStatDown->SetBinContent(i,binContent-sqrt(binContent));
         hist_mc_expDataStatUp->SetBinContent(i,binContent+sqrt(binContent));
      }
   }
   
   hist_mc_totalDown->Add(hist_mc,-1.);
   hist_mc_totalDown->Scale(-1.);
   hist_mc_totalUp->Add(hist_mc);
            
   TH1F ratio_mc_totalDown = hist::getRatio(*hist_mc_totalDown,*hist_mc,"data/MC",hist::ONLY1);
   TH1F ratio_mc_totalUp = hist::getRatio(*hist_mc_totalUp,*hist_mc,"data/MC",hist::ONLY1);
   
   TH1F ratio_mc_expDataStatDown = hist::getRatio(*hist_mc_expDataStatDown,*hist_mc,"data/MC",hist::ONLY1);
   TH1F ratio_mc_expDataStatUp = hist::getRatio(*hist_mc_expDataStatUp,*hist_mc,"data/MC",hist::ONLY1);
   
   TGraphAsymmErrors totalUncGraphRatio = hist::getErrorGraph(&ratio_mc_totalUp,&ratio_mc_totalDown,&ratio_mc,false);
   
   TGraphAsymmErrors expDataStatRatio = hist::getErrorGraph(&ratio_mc_expDataStatUp,&ratio_mc_expDataStatDown,&ratio_mc,false);
   
   if(sPresel.Contains("cutflow")){    // set cutflow specific axis labels
      ratio_mc.GetXaxis()->SetBinLabel(1,"diLepton");
      ratio_mc.GetXaxis()->SetBinLabel(2,"mll");
      ratio_mc.GetXaxis()->SetBinLabel(3,"jets");
      ratio_mc.GetXaxis()->SetBinLabel(4,"btag");
      ratio_mc.GetXaxis()->SetBinLabel(5,"DNN MET");
      // ~ratio_mc.GetXaxis()->SetBinLabel(7,"(addLepton veto)");
      ratio_mc.GetXaxis()->SetRangeUser(0.5,6.5);
      ratio_mc.GetXaxis()->SetLabelOffset(0.03);
   }
   
   if(is2D){
      ratio_mc.GetYaxis()->SetTitleOffset(0.25);
      ratio_mc.SetMaximum(1.35);
      ratio_mc.SetMinimum(0.65);
   }
   else{
      ratio_mc.GetYaxis()->SetTitleOffset(0.45);
      ratio_mc.SetMaximum(1.35);
      ratio_mc.SetMinimum(0.65);
   }
   ratio_mc.SetStats(0);
   ratio.SetLineColor(kBlack);
   totalUncGraphRatio.SetFillColor(kGray);
   totalUncGraphRatio.SetLineWidth(0);
   expDataStatRatio.SetFillColor(kBlue);
   systGraphRatio.SetFillColor(kGray+2);
   systGraphRatio.SetFillStyle(3004);
   expDataStatRatio.SetFillStyle(3005);
   expDataStatRatio.SetLineWidth(0);
   systGraphRatio.SetLineWidth(0);
   systGraphRatio.SetMarkerSize(0);
   ratio_mc.SetMarkerSize(0);
   ratio_mc.Draw("e2");    // only for axis
   totalUncGraphRatio.Draw("same e2");
   systGraphRatio.Draw("same e2");
   if(plotStatUncExpData) expDataStatRatio.Draw("same e2");
   ratio_mc.Draw("axis same");
   if(is2D){
      ratio.SetMarkerSize(0.2);
      ratio.SetLineWidth(1.5);
   }
   gPad->RedrawAxis("g");
   if(plotData) ratio.Draw("pe1 same");
   
   if(is2D){
      fixAxis2D(ratio_mc.GetXaxis(),binEdgesY.size()-1);
      drawVertLines2D(binEdgesY,ratio_mc.GetMinimum(),ratio_mc.GetMaximum(),false);
   }
   
   gfx::LegendEntries le_low;
   le_low.append(totalUncGraphRatio,"#sigma_{tot.}","f");
   le_low.append(systGraphRatio,"#sigma_{syst.}","f");
   if(plotStatUncExpData) le_low.append(expDataStatRatio,"#sigma_{stat.(exp. data)}","f");
   TLegend leg_low=le_low.buildLegend(.45,.8,0.75,0.95,2);
   leg_low.Draw();
      
   saver.save(sp_can,loc,false,true);
   
   //Draw S over sqrt(B) plots
   if (systHists_vec.size() == 1 && cfgDistr.year_int==3){  // single period plot and year 2018
      TCanvas can;
      can.cd();
      if(is2D) can.Size(1200,600);
      
      S_sqrtB_hist->GetYaxis()->SetTitle("S/#sqrt{B}");
      S_sqrtB_hist->SetStats(0);
      S_sqrtB_hist->Draw("hist");
      
      le.clear();
      le.append(*S_sqrtB_hist,"T2tt_525_350","l");
      TLegend leg_sensitivity=le.buildLegend(.5,.8,0.7,0.9,2);
      leg_sensitivity.SetFillStyle(1001);
      leg_sensitivity.Draw();
      label.Draw();
      
      if(is2D){
         fixAxis2D(S_sqrtB_hist->GetXaxis(),binEdgesY.size()-1);
         S_sqrtB_hist->SetMaximum(1.5);
         drawVertLines2D(binEdgesY,S_sqrtB_hist->GetMinimum(),S_sqrtB_hist->GetMaximum(),true);
      }
      
      saver.save(can,"S_sqrtB/"+loc,true,true,true);
   }
   
   
   //normalized distributions
   /*
   can.cd();
   auto ttbar_hist=hs->getHistogram(loc,"TTbar_diLepton");
   ttbar_hist->Scale(1.0/(ttbar_hist->Integral()));
   ttbar_hist->SetFillColor(ttbar_hist->GetLineColor());
   // ~ttbar_hist->SetFillStyle(1001);
   auto SMbkg_hist=hs->getHistogram(loc,"SM bkg.");
   SMbkg_hist->Scale(1.0/(SMbkg_hist->Integral()));
   SMbkg_hist->SetFillColor(kGray);
   // ~SMbkg_hist->SetFillStyle(1001);
   THStack st_norm;
   st_norm.Add(SMbkg_hist);
   st_norm.Add(ttbar_hist);
   axis.SetStats(0);
   axis.SetMaximum(1.0);
   axis.SetMinimum(1e-3);
   axis.GetYaxis()->SetTitle("normalized distribution");
   axis.Draw("axis");
   // ~st_norm.Draw("hist same");
   ttbar_hist->SetMarkerSize(0);
   SMbkg_hist->SetMarkerSize(0);
   ttbar_hist->Draw("hist e same");
   SMbkg_hist->Draw("hist e same");
   if (cfg.year_int==1){    //Plot BSM in 2016
      for (auto const &h: hists_BSM) {
         h->Scale(1.0/(h->Integral()));
         h->SetMarkerSize(0);
         h->Draw("same hist e");
      }
   }
   axis.Draw("axis same");
   le.clear();
   le.prepend(*ttbar_hist,"tt ll","l");
   le.append(*SMbkg_hist,"SM bkg.","l");
   if (cfg.year_int==1) le+=BSM_legend;
   TLegend leg2=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
   leg2.Draw();
   label.Draw();
   saver.save(can,"normalized/"+loc,true,false,true);
   */
}

void distributionsplotting::getCombinedDistributions(hist::Histograms<TH1F> &hs_combined, std::vector<hist::Histograms<TH1F>*> const &hs_vec, std::vector<TString> const &mcSamples_merged){
   for (auto const &loc : hs_vec[0]->getVariableNames()) {
      for (auto const &sample : mcSamples_merged) {
         TH1F* temp = hs_vec[0]->getHistogram(loc,sample);
         temp->Add(hs_vec[1]->getHistogram(loc,sample));
         temp->Add(hs_vec[2]->getHistogram(loc,sample));
         temp->Add(hs_vec[3]->getHistogram(loc,sample));
         hs_combined.addFilledHist(loc,sample,*temp);
      }
   }
}

// get up and down shift for set of systematics for all periods combined (partly copied from getCombinedUnc in tunfoldPlottingHelper)
std::pair<TH1F*,TH1F*> distributionsplotting::getTotalSystCombined(std::vector<std::vector<systHists*>> const &systHists_vec_all, TString const loc){
   
   // ~auto start = high_resolution_clock::now();
   
   // get nominal histograms per period
   std::vector<TH1F*> nominals(4);
   for (int i=0; i<nominals.size(); i++){
      nominals[i] = systHists_vec_all[i][0]->hists_.getHistogram(loc,{"MC"});
   }
   
   // define empty histograms to store total shifts
   TH1F* hist_TotalShiftUP = (TH1F*)nominals[0]->Clone();
   TH1F* hist_TotalShiftDOWN = (TH1F*)nominals[0]->Clone();
   hist_TotalShiftUP->Reset();
   hist_TotalShiftDOWN->Reset();
   TH1F* zeroes = (TH1F*)nominals[0]->Clone();
   zeroes->Reset();
   
   // define container for individual and combined shifts
   std::vector<std::map<TString,TH1F>> vec_systShifts(4);
   std::map<TString,TH1F> map_combinedShifts;
   TH1F tempHist;
   
   // covariance matrix for lumi uncertainty
   float lumiCorr[5][3] = {
      {1.0,0.0,0.0},
      {0.0,2.0,0.0},
      {0.0,0.0,1.5},
      {0.6,0.9,2.0},
      {0.0,0.6,0.2}
   };
   
   // booleans to check if lumi or envelope uncertainty was already produced
   bool lumiDone = false;
   bool CRenvelopeDone = false;
   bool MEenvelopeDone = false;
   bool MTopDone = false;
   
   // loop over systematics to store individual shifts
   for (int i=1; i<(systHists_vec_all[3].size()); i++){
      
      for (int j=0; j<nominals.size(); j++){ // loop over years
         
         if (i>=systHists_vec_all[j].size()) continue;   // handle different number of systematics per year (mainly for HEM)
         
         TString syst = systHists_vec_all[j][i]->systematicName_;
                  
         if (syst.BeginsWith("LUMI")) continue;    // Lumi derived differently

         if (syst == "JESUserDefinedHEM1516_DOWN" && j!=3) continue; // HEM only used for 2018
                  
         std::vector<systHists*> tempVec = {systHists_vec_all[j][0],systHists_vec_all[j][i]};   // temp vec with only one source
         std::pair<TH1F*,TH1F*> tempPair = getTotalSyst(nominals[j],tempVec,loc,"",true);

         tempPair.second->Add(tempPair.first,-1);   // get one histogram with up and down shifts (no cancelation since only one syst at a time is considered)
         vec_systShifts[j][syst] = *tempPair.second;
         
      }
   }
   
   for (int i=1; i<(systHists_vec_all[3].size()); i++){
      TString syst = systHists_vec_all[3][i]->systematicName_;
      
      if (syst.BeginsWith("LUMI")) { // Lumi unc correlation more complex then just (un)corr.
         if (!lumiDone){
            nominals[1]->Add(nominals[0]);   // add preVFP and post VFP
            
            // uncorrelated part
            TH1F tempHist_uncorr = *(TH1F*)nominals[1]->Clone();
            tempHist_uncorr.Reset();
            tempHist_uncorr.Add(nominals[1],lumiCorr[0][0]*1e-2);
            
            for (int i=1; i<3; i++){
               hist::addQuadr(tempHist_uncorr,*nominals[i+1],lumiCorr[i][i]*1e-2);
            }
            
            // correlated part
            TH1F tempHist_corr = *(TH1F*)nominals[1]->Clone();
            tempHist_corr.Reset();
            tempHist_corr.Add(nominals[1],lumiCorr[3][0]*1e-2);
            for (int i=1; i<3; i++){
               tempHist_corr.Add(nominals[i+1],lumiCorr[3][i]*1e-2);
            }
            
            // correlated part (17/18)
            TH1F tempHist_corr1718 = *(TH1F*)nominals[2]->Clone();
            tempHist_corr1718.Reset();
            tempHist_corr1718.Add(nominals[2],lumiCorr[4][1]*1e-2);
            tempHist_corr1718.Add(nominals[3],lumiCorr[4][2]*1e-2);
            
            // combine parts
            hist::addQuadr(tempHist_uncorr,tempHist_corr);
            hist::addQuadr(tempHist_uncorr,tempHist_corr1718);
            
            map_combinedShifts["LUMI_DOWN"] = tempHist_uncorr;
            tempHist_uncorr.Scale(-1.);
            map_combinedShifts["LUMI_UP"] = tempHist_uncorr;
            
            lumiDone = true;
         }
         
         continue;
      }
      
      // derive CR envelope
      if ((std::find(Systematic::crTypes.begin(), Systematic::crTypes.end(), Systematic::convertType(syst)) != Systematic::crTypes.end()) && CRenvelopeDone == false) {
         map_combinedShifts["CR_ENVELOPE_DOWN"] = tunfoldplotting::getCRenvelopeCombined(vec_systShifts,false);
         map_combinedShifts["CR_ENVELOPE_UP"] = tunfoldplotting::getCRenvelopeCombined(vec_systShifts,true);
         CRenvelopeDone = true;
         continue;
      }
      else if ((std::find(Systematic::crTypes.begin(), Systematic::crTypes.end(), Systematic::convertType(syst)) != Systematic::crTypes.end())) continue;  //ignore shifts already used in envelope
      
      // derive ME envelope
      if ((std::find(Systematic::meTypes.begin(), Systematic::meTypes.end(), Systematic::convertType(syst)) != Systematic::meTypes.end()) && MEenvelopeDone == false) {
         map_combinedShifts["MESCALE_ENVELOPE_DOWN"] = tunfoldplotting::getMESCALEenvelopeCombined(vec_systShifts,false);
         map_combinedShifts["MESCALE_ENVELOPE_UP"] = tunfoldplotting::getMESCALEenvelopeCombined(vec_systShifts,true);
         MEenvelopeDone = true;
         continue;
      }
      else if ((std::find(Systematic::meTypes.begin(), Systematic::meTypes.end(), Systematic::convertType(syst)) != Systematic::meTypes.end())) continue;  //ignore shifts already used in envelope
      
      // derive mTop uncertainty
      if ((std::find(Systematic::mTopTypes.begin(), Systematic::mTopTypes.end(), Systematic::convertType(syst)) != Systematic::mTopTypes.end()) && MTopDone == false) {
         map_combinedShifts["MTOP_DOWN"] = tunfoldplotting::getMTOPuncCombined(vec_systShifts,false);
         map_combinedShifts["MTOP_UP"] = tunfoldplotting::getMTOPuncCombined(vec_systShifts,true);
         MTopDone = true;
         continue;
      }
      else if ((std::find(Systematic::mTopTypes.begin(), Systematic::mTopTypes.end(), Systematic::convertType(syst)) != Systematic::mTopTypes.end())) continue;  //ignore shifts already used in envelope
      
      if (syst == "JESUserDefinedHEM1516_DOWN") { // HEM only used for 2018
         map_combinedShifts[syst] = vec_systShifts[3][syst];
         continue;
      }
      
      if (Systematic::isCorrelated(syst)){   // Add correlated parts
         tempHist = vec_systShifts[0][syst];
         for (int i=1; i<vec_systShifts.size(); i++){
            tempHist.Add(&vec_systShifts[i][syst]);
         }
         map_combinedShifts[syst] = tempHist;
      }
      else {   // Add non-correlated parts
         if (std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), Systematic::convertType(syst)) == Systematic::upDownTypes.end()){
            std::cout<<"Error: Combination of non-correlated uncertainties between years not supported of uncertainty not up_down type"<<std::endl;
            exit(301);
         }
         
         TString systName_down;
         TString systName_up;
         if (Systematic::convertVariation(syst) == Systematic::up){
            systName_up = syst;
            systName_down = syst;
            systName_down.ReplaceAll("_UP","_DOWN");
         }
         else if (Systematic::convertVariation(syst) == Systematic::down){
            systName_down = syst;
            systName_up = syst;
            systName_up.ReplaceAll("_DOWN","_UP");
         }
         
         TH1F down(vec_systShifts[0][syst]);
         TH1F up(vec_systShifts[0][syst]);
         down.Reset();
         up.Reset();
         
         for (int i=0; i<vec_systShifts.size(); i++){
            std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{vec_systShifts[i][systName_down],vec_systShifts[i][systName_up]});
            hist::addQuadr(down,*envelopes.first);
            hist::addQuadr(up,*envelopes.second);
         }
         down.Scale(-1.);
         map_combinedShifts[systName_down] = down;
         map_combinedShifts[systName_up] = up;
      }
      
   }
   
   //Store for total uncertainty
   for (const auto &[key, value]: map_combinedShifts){
      for (int i=0; i<=value.GetNbinsX(); i++){
         float content = value.GetBinContent(i);
         if (content>0) hist_TotalShiftUP->SetBinContent(i,hist_TotalShiftUP->GetBinContent(i)+content*content);
         else hist_TotalShiftDOWN->SetBinContent(i,hist_TotalShiftDOWN->GetBinContent(i)+content*content);
      }
   }
   hist::sqrtHist(*hist_TotalShiftUP);
   hist::sqrtHist(*hist_TotalShiftDOWN);
   
   // ~auto stop = high_resolution_clock::now();
   // ~auto duration = duration_cast<microseconds>(stop - start);
   // ~std::cout<<duration.count()<<std::endl;
   
   return std::make_pair(hist_TotalShiftDOWN,hist_TotalShiftUP);
}
