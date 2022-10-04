#include "distributionsPlottingHelper.hpp"
#include "Config.hpp"
#include <iostream>

Config const &cfgDistr=Config::get();

systHists::systHists(TString const &systematicName, TString filePath, TString const &histPath, std::vector<TString> const &datasets, std::vector<TString> const &datasets_ttBar, std::vector<TString> const &datasets_ttBar2L):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(hist::Histograms<TH1F>(datasets)),
      systematicName_(systematicName),
      datasets_(datasets),
      datasets_ttBar_(datasets_ttBar),
      datasets_ttBar2L_(datasets_ttBar2L),
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
            datasets_ = datasets_ttBar;
         }
         
         // Adapt datasets if syst. only changes ttbar (currently used for bFrag and bSemi) and scale unc. since not working for MadGraph Samples
         if (systematic_.type() == Systematic::bFrag || systematic_.type() == Systematic::bSemilep || systematic_.type() == Systematic::meFacScale || systematic_.type() == Systematic::meRenScale || systematic_.type() == Systematic::meScale || systematic_.type() == Systematic::meScale_envelope_ind){
         // ~if (systematic_.type() == Systematic::bFrag || systematic_.type() == Systematic::bSemilep){
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
void distributionsplotting::getSampleVectors(int const &year_int, std::vector<TString> &mcSamples, std::vector<TString> &dataSamples, std::vector<TString> &ttbarSamples, std::vector<TString> &signalSamples){
   switch(year_int){
      case(3): //2018
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      break;
      case(2): //2017
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      break;
      case(1): //2016_postVFP 
      case(0): //2016_preVFP
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      break;
   }
}

//Import hists for Nominal and systematics
void distributionsplotting::importHists(std::vector<systHists*> &systHists_vec, std::vector<TString> const &samplesToPlot, std::vector<TString> const &mcSamples,
                     std::vector<distr> const &vecDistr, std::vector<distr2D> const &vecDistr2D, bool const &standardDict) {
   for (systHists* &current : systHists_vec){
      std::vector<TString> inputSamples = current->datasets_;
      current->openFile(standardDict);
      
      for (TString sSample : inputSamples){
         // ~auto start = high_resolution_clock::now();
         current->hists_.setCurrentSample(sSample);
         TString sSample_alt = sSample;
         if (current->altSampleType) sSample_alt = sSample+"_"+current->systematicName_;    // get correct Sample name if alternative Sample Type
         //import 1D Hists
         for (auto const &distr_:vecDistr){
            TString loc;
            loc=distr_.path+distr_.name;
            TH1F* tempHist=current->histReader_->read<TH1F>(loc+"/"+sSample_alt);
            
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
            
            TH2F rebinHist;   //rebin hist
            if (distr_.binEdgesX.empty()) {   // use same bin width
               rebinHist = hist::rebinned(*tempHist,distr_.xMin,distr_.xMax,distr_.nBinsX,distr_.yMin,distr_.yMax,distr_.nBinsY);
            }
            else{
               rebinHist = hist::rebinned(*tempHist,distr_.binEdgesX,distr_.binEdgesY);
            }
            
            TH1F trafoHist = hist::histTrafo_2D(&rebinHist);
            
            if (!current->hasRootFile_) trafoHist.Scale(current->sf_);
            current->hists_.addFilledHist(loc,sSample,trafoHist);
         }
         // ~auto stop = high_resolution_clock::now();
         // ~auto duration = duration_cast<milliseconds>(stop - start);
         // ~std::cout<<sSample<<"   "<<duration.count()<<std::endl;
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

// get up and down shift for set of systematics
std::pair<TH1F*,TH1F*> distributionsplotting::getTotalSyst(TH1F* const &nominal, std::vector<systHists*> const &systHists_vec, TString const loc, TString const sample, bool envelope){
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   TH1F* tempSys;
   TH1F tempShift;
   TH1F* nominal_ttbar;
   TH1F* nominal_ttbar2L;
   for (auto &current : systHists_vec){
      if (current->systematic_.type() == Systematic::nominal || current->systematic_.type() == Systematic::met40Cut){  //Store sum of ttBar for syst. with alt. samples
         if (sample == "") {
            nominal_ttbar = current->hists_.getSummedHist(loc,current->datasets_ttBar_);
            nominal_ttbar2L = current->hists_.getSummedHist(loc,current->datasets_ttBar2L_);   //for PDF unc.
         }
         else {
            nominal_ttbar = current->hists_.getHistogram(loc,sample);
            nominal_ttbar2L = current->hists_.getHistogram(loc,sample);
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
            
      if (current->altSampleType || current->onlyTTbar) tempShift= phys::getSystShift(*nominal_ttbar,*tempSys);  // Choose correct reference
      else if (current->onlyTTbar2L) tempShift= phys::getSystShift(*nominal_ttbar2L,*tempSys);  // for PDF unc.
      else tempShift= phys::getSystShift(*nominal,*tempSys);
      
      if(envelope){
         for (int i=0; i<=tempShift.GetNbinsX(); i++){
            float content = tempShift.GetBinContent(i);
            if (content>0 && (content>hist_shiftUP->GetBinContent(i))) hist_shiftUP->SetBinContent(i,content);
            else if (content<0 && content<hist_shiftDOWN->GetBinContent(i)) hist_shiftDOWN->SetBinContent(i,abs(content));
         }
      }
      else{
         for (int i=0; i<=tempShift.GetNbinsX(); i++){
            float content = tempShift.GetBinContent(i);
            if (content>0) hist_shiftUP->SetBinContent(i,hist_shiftUP->GetBinContent(i)+content*content);
            else hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
         }
      }
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
void distributionsplotting::printTotalYields(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples_merged){
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
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat);
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      float mcYield = mc_total->GetBinContent(6);
      float mcYield_down = syst.first->GetBinContent(6);
      float mcYield_up = syst.second->GetBinContent(6);
      mcYield_total += mcYield;
      mcYield_total_down += mcYield_down;
      mcYield_total_up += mcYield_up;
      bool isMCtotal = false;
      for (TString sample:outputSamples){
         isMCtotal = (sample=="MC");
         TH1F* temp_hist=hs->getHistogram("cutflow/"+cat,sample);
         // ~TH1F* temp_hist=hs->getHistogram("baseline/"+cat+"/nBjets",sample);
         float sampleYield = temp_hist->GetBinContent(6);
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
   for (TString cat:{"ee","emu","mumu"}){
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      // ~TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"DrellYan_NLO");
      for (int i=1; i<(systHists_vec.size()); i++){
         if ((i+1)<systHists_vec.size()){    // check if still room for up and down shift
            if (systHists_vec[i]->systematic_.type_str()==systHists_vec[i+1]->systematic_.type_str()){     // check if up and down shift
               std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i],systHists_vec[i+1]};
               std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat);
               // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"DrellYan_NLO");
               printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(6),syst.second->GetBinContent(6),mc_total->GetBinContent(6));
               i++;
               continue;
            }
         }
         std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i]};   // print unc. without up and down shift
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat);
         // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat,"DrellYan_NLO");
         printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(6),syst.second->GetBinContent(6),mc_total->GetBinContent(6));
      }
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat);     // print total uncertainty
      // ~std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat,"DrellYan_NLO");     // print total uncertainty
      std::cout<<"\\hline"<<std::endl;
      printUnc("TOTAL",syst.first->GetBinContent(6),syst.second->GetBinContent(6),mc_total->GetBinContent(6));
   }
}

// print uncertainty shift for individual samples
void distributionsplotting::printShiftBySample(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples){
   std::cout<<std::endl<<"-----------------Uncertainties per sample-------------------"<<std::endl;
   for (TString cat:{"ee","emu","mumu"}){
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      float mcYield = mc_total->GetBinContent(6);
      for (TString sample:mcSamples){
         TH1F* temp_hist=hs->getHistogram("cutflow/"+cat,sample);
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(temp_hist,systHists_vec,"cutflow/"+cat,sample);
         float sampleYield = temp_hist->GetBinContent(6);
         float sampleYield_down = syst.first->GetBinContent(6);
         float sampleYield_up = syst.second->GetBinContent(6);
         std::cout<<sample<<"   "<<sampleYield_down<<"   "<<sampleYield_up<<std::endl;
         std::cout<<sample<<"   "<<sampleYield_down/sampleYield*100<<"   "<<sampleYield_up/sampleYield*100<<std::endl;
      }
   }
}

void distributionsplotting::fixAxis2D(TAxis* axis){     // currently only works for 2D SR with 6x11 bins
   axis->SetNdivisions(12);
   axis->ChangeLabel(13,-1,-1,-1,-1,-1," ");
   for (int i=0; i<3; i++){
      axis->ChangeLabel(i*4+1,-1,-1,-1,-1,-1,"  ");
      axis->ChangeLabel(i*4+2,-1,-1,-1,-1,-1,"100");
      axis->ChangeLabel(i*4+3,-1,-1,-1,-1,-1," ");
      axis->ChangeLabel(i*4+4,-1,-1,-1,-1,-1,"300");
   }
   // ~axis->SetNdivisions(24);
   // ~axis->ChangeLabel(25,-1,-1,-1,-1,-1," ");
   // ~for (int i=0; i<6; i++){
      // ~axis->ChangeLabel(i*4+1,-1,-1,-1,-1,-1,"  ");
      // ~axis->ChangeLabel(i*4+2,-1,-1,-1,-1,-1,"100");
      // ~axis->ChangeLabel(i*4+3,-1,-1,-1,-1,-1," ");
      // ~axis->ChangeLabel(i*4+4,-1,-1,-1,-1,-1,"300");
   // ~}
}

void distributionsplotting::drawVertLines2D(float const &lowerEnd, float const &upperEnd,bool text){
   TLine * aline = new TLine();
   TLatex * atext = new TLatex();
   atext->SetTextSize(0.025);
   aline->SetLineWidth(2);
   for (int i=1; i<3; i++){
      aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
   }
   if (text){
      atext->DrawLatex(150,1.7*lowerEnd,"|#Delta#phi|<0.35");
      atext->DrawLatex(500,1.7*lowerEnd,"0.7<|#Delta#phi|<1.4");
      atext->DrawLatex(900,1.7*lowerEnd,"1.4<|#Delta#phi|");
   }
   // ~for (int i=1; i<6; i++){
      // ~aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
   // ~}
   // ~if (text){
      // ~atext->DrawLatex(100,1.7*lowerEnd,"|#Delta#phi|<0.35");
      // ~atext->DrawLatex(450,1.7*lowerEnd,"0.35<|#Delta#phi|<0.7");
      // ~atext->DrawLatex(850,1.7*lowerEnd,"0.7<|#Delta#phi|<1.05");
      // ~atext->DrawLatex(1250,1.7*lowerEnd,"1.05<|#Delta#phi|<1.4");
      // ~atext->DrawLatex(1650,1.7*lowerEnd,"1.4<|#Delta#phi|<2.27");
      // ~atext->DrawLatex(2050,1.7*lowerEnd,"2.27<|#Delta#phi|");
   // ~}
}

void distributionsplotting::plotHistograms(TString const &sPresel, TString const &sVar, hist::Histograms<TH1F>* const &hs, std::vector<TString> const &mcSamples_merged, 
                     std::map<const TString,Color_t> const & colormap, std::vector<systHists*> const &systHists_vec, io::RootFileSaver const &saver, bool plotStatUncExpData, bool is2D){
   TCanvas can;
   can.SetLogy();
   gfx::SplitCan sp_can;
   sp_can.pU_.SetLogy();
   sp_can.pU_.cd();
   TString loc;
   loc=sPresel+sVar;
   THStack st_mc=hs->getStack(loc,mcSamples_merged,colormap);
   gfx::LegendEntries le=hs->getLegendEntries();
   
   //systematics
   TH1F* hist_mc = hs->getHistogram(loc,{"MC"});
   std::pair<TH1F*,TH1F*> syst = getTotalSyst(hist_mc,systHists_vec,loc);
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
   
   systGraph.SetFillStyle(3001);
   systGraph.Draw("same 2");    // draw syst.
   
   if (is2D){
      fixAxis2D(st_mc.GetXaxis());
      drawVertLines2D(st_mc.GetHistogram()->GetMinimum(),st_mc.GetHistogram()->GetMaximum(),true);
   }
   
   // data plotting part
   auto hist_data = hs->getHistogram(loc,{"data"});
   hist_data->SetLineColor(kBlack);
   hist_data->SetMarkerSize(0.5);
   bool plotData = false;
   if(sPresel.Contains("GOF2D") || !(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep")|| sVar.Contains("C_em")|| sVar.Contains("Met"))) {
      hist_data->Draw("same");
      le.append(*hist_data,"data","lep");
      plotData = true;
   }
   
   //auto hists_BSM=hs->getHistograms(loc,{"TTbar_diLepton"});     //placeholder
   //auto BSM_legend=hs->getLegendEntries();    //placeholder
   //if (cfg.year_int==1){    //Plot BSM in 2016
   //   hists_BSM=hs->getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
   //   for (auto const &h: hists_BSM) {
   //      h->Draw("same hist");
   //      // ~h->SetLineWidth(4);
   //   }
   //   le+=hs->getLegendEntries();
   //   BSM_legend=hs->getLegendEntries();
   //}
   
   // redraw axis, draw label and legend
   auto hists_SM=hs->getHistograms(loc,{"TTbar_diLepton"});
   TH1F axis=*(hists_SM[0]);
   axis.SetStats(0);
   axis.SetYTitle("Events/Bin");
   axis.Draw("same axis");
   TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
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
      ratio_mc.GetXaxis()->SetBinLabel(4,"met");
      ratio_mc.GetXaxis()->SetBinLabel(5,"btag");
      ratio_mc.GetXaxis()->SetBinLabel(6,"triggerSF");
      // ~ratio_mc.GetXaxis()->SetBinLabel(7,"(addLepton veto)");
      ratio_mc.GetXaxis()->SetRangeUser(0.5,6.5);
      ratio_mc.GetXaxis()->SetLabelOffset(0.03);
   }
   
   ratio_mc.GetYaxis()->SetTitleOffset(0.45);
   ratio_mc.SetStats(0);
   ratio.SetLineColor(kBlack);
   // ~ratio_mc.SetMaximum(1.04);
   // ~ratio_mc.SetMinimum(0.9);
   ratio_mc.SetMaximum(1.35);
   ratio_mc.SetMinimum(0.65);
   totalUncGraphRatio.SetFillColor(kGray);
   totalUncGraphRatio.SetLineWidth(0);
   expDataStatRatio.SetFillColor(kBlue);
   systGraphRatio.SetFillColor(kGray+3);
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
   // ~ratio.SetMarkerSize(4);
   gPad->RedrawAxis("g");
   if(plotData) ratio.Draw("pe1 same");
   
   if(is2D){
      fixAxis2D(ratio_mc.GetXaxis());
      drawVertLines2D(ratio_mc.GetMinimum(),ratio_mc.GetMaximum(),false);
   }
   
   gfx::LegendEntries le_low;
   le_low.append(totalUncGraphRatio,"#sigma_{tot.}","f");
   le_low.append(systGraphRatio,"#sigma_{syst.}","f");
   if(plotStatUncExpData) le_low.append(expDataStatRatio,"#sigma_{stat.(exp. data)}","f");
   TLegend leg_low=le_low.buildLegend(.5,.8,0.8,0.95,2);
   leg_low.Draw();
      
   saver.save(sp_can,loc,false,true);
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
