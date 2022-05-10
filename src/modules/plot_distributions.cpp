//Script to (re-)plot distributions from distributions.cpp 

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/systematics.hpp"
#include "tools/util.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraphAsymmErrors.h>
#include <iomanip>
#include <chrono>
#include <algorithm>
using namespace std::chrono;

Config const &cfg=Config::get();

// struct to combine name of distributions and intended binning
struct distr {
   TString path;
   TString name;
   float xMin;
   float xMax;
   int nBins;
   std::vector<double> binEdges = {};
};

struct distr2D {
   TString path;
   TString name;
   float xMin;
   float xMax;
   int nBinsX;
   float yMin;
   float yMax;
   int nBinsY;
   std::vector<float> binEdgesX = {};
   std::vector<float> binEdgesY = {};
};

// small class for connecting syst., reader and histograms
class systHists
{
   public:
      systHists(TString const &systematicName, TString filePath, TString const &histPath, std::vector<TString> const &datasets, std::vector<TString> const &datasets_ttBar, std::vector<TString> const &datasets_ttBar2L):
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
            if(cfg.systUncFactor.find(systematic_.type_str()) != cfg.systUncFactor.end()){
               float unc = cfg.systUncFactor.at(systematic_.type_str()).first;
               sf_ = (systematic_.variation() == Systematic::up)? 1.0+unc : 1.0-unc;
               datasets_SF = cfg.systUncFactor.at(systematic_.type_str()).second;
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
         if (systematic_.type() == Systematic::bFrag || systematic_.type() == Systematic::bSemilep || systematic_.type() == Systematic::meFacScale || systematic_.type() == Systematic::meRenScale){
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
      
      void openFile(){
         histReader_ = new io::RootFileReader(filePath_,histPath_,true,true);
      }
      
      void combineChannel(){
         hists_.combineChannel("/all",{"/ee","/mumu","/emu"});
      }
      
      Systematic::Systematic systematic_;
      TString filePath_;
      TString histPath_;
      io::RootFileReader* histReader_;
      hist::Histograms<TH1F> hists_;
      std::vector<TString> datasets_;
      std::vector<TString> datasets_ttBar_;
      std::vector<TString> datasets_ttBar2L_;
      TString const systematicName_;
      float sf_ = 1.0;
      std::vector<std::string> datasets_SF;   //datasets for which the SF is applied when using the syst.
      bool hasRootFile_ = true;
      bool altSampleType = false;
      bool onlyTTbar = false;
      bool onlyTTbar2L = false;
};

//Import hists for Nominal and systematics
void importHists(std::vector<systHists*> &systHists_vec, std::vector<TString> const &samplesToPlot, std::vector<TString> const &mcSamples,
                     std::vector<distr> const &vecDistr, std::vector<distr2D> const &vecDistr2D) {
   for (systHists* &current : systHists_vec){
      std::vector<TString> inputSamples = current->datasets_;
      current->openFile();
      
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

/*
//Import 2D hists for Nominal and systematics (onyl signal variable currently)
void importHists2D(std::vector<systHists*> &systHists_vec, std::vector<TString> const &samplesToPlot, std::vector<TString> const &mcSamples,
                     std::map<TString,std::vector<TString>> const &msPresel_vVars) {
                        
   std::vector<float> binedges_x = {0, 20, 40, 60, 80, 100, 120, 140, 160, 195, 230, 400};
   std::vector<float> binedges_y = {0, 0.35, 0.7, 1.05, 1.4, 2.27, 3.14};
   for (systHists* &current : systHists_vec){
      std::vector<TString> inputSamples;
      if (current->systematic_.type() == Systematic::nominal) inputSamples = samplesToPlot;
      else inputSamples = mcSamples;
      
      for (TString sSample : inputSamples){
         current->hists_.setCurrentSample(sSample);
         for (auto const &sPresel_vVars:msPresel_vVars){
            TString const &sPresel=sPresel_vVars.first;
            for (TString sVar:sPresel_vVars.second){
               TString loc;
               loc=sPresel+sVar;
               TH2F* tempHist=current->histReader_->read<TH2F>(loc+"/"+sSample);
               TH1F tempHist1D=HistTrafo_2D(tempHist,binedges_x,binedges_y);
               current->hists_.addFilledHist(loc,sSample,tempHist1D);
            }
            // ~TH1F* tempHist=current->histReader_->read<TH1F>(sPresel+"ee/"+sSample);
            // ~tempHist->Add(current->histReader_->read<TH1F>(sPresel+"emu/"+sSample));
            // ~tempHist->Add(current->histReader_->read<TH1F>(sPresel+"mumu/"+sSample));
            // ~current->hists_.addFilledHist(sPresel+"all",sSample,*(tempHist));
         }
      }
   }
}
*/

void add_Categories(TString const path, io::RootFileReader const &reader_hist, TH1F &out_hist) {   //Function to add the three different categories
   TH1F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH1F*) reader_hist.read<TH1F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH1F) *(reader_hist.read<TH1F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

// get up and down shift for set of systematics
std::pair<TH1F*,TH1F*> getTotalSyst(TH1F* const &nominal, std::vector<systHists*> const &systHists_vec, TString const loc, TString const sample="", bool envelope=false){
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

void printUnc(TString name, const float &down, const float &up, const float &nominal){
   name = name.ReplaceAll("ELECTRON_SCALESMEARING","EL_SCALESMEARING");
   TString out = TString::Format(" %s & $%.1f(%.1f)$ & $%.1f(%.1f)$\\\\\n",name.ReplaceAll("_","\\_").Data(),down,down/nominal*100,up,up/nominal*100);
   std::cout<<out;
}

// print total yields per contribution
void printTotalYields(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples_merged){
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
         float sampleYield = temp_hist->GetBinContent(6);
         totalMap[sample.ReplaceAll("_","\\_")] += sampleYield;
         if(!isMCtotal){
            std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<sampleYield<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
         }
         else{
            TString yieldUnc = TString::Format("$%.1f^{+%.1f}_{-%.1f}$",sampleYield,mcYield_up,mcYield_down);
            std::cout<<"\\hline"<<std::endl;
            std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<yieldUnc<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
         }
         // ~if (sample=="MC") std::cout<<mcYield_down<<"   "<<mcYield_up<<std::endl;
         // ~if (sample=="MC") std::cout<<mcYield_down/mcYield*100<<"   "<<mcYield_up/mcYield*100<<std::endl;
      }
   }
   // ~//print channels combined
   // ~std::cout<<"----------------"<<"combined"<<"-----------------------"<<std::endl;
   // ~for (TString sample:outputSamples){
      // ~if(sample=="MC"){
         // ~TString yieldUnc = TString::Format("$%.1f^{+%.1f}_{-%.1f}$",mcYield_total,mcYield_total_up,mcYield_total_down);
         // ~std::cout<<"\\hline"<<std::endl;
         // ~std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<yieldUnc<<"&"<<std::setprecision(1)<<mcYield_total/mcYield_total*100<<"\\\\"<<std::endl;
      // ~}
      // ~else{
         // ~std::cout<<std::fixed<<sample.ReplaceAll("_","\\_")<<"&"<<totalMap[sample]<<"&"<<std::setprecision(1)<<totalMap[sample]/mcYield_total*100<<"\\\\"<<std::endl;
      // ~}
   // ~}
   
}

// print breakdown of syst uncertainties
void printUncBreakDown(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples) {
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
void printShiftBySample(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples){
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

void fixAxis2D(TAxis* axis){     // currently only works for 2D SR with 6x11 bins
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

void drawVertLines2D(float const &lowerEnd, float const &upperEnd,bool text){
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

void plotHistograms(TString const &sPresel, TString const &sVar, hist::Histograms<TH1F>* const &hs, std::vector<TString> const &mcSamples_merged, 
                     std::map<const TString,Color_t> const & colormap, std::vector<systHists*> const &systHists_vec, io::RootFileSaver const &saver, bool plotStatUncExpData=false, bool is2D=false){
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
      ratio_mc.GetXaxis()->SetBinLabel(1,"DiLepton");
      ratio_mc.GetXaxis()->SetBinLabel(2,"mll");
      ratio_mc.GetXaxis()->SetBinLabel(3,"jets");
      ratio_mc.GetXaxis()->SetBinLabel(4,"met");
      ratio_mc.GetXaxis()->SetBinLabel(5,"btag");
      ratio_mc.GetXaxis()->SetBinLabel(6,"ScaleFactors");
      // ~ratio_mc.GetXaxis()->SetBinLabel(7,"(addLepton veto)");
      ratio_mc.GetXaxis()->SetRangeUser(0.5,6.5);
   }
   
   ratio_mc.GetYaxis()->SetTitleOffset(0.45);
   ratio_mc.SetStats(0);
   ratio.SetLineColor(kBlack);
   // ~ratio_mc.SetMaximum(1.04);
   // ~ratio_mc.SetMinimum(0.9);
   ratio_mc.SetMaximum(1.35);
   ratio_mc.SetMinimum(0.65);
   totalUncGraphRatio.SetFillColor(kGray);
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
   TLegend leg_low=le_low.buildLegend(.2,.8,0.5,0.95,2);
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
   
extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions");
   
   std::vector<TString> mcSamples={};
   std::vector<TString> dataSamples={};
   std::vector<TString> ttbarSamples={};
   std::vector<TString> signalSamples={};
   switch(cfg.year_int){
      case(3): //2018
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_ee","DrellYan_mumu","DrellYan_tautau","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
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
   std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
   
   // JESTotal
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_IND_UP","MTOP_IND_DOWN","CR_ENVELOPE_IND_UP","CR_ENVELOPE_IND_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // JESTotal without PDF
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_IND_UP","MTOP_IND_DOWN","CR_ENVELOPE_IND_UP","CR_ENVELOPE_IND_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   
   // Split JES
   // ~std::vector<TString> systToPlot = {"Nominal","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_IND_UP","MTOP_IND_DOWN","CR_ENVELOPE_IND_UP","CR_ENVELOPE_IND_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorQCD_UP","JESFlavorQCD_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN"};
   
   // Regrouped JES
   std::vector<TString> systToPlot = {"Nominal","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_IND_UP","MTOP_IND_DOWN","CR_ENVELOPE_IND_UP","CR_ENVELOPE_IND_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN","L1PREFIRING_UP","L1PREFIRING_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESBBEC1_UP","JESBBEC1_DOWN"};
   
   // Only split JES
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN","JESAbsoluteScale_UP","JESAbsoluteScale_DOWN","JESAbsoluteStat_UP","JESAbsoluteStat_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESFragmentation_UP","JESFragmentation_DOWN","JESPileUpDataMC_UP","JESPileUpDataMC_DOWN","JESPileUpPtBB_UP","JESPileUpPtBB_DOWN","JESPileUpPtEC1_UP","JESPileUpPtEC1_DOWN","JESPileUpPtRef_UP","JESPileUpPtRef_DOWN","JESRelativeBal_UP","JESRelativeBal_DOWN","JESRelativeFSR_UP","JESRelativeFSR_DOWN","JESRelativeJEREC1_UP","JESRelativeJEREC1_DOWN","JESRelativePtBB_UP","JESRelativePtBB_DOWN","JESRelativePtEC1_UP","JESRelativePtEC1_DOWN","JESRelativeSample_UP","JESRelativeSample_DOWN","JESRelativeStatEC_UP","JESRelativeStatEC_DOWN","JESRelativeStatFSR_UP","JESRelativeStatFSR_DOWN","JESSinglePionECAL_UP","JESSinglePionECAL_DOWN","JESSinglePionHCAL_UP","JESSinglePionHCAL_DOWN","JESTimePtEta_UP","JESTimePtEta_DOWN"};
   
   // Only regrouped JES
   // ~std::vector<TString> systToPlot = {"Nominal","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESFlavorQCDreg_UP","JESFlavorQCDreg_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESBBEC1_UP","JESBBEC1_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESRelativeBalreg_UP","JESRelativeBalreg_DOWN","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN","JESRelativeSampleYear_UP","JESRelativeSampleYear_DOWN","JESAbsoluteYear_UP","JESAbsoluteYear_DOWN","JESAbsolute_UP","JESAbsolute_DOWN","JESBBEC1Year_UP","JESBBEC1Year_DOWN","JESBBEC1_UP","JESBBEC1_DOWN"};
   
   // ~std::vector<TString> systToPlot = {"Nominal","JESFlavorPureGluon_UP","JESFlavorPureGluon_DOWN","JESFlavorPureQuark_UP","JESFlavorPureQuark_DOWN","JESFlavorPureCharm_UP","JESFlavorPureCharm_DOWN","JESFlavorPureBottom_UP","JESFlavorPureBottom_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESFlavorRealistic_UP","JESFlavorRealistic_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","ELECTRON_ID_UP","ELECTRON_ID_DOWN","MUON_ID_UP","MUON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ISO_UP","MUON_ISO_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MUON_ISO_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP","PU_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PU_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","UNCLUSTERED_UP","UNCLUSTERED_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","UNCLUSTERED_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","MEFACSCALE_UP","MEFACSCALE_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","MERENSCALE_UP","MERENSCALE_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP","PSFSRSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_UP","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BFRAG_UP","BFRAG_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BFRAG_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR_ENVELOPE_IND_UP","CR_ENVELOPE_IND_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR_ENVELOPE_IND_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","CR1","CR2","ERDON"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_IND_UP","MTOP_IND_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MTOP_IND_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_UP","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","UETUNE_UP","UETUNE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","UETUNE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PDF_50_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN"};
   //std::vector<TString> systToPlot = {"Nominal","LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteMPFBias_UP","JESAbsoluteMPFBias_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESRelativeSample_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESAbsoluteYear_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal"};
   
   // ~for (int i=1;i<=50;i++){
       // ~systToPlot.push_back(TString::Format("PDF_%i_UP",i));
       // ~systToPlot.push_back(TString::Format("PDF_%i_DOWN",i));
   // ~}
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         /*
         // ~vecDistr.push_back({selection+channel,"Lep_e_pt",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep_mu_pt",0.,600.,50});
         vecDistr.push_back({selection+channel,"Lep1_pt",0.,360.,30});
         vecDistr.push_back({selection+channel,"Lep2_pt",0.,300.,25});
         // ~vecDistr.push_back({selection+channel,"MET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET",0.,500.,7,{0,40,70,110,170,260,370,500}});
         // ~vecDistr.push_back({selection+channel,"PuppiMET",0.,500.,7,{0,20,40,55,70,90,110,140,170,215,260,315,370,435,500}});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0,3.2,40});
         // ~vecDistr.push_back({selection+channel,"met1000",0.,1000.,50});
         vecDistr.push_back({selection+channel,"mLL",0,200,25});
         // ~vecDistr.push_back({selection+channel,"pTsumlep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"sumpTlep",0.,600.,30});
         vecDistr.push_back({selection+channel,"pTbJet",0.,600.,30});
         vecDistr.push_back({selection+channel,"Jet1_pt",0.,600.,30});
         vecDistr.push_back({selection+channel,"Jet2_pt",0.,400.,20});
         // ~vecDistr.push_back({selection+channel,"dPhiMETnearJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETleadJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETlead2Jet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi",0.,3.2,8});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi",0.,3.2,16});
         // ~vecDistr.push_back({selection+channel,"COSdphi_metNearLep",-1.,1,50});
         // ~vecDistr.push_back({selection+channel,"SINdphi_metNearLep",0.,1,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETbJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1bJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_bJetLep2",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_bJetnearLep",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_b1b2",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dR_b1b2",0.,5.,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLep1",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLep2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dphi_metLepsum",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1Lep2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dR_Lep1Lep2",0.,5.,100});
         vecDistr.push_back({selection+channel,"nJets",1.5,9.5,8});
         vecDistr.push_back({selection+channel,"nBjets",-0.5,4.5,5});
         // ~vecDistr.push_back({selection+channel,"MT2",0.,200.,50});
         // ~vecDistr.push_back({selection+channel,"MT",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"mt_MetLep2",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"mt_MetNextLep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"conMt_Lep1Lep2",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"ST",0.,1500.,50});
         // ~vecDistr.push_back({selection+channel,"HT",0.,2500.,50});
         // ~vecDistr.push_back({selection+channel,"sum_STHT",0.,4000.,50});
         // ~vecDistr.push_back({selection+channel,"sum_mlb",0.,3000.,50});
         // ~vecDistr.push_back({selection+channel,"METunc_Puppi",0.,50.,50});
         // ~vecDistr.push_back({selection+channel,"n_Interactions",0.,100.,100});
         // ~vecDistr.push_back({selection+channel,"Lep1_flavor",0.5,2.5,2});
         // ~vecDistr.push_back({selection+channel,"Lep2_flavor",0.5,2.5,2});
         // ~vecDistr.push_back({selection+channel,"Lep1_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"Lep2_phi",-3.2,3.2,50});
         vecDistr.push_back({selection+channel,"Lep1_eta",-2.5,2.5,25});
         vecDistr.push_back({selection+channel,"Lep2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Lep1_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep2_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Jet1_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_phi",-3.2,3.2,50});
         vecDistr.push_back({selection+channel,"Jet1_eta",-2.5,2.5,25});
         vecDistr.push_back({selection+channel,"Jet2_eta",-2.5,2.5,25});
         // ~vecDistr.push_back({selection+channel,"Jet1_E",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_E",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETfarJet",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiJet1Jet2",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"METsig",0.,200.,50});
         // ~vecDistr.push_back({selection+channel,"MHT",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"looseLeptonVeto",-0.5,1.5,2});
         // ~vecDistr.push_back({selection+channel,"dPhiMETnearJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETfarJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETleadJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETlead2Jet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiMETbJet_Puppi",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"dPhiLep1Jet1",0.,3.2,50});
         // ~vecDistr.push_back({selection+channel,"PFMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"CaloMET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"CaloMET_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_allJet",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"vecsum_pT_l1l2_allJet",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"mass_l1l2_allJet",0.,3000.,50});
         // ~vecDistr.push_back({selection+channel,"ratio_vecsumpTlep_vecsumpTjet",0.,20.,50});
         // ~vecDistr.push_back({selection+channel,"mjj",0.,2000.,50});
         // ~vecDistr.push_back({selection+channel,"Lep1_pt*cos(Lep1_phi)",-250,250,25});
         // ~vecDistr.push_back({"baseline_GOF2D"+channel,"PuppiMET_xy*sin(PuppiMET_xy_phi)_VS_MET_xy*sin(MET_xy_phi)",0.5,36.5,36});
         */
      }
   }
   
   // 2D plots
   std::vector<distr2D> vecDistr2D;
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_Puppi",0.,400.,6,0.,3.14,3,{0,40,80,120,160,230,400},{0,0.7,1.4,3.14}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_Puppi",0.,400.,6,0.,3.14,3,{0,20,40,60,80,100,120,140,160,195,230,400},{0,0.35,0.7,1.05,1.4,2.27,3.14}});
         // ~vecDistr2D.push_back({selection+channel,"2d_MetVSdPhiMetNearLep_DNN",0.,400.,10,0.,3.2,8});
      }
   }
   
   // Setup systematics
   std::vector<systHists*> systHists_vec;
   for (TString syst : systToPlot){
      systHists* temp = new systHists(syst,TString::Format("multiHists/%s/histograms_merged_%s.root",syst.Data(),cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal" || syst=="met40Cut")? samplesToPlot : mcSamples,(syst=="CR_ENVELOPE_UP" || syst=="CR_ENVELOPE_DOWN" || syst=="MTOP_DOWN" || syst=="MTOP_UP")? signalSamples : ttbarSamples, signalSamples);
      systHists_vec.push_back(temp);
   }
   
   std::vector<TString> mcSamples_merged={};
   
   for (auto const & distr_ : vecDistr){
      // Import hists
      importHists(systHists_vec,samplesToPlot,mcSamples,{distr_},{});
      
      // Define hist collection for nominal
      hist::Histograms<TH1F>* hs;
      hs = &(systHists_vec[0]->hists_);
      
      // combine different samples to improve readability
      hs->combineSamples("Diboson",{"WW","WZ","ZZ"});
      hs->combineSamples("DrellYan_comb",{"DrellYan_NLO","DrellYan_M10to50"});
      // ~hs->combineSamples("DrellYan_comb",{"DrellYan_ee","DrellYan_mumu","DrellYan_tautau","DrellYan_M10to50"});
      hs->combineSamples("ttZ",{"ttZ_2L","ttZ_QQ"});
      hs->combineSamples("ttW/Z",{"ttW","ttZ"});
      hs->combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
      switch(cfg.year_int){
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
      std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
      
      plotHistograms(distr_.path,distr_.name,hs,mcSamples_merged,colormap,systHists_vec,saver,false,false);
      
      if (distr_.path.Contains("/emu") || distr_.name == "emu"){  //plot all channels combined
         for (auto &current : systHists_vec){
            current->combineChannel();
         }
         TString combinedPath(distr_.path);
         combinedPath.ReplaceAll("/emu","/all");
         if (distr_.name == "emu") {   // fot cutflow plot
            plotHistograms(combinedPath,"all",hs,mcSamples_merged,colormap,systHists_vec,saver,false,false);
         }
         else {
            plotHistograms(combinedPath,distr_.name,hs,mcSamples_merged,colormap,systHists_vec,saver,false,false);
         }
      }
   }
         
   // ~}
   
   hist::Histograms<TH1F>* hs;
   hs = &(systHists_vec[0]->hists_);
   // Print total yields
   printTotalYields(hs,systHists_vec,mcSamples_merged);
   
   // Print uncertainties
   printUncBreakDown(hs,systHists_vec,mcSamples);
      
   // Print shift per sample
   printShiftBySample(hs,systHists_vec,mcSamples);
   
   /*
   //Plot 2D SignalRegion Plot (reco level)
   for (auto &sPresel_vVars:msPresel_vVars_2D){
      TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         TString loc=sPresel+sVar;
         THStack st_mc=hs->getStack(loc,MCsamples,colormap);
         st_mc.SetMaximum(st_mc.GetMaximum()*3);
         st_mc.SetTitle(";p_{#scale[.8]{T}}^{#scale[.8]{miss}}(GeV);Events/Bin");
         st_mc.Draw();
         gfx::LegendEntries le=hs->getLegendEntries();
         
         //systematics
         TH1F* hist_mc = hs->getHistogram(loc,{"MC"});
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(hist_mc,systHists_vec,loc);
         
         TGraphAsymmErrors systGraph = hist::getErrorGraph(syst.first,syst.second,hist_mc,true);
         
         systGraph.SetFillStyle(3001);
         systGraph.Draw("same 2");
         
         st_mc.GetXaxis()->SetNdivisions(24);
         st_mc.GetXaxis()->ChangeLabel(25,-1,-1,-1,-1,-1," ");
         for (int i=0; i<=5; i++){
            st_mc.GetXaxis()->ChangeLabel(i*4+1,-1,-1,-1,-1,-1,"  ");
            st_mc.GetXaxis()->ChangeLabel(i*4+2,-1,-1,-1,-1,-1,"100");
            st_mc.GetXaxis()->ChangeLabel(i*4+3,-1,-1,-1,-1,-1," ");
            st_mc.GetXaxis()->ChangeLabel(i*4+4,-1,-1,-1,-1,-1,"300");
         }
         
         // ~auto hists=hs->getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
         // ~for (auto const &h: hists) {
            // ~h->SetLineWidth(4);
            // ~h->Draw("same hist");
         // ~}
         // ~le+=hs->getLegendEntries();
         
         TLine * aline = new TLine();
         TLatex * atext = new TLatex();
         atext->SetTextSize(0.03);
         aline->SetLineWidth(2);
         float upperEnd=st_mc.GetHistogram()->GetMaximum();
         float lowerEnd=st_mc.GetHistogram()->GetMinimum();
         for (int i=1; i<=5; i++){
            aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
         }
         atext->DrawLatex(100,0.2*upperEnd,"|#Delta#phi|<0.35");
         atext->DrawLatex(475,0.2*upperEnd,"0.35<|#Delta#phi|<0.7");
         atext->DrawLatex(875,0.2*upperEnd,"0.7<|#Delta#phi|<1.05");
         atext->DrawLatex(1275,0.2*upperEnd,"1.05<|#Delta#phi|<1.4");
         atext->DrawLatex(1675,0.2*upperEnd,"1.4<|#Delta#phi|<2.27");
         atext->DrawLatex(2075,0.2*upperEnd,"2.27<|#Delta#phi|");
         TString cat;
         if (loc.Contains("ee")) cat="ee";
         else if (loc.Contains("emu")) cat="e#mu";
         else if (loc.Contains("mumu")) cat="#mu#mu";
         else if (loc.Contains("all")) cat="all";
         atext->DrawLatex(50,0.3*upperEnd,cat);
         
         TLegend leg=le.buildLegend(gPad->GetLeftMargin(),.9,1-gPad->GetRightMargin(),1-gPad->GetTopMargin(),10);
         leg.SetTextSize(0.025);
         leg.SetBorderSize(2);
         leg.SetShadowColor(0);
         leg.SetFillColor(kWhite);
         leg.SetFillStyle(1001);
         leg.Draw();
      
         saver.save(can,loc,true,false,true);
      }
   }
   */
}
