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
};

// small class for connecting syst., reader and histograms
class systHists
{
   public:
      systHists(TString const &systematicName, TString filePath, TString const &histPath, std::vector<TString> const &datasets, std::vector<TString> const &datasets_ttBar):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(hist::Histograms<TH1F>(datasets)),
      systematicName_(systematicName),
      datasets_(datasets),
      datasets_ttBar_(datasets_ttBar)
      {
         const std::vector<Systematic::Type> typeVec = Systematic::fileIndependentTypes;
         if (std::find(typeVec.begin(), typeVec.end(), systematic_.type()) == typeVec.end()){
            histReader_ = new io::RootFileReader(filePath,histPath);
         }
         else {      //currently only used for lumi unc.
            hasRootFile_ = false;
            filePath.ReplaceAll(systematic_.name(),"Nominal");
            histReader_ = new io::RootFileReader(filePath,histPath);
            if(cfg.systUncFactor.find(systematic_.type_str()) != cfg.systUncFactor.end()){
               float unc = cfg.systUncFactor.at(systematic_.type_str());
               sf_ = (systematic_.variation() == Systematic::up)? 1.0+unc : 1.0-unc;
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
         
      }
      
      Systematic::Systematic systematic_;
      io::RootFileReader* histReader_;
      hist::Histograms<TH1F> hists_;
      std::vector<TString> datasets_;
      std::vector<TString> datasets_ttBar_;
      TString const systematicName_;
      float sf_ = 1.0;
      bool hasRootFile_ = true;
      bool altSampleType = false;
      bool onlyTTbar = false;
};

// transform 2D hist to 1D hist
TH1F HistTrafo_2D(TH2F* const &hist2D, std::vector<float> binedges_x, std::vector<float> binedges_y){   
   *hist2D=hist::rebinned(*hist2D, binedges_x, binedges_y);
   
   int numBins_x=binedges_x.size()-1;
   int numBins_y=binedges_y.size()-1;
   int numBins=numBins_x*numBins_y;
   double binedges_1d[numBins+1];
	binedges_1d[0]=0;
   int phi_bin = 0;
   for (int i=0; i<(numBins); i++)   {
      binedges_1d[i+1] = binedges_x[i%numBins_x+1]+phi_bin*binedges_x[numBins_x];
      if (i%numBins_x==numBins_x-1) phi_bin++;
   }
   
   
   TH1F* tempHist= new TH1F("", "", numBins_x*numBins_y, binedges_1d);
   int binNew = 1;
   for (int j=1; j<=numBins_y; j++){
      for (int i=1; i<=numBins_x; i++){
         tempHist->SetBinContent(binNew, hist2D->GetBinContent(i,j));
         tempHist->SetBinError(binNew, hist2D->GetBinError(i,j));
         binNew++;
      }
   }
   
	return *tempHist;
}

//Import 1D hists for Nominal and systematics
void importHists1D(std::vector<systHists*> &systHists_vec, std::vector<TString> const &samplesToPlot, std::vector<TString> const &mcSamples,
                     std::vector<distr> const &vecDistr) {
   for (systHists* &current : systHists_vec){
      std::cout<<"Importing "<<current->systematicName_<<std::endl;
      std::vector<TString> inputSamples = current->datasets_;
      
      for (TString sSample : inputSamples){
         // ~auto start = high_resolution_clock::now();
         current->hists_.setCurrentSample(sSample);
         TString sSample_alt = sSample;
         if (current->altSampleType) sSample_alt = sSample+"_"+current->systematicName_;    // get correct Sample name if alternative Sample Type
         for (auto const &distr_:vecDistr){
            TString loc;
            loc=distr_.path+distr_.name;
            TH1F* tempHist=current->histReader_->read<TH1F>(loc+"/"+sSample_alt);
            // ~if (tempHist->GetNbinsX()>25) tempHist->Rebin(2);
            // ~if (distr_.name=="dphi_metNearLep" or distr_.name=="dphi_metNearLep_puppi") tempHist->Rebin(4);
            TH1F rebinHist = hist::rebinned(*tempHist,distr_.xMin,distr_.xMax,distr_.nBins);
            if (!current->hasRootFile_) rebinHist.Scale(current->sf_);
            current->hists_.addFilledHist(loc,sSample,rebinHist);
         }
         // ~auto stop = high_resolution_clock::now();
         // ~auto duration = duration_cast<milliseconds>(stop - start);
         // ~std::cout<<sSample<<"   "<<duration.count()<<std::endl;
      }
   }
}

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

void add_Categories(TString const path, io::RootFileReader const &reader_hist, TH1F &out_hist) {   //Function to add the three different categories
   TH1F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH1F*) reader_hist.read<TH1F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH1F) *(reader_hist.read<TH1F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

// get up and down shift for set of systematics
std::pair<TH1F*,TH1F*> getTotalSyst(TH1F* const &nominal, std::vector<systHists*> &systHists_vec, TString const loc, TString const sample="", bool envelope=false){
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   TH1F* tempSys;
   TH1F tempShift;
   TH1F* nominal_ttbar;
   for (auto &current : systHists_vec){
      if (current->systematic_.type() == Systematic::nominal || current->systematic_.type() == Systematic::met40Cut){  //Store sum of ttBar for syst. with alt. samples
         if (sample == "") nominal_ttbar = current->hists_.getSummedHist(loc,current->datasets_ttBar_);
         else nominal_ttbar = current->hists_.getHistogram(loc,sample);
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

// get graph with asym. errors from three histograms (shift=true if only shift and not shift+nominal is given)
TGraphAsymmErrors getErrorGraph(TH1F* const eDOWN, TH1F* const eUP, TH1F* const nominal, bool const shift){
   TGraphAsymmErrors asymmerrors(nominal);
   if (shift) {
      for (int i=0; i<=eUP->GetNbinsX(); i++){
         asymmerrors.SetPointEYhigh(i,eUP->GetBinContent(i+1));
         asymmerrors.SetPointEYlow(i,eDOWN->GetBinContent(i+1));
      }
   }
   else {
      for (int i=0; i<=eUP->GetNbinsX(); i++){
         asymmerrors.SetPointEYhigh(i,abs(eUP->GetBinContent(i+1)-nominal->GetBinContent(i+1)));
         asymmerrors.SetPointEYlow(i,abs(eDOWN->GetBinContent(i+1)-nominal->GetBinContent(i+1)));
      }
   }
   return asymmerrors;
}

void printUnc(TString name, const float &down, const float &up, const float &nominal){
   TString out = TString::Format(" %s & $%.1f(%.1f)$ & $%.1f(%.1f)$\\\\\n",name.ReplaceAll("_","\\_").Data(),down,down/nominal*100,up,up/nominal*100);
   std::cout<<out;
}

// print total yields per contribution
void printTotalYields(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples_merged){
   std::vector<TString> outputSamples(mcSamples_merged);
   std::reverse(outputSamples.begin(),outputSamples.end());
   outputSamples.push_back("MC");
   outputSamples.push_back("data");
   for (TString cat:{"ee","emu","mumu"}){    //Get the number of events per category
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat);
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      float mcYield = mc_total->GetBinContent(6);
      float mcYield_down = syst.first->GetBinContent(6);
      float mcYield_up = syst.second->GetBinContent(6);
      for (TString sample:outputSamples){
         TH1F* temp_hist=hs->getHistogram("cutflow/"+cat,sample);
         float sampleYield = temp_hist->GetBinContent(6);
         std::cout<<std::fixed<<sample<<"&"<<sampleYield<<"&"<<std::setprecision(1)<<sampleYield/mcYield*100<<"\\\\"<<std::endl;
         if (sample=="MC") std::cout<<mcYield_down<<"   "<<mcYield_up<<std::endl;
         if (sample=="MC") std::cout<<mcYield_down/mcYield*100<<"   "<<mcYield_up/mcYield*100<<std::endl;
      }
   }
}

// print brackdown of syst uncertainties
void printUncBreakDown(hist::Histograms<TH1F>* hs, std::vector<systHists*> &systHists_vec, const std::vector<TString> &mcSamples) {
   std::cout<<std::endl<<"-----------------Uncertainties-------------------"<<std::endl;
   for (TString cat:{"ee","emu","mumu"}){
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      TH1F* mc_total=hs->getHistogram("cutflow/"+cat,"MC");
      for (int i=1; i<(systHists_vec.size()); i++){
         if ((i+1)<systHists_vec.size()){    // check if still room for up and down shift
            if (systHists_vec[i]->systematic_.type_str()==systHists_vec[i+1]->systematic_.type_str()){     // check if up and down shift
               std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i],systHists_vec[i+1]};
               std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat);
               printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(6),syst.second->GetBinContent(6),mc_total->GetBinContent(6));
               i++;
               continue;
            }
         }
         std::vector<systHists*> tempVec = {systHists_vec[0],systHists_vec[i]};   // print unc. without up and down shift
         std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,tempVec,"cutflow/"+cat);
         printUnc(systHists_vec[i]->systematic_.type_str(),syst.first->GetBinContent(6),syst.second->GetBinContent(6),mc_total->GetBinContent(6));
      }
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"cutflow/"+cat);     // print total uncertainty
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
         std::cout<<sample<<"   "<<sampleYield_down/mcYield*100<<"   "<<sampleYield_up/mcYield*100<<std::endl;
      }
   }
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
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      // ~mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","DrellYan_M10to50","WW","WZ","ZZ","ttZ_2L","ttZ_QQ","ttW"};
      dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      ttbarSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"};
      signalSamples = {"TTbar_diLepton"};
      break;
      case(2): //2017
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      signalSamples = {"TTbar_diLepton"};
      break;
      case(1): //2016
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      signalSamples = {"TTbar_diLepton"};
      dataSamples = {};
      break;
   }
   std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
   
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_UP","MTOP_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN",};
   std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","MTOP_UP","MTOP_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","ELECTRON_ID_UP","ELECTRON_ID_DOWN","MUON_ID_UP","MUON_ID_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MERENSCALE_UP","MERENSCALE_DOWN",};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_DOWN","PSFSRSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","PSISRSCALE_UP","PSFSRSCALE_UP"};
   // ~std::vector<TString> systToPlot = {"Nominal","BFRAG_UP","BFRAG_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","BSEMILEP_UP","BSEMILEP_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_UP","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MATCH_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","MERENSCALE_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal"};
   
   // 1D plots
   std::vector<distr> vecDistr;
   for(TString channel:{"ee","mumu","emu"}){
      vecDistr.push_back({"cutflow/",channel,0.5,9.5,9});
   }
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         // ~vecDistr.push_back({selection+channel,"Lep_e_pt",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep_mu_pt",0.,600.,50});
         vecDistr.push_back({selection+channel,"Lep1_pt",0.,480.,40});
         vecDistr.push_back({selection+channel,"Lep2_pt",0.,480.,40});
         // ~vecDistr.push_back({selection+channel,"MET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"PuppiMET",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_pT",0.,500.,50});
         // ~vecDistr.push_back({selection+channel,"DNN_MET_dPhi_nextLep",0,3.2,40});
         // ~vecDistr.push_back({selection+channel,"met1000",0.,1000.,50});
         vecDistr.push_back({selection+channel,"mLL",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"pTsumlep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"sumpTlep",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"pTbJet",0.,600.,30});
         vecDistr.push_back({selection+channel,"Jet1_pt",0.,600.,30});
         vecDistr.push_back({selection+channel,"Jet2_pt",0.,600.,30});
         // ~vecDistr.push_back({selection+channel,"dPhiMETnearJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETleadJet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dPhiMETlead2Jet",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep",0.,3.2,32});
         // ~vecDistr.push_back({selection+channel,"dphi_metNearLep_puppi",0.,3.2,32});
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
         vecDistr.push_back({selection+channel,"nJets",-0.5,10.5,11});
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
         vecDistr.push_back({selection+channel,"Lep1_eta",-2.5,2.5,50});
         vecDistr.push_back({selection+channel,"Lep2_eta",-2.5,2.5,50});
         // ~vecDistr.push_back({selection+channel,"Lep1_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Lep2_E",0.,600.,50});
         // ~vecDistr.push_back({selection+channel,"Jet1_phi",-3.2,3.2,50});
         // ~vecDistr.push_back({selection+channel,"Jet2_phi",-3.2,3.2,50});
         vecDistr.push_back({selection+channel,"Jet1_eta",-2.5,2.5,50});
         vecDistr.push_back({selection+channel,"Jet2_eta",-2.5,2.5,50});
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
      }
   }
   
   // 2D plots
   std::map<TString,std::vector<TString>> msPresel_vVars_2D;
      for(TString selection:{"baseline"}){
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         msPresel_vVars_2D.insert(std::pair<TString,std::vector<TString>>(selection+channel,
         {"2d_MetVSdPhiMetNearLep_Puppi",
         "2d_MetVSdPhiMetNearLep_DNN",
         }));
      }
   }
   
   // Setup systematics
   std::vector<systHists*> systHists_vec;
   for (TString syst : systToPlot){
      systHists* temp = new systHists(syst,TString::Format("multiHists/%s/histograms_merged_%s.root",syst.Data(),cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal" || syst=="met40Cut")? samplesToPlot : mcSamples,(syst=="CR_ENVELOPE_UP" || syst=="CR_ENVELOPE_DOWN" || syst=="MTOP_DOWN" || syst=="MTOP_UP")? signalSamples : ttbarSamples);
      systHists_vec.push_back(temp);
   }
   
   // Import 1D hists
   importHists1D(systHists_vec,samplesToPlot,mcSamples,vecDistr);
   
   // Import 2D hists
   // ~importHists2D(systHists_vec,samplesToPlot,mcSamples,msPresel_vVars_2D);
   
   // Define hist collection for nominal
   hist::Histograms<TH1F>* hs;
   hs = &(systHists_vec[0]->hists_);
   
   /*
   //Import 2D hist (trafo to 1D)
   std::vector<float> binedges_x = {0, 20, 40, 60, 80, 100, 120, 140, 160, 195, 230, 400};
   std::vector<float> binedges_y = {0, 0.35, 0.7, 1.05, 1.4, 2.27, 3.14};
   
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      TH1F sumHist;
      for (TString channel:{"ee","mumu","emu"}){
         TString loc="baseline/"+channel+"/2d_MetVSdPhiMetNearLep_Puppi/";
         TH2F* tempHist=histReader.read<TH2F>(loc+sSample);
         TH1F tempHist1D=HistTrafo_2D(tempHist,binedges_x,binedges_y);
         
         for(int i=1; i<=tempHist1D.GetNbinsX(); i++){
            if(tempHist1D.GetBinContent(i)<0) tempHist1D.SetBinContent(i,0);
         }
         
         hs.addFilledHist(loc,sSample,tempHist1D);
         if(channel=="ee") sumHist=tempHist1D;
         else sumHist.Add(&tempHist1D);
      }
      hs.addFilledHist("baseline/all/2d_MetVSdPhiMetNearLep_Puppi/",sSample,sumHist);
   }
   */
   
   // combine different samples to improve readability
   hs->combineSamples("Diboson",{"WW","WZ","ZZ"});
   hs->combineSamples("DrellYan_comb",{"DrellYan","DrellYan_M10to50"});
   // ~hs->combineSamples("DrellYan_comb",{"DrellYan_NLO","DrellYan_M10to50"});
   hs->combineSamples("ttZ",{"ttZ_2L","ttZ_QQ"});
   hs->combineSamples("ttW/Z",{"ttW","ttZ"});
   hs->combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
   std::vector<TString> mcSamples_merged={};
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
      case(1): //2016
      mcSamples_merged = {"ttW/Z","WJetsToLNu","Diboson","DrellYan_NLO","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      hs->combineSamples("MC",mcSamples_merged);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW"});
      break;
   }
   std::map<const TString,Color_t> colormap = {{"TTbar_diLepton",kRed-6},{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   // start 1D plotting
   TCanvas can;
   can.SetLogy();
   gfx::SplitCan sp_can;
   sp_can.pU_.SetLogy();
   for (auto &distr_:vecDistr){    // loop over all histograms defined above
      TString const &sPresel=distr_.path;
      TString const &sVar=distr_.name;
      sp_can.pU_.cd();
      TString loc;
      loc=sPresel+sVar;
      THStack st_mc=hs->getStack(loc,mcSamples_merged,colormap);
      gfx::LegendEntries le=hs->getLegendEntries();
      
      //systematics
      TH1F* hist_mc = hs->getHistogram(loc,{"MC"});
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(hist_mc,systHists_vec,loc);
      TGraphAsymmErrors systGraph = getErrorGraph(syst.first,syst.second,hist_mc,true);
      
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
      st_mc.SetMaximum(1e3*st_mc.GetMaximum());
      st_mc.Draw();     // draw stack
      if(sPresel.Contains("cutflow")) st_mc.GetXaxis()->SetRangeUser(0.5,6.5);
      
      systGraph.SetFillStyle(3001);
      systGraph.Draw("same 2");    // draw syst.
      
      // data plotting part
      auto hist_data = hs->getHistogram(loc,{"data"});
      hist_data->SetLineColor(kBlack);
      hist_data->SetMarkerSize(0.5);
      if(!(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep")|| sVar.Contains("C_em"))) {
         hist_data->Draw("same");
         le.append(*hist_data,"data","lep");
      }
      
      auto hists_BSM=hs->getHistograms(loc,{"TTbar_diLepton"});     //placeholder
      auto BSM_legend=hs->getLegendEntries();    //placeholder
      if (cfg.year_int==1){    //Plot BSM in 2016
         hists_BSM=hs->getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
         for (auto const &h: hists_BSM) {
            h->Draw("same hist");
            // ~h->SetLineWidth(4);
         }
         le+=hs->getLegendEntries();
         BSM_legend=hs->getLegendEntries();
      }
      
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
      
      TGraphAsymmErrors systGraphRatio = getErrorGraph(&ratio_mc_systUp,&ratio_mc_systDown,&ratio_mc,false);
      
      // total unc.
      TH1F* hist_mc_totalDown = (TH1F*)syst.second->Clone();
      TH1F* hist_mc_totalUp = (TH1F*)syst.first->Clone();
      
      for (int i=0; i<=hist_mc_totalUp->GetNbinsX(); i++){
         float stat2 = hist_mc->GetBinError(i)*hist_mc->GetBinError(i);
         float down2 = hist_mc_totalDown->GetBinContent(i)*hist_mc_totalDown->GetBinContent(i);
         float up2 = hist_mc_totalUp->GetBinContent(i)*hist_mc_totalUp->GetBinContent(i);
         hist_mc_totalDown->SetBinContent(i,sqrt(stat2+down2));
         hist_mc_totalUp->SetBinContent(i,sqrt(stat2+up2));
      }
      
      hist_mc_totalDown->Add(hist_mc,-1.);
      hist_mc_totalDown->Scale(-1.);
      hist_mc_totalUp->Add(hist_mc);
               
      TH1F ratio_mc_totalDown = hist::getRatio(*hist_mc_totalDown,*hist_mc,"data/MC",hist::ONLY1);
      TH1F ratio_mc_totalUp = hist::getRatio(*hist_mc_totalUp,*hist_mc,"data/MC",hist::ONLY1);
      
      TGraphAsymmErrors totalUncGraphRatio = getErrorGraph(&ratio_mc_totalUp,&ratio_mc_totalDown,&ratio_mc,false);
      
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
      ratio_mc.SetMaximum(1.25);
      ratio_mc.SetMinimum(0.75);
      totalUncGraphRatio.SetFillColor(kGray);
      systGraphRatio.SetFillColor(kGray+3);
      systGraphRatio.SetFillStyle(3004);
      systGraphRatio.SetLineWidth(0);
      systGraphRatio.SetMarkerSize(0);
      ratio_mc.SetMarkerSize(0);
      ratio_mc.Draw("e2");    // only for axis
      totalUncGraphRatio.Draw("same e2");
      systGraphRatio.Draw("same e2");
      ratio_mc.Draw("axis same");
      // ~ratio.SetMarkerSize(4);
      if(!(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep")|| sVar.Contains("C_em"))) ratio.Draw("pe1 same");
      
      gfx::LegendEntries le_low;
      le_low.append(totalUncGraphRatio,"#sigma_{tot.}","f");
      le_low.append(systGraphRatio,"#sigma_{syst.}","f");
      TLegend leg_low=le_low.buildLegend(.2,.8,0.5,0.95,2);
      leg_low.Draw();
      
      saver.save(sp_can,loc,false,true);
      // ~saver.save(can,loc);
      
      
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
         
         TGraphAsymmErrors systGraph = getErrorGraph(syst.first,syst.second,hist_mc,true);
         
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
