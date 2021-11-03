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

Config const &cfg=Config::get();

// small class for connecting syst., reader and histograms
class systHists
{
   public:
      systHists(TString const systematicName, TString filePath, TString const histPath, std::vector<TString> const datasets):
      systematic_(Systematic::Systematic(systematicName)),
      hists_(hist::Histograms<TH1F>(datasets)),
      datasets_(datasets){
         const std::vector<Systematic::Type> typeVec = Systematic::fileIndependentTypes;
         if (std::find(typeVec.begin(), typeVec.end(), systematic_.type()) == typeVec.end()){
            histReader_ = new io::RootFileReader(filePath,histPath);
         }
         else {
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
      }
      
      Systematic::Systematic systematic_;
      io::RootFileReader* histReader_;
      hist::Histograms<TH1F> hists_;
      std::vector<TString> datasets_;
      float sf_ = 1.0;
      bool hasRootFile_ = true;
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
void importHists1D(std::vector<systHists> &systHists_vec, std::vector<TString> const samplesToPlot, std::vector<TString> const mcSamples,
                     std::map<TString,std::vector<TString>> const msPresel_vVars) {
   for (systHists& current : systHists_vec){
      std::vector<TString> inputSamples = current.datasets_;
      
      for (TString sSample : inputSamples){
         current.hists_.setCurrentSample(sSample);
         for (auto const &sPresel_vVars:msPresel_vVars){
            TString const &sPresel=sPresel_vVars.first;
            for (TString sVar:sPresel_vVars.second){
               TString loc;
               loc=sPresel+sVar;
               TH1F* tempHist=current.histReader_->read<TH1F>(loc+"/"+sSample);
               if (tempHist->GetNbinsX()>25) tempHist->Rebin(2);
               if (sVar=="dphi_metNearLep" or sVar=="dphi_metNearLep_puppi") tempHist->Rebin(4);
               if (!current.hasRootFile_) tempHist->Scale(current.sf_);
               current.hists_.addFilledHist(loc,sSample,*(tempHist));
            }
            if(sPresel.Contains("cutflow")){
               TH1F* tempHist=current.histReader_->read<TH1F>(sPresel+"ee/"+sSample);
               tempHist->Add(current.histReader_->read<TH1F>(sPresel+"emu/"+sSample));
               tempHist->Add(current.histReader_->read<TH1F>(sPresel+"mumu/"+sSample));
               if (!current.hasRootFile_) tempHist->Scale(current.sf_);
               current.hists_.addFilledHist(sPresel+"all",sSample,*(tempHist));
            }
         }
      }
   }
}

//Import 2D hists for Nominal and systematics (onyl signal variable currently)
void importHists2D(std::vector<systHists> &systHists_vec, std::vector<TString> const samplesToPlot, std::vector<TString> const mcSamples,
                     std::map<TString,std::vector<TString>> const msPresel_vVars) {
                        
   std::vector<float> binedges_x = {0, 20, 40, 60, 80, 100, 120, 140, 160, 195, 230, 400};
   std::vector<float> binedges_y = {0, 0.35, 0.7, 1.05, 1.4, 2.27, 3.14};
   for (systHists& current : systHists_vec){
      std::vector<TString> inputSamples;
      if (current.systematic_.type() == Systematic::nominal) inputSamples = samplesToPlot;
      else inputSamples = mcSamples;
      
      for (TString sSample : inputSamples){
         current.hists_.setCurrentSample(sSample);
         for (auto const &sPresel_vVars:msPresel_vVars){
            TString const &sPresel=sPresel_vVars.first;
            for (TString sVar:sPresel_vVars.second){
               TString loc;
               loc=sPresel+sVar;
               TH2F* tempHist=current.histReader_->read<TH2F>(loc+"/"+sSample);
               TH1F tempHist1D=HistTrafo_2D(tempHist,binedges_x,binedges_y);
               current.hists_.addFilledHist(loc,sSample,tempHist1D);
            }
            // ~TH1F* tempHist=current.histReader_->read<TH1F>(sPresel+"ee/"+sSample);
            // ~tempHist->Add(current.histReader_->read<TH1F>(sPresel+"emu/"+sSample));
            // ~tempHist->Add(current.histReader_->read<TH1F>(sPresel+"mumu/"+sSample));
            // ~current.hists_.addFilledHist(sPresel+"all",sSample,*(tempHist));
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
std::pair<TH1F*,TH1F*> getTotalSyst(TH1F* const nominal, std::vector<systHists> &systHists_vec, TString const loc){
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   for (auto &current : systHists_vec){
      if (current.systematic_.type() == Systematic::nominal) continue;
      TH1F* tempSys = current.hists_.getSummedHist(loc);
      TH1F tempShift = phys::getSystShift(*nominal,*tempSys);
      for (int i=0; i<=tempShift.GetNbinsX(); i++){
         float content = tempShift.GetBinContent(i);
         if (content>0) hist_shiftUP->SetBinContent(i,hist_shiftUP->GetBinContent(i)+content*content);
         else hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
      }
   }
   
   hist::sqrtHist(*hist_shiftUP);
   hist::sqrtHist(*hist_shiftDOWN);
   
   return std::make_pair(hist_shiftUP,hist_shiftDOWN);
}

// get graph with asym. errors from three histograms (shift=true if only shift and not shift+nominal is given)
TGraphAsymmErrors getErrorGraph(TH1F* const eUP, TH1F* const eDOWN, TH1F* const nominal, bool const shift){
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
   
extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions");
   
   std::vector<TString> mcSamples={};
   std::vector<TString> dataSamples={};
   switch(cfg.year_int){
      case(3): //2018
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW"};
      dataSamples = {"DoubleMuon","EGamma","MuonEG","SingleMuon"};
      break;
      case(2): //2017
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW"};
      dataSamples = {"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      break;
      case(1): //2016
      mcSamples = {"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW","T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      dataSamples = {};
      break;
   }
   std::vector<TString> samplesToPlot = util::addVectors(mcSamples,dataSamples);
   
   std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","LUMI_UP","LUMI_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN"};
   // ~std::vector<TString> systToPlot = {"Nominal","JESTotal_UP","JESTotal_DOWN"};
   
   // 1D plots
   std::map<TString,std::vector<TString>> msPresel_vVars={
   {"cutflow/",{"ee","emu","mumu"}}
   };
   // ~for(TString selection:{"baseline","baseline_Met200"}){ //Reco 1D Histograms
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         msPresel_vVars.insert(std::pair<TString,std::vector<TString>>(selection+channel,
         {"MET"
         ,"PuppiMET"
         ,"DNN_MET_pT"
         ,"DNN_MET_dPhi_nextLep"
         ,"met1000"
         ,"mLL"
         ,"Lep1_pt"
         ,"Lep2_pt"
         ,"pTsumlep"
         ,"sumpTlep"
         ,"pTbJet"
         ,"Jet1_pt"
         ,"Jet2_pt"
         ,"dPhiMETnearJet"
         ,"dPhiMETleadJet"
         ,"dPhiMETlead2Jet"
         ,"dphi_metNearLep"
         ,"dphi_metNearLep_puppi"
         ,"COSdphi_metNearLep"
         ,"SINdphi_metNearLep"
         ,"dPhiMETbJet"
         ,"dPhiLep1bJet"
         ,"dR_bJetLep1"
         ,"dphi_bJetLep2"
         ,"dphi_bJetnearLep"
         ,"dphi_b1b2"
         ,"dR_b1b2"
         ,"dphi_metLep1"
         ,"dphi_metLep2"
         ,"dphi_metLepsum"
         ,"dPhiLep1Lep2"
         ,"dR_Lep1Lep2"
         ,"nJets"
         ,"nBjets"
         ,"MT2"
         ,"MT"
         ,"mt_MetLep2"
         ,"mt_MetNextLep"
         ,"conMt_Lep1Lep2"
         ,"ST"
         ,"HT"
         ,"sum_STHT"
         ,"sum_mlb"
         ,"METunc_Puppi"
         ,"n_Interactions"
         ,"Lep1_flavor"
         ,"Lep2_flavor"
         ,"Lep1_phi"
         ,"Lep2_phi"
         ,"Lep1_eta"
         ,"Lep2_eta"
         ,"Lep1_E"
         ,"Lep2_E"
         ,"Jet1_phi"
         ,"Jet2_phi"
         ,"Jet1_eta"
         ,"Jet2_eta"
         ,"Jet1_E"
         ,"Jet2_E"
         ,"dPhiMETfarJet"
         ,"dPhiJet1Jet2"
         ,"METsig"
         ,"MHT"
         ,"looseLeptonVeto"
         ,"dPhiMETnearJet_Puppi"
         ,"dPhiMETfarJet_Puppi"
         ,"dPhiMETleadJet_Puppi"
         ,"dPhiMETlead2Jet_Puppi"
         ,"dPhiMETbJet_Puppi"
         ,"dPhiLep1Jet1"
         ,"PFMET_phi"
         ,"PuppiMET_phi"
         ,"CaloMET"
         ,"CaloMET_phi"
         ,"vecsum_pT_allJet"
         ,"vecsum_pT_l1l2_allJet"
         ,"mass_l1l2_allJet"
         ,"ratio_vecsumpTlep_vecsumpTjet"
         ,"mjj"
         ,"C_em_W_p"
         ,"C_em_W_m"
         }));
         // ~{"Jet1_pt","looseLeptonVeto"}));
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
   std::vector<systHists> systHists_vec;
   for (TString syst : systToPlot){
      systHists temp(syst,TString::Format("multiHists/%s/histograms_merged_%s.root",syst.Data(),cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100),(syst=="Nominal")? samplesToPlot : mcSamples);
      systHists_vec.push_back(temp);
   }
   
   // Import 1D hists
   importHists1D(systHists_vec,samplesToPlot,mcSamples,msPresel_vVars);
   
   // Import 2D hists
   importHists2D(systHists_vec,samplesToPlot,mcSamples,msPresel_vVars_2D);
   
   // Define hist collection for nominal
   hist::Histograms<TH1F>* hs;
   hs = &(systHists_vec[0].hists_);
   
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
   hs->combineSamples("ttW/Z",{"ttW","ttZ"});
   hs->combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
   std::vector<TString> MCsamples={};
   switch(cfg.year_int){
      case(3): //2018
      MCsamples = {"ttW/Z","WJetsToLNu","Diboson","DrellYan","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","EGamma","MuonEG","SingleMuon"});
      hs->combineSamples("MC",MCsamples);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttZ","ttW"});
      break;
      case(2): //2017
      MCsamples = {"ttW/Z","WJetsToLNu","Diboson","DrellYan","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      hs->combineSamples("MC",MCsamples);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttZ","ttW"});
      break;
      case(1): //2016
      MCsamples = {"ttW/Z","WJetsToLNu","Diboson","DrellYan_NLO","SingleTop","tt other","TTbar_diLepton"};
      hs->combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      hs->combineSamples("MC",MCsamples);
      hs->combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW"});
      break;
   }
   std::map<const TString,Color_t> colormap = {{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   // start 1D plotting
   TCanvas can;
   can.SetLogy();
   gfx::SplitCan sp_can;
   sp_can.pU_.SetLogy();
   for (auto &sPresel_vVars:msPresel_vVars){    // loop over all histograms defined above
   TString const &sPresel=sPresel_vVars.first;
      if(sPresel.Contains("cutflow")) sPresel_vVars.second.push_back("all");
      for (TString sVar:sPresel_vVars.second){
         sp_can.pU_.cd();
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs->getStack(loc,MCsamples,colormap);
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
         
      }
   }
   

   std::vector<TString> outputSamples(MCsamples);
   std::reverse(outputSamples.begin(),outputSamples.end());
   outputSamples.push_back("MC");
   outputSamples.push_back("data");
   for (TString cat:{"ee","emu","mumu"}){    //Get the number of events per category
      TH1F* mc_total=hs->getHistogram("baseline/"+cat+"/looseLeptonVeto","MC");
      std::pair<TH1F*,TH1F*> syst = getTotalSyst(mc_total,systHists_vec,"baseline/"+cat+"/looseLeptonVeto");
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      // ~for (TString sample:{"TTbar_diLepton","tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttW/Z","MC","data"}){
      for (TString sample:outputSamples){
         TH1F* temp_hist=hs->getHistogram("baseline/"+cat+"/looseLeptonVeto",sample);
         std::cout<<std::fixed<<sample<<"&"<<temp_hist->Integral()<<"&"<<std::setprecision(1)<<temp_hist->Integral()/mc_total->Integral()*100<<"\\\\"<<std::endl;
         // ~std::cout<<sample<<"   "<<temp_hist->Integral()<<"   "<<temp_hist->Integral()/mc_total->Integral()*100<<std::endl;
         if (sample=="MC") std::cout<<syst.first->Integral()<<"   "<<syst.second->Integral()<<std::endl;
         if (sample=="MC") std::cout<<syst.first->Integral()/temp_hist->Integral()*100<<"   "<<syst.second->Integral()/temp_hist->Integral()*100<<std::endl;
      }
   }
   
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
}
