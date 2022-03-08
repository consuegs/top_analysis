//Script to plot the result of TUnfold_unfolding
#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include "TUnfoldDensity.h"
#include <TRandom3.h>
#include <TProfile.h>
#include <TVectorD.h>

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

using namespace std;
using namespace Systematic;

Config const &cfg=Config::get();

struct distrUnfold {
   TString varName;
   float xMin;
   float xMax;
   TString title;
   TString labelFormat;
   bool is2D;
};

std::pair<float,int> getChi2NDF(TH1F* hist_res, TH1F* hist_true) {
   if (hist_res->GetNbinsX()!=hist_true->GetNbinsX()){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=1; i<=hist_res->GetNbinsX(); i++) {
         float diff = hist_res->GetBinContent(i)-hist_true->GetBinContent(i);
         float err = hist_res->GetBinError(i)*1.;
         chi2+= diff*diff/(err*err);
      }
      std::pair<float,int>result(chi2,hist_res->GetNbinsX());
      return result;
   }
}

std::pair<float,int> getChi2NDF_withCorr(TH1F* hist_res, TH1F* hist_true, TH2F* corr_res) {
   if (hist_res->GetNbinsX()!=hist_true->GetNbinsX()){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      
      TMatrixD diff(hist_res->GetNbinsX(),1);
      for (int i=1; i<=hist_res->GetNbinsX();i++){
         diff[i-1][0]=hist_res->GetBinContent(i)-hist_true->GetBinContent(i);
      }
      
      TMatrixD corr(hist_res->GetNbinsX(),hist_res->GetNbinsX());
      for (int i=1; i<=corr_res->GetNbinsX();i++){
         for (int j=1; j<=corr_res->GetNbinsY();j++){
            corr[i-1][j-1]=corr_res->GetBinContent(i,j);
         }
      }
      corr.Invert();
      
      TMatrixD resultMatrix(diff, TMatrixD::kTransposeMult,corr*diff);
      
      float chi2 = resultMatrix[0][0];
      std::pair<float,int>result(chi2,hist_res->GetNbinsX());
      return result;
   }
}

void plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver) {
   
   TH2F* tempHist;
   for (TString norm :{"","column"}){
      float sum=0;
      if (norm=="column"){ //Normalize each individual column of diagram
         tempHist=(TH2F*)responseHist->Clone();
         for (int x=1; x<=tempHist->GetNbinsX(); x++){
            sum=tempHist->Integral(x,x,1,tempHist->GetNbinsY());
            if (sum==0) continue;
            for (int y=1; y<=tempHist->GetNbinsY(); y++){
               if (tempHist->GetBinContent(x,y)!=0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
               else tempHist->SetBinContent(x,y,0.000002);
            }
         }
      }
      else { //Normalize each individual line of diagram
         tempHist=(TH2F*)responseHist->Clone();
         for (int y=1; y<=tempHist->GetNbinsY(); y++){
            sum=tempHist->Integral(1,tempHist->GetNbinsX(),y,y);
            if (sum==0) continue;
            for (int x=1; x<=tempHist->GetNbinsY(); x++){
               if (tempHist->GetBinContent(x,y)!=0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
               else tempHist->SetBinContent(x,y,0.000002);
            }
         }
      }
            
      TCanvas can;
      can.cd();
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.11);
      
      tempHist->GetYaxis()->SetTitleOffset(1.3);
      tempHist->GetXaxis()->SetTitleOffset(0.9);
      tempHist->GetZaxis()->SetTitleOffset(1.3);
      tempHist->GetYaxis()->SetTitleSize(0.05);
      tempHist->GetXaxis()->SetTitleSize(0.05);
      tempHist->GetZaxis()->SetTitleSize(0.05);
      tempHist->GetYaxis()->SetLabelSize(0.04);
      tempHist->GetXaxis()->SetLabelSize(0.04);
      tempHist->GetZaxis()->SetLabelSize(0.04);
      
      tempHist->SetMarkerColor(kRed);
      tempHist->SetMarkerSize(1.5);
      tempHist->SetTitle("");
      tempHist->GetZaxis()->SetTitle((norm=="column")?"column normalized distribution":"line normalized distribution");
      tempHist->SetMinimum(0.000001);
      tempHist->SetMaximum(1.);
      
      tempHist->SetStats(false);
      tempHist->Draw("colz text");
      
      tempHist->GetYaxis()->SetTitle("reco binNumber");
      tempHist->GetXaxis()->SetTitle("gen binNumber");
      tempHist->GetZaxis()->SetTitleOffset(0.55);
      tempHist->GetZaxis()->SetLabelOffset(0.0015);
      
      TLine line;
      line.SetLineColor(kGreen);
      line.SetLineWidth(2);
      line.SetLineStyle(2);
      line.DrawLine(6.5,0.5,6.5,18.5);
      line.DrawLine(12.5,0.5,12.5,18.5);
      
      line.DrawLine(0.5,6.5,18.5,6.5);
      line.DrawLine(0.5,12.5,18.5,12.5);
      
      can.RedrawAxis();
      saver->save(can,"response"+norm+"/"+name,true,true);
   }
}

void plot_correlation(TH2F* corrMatrix, TString name, io::RootFileSaver* saver){
   TCanvas can2D;
   can2D.cd();
   
   gPad->SetRightMargin(0.2);
   corrMatrix->SetStats(0);
   corrMatrix->SetTitle("");
   corrMatrix->GetYaxis()->SetTitleOffset(0.6);
   corrMatrix->GetYaxis()->SetTitle("BinNo");
   corrMatrix->GetXaxis()->SetTitle("BinNo");
   corrMatrix->GetZaxis()->SetTitle("Correlation");
   corrMatrix->GetZaxis()->SetTitleOffset(1.2);
   corrMatrix->GetZaxis()->SetLabelOffset(0.01);
   corrMatrix->SetMarkerColor(kRed);
   corrMatrix->Draw("hcolz text");
   saver->save(can2D,name,true,true);
}

void plot_systBreakdown(std::map<TString,TH1F> const &indShifts, io::RootFileSaver* saver, TString const &name, TString const &method, TString var,
                        std::map<TString,std::vector<TString>> const &systCombinations={}){
   TH1F* zeroes = (TH1F*)indShifts.begin()->second.Clone();
   zeroes->Reset();
   
   std::map<TString,std::pair<TH1F,TH1F>> shiftsMap;
   std::pair<TH1F,TH1F> totalShift;
   std::pair<TH1F,TH1F> statShift;
      
   for (auto const &shift : indShifts){      // connect up and down shifts for plotting
      
      if (Systematic::convertType(shift.first,true) == Systematic::undefinedType){
         if(shift.first == "TOTAL_DOWN") totalShift.first = shift.second;
         else if(shift.first == "TOTAL_UP") totalShift.second = shift.second;
         else if(shift.first == "STAT_DOWN") statShift.first = shift.second;
         else if(shift.first == "STAT_UP") statShift.second = shift.second;
         else std::cout<<"Syst not valid:"<<shift.first<<std::endl;
         continue;
      }
      
      const std::vector<Systematic::Type> up_downTypes = Systematic::upDownTypes;
      if (std::find(up_downTypes.begin(), up_downTypes.end(), Systematic::convertType(shift.first)) != up_downTypes.end()){
         if (Systematic::convertVariation(shift.first) == Systematic::down) shiftsMap[Systematic::convertTypeString(shift.first)].first = shift.second;
         else shiftsMap[Systematic::convertTypeString(shift.first)].second = shift.second;
      }
      else{
         shiftsMap[Systematic::convertTypeString(shift.first)].first = shift.second;
         shiftsMap[Systematic::convertTypeString(shift.first)].second = shift.second;
      }
   }
   
   //combine uncertainties to reduce number of lines in plot
   for (auto const & combi : systCombinations){
      if (shiftsMap.find(combi.second[0]) != shiftsMap.end()){
         TH1F down(shiftsMap[combi.second[0]].first);
         TH1F up(shiftsMap[combi.second[0]].second);
         down.Reset();
         up.Reset();
         
         for (TString const &syst : combi.second){
            std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{shiftsMap[syst].first,shiftsMap[syst].second});
            hist::addQuadr(down,*envelopes.first);
            hist::addQuadr(up,*envelopes.second);
            shiftsMap.erase(syst);
         }
         
         down.Scale(-1.);
         shiftsMap[combi.first].first = down;
         shiftsMap[combi.first].second = up;
      }
   }
      
   
   TCanvas can;
   can.cd();
   gPad->SetRightMargin(0.2);
   
   totalShift.first.Scale(-100.);
   totalShift.second.Scale(100.);
   totalShift.first.SetMaximum(1.5*abs(totalShift.first.GetMaximum()));
   totalShift.first.SetMinimum(-1*totalShift.first.GetMaximum());
   totalShift.first.SetStats(false);
   if (!saver->getInternalPath().Contains("2D")){
      totalShift.first.SetTitle(TString::Format(";%s;Uncertainty(%)",var.Data()));
   }
   else totalShift.first.SetTitle(";Signal Bin;Uncertainty(%)");
   
   totalShift.first.SetFillColor(kOrange);
   totalShift.second.SetFillColor(kOrange);
   totalShift.first.SetLineWidth(0.);
   totalShift.second.SetLineWidth(0.);
   totalShift.first.SetFillStyle(1001);
   totalShift.second.SetFillStyle(1001);
   
   statShift.first.Scale(-100.);
   statShift.second.Scale(100.);
   statShift.first.SetFillColor(kGray);
   statShift.second.SetFillColor(kGray);
   statShift.first.SetLineWidth(0.);
   statShift.second.SetLineWidth(0.);
   statShift.first.SetFillStyle(1001);
   statShift.second.SetFillStyle(1001);
   
   totalShift.first.Draw("hist");
   totalShift.second.Draw("hist same");
   statShift.first.Draw("hist same");
   statShift.second.Draw("hist same");
   
   totalShift.first.GetYaxis()->SetRangeUser(-30,30);
   
   gfx::LegendEntries legE;
   legE.append(totalShift.first,"Total","f");
   legE.append(statShift.first,"Stat.","f");
   
   int currentColor = 2;
   int currentLineStyle = 0;
   
   std::vector<int> lineStyles = {1,7,2,9};
   
   for (auto const &shift : shiftsMap){   // loop over shifts to plot unc.      
      std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{shift.second.first,shift.second.second});
      
      // unc. in %
      envelopes.first->Scale(-100.);  // getEnvelope return abs. shifts
      envelopes.second->Scale(100.);
      
      envelopes.first->SetLineColor(currentColor);
      envelopes.second->SetLineColor(currentColor);
      envelopes.first->SetLineWidth(1);
      envelopes.second->SetLineWidth(1);
      envelopes.first->SetLineStyle(lineStyles[currentLineStyle]);
      envelopes.second->SetLineStyle(lineStyles[currentLineStyle]);
      
      envelopes.first->Draw("hist same");
      envelopes.second->Draw("hist same");
      
      legE.append(*envelopes.first,shift.first,"l");
      
      currentColor++;
      if (currentColor > 9) {
         currentColor = 2;
         currentLineStyle += 1;
      }
      
   }
   
   var.ReplaceAll("(GeV)","");
   
   TLatex label=gfx::cornerLabel(var,1);
   label.Draw();
   TLatex label2=gfx::cornerLabel(method,3);
   label2.Draw();
   
   TLegend leg=legE.buildLegend(.80,.15,0.95,.95,1);
   leg.SetTextSize(0.04);
   leg.Draw();
   
   totalShift.first.Draw("axis same");
   saver->save(can,"systBreakdown/"+name+"_"+method,true,true);
   
}

TH1F getCRenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up){
   io::RootFileReader histReader_CR1(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"CR1",cfg.processFraction*100));
   io::RootFileReader histReader_CR2(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"CR2",cfg.processFraction*100));
   io::RootFileReader histReader_ERDON(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"ERDON",cfg.processFraction*100));
   
   TH1F* hist_CR1 = histReader_CR1.read<TH1F>(path);
   TH1F* hist_CR2 = histReader_CR2.read<TH1F>(path);
   TH1F* hist_ERDON = histReader_ERDON.read<TH1F>(path);
   
   hist_CR1->Scale(1./cfg.lumi);
   hist_CR2->Scale(1./cfg.lumi);
   hist_ERDON->Scale(1./cfg.lumi);
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(nominal,{hist_CR1,hist_CR2,hist_ERDON});
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F getPDFenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up){
   
   std::vector<TH1F> histVec;
   for (int i=1; i<=50; i++){    // create reader for each shift
      std::cout<<"Reading PDF variation "+std::to_string(i)<<std::endl;
      io::RootFileReader histReader_up(TString::Format("TUnfold/%s/TUnfold_PDF_%i_UP_%.1f.root",filePath.Data(),i,cfg.processFraction*100));
      io::RootFileReader histReader_down(TString::Format("TUnfold/%s/TUnfold_PDF_%i_DOWN_%.1f.root",filePath.Data(),i,cfg.processFraction*100));
      
      TH1F temp_up = *(histReader_up.read<TH1F>(path));
      TH1F temp_down = *(histReader_down.read<TH1F>(path));
      
      temp_up.Scale(1./cfg.lumi);
      temp_down.Scale(1./cfg.lumi);
      
      histVec.push_back(temp_up);
      histVec.push_back(temp_down);
      
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(nominal,histVec);
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F getMTOPunc(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up){
   io::RootFileReader histReader_shift(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),(up)? "MTOP175p5" : "MTOP169p5",cfg.processFraction*100));
   TH1F* hist_shift = histReader_shift.read<TH1F>(path);
   hist_shift->Scale(1./cfg.lumi);
   
   hist_shift->Add(nominal,-1.);
   hist_shift->Scale(1/3.);
   hist_shift->Add(nominal);
   
   return *hist_shift;
}

TH1F getLUMIunc(TH1F* const &nominal, bool const &up){
   TH1F hist_shift(*nominal);
   float unc = cfg.systUncFactor.at("LUMI").first;
   float sf = (up)? 1.0-unc : 1.0+unc;    //larger lumi leads to smaller xsec!!
   
   hist_shift.Scale(sf);
   
   return hist_shift;
}

std::pair<TH1F*,TH1F*> getSystUnc(TH1F* const &nominal, TString const &path, TString const &filePath, std::vector<TString> const &systVec, 
                                    bool const &addStat, bool const &withScaleFactor, std::map<TString,TH1F> &indShifts){
   
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   
   TH1F tempSys;
   TH1F tempShift;
   
   for (auto &syst : systVec){
            
      if (Systematic::convertType(syst) == Systematic::CR_envelope) {   // derive CR envelope
         tempSys = getCRenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up);
      }
      else if (Systematic::convertType(syst) == Systematic::mtop){      // derive mtop uncertainty
         tempSys = getMTOPunc(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up);
      }
      else if (Systematic::convertType(syst) == Systematic::pdf_envelope){      // derive mtop uncertainty
         tempSys = getPDFenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up);
      }
      else if (Systematic::convertType(syst) == Systematic::lumi){      // derive lumi unc.
         tempSys = getLUMIunc(nominal,Systematic::convertVariation(syst) == Systematic::up);
      }
      else {   // load "nominal" syst., which need no further handling
         io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",filePath.Data(),syst.Data(),cfg.processFraction*100));
         tempSys = *(histReader.read<TH1F>(path));
         tempSys.Scale(1./cfg.lumi);
      }
      
      // derive shift
      tempShift = phys::getSystShift(*nominal,tempSys);
      
       
      for (int i=0; i<=tempShift.GetNbinsX(); i++){
         float content = tempShift.GetBinContent(i);
         if (content>0) hist_shiftUP->SetBinContent(i,hist_shiftUP->GetBinContent(i)+content*content);
         else hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
      }
      
      // store individual shift
      tempShift.Divide(nominal);
      indShifts.insert(std::pair<TString,TH1F>(syst,tempShift));
      
   }
   
   // store stat shift for unc. plotting
   TH1F statShift(*nominal);
   statShift.Reset();
   
   if (addStat){  //add stat. unc.
      for (int i=0; i<=nominal->GetNbinsX(); i++){
         float content = nominal->GetBinError(i);
         hist_shiftUP->SetBinContent(i,hist_shiftUP->GetBinContent(i)+content*content);
         hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
         
         statShift.SetBinContent(i,content/nominal->GetBinContent(i));
      }
   }
   
   indShifts.insert(std::pair<TString,TH1F>("STAT_UP",statShift));
   indShifts.insert(std::pair<TString,TH1F>("STAT_DOWN",statShift));
   
   hist::sqrtHist(*hist_shiftUP);
   hist::sqrtHist(*hist_shiftDOWN);
   
   // store total shift for unc. plotting
   TH1F hist_shiftUP_rel(*hist_shiftUP);
   TH1F hist_shiftDOWN_rel(*hist_shiftDOWN);
   hist_shiftUP_rel.Divide(nominal);
   hist_shiftDOWN_rel.Divide(nominal);
   indShifts.insert(std::pair<TString,TH1F>("TOTAL_UP",hist_shiftUP_rel));
   indShifts.insert(std::pair<TString,TH1F>("TOTAL_DOWN",hist_shiftDOWN_rel));
   
   return std::make_pair(hist_shiftDOWN,hist_shiftUP);
}

extern "C"
void run()
{
   TString sample = cfg.tunfold_InputSamples[0];
   TString sample_response = cfg.tunfold_ResponseSample;
   
   // Use pT reweighted
   bool withPTreweight = cfg.tunfold_withPTreweight;
   TString scale = cfg.tunfold_scalePTreweight;
   
   // Use DNN instead of pfMET
   bool withDNN = cfg.tunfold_withDNN;
   
   // Use pf instead of PuppiMET
   bool withPF = cfg.tunfold_withPF;
   
   // Use puppi instead of pfMET
   bool withPuppi = !withDNN && !withPF;
   
   // Use same bin numbers for gen/true
   bool withSameBins = cfg.tunfold_withSameBins;
   
   // include signal to pseudo data
   bool withBSM = cfg.tunfold_withBSM;
   
   //Use scale factor
   bool withScaleFactor = cfg.tunfold_withScaleFactor;
   
   //Plot comparison
   bool plotComparison = cfg.tunfold_plotComparison;
   
   //Plot toy studies
   bool plotToyStudies = cfg.tunfold_plotToyStudies;
   
   //Define distributions
   std::vector<distrUnfold> vecDistr;
   vecDistr.push_back({"2D_dPhi_pTnunu",0,400.,";p_{T}^{#nu#nu} (GeV);#sigma (pb)","%.0f",true});
   vecDistr.push_back({"2D_dPhi_pTnunu_new",0,400.,";p_{T}^{#nu#nu} (GeV);#sigma (pb)","%.0f",true});
   vecDistr.push_back({"pTnunu",0,600.,";p_{T}^{#nu#nu} (GeV);#sigma (pb)","%.0f",false});
   vecDistr.push_back({"dPhi",0,3.2,";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);#sigma (pb)","%.1f",false});
   
   //Define syst. unc. to plot
   std::cout<<"!!!!!!!!!!!!!!!!!!!!Uncl. Energy is removed!!!!!!!!!!!!!!!!!!!!!!!!"<<std::endl;
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN","LUMI_UP","LUMI_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN"};
   std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN","BTAGBC_UP","BTAGBC_DOWN","BTAGL_UP","BTAGL_DOWN","ELECTRON_ID_UP","ELECTRON_ID_DOWN","ELECTRON_RECO_UP","ELECTRON_RECO_DOWN","ELECTRON_SCALESMEARING_UP","ELECTRON_SCALESMEARING_DOWN","MUON_ID_UP","MUON_ID_DOWN","MUON_ISO_UP","MUON_ISO_DOWN","MUON_SCALE_UP","MUON_SCALE_DOWN","PU_UP","PU_DOWN","UETUNE_UP","UETUNE_DOWN","MATCH_UP","MATCH_DOWN","TRIG_UP","TRIG_DOWN","MERENSCALE_UP","MERENSCALE_DOWN","MEFACSCALE_UP","MEFACSCALE_DOWN","PSISRSCALE_UP","PSISRSCALE_DOWN","PSFSRSCALE_UP","PSFSRSCALE_DOWN","BFRAG_UP","BFRAG_DOWN","BSEMILEP_UP","BSEMILEP_DOWN","PDF_ALPHAS_UP","PDF_ALPHAS_DOWN","TOP_PT","XSEC_TTOTHER_UP","XSEC_TTOTHER_DOWN","XSEC_DY_UP","XSEC_DY_DOWN","XSEC_ST_UP","XSEC_ST_DOWN","XSEC_OTHER_UP","XSEC_OTHER_DOWN","LUMI_UP","LUMI_DOWN","CR_ENVELOPE_UP","CR_ENVELOPE_DOWN","MTOP_UP","MTOP_DOWN","PDF_ENVELOPE_UP","PDF_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","JER_UP","JER_DOWN"};
   // ~std::vector<TString> systVec = {"JESTotal_UP","JESTotal_DOWN","UNCLUSTERED_UP","UNCLUSTERED_DOWN"};
   // ~std::vector<TString> systVec = {"CR_ENVELOPE_UP","CR_ENVELOPE_DOWN"};
   // ~std::vector<TString> systVec = {"MTOP_UP","MTOP_DOWN"};
   // ~std::vector<TString> systVec = {"LUMI_UP","LUMI_DOWN"};
   // ~std::vector<TString> systVec = {};
   
   //Define combination of syst. unc. for plotting
   // ~std::map<TString,std::vector<TString>> systCombinations = {
      // ~{"JES/JER",{"JESTotal","JER"}},
      // ~{"BTAG",{"BTAGBC","BTAGL"}},
      // ~{"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID","MUON_ISO","MUON_SCALE"}},
      // ~{"PS",{"PSISRSCALE","PSFSRSCALE"}},
      // ~{"XSEC BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
      // ~{"MESCALE",{"MERENSCALE","MEFACSCALE"}},
   // ~};
   std::map<TString,std::vector<TString>> systCombinations = {
      {"JES/JER",{"JESTotal","JER"}},
      {"BTAG",{"BTAGBC","BTAGL"}},
      {"LEPTON",{"ELECTRON_ID","ELECTRON_RECO","ELECTRON_SCALESMEARING","MUON_ID","MUON_ISO","MUON_SCALE"}},
      {"BKG",{"XSEC_TTOTHER","XSEC_DY","XSEC_ST","XSEC_OTHER"}},
      {"THEORY",{"PSISRSCALE","PSFSRSCALE","MERENSCALE","MEFACSCALE","UETUNE","MATCH","BFRAG","BSEMILEP","PDF_ALPHAS","TOP_PT","CR_ENVELOPE","MTOP","PDF_ENVELOPE",}},
   };
   
   for (distrUnfold &dist : vecDistr){
   
      //==============================================
      // step 1 : open output file
      io::RootFileSaver saver(TString::Format("TUnfold/plots%.1f.root",cfg.processFraction*100),TString::Format(!withScaleFactor ? "TUnfold_plotting%.1f/%s" : "TUnfold_plotting_SF_%.1f/%s",cfg.processFraction*100,dist.varName.Data()));

      //==============================================
      // step 2 : read binning schemes and input histograms
      TString input_loc="TUnfold_binning_"+sample+"_"+sample_response;
      TString input_loc_result="TUnfold_results_"+sample+"_"+sample_response;
      if (withBSM) {
         input_loc+="_BSM";
         input_loc_result+="_BSM";
      }
      if (withPF) {
         input_loc+="_PF";
         input_loc_result+="_PF";
      }
      if (withPuppi) {
         input_loc+="_Puppi";
         input_loc_result+="_Puppi";
      }
      if (withDNN) {
         input_loc+="_DNN";
         input_loc_result+="_DNN";
      }
      if (withSameBins) {
         input_loc+="_SameBins";
         input_loc_result+="_SameBins";
      }
      if (withPTreweight) {
         input_loc+="_PTreweight"+scale;
         input_loc_result+="_PTreweight"+scale;
      }
      io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold/%s/TUnfold_Nominal_%.1f.root" : "TUnfold/%s/TUnfold_SF91_Nominal_%.1f.root",input_loc.Data(),cfg.processFraction*100));
      
      TString input_loc_old = input_loc;
      input_loc += "/"+dist.varName;
      input_loc_result += "/"+dist.varName;

      TUnfoldBinning *detectorBinning=histReader.read<TUnfoldBinning>(input_loc+"/detector_binning");
      TUnfoldBinning *generatorBinning=histReader.read<TUnfoldBinning>(input_loc+"/generator_binning");

      if((!detectorBinning)||(!generatorBinning)) {
         cout<<"problem to read binning schemes\n";
      }

      // read histograms
      TH1F *realDis=histReader.read<TH1F>(input_loc+"/histDataTruth");
      TH1F *realDis_response=histReader.read<TH1F>(input_loc+"/histMCGenRec_projX");
      TH1F *unfolded=histReader.read<TH1F>(input_loc_result+"/hist_unfoldedResult");
      TH1F *unfolded_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_unfoldedResult");
      TH1F *unfolded_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_unfoldedResult");
      TH2F *corrMatrix=histReader.read<TH2F>(input_loc_result+"/corr_matrix");
      TH2F *corrMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/corr_matrix");
      TH2F *covMatrix=histReader.read<TH2F>(input_loc_result+"/cov_output");
      TH2F *covMatrix_reg=histReader.read<TH2F>(input_loc_result+"/reg/cov_output");
      TH2F *covMatrix_bbb=histReader.read<TH2F>(input_loc_result+"/BBB/cov_output");
      TH2F *responseMatrix=histReader.read<TH2F>(input_loc+"/histMCGenRec_sameBins");

      if((!realDis)||(!realDis_response)||(!unfolded)) {
         cout<<"problem to read input histograms\n";
      }
      
      // divide by lumi to get xsec
      realDis->Scale(1./cfg.lumi);
      realDis_response->Scale(1./cfg.lumi);
      unfolded->Scale(1./cfg.lumi);
      unfolded_reg->Scale(1./cfg.lumi);
      unfolded_bbb->Scale(1./cfg.lumi);
      
      // Get syst unc
      std::map<TString,TH1F> indShifts;
      std::map<TString,TH1F> indShifts_reg;
      std::map<TString,TH1F> indShifts_bbb;
      std::pair<TH1F*,TH1F*> unfolded_total = getSystUnc(unfolded,input_loc_result+"/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts);
      std::pair<TH1F*,TH1F*> unfolded_reg_total = getSystUnc(unfolded_reg,input_loc_result+"/reg/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_reg);
      std::pair<TH1F*,TH1F*> unfolded_bbb_total = getSystUnc(unfolded_bbb,input_loc_result+"/BBB/hist_unfoldedResult",input_loc_old,systVec,true,withScaleFactor,indShifts_bbb);
      
      //========================
      // Step 3: plotting
      
      gfx::SplitCan can;
      can.can_.SetWindowSize(1800,600);
      can.pU_.SetLogy();
      can.pU_.cd();  //Use upper part of canvas
      
      //Initialize proper binning for plotting
      TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
      
      int num_bins = binning_met.GetNoElements();
      int num_bins_met = binning_met.GetNoElements();
      
      if (dist.is2D){
         TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning(1)));
         num_bins *= (binning_phi.GetNoElements()-1);
      }
      
      binning_met.ResizeTo(binning_met.GetNoElements()+1);
      binning_met[binning_met.GetNoElements()-1] = dist.xMax;  //Plotting end for overflow bin
      
      Double_t xbins[num_bins+1];
      Double_t xbins_minus[num_bins+1];
      Double_t xbins_plus[num_bins+1];
      xbins[0] = 0;   //Plotting start
      xbins_minus[0] = -1.*(binning_met[1]/4.);
      xbins_plus[0] = binning_met[1]/4.;
      
      
      int phi_bin = 0;
      for (int i=0; i<(num_bins); i++)   {
         xbins[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax;
         xbins_minus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-binning_met[1]/4.;
         xbins_plus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+binning_met[1]/4.;
         if (i%num_bins_met==num_bins_met-1) phi_bin++;
      }
      
      unfolded->SetBins(num_bins,xbins);
      unfolded_total.first->SetBins(num_bins,xbins);
      unfolded_total.second->SetBins(num_bins,xbins);
      unfolded_reg->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.first->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.second->SetBins(num_bins,xbins_minus);
      unfolded_bbb->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.first->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.second->SetBins(num_bins,xbins_plus);
      
      if (!dist.is2D){ // set correct binning for plotting syst breakdown (currently only 1D)
         for(std::pair<TString,TH1F> const &shift : indShifts){
            indShifts[shift.first].SetBins(num_bins,xbins);
            indShifts_reg[shift.first].SetBins(num_bins,xbins);
            indShifts_bbb[shift.first].SetBins(num_bins,xbins);
         }
      }
      
      realDis->SetBins(num_bins,xbins);
      realDis_response->SetBins(num_bins,xbins);
      for (int i=1; i<=num_bins; i++) {  //Set proper label for x axis
         int bin_label_no = (i-1)%num_bins_met+1;
         TString label;
         if (bin_label_no == num_bins_met) label = TString::Format(">"+dist.labelFormat,binning_met[bin_label_no-1]);
         else label = TString::Format(dist.labelFormat+"-"+dist.labelFormat,binning_met[bin_label_no-1],binning_met[bin_label_no]);
         unfolded->GetXaxis()->SetBinLabel(i,label);
         realDis->GetXaxis()->SetBinLabel(i,label);
      }
      unfolded->GetXaxis()->SetTickLength(0.);
      unfolded->GetYaxis()->SetTickLength(0.008);
      unfolded->GetXaxis()->SetTitleOffset(1.5);
      unfolded->GetYaxis()->SetTitleOffset(0.8);
      unfolded->GetXaxis()->CenterLabels(false);
         
      unfolded->LabelsOption("v");
      realDis->LabelsOption("v");
      unfolded->SetMaximum(4*unfolded->GetMaximum());
      unfolded->SetMinimum(0.5*unfolded->GetMinimum());
      unfolded->SetLineColor(kBlack);
      unfolded->SetTitle(dist.title);
      unfolded->SetStats(false);
      realDis->SetTitle(dist.title);
      realDis->SetStats(false);
      
      
      unfolded_reg->SetLineColor(kGreen+2);
      unfolded_bbb->SetLineColor(kViolet);
      unfolded_reg->SetMarkerColor(kGreen+2);
      unfolded_bbb->SetMarkerColor(kViolet);
      
      // Setup unc. plotting
      TGraphAsymmErrors unfolded_totalGraph = hist::getErrorGraph(unfolded_total.first,unfolded_total.second,unfolded,true,true);
      TGraphAsymmErrors unfolded_reg_totalGraph = hist::getErrorGraph(unfolded_reg_total.first,unfolded_reg_total.second,unfolded_reg,true,true);
      TGraphAsymmErrors unfolded_bbb_totalGraph = hist::getErrorGraph(unfolded_bbb_total.first,unfolded_bbb_total.second,unfolded_bbb,true,true);
      unfolded_totalGraph.SetLineColor(kBlack);
      unfolded_reg_totalGraph.SetLineColor(kGreen+2);
      unfolded_bbb_totalGraph.SetLineColor(kViolet);
      TGraphErrors unfolded_graph(unfolded);
      TGraphErrors unfolded_reg_graph(unfolded_reg);
      TGraphErrors unfolded_bbb_graph(unfolded_bbb);
      for (int i=0; i<unfolded->GetNbinsX(); i++){    //change x error for plotting
         unfolded_graph.SetPointError(i,binning_met[1]/10.,unfolded->GetBinError(i+1));
         unfolded_reg_graph.SetPointError(i,binning_met[1]/10.,unfolded_reg->GetBinError(i+1));
         unfolded_bbb_graph.SetPointError(i,binning_met[1]/10.,unfolded_bbb->GetBinError(i+1));
      }
      unfolded_graph.SetFillStyle(1001);
      unfolded_reg_graph.SetFillStyle(1001);
      unfolded_bbb_graph.SetFillStyle(1001);
      unfolded_graph.SetLineWidth(0);
      unfolded_reg_graph.SetLineWidth(0);
      unfolded_bbb_graph.SetLineWidth(0);
      unfolded_graph.SetFillColor(kGray);
      unfolded_reg_graph.SetFillColor(kGreen-9);
      unfolded_bbb_graph.SetFillColor(kMagenta-9);
      
      unfolded_totalGraph.SetLineWidth(1.);
      unfolded_reg_totalGraph.SetLineWidth(1.);
      unfolded_bbb_totalGraph.SetLineWidth(1.);
      unfolded_totalGraph.SetMarkerSize(2.);
      unfolded_reg_totalGraph.SetMarkerSize(2.);
      unfolded_bbb_totalGraph.SetMarkerSize(2.);
      
      realDis->SetLineColor(kRed-6);
      realDis->SetFillColor(kRed-6);
      
      if (plotComparison) {
         
         // first draw used for axis and stuff
         unfolded->Draw("px0");
         
         realDis->Draw("hist same");   //Draw into current canvas (axis are not drawn again due to "same")
         
         //draw stat. unc.
         unfolded_graph.Draw("pe2 same");
         unfolded_reg_graph.Draw("pe2 same");
         unfolded_bbb_graph.Draw("pe2 same");
         
         // draw total unc.
         unfolded_totalGraph.Draw("p same");
         unfolded_reg_totalGraph.Draw("p same");
         unfolded_bbb_totalGraph.Draw("p same");
      }
      else unfolded->Draw("pe1");  //Draw into current canvas
      
      realDis_response->SetLineColor(kBlue);
      realDis_response->SetFillColor(kBlue);
      // ~realDis_response->Draw("hist same");
      
      TLatex * atext = new TLatex();
      TLine * aline = new TLine();
      if (dist.is2D){
         //Draw vertical lines and binning ranges for deltaPhi
         atext->SetTextSize(0.03);
         aline->SetLineWidth(2);
         aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
         aline->DrawLine(400,unfolded->GetMinimum(),400,unfolded->GetMaximum());
         aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
         atext->DrawLatex(75,0.5*unfolded->GetMaximum(),"0<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<0.7");
         atext->DrawLatex(475,0.5*unfolded->GetMaximum(),"0.7<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<1.4");
         atext->DrawLatex(875,0.5*unfolded->GetMaximum(),"1.4<|#Delta#phi|(p_{T}^{#nu#nu},nearest l)<3.14");
      }
      
      //Get Chi2 and NDF
      auto Chi2Pair = getChi2NDF(unfolded,realDis);
      auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
      
      //Draw legend
      gfx::LegendEntries legE;
      if (plotComparison) {
         auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
         auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
         auto Chi2Pair_bbb = getChi2NDF(unfolded_bbb,realDis);
         auto Chi2Pair_corr_bbb = getChi2NDF_withCorr(unfolded_bbb,realDis,covMatrix_bbb);
         legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair.first,Chi2Pair.second,Chi2Pair_corr.first,Chi2Pair_corr.second),"pe");
         legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
         legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_bbb.first,Chi2Pair_bbb.second,Chi2Pair_corr_bbb.first,Chi2Pair_corr_bbb.second),"pe");
      }
      else {
         legE.append(*unfolded,"Unfolded","pe");
         atext->DrawLatex(30,10,TString::Format("#chi^{2}/NDF=%.1f/%i",Chi2Pair.first,Chi2Pair.second));
         atext->DrawLatex(30,3,TString::Format("#chi^{2}/NDF(corr.)=%.1f/%i",Chi2Pair_corr.first,Chi2Pair_corr.second));
      }
      legE.append(*realDis,"MC true ttbar","l");
      // ~legE.append(*realDis_response,"MC true ttbar (response)","l");
      TLegend leg=legE.buildLegend(.10,.45,0.25,.65,1);
      leg.SetTextSize(0.03);
      leg.Draw();
      
      unfolded->Draw("axis same");
      
      if (withScaleFactor) {
         atext->DrawLatex(75,10,"With ScaleFactor");
      }
      
      //Change to lower part of canvas
      can.pL_.cd();
      can.pL_.SetBottomMargin(0.45);
      can.pL_.SetTickx(0);
      TH1F ratio;
      TH1F ratio_response;
      TH1F ratio_unfolded;
      TH1F ratio_unfolded_reg;
      TH1F ratio_unfolded_bbb;
      
      // derive ratios
      ratio=hist::getRatio(*realDis,*realDis,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
      ratio_response=hist::getRatio(*realDis_response,*realDis,"ratio",hist::NOERR);
      ratio_unfolded=hist::getRatio(*unfolded,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_reg=hist::getRatio(*unfolded_reg,*realDis,"ratio",hist::ONLY1);
      ratio_unfolded_bbb=hist::getRatio(*unfolded_bbb,*realDis,"ratio",hist::ONLY1);
      
      // derive syst. ratios
      TGraphAsymmErrors ratio_totalGraph = hist::getRatioAsymmGraph(*unfolded_total.first,*unfolded_total.second,*unfolded,*realDis);
      TGraphAsymmErrors ratio_reg_totalGraph = hist::getRatioAsymmGraph(*unfolded_reg_total.first,*unfolded_reg_total.second,*unfolded_reg,*realDis);
      TGraphAsymmErrors ratio_bbb_totalGraph = hist::getRatioAsymmGraph(*unfolded_bbb_total.first,*unfolded_bbb_total.second,*unfolded_bbb,*realDis);
      
      // setup stat. unc. plotting
      TGraphErrors ratio_graph(&ratio_unfolded);
      TGraphErrors ratio_reg_graph(&ratio_unfolded_reg);
      TGraphErrors ratio_bbb_graph(&ratio_unfolded_bbb);
      for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){    //change x error for plotting
         ratio_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded.GetBinError(i+1));
         ratio_reg_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_reg.GetBinError(i+1));
         ratio_bbb_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_bbb.GetBinError(i+1));
      }
      
      ratio_graph.SetFillStyle(1001);
      ratio_reg_graph.SetFillStyle(1001);
      ratio_bbb_graph.SetFillStyle(1001);
      ratio_graph.SetLineWidth(0);
      ratio_reg_graph.SetLineWidth(0);
      ratio_bbb_graph.SetLineWidth(0);
      ratio_graph.SetFillColor(kGray);
      ratio_reg_graph.SetFillColor(kGreen-9);
      ratio_bbb_graph.SetFillColor(kMagenta-9);
      
      // setup axis
      ratio.SetMaximum(1.5);
      ratio.SetMinimum(0.5);
      ratio.SetLineColor(kRed-6);
      ratio.SetMarkerColor(kRed-6);
      ratio.GetYaxis()->SetTitleOffset(0.3);
      ratio.GetXaxis()->SetTitleOffset(1.3);
      ratio.GetXaxis()->SetTickLength(0.);
      if (plotComparison) ratio.Draw("hist");
      else ratio.Draw();
      
      ratio_response.SetLineColor(kBlue);
      ratio_response.SetMarkerColor(kBlue);
      // ~ratio_response.Draw("same");
      
      ratio_totalGraph.SetLineWidth(1.);
      ratio_reg_totalGraph.SetLineWidth(1.);
      ratio_bbb_totalGraph.SetLineWidth(1.);
      ratio_totalGraph.SetMarkerSize(2.);
      ratio_reg_totalGraph.SetMarkerSize(2.);
      ratio_bbb_totalGraph.SetMarkerSize(2.);
      
      if (plotComparison) {
         
         // first draw used for axis and stuff
         ratio_unfolded.Draw("pex0 same");
         
         //draw stat. unc
         ratio_graph.Draw("pe2 same");
         ratio_reg_graph.Draw("pe2 same");
         ratio_bbb_graph.Draw("pe2 same");
         
         //draw tot. unc.
         ratio_totalGraph.Draw("p same");
         ratio_reg_totalGraph.Draw("p same");
         ratio_bbb_totalGraph.Draw("p same");
         
      }
      else {
         ratio_unfolded.SetMarkerSize(0);
         ratio_unfolded.Draw("same");
      }
      
      if (dist.is2D){
         aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
         aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
         aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
      }
      
      //Print rel. uncertainties:
      for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){
         std::cout<<roundf(ratio_totalGraph.GetErrorYlow(i)*100 * 100) / 100<<"   "<<roundf(ratio_totalGraph.GetErrorYhigh(i)*100 * 100) / 100<<std::endl;
      }
      std::cout<<"-----------------------"<<std::endl;
      for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
         std::cout<<roundf(ratio_reg_totalGraph.GetErrorYlow(i)*100 * 100) / 100<<"   "<<roundf(ratio_reg_totalGraph.GetErrorYhigh(i)*100 * 100) / 100<<std::endl;
      }
      std::cout<<"-----------------------"<<std::endl;
      for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
         std::cout<<roundf(ratio_bbb_totalGraph.GetErrorYlow(i)*100 * 100) / 100<<"   "<<roundf(ratio_bbb_totalGraph.GetErrorYhigh(i)*100 * 100) / 100<<std::endl;
      }
      
      //===========================
      // Step 4: save plot
      TString saveName=sample+"_"+sample_response;
      TString saveName2D="correlations/"+sample+"_"+sample_response;
      if (withBSM) {
         saveName+="_BSM";
         saveName2D+="_BSM";
      }
      if (withPF)  {
         saveName+="_PF";
         saveName2D+="_PF";
      }
      if (withPuppi)  {
         saveName+="_Puppi";
         saveName2D+="_Puppi";
      }
      if (withDNN)  {
         saveName+="_DNN";
         saveName2D+="_DNN";
      }
      if (withSameBins) {
         saveName+="_SameBins";
         saveName2D+="_SameBins";
      }
      if (withPTreweight) {
         saveName+="_PTreweight"+scale;
         saveName2D+="_PTreweight"+scale;
      }
      if (plotComparison) saveName="compareMethods/"+saveName;
      saver.save(can,saveName);
      
      //Plot response matrix
      plot_response(responseMatrix,saveName,&saver);
      plot_correlation(corrMatrix,saveName2D,&saver);
      plot_correlation(corrMatrix_reg,saveName2D+"_reg",&saver);
      
      //Plot syst. unc. breakdown
      plot_systBreakdown(indShifts,&saver,saveName,"Nominal",ratio.GetXaxis()->GetTitle(),systCombinations);
      plot_systBreakdown(indShifts_reg,&saver,saveName,"Reg",ratio.GetXaxis()->GetTitle(),systCombinations);
      plot_systBreakdown(indShifts_bbb,&saver,saveName,"BBB",ratio.GetXaxis()->GetTitle(),systCombinations);
      
      //Plot toy studies
      if (plotComparison) saveName=sample+"_"+sample_response;
      if (plotToyStudies) {
         TProfile *profPull=histReader.read<TProfile>(input_loc_result+"/prof_pull");
         TProfile *profPull_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_pull");
         TProfile *profPull_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_pull");
         TProfile *profRes=histReader.read<TProfile>(input_loc_result+"/prof_res");
         TProfile *profRes_reg=histReader.read<TProfile>(input_loc_result+"/reg/prof_res");
         TProfile *profRes_bbb=histReader.read<TProfile>(input_loc_result+"/BBB/prof_res");
         TH1F *histPull=histReader.read<TH1F>(input_loc_result+"/hist_pull");
         TH1F *histPull_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_pull");
         TH1F *histPull_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_pull");
         TH1F *histRes=histReader.read<TH1F>(input_loc_result+"/hist_res");
         TH1F *histRes_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_res");
         TH1F *histRes_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_res");
         TH1F *histChi=histReader.read<TH1F>(input_loc_result+"/hist_chi");
         TH1F *histChi_reg=histReader.read<TH1F>(input_loc_result+"/reg/hist_chi");
         TH1F *histChi_bbb=histReader.read<TH1F>(input_loc_result+"/BBB/hist_chi");
         for (TString type : {"Pull","Residual"}) {   //Plot bin-wise pull and residual
            TCanvas canToy;
            canToy.cd();
            
            //Set temp profiles/hists
            TProfile* temp=0;
            TProfile* temp_reg=0;
            TProfile* temp_bbb=0;
            if (type=="Pull") {
               temp=profPull;
               temp_reg=profPull_reg;
               temp_bbb=profPull_bbb;
            }
            else {
               temp=profRes;
               temp_reg=profRes_reg;
               temp_bbb=profRes_bbb;
            }
            
            //Get RMS from profile
            TH1D tempRMS=*(temp->ProjectionX());
            TH1D tempRMS_reg=*(TH1D*)tempRMS.Clone();
            TH1D tempRMS_bbb=*(TH1D*)tempRMS.Clone();
            for (int i=1; i<=profPull->GetNbinsX(); i++){
               tempRMS.SetBinContent(i,temp->GetBinError(i)*sqrt(temp->GetBinEntries(i)));
               tempRMS_reg.SetBinContent(i,temp_reg->GetBinError(i)*sqrt(temp_reg->GetBinEntries(i)));
               tempRMS_bbb.SetBinContent(i,temp_bbb->GetBinError(i)*sqrt(temp_bbb->GetBinEntries(i)));
            }
            
            temp->SetStats(0);
            temp->GetXaxis()->SetTitle("Bin");
            if (type=="Pull") temp->SetMaximum(1.2);
            else temp->SetMaximum(1.5);
            
            temp->SetMarkerColor(kBlue);
            temp_reg->SetMarkerColor(kRed);
            temp_bbb->SetMarkerColor(kBlack);
            tempRMS.SetMarkerColor(kBlue);
            tempRMS_reg.SetMarkerColor(kRed);
            tempRMS_bbb.SetMarkerColor(kBlack);
            
            temp->SetLineWidth(0);
            temp_reg->SetLineWidth(0);
            temp_bbb->SetLineWidth(0);
            tempRMS.SetLineWidth(0);
            tempRMS_reg.SetLineWidth(0);
            tempRMS_bbb.SetLineWidth(0);
            
            temp->SetMarkerSize(1.3);
            temp_reg->SetMarkerSize(1.3);
            temp_bbb->SetMarkerSize(1.3);
            tempRMS.SetMarkerSize(1.3);
            tempRMS_reg.SetMarkerSize(1.3);
            tempRMS_bbb.SetMarkerSize(1.3);
            
            tempRMS.SetMarkerStyle(4);
            tempRMS_reg.SetMarkerStyle(4);
            tempRMS_bbb.SetMarkerStyle(4);
            
            temp->Draw("");
            temp_reg->Draw("same");
            temp_bbb->Draw("same");
            tempRMS.Draw("same");
            tempRMS_reg.Draw("same");
            tempRMS_bbb.Draw("same");
            
            // ~aline->DrawLine(0.5,0,num_bins+0.5,0);
            // ~if (type=="Pull") aline->DrawLine(0.5,1,num_bins+0.5,1);
            
            gfx::LegendEntries legE_pull;
            legE_pull.append(*temp,"NoReg "+type+" Mean","p");
            legE_pull.append(tempRMS,"NoReg "+type+" RMS","p");
            legE_pull.append(*temp_reg,"Reg "+type+" Mean","p");
            legE_pull.append(tempRMS_reg,"Reg "+type+" RMS","p");
            legE_pull.append(*temp_bbb,"BBB "+type+" Mean","p");
            legE_pull.append(tempRMS_bbb,"BBB "+type+" RMS","p");
            TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.8,.65,2);
            leg_pull.SetTextSize(0.04);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/prof"+type);
         }
         
         for (TString type : {"Pull","Residual","Chi"}) {
            TCanvas canToy;
            canToy.cd();
            
            //Set temp profiles/hists
            TH1* temp=0;
            TH1* temp_reg=0;
            TH1* temp_bbb=0;
            if (type=="Pull") {
               temp=histPull;
               temp_reg=histPull_reg;
               temp_bbb=histPull_bbb;
            }
            else if (type=="Residual") {
               temp=histRes;
               temp_reg=histRes_reg;
               temp_bbb=histRes_bbb;
            }
            else {
               temp=histChi;
               temp_reg=histChi_reg;
               temp_bbb=histChi_bbb;
            }
            
            temp->SetLineColor(kBlue);
            temp_reg->SetLineColor(kRed);
            temp_bbb->SetLineColor(kBlack);
            temp->SetLineWidth(2);
            temp_reg->SetLineWidth(2);
            temp_bbb->SetLineWidth(2);
            
            temp_bbb->SetStats(0);
            temp_bbb->Draw("hist");
            temp->Draw("hist same");
            temp_reg->Draw("hist same");
            gfx::LegendEntries legE_pull;
            legE_pull.append(*temp,TString::Format("NoReg [#mu=%.3f #sigma=%.3f]",temp->GetMean(),temp->GetRMS()),"l");
            legE_pull.append(*temp_reg,TString::Format("Reg [#mu=%.3f #sigma=%.3f]",temp_reg->GetMean(),temp_reg->GetRMS()),"l");
            legE_pull.append(*temp_bbb,TString::Format("BBB [#mu=%.3f #sigma=%.3f]",temp_bbb->GetMean(),temp_bbb->GetRMS()),"l");
            TLegend leg_pull=legE_pull.buildLegend(.10,.45,0.4,.65,1);
            leg_pull.SetTextSize(0.04);
            leg_pull.Draw();
            
            saver.save(canToy,"toyStudies/"+saveName+"/hist"+type);
         }
            
      }
   }
      
}
