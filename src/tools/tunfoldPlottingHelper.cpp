#include "tunfoldPlottingHelper.hpp"
#include "Config.hpp"
#include <iostream>

using namespace Systematic;

Config const &cfgTemp=Config::get();

std::pair<float,int> tunfoldplotting::getChi2NDF(TH1F* hist_res, TH1F* hist_true) {
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

std::pair<float,int> tunfoldplotting::getChi2NDF_withCorr(TH1F* hist_res, TH1F* hist_true, TH2F* corr_res) {
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

TH2F* tunfoldplotting::get_response(TH2F* responseHist,bool columnNormalized) {
   TH2F* tempHist;
   float sum=0;
   if (columnNormalized){ //Normalize each individual column of diagram
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
   return tempHist;
}

void tunfoldplotting::plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver, const bool is2D) {
   
   for (TString norm :{"","column"}){
      
      TH2F* tempHist = get_response(responseHist, norm == "column");
            
      TCanvas can;
      can.cd();
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.15);
      
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
      
      if (is2D){
         TLine line;
         line.SetLineColor(kGreen);
         line.SetLineWidth(2);
         line.SetLineStyle(2);
         line.DrawLine(7.5,0.5,7.5,21.5);
         line.DrawLine(14.5,0.5,14.5,21.5);
         
         line.DrawLine(0.5,7.5,21.5,7.5);
         line.DrawLine(0.5,14.5,21.5,14.5);
      }
      
      can.RedrawAxis();
      saver->save(can,"response"+norm+"/"+name,true,true,true);
   }
}

void tunfoldplotting::plot_correlation(TH2F* corrMatrix, TString name, io::RootFileSaver* saver){
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

void tunfoldplotting::plot_systBreakdown(std::map<TString,TH1F> const &indShifts, io::RootFileSaver* saver, TString const &name, TString const &method, TString var,
                        std::map<TString,std::vector<TString>> const &systCombinations){
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
   
   if (var == "inclusive") {    //only print unc. for inclusive cross section
      for (auto const& [key, val] : shiftsMap){
         float min = std::min(val.first.GetBinContent(1)*100,val.second.GetBinContent(1)*100);
         float max = std::max(val.first.GetBinContent(1)*100,val.second.GetBinContent(1)*100);
         TString systPrint(key);
         
         // one sided unc.
         if (max==min){
            if (max<0){
               max = 0;
            }
            else{
               min = 0;
            }
         }
         std::cout<<std::fixed<<systPrint.ReplaceAll("_","\\_")<<"&"<<TString::Format("%.2f",abs(min))<<"&"<<TString::Format("%.2f",abs(max))<<"\\\\"<<std::endl;
      }
      
      std::cout<<std::fixed<<"STAT"<<"&"<<TString::Format("%.2f",abs(statShift.first.GetBinContent(1)*100))<<"&"<<TString::Format("%.2f",abs(statShift.second.GetBinContent(1)*100))<<"\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<std::fixed<<"TOTAL"<<"&"<<TString::Format("%.2f",abs(totalShift.first.GetBinContent(1)*100))<<"&"<<TString::Format("%.2f",abs(totalShift.second.GetBinContent(1)*100))<<"\\\\"<<std::endl;
      return;
   }
      
   
   TCanvas can;
   can.cd();
   gPad->SetRightMargin(0.2);
   bool is2D = true;
   bool isNorm = false;
   
   totalShift.first.Scale(-100.);
   totalShift.second.Scale(100.);
   totalShift.first.SetMaximum(1.5*abs(totalShift.first.GetMaximum()));
   totalShift.first.SetMinimum(-1*totalShift.first.GetMaximum());
   totalShift.first.SetStats(false);
   if (!saver->getInternalPath().Contains("2D")){
      totalShift.first.SetTitle(TString::Format(";%s;Uncertainty(%)",var.Data()));
      is2D = false;
   }
   else totalShift.first.SetTitle(";Signal Bin;Uncertainty(%)");
   
   if (saver->getInternalPath().Contains("_norm")){
      isNorm = true;
   }
   
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
   
   totalShift.first.Draw("hist axis");
   
   totalShift.first.Draw("hist same");
   totalShift.second.Draw("hist same");
   statShift.first.Draw("hist same");
   statShift.second.Draw("hist same");
   
   if (is2D){
      // ~totalShift.first.GetYaxis()->SetRangeUser(-100,100);
      totalShift.first.GetYaxis()->SetRangeUser(-60,60);
   }
   else {
      totalShift.first.GetYaxis()->SetRangeUser(-30,30);
   }
   
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
   
   if (is2D){
      var = "2D";
   }
   
   if (isNorm){
      var +=", normalized";
   }
   
   TLatex label=gfx::cornerLabel(var,1);
   label.Draw();
   TLatex label2=gfx::cornerLabel(method,3);
   label2.Draw();
   
   TLegend leg=legE.buildLegend(.80,.15,0.95,.95,1);
   leg.SetTextSize(0.027);
   leg.Draw();
   
   totalShift.first.Draw("axis same");
   saver->save(can,"systBreakdown/"+name+"_"+method,true,true,true);
   
}

void tunfoldplotting::plot_response_diff(std::map<TString,TH2F> const &indResponse, TH2F* const &nomResponse, io::RootFileSaver* saver, TString const &name){
   for (auto const &response : indResponse){
      TH2F* tempHist = (TH2F*)response.second.Clone();
      tempHist->Add(nomResponse,-1);
      tempHist->Divide(nomResponse);
      
      for (int i=1; i<=tempHist->GetNbinsX(); i++){
         for (int j=1; j<=tempHist->GetNbinsX(); j++){
            tempHist->SetBinContent(i,j,abs(tempHist->GetBinContent(i,j)));
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
      tempHist->SetTitle("");
      tempHist->SetMaximum(0.5);
      
      tempHist->SetStats(false);
      tempHist->Draw("colz text");
      
      tempHist->GetYaxis()->SetTitle("reco binNumber");
      tempHist->GetXaxis()->SetTitle("gen binNumber");
      tempHist->GetZaxis()->SetTitleOffset(0.55);
      tempHist->GetZaxis()->SetLabelOffset(0.0015);
      
      saver->save(can,"response_diff/"+name+"/"+response.first,true,true);
   }
}

TH1F tunfoldplotting::getCRenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm){
   io::RootFileReader histReader_CR1(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"CR1",cfgTemp.processFraction*100));
   io::RootFileReader histReader_CR2(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"CR2",cfgTemp.processFraction*100));
   io::RootFileReader histReader_ERDON(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"ERDON",cfgTemp.processFraction*100));
   
   TH1F* hist_CR1 = histReader_CR1.read<TH1F>(path);
   TH1F* hist_CR2 = histReader_CR2.read<TH1F>(path);
   TH1F* hist_ERDON = histReader_ERDON.read<TH1F>(path);
   
   if(norm){
      hist_CR1->Scale(1./hist_CR1->Integral());
      hist_CR2->Scale(1./hist_CR2->Integral());
      hist_ERDON->Scale(1./hist_ERDON->Integral());
   }
   else{
      hist_CR1->Scale(1./cfgTemp.lumi);
      hist_CR2->Scale(1./cfgTemp.lumi);
      hist_ERDON->Scale(1./cfgTemp.lumi);
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(nominal,{hist_CR1,hist_CR2,hist_ERDON});
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getMESCALEenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm){
   io::RootFileReader histReader_MEFACSCALE_UP(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MEFACSCALE_UP",cfgTemp.processFraction*100));
   io::RootFileReader histReader_MEFACSCALE_DOWN(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MEFACSCALE_DOWN",cfgTemp.processFraction*100));
   io::RootFileReader histReader_MERENSCALE_UP(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MERENSCALE_UP",cfgTemp.processFraction*100));
   io::RootFileReader histReader_MERENSCALE_DOWN(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MERENSCALE_DOWN",cfgTemp.processFraction*100));
   io::RootFileReader histReader_MESCALE_UP(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MESCALE_UP",cfgTemp.processFraction*100));
   io::RootFileReader histReader_MESCALE_DOWN(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),"MESCALE_DOWN",cfgTemp.processFraction*100));

   
   TH1F* hist_MEFACSCALE_UP = histReader_MEFACSCALE_UP.read<TH1F>(path);
   TH1F* hist_MEFACSCALE_DOWN = histReader_MEFACSCALE_DOWN.read<TH1F>(path);
   TH1F* hist_MERENSCALE_UP = histReader_MERENSCALE_UP.read<TH1F>(path);
   TH1F* hist_MERENSCALE_DOWN = histReader_MERENSCALE_DOWN.read<TH1F>(path);
   TH1F* hist_MESCALE_UP = histReader_MESCALE_UP.read<TH1F>(path);
   TH1F* hist_MESCALE_DOWN = histReader_MESCALE_DOWN.read<TH1F>(path);
   
   if(norm){
      hist_MEFACSCALE_UP->Scale(1./hist_MEFACSCALE_UP->Integral());
      hist_MEFACSCALE_DOWN->Scale(1./hist_MEFACSCALE_DOWN->Integral());
      hist_MERENSCALE_UP->Scale(1./hist_MERENSCALE_UP->Integral());
      hist_MERENSCALE_DOWN->Scale(1./hist_MERENSCALE_DOWN->Integral());
      hist_MESCALE_UP->Scale(1./hist_MESCALE_UP->Integral());
      hist_MESCALE_DOWN->Scale(1./hist_MESCALE_DOWN->Integral());
   }
   else{
      hist_MEFACSCALE_UP->Scale(1./cfgTemp.lumi);
      hist_MEFACSCALE_DOWN->Scale(1./cfgTemp.lumi);
      hist_MERENSCALE_UP->Scale(1./cfgTemp.lumi);
      hist_MERENSCALE_DOWN->Scale(1./cfgTemp.lumi);
      hist_MESCALE_UP->Scale(1./cfgTemp.lumi);
      hist_MESCALE_DOWN->Scale(1./cfgTemp.lumi);
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(nominal,{hist_MEFACSCALE_UP,hist_MEFACSCALE_DOWN,hist_MERENSCALE_UP,hist_MESCALE_UP,hist_MESCALE_DOWN});
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getPDFenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm){
   
   std::vector<TH1F> histVec;
   for (int i=1; i<=50; i++){    // create reader for each shift
      io::RootFileReader histReader_up(TString::Format("TUnfold/%s/TUnfold_PDF_%i_UP_%.1f.root",filePath.Data(),i,cfgTemp.processFraction*100));
      io::RootFileReader histReader_down(TString::Format("TUnfold/%s/TUnfold_PDF_%i_DOWN_%.1f.root",filePath.Data(),i,cfgTemp.processFraction*100));
      
      TH1F temp_up = *(histReader_up.read<TH1F>(path));
      TH1F temp_down = *(histReader_down.read<TH1F>(path));
      
      if(norm){
         temp_up.Scale(1./temp_up.Integral());
         temp_down.Scale(1./temp_down.Integral());
      }
      else{
         temp_up.Scale(1./cfgTemp.lumi);
         temp_down.Scale(1./cfgTemp.lumi);
      }
      
      histVec.push_back(temp_up);
      histVec.push_back(temp_down);
      
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(nominal,histVec);
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getMTOPunc(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm){
   io::RootFileReader histReader_shift(TString::Format("TUnfold/%s/TUnfold_%s_%.1f.root",filePath.Data(),(up)? "MTOP175p5" : "MTOP169p5",cfgTemp.processFraction*100));
   TH1F* hist_shift = histReader_shift.read<TH1F>(path);
   if(norm) hist_shift->Scale(1./hist_shift->Integral());
   else hist_shift->Scale(1./cfgTemp.lumi);
   
   hist_shift->Add(nominal,-1.);
   hist_shift->Scale(1/3.);
   hist_shift->Add(nominal);
   
   return *hist_shift;
}

TH1F tunfoldplotting::getLUMIunc(TH1F* const &nominal, bool const &up, bool const &norm){
   TH1F hist_shift(*nominal);
   float unc = cfgTemp.systUncFactor.at("LUMI").first;
   float sf = (up)? 1.0-unc : 1.0+unc;    //larger lumi leads to smaller xsec!!
   
   if(norm) hist_shift.Scale(1./hist_shift.Integral());
   else hist_shift.Scale(sf);
   
   return hist_shift;
}

std::pair<TH1F*,TH1F*> tunfoldplotting::getSystUnc(TH1F* const &nominal, TString const &path, TString const &filePath, std::vector<TString> const &systVec, 
                                    bool const &addStat, bool const &withScaleFactor, std::map<TString,TH1F> &indShifts, std::map<TString,TH2F> &indResponse, bool const &norm){
   
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   
   TH1F tempSys;
   TH1F tempShift;
   TH2F tempResponse;
   
   TString responsePath = path;
   TString tauPath = path;
   responsePath.ReplaceAll("_results_","_binning_");
   responsePath.ReplaceAll("hist_unfoldedResult","histMCGenRec");
   tauPath.ReplaceAll("hist_unfoldedResult","tau");
   bool saveResponse = !(path.Contains("/reg") || path.Contains("/BBB"));     // Response has to be read only once
   bool printTau = path.Contains("/reg");     // Tau parameter only available for reg
   
   for (auto &syst : systVec){
            
      if (Systematic::convertType(syst) == Systematic::CR_envelope) {   // derive CR envelope
         tempSys = getCRenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::mtop){      // derive mtop uncertainty
         tempSys = getMTOPunc(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::meScale_envelope){      // derive mtop uncertainty
         tempSys = getMESCALEenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::pdf_envelope){      // derive mtop uncertainty
         tempSys = getPDFenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::lumi){      // derive lumi unc.
         tempSys = getLUMIunc(nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else {   // load "nominal" syst., which need no further handling
         io::RootFileReader histReader(TString::Format(!withScaleFactor ? "TUnfold/%s/TUnfold_%s_%.1f.root" : "TUnfold/%s/TUnfold_SF91_%s_%.1f.root",filePath.Data(),syst.Data(),cfgTemp.processFraction*100));
         tempSys = *(histReader.read<TH1F>(path));
         if (norm) tempSys.Scale(1./tempSys.Integral());
         else tempSys.Scale(1./cfgTemp.lumi);
         
         if (saveResponse) {
            // ~tempResponse = *(get_response(histReader.read<TH2F>(responsePath)));
            tempResponse = *(histReader.read<TH2F>(responsePath));
            indResponse.insert(std::pair<TString,TH2F>(syst,tempResponse));
         }
         // ~else if (printTau) {// print tau parameter
            // ~TParameter<float> *tau_par = histReader.read<TParameter<float>>(tauPath);
            // ~std::cout<<syst<<"   "<<tau_par->GetVal()<<std::endl;
         // ~}
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
