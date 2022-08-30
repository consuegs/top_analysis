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
      
      if(is2D){
         can.Size(1200,600);
         gPad->SetRightMargin(0.13);
         gPad->SetLeftMargin(0.08);
         gPad->SetBottomMargin(0.15);
         tempHist->GetYaxis()->SetTitleOffset(0.5);
         tempHist->GetZaxis()->SetTitleOffset(0.55);
         tempHist->GetZaxis()->SetLabelOffset(0.0015);
      }
      else{
         gPad->SetRightMargin(0.2);
         gPad->SetLeftMargin(0.13);
         gPad->SetBottomMargin(0.15);
         tempHist->GetYaxis()->SetTitleOffset(0.9);
         tempHist->GetXaxis()->SetTitleOffset(0.9);
         tempHist->GetZaxis()->SetTitleOffset(0.8);
         tempHist->GetZaxis()->SetLabelOffset(0.0015);
      }
      
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

void tunfoldplotting::plot_pur_stab_eff(TH1F* purity, TH1F* stability, TH1F* efficiency, TString const &name, TString const &xTitle,
                                       const distrUnfold &dist, TUnfoldBinning* const generatorBinning, io::RootFileSaver* saver){
   
   bool is2D = true;
   if (!saver->getInternalPath().Contains("2D")){
      is2D = false;
   }
   
   // Combined 1D Histogram
   TCanvas can;
   can.cd();
   
   purity->GetYaxis()->SetRangeUser(0.,1.);
   
   if (is2D) purity->SetXTitle(xTitle);
   else purity->SetXTitle("Signal Bin");
   
   purity->SetLineColor(kBlue);
   stability->SetLineColor(kRed);
   efficiency->SetLineColor(kGreen);
   
   purity->Draw("hist");
   stability->Draw("hist same");
   efficiency->Draw("hist same");
   
   gfx::LegendEntries legE;
   legE.append(*purity,"Purity","l");
   legE.append(*stability,"Stability","l");
   legE.append(*efficiency,"Efficiency","l");
   
   TLegend leg=legE.buildLegend(.75,.7,0.95,.9,1);
   leg.SetTextSize(0.03);
   leg.Draw();
   
   can.RedrawAxis();
   saver->save(can,"migrations_efficiency/"+name,true,true,true);
   
   // Single 2D Histogram (only for 2D measurement)
   if(is2D){
      TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
      TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning((dist.is2D)? 1 : 0)));
      
      binning_met.ResizeTo(binning_met.GetNoElements()+1);
      binning_met[binning_met.GetNoElements()-1] = dist.xMax;  //Plotting end for overflow bin
      
      Double_t* xbins = &binning_met(0);
      Double_t* ybins = &binning_phi(0);
      TH2F temp2D("",dist.title,binning_met.GetNoElements()-1,xbins,binning_phi.GetNoElements()-1,ybins);
      
      temp2D.SetYTitle("min[#Delta#phi(p_{T}^{#nu#nu},l)]");
      
      for (const TH1F* tempHist : {purity,stability,efficiency}){
         
         // Fill from 1D hist to 2D
         for (int i=1; i<=temp2D.GetNbinsX(); i++){
            for (int j=1; j<=temp2D.GetNbinsY(); j++){
               temp2D.SetBinContent(i,j,tempHist->GetBinContent((j-1)*temp2D.GetNbinsX()+i));
            }
         }
      
         TCanvas can2D;
         can2D.cd();
         
         std::string zTitle = tempHist->GetName();
         zTitle[0] = std::toupper(zTitle[0]);
         
         gPad->SetRightMargin(0.2);
         gPad->SetLeftMargin(0.13);
         temp2D.GetYaxis()->SetTitleOffset(0.8);
         temp2D.GetZaxis()->SetLabelOffset(0.015);
         temp2D.GetZaxis()->SetTitle(zTitle.data());
         temp2D.SetStats(0);
         temp2D.SetMarkerColor(kRed);
         temp2D.SetMarkerSize(1.6);
         
         temp2D.Draw("colz text");
            
         saver->save(can2D,"migrations_efficiency/"+name+"_"+tempHist->GetName()+"_2D",true,true,true);
      }
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
                        std::map<TString,std::vector<TString>> const &systCombinations, const bool jesComparison){
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
            if(!jesComparison) shiftsMap.erase(syst);
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
         TString systPrint(Systematic::getPrintName(key));
         
         // one sided unc.
         if (max==min){
            if (max<0){
               max = 0;
            }
            else{
               min = 0;
            }
         }
         std::cout<<std::fixed<<systPrint<<"&"<<TString::Format("%.2f",abs(min))<<"&"<<TString::Format("%.2f",abs(max))<<"\\\\"<<std::endl;
      }
      
      std::cout<<std::fixed<<"Statistics"<<"&"<<TString::Format("%.2f",abs(statShift.first.GetBinContent(1)*100))<<"&"<<TString::Format("%.2f",abs(statShift.second.GetBinContent(1)*100))<<"\\\\"<<std::endl;
      std::cout<<"\\hline"<<std::endl;
      std::cout<<std::fixed<<"Total"<<"&"<<TString::Format("%.2f",abs(totalShift.first.GetBinContent(1)*100))<<"&"<<TString::Format("%.2f",abs(totalShift.second.GetBinContent(1)*100))<<"\\\\"<<std::endl;
      return;
   }
      
   
   TCanvas can;
   can.cd();
   can.Size(1200,600);
   gPad->SetRightMargin(0.25);
   gPad->SetLeftMargin(0.1);
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
   totalShift.first.GetYaxis()->SetTitleOffset(0.6);
   
   statShift.first.Scale(-100.);
   statShift.second.Scale(100.);
   statShift.first.SetFillColor(kGray);
   statShift.second.SetFillColor(kGray);
   statShift.first.SetLineWidth(0.);
   statShift.second.SetLineWidth(0.);
   statShift.first.SetFillStyle(1001);
   statShift.second.SetFillStyle(1001);
   
   totalShift.first.Draw("hist axis");
   
   gfx::LegendEntries legE;
   
   if(!jesComparison){
      totalShift.first.Draw("hist same");
      totalShift.second.Draw("hist same");
      statShift.first.Draw("hist same");
      statShift.second.Draw("hist same");
      legE.append(totalShift.first,"Total","f");
      legE.append(statShift.first,"Stat.","f");
   }
   
   if (is2D){
      totalShift.first.GetYaxis()->SetRangeUser(-55,55);
   }
   else if(saver->getInternalPath().Contains("dPhi")){
      totalShift.first.GetYaxis()->SetRangeUser(-10,10);
   }
   else {
      totalShift.first.GetYaxis()->SetRangeUser(-30,30);
   }
   
   int currentColor = 2;
   int currentLineStyle = 0;
   
   std::vector<int> lineStyles = {1,7,2,9};
   
   for (auto const &shift : shiftsMap){   // loop over shifts to plot unc.
      
      if (shift.first!="JES Regrouped" && shift.first!="JES Split" && jesComparison) continue;
            
      std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{shift.second.first,shift.second.second});
      
      // unc. in %
      envelopes.first->Scale(100.);
      envelopes.second->Scale(100.);
      
      envelopes.first->SetLineColor(currentColor);
      envelopes.second->SetLineColor(currentColor);
      envelopes.first->SetLineWidth(1);
      envelopes.second->SetLineWidth(1);
      envelopes.first->SetLineStyle(lineStyles[currentLineStyle]);
      envelopes.second->SetLineStyle(lineStyles[currentLineStyle]);
      
      envelopes.first->Draw("hist same");
      envelopes.second->Draw("hist same");
      
      legE.append(*envelopes.first,Systematic::getPrintName(shift.first),"l");
      
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
   
   TLegend leg=legE.buildLegend(.75,.12,0.95,.95,1);
   leg.SetTextSize(0.035);
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


std::pair<TH1F*,TH1F*> tunfoldplotting::getTotalShifts(const std::map<TString,TH1F> &map_combinedShifts, const TH1F &nominal, const bool isNorm, const float &scale){
   
   TH1F* hist_shiftUP = (TH1F*)nominal.Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal.Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   
   for (auto &[key, value]: map_combinedShifts){
      if (key.BeginsWith("TOTAL")) continue;
      
      TH1F tempShift = value;
      tempShift.Scale(1./scale);
      if(isNorm) tempShift.Scale(1./value.Integral());
      
      if (key.BeginsWith("STAT_DOWN")){
         for (int i=0; i<=tempShift.GetNbinsX(); i++){
            float content = tempShift.GetBinContent(i);
            hist_shiftDOWN->SetBinContent(i,hist_shiftDOWN->GetBinContent(i)+content*content);
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
         
   hist::sqrtHist(*hist_shiftUP);
   hist::sqrtHist(*hist_shiftDOWN);
   
   return std::make_pair(hist_shiftDOWN,hist_shiftUP);
}


std::vector<double> tunfoldplotting::plot_UnfoldedResult(TUnfoldBinning* generatorBinning, TH1F* unfolded, TH1F* unfolded_reg, TH1F* unfolded_bbb,
                                                         std::pair<TH1F*,TH1F*> &unfolded_total, std::pair<TH1F*,TH1F*> &unfolded_reg_total, std::pair<TH1F*,TH1F*> &unfolded_bbb_total,
                                                         const float &tau_par, TH1F* realDis, TH1F* realDisAlt,const distrUnfold &dist, const bool plotComparison,
                                                         const TString &saveName, io::RootFileSaver* saver, int &num_bins, const bool rewStudy, const bool divideBinWidth, const bool adaptYaxis){
   //========================
   // Step 3: plotting
   
   // ~gfx::SplitCan can;
   gfx::SplitCan can(0.55);
   can.can_.Size(1000,600);
   can.pU_.SetLogy();
   can.pU_.cd();  //Use upper part of canvas
         
   //Initialize proper binning for plotting
   TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
   TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning((dist.is2D)? 1 : 0)));
   
   num_bins = binning_met.GetNoElements();
   int num_bins_met = binning_met.GetNoElements();
   
   if (dist.is2D){
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
   
   std::vector<double> xbins_vec(xbins, xbins + sizeof xbins / sizeof xbins[0]);
   
   unfolded->SetBins(num_bins,xbins);
   unfolded_total.first->SetBins(num_bins,xbins);
   unfolded_total.second->SetBins(num_bins,xbins);
   if (plotComparison){
      unfolded_reg->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.first->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.second->SetBins(num_bins,xbins_minus);
      unfolded_bbb->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.first->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.second->SetBins(num_bins,xbins_plus);
   }
   else{
      unfolded_reg->SetBins(num_bins,xbins);
      unfolded_reg_total.first->SetBins(num_bins,xbins);
      unfolded_reg_total.second->SetBins(num_bins,xbins);
   }
   
   realDis->SetBins(num_bins,xbins);
   realDisAlt->SetBins(num_bins,xbins);
   std::vector<TString> labelVec;
   for (int i=1; i<=num_bins; i++) {  //Set proper label for x axis
      int bin_label_no = (i-1)%num_bins_met+1;
      TString label;
      if (bin_label_no == num_bins_met) label = TString::Format(">"+dist.labelFormat,binning_met[bin_label_no-1]);
      else label = TString::Format(dist.labelFormat+"-"+dist.labelFormat,binning_met[bin_label_no-1],binning_met[bin_label_no]);
      unfolded->GetXaxis()->SetBinLabel(i,label);
      unfolded_reg->GetXaxis()->SetBinLabel(i,label);
      unfolded_bbb->GetXaxis()->SetBinLabel(i,label);
      realDis->GetXaxis()->SetBinLabel(i,label);
      realDisAlt->GetXaxis()->SetBinLabel(i,label);
      labelVec.push_back(label);
   }
   
   // Divide by bin width
   if(divideBinWidth){
      hist::divideByBinWidth(*realDis);
      hist::divideByBinWidth(*realDisAlt);
      hist::divideByBinWidth(*unfolded_reg);
      hist::divideByBinWidth(*unfolded_reg_total.first);
      hist::divideByBinWidth(*unfolded_reg_total.second);
      hist::divideByBinWidth(*unfolded_bbb);
      hist::divideByBinWidth(*unfolded_bbb_total.first);
      hist::divideByBinWidth(*unfolded_bbb_total.second);
      hist::divideByBinWidth(*unfolded);
      hist::divideByBinWidth(*unfolded_total.first);
      hist::divideByBinWidth(*unfolded_total.second);
   }
   
   // Plotting options
   unfolded->GetXaxis()->SetTickLength(0.);
   unfolded->GetYaxis()->SetTickLength(0.008);
   unfolded->GetXaxis()->SetTitleOffset(1.5);
   unfolded->GetYaxis()->SetTitleOffset(1.5);
   unfolded->GetYaxis()->SetTitleSize(0.04);
   unfolded->GetXaxis()->CenterLabels(false);
      
   unfolded->LabelsOption("v");
   realDis->LabelsOption("v");
   if(adaptYaxis){
      unfolded->SetMaximum(2.2*unfolded->GetMaximum());
      unfolded->SetMinimum(0.4*unfolded->GetMinimum());
   }
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
   unfolded_graph.SetMarkerSize(0.4);
   unfolded_reg_graph.SetMarkerSize(0.4);
   unfolded_bbb_graph.SetMarkerSize(0.4);
   
   unfolded_totalGraph.SetLineWidth(1.);
   unfolded_reg_totalGraph.SetLineWidth(1.);
   unfolded_bbb_totalGraph.SetLineWidth(1.);
   unfolded_totalGraph.SetMarkerSize(0.4);
   unfolded_reg_totalGraph.SetMarkerSize(0.4);
   unfolded_bbb_totalGraph.SetMarkerSize(0.4);
   
   realDis->SetLineColor(kRed-6);
   realDis->SetFillColor(kRed-6);
   realDisAlt->SetLineColor(kBlue-6);
   realDisAlt->SetFillColor(kBlue-6);
   
   if (plotComparison) {
      
      // first draw used for axis and stuff
      // ~unfolded->Draw("px0");
      unfolded->Draw("axis");
      
      realDis->Draw("hist same");   //Draw into current canvas (axis are not drawn again due to "same")
      realDisAlt->Draw("hist same");
      
      //draw stat. unc.
      unfolded_graph.Draw("pe2 same");
      unfolded_reg_graph.Draw("pe2 same");
      unfolded_bbb_graph.Draw("pe2 same");
      
      // draw total unc.
      unfolded_totalGraph.Draw("p same");
      unfolded_reg_totalGraph.Draw("p same");
      unfolded_bbb_totalGraph.Draw("p same");
   }
   else {
      unfolded->Draw("axis");  //Draw into current canvas
      
      realDis->Draw("hist same");
      realDisAlt->Draw("hist same");
      
      unfolded_reg_graph.Draw("pe2 same");
      unfolded_reg_totalGraph.Draw("p same");
   }
   
   
   TLatex * atext = new TLatex();
   TLine * aline = new TLine();
   if (dist.is2D){
      //Draw vertical lines and binning ranges for deltaPhi
      atext->SetTextSize(0.03);
      aline->SetLineWidth(2);
      aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
      aline->DrawLine(400,unfolded->GetMinimum(),400,unfolded->GetMaximum());
      aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
      atext->DrawLatex(75,0.5*unfolded->GetMaximum(),TString::Format("0<|#Delta#phi|<%.2f",binning_phi(1)));
      atext->DrawLatex(475,0.5*unfolded->GetMaximum(),TString::Format("%.2f<|#Delta#phi|<%.2f",binning_phi(1),binning_phi(2)));
      atext->DrawLatex(875,0.5*unfolded->GetMaximum(),TString::Format("%.2f<|#Delta#phi|<%.2f",binning_phi(2),binning_phi(3)));
   }
   
   //Get Chi2 and NDF
   auto Chi2Pair = getChi2NDF(unfolded,realDis);
   // ~auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
   
   //Draw legend
   gfx::LegendEntries legE;
   if (plotComparison) {
      auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
      // ~auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
      auto Chi2Pair_bbb = getChi2NDF(unfolded_bbb,realDis);
      // ~auto Chi2Pair_corr_bbb = getChi2NDF_withCorr(unfolded_bbb,realDis,covMatrix_bbb);
      // ~legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair.first,Chi2Pair.second,Chi2Pair_corr.first,Chi2Pair_corr.second),"pe");
      // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
      // ~legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_bbb.first,Chi2Pair_bbb.second,Chi2Pair_corr_bbb.first,Chi2Pair_corr_bbb.second),"pe");
      // ~legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i]",Chi2Pair.first,Chi2Pair.second),"pe");
      // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_reg.first,Chi2Pair_reg.second),"pe");
      // ~legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_bbb.first,Chi2Pair_bbb.second),"pe");
      legE.append(unfolded_totalGraph,"NoReg","pe");
      legE.append(unfolded_reg_totalGraph,"Reg","pe");
      legE.append(unfolded_bbb_totalGraph,"BBB","pe");
   }
   else {
      // ~legE.append(*unfolded,"Unfolded","pe");
      // ~atext->DrawLatex(30,10,TString::Format("#chi^{2}/NDF=%.1f/%i",Chi2Pair.first,Chi2Pair.second));
      // ~atext->DrawLatex(30,3,TString::Format("#chi^{2}/NDF(corr.)=%.1f/%i",Chi2Pair_corr.first,Chi2Pair_corr.second));
      
      auto Chi2Pair_reg = getChi2NDF(unfolded_reg,realDis);
      // ~auto Chi2Pair_corr_reg = getChi2NDF_withCorr(unfolded_reg,realDis,covMatrix_reg);
      // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i(%.1f/%i)]",Chi2Pair_reg.first,Chi2Pair_reg.second,Chi2Pair_corr_reg.first,Chi2Pair_corr_reg.second),"pe");
      if(tau_par>0) legE.append(unfolded_reg_totalGraph,TString::Format("Unfolded (#tau=%.5f)",tau_par),"pe");
      else legE.append(unfolded_reg_totalGraph,"Unfolded","pe");
   }
   
   if(!rewStudy){
      legE.append(*realDis,"Powheg","l");
      legE.append(*realDisAlt,"MadGraph","l");
   }
   else{
      std::cout<<"------------------------------------------------"<<std::endl;
      legE.append(*realDis,"Truth","l");
      legE.append(*realDisAlt,"Response","l");
   }
   TLegend leg=legE.buildLegend(.16,.5,0.35,.67,1);
   leg.SetTextSize(0.03);
   leg.Draw();
   
   unfolded->Draw("axis same");
   
   //Change to lower part of canvas
   can.pL_.cd();
   can.pL_.SetBottomMargin(0.45);
   can.pL_.SetTickx(0);
   TH1F ratio;
   TH1F ratio_alt;
   TH1F ratio_unfolded;
   TH1F ratio_unfolded_reg;
   TH1F ratio_unfolded_bbb;
   
   // derive ratios
   ratio=hist::getRatio(*realDis,*realDis,"ratio",hist::NOERR);   //Get Ratio between unfolded and true hists
   ratio_alt=hist::getRatio(*realDisAlt,*realDis,"ratio",hist::NOERR);
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
   ratio_graph.SetMarkerSize(0.4);
   ratio_reg_graph.SetMarkerSize(0.4);
   ratio_bbb_graph.SetMarkerSize(0.4);
      
   // setup axis
   if(dist.is2D){
      ratio.SetMaximum(1.49);
      ratio.SetMinimum(0.51);
   }
   else {
      ratio.SetMaximum((rewStudy)? 1.49 : 1.25);
      ratio.SetMinimum((rewStudy)? 0.51 : 0.75);
   }
   ratio.SetLineColor(kRed-6);
   ratio.SetMarkerColor(kRed-6);
   ratio.GetYaxis()->SetTitleOffset(0.4);
   ratio.GetXaxis()->SetTitleOffset(1.7);
   ratio.GetXaxis()->SetLabelOffset(0.015);
   ratio.GetXaxis()->SetTickLength(0.);
   // ~if (plotComparison) ratio.Draw("hist");
   // ~else ratio.Draw();
   ratio.Draw("hist");
   
   ratio_alt.SetLineColor(kBlue-6);
   ratio_alt.SetMarkerColor(kBlue-6);
   ratio_alt.Draw("same");
   
   ratio_totalGraph.SetLineWidth(1.);
   ratio_reg_totalGraph.SetLineWidth(1.);
   ratio_bbb_totalGraph.SetLineWidth(1.);
   ratio_totalGraph.SetMarkerSize(0.4);
   ratio_reg_totalGraph.SetMarkerSize(0.4);
   ratio_bbb_totalGraph.SetMarkerSize(0.4);
   
   gfx::LegendEntries legE_ratio;
   
   if (plotComparison) {
      
      // first draw used for axis and stuff
      // ~ratio_unfolded.Draw("pex0 same");
      
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
      ratio_reg_graph.Draw("pe2 same");
      ratio_reg_totalGraph.Draw("p same");
      ratio.Draw("axis same");
      
      legE_ratio.append(ratio_reg_totalGraph,"#sigma_{tot.}","e");
      legE_ratio.append(ratio_reg_graph,"#sigma_{stat.}","f");
   }
   
   TLegend leg_ratio=legE_ratio.buildLegend(.16,.85,0.45,.95,2);
   leg_ratio.Draw();
   
   if (dist.is2D){
      aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
      aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
      aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   }
   
   
   /*
   //Print rel. uncertainties:
   std::cout<<"---------------------Uncertainties for "<<dist.varName<<"-------------------------------"<<std::endl;
   float meanRelError = 0;
   for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){
      meanRelError = roundf((ratio_totalGraph.GetErrorYlow(i)+ratio_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
      std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
   }
   std::cout<<"-----------------------"<<std::endl;
   for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
      meanRelError = roundf((ratio_reg_totalGraph.GetErrorYlow(i)+ratio_reg_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
      std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
   }
   std::cout<<"-----------------------"<<std::endl;
   for (int i=0; i<ratio_unfolded_reg.GetNbinsX(); i++){
      meanRelError = roundf((ratio_bbb_totalGraph.GetErrorYlow(i)+ratio_bbb_totalGraph.GetErrorYhigh(i))*100 * 100) / 200;
      std::cout<<TString::Format(" %i & $%.1f$ \\\\\n",i+1,meanRelError);
   }
   */
   
   //===========================
   // Step 4: save plot
   saver->save(can,saveName,true,true);
   
   return xbins_vec;
}

TH1F tunfoldplotting::getCRenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, bool const &up, bool const &norm){
   TH1F hist_CR1 = vec_systShifts[0]["CR1"];
   TH1F hist_CR2 = vec_systShifts[0]["CR2"];
   TH1F hist_ERDON = vec_systShifts[0]["ERDON"];
   
   TH1F zeroes = hist_CR1;
   zeroes.Reset();
   
   for (int i=1; i<vec_systShifts.size(); i++){
      hist_CR1.Add(&vec_systShifts[i]["CR1"]);
      hist_CR2.Add(&vec_systShifts[i]["CR2"]);
      hist_ERDON.Add(&vec_systShifts[i]["ERDON"]);
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(&zeroes,{&hist_CR1,&hist_CR2,&hist_ERDON});
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getMESCALEenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, bool const &up, bool const &norm){
   TH1F hist_MESCALE_UP = vec_systShifts[0]["MESCALE_UP"];
   TH1F hist_MESCALE_DOWN = vec_systShifts[0]["MESCALE_DOWN"];
   TH1F hist_MERENSCALE_UP = vec_systShifts[0]["MERENSCALE_UP"];
   TH1F hist_MERENSCALE_DOWN = vec_systShifts[0]["MERENSCALE_DOWN"];
   TH1F hist_MEFACSCALE_UP = vec_systShifts[0]["MEFACSCALE_UP"];
   TH1F hist_MEFACSCALE_DOWN = vec_systShifts[0]["MEFACSCALE_DOWN"];
   
   TH1F zeroes = hist_MESCALE_UP;
   zeroes.Reset();
   
   for (int i=1; i<vec_systShifts.size(); i++){
      hist_MESCALE_UP.Add(&vec_systShifts[i]["MESCALE_UP"]);
      hist_MESCALE_DOWN.Add(&vec_systShifts[i]["MESCALE_DOWN"]);
      hist_MERENSCALE_UP.Add(&vec_systShifts[i]["MERENSCALE_UP"]);
      hist_MERENSCALE_DOWN.Add(&vec_systShifts[i]["MERENSCALE_DOWN"]);
      hist_MEFACSCALE_UP.Add(&vec_systShifts[i]["MEFACSCALE_UP"]);
      hist_MEFACSCALE_DOWN.Add(&vec_systShifts[i]["MEFACSCALE_DOWN"]);
   }
   
   std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(&zeroes,{&hist_MESCALE_UP,&hist_MESCALE_DOWN,&hist_MERENSCALE_UP,&hist_MERENSCALE_DOWN,&hist_MEFACSCALE_UP,&hist_MEFACSCALE_DOWN});
   
   TH1F env_down = *envelopes.first;
   TH1F env_up = *envelopes.second;
               
   if(up) return env_up;
   else return env_down;
}
      


std::map<TString,TH1F> tunfoldplotting::getCombinedUnc(std::vector<std::map<TString,TH1F>>& vec_systShifts, const std::vector<TString>& systVec, const TH1F& combinedResult,
                                                       std::vector<TH1F> nominalResults){
   std::map<TString,TH1F>map_combinedShifts;
   
   TH1F* hist_TotalShiftUP = (TH1F*)combinedResult.Clone();
   TH1F* hist_TotalShiftDOWN = (TH1F*)combinedResult.Clone();
   hist_TotalShiftUP->Reset();
   hist_TotalShiftDOWN->Reset();
   
   TH1F* zeroes = (TH1F*)combinedResult.Clone();
   zeroes->Reset();
   
   float lumiCorr[5][3] = {
      {1.0,0.0,0.0},
      {0.0,2.0,0.0},
      {0.0,0.0,1.5},
      {0.6,0.9,2.0},
      {0.0,0.6,0.2}
   };
   
   bool lumiDone = false;
   bool useCRenvelope = false;
   bool useMEenvelope = false;
   
   std::vector<TString> doneVec;    // vector to store syst which are already processed (used for uncorrelated combinations)
   
   std::cout<<"combined:"<<combinedResult.GetBinContent(1)<<std::endl;
   
   for (const TString& syst : systVec){
      if (syst == "JESUserDefinedHEM1516_DOWN") { // HEM only used for 2018
         map_combinedShifts[syst] = vec_systShifts[3][syst];
         continue;
      }
      
      if (syst.BeginsWith("LUMI")) { // Lumi unc correlation more complex then just (un)corr.
         if (!lumiDone){
            nominalResults[1].Add(&nominalResults[0]);   // add preVFP and post VFP
            
            // uncorrelated part
            TH1F tempHist_uncorr = *(TH1F*)nominalResults[1].Clone();
            tempHist_uncorr.Reset();
            tempHist_uncorr.Add(&nominalResults[1],lumiCorr[0][0]*1e-2);
            
            for (int i=1; i<3; i++){
               hist::addQuadr(tempHist_uncorr,nominalResults[i+1],lumiCorr[i][i]*1e-2);
            }
            
            // correlated part
            TH1F tempHist_corr = *(TH1F*)nominalResults[1].Clone();
            tempHist_corr.Reset();
            tempHist_corr.Add(&nominalResults[1],lumiCorr[3][0]*1e-2);
            for (int i=1; i<3; i++){
               tempHist_corr.Add(&nominalResults[i+1],lumiCorr[3][i]*1e-2);
            }
            
            // correlated part (17/18)
            TH1F tempHist_corr1718 = *(TH1F*)nominalResults[2].Clone();
            tempHist_corr1718.Reset();
            tempHist_corr1718.Add(&nominalResults[2],lumiCorr[4][1]*1e-2);
            tempHist_corr1718.Add(&nominalResults[3],lumiCorr[4][2]*1e-2);
            
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
      
      TH1F tempHist;
      
      // Get CR envelope
      if (Systematic::convertType(syst) == Systematic::CR_envelope) {
         tempHist = getCRenvelopeCombined(vec_systShifts,Systematic::convertVariation(syst) == Systematic::up);
         map_combinedShifts[syst] = tempHist;
         useCRenvelope = true;
         continue;
      }
      else if ((std::find(Systematic::crTypes.begin(), Systematic::crTypes.end(), Systematic::convertType(syst)) != Systematic::crTypes.end()) && useCRenvelope) continue;  //ignore shifts already used in envelope
      
      // Get MEscale envelope
      if (Systematic::convertType(syst) == Systematic::meScale_envelope) {
         tempHist = getMESCALEenvelopeCombined(vec_systShifts,Systematic::convertVariation(syst) == Systematic::up);
         map_combinedShifts[syst] = tempHist;
         useMEenvelope = true;
         continue;
      }
      else if ((std::find(Systematic::meTypes.begin(), Systematic::meTypes.end(), Systematic::convertType(syst)) != Systematic::meTypes.end()) && useMEenvelope) continue;  //ignore shifts already used in envelope
      
      // Combine all other systematics
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
         if (std::find(doneVec.begin(), doneVec.end(), syst) != doneVec.end()) continue;     //ignore syst which are already processed
         
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
   
   // Add statistic uncertainty
   map_combinedShifts["STAT_UP"] = combinedResult;
   map_combinedShifts["STAT_DOWN"] = combinedResult;
   for (int i=1; i<=map_combinedShifts["STAT_UP"].GetNbinsX(); i++){
      
      float content = map_combinedShifts["STAT_UP"].GetBinError(i);
      
      hist_TotalShiftUP->SetBinContent(i,hist_TotalShiftUP->GetBinContent(i)+content*content);
      hist_TotalShiftDOWN->SetBinContent(i,hist_TotalShiftDOWN->GetBinContent(i)+content*content);
      
      map_combinedShifts["STAT_UP"].SetBinContent(i,content);
      map_combinedShifts["STAT_DOWN"].SetBinContent(i,content);
   }
   
   // Add total uncertainty
   hist::sqrtHist(*hist_TotalShiftUP);
   hist::sqrtHist(*hist_TotalShiftDOWN);
   map_combinedShifts["TOTAL_UP"] = *hist_TotalShiftUP;
   map_combinedShifts["TOTAL_DOWN"] = *hist_TotalShiftDOWN;
   
   return map_combinedShifts;
}
            
