#include "tunfoldPlottingHelper.hpp"
#include "Config.hpp"
#include <iostream>
#include <fstream>
#include <algorithm>
#include <iterator>

using namespace Systematic;

Config const &cfgTemp=Config::get();

std::pair<float,int> tunfoldplotting::getChi2NDF(TH1F* hist_res, TH1F* hist_true, const bool includeTrueError) {
   if (hist_res->GetNbinsX()!=hist_true->GetNbinsX()){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=1; i<=hist_res->GetNbinsX(); i++) {
         float diff = hist_res->GetBinContent(i)-hist_true->GetBinContent(i);
         float err = hist_res->GetBinError(i)*1.;
         float errTrue = (includeTrueError)? hist_true->GetBinError(i)*1 : 0.;
         chi2+= diff*diff/(err*err+errTrue*errTrue);
      }
      std::pair<float,int>result(chi2,hist_res->GetNbinsX());
      return result;
   }
}

std::pair<float,int> tunfoldplotting::getChi2NDF(TGraphAsymmErrors* graph_res, TH1F* hist_true, const bool includeTrueError) {
   if ((graph_res->GetN()-1)!=hist_true->GetNbinsX()){
      std::cout<<"Histogram and Graph in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=1; i<=hist_true->GetNbinsX(); i++) {
         double tempX;
         double tempY;
         graph_res->GetPoint(i-1,tempX,tempY);
         float diff = tempY-hist_true->GetBinContent(i);
         float err = graph_res->GetErrorY(i-1)*1.;   //symmetrizes error in y for graph
         float errTrue = (includeTrueError)? hist_true->GetBinError(i)*1 : 0.;
         chi2+= diff*diff/(err*err+errTrue*errTrue);
         
      }
      std::pair<float,int>result(chi2,hist_true->GetNbinsX());
      return result;
   }
}

std::pair<float,int> tunfoldplotting::getChi2NDF(TGraph* graph_res, TH1F* hist_true, const bool includeTrueError) {
   if (graph_res->GetN()!=hist_true->GetNbinsX()){
      std::cout<<"Histogram and Graph in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=1; i<=hist_true->GetNbinsX(); i++) {
         double tempX;
         double tempY;
         graph_res->GetPoint(i-1,tempX,tempY);
         float diff = tempY-hist_true->GetBinContent(i);
         float err = graph_res->GetErrorY(i-1)*1.;   //symmetrizes error in y for graph
         float errTrue = (includeTrueError)? hist_true->GetBinError(i)*1 : 0.;
         chi2+= diff*diff/(err*err+errTrue*errTrue);
      }
      std::pair<float,int>result(chi2,hist_true->GetNbinsX());
      return result;
   }
}

std::pair<float,int> tunfoldplotting::getChi2NDF(TGraphAsymmErrors* graph_res, TGraphAsymmErrors* graph_true, const bool includeTrueError) {
   if (graph_res->GetN()!=graph_true->GetN()){
      std::cout<<"Graphs in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      float chi2 = 0;
      for (int i=0; i<(graph_res->GetN()-1); i++) {
         double tempX;
         double tempY;
         graph_res->GetPoint(i,tempX,tempY);
         float diff = tempY;
         graph_true->GetPoint(i,tempX,tempY);
         diff -= tempY;
         float err = graph_res->GetErrorY(i)*1.;   //symmetrizes error in y for graph
         float errTrue = (includeTrueError)? graph_true->GetErrorY(i)*1. : 0.;
         chi2+= diff*diff/(err*err+errTrue*errTrue);
      }
      std::pair<float,int>result(chi2,graph_res->GetN()-1);
      return result;
   }
}

std::pair<float,int> tunfoldplotting::getChi2NDF_withCorr(TH1F* hist_res, TGraphAsymmErrors* graph_true, TH2D* cov_res, bool const norm, const bool includeTrueError) {
   if (hist_res->GetNbinsX()!=(graph_true->GetN()-1)){
      std::cout<<"Histograms in Chi2 function have different number of bins"<<std::endl;
      throw;
   }
   else {
      
      int nBins = hist_res->GetNbinsX();
      
      if (norm) nBins-=1;     //remove last bin in normalized distributions due to loosing one degree of freedom
      
      TMatrixD diff(nBins,1);
      for (int i=1; i<=nBins;i++){
         double tempX;
         double tempY;
         graph_true->GetPoint(i-1,tempX,tempY);
         // ~graph_true->GetPoint(i,tempX,tempY);
         diff[i-1][0]=hist_res->GetBinContent(i)-tempY;
         // ~diff[i-1][0]=hist_res->GetBinContent(i+1)-tempY;
      }
      
      TMatrixD cov(nBins,nBins);
      TMatrixD corr(nBins,nBins);
      for (int i=1; i<=nBins;i++){
         // ~cov[i-1][i-1]=cov_res->GetBinContent(i,i);  //only for testing diagonal elements
         for (int j=1; j<=nBins;j++){
            cov[i-1][j-1]=cov_res->GetBinContent(i,j);
            // ~cov[i-1][j-1]=cov_res->GetBinContent(i+1,j+1);
            if(i==j && includeTrueError) cov[i-1][j-1] += graph_true->GetErrorY(i-1)*graph_true->GetErrorY(i-1);
            corr[i-1][j-1]=cov_res->GetBinContent(i,j)/(sqrt(cov_res->GetBinContent(i,i)*cov_res->GetBinContent(j,j)));
         }
      }
      
      if (!includeTrueError){
         diff.Print();
         // ~cov.Print();
         // ~std::cout<<cov.Determinant()<<std::endl;
         corr.Print();
      }
      
      cov.Invert();
      
      // ~cov.Print();
      
      TMatrixD resultMatrix(diff, TMatrixD::kTransposeMult,cov*diff);
            
      float chi2 = resultMatrix[0][0];
      std::pair<float,int>result(chi2,nBins);
      return result;
   }
}

TH2F* tunfoldplotting::get_response(TH2F* responseHist,bool columnNormalized, bool includeUnderflow) {
   TH2F* tempHist;
   float sum=0;
   int startBin = (includeUnderflow)? 0 : 1;
   if (columnNormalized){ //Normalize each individual column of diagram
      tempHist=(TH2F*)responseHist->Clone();
      for (int x=startBin; x<=tempHist->GetNbinsX(); x++){
         sum=tempHist->Integral(x,x,startBin,tempHist->GetNbinsY());
         if (sum==0) continue;
         for (int y=startBin; y<=tempHist->GetNbinsY(); y++){
            if (tempHist->GetBinContent(x,y)>0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
            else tempHist->SetBinContent(x,y,0.0000002);
         }
      }
   }
   else { //Normalize each individual line of diagram
      tempHist=(TH2F*)responseHist->Clone();
      for (int y=startBin; y<=tempHist->GetNbinsY(); y++){
         sum=tempHist->Integral(startBin,tempHist->GetNbinsX(),y,y);
         if (sum==0) continue;
         for (int x=startBin; x<=tempHist->GetNbinsY(); x++){
            if (tempHist->GetBinContent(x,y)>0)tempHist->SetBinContent(x,y,tempHist->GetBinContent(x,y)/sum);
            else tempHist->SetBinContent(x,y,0.0000002);
         }
      }
   }
   return tempHist;
}

void tunfoldplotting::plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver, const bool is2D) {
   
   for (TString norm :{"","column"}){
      for (TString sUnderFlow : {"","includeUnderflow"}){
      
         TH2F* tempHist = get_response(responseHist, norm == "column", sUnderFlow == "includeUnderflow");
               
         TCanvas can;
         can.cd();
         
         if(is2D){
            can.Size(1200,600);
            gPad->SetRightMargin(0.13);
            gPad->SetLeftMargin(0.09);
            gPad->SetBottomMargin(0.15);
            tempHist->GetYaxis()->SetTitleOffset(0.5);
            tempHist->GetZaxis()->SetTitleOffset(0.55);
            tempHist->GetZaxis()->SetLabelOffset(0.003);
         }
         else{
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.14);
            gPad->SetBottomMargin(0.15);
            tempHist->GetYaxis()->SetTitleOffset(0.9);
            tempHist->GetXaxis()->SetTitleOffset(0.9);
            tempHist->GetZaxis()->SetTitleOffset(0.8);
            tempHist->GetZaxis()->SetLabelOffset(0.003);
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
         tempHist->GetZaxis()->SetTitle((norm=="column")?"Column normalized distribution":"Line-normalized distribution");
         //tempHist->SetMinimum(0.000000000001);
         tempHist->SetMaximum(1.);
         
         if (sUnderFlow=="includeUnderflow") {
            tempHist->GetYaxis()->SetRangeUser(-1, tempHist->GetNbinsY());
            tempHist->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1,"#notin");
            
            // ~std::cout<<(TObjString)*tempHist->GetYaxis()->GetLabels()->At(0)<<std::endl;
         } 
         
         tempHist->SetStats(false);
         tempHist->Draw("colz text");
         
         tempHist->GetYaxis()->SetTitle("Bin(Detector level)");
         tempHist->GetXaxis()->SetTitle("Bin(Particle level)");

         
         if (is2D){
            TLine line;
            line.SetLineColor(kGreen);
            line.SetLineWidth(2);
            line.SetLineStyle(2);
            if (sUnderFlow=="includeUnderflow"){
               line.DrawLine(4.5,-0.5,4.5,12.5);
               line.DrawLine(8.5,-0.5,8.5,12.5);
               line.DrawLine(0.5,0.5,12.5,0.5);
            }
            else{
               line.DrawLine(4.5,0.5,4.5,12.5);
               line.DrawLine(8.5,0.5,8.5,12.5);
            }
            
            line.DrawLine(0.5,4.5,12.5,4.5);
            line.DrawLine(0.5,8.5,12.5,8.5);
         }
         
         can.RedrawAxis();
         saver->save(can,"response"+norm+sUnderFlow+"/"+name,true,true,true);
      }
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
   
   purity->GetYaxis()->SetRangeUser(0.,0.95);
   
   if (!is2D) purity->SetXTitle(xTitle);
   else purity->SetXTitle("Particle level bin");
   
   purity->SetLineColor(kBlue);
   stability->SetLineColor(kRed);
   efficiency->SetLineColor(kGreen);
   purity->SetMarkerColor(kBlue);
   stability->SetMarkerColor(kRed);
   efficiency->SetMarkerColor(kGreen);
   purity->SetMarkerSize(0.7);
   stability->SetMarkerSize(0.7);
   efficiency->SetMarkerSize(0.7);
   
   purity->SetBinError(1,1e-6);    // required to plot bin width
   stability->SetBinError(1,1e-6);
   efficiency->SetBinError(1,1e-6);
   purity->Draw("axis");
   
   // Draw horizontal lines for min. requirements
   TLine * requLine = new TLine();
   float lim_pur_stab = (is2D)? 0.3 : 0.5;
   requLine->SetLineWidth(2);
   gStyle->SetLineStyleString(25,"12 36");
   requLine->SetLineStyle(2);
   requLine->SetLineColor(kGreen);
   requLine->DrawLine(purity->GetXaxis()->GetXmin(),0.2,purity->GetXaxis()->GetXmax(),0.2);
   requLine->SetLineStyle(2);
   requLine->SetLineColor(kBlue);
   requLine->DrawLine(purity->GetXaxis()->GetXmin(),lim_pur_stab,purity->GetXaxis()->GetXmax(),lim_pur_stab);
   requLine->SetLineStyle(25);
   requLine->SetLineColor(kRed);
   requLine->DrawLine(purity->GetXaxis()->GetXmin(),lim_pur_stab,purity->GetXaxis()->GetXmax(),lim_pur_stab);
   
   purity->Draw("pe same");
   stability->Draw("pe same");
   efficiency->Draw("pe same");
   
   gfx::LegendEntries legE;
   legE.append(*purity,"Purity","pl");
   legE.append(*stability,"Stability","pl");
   legE.append(*efficiency,"Efficiency","pl");
   
   // Draw vertical lines for phi bins in 2D
   TVectorD binning_met(*(generatorBinning->FindNode("signal")->GetDistributionBinning(0)));
   TVectorD binning_phi(*(generatorBinning->FindNode("signal")->GetDistributionBinning((dist.is2D)? 1 : 0)));
   if(is2D){
      TLine * aline = new TLine();
      int nBinsx = binning_met.GetNoElements();
      for (int i=1; i<(binning_phi.GetNoElements()-1); i++){
         aline->DrawLine(i*nBinsx+0.5,purity->GetMinimum(),i*nBinsx+0.5,purity->GetMaximum());
      }
   }
   
   can.RedrawAxis();
   
   // ~TLegend leg=legE.buildLegend(.2,.88,.95,.92,3);
   TLegend leg=legE.buildLegend(gPad->GetLeftMargin(),.88,1.-gPad->GetRightMargin(),1.-gPad->GetTopMargin(),3);
   leg.SetTextSize(0.04);
   leg.SetBorderSize(1);
   leg.SetFillStyle(1001);
   leg.Draw();
   
   saver->save(can,"migrations_efficiency/"+name,true,true,true);
   
   // Single 2D Histogram (only for 2D measurement)
   if(is2D){
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
                                          const bool & isData, std::map<TString,std::vector<TString>> const &systCombinations, const bool jesComparison){
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
   can.Size(1000,600);
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
      totalShift.first.GetYaxis()->SetRangeUser(-30,30);
   }
   else if(saver->getInternalPath().Contains("dPhi")){
      totalShift.first.GetYaxis()->SetRangeUser(-10,10);
   }
   else {
      // ~totalShift.first.GetYaxis()->SetRangeUser(-30,30);
      totalShift.first.GetYaxis()->SetRangeUser(-45,45);
   }
   
   int currentColor = 0;
   int currentLineStyle = 0;
   
   std::vector<int> lineStyles = {1,7,2,9};
   // ~std::vector<int> colors = {1,2,3,4,6,7,12,30,49};
   std::vector<int> colors = {1,2,3,4,6,7,45};
   // ~std::vector<int> colors = {2,3,4};
   // ~std::vector<int> colors = {1,2,4};
   
   for (auto const &shift : shiftsMap){   // loop over shifts to plot unc.
      
      if (shift.first!="JES Regrouped" && shift.first!="JES Split" && jesComparison) continue;
            
      std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{shift.second.first,shift.second.second});
      
      // unc. in %
      envelopes.first->Scale(100.);
      envelopes.second->Scale(100.);
      
      envelopes.first->SetLineColor(colors[currentColor]);
      envelopes.second->SetLineColor(colors[currentColor]);
      // ~envelopes.first->SetLineWidth(1.5);
      // ~envelopes.second->SetLineWidth(1.5);
      envelopes.first->SetLineWidth(2);
      envelopes.second->SetLineWidth(2);
      envelopes.first->SetLineStyle(lineStyles[currentLineStyle]);
      envelopes.second->SetLineStyle(lineStyles[currentLineStyle]);
      
      envelopes.first->Draw("hist same");
      envelopes.second->Draw("hist same");
      
      legE.append(*envelopes.first,Systematic::getPrintName(shift.first),"l");
      
      currentColor++;
      if (currentColor == colors.size()) {
         currentColor = 0;
         currentLineStyle += 1;
      }
      
   }
   
   var.ReplaceAll("(GeV)","");
   
   if (is2D){
      var = "2D";
      TLine * aline = new TLine();
      aline->SetLineStyle(7);
      aline->DrawLine(4.5,-30.,4.5,30.);
      aline->DrawLine(8.5,-30.,8.5,30.);
   }
   
   if (isNorm){
      var +=", normalized";
   }
   
   TLatex label=gfx::cornerLabel(var,1);
   label.Draw();
   TLatex label2=gfx::cornerLabel(method,3);
   if(method != "NoReg") label2.Draw();
   
   //Draw corner label for TWDS option
   bool isTWDS = saver->getInternalPath().Contains("SingleTopDS");
   TLatex label_TW=gfx::cornerLabel((isTWDS)? "tW DS" : "tW DR",2);
   label_TW.SetTextSize(0.045);
   label_TW.Draw();
   
   TLegend leg=legE.buildLegend(.75,.12,0.95,.95,1);
   leg.SetTextSize(0.035);
   leg.Draw();
   
   totalShift.first.Draw("axis same");
   saver->save(can,"systBreakdown/"+name+"_"+method,true,!isData,true);
   
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

TH1F tunfoldplotting::getPDFenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm, bool const &scaleLumi){
   
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
      else if (scaleLumi){
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
                                    bool const &addStat, bool const &withScaleFactor, std::map<TString,TH1F> &indShifts, std::map<TString,TH2F> &indResponse, std::map<TString,TH2F> &indResponseAlt, bool const &norm){
   
   TH1F* hist_shiftUP = (TH1F*)nominal->Clone();
   TH1F* hist_shiftDOWN = (TH1F*)nominal->Clone();
   hist_shiftUP->Reset();
   hist_shiftDOWN->Reset();
   
   TH1F tempSys;
   TH1F tempShift;
   TH2F tempResponse;
   TH2F tempResponseAlt;
   
   TString responsePath = path;
   TString tauPath = path;
   responsePath.ReplaceAll("_results_","_binning_");
   responsePath.ReplaceAll("hist_unfoldedResult","histMCGenRec");
   responsePath.ReplaceAll("realData","TTbar_diLepton");
   responsePath.ReplaceAll("/BBB","");
   responsePath.ReplaceAll("_SingleTopDS","");
   tauPath.ReplaceAll("hist_unfoldedResult","tau");
   bool saveResponse = !(path.Contains("/reg") || path.Contains("/BBB")) || responsePath.Contains("inclusive");     // Response has to be read only once
   bool printTau = path.Contains("/reg");     // Tau parameter only available for reg
      
   for (auto &syst : systVec){
            
      if (Systematic::convertType(syst) == Systematic::CR_envelope) {   // derive CR envelope
         tempSys = getCRenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::mtop){      // derive mtop uncertainty
         tempSys = getMTOPunc(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::meScale_envelope){      // derive meScale uncertainty
         tempSys = getMESCALEenvelope(path,filePath,nominal,Systematic::convertVariation(syst) == Systematic::up,norm);
      }
      else if (Systematic::convertType(syst) == Systematic::pdf_envelope){      // derive pdf uncertainty
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
            tempResponseAlt = *(histReader.read<TH2F>(responsePath+"Alt"));
            indResponse.insert(std::pair<TString,TH2F>(syst,tempResponse));
            indResponseAlt.insert(std::pair<TString,TH2F>(syst,tempResponseAlt));
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
      // ~if (key.BeginsWith("STAT")) continue;
      
      TH1F tempShift = value;
      if(!isNorm) tempShift.Scale(1./scale);
            
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
                                                         std::pair<TH1F*,TH1F*> &realDis_syst_total, std::pair<TH1F*,TH1F*> &realDisAlt_syst_total, std::pair<TH1F*,TH1F*> &realDisHerwig_syst_total, const float &tau_par, TH1F* realDis, TH1F* realDisAlt, TH1F* realDisHerwig, TH2D* cov,const distrUnfold &dist, const bool plotComparison,const TString &saveName, io::RootFileSaver* saver, int &num_bins, const bool rewStudy, const bool onlyTheo, const bool plotTheo, std::ofstream &chi2_file, const bool divideBinWidth, const bool adaptYaxis){
   //========================
   // Step 3: plotting
   
   // ~gfx::SplitCan can;
   gfx::SplitCan can(0.55);
   can.can_.Size(1000,600);
   can.pU_.SetLogy();
   can.pU_.cd();  //Use upper part of canvas
   
   //Check if dPhi distributions and if realData
   bool is_dPhi_1d = (!dist.is2D && dist.varName.Contains("dPhi"));
   bool isRealData = saver->getInternalPath().Contains("realData");
   bool isTWDS = saver->getInternalPath().Contains("SingleTopDS");
   bool isAltMC = saveName.Contains("amcatnlo");
   bool isBSM = saveName.Contains("BSM");
         
   //Get fixed order prediction
   std::pair<TH1F,std::pair<TH1F,TH1F>> fixedOrderNNLOpair;
   std::pair<TH1F,std::pair<TH1F,TH1F>> fixedOrderNLOpair;
   std::pair<TH1F,std::pair<TH1F,TH1F>> fixedOrderLOpair;
   
   if(plotTheo){
      if(dist.is2D){
         fixedOrderNNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NNLO_mt172.5_HT4_NNPDF31_pTnnbar_x_dphilnnbar.dat",true,dist.norm,true);
         fixedOrderNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NLO_mt172.5_HT4_NNPDF31_pTnnbar_x_dphilnnbar.dat",true,dist.norm,false);
         fixedOrderLOpair = readFixedOrderPred("../data/theoryPredictions/fid_LO_mt172.5_HT4_NNPDF31_pTnnbar_x_dphilnnbar.dat",true,dist.norm,false);
      }
      else if (dist.varName.Contains("pTnunu")){
         fixedOrderNNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NNLO_mt172.5_HT4_NNPDF31_pTnnbar.dat",false,dist.norm,true);
         fixedOrderNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NLO_mt172.5_HT4_NNPDF31_pTnnbar.dat",false,dist.norm,false);
         fixedOrderLOpair = readFixedOrderPred("../data/theoryPredictions/fid_LO_mt172.5_HT4_NNPDF31_pTnnbar.dat",false,dist.norm,false);
      }
      else if (dist.varName.Contains("dPhi")){
         fixedOrderNNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NNLO_mt172.5_HT4_NNPDF31_dphilnnbar.dat",false,dist.norm,true,true);
         fixedOrderNLOpair = readFixedOrderPred("../data/theoryPredictions/fid_NLO_mt172.5_HT4_NNPDF31_dphilnnbar.dat",false,dist.norm,false,true);
         fixedOrderLOpair = readFixedOrderPred("../data/theoryPredictions/fid_LO_mt172.5_HT4_NNPDF31_dphilnnbar.dat",false,dist.norm,false,true);
      }
   }
         
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
   Double_t xbins_minus_extra[num_bins+1];
   Double_t xbins_plus_extra[num_bins+1];
   Double_t xbins_plus_extra_2[num_bins+1];
   double minBinWidth = binning_met[2]-binning_met[1];
   if(is_dPhi_1d) minBinWidth = binning_met[1]-binning_met[0];
   xbins[0] = 0;   //Plotting start
   
   //Without Herwig
   // ~xbins_minus[0] = -1.*(minBinWidth/5.);
   // ~xbins_plus[0] = minBinWidth/5.;
   // ~xbins_minus_extra[0] = -1.*(2*minBinWidth/5.);
   // ~xbins_plus_extra[0] = 2*minBinWidth/5.;
   
   //With Herwig
   xbins_minus[0] = -1.*(minBinWidth/7.);
   xbins_plus[0] = minBinWidth/7.;
   xbins_minus_extra[0] = -1.*(2*minBinWidth/7.);
   xbins_plus_extra[0] = 2*minBinWidth/7.;
   xbins_plus_extra_2[0] = 3*minBinWidth/7.;
   
   int phi_bin = 0;
   for (int i=0; i<(num_bins); i++)   {
      xbins[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax;
      //Without Herwig
      // ~xbins_minus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-minBinWidth/5.;
      // ~xbins_plus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+minBinWidth/5.;
      // ~xbins_minus_extra[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-2*minBinWidth/5.;
      // ~xbins_plus_extra[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+2*minBinWidth/5.;
      
      //With Herwig
      xbins_minus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-minBinWidth/6.;
      xbins_plus[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+minBinWidth/6.;
      xbins_minus_extra[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax-2*minBinWidth/7.;
      xbins_plus_extra[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+2*minBinWidth/7.;
      xbins_plus_extra_2[i+1] = binning_met[i%num_bins_met+1]+phi_bin*dist.xMax+3*minBinWidth/7.;
      if (i%num_bins_met==num_bins_met-1) phi_bin++;
   }
   
   std::vector<double> xbins_vec(xbins, xbins + sizeof xbins / sizeof xbins[0]);
   
   unfolded->SetBins(num_bins,xbins);
   unfolded_total.first->SetBins(num_bins,xbins);
   unfolded_total.second->SetBins(num_bins,xbins);
   realDis_syst_total.first->SetBins(num_bins,xbins);
   realDis_syst_total.second->SetBins(num_bins,xbins);
   realDisAlt_syst_total.first->SetBins(num_bins,xbins);
   realDisAlt_syst_total.second->SetBins(num_bins,xbins);
   realDisHerwig_syst_total.first->SetBins(num_bins,xbins);
   realDisHerwig_syst_total.second->SetBins(num_bins,xbins);
   realDis->SetBins(num_bins,xbins);
   realDisAlt->SetBins(num_bins,xbins);
   realDisHerwig->SetBins(num_bins,xbins);
   cov->SetBins(num_bins,xbins,num_bins,xbins);
   if(dist.is2D && plotTheo){
      fixedOrderNNLOpair.first.SetBins(num_bins,xbins);
      fixedOrderNLOpair.first.SetBins(num_bins,xbins);
      fixedOrderLOpair.first.SetBins(num_bins,xbins);
      fixedOrderNNLOpair.second.first.SetBins(num_bins,xbins);
      fixedOrderNLOpair.second.first.SetBins(num_bins,xbins);
      fixedOrderLOpair.second.first.SetBins(num_bins,xbins);
      fixedOrderNNLOpair.second.second.SetBins(num_bins,xbins);
      fixedOrderNLOpair.second.second.SetBins(num_bins,xbins);
      fixedOrderLOpair.second.second.SetBins(num_bins,xbins);
   }
   if (plotComparison){
      unfolded_reg->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.first->SetBins(num_bins,xbins_minus);
      unfolded_reg_total.second->SetBins(num_bins,xbins_minus);
      unfolded_bbb->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.first->SetBins(num_bins,xbins_plus);
      unfolded_bbb_total.second->SetBins(num_bins,xbins_plus);
   }
   else if (!onlyTheo){ //set offset for proper theory comparison
      realDis->SetBins(num_bins,xbins_minus);
      realDis_syst_total.first->SetBins(num_bins,xbins_minus);
      realDis_syst_total.second->SetBins(num_bins,xbins_minus);
      realDisAlt_syst_total.first->SetBins(num_bins,xbins_minus_extra);
      realDisAlt_syst_total.second->SetBins(num_bins,xbins_minus_extra);
      realDisAlt->SetBins(num_bins,xbins_minus_extra);
      realDisHerwig_syst_total.first->SetBins(num_bins,xbins_plus_extra_2);
      realDisHerwig_syst_total.second->SetBins(num_bins,xbins_plus_extra_2);
      realDisHerwig->SetBins(num_bins,xbins_plus_extra_2);
      fixedOrderNNLOpair.first.SetBins(num_bins,xbins_plus_extra);
      fixedOrderNLOpair.first.SetBins(num_bins,xbins_plus);
      fixedOrderNNLOpair.second.first.SetBins(num_bins,xbins_plus_extra);
      fixedOrderNLOpair.second.first.SetBins(num_bins,xbins_plus);
      fixedOrderNNLOpair.second.second.SetBins(num_bins,xbins_plus_extra);
      fixedOrderNLOpair.second.second.SetBins(num_bins,xbins_plus);
   }

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
      realDis_syst_total.first->GetXaxis()->SetBinLabel(i,label);
      realDis_syst_total.second->GetXaxis()->SetBinLabel(i,label);
      realDisAlt_syst_total.first->GetXaxis()->SetBinLabel(i,label);
      realDisAlt_syst_total.second->GetXaxis()->SetBinLabel(i,label);
      realDisAlt->GetXaxis()->SetBinLabel(i,label);
      realDisHerwig_syst_total.first->GetXaxis()->SetBinLabel(i,label);
      realDisHerwig_syst_total.second->GetXaxis()->SetBinLabel(i,label);
      realDisHerwig->GetXaxis()->SetBinLabel(i,label);
      labelVec.push_back(label);
   }
   
   // Divide by bin width
   if(divideBinWidth){
      hist::divideByBinWidth(*realDis);
      hist::divideByBinWidth(*realDisAlt);
      hist::divideByBinWidth(*realDisHerwig);
      hist::divideByBinWidth(*unfolded_reg);
      hist::divideByBinWidth(*unfolded_reg_total.first);
      hist::divideByBinWidth(*unfolded_reg_total.second);
      hist::divideByBinWidth(*unfolded_bbb);
      hist::divideByBinWidth(*unfolded_bbb_total.first);
      hist::divideByBinWidth(*unfolded_bbb_total.second);
      hist::divideByBinWidth(*unfolded);
      hist::divideByBinWidth(*unfolded_total.first);
      hist::divideByBinWidth(*unfolded_total.second);
      hist::divideByBinWidth(*realDis_syst_total.first);
      hist::divideByBinWidth(*realDis_syst_total.second);
      hist::divideByBinWidth(*realDisAlt_syst_total.first);
      hist::divideByBinWidth(*realDisAlt_syst_total.second);
      hist::divideByBinWidth(*realDisHerwig_syst_total.first);
      hist::divideByBinWidth(*realDisHerwig_syst_total.second);
      hist::divideByBinWidth(*cov);
      if(dist.is2D && plotTheo){
         hist::divideByBinWidth(fixedOrderNNLOpair.first);
         hist::divideByBinWidth(fixedOrderNLOpair.first);
         hist::divideByBinWidth(fixedOrderLOpair.first);
         hist::divideByBinWidth(fixedOrderNNLOpair.second.first);
         hist::divideByBinWidth(fixedOrderNLOpair.second.first);
         hist::divideByBinWidth(fixedOrderLOpair.second.first);
         hist::divideByBinWidth(fixedOrderNNLOpair.second.second);
         hist::divideByBinWidth(fixedOrderNLOpair.second.second);
         hist::divideByBinWidth(fixedOrderLOpair.second.second);
      }
   }
   if(dist.is2D && plotTheo && saveName.Contains("Combined")){    //Needed for combined 2D results
      hist::divideByBinWidth(fixedOrderNNLOpair.first);
      hist::divideByBinWidth(fixedOrderNLOpair.first);
      hist::divideByBinWidth(fixedOrderLOpair.first);
      hist::divideByBinWidth(fixedOrderNNLOpair.second.first);
      hist::divideByBinWidth(fixedOrderNLOpair.second.first);
      hist::divideByBinWidth(fixedOrderLOpair.second.first);
      hist::divideByBinWidth(fixedOrderNNLOpair.second.second);
      hist::divideByBinWidth(fixedOrderNLOpair.second.second);
      hist::divideByBinWidth(fixedOrderLOpair.second.second);
   }
   
   // Derive Graph for fixed order uncertainties
   TGraphAsymmErrors fixedOrderNNLOgraph = hist::getErrorGraph(&fixedOrderNNLOpair.second.first,&fixedOrderNNLOpair.second.second,&fixedOrderNNLOpair.first,true,!onlyTheo);
   TGraphAsymmErrors fixedOrderNLOgraph = hist::getErrorGraph(&fixedOrderNLOpair.second.first,&fixedOrderNLOpair.second.second,&fixedOrderNLOpair.first,true,!onlyTheo);
   TGraphAsymmErrors fixedOrderLOgraph = hist::getErrorGraph(&fixedOrderLOpair.second.first,&fixedOrderLOpair.second.second,&fixedOrderLOpair.first,true,!onlyTheo);
   
   // Plotting options
   unfolded->GetXaxis()->SetTickLength(0.);
   unfolded->GetYaxis()->SetTickLength(0.008);
   unfolded->GetXaxis()->SetTitleOffset(1.5);
   unfolded->GetYaxis()->SetTitleSize(0.045);
   unfolded->GetXaxis()->CenterLabels(false);
      
   unfolded->LabelsOption("v");
   realDis->LabelsOption("v");
   realDisAlt->LabelsOption("v");
   if(adaptYaxis){
      unfolded->SetMaximum(2.5*unfolded->GetMaximum());
      // ~unfolded->SetMinimum(0.4*unfolded->GetMinimum());
      if(!is_dPhi_1d && !dist.is2D) unfolded->SetMinimum(0.15*unfolded->GetMinimum());
      else unfolded->SetMinimum((isBSM)? 0.1*unfolded->GetMinimum() : 0.3*unfolded->GetMinimum());
   }
   unfolded->SetLineColor(kBlack);
   unfolded->SetTitle(dist.title);
   unfolded->SetStats(false);
   realDis->SetTitle(dist.title);
   realDisAlt->SetTitle(dist.title);
   realDis->SetStats(false);
   
   unfolded_reg->SetLineColor(kGreen+2);
   unfolded_bbb->SetLineColor(kViolet);
   unfolded_reg->SetMarkerColor(kGreen+2);
   unfolded_bbb->SetMarkerColor(kViolet);
   
   // Setup unc. plotting
   // ~TGraphAsymmErrors unfolded_totalGraph = hist::getErrorGraph(unfolded_total.first,unfolded_total.second,unfolded,true,false);
   TGraphAsymmErrors unfolded_totalGraph = hist::getErrorGraph(unfolded_total.first,unfolded_total.second,unfolded,true,plotComparison);
   TGraphAsymmErrors unfolded_reg_totalGraph = hist::getErrorGraph(unfolded_reg_total.first,unfolded_reg_total.second,unfolded_reg,true,true);
   TGraphAsymmErrors unfolded_bbb_totalGraph = hist::getErrorGraph(unfolded_bbb_total.first,unfolded_bbb_total.second,unfolded_bbb,true,true);
   TGraphAsymmErrors realDis_totalGraph = hist::getErrorGraph(realDis_syst_total.first,realDis_syst_total.second,realDis,true,true);
   TGraphAsymmErrors realDisAlt_totalGraph = hist::getErrorGraph(realDisAlt_syst_total.first,realDisAlt_syst_total.second,realDisAlt,true,true);
   TGraphAsymmErrors realDisHerwig_totalGraph = hist::getErrorGraph(realDisHerwig_syst_total.first,realDisHerwig_syst_total.second,realDisHerwig,true,true);
   unfolded_totalGraph.SetLineColor(kBlack);
   unfolded_reg_totalGraph.SetLineColor(kGreen+2);
   unfolded_bbb_totalGraph.SetLineColor(kViolet);
   TGraphErrors unfolded_graph(unfolded);
   TGraphErrors unfolded_reg_graph(unfolded_reg);
   TGraphErrors unfolded_bbb_graph(unfolded_bbb);
   if (plotComparison){
      for (int i=0; i<unfolded->GetNbinsX(); i++){    //change x error for plotting
         unfolded_graph.SetPointError(i,binning_met[1]/10.,unfolded->GetBinError(i+1));
         unfolded_reg_graph.SetPointError(i,binning_met[1]/10.,unfolded_reg->GetBinError(i+1));
         unfolded_bbb_graph.SetPointError(i,binning_met[1]/10.,unfolded_bbb->GetBinError(i+1));
      }
   }
   unfolded_graph.SetFillStyle(1001);
   unfolded_reg_graph.SetFillStyle(1001);
   unfolded_bbb_graph.SetFillStyle(1001);
   unfolded_graph.SetLineWidth(0);
   unfolded_reg_graph.SetLineWidth(0);
   unfolded_bbb_graph.SetLineWidth(0);
   unfolded_graph.SetFillColor(kGray+1);
   unfolded_reg_graph.SetFillColor(kGreen-9);
   unfolded_bbb_graph.SetFillColor(kMagenta-9);
   unfolded_graph.SetMarkerSize(0.4);
   unfolded_reg_graph.SetMarkerSize(0.4);
   unfolded_bbb_graph.SetMarkerSize(0.4);
   
   unfolded_totalGraph.SetFillStyle(1001);
   unfolded_totalGraph.SetFillColor(kOrange+1);
   
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
   realDisHerwig->SetLineColor(kCyan-6);
   realDisHerwig->SetFillColor(kCyan-6);
   
   fixedOrderNNLOpair.first.SetLineColor(kRed);
   fixedOrderNLOpair.first.SetLineColor(kBlue);
   fixedOrderLOpair.first.SetLineColor(kGreen);
   
   fixedOrderNNLOgraph.SetFillStyle(1001);
   fixedOrderNNLOgraph.SetFillColorAlpha(kRed, 0.35);
   fixedOrderNLOgraph.SetFillStyle(1001);
   fixedOrderNLOgraph.SetFillColorAlpha(kBlue, 0.35);
   fixedOrderLOgraph.SetFillStyle(1001);
   fixedOrderLOgraph.SetFillColorAlpha(kGreen, 0.35);
      
   if (onlyTheo){
      unfolded->Draw("axis");
      fixedOrderLOgraph.Draw("same 2");
      fixedOrderLOpair.first.Draw("hist same");
      fixedOrderNLOgraph.Draw("same 2");
      fixedOrderNLOpair.first.Draw("hist same");
      fixedOrderNNLOgraph.Draw("same 2");
      fixedOrderNNLOpair.first.Draw("hist same");
      
      realDis->SetLineStyle(7);
      realDisAlt->SetLineStyle(7);
      realDis->Draw("hist same");
      realDisAlt->Draw("hist same");
   }
   else if (plotComparison) {
      
      // first draw used for axis and stuff
      // ~unfolded->Draw("px0");
      unfolded->Draw("axis");
      
      
      if(!isAltMC && !isBSM){
         realDis->Draw("hist same");   //Draw into current canvas (axis are not drawn again due to "same")
         if(rewStudy || isBSM) realDisAlt->Draw("hist same");
      }
      else{
         realDis->SetLineColor(kBlue-6);
         realDisAlt->SetLineColor(kRed-6);
         realDis->Draw("hist same");
         realDisAlt->Draw("hist same");
      }
      if(plotTheo){
         fixedOrderNNLOpair.first.Draw("hist same");
         fixedOrderNLOpair.first.Draw("hist same");
         fixedOrderLOpair.first.Draw("hist same");
      }
      
      //draw stat. unc.
      // ~unfolded_graph.Draw("pe2 same");
      // ~unfolded_reg_graph.Draw("pe2 same");
      // ~unfolded_bbb_graph.Draw("pe2 same");
      
      // draw total unc.
      unfolded_totalGraph.Draw("pe same");
      unfolded_reg_totalGraph.Draw("pe same");
      unfolded_bbb_totalGraph.Draw("pe same");
   }
   else {   //nominal with real data
      
      
      realDis_totalGraph.SetMarkerSize(0.5);
      realDisAlt_totalGraph.SetMarkerSize(0.5);
      realDisHerwig_totalGraph.SetMarkerSize(0.5);
      realDis_totalGraph.SetMarkerColor(kGreen+2);
      realDis_totalGraph.SetLineColor(kGreen+2);
      realDisAlt_totalGraph.SetMarkerColor(kViolet);
      realDisAlt_totalGraph.SetLineColor(kViolet);
      realDisHerwig_totalGraph.SetMarkerColor(kCyan+2);
      realDisHerwig_totalGraph.SetLineColor(kCyan+2);
      realDis_totalGraph.SetMarkerStyle(22);
      realDisAlt_totalGraph.SetMarkerStyle(23);
      realDisHerwig_totalGraph.SetMarkerStyle(43);
      realDis_totalGraph.SetLineWidth(1);
      realDisAlt_totalGraph.SetLineWidth(1);
      realDisHerwig_totalGraph.SetLineWidth(1);
            
      unfolded->Draw("axis");  //Draw into current canvas
      
      unfolded_totalGraph.SetLineWidth(0);
      
      unfolded_totalGraph.Draw("2 same");
      unfolded_graph.Draw("2 same");
      unfolded_totalGraph.Draw("pX same");
      
      realDis_totalGraph.Draw("pe same");
      realDisAlt_totalGraph.Draw("pe same");
      realDisHerwig_totalGraph.Draw("pe same");
      
      if(plotTheo){
         fixedOrderNNLOgraph.SetMarkerSize(0.5);
         fixedOrderNLOgraph.SetMarkerSize(0.5);
         fixedOrderNNLOgraph.SetMarkerColor(kRed);
         fixedOrderNLOgraph.SetMarkerColor(kBlue);
         fixedOrderNNLOgraph.SetLineColor(kRed);
         fixedOrderNLOgraph.SetLineColor(kBlue);
         fixedOrderNNLOgraph.SetMarkerStyle(34);
         fixedOrderNLOgraph.SetMarkerStyle(33);
         fixedOrderNNLOgraph.SetLineWidth(1);
         fixedOrderNLOgraph.SetLineWidth(1);
         
         gStyle->SetEndErrorSize(0);      // do not plot small lines at the end of errorbar
      
         fixedOrderNNLOgraph.Draw("pe same");
         fixedOrderNLOgraph.Draw("pe same");
      }
      
   }
   
   
   TLatex * atext = new TLatex();
   TLine * aline = new TLine();
   if (dist.is2D){
      //Draw vertical lines and binning ranges for deltaPhi
      atext->SetTextSize(0.035);
      aline->SetLineWidth(2);
      aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
      aline->DrawLine(400,unfolded->GetMinimum(),400,unfolded->GetMaximum());
      aline->DrawLine(800,unfolded->GetMinimum(),800,unfolded->GetMaximum());
      atext->DrawLatex(50,0.5*unfolded->GetMaximum(),TString::Format("0<|#Delta#phi|<%.2f",binning_phi(1)));
      atext->DrawLatex(450,0.5*unfolded->GetMaximum(),TString::Format("%.2f<|#Delta#phi|<%.2f",binning_phi(1),binning_phi(2)));
      atext->DrawLatex(850,0.5*unfolded->GetMaximum(),TString::Format("%.2f<|#Delta#phi|<%.2f",binning_phi(2),binning_phi(3)));
   }
   
   //Get Chi2 and NDF
   auto Chi2Pair = getChi2NDF(unfolded,(isAltMC || isBSM)? realDisAlt:realDis);
   // ~auto Chi2Pair_corr = getChi2NDF_withCorr(unfolded,realDis,covMatrix);
   
   //Draw legend
   gfx::LegendEntries legE;
   if (plotComparison && !onlyTheo) {
      auto Chi2Pair_reg = getChi2NDF(unfolded_reg,(isAltMC || isBSM)? realDisAlt:realDis);
      auto Chi2Pair_bbb = getChi2NDF(unfolded_bbb,(isAltMC || isBSM)? realDisAlt:realDis);
      // ~legE.append(*unfolded,TString::Format("NoReg [#chi^{2}/NDF=%.1f/%i]",Chi2Pair.first,Chi2Pair.second),"pe");
      // ~legE.append(*unfolded_reg,TString::Format("Reg [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_reg.first,Chi2Pair_reg.second),"pe");
      // ~legE.append(*unfolded_bbb,TString::Format("BBB [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_bbb.first,Chi2Pair_bbb.second),"pe");
      if(isAltMC and !isBSM){
         legE.append(*realDisAlt,"Truth","l");
         legE.append(*realDis,"Response","l");
         // ~legE.append(unfolded_totalGraph,"Nominal","pe");
         // ~legE.append(unfolded_reg_totalGraph,"Regularized","pe");
         // ~legE.append(unfolded_bbb_totalGraph,"BBB","pe");
         legE.append(unfolded_totalGraph,TString::Format("Nominal (#chi^{2}/dof=%.1f/%i)",Chi2Pair.first,Chi2Pair.second),"pe");
         legE.append(unfolded_reg_totalGraph,TString::Format("Regularized (#chi^{2}/dof=%.1f/%i)",Chi2Pair_reg.first,Chi2Pair_reg.second),"pe");
         legE.append(unfolded_bbb_totalGraph,TString::Format("BBB (#chi^{2}/dof=%.1f/%i)",Chi2Pair_bbb.first,Chi2Pair_bbb.second),"pe");
      }
      else if (rewStudy){
         legE.append(*realDis,"Truth","l");
         legE.append(*realDisAlt,"Response","l");
         legE.append(unfolded_totalGraph,TString::Format("Nominal (#chi^{2}/dof=%.1f/%i)",Chi2Pair.first,Chi2Pair.second),"pe");
         legE.append(unfolded_reg_totalGraph,TString::Format("Regularized (#chi^{2}/dof=%.1f/%i)",Chi2Pair_reg.first,Chi2Pair_reg.second),"pe");
         legE.append(unfolded_bbb_totalGraph,TString::Format("BBB (#chi^{2}/dof=%.1f/%i)",Chi2Pair_bbb.first,Chi2Pair_bbb.second),"pe");
      }
      else if (isBSM){
         legE.append(*realDis,"Response","l");
         legE.append(*realDisAlt,"tt+BSM (p_{T}^{#nu#nu}=p_{T}^{#nu#nu+BSM})","l");
         legE.append(unfolded_totalGraph,TString::Format("Nominal [%.1f/%i]",Chi2Pair.first,Chi2Pair.second),"pe");
         legE.append(unfolded_reg_totalGraph,TString::Format("Regularized [%.1f/%i]",Chi2Pair_reg.first,Chi2Pair_reg.second),"pe");
         legE.append(unfolded_bbb_totalGraph,TString::Format("BBB [%.1f/%i]",Chi2Pair_bbb.first,Chi2Pair_bbb.second),"pe");
      }
      else {
         legE.append(*realDis,"Truth","l");
         legE.append(unfolded_totalGraph,"Nominal","pe");
         legE.append(unfolded_reg_totalGraph,"Regularized","pe");
         legE.append(unfolded_bbb_totalGraph,"BBB","pe");
      }
   }
   else if (!onlyTheo) { //nominal with real data (save chi2 values to file)
      
      auto Chi2Pair_powheg = getChi2NDF(&unfolded_totalGraph,realDis);
      auto Chi2Pair_madgraph = getChi2NDF(&unfolded_totalGraph,realDisAlt);
      auto Chi2Pair_herwig = getChi2NDF(&unfolded_totalGraph,realDisHerwig);
      auto Chi2Pair_NNLO = getChi2NDF(&unfolded_totalGraph,&fixedOrderNNLOgraph);
      auto Chi2Pair_NLO = getChi2NDF(&unfolded_totalGraph,&fixedOrderNLOgraph);
      
      auto Chi2Pair_powheg_theo = getChi2NDF(&unfolded_totalGraph,&realDis_totalGraph,true);
      auto Chi2Pair_madgraph_theo = getChi2NDF(&unfolded_totalGraph,&realDisAlt_totalGraph,true);
      auto Chi2Pair_herwig_theo = getChi2NDF(&unfolded_totalGraph,&realDisHerwig_totalGraph,true);
      auto Chi2Pair_NNLO_theo = getChi2NDF(&unfolded_totalGraph,&fixedOrderNNLOgraph,true);
      auto Chi2Pair_NLO_theo = getChi2NDF(&unfolded_totalGraph,&fixedOrderNLOgraph,true);
      
      auto Chi2Pair_powheg_cov = getChi2NDF_withCorr(unfolded,&realDis_totalGraph,cov,dist.norm);
      auto Chi2Pair_madgraph_cov = getChi2NDF_withCorr(unfolded,&realDisAlt_totalGraph,cov,dist.norm);
      auto Chi2Pair_herwig_cov = getChi2NDF_withCorr(unfolded,&realDisHerwig_totalGraph,cov,dist.norm);
      auto Chi2Pair_NNLO_cov = getChi2NDF_withCorr(unfolded,&fixedOrderNNLOgraph,cov,dist.norm);
      auto Chi2Pair_NLO_cov = getChi2NDF_withCorr(unfolded,&fixedOrderNLOgraph,cov,dist.norm);
      
      auto Chi2Pair_powheg_cov_theo = getChi2NDF_withCorr(unfolded,&realDis_totalGraph,cov,dist.norm,true);
      auto Chi2Pair_madgraph_cov_theo = getChi2NDF_withCorr(unfolded,&realDisAlt_totalGraph,cov,dist.norm,true);
      auto Chi2Pair_herwig_cov_theo = getChi2NDF_withCorr(unfolded,&realDisHerwig_totalGraph,cov,dist.norm,true);
      auto Chi2Pair_NNLO_cov_theo = getChi2NDF_withCorr(unfolded,&fixedOrderNNLOgraph,cov,dist.norm,true);
      auto Chi2Pair_NLO_cov_theo = getChi2NDF_withCorr(unfolded,&fixedOrderNLOgraph,cov,dist.norm,true);
      
      // store chi2 values
      chi2_file<<"--------------------------------------------------------------------------"<<std::endl;
      chi2_file<<dist.varName<<" chi2";
      if(dist.norm) chi2_file<<" normalized";
      chi2_file<<std::endl;
      chi2_file<<std::setprecision(1)<<std::fixed;
      chi2_file<<"NDOF  "<<Chi2Pair_powheg_cov.second<<std::endl;
      chi2_file<<"POWHEG  "<<Chi2Pair_powheg.first<<"   "<<Chi2Pair_powheg_theo.first<<"   "<<Chi2Pair_powheg_cov.first<<"   "<<Chi2Pair_powheg_cov_theo.first<<std::endl;
      chi2_file<<"MADGRAPH  "<<Chi2Pair_madgraph.first<<"   "<<Chi2Pair_madgraph_theo.first<<"   "<<Chi2Pair_madgraph_cov.first<<"   "<<Chi2Pair_madgraph_cov_theo.first<<std::endl;
      chi2_file<<"HERWIG  "<<Chi2Pair_herwig.first<<"   "<<Chi2Pair_herwig_theo.first<<"   "<<Chi2Pair_herwig_cov.first<<"   "<<Chi2Pair_herwig_cov_theo.first<<std::endl;
      chi2_file<<"NNLO  "<<Chi2Pair_NNLO.first<<"   "<<Chi2Pair_NNLO_theo.first<<"   "<<Chi2Pair_NNLO_cov.first<<"   "<<Chi2Pair_NNLO_cov_theo.first<<std::endl;
      chi2_file<<"NLO  "<<Chi2Pair_NLO.first<<"   "<<Chi2Pair_NLO_theo.first<<"   "<<Chi2Pair_NLO_cov.first<<"   "<<Chi2Pair_NLO_cov_theo.first<<std::endl<<std::endl;
      
      // store p values values
      int nBins = Chi2Pair_powheg_cov.second;
      chi2_file<<"--------------------------------------------------------------------------"<<std::endl;
      chi2_file<<dist.varName<<" p-value";
      if(dist.norm) chi2_file<<" normalized";
      chi2_file<<std::endl;
      chi2_file<<std::setprecision(3)<<std::fixed;
      chi2_file<<"NDOF  "<<nBins<<std::endl;
      chi2_file<<"POWHEG  "<<TMath::Prob(Chi2Pair_powheg.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_powheg_theo.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_powheg_cov.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_powheg_cov_theo.first,nBins)<<std::endl;
      chi2_file<<"MADGRAPH  "<<TMath::Prob(Chi2Pair_madgraph.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_madgraph_theo.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_madgraph_cov.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_madgraph_cov_theo.first,nBins)<<std::endl;
      chi2_file<<"HERWIG  "<<TMath::Prob(Chi2Pair_herwig.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_herwig_theo.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_herwig_cov.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_herwig_cov_theo.first,nBins)<<std::endl;
      chi2_file<<"NNLO  "<<TMath::Prob(Chi2Pair_NNLO.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NNLO_theo.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NNLO_cov.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NNLO_cov_theo.first,nBins)<<std::endl;
      chi2_file<<"NLO  "<<TMath::Prob(Chi2Pair_NLO.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NLO_theo.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NLO_cov.first,nBins)<<"   "<<TMath::Prob(Chi2Pair_NLO_cov_theo.first,nBins)<<std::endl<<std::endl;
            
      // ~for (int i=0; i<unfolded_totalGraph.GetN()-1; i++) std::cout<<unfolded_totalGraph.GetErrorY(i)*unfolded_totalGraph.GetErrorY(i)<<std::endl;
                  
      legE.append(unfolded_totalGraph,"Data","p");
      legE.append(unfolded_totalGraph,"#sigma_{tot.}","f");
      legE.append(unfolded_graph,"#sigma_{stat.}","f");
      // ~legE.append(realDis_totalGraph,"Powheg","p");
      // ~legE.append(realDisAlt_totalGraph,"MadGraph","p");
      // ~legE.append(fixedOrderNNLOgraph,"NNLO","p");
      // ~legE.append(fixedOrderNLOgraph,"NLO","p");
      legE.append(realDis_totalGraph,TString::Format("Powheg [%.1f/%i]",Chi2Pair_powheg_cov.first,Chi2Pair_powheg_cov.second),"p");
      legE.append(realDisAlt_totalGraph,TString::Format("MC@NLO [%.1f/%i]",Chi2Pair_madgraph_cov.first,Chi2Pair_madgraph_cov.second),"p");
      legE.append(realDisHerwig_totalGraph,TString::Format("HERWIG [%.1f/%i]",Chi2Pair_herwig_cov.first,Chi2Pair_herwig_cov.second),"p");
      legE.append(fixedOrderNNLOgraph,TString::Format("NNLO [%.1f/%i]",Chi2Pair_NNLO_cov.first,Chi2Pair_NNLO_cov.second),"p");
      legE.append(fixedOrderNLOgraph,TString::Format("NLO [%.1f/%i]",Chi2Pair_NLO_cov.first,Chi2Pair_NLO_cov.second),"p");
      // ~legE.append(*realDis,TString::Format("Powheg [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_powheg.first,Chi2Pair_powheg.second),"p");
      // ~legE.append(*realDisAlt,TString::Format("MadGraph [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_madgraph.first,Chi2Pair_madgraph.second),"p");
      // ~legE.append(fixedOrderNNLOgraph,TString::Format("NNLO [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_NNLO.first,Chi2Pair_NNLO.second),"p");
      // ~legE.append(fixedOrderNLOgraph,TString::Format("NLO [#chi^{2}/NDF=%.1f/%i]",Chi2Pair_NLO.first,Chi2Pair_NLO.second),"p");
   }
   
   if(!rewStudy && onlyTheo){
      legE.append(*realDis,"Powheg","l");
      legE.append(*realDisAlt,"MadGraph","l");
      legE.append(fixedOrderNNLOpair.first,"NNLO","l");
      legE.append(fixedOrderNLOpair.first,"NLO","l");
      legE.append(fixedOrderLOpair.first,"LO","l");
   }
   
   TLegend leg=legE.buildLegend(.149,.45,0.339,(isBSM)? 0.68 : .72,1);
   leg.SetTextSize(0.03);
   leg.Draw();
   
   unfolded->Draw("axis same");
   
   //Draw corner label for TWDS option
   TLatex label=gfx::cornerLabel((isTWDS)? "tW DS" : "tW DR",2);
   label.SetTextSize(0.045);
   label.Draw();
   
   //Change to lower part of canvas
   can.pL_.cd();
   can.pL_.SetBottomMargin(0.45);
   can.pL_.SetTickx(0);
   TH1F ratio;
   TH1F ratio_alt;
   TH1F ratio_herwig;
   TH1F ratio_data;
   TH1F ratio_data_axis;
   TH1F ratio_alt_data;
   TH1F ratio_herwig_data;
   TH1F ratio_NNLO;
   TH1F ratio_NLO;
   TH1F ratio_LO;
   TH1F ratio_NNLO_data;
   TH1F ratio_NLO_data;
   TH1F ratio_LO_data;
   TH1F ratio_unfolded;
   TH1F ratio_unfolded_data;
   TH1F ratio_unfolded_reg;
   TH1F ratio_unfolded_bbb;
   TGraphAsymmErrors ratio_NNLO_graph;
   TGraphAsymmErrors ratio_NLO_graph;
   TGraphAsymmErrors ratio_LO_graph;
   TGraphAsymmErrors ratio_graph_asym;
   TGraphAsymmErrors ratio_alt_graph_asym;
   TGraphAsymmErrors ratio_herwig_graph_asym;
   TGraphAsymmErrors ratio_NNLO_data_graph;
   TGraphAsymmErrors ratio_NLO_data_graph;
   TGraphAsymmErrors ratio_LO_data_graph;
   TGraphAsymmErrors ratio_data_graph_asym;
   TGraphAsymmErrors ratio_alt_data_graph_asym;
   TGraphAsymmErrors ratio_herwig_data_graph_asym;
   
   //change realDis if alt MC is used (easier than using multiple if statements)
   if(isAltMC || isBSM){
      TH1F* dummyHist = realDis;
      realDis = realDisAlt;
      realDisAlt = dummyHist;
   }
   
   // derive ratios
   ratio=hist::getRatio(*realDis,*realDis,"Result/Truth",hist::NOERR);   //Get Ratio between unfolded and true hists
   ratio_data_axis=hist::getRatio(*unfolded,*unfolded,"Pred./Data",hist::NOERR);   
   ratio_data=hist::getRatio(*realDis,*unfolded,"Pred./Data",hist::NOERR);   
   ratio_alt=hist::getRatio(*realDisAlt,*realDis,"Result/POWHEG",hist::NOERR);
   ratio_alt_data=hist::getRatio(*realDisAlt,*unfolded,"Pred./Data",hist::NOERR);
   ratio_herwig=hist::getRatio(*realDisHerwig,*realDis,"Result/POWHEG",hist::NOERR);
   ratio_herwig_data=hist::getRatio(*realDisHerwig,*unfolded,"Pred./Data",hist::NOERR);
   ratio_graph_asym = hist::getRatioAsymmGraph(*realDis_syst_total.first,*realDis_syst_total.second,*realDis,*realDis,false);
   ratio_data_graph_asym = hist::getRatioAsymmGraph(*realDis_syst_total.first,*realDis_syst_total.second,*realDis,*unfolded,!onlyTheo);
   ratio_alt_graph_asym = hist::getRatioAsymmGraph(*realDisAlt_syst_total.first,*realDisAlt_syst_total.second,*realDisAlt,*realDis,false);
   ratio_alt_data_graph_asym = hist::getRatioAsymmGraph(*realDisAlt_syst_total.first,*realDisAlt_syst_total.second,*realDisAlt,*unfolded,!onlyTheo);
   ratio_herwig_graph_asym = hist::getRatioAsymmGraph(*realDisHerwig_syst_total.first,*realDisHerwig_syst_total.second,*realDisHerwig,*realDis,false);
   ratio_herwig_data_graph_asym = hist::getRatioAsymmGraph(*realDisHerwig_syst_total.first,*realDisHerwig_syst_total.second,*realDisHerwig,*unfolded,!onlyTheo);
   if(plotTheo){
      ratio_NNLO=hist::getRatio(fixedOrderNNLOpair.first,*realDis,"Result/POWHEG",hist::ONLY1);
      ratio_NLO=hist::getRatio(fixedOrderNLOpair.first,*realDis,"Result/POWHEG",hist::ONLY1);
      ratio_LO=hist::getRatio(fixedOrderLOpair.first,*realDis,"Result/POWHEG",hist::ONLY1);
      ratio_NNLO_data=hist::getRatio(fixedOrderNNLOpair.first,*unfolded,"Pred./Data",hist::ONLY1);
      ratio_NLO_data=hist::getRatio(fixedOrderNLOpair.first,*unfolded,"Pred./Data",hist::ONLY1);
      ratio_LO_data=hist::getRatio(fixedOrderLOpair.first,*unfolded,"Pred./Data",hist::ONLY1);
      
      ratio_NNLO_graph = hist::getRatioAsymmGraph(fixedOrderNNLOpair.second.first,fixedOrderNNLOpair.second.second,fixedOrderNNLOpair.first,*realDis,false);
      ratio_NLO_graph = hist::getRatioAsymmGraph(fixedOrderNLOpair.second.first,fixedOrderNLOpair.second.second,fixedOrderNLOpair.first,*realDis,false);
      ratio_LO_graph = hist::getRatioAsymmGraph(fixedOrderLOpair.second.first,fixedOrderLOpair.second.second,fixedOrderLOpair.first,*realDis,false);
      ratio_NNLO_data_graph = hist::getRatioAsymmGraph(fixedOrderNNLOpair.second.first,fixedOrderNNLOpair.second.second,fixedOrderNNLOpair.first,*unfolded,!onlyTheo);
      ratio_NLO_data_graph = hist::getRatioAsymmGraph(fixedOrderNLOpair.second.first,fixedOrderNLOpair.second.second,fixedOrderNLOpair.first,*unfolded,!onlyTheo);
      ratio_LO_data_graph = hist::getRatioAsymmGraph(fixedOrderLOpair.second.first,fixedOrderLOpair.second.second,fixedOrderLOpair.first,*unfolded,!onlyTheo);
   }
   ratio_unfolded=hist::getRatio(*unfolded,*realDis,"Result/POWHEG",hist::ONLY1);
   ratio_unfolded_data=hist::getRatio(*unfolded,*unfolded,"Pred./Data",hist::ONLY1);
   ratio_unfolded_reg=hist::getRatio(*unfolded_reg,*realDis,"Result/POWHEG",hist::ONLY1);
   ratio_unfolded_bbb=hist::getRatio(*unfolded_bbb,*realDis,"Result/POWHEG",hist::ONLY1);
      
   // derive syst. ratios
   // ~TGraphAsymmErrors ratio_totalGraph = hist::getRatioAsymmGraph(*unfolded_total.first,*unfolded_total.second,*unfolded,*realDis,false);
   TGraphAsymmErrors ratio_totalGraph = hist::getRatioAsymmGraph(*unfolded_total.first,*unfolded_total.second,*unfolded,*realDis,true);
   TGraphAsymmErrors ratio_totalGraph_data = hist::getRatioAsymmGraph(*unfolded_total.first,*unfolded_total.second,*unfolded,*unfolded,false);
   TGraphAsymmErrors ratio_reg_totalGraph = hist::getRatioAsymmGraph(*unfolded_reg_total.first,*unfolded_reg_total.second,*unfolded_reg,*realDis);
   TGraphAsymmErrors ratio_bbb_totalGraph = hist::getRatioAsymmGraph(*unfolded_bbb_total.first,*unfolded_bbb_total.second,*unfolded_bbb,*realDis);
      
   // setup stat. unc. plotting
   TGraphErrors ratio_graph(&ratio_unfolded);
   TGraphErrors ratio_graph_data(&ratio_unfolded_data);
   TGraphErrors ratio_reg_graph(&ratio_unfolded_reg);
   TGraphErrors ratio_bbb_graph(&ratio_unfolded_bbb);
   if (plotComparison){
      for (int i=0; i<ratio_unfolded.GetNbinsX(); i++){    //change x error for plotting
         ratio_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded.GetBinError(i+1));
         ratio_graph_data.SetPointError(i,binning_met[1]/10.,ratio_unfolded_data.GetBinError(i+1));
         ratio_reg_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_reg.GetBinError(i+1));
         ratio_bbb_graph.SetPointError(i,binning_met[1]/10.,ratio_unfolded_bbb.GetBinError(i+1));
      }
   }
   ratio_graph.SetFillStyle(1001);
   ratio_graph_data.SetFillStyle(1001);
   ratio_reg_graph.SetFillStyle(1001);
   ratio_bbb_graph.SetFillStyle(1001);
   ratio_graph.SetLineWidth(0);
   ratio_graph_data.SetLineWidth(0);
   ratio_reg_graph.SetLineWidth(0);
   ratio_bbb_graph.SetLineWidth(0);
   ratio_graph.SetFillColor(kGray);
   ratio_graph_data.SetFillColor(kGray);
   ratio_reg_graph.SetFillColor(kGreen-9);
   ratio_bbb_graph.SetFillColor(kMagenta-9);
   ratio_graph.SetMarkerSize(0.4);
   ratio_graph_data.SetMarkerSize(0.4);
   ratio_reg_graph.SetMarkerSize(0.4);
   ratio_bbb_graph.SetMarkerSize(0.4);
      
   // setup axis
   if(dist.is2D){
      ratio.SetMaximum(1.49);
      ratio.SetMinimum(0.51);
      ratio_data_axis.SetMaximum(1.49);
      ratio_data_axis.SetMinimum(0.51);
      if (onlyTheo) {
         ratio.SetMaximum(1.55);
         ratio_data_axis.SetMaximum(1.55);
         ratio.SetMinimum(0.21);
         ratio_data_axis.SetMinimum(0.21);
      }
   }
   else if (!is_dPhi_1d) {
      ratio.SetMaximum((rewStudy)? 1.49 : 1.89);
      ratio.SetMinimum((rewStudy)? 0.51 : 0.11);
      ratio_data_axis.SetMaximum((rewStudy)? 1.49 : 1.89);
      ratio_data_axis.SetMinimum((rewStudy)? 0.51 : 0.11);
   }
   else {
      ratio.SetMaximum((rewStudy)? 1.49 : 1.25);
      ratio.SetMinimum((rewStudy)? 0.51 : 0.75);
      ratio_data_axis.SetMaximum((rewStudy)? 1.49 : 1.25);
      ratio_data_axis.SetMinimum((rewStudy)? 0.51 : 0.75);
   }
   ratio.SetLineColor(kRed-6);
   ratio.SetMarkerColor(kRed-6);
   ratio.GetYaxis()->SetTitleOffset(0.4);
   ratio.GetXaxis()->SetTitleOffset(1.7);
   ratio.GetXaxis()->SetLabelOffset(0.015);
   ratio.GetXaxis()->SetTickLength(0.);
   ratio_data.SetLineColor(kRed-6);
   ratio_data.SetMarkerColor(kRed-6);
   ratio_data_axis.GetYaxis()->SetTitleOffset(0.4);
   ratio_data_axis.GetXaxis()->SetTitleOffset(1.7);
   ratio_data_axis.GetXaxis()->SetLabelOffset(0.015);
   ratio_data_axis.GetXaxis()->SetTickLength(0.);
   
   ratio_NNLO_graph.SetFillStyle(1001);
   ratio_NNLO_graph.SetFillColorAlpha(kRed, 0.35);
   ratio_NNLO_data_graph.SetFillStyle(1001);
   ratio_NNLO_data_graph.SetFillColorAlpha(kRed, 0.35);
   ratio_NLO_graph.SetFillStyle(1001);
   ratio_NLO_graph.SetFillColorAlpha(kBlue, 0.35);
   ratio_NLO_data_graph.SetFillStyle(1001);
   ratio_NLO_data_graph.SetFillColorAlpha(kBlue, 0.35);
   ratio_LO_graph.SetFillStyle(1001);
   ratio_LO_graph.SetFillColorAlpha(kGreen, 0.35);
   ratio_LO_data_graph.SetFillStyle(1001);
   ratio_LO_data_graph.SetFillColorAlpha(kGreen, 0.35);
   
   if (plotComparison || onlyTheo) {   //first ratio draw for axis
      ratio.Draw("hist");
   }
   else {
      ratio_data_axis.Draw("hist");
   }
   
   ratio_totalGraph.SetLineWidth(1.);
   ratio_totalGraph_data.SetLineWidth(0.);
   ratio_reg_totalGraph.SetLineWidth(1.);
   ratio_bbb_totalGraph.SetLineWidth(1.);
   ratio_totalGraph.SetMarkerSize(0.4);
   ratio_totalGraph_data.SetMarkerSize(0.4);
   ratio_reg_totalGraph.SetMarkerSize(0.4);
   ratio_bbb_totalGraph.SetMarkerSize(0.4);
   
   ratio_totalGraph_data.SetFillStyle(1001);
   ratio_totalGraph_data.SetFillColor(kOrange+1);
   
   gfx::LegendEntries legE_ratio;
   
   if (plotComparison || onlyTheo) {   // draw MC
      ratio_alt.SetLineColor(kBlue-6);
      ratio_alt.SetMarkerColor(kBlue-6);
      if(isAltMC || rewStudy || isBSM) ratio_alt.Draw("same");
   }
   else {
      ratio_data.SetLineColor(kGreen+2);
      ratio_data.SetMarkerColor(kGreen+2);
      ratio_data.SetMarkerSize(0.5);
      ratio_data.SetLineWidth(1);
      ratio_data.SetMarkerStyle(22);
      ratio_data.Draw("p same");
   }
   
   if (plotComparison && !onlyTheo) {     //draw exp results
      
      // first draw used for axis and stuff
      // ~ratio_unfolded.Draw("pex0 same");
      
      //draw stat. unc
      // ~ratio_graph.Draw("pe20 same");
      // ~ratio_reg_graph.Draw("pe20 same");
      // ~ratio_bbb_graph.Draw("pe20 same");
      
      //draw tot. unc.
      ratio_totalGraph.Draw("p0 same");
      ratio_reg_totalGraph.Draw("p0 same");
      ratio_bbb_totalGraph.Draw("p0 same");
      
   }
   else if(!onlyTheo){
      ratio_totalGraph_data.Draw("2 same");
      ratio_graph_data.Draw("2 same");
      
      legE_ratio.append(ratio_totalGraph_data,"#sigma_{tot.}","f");
      legE_ratio.append(ratio_graph_data,"#sigma_{stat.}","f");
   }
   
   if(!plotComparison){
      aline->SetLineWidth(1);
      aline->DrawLine(ratio_data_axis.GetXaxis()->GetXmin(),1.,ratio_data_axis.GetXaxis()->GetXmax(),1.);   //Line at one for ratio
   }
      
   if(plotTheo){  //draw theory
      if (plotComparison || onlyTheo) {
         ratio_NNLO.SetMarkerSize(0);
         ratio_NLO.SetMarkerSize(0);
         ratio_LO.SetMarkerSize(0);
         ratio_NNLO.SetLineColor(kRed);
         ratio_NLO.SetLineColor(kBlue);
         ratio_LO.SetLineColor(kGreen);
         ratio_LO_graph.Draw("same 2");
         ratio_LO.Draw("same hist");
         ratio_NLO_graph.Draw("same 2");
         ratio_NLO.Draw("same hist");
         ratio_NNLO_graph.Draw("same 2");
         ratio_NNLO.Draw("same hist");
         
         ratio.Draw("hist same");
         ratio_alt.Draw("same");
         
      }
      else {
         ratio_NNLO_data_graph.SetMarkerSize(0.5);
         ratio_NLO_data_graph.SetMarkerSize(0.5);
         ratio_data_graph_asym.SetMarkerSize(0.5);
         ratio_alt_data_graph_asym.SetMarkerSize(0.5);
         ratio_herwig_data_graph_asym.SetMarkerSize(0.5);
         ratio_NNLO_data_graph.SetMarkerColor(kRed);
         ratio_NLO_data_graph.SetMarkerColor(kBlue);
         ratio_data_graph_asym.SetMarkerColor(kGreen+2);
         ratio_alt_data_graph_asym.SetMarkerColor(kViolet);
         ratio_herwig_data_graph_asym.SetMarkerColor(kCyan+2);
         ratio_NNLO_data_graph.SetLineColor(kRed);
         ratio_NLO_data_graph.SetLineColor(kBlue);
         ratio_data_graph_asym.SetLineColor(kGreen+2);
         ratio_alt_data_graph_asym.SetLineColor(kViolet);
         ratio_herwig_data_graph_asym.SetLineColor(kCyan+2);
         ratio_NNLO_data_graph.SetMarkerStyle(34);
         ratio_NLO_data_graph.SetMarkerStyle(33);
         ratio_data_graph_asym.SetMarkerStyle(22);
         ratio_alt_data_graph_asym.SetMarkerStyle(23);
         ratio_herwig_data_graph_asym.SetMarkerStyle(43);
         ratio_NNLO_data_graph.SetLineWidth(1);
         ratio_NLO_data_graph.SetLineWidth(1);
         ratio_data_graph_asym.SetLineWidth(1);
         ratio_alt_data_graph_asym.SetLineWidth(1);
         ratio_herwig_data_graph_asym.SetLineWidth(1);
         
         ratio_NNLO_data_graph.Draw("pe0 same");
         ratio_NLO_data_graph.Draw("pe0 same");
         ratio_data_graph_asym.Draw("pe0 same");
         ratio_alt_data_graph_asym.Draw("pe0 same");
         ratio_herwig_data_graph_asym.Draw("pe0 same");
         
         // ~ratio_alt_data.SetLineColor(kViolet);
         // ~ratio_alt_data.SetMarkerColor(kViolet);
         // ~ratio_alt_data.Draw("p same");
      }
   }
   
   if (plotComparison) ratio.Draw("axis same");
   else ratio_data_axis.Draw("axis same");
   
   TLegend leg_ratio=legE_ratio.buildLegend(.16,.85,0.45,.95,2);
   if (plotComparison || onlyTheo) leg_ratio.Draw();
   
   //Draw vertical lines as bins border for ratio
   aline->SetLineWidth(1);
   // ~aline->SetLineStyle(1);
   // ~aline->SetLineColor(kWhite);
   aline->SetLineStyle(7);
   aline->SetLineColor(kBlack);
   for (int i=0; i<(num_bins); i++){
      aline->DrawLine(xbins[i+1],ratio.GetMinimum(),xbins[i+1],ratio.GetMaximum());
   }
   
   if (dist.is2D){
      aline->SetLineStyle(1);
      aline->SetLineWidth(2);
      aline->SetLineColor(kBlack);
      aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
      aline->DrawLine(400,ratio.GetMinimum(),400,ratio.GetMaximum());
      aline->DrawLine(800,ratio.GetMinimum(),800,ratio.GetMaximum());
   }
   
   gStyle->SetTickLength(0.,"x");
   gStyle->SetLabelSize(0.06,"x");
   
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
   saver->save(can,saveName,(!isRealData || onlyTheo),true);
   
   return xbins_vec;
}

TH1F tunfoldplotting::getCRenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm){
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
   
   // Derive contributions to cov matrix (contribution only derived once, therefore only when up variation is derived)
   if (up) {
      int N_source = 3;
      if (norm){
         TH1F nominalResult_norm = nominalResult;
         nominalResult_norm.Scale(1./nominalResult.Integral());
         hist_CR1.Add(&nominalResult);
         hist_CR1.Scale(1./hist_CR1.Integral());
         hist_CR1 = phys::getSystShift(nominalResult_norm,hist_CR1);
         hist_CR2.Add(&nominalResult);
         hist_CR2.Scale(1./hist_CR2.Integral());
         hist_CR2 = phys::getSystShift(nominalResult_norm,hist_CR2);
         hist_ERDON.Add(&nominalResult);
         hist_ERDON.Scale(1./hist_ERDON.Integral());
         hist_ERDON = phys::getSystShift(nominalResult_norm,hist_ERDON);
      }
      for (int i=1; i<=hist_CR1.GetNbinsX(); i++){
         for (int j=1; j<=hist_CR1.GetNbinsX(); j++){
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_CR1.GetBinContent(i)*hist_CR1.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_CR2.GetBinContent(i)*hist_CR2.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_ERDON.GetBinContent(i)*hist_ERDON.GetBinContent(j));
         }
      }
   }
                  
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getMESCALEenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm){
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

   // Derive contributions to cov matrix (contribution only derived once, therefore only when up variation is derived)
   if (up) {
      int N_source = 6;
      if (norm){
         TH1F nominalResult_norm = nominalResult;
         nominalResult_norm.Scale(1./nominalResult.Integral());
         hist_MESCALE_UP.Add(&nominalResult);
         hist_MESCALE_UP.Scale(1./hist_MESCALE_UP.Integral());
         hist_MESCALE_UP = phys::getSystShift(nominalResult_norm,hist_MESCALE_UP);
         hist_MESCALE_DOWN.Add(&nominalResult);
         hist_MESCALE_DOWN.Scale(1./hist_MESCALE_DOWN.Integral());
         hist_MESCALE_DOWN = phys::getSystShift(nominalResult_norm,hist_MESCALE_DOWN);
         hist_MERENSCALE_UP.Add(&nominalResult);
         hist_MERENSCALE_UP.Scale(1./hist_MERENSCALE_UP.Integral());
         hist_MERENSCALE_UP = phys::getSystShift(nominalResult_norm,hist_MERENSCALE_UP);
         hist_MERENSCALE_DOWN.Add(&nominalResult);
         hist_MERENSCALE_DOWN.Scale(1./hist_MERENSCALE_DOWN.Integral());
         hist_MERENSCALE_DOWN = phys::getSystShift(nominalResult_norm,hist_MERENSCALE_DOWN);
         hist_MEFACSCALE_UP.Add(&nominalResult);
         hist_MEFACSCALE_UP.Scale(1./hist_MEFACSCALE_UP.Integral());
         hist_MEFACSCALE_UP = phys::getSystShift(nominalResult_norm,hist_MEFACSCALE_UP);
         hist_MEFACSCALE_DOWN.Add(&nominalResult);
         hist_MEFACSCALE_DOWN.Scale(1./hist_MEFACSCALE_DOWN.Integral());
         hist_MEFACSCALE_DOWN = phys::getSystShift(nominalResult_norm,hist_MEFACSCALE_DOWN);
      }
      for (int i=1; i<=hist_MESCALE_UP.GetNbinsX(); i++){
         for (int j=1; j<=hist_MESCALE_UP.GetNbinsX(); j++){
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MESCALE_UP.GetBinContent(i)*hist_MESCALE_UP.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MESCALE_DOWN.GetBinContent(i)*hist_MESCALE_DOWN.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MERENSCALE_UP.GetBinContent(i)*hist_MERENSCALE_UP.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MEFACSCALE_DOWN.GetBinContent(i)*hist_MEFACSCALE_DOWN.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MEFACSCALE_UP.GetBinContent(i)*hist_MEFACSCALE_UP.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*hist_MEFACSCALE_DOWN.GetBinContent(i)*hist_MEFACSCALE_DOWN.GetBinContent(j));
         }
      }
   }
               
   if(up) return env_up;
   else return env_down;
}

TH1F tunfoldplotting::getMTOPuncCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm){
   TH1F hist_MTOP_UP = vec_systShifts[0]["MTOP175p5"];
   TH1F hist_MTOP_DOWN = vec_systShifts[0]["MTOP169p5"];
   
   TH1F zeroes = hist_MTOP_UP;
   zeroes.Reset();
   
   for (int i=1; i<vec_systShifts.size(); i++){
      hist_MTOP_UP.Add(&vec_systShifts[i]["MTOP175p5"]);
      hist_MTOP_DOWN.Add(&vec_systShifts[i]["MTOP169p5"]);
   }
   
   hist_MTOP_UP.Scale(1/3.);
   hist_MTOP_DOWN.Scale(1/3.);
   
   // Derive contributions to cov matrix (contribution only derived once, therefore only when up variation is derived)
   if (up) {
      int N_source = 2;
      TH1F temp_UP = hist_MTOP_UP;
      TH1F temp_DOWN = hist_MTOP_DOWN;
      if (norm){
         TH1F nominalResult_norm = nominalResult;
         nominalResult_norm.Scale(1./nominalResult.Integral());
         temp_UP.Add(&nominalResult);
         temp_UP.Scale(1./temp_UP.Integral());
         temp_UP = phys::getSystShift(nominalResult_norm,temp_UP);
         temp_DOWN.Add(&nominalResult);
         temp_DOWN.Scale(1./temp_DOWN.Integral());
         temp_DOWN = phys::getSystShift(nominalResult_norm,temp_DOWN);
      }
      for (int i=1; i<=hist_MTOP_UP.GetNbinsX(); i++){
         for (int j=1; j<=hist_MTOP_UP.GetNbinsX(); j++){
            cov_syst.Fill(i,j,1/(1.*N_source)*temp_UP.GetBinContent(i)*temp_UP.GetBinContent(j));
            cov_syst.Fill(i,j,1/(1.*N_source)*temp_DOWN.GetBinContent(i)*temp_DOWN.GetBinContent(j));
         }
      }
   }
               
   if(up) return hist_MTOP_UP;
   else return hist_MTOP_DOWN;
}

TH2D tunfoldplotting::getNormCov(TH2D const &input_cov, TH1F const &input_hist){
   
   int nBins = input_hist.GetNbinsX();
   float integral = input_hist.Integral();
   
   TMatrixD cov(nBins,nBins);
   TMatrixD jac(nBins,nBins);
   for (int i=1; i<=nBins;i++){
      for (int j=1; j<=nBins;j++){
         cov[i-1][j-1] = input_cov.GetBinContent(i,j);
         if (i==j) jac[i-1][j-1] = (integral-input_hist.GetBinContent(i))/(integral*integral);
         else jac[i-1][j-1] = (-1.*input_hist.GetBinContent(i))/(integral*integral);
      }
   }
      
   TMatrixD resultCov(cov,TMatrixD::kMultTranspose,jac);
   resultCov = jac*resultCov;
   
   TH2D resultCov_hist = input_cov;
   resultCov_hist.Reset();
   
   for (int i=1; i<=nBins;i++){
      for (int j=1; j<=nBins;j++){
         resultCov_hist.SetBinContent(i,j,resultCov[i-1][j-1]);
      }
   }
   
   return resultCov_hist;
}



std::map<TString,TH1F> tunfoldplotting::getCombinedUnc(std::vector<std::map<TString,TH1F>>& vec_systShifts, const std::vector<TString>& systVec, const TH1F& combinedResult,
                                                       std::vector<TH1F> nominalResults, TH2D& cov_syst, bool const &norm){
   std::map<TString,TH1F>map_combinedShifts;
   
   TH1F* hist_TotalShiftUP = (TH1F*)combinedResult.Clone();
   TH1F* hist_TotalShiftDOWN = (TH1F*)combinedResult.Clone();
   hist_TotalShiftUP->Reset();
   hist_TotalShiftDOWN->Reset();
   
   TH1F* zeroes = (TH1F*)combinedResult.Clone();
   zeroes->Reset();
   
   TH1F combinedResult_norm = hist::getNormalizedHist(combinedResult);
   
   std::vector<TH1F> nominalResults_norm;
   std::vector<TH1F> nominalResults_clone;
   for (int y=0; y<=3; y++){
      nominalResults_clone.push_back(nominalResults[y]);
      nominalResults_norm.push_back(nominalResults[y]);
      nominalResults_norm[y].Scale(1./nominalResults_norm[y].Integral());
   }
   
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
   bool useMTop = false;
         
   for (const TString& syst : systVec){
      if (syst == "JESUserDefinedHEM1516_DOWN") { // HEM only used for 2018
         map_combinedShifts[syst] = vec_systShifts[3][syst];
         continue;
      }
      
      if (syst.BeginsWith("LUMI")) { // Lumi unc correlation more complex then just (un)corr.
         if (!lumiDone){
            nominalResults[1].Add(&nominalResults[0]);   // add preVFP and post VFP
            
            std::vector<TH1F> nominalResults_norm_lumi;    //normalized with pre and post VFP added
            for (int y=0; y<=3; y++){
               nominalResults_norm_lumi.push_back(nominalResults[y]);
               nominalResults_norm_lumi[y].Scale(1./nominalResults_norm_lumi[y].Integral());
            }
            
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
            
            // Derive contributions to cov matrix (no separation between up and down due the cancellation in minus sign! effectively it is N_source=1
            // Uncorrelated part should be zero in case of normalized distributions!
            if (norm){
               tempHist_corr.Add(&combinedResult);
               tempHist_corr.Scale(1./tempHist_corr.Integral());
               tempHist_corr = phys::getSystShift(combinedResult_norm,tempHist_corr);
               tempHist_corr1718.Add(&combinedResult);
               tempHist_corr1718.Scale(1./tempHist_corr1718.Integral());
               tempHist_corr1718 = phys::getSystShift(combinedResult_norm,tempHist_corr1718);               
            }
            
            for (int i=1; i<=vec_systShifts[0][syst].GetNbinsX(); i++){
               for (int j=1; j<=vec_systShifts[0][syst].GetNbinsX(); j++){
                  for (int y=1; y<=3; y++){ //uncorrelated part
                     if (!norm) cov_syst.Fill(i,j,(nominalResults[y].GetBinContent(i)*lumiCorr[y-1][y-1]*1e-2*nominalResults[y].GetBinContent(j)*lumiCorr[y-1][y-1]*1e-2));
                  }
                  cov_syst.Fill(i,j,tempHist_corr.GetBinContent(i)*tempHist_corr.GetBinContent(j));      //correlated parts
                  cov_syst.Fill(i,j,tempHist_corr1718.GetBinContent(i)*tempHist_corr1718.GetBinContent(j));
               }
            }
         }
         continue;
      }
      
      // TTBAR norm. uncertainty for MC (5%)
      if (syst.BeginsWith("XSEC_TTSIGNAL")) {
         TH1F combinedResult_scaled = combinedResult;
         combinedResult_scaled.Scale(0.05);
         map_combinedShifts["XSEC_TTSIGNAL_UP"] = combinedResult_scaled;
         combinedResult_scaled.Scale(-1.);
         map_combinedShifts["XSEC_TTSIGNAL_DOWN"] = combinedResult_scaled;
         continue;
      }
      
      TH1F tempHist;
      
      // Get CR envelope
      if (Systematic::convertType(syst) == Systematic::CR_envelope) {
         tempHist = getCRenvelopeCombined(vec_systShifts,combinedResult,cov_syst,Systematic::convertVariation(syst) == Systematic::up,norm);
         map_combinedShifts[syst] = tempHist;
         useCRenvelope = true;
         continue;
      }
      else if ((std::find(Systematic::crTypes.begin(), Systematic::crTypes.end(), Systematic::convertType(syst)) != Systematic::crTypes.end()) && useCRenvelope) continue;  //ignore shifts already used in envelope
      
      // Get MEscale envelope
      if (Systematic::convertType(syst) == Systematic::meScale_envelope) {
         tempHist = getMESCALEenvelopeCombined(vec_systShifts,combinedResult,cov_syst,Systematic::convertVariation(syst) == Systematic::up,norm);
         map_combinedShifts[syst] = tempHist;
         useMEenvelope = true;
         continue;
      }
      else if ((std::find(Systematic::meTypes.begin(), Systematic::meTypes.end(), Systematic::convertType(syst)) != Systematic::meTypes.end()) && useMEenvelope) continue;  //ignore shifts already used in envelope
      
      // Get MTop uncertainty
      if (Systematic::convertType(syst) == Systematic::mtop) {
         tempHist = getMTOPuncCombined(vec_systShifts,combinedResult,cov_syst,Systematic::convertVariation(syst) == Systematic::up,norm);
         map_combinedShifts[syst] = tempHist;
         useMTop = true;
         continue;
      }
      else if ((std::find(Systematic::mTopTypes.begin(), Systematic::mTopTypes.end(), Systematic::convertType(syst)) != Systematic::mTopTypes.end()) && useMTop) continue;  //ignore shifts already used in envelope
      
      // Combine all other systematics
      if (Systematic::isCorrelated(syst)){   // Add correlated parts
         tempHist = vec_systShifts[0][syst];
         for (int i=1; i<vec_systShifts.size(); i++){
            tempHist.Add(&vec_systShifts[i][syst]);
         }
         map_combinedShifts[syst] = tempHist;
         
         // Derive contributions to cov matrix
         int N_source = 2;
         if (std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), Systematic::convertType(syst)) == Systematic::upDownTypes.end()){
            N_source = 1;
         }
         if (norm){
            tempHist.Add(&combinedResult);
            tempHist.Scale(1./tempHist.Integral());
            tempHist = phys::getSystShift(combinedResult_norm,tempHist);
         }
         for (int i=1; i<=tempHist.GetNbinsX(); i++){
            for (int j=1; j<=tempHist.GetNbinsX(); j++){
               cov_syst.Fill(i,j,1/(1.*N_source)*tempHist.GetBinContent(i)*tempHist.GetBinContent(j));
            }
         }
      }
      else {   // Add non-correlated parts
         if (std::find(Systematic::upDownTypes.begin(), Systematic::upDownTypes.end(), Systematic::convertType(syst)) == Systematic::upDownTypes.end()){
            std::cout<<"Error: Combination of non-correlated uncertainties between years not supported of uncertainty not up_down type: "<<syst<<std::endl;
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
         
         for (int i=0; i<vec_systShifts.size(); i++){    //derive envelope to get up and down shifts (removes correct correlations!!)
            std::pair<TH1F*,TH1F*> envelopes =  hist::getEnvelope(zeroes,{vec_systShifts[i][systName_down],vec_systShifts[i][systName_up]});
            if (norm){
               TH1F temp_down = vec_systShifts[i][systName_down];
               TH1F temp_up = vec_systShifts[i][systName_up];
               temp_down.Add(&combinedResult);
               temp_up.Add(&combinedResult);
               temp_down.Scale(1./temp_down.Integral());
               temp_up.Scale(1./temp_up.Integral());
               temp_up = phys::getSystShift(combinedResult_norm,temp_up);
               temp_down = phys::getSystShift(combinedResult_norm,temp_down);
               envelopes =  hist::getEnvelope(zeroes,{temp_down,temp_up});
            }
            hist::addQuadr(down,*envelopes.first);
            hist::addQuadr(up,*envelopes.second);
         }
         down.Scale(-1.);
         map_combinedShifts[systName_down] = down;
         map_combinedShifts[systName_up] = up;
                  
         // Derive contributions to cov matrix
         int N_source = 2;    //only up and down uncertainties in this part
         for (int y=0; y<=3; y++){  //Uncorrelated contributions have to be treated individually per year
            TH1F temp = vec_systShifts[y][syst];
            if (norm){
               temp.Add(&combinedResult);
               temp.Scale(1./temp.Integral());
               temp = phys::getSystShift(combinedResult_norm,temp);
            }
            for (int i=1; i<=vec_systShifts[0][syst].GetNbinsX(); i++){
               for (int j=1; j<=vec_systShifts[0][syst].GetNbinsX(); j++){
                  cov_syst.Fill(i,j,1/(1.*N_source)*temp.GetBinContent(i)*temp.GetBinContent(j));
               }
            }
         }
      }
   }
   
   //Apply correction for normalized distributions
   if (norm){      
      for (auto &[key, value]: map_combinedShifts){
         
         TString tempKey = key;  //temp key to get proper syst name
         tempKey.ReplaceAll("_UP","");
         tempKey.ReplaceAll("_DOWN","");
         if (!Systematic::isCorrelated(tempKey) && tempKey!="LUMI") continue;
         
         value.Add(&combinedResult);
         value.Scale(1./value.Integral());
                  
         map_combinedShifts[key] = phys::getSystShift(combinedResult_norm,value);
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
   if (norm){
      map_combinedShifts["STAT_UP"] = combinedResult_norm;
      map_combinedShifts["STAT_DOWN"] = combinedResult_norm;
   }
   else{
      map_combinedShifts["STAT_UP"] = combinedResult;
      map_combinedShifts["STAT_DOWN"] = combinedResult;
   }
   
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

std::pair<TH1F,std::pair<TH1F,TH1F>> tunfoldplotting::readFixedOrderPred(TString const &filePath,bool const &is2D,bool const &norm,bool const &isNNLO, bool const &isdPhi){
   
   ///////////////////////////////////////////////////////
   //Uncertainty for normalized result not yet correct!!//
   ///////////////////////////////////////////////////////
   std::ifstream inputFile(filePath.Data());
   if(!inputFile.is_open()){
      std::cout<<"File not found: "<<filePath.Data()<<std::endl;
      exit(400);
   }
   std::string line;
   std::vector<std::vector<float>> inputData((is2D&&isNNLO)? 10: ((is2D||isNNLO )? 8: 6));
   
   int numBins = 0;
   while (std::getline(inputFile, line)){    //loop over inputs and write to array
      std::istringstream iss(line);
      std::vector<std::string> stringSplit((std::istream_iterator<std::string>(iss)), std::istream_iterator<std::string>());   //split by space
      for(int i=0; i<stringSplit.size();i++){
         inputData[i].push_back(std::stof(stringSplit[i]));
      }
      numBins++;
   }
   inputFile.close();
   
   TH1F outHist;
   TH1F outShiftUp;
   TH1F outShiftDown;
   int binOffset = 1;
      
   if(is2D){
      TH1F temp("","",numBins,0.5,numBins+0.5);
      TH1F temp_UP(temp);
      TH1F temp_DOWN(temp);
      float currentBinWidth = 1.;
      for (int i=1; i<=temp.GetNbinsX(); i++){
         // ~if(i%7==0){    //handle missing overflow bins
            // ~binOffset++;
            // ~continue;
         // ~}
         currentBinWidth = (inputData[3][i-binOffset]-inputData[2][i-binOffset])*(inputData[1][i-binOffset]-inputData[0][i-binOffset]);
         temp.SetBinContent(i,inputData[4][i-binOffset]*currentBinWidth);
         if(isNNLO){
            temp_UP.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[5][i-binOffset])*currentBinWidth,2)+std::pow(abs(inputData[4][i-binOffset]-inputData[8][i-binOffset])*currentBinWidth,2)+std::pow(inputData[7][i-binOffset]*currentBinWidth,2)));
            temp_DOWN.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[6][i-binOffset])*currentBinWidth,2)+std::pow(abs(inputData[4][i-binOffset]-inputData[9][i-binOffset])*currentBinWidth,2)+std::pow(inputData[7][i-binOffset]*currentBinWidth,2)));
         }
         // ~else {
            // ~temp_UP.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[5][i-binOffset])*currentBinWidth,2)+std::pow(inputData[7][i-binOffset]*currentBinWidth,2)));
            // ~temp_DOWN.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[6][i-binOffset])*currentBinWidth,2)+std::pow(inputData[7][i-binOffset]*currentBinWidth,2)));
         // ~}
         else {
            temp_UP.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[5][i-binOffset])*currentBinWidth,2)));
            temp_DOWN.SetBinContent(i,std::sqrt(std::pow(abs(inputData[4][i-binOffset]-inputData[6][i-binOffset])*currentBinWidth,2)));
         }
         temp.SetBinError(i,0.5*(temp_UP.GetBinContent(i)+temp_DOWN.GetBinContent(i)));
      }
      
      if(norm){   //for normalized distributions
         temp_UP.Add(&temp);  //First reproduce shifted values (not only differences)
         temp_DOWN.Add(&temp,-1);
         temp_DOWN.Scale(-1);
         
         temp.Scale(1./temp.Integral());  //Scale by individal integral
         temp_UP.Scale(1./temp_UP.Integral());
         temp_DOWN.Scale(1./temp_DOWN.Integral());
         
         std::pair<TH1F*,TH1F*> shiftPair = hist::getEnvelope(&temp,{temp_UP,temp_DOWN});    //get envelope to correctly assign up and down shifts
         temp_DOWN = *shiftPair.first;
         temp_UP = *shiftPair.second;
         
         temp_UP.Add(&temp,-1);  //Take difference to normalized
         temp_DOWN.Add(&temp,-1);
         temp_DOWN.Scale(-1.);
                           
         for (int i=1; i<=temp.GetNbinsX(); i++){
            temp.SetBinError(i,0.5*(temp_UP.GetBinContent(i)+temp_DOWN.GetBinContent(i)));
         }
      }
      
      outHist = temp;
      outShiftUp = temp_UP;
      outShiftDown = temp_DOWN;
      
   }
   else{
      std::vector<float> edges(inputData[0]);
      if (!isdPhi) edges.push_back(500);   // required because last bin edge in theo pred. is taken as 13000 GeV
      else edges.push_back(3.2);
            
      TH1F temp("","",numBins,&edges[0]);
      TH1F temp_UP(temp);
      TH1F temp_DOWN(temp);
      for (int i=1; i<=temp.GetNbinsX(); i++){
         
         if(i==temp.GetNbinsX()){   //fix width of overflow bin to plotted range
            for(int j=2; j<inputData.size();j++){
               inputData[j][i-1] = inputData[j][i-1]*(inputData[1][i-1]-inputData[0][i-1])/(edges[i]-edges[i-1]);
            }
         }
         
         temp.SetBinContent(i,inputData[2][i-1]);
         if(isNNLO){
            temp_UP.SetBinContent(i,std::sqrt(std::pow(abs(inputData[2][i-1]-inputData[3][i-1]),2)+std::pow(abs(inputData[2][i-1]-inputData[6][i-1]),2)+std::pow(inputData[5][i-1],2)));
            temp_DOWN.SetBinContent(i,std::sqrt(std::pow(abs(inputData[2][i-1]-inputData[4][i-1]),2)+std::pow(abs(inputData[2][i-1]-inputData[7][i-1]),2)+std::pow(inputData[5][i-1],2)));
         }
         else {
            temp_UP.SetBinContent(i,std::sqrt(std::pow(abs(inputData[2][i-1]-inputData[3][i-1]),2)+std::pow(inputData[5][i-1],2)));
            temp_DOWN.SetBinContent(i,std::sqrt(std::pow(abs(inputData[2][i-1]-inputData[4][i-1]),2)+std::pow(inputData[5][i-1],2)));
         }
         temp.SetBinError(i,0.5*(temp_UP.GetBinContent(i)+temp_DOWN.GetBinContent(i)));
      }
      
      if(norm){   //for normalized distributions
         temp_UP.Add(&temp);  //First reproduce shifted values (not only differences)
         temp_DOWN.Add(&temp,-1);
         temp_DOWN.Scale(-1);
         
         temp.Scale(1./temp.Integral("width"));  //Scale by individal integral
         temp_UP.Scale(1./temp_UP.Integral("width"));
         temp_DOWN.Scale(1./temp_DOWN.Integral("width"));
         
         std::pair<TH1F*,TH1F*> shiftPair = hist::getEnvelope(&temp,{temp_UP,temp_DOWN});    //get envelope to correctly assign up and down shifts
         temp_DOWN = *shiftPair.first;
         temp_UP = *shiftPair.second;
         
         temp_UP.Add(&temp,-1);  //Take difference to normalized
         temp_DOWN.Add(&temp,-1);
         temp_DOWN.Scale(-1.);
                           
         for (int i=1; i<=temp.GetNbinsX(); i++){
            temp.SetBinError(i,0.5*(temp_UP.GetBinContent(i)+temp_DOWN.GetBinContent(i)));
         }
      }
      
      outHist = temp;
      outShiftUp = temp_UP;
      outShiftDown = temp_DOWN;
   }
   
   return std::make_pair(outHist,std::make_pair(outShiftDown,outShiftUp));
}
            
