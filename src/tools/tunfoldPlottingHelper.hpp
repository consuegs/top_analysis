#ifndef TUNFOLDPLOTTINGHELPER_HPP__
#define TUNFOLDPLOTTINGHELPER_HPP__

#include <iostream>
#include <iomanip>
#include <cmath>
#include <map>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TH1.h>
#include "TUnfoldDensity.h"
#include <TRandom3.h>
#include <TProfile.h>
#include <TVectorD.h>
#include <TParameter.h>
#include <THashList.h>

#include "tools/hist.hpp"
#include "tools/io.hpp"
#include "tools/physics.hpp"

namespace tunfoldplotting
{
   struct distrUnfold {
   TString varName;
   float xMin;
   float xMax;
   TString title;
   TString labelFormat;
   bool is2D;
   bool norm = false;
   };
   
   std::pair<float,int> getChi2NDF(TH1F* hist_res, TH1F* hist_true, const bool includeTrueError = false);
   std::pair<float,int> getChi2NDF(TGraphAsymmErrors* graph_res, TH1F* hist_true, const bool includeTrueError = false);
   std::pair<float,int> getChi2NDF(TGraph* graph_res, TH1F* hist_true, const bool includeTrueError = false);
   std::pair<float,int> getChi2NDF(TGraphAsymmErrors* graph_res, TGraphAsymmErrors* graph_true, const bool includeTrueError = false);
   std::pair<float,int> getChi2NDF_withCorr(TH1F* hist_res, TGraphAsymmErrors* graph_true, TH2D* corr_res, bool const norm=false, const bool includeTrueError = false);
   TH2F* get_response(TH2F* responseHist,bool columnNormalized = false,bool includeUnderflow = false);
   void plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver, const bool is2D);
   void plot_correlation(TH2F* corrMatrix, TString name, io::RootFileSaver* saver);
   void plot_systBreakdown(std::map<TString,TH1F> const &indShifts, io::RootFileSaver* saver, TString const &name, TString const &method, TString var, const bool & isData,
                        std::map<TString,std::vector<TString>> const &systCombinations={}, const bool jesComparison=false);
   void plot_response_diff(std::map<TString,TH2F> const &indResponse, TH2F* const &nomResponse, io::RootFileSaver* saver, TString const &name);
   void plot_pur_stab_eff(TH1F* purity, TH1F* stability, TH1F* efficiency, TString const &name, TString const &xTitle, const distrUnfold &dist, TUnfoldBinning* const generatorBinning,io::RootFileSaver* saver);
   TH1F getCRenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getMESCALEenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getPDFenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false, bool const &scaleLumi=true);
   TH1F getMTOPunc(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getLUMIunc(TH1F* const &nominal, bool const &up, bool const &norm=false);
   std::pair<TH1F*,TH1F*> getSystUnc(TH1F* const &nominal, TString const &path, TString const &filePath, std::vector<TString> const &systVec, 
                                     bool const &addStat, bool const &withScaleFactor, std::map<TString,TH1F> &indShifts, std::map<TString,TH2F> &indResponse, std::map<TString,TH2F> &indResponseAlt, bool const &norm=false);
                                    
   std::map<TString,TH1F> getCombinedUnc(std::vector<std::map<TString,TH1F>>& vec_systShifts,const std::vector<TString>& systVec, const TH1F& combinedResult, std::vector<TH1F> nominalResults, TH2D& cov_syst, bool const &norm=false);
   TH1F getCRenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm=false);
   TH1F getMESCALEenvelopeCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm=false);
   TH1F getMTOPuncCombined(std::vector<std::map<TString,TH1F>>& vec_systShifts, const TH1F& nominalResult, TH2D& cov_syst, bool const &up, bool const &norm=false);
   std::pair<TH1F*,TH1F*> getTotalShifts(const std::map<TString,TH1F> &map_combinedShifts, const TH1F &nominal, const bool isNorm, const float &scale = 1.);
   
   std::vector<double> plot_UnfoldedResult(TUnfoldBinning* generatorBinning, TH1F* unfolded, TH1F* unfolded_reg, TH1F* unfolded_bbb,
                            std::pair<TH1F*,TH1F*> &unfolded_total, std::pair<TH1F*,TH1F*> &unfolded_reg_total, std::pair<TH1F*,TH1F*> &unfolded_bbb_total, std::pair<TH1F*,TH1F*> &realDis_syst_total, std::pair<TH1F*,TH1F*> &realDisAlt_syst_total, std::pair<TH1F*,TH1F*> &realDisHerwig_syst_total, const float &tau_par, TH1F* realDis, TH1F* realDisAlt, TH1F* realDisHerwig, TH2D* cov,const distrUnfold &dist, const bool plotComparison, const TString &saveName,io::RootFileSaver* saver, int &num_bins, const bool rewStudy, const bool onlyTheo, const bool plotTheo, std::ofstream &chi2_file, const bool divideBinWidth=true, const bool adaptYaxis=true);
   
   std::pair<TH1F,std::pair<TH1F,TH1F>> readFixedOrderPred(TString const &filePath,bool const &is2D,bool const &norm,bool const &isNNLO,bool const &isdPhi=false);
   
   TH2D getNormCov(TH2D const &input_cov, TH1F const &input_hist);
} // namespace tunfoldplotting

#endif /* TUNFOLDPLOTTINGHELPER_HPP__ */
