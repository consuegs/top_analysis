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
#include <TFile.h>
#include <TH1.h>
#include "TUnfoldDensity.h"
#include <TRandom3.h>
#include <TProfile.h>
#include <TVectorD.h>
#include <TParameter.h>

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
   
   std::pair<float,int> getChi2NDF(TH1F* hist_res, TH1F* hist_true);
   std::pair<float,int> getChi2NDF_withCorr(TH1F* hist_res, TH1F* hist_true, TH2F* corr_res);
   TH2F* get_response(TH2F* responseHist,bool columnNormalized = false);
   void plot_response(TH2F* responseHist, TString name, io::RootFileSaver* saver, const bool is2D);
   void plot_correlation(TH2F* corrMatrix, TString name, io::RootFileSaver* saver);
   void plot_systBreakdown(std::map<TString,TH1F> const &indShifts, io::RootFileSaver* saver, TString const &name, TString const &method, TString var,
                        std::map<TString,std::vector<TString>> const &systCombinations={}, const bool jesComparison=false);
   void plot_response_diff(std::map<TString,TH2F> const &indResponse, TH2F* const &nomResponse, io::RootFileSaver* saver, TString const &name);
   TH1F getCRenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getMESCALEenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getPDFenvelope(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getMTOPunc(TString const &path, TString const &filePath, TH1F* const &nominal, bool const &up, bool const &norm=false);
   TH1F getLUMIunc(TH1F* const &nominal, bool const &up, bool const &norm=false);
   std::pair<TH1F*,TH1F*> getSystUnc(TH1F* const &nominal, TString const &path, TString const &filePath, std::vector<TString> const &systVec, 
                                     bool const &addStat, bool const &withScaleFactor, std::map<TString,TH1F> &indShifts, std::map<TString,TH2F> &indResponse, bool const &norm=false);
                                    
   std::map<TString,TH1F> getCombinedUnc(std::vector<std::map<TString,TH1F>>& vec_systShifts,const std::vector<TString>& systVec, const TH1F& combinedResult, std::vector<TH1F> nominalResults);
   std::pair<TH1F*,TH1F*> getTotalShifts(const std::map<TString,TH1F> &map_combinedShifts, const TH1F &nominal, const bool isNorm, const float &scale = 1.);
   
   std::vector<double> plot_UnfoldedResult(TUnfoldBinning* generatorBinning, TH1F* unfolded, TH1F* unfolded_reg, TH1F* unfolded_bbb,
                            std::pair<TH1F*,TH1F*> &unfolded_total, std::pair<TH1F*,TH1F*> &unfolded_reg_total, std::pair<TH1F*,TH1F*> &unfolded_bbb_total,
                            const float &tau_par, TH1F* realDis, TH1F* realDisAlt,const distrUnfold &dist, const bool plotComparison, const TString &saveName,
                            io::RootFileSaver* saver, int &num_bins, const bool withScaleFactor, const bool divideBinWidth=true, const bool adaptYaxis=true);
} // namespace tunfoldplotting

#endif /* TUNFOLDPLOTTINGHELPER_HPP__ */
