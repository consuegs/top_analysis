#include "hist.hpp"

#include "tools/io.hpp"
#include "tools/gfx.hpp"
#include "tools/util.hpp"

#include <TMath.h>
#include <iomanip> 

static Config const &cfg=Config::get();

/*******************************************************************************
 * class Histograms
 ******************************************************************************/
template <class HIST>
int hist::Histograms<HIST>::s_iInstances = 0;

template <class HIST>
hist::Histograms<HIST>::Histograms(std::vector<TString> const &samples,DatasetCollection const &datasets)
   : iInstance_(++s_iInstances)
   , datasets(datasets)
   , vsSamples_(samples)
{}

template <class HIST>
hist::Histograms<HIST>::Histograms(std::vector<std::string> const &samples,DatasetCollection const &datasets)
   : iInstance_(++s_iInstances)
   , datasets(datasets)
{
   vsSamples_.clear();
   for (auto const &s:samples){
      vsSamples_.push_back(s);
   }
}

template <class HIST>
void hist::Histograms<HIST>::addHist(TString const &varName,TString const &title, Int_t nbinsx, Double_t xlow, Double_t xup)
{
   TH1::SetDefaultSumw2();
   mmH_[varName]=std::map<TString,HIST>();
   for (TString const &s: vsSamples_)
      mmH_[varName][s]=HIST(varName+"_"+s+TString::Format("_%d",iInstance_),title,nbinsx,xlow,xup);
}

template <class HIST>
void hist::Histograms<HIST>::addHist(TString const &varName,TString const &title, std::vector<float> edges, std::vector<float> widths)
{
   TH1::SetDefaultSumw2();
   mmH_[varName]=std::map<TString,HIST>();
   for (TString const &s: vsSamples_)
      mmH_[varName][s]=fromWidths(varName+"_"+s+TString::Format("_%d",iInstance_),title,edges,widths);
}

template <class HIST>
void hist::Histograms<HIST>::addFilledHist(TString const &varName,TString const &s,TH1F const &filledHist)
{
   TH1::SetDefaultSumw2();
   if (mmH_.find(varName)==mmH_.end()) mmH_[varName]=std::map<TString,HIST>();
   mmH_[varName][s]=filledHist;
}

template <class HIST>
void hist::Histograms<HIST>::addFilledHist(TString const &varName,TString const &s,TH2F const &filledHist)
{
   TH1::SetDefaultSumw2();
   if (mmH_.find(varName)==mmH_.end()) mmH_[varName]=std::map<TString,HIST>();
   mmH_[varName][s]=filledHist;
}

template <class HIST>
void hist::Histograms<HIST>::addHist(TString const &varName,TString const &title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
{
   TH1::SetDefaultSumw2();
   mmH_[varName]=std::map<TString,HIST>();
   for (TString const &s: vsSamples_)
      mmH_[varName][s]=HIST(varName+"_"+s+TString::Format("_%d",iInstance_),title,nbinsx,xlow,xup,nbinsy,ylow,yup);
}

template <class HIST>
void hist::Histograms<HIST>::addHist(TString const &varName,TString const &title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y)
{
   TH1::SetDefaultSumw2();
   mmH_[varName]=std::map<TString,HIST>();
   for (TString const &s: vsSamples_)
      mmH_[varName][s]=fromWidths_2d(varName+"_"+s+TString::Format("_%d",iInstance_),title,edges_x,widths_x,edges_y,widths_y);
}

template <class HIST>
void hist::Histograms<HIST>::addCounter(TString const &varName)
{
   mCount_[varName]=std::map<TString,float>();
   for (TString const &s: vsSamples_)
      mCount_[varName][s]=0.0;
}

template <class HIST>
void hist::Histograms<HIST>::setCurrentSample(TString const &current)
{
   sCurrentSample_=current;
}

template <class HIST>
void hist::Histograms<HIST>::setFillWeight(float w)
{
   fWeight_=w;
}

template <class HIST>
void hist::Histograms<HIST>::fill(TString const &varName,float x)
{
   if (mmH_.count(varName)<1) {debug_io*varName>>"unkown"; throw;}
   mmH_[varName][sCurrentSample_].Fill(x,fWeight_);
}

template <class HIST>
void hist::Histograms<HIST>::fillbin(TString const &varName,TString const &binName)
{
   if (mmH_.count(varName)<1) {debug_io*varName>>"unkown"; throw;}
   mmH_[varName][sCurrentSample_].Fill(binName,fWeight_);
}

template <class HIST>
void hist::Histograms<HIST>::fillbinFake(TString const &varName,TString const &binName)
{
   if (mmH_.count(varName)<1) {debug_io*varName>>"unkown"; throw;}
   mmH_[varName][sCurrentSample_].Fill(binName,0);
}

template <class HIST>
void hist::Histograms<HIST>::fillweight(TString const &varName,float x,float w)
{
   if (mmH_.count(varName)<1) {debug_io*varName>>"unkown"; throw;}
   mmH_[varName][sCurrentSample_].Fill(x,w*fWeight_);
}

template <class HIST>
void hist::Histograms<HIST>::fill(TString const &varName,float x,float y)
{
   if (mmH_.count(varName)<1) {debug_io*varName>>"unkown"; throw;}
   mmH_[varName][sCurrentSample_].Fill(x,y,fWeight_);
}

template <class HIST>
void hist::Histograms<HIST>::count(TString const &varName)
{
   mCount_[varName][sCurrentSample_]+=fWeight_;
   mCountError2_[varName][sCurrentSample_]+=fWeight_*fWeight_;
}

template <class HIST>
void hist::Histograms<HIST>::scaleLumi(const Systematic::Systematic& systematic)
{
   Datasubset curS=datasets.getDatasubset(sCurrentSample_);
   if (curS.isData) return; // don't scale data
   // Lumi weight
   float w=curS.xsec/curS.getNgen_syst(systematic)*cfg.lumi;
   for (auto &mH:mmH_){
      mH.second[sCurrentSample_].Scale(w);
   }
   for (auto &mC:mCount_){
      mC.second[sCurrentSample_]*=w;
   }
   for (auto &mCe:mCountError2_){
      mCe.second[sCurrentSample_]*=(w*w);
   }
}

template <class HIST>
void hist::Histograms<HIST>::mergeOverflow(bool includeUnderflow)
{
   // ~Datasubset curS=datasets.getDatasubset(sCurrentSample_);
   for (auto &mH:mmH_){
      hist::mergeOverflow(mH.second[sCurrentSample_],includeUnderflow);
   }
}

template <class HIST>
void hist::Histograms<HIST>::normHists()
{
   for (auto &mH:mmH_){
      mH.second[sCurrentSample_].Scale(1/mH.second[sCurrentSample_].Integral());
      mH.second[sCurrentSample_].GetYaxis()->SetTitle("Normalized Distributions");
   }
}

/*
 * Build histograms for whole processes from the corresponding
 * subsamples (e.g. QCD HT-bins). Counters are also combined.
 * Careful: don't combine samples with only one subsample having the same name,
 * it will simply be deleted (should work for histograms now, counters maybe not)
 */
template <class HIST>
void hist::Histograms<HIST>::combineFromSubsamples(std::vector<TString> const &samples)
{
   for (auto const &mH: mmH_){
      for (auto const &s: samples){
         auto const &subsets=datasets.getDataset(s).subsets;
         if (subsets.size()==1 && subsets[0].name==s) continue;
         mmH_[mH.first][s]=*(HIST*)mH.second.begin()->second.Clone();
         HIST &newH=mmH_[mH.first][s];
         newH.Reset();
         for (Datasubset const &dss: subsets){
            newH.Add(&mH.second.at(dss.name));
         }
      }
   }
   for (auto const &mC: mCount_){
      for (auto const &s: samples){
         mCount_[mC.first][s]=0.0;
         mCountError2_[mC.first][s]=0.0;
         for (Datasubset const &dss: datasets.getDataset(s).subsets){
            mCount_[mC.first][s]+=mC.second.at(dss.name);
            mCountError2_[mC.first][s]+=mCountError2_[mC.first][dss.name];
         }
      }
   }
}

template <class HIST>
void hist::Histograms<HIST>::combineSamples(TString const &sampleCombined, std::vector<TString> const &samples)
{
   for (auto const &mH: mmH_){
      mmH_[mH.first][sampleCombined]=*(HIST*)mH.second.begin()->second.Clone();
      HIST &newH=mmH_[mH.first][sampleCombined];
      newH.Reset();
      for (auto const &s: samples){
         newH.Add(&mH.second.at(s));
      }
   }
}

template <class HIST>
void hist::Histograms<HIST>::combineChannel(TString const &combinedName, std::vector<TString> const &channels)
{
   for (auto const &mH: mmH_){
      if (!mH.first.Contains(channels.back())) continue;
      TString combinedPath(mH.first);
      combinedPath.ReplaceAll(channels.back(),combinedName);
      if (mmH_.find(combinedPath) != mmH_.end()) continue;
      for (auto const &s: mH.second){
         mmH_[combinedPath][s.first]=*(HIST*)mH.second.begin()->second.Clone();
         HIST &newH=mmH_[combinedPath][s.first];
         newH.Reset();
         for (auto const &channel: channels){
            TString channelPath(mH.first);
            channelPath.ReplaceAll(channels.back(),channel);
            newH.Add(&mmH_[channelPath][s.first]);
         }
      }
   }
}

template <class HIST>
std::vector<TString> hist::Histograms<HIST>::getVariableNames()
{
   std::vector<TString> vars;
   for (auto const &mv:mmH_) vars.push_back(mv.first);
   return vars;
}

template <class HIST>
std::vector<HIST*> hist::Histograms<HIST>::getHistograms(TString const &varName,std::vector<TString> const &samples,bool divideByBinWidth)
{
   le_.clear();
   std::vector<HIST*> v;
   if (samples.size()>1) Color::reset();
   for (TString const &s: samples){
      HIST *h=(HIST*)mmH_.at(varName).at(s).Clone();
      if (divideByBinWidth) hist::divideByBinWidth(*h);
      v.push_back(h);
      h->SetLineWidth(2);
      try {
         Dataset const &ds=datasets.getDataset(s);
         h->SetLineColor(ds.color);
         if (ds.isData) h->SetLineWidth(1);
         le_.append(*h,datasets.getDataset(s).label,ds.isData?"pe":"l");
      } catch (const std::out_of_range& exc) {
         // no full dataset, use default colors
         h->SetLineColor(Color::next());
         le_.append(*h,s,"l");
      }
   }
   return v;
}

template <class HIST>
HIST* hist::Histograms<HIST>::getHistogram(TString const &varName,TString const &sample,bool divideByBinWidth)
{
   return getHistograms(varName,{sample},divideByBinWidth)[0];
}


template <class HIST>
THStack hist::Histograms<HIST>::getStack(TString const &varName,std::vector<TString> const& samples,std::map<const TString,Color_t> const& colormap,bool divideByBinWidth, bool includeData)
{
   le_.clear();
   THStack st;
   TAxis xax(*mmH_[varName].begin()->second.GetXaxis());
   TAxis yax(*mmH_[varName].begin()->second.GetYaxis());
   if (divideByBinWidth) {
      TString yt = yax.GetTitle();
      yt.ReplaceAll("BIN","BINW");
      yt.ReplaceAll(" / bin","BINW");
      yax.SetTitle(yt);
   }
   gfx::setupAxes(xax,yax);
   st.SetTitle(TString::Format(";%s;%s",xax.GetTitle(),yax.GetTitle()));
   Color::reset();
   for (TString const &s: samples){
      HIST *h=(HIST*)mmH_[varName][s].Clone();
      if (divideByBinWidth) hist::divideByBinWidth(*h);
      h->SetLineColor(kBlack);
      h->SetLineWidth(1);
      try {
         if (!includeData && datasets.getDataset(s).isData) continue;
         h->SetFillColor(datasets.getDataset(s).color);
         le_.prepend(*h,datasets.getDataset(s).label,"f");
      } catch (const std::out_of_range& exc) {
         // no full dataset, use default colors or color defined in colormap
         if (colormap.find(s)!=colormap.end()) h->SetFillColor(colormap.at(s));
         else h->SetFillColor(Color::next());
         le_.prepend(*h,s,"f");
      }
      h->SetFillStyle(1001);
      st.Add(h,"hist");
   }
   st.SetMinimum(1);
   return st;
}

template <class HIST>
HIST* hist::Histograms<HIST>::getSummedHist(TString const &varName,bool divideByBinWidth)
{
   HIST *h=(HIST*)mmH_[varName][vsSamples_[0]].Clone();
   h->Reset();
   for (TString const &s: vsSamples_){
      h->Add(&mmH_[varName][s]);
      if (divideByBinWidth) hist::divideByBinWidth(*h);
   }
   return h;
}

template <class HIST>
HIST* hist::Histograms<HIST>::getSummedHist(TString const &varName,std::vector<TString> const &samples,bool divideByBinWidth)
{
   HIST *h=(HIST*)mmH_[varName][samples[0]].Clone();
   h->Reset();
   for (TString const &s: samples){
      h->Add(&mmH_[varName][s]);
      if (divideByBinWidth) hist::divideByBinWidth(*h);
   }
   return h;
}

template <class HIST>
float hist::Histograms<HIST>::getCount(TString const &varName, TString const &sample)
{
   return mCount_[varName][sample];
}

template <class HIST>
float hist::Histograms<HIST>::getCountError(TString const &varName, TString const &sample)
{
   return TMath::Sqrt(mCountError2_[varName][sample]);
}

/* return the list of legend entries created by the last method call */
template <class HIST>
gfx::LegendEntries hist::Histograms<HIST>::getLegendEntries()
{
   return le_;
}

template <class HIST>
void hist::Histograms<HIST>::saveHistograms(io::RootFileSaver const &saver_hist, std::vector<TString> const &Samples)
{   
   for (auto const &mv:mmH_){
      for (TString sSample: Samples){
         auto temp=getHistogram(mv.first,sSample);
         temp->SetName(sSample);
         // ~saver_hist.save(*getHistogram(mv.first,sSample),mv.first+"/"+sSample);
         saver_hist.save(*temp,mv.first+"/"+sSample);
      }
   }
}

template <class HIST>
void hist::Histograms<HIST>::saveHistograms(io::RootFileSaver const &saver_hist, std::vector<std::string> const &Samples)
{   
   for (auto const &mv:mmH_){
      for (std::string sSample: Samples){
         TString tempString(sSample);
         auto temp=getHistogram(mv.first,tempString);
         temp->SetName(tempString);
         // ~saver_hist.save(*getHistogram(mv.first,sSample),mv.first+"/"+sSample);
         saver_hist.save(*temp,mv.first+"/"+tempString);
      }
   }
}

template <class HIST>
void hist::Histograms<HIST>::saveHistograms2D_as1D(io::RootFileSaver const &saver_hist, std::vector<TString> const &Samples)
{   
   for (auto const &mv:mmH_){
      for (TString sSample: Samples){
         auto temp=(TH2F*)getHistogram(mv.first,sSample);
         temp->SetName(sSample);
         saver_hist.save(histTrafo_2D(temp),mv.first+"/"+sSample);
      }
   }
}

/*******************************************************************************
 * end class Histograms
 ******************************************************************************/

std::vector<double> hist::getBinVector(std::vector<float> edges, std::vector<float> widths)
{
   assert(edges.size()==widths.size()+1);
   std::vector<double> xbins={edges[0]};
   for (uint i=0; i<edges.size()-1; i++){
      float x=edges[i];
      while ((edges[i+1]-x)>1e-5){
         x+=widths[i];
         xbins.push_back(x);
      }
   }
   return xbins;
}

std::vector<double> hist::getBinVector(const double Xmin, const double Xmax, const int nBins)
{
   
   std::vector<double> xbins={Xmin};
   float binWidth = (Xmax-Xmin)/(1.0*nBins);
   for (uint i=1; i<=nBins; i++){
      xbins.push_back(xbins[i-1]+binWidth);
   }
   return xbins;
}

TH1F hist::fromWidths(const char *name, const char *title,std::vector<float> edges, std::vector<float> widths)
{
   std::vector<double> xbins=getBinVector(edges, widths);
   return TH1F(name,title,xbins.size()-1,&xbins[0]);
}

TH2F hist::fromWidths_2d(const char *name, const char *title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y)
{
   std::vector<double> xbins=getBinVector(edges_x, widths_x);
   std::vector<double> ybins=getBinVector(edges_y, widths_y);
   return TH2F(name,title,xbins.size()-1,&xbins[0],ybins.size()-1,&ybins[0]);
}

TProfile2D hist::ProfilefromWidths_2d(const char *name, const char *title, std::vector<float> edges_x, std::vector<float> widths_x, std::vector<float> edges_y, std::vector<float> widths_y)
{
   std::vector<double> xbins=getBinVector(edges_x, widths_x);
   std::vector<double> ybins=getBinVector(edges_y, widths_y);
   return TProfile2D(name,title,xbins.size()-1,&xbins[0],ybins.size()-1,&ybins[0]);
}

std::vector<float> hist::getWidths(std::vector<float> const &bins){
   std::vector<float> widths;
   int nBins=bins.size();
   for(int i=0;i<(nBins-1);i++){
      widths.push_back(bins[i+1]-bins[i]);
   }
   return widths;
}

bool checkRebinningConistency(const TAxis* axis, std::vector<double> const &newBinning)
{
   int nBins_old = axis->GetNbins();
   std::vector<double> oldBinning;
   for(int binId = 0; binId <= nBins_old; ++binId) {
     oldBinning.push_back(round( axis->GetBinLowEdge(binId+1) * 1000.0 ) / 1000.0);
   }
   
   for(auto binEdge : newBinning){
      binEdge = round( binEdge * 1000.0 ) / 1000.0;
      if (std::find(oldBinning.begin(), oldBinning.end(), binEdge) == oldBinning.end()){
         return false;
      }
   }
   return true;
}

TH1F hist::rebinned(TH1F const &h, std::vector<float> const &edges, std::vector<float> const &widths,bool mergeOverflow,bool mergeUnderflow)
{
   std::vector<double> binedges=getBinVector(edges, widths);
   return rebinned(h,binedges,mergeOverflow,mergeUnderflow);
}

TH1F hist::rebinned(TH1F const &h, float const &Xmin, float const &Xmax, int const &nBins,bool mergeOverflow,bool mergeUnderflow)
{
   std::vector<double> binedges=getBinVector(Xmin, Xmax, nBins);
   return rebinned(h,binedges,mergeOverflow,mergeUnderflow);
}

TH1F hist::rebinned(TH1F const &h, std::vector<double> const &binedges,bool mergeOverflow,bool mergeUnderflow)
{
   TH1F hClone(h);
   std::string name(hClone.GetName());
   name+="_rebinned";
   if (checkRebinningConistency(h.GetXaxis(),binedges)==false){
      std::cout<<"Warning: Binning used for rebinning is not compatible:"<<hClone.GetXaxis()->GetTitle()<<std::endl;
   }
   TH1F *hnew=(TH1F*)hClone.Rebin(binedges.size()-1,name.c_str(),&binedges[0]);
   if (mergeOverflow) hist::mergeOverflow(*hnew,mergeUnderflow);
   TString yTitle=hClone.GetYaxis()->GetTitle();
   // ~yTitle.ReplaceAll("BIN"," / bin");
   hnew->GetYaxis()->SetTitle(yTitle);
   return *hnew;
   // ~return h;
}

TH2F hist::rebinned(TH2F const &h, float const &Xmin, float const &Xmax, int const &nBinsX, float const &Ymin, float const &Ymax, int const &nBinsY, bool mergeOverflow,bool mergeUnderflow)
{
   std::vector<double> binedgesX_d = getBinVector(Xmin, Xmax, nBinsX);
   std::vector<double> binedgesY_d = getBinVector(Ymin, Ymax, nBinsY);
   
   std::vector<float> binedgesX(binedgesX_d.begin(),binedgesX_d.end());
   std::vector<float> binedgesY(binedgesY_d.begin(),binedgesY_d.end());
   
   return rebinned(h,binedgesX,binedgesY,mergeOverflow,mergeUnderflow);
}

TH2F hist::rebinned(TH2F const &h, std::vector<float> const &binedges_x, std::vector<float> const &binedges_y,bool mergeOverflow,bool mergeUnderflow)
{
   TH2F hClone(h);
   int entries = h.GetEntries();
   std::string name(hClone.GetName());
   name+="_rebinned";
   TH2F hnew=(TH2F)fromWidths_2d("","",binedges_x,getWidths(binedges_x),binedges_y,getWidths(binedges_y));
   float_t errors[hnew.GetNbinsX()+2][hnew.GetNbinsY()+2] = {};
   TAxis *xaxis = hClone.GetXaxis();
   TAxis *yaxis = hClone.GetYaxis();
   TAxis *xaxis_new = hnew.GetXaxis();
   TAxis *yaxis_new = hnew.GetYaxis();
   for (int j=1; j<=yaxis->GetNbins();j++) {
      for (int i=1; i<=xaxis->GetNbins();i++) {
         hnew.Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j),h.GetBinContent(i,j));
         errors[xaxis_new->FindBin(xaxis->GetBinCenter(i))][yaxis_new->FindBin(yaxis->GetBinCenter(j))] += h.GetBinError(i,j)*h.GetBinError(i,j);
      } 
   }
   for (int i=0; i<hnew.GetNbinsX()+2; i++) {
      for (int j=0; j<hnew.GetNbinsY()+2; j++) {
         hnew.SetBinError(i,j,sqrt(errors[i][j]));
      }
   }
   
   if (mergeOverflow) hist::mergeOverflow(hnew,mergeUnderflow);
   TString yTitle=hClone.GetYaxis()->GetTitle();
   TString xTitle=hClone.GetXaxis()->GetTitle();
   yTitle.ReplaceAll("BIN"," / bin");
   hnew.GetYaxis()->SetTitle(yTitle);
   hnew.GetXaxis()->SetTitle(xTitle);
   hnew.SetEntries(entries);
   return hnew;
}

TH1F hist::histTrafo_2D(TH2F* const &hist2D){      
   
   int numBins_x = hist2D->GetNbinsX();
   int numBins_y = hist2D->GetNbinsY();
   int numBins = numBins_x*numBins_y;
   const TArrayD* binedges_x = hist2D->GetXaxis()->GetXbins();
   double binedges_1d[numBins+1];
   int phi_bin = 0;
   if (binedges_x->GetSize()==0){   //For Histograms with constant bin width (only use bin numbers)
      binedges_1d[0]=0.5;
      for (int i=1; i<=numBins; i++){
         binedges_1d[i]=i+0.5;
      }
   }
   else{
      binedges_1d[0]=hist2D->GetXaxis()->GetXmin();
      for (int i=0; i<(numBins); i++){
         binedges_1d[i+1] = binedges_x->GetAt(i%numBins_x+1)+phi_bin*binedges_x->GetAt(numBins_x);
         if (i%numBins_x==numBins_x-1) phi_bin++;
      }
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


void hist::divideByBinWidth(TH1& h,bool divideLastBin)
{
   int N=h.GetNbinsX();
   if (h.GetBinContent(N+1) != 0) {
      debug_io<<"non-emtpy overflow. merge first!";
      throw;
   }
   float w;
   if (divideLastBin) N++;
   //else: not dividing the last bin=merged overflow
   for (int i=0; i<N; i++) {
      w=h.GetBinWidth(i);
      h.SetBinContent(i,h.GetBinContent(i)/w);
      h.SetBinError(i,h.GetBinError(i)/w);
   }
   // labels
   TString yt=h.GetYaxis()->GetTitle();
   yt.ReplaceAll("BIN","BINW");
   yt.ReplaceAll(" / bin","BINW");
   h.GetYaxis()->SetTitle(yt);
}

void hist::mergeOverflow(TH1& h, bool includeUnderflow)
{
   int N=h.GetNbinsX();
   int entries=h.GetEntries();
   // -- overflow
   float cont=h.GetBinContent(N)+h.GetBinContent(N+1);
   float err2=util::quadSum<double>({h.GetBinError(N),h.GetBinError(N+1)});
   // set content+error of last bin
   h.SetBinContent(N,cont);
   h.SetBinError(N,TMath::Sqrt(err2));
   // clear overflow
   h.SetBinContent(N+1,0);
   h.SetBinError(N+1,0);
   if (!includeUnderflow) return;
   // -- underflow
   cont=h.GetBinContent(0)+h.GetBinContent(1);
   err2=util::quadSum<double>({h.GetBinError(0),h.GetBinError(1)});
   // set content+error of first bin
   h.SetBinContent(1,cont);
   h.SetBinError(1,TMath::Sqrt(err2));
   // clear overflow
   h.SetBinContent(0,0);
   h.SetBinError(0,0);
   // restore correct number of entries
   h.SetEntries(entries);
}

void hist::mergeOverflow(TH2& h, bool includeUnderflow)
{
   int N_X=h.GetNbinsX();
   int N_Y=h.GetNbinsY();
   int entries=h.GetEntries();
   includeUnderflow = true;
   //loop over Y-bins
   for (int i=0; i<=N_Y+1; i++){
      // -- overflow
      float cont=h.GetBinContent(N_X,i)+h.GetBinContent(N_X+1,i);
      float err2=util::quadSum<double>({h.GetBinError(N_X,i),h.GetBinError(N_X+1,i)});
      // set content+error of last bin
      h.SetBinContent(N_X,i,cont);
      h.SetBinError(N_X,i,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(N_X+1,i,0);
      h.SetBinError(N_X+1,i,0);
      if (!includeUnderflow) return;
      // -- underflow
      cont=h.GetBinContent(0,i)+h.GetBinContent(1,i);
      err2=util::quadSum<double>({h.GetBinError(0,i),h.GetBinError(1,i)});
      // set content+error of first bin
      h.SetBinContent(1,i,cont);
      h.SetBinError(1,i,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(0,i,0);
      h.SetBinError(0,i,0);
   }
   //Loop over X-bins
   for (int i=0; i<=N_X+1; i++){
      // -- overflow
      float cont=h.GetBinContent(i,N_Y)+h.GetBinContent(i,N_Y+1);
      float err2=util::quadSum<double>({h.GetBinError(i,N_Y),h.GetBinError(i,N_Y+1)});
      // set content+error of last bin
      h.SetBinContent(i,N_Y,cont);
      h.SetBinError(i,N_Y,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(i,N_Y+1,0);
      h.SetBinError(i,N_Y+1,0);
      if (!includeUnderflow) return;
      // -- underflow
      cont=h.GetBinContent(i,0)+h.GetBinContent(i,1);
      err2=util::quadSum<double>({h.GetBinError(i,0),h.GetBinError(i,1)});
      // set content+error of first bin
      h.SetBinContent(i,1,cont);
      h.SetBinError(i,1,TMath::Sqrt(err2));
      // clear overflow
      h.SetBinContent(i,0,0);
      h.SetBinError(i,0,0);
   }
   // restore correct number of entries
   h.SetEntries(entries);
}


void hist::setMaximum(TH1& h,std::vector<TH1F> hists,float multiplier)
{
   double max=h.GetMaximum();
   for (TH1 const &hh: hists) max=TMath::Max(max,hh.GetMaximum());
   h.SetMaximum(max*multiplier);
}

void hist::setMinimum(TH1& h,std::vector<TH1F> hists,float multiplier, bool allowNegative)
{
   double min=h.GetMaximum();
   hists.push_back(TH1F(*(TH1F*)&h));
   double v;
   for (TH1 const &hh: hists) {
      for (int i=0; i<=hh.GetNbinsX()+1; i++){
         v=hh.GetBinContent(i);
         if (allowNegative || v>0) min=TMath::Min(min,v);
      }
   }
   h.SetMinimum(min*multiplier);
}

void hist::sqrtHist(TH1& h)
{
   for (int bin=0; bin<=h.GetNcells(); ++bin) {
      float content = h.GetBinContent(bin);
      if (content >= 0){
         h.SetBinContent(bin, sqrt(content));
         if (content != 0){
            h.SetBinError(bin, 1./sqrt(content)*h.GetBinError(bin));
         }
      }
      else {
         std::cout<<"Error in sqrtHist: Bin "<<bin<<" has a negative content"<<std::endl;
      }
   }
}

void hist::addQuadr(TH1F &h1, TH1F const &h2)
{
   TH1F add(h2);
   h1.Multiply(&h1);
   add.Multiply(&add);
   h1.Add(&add);
   sqrtHist(h1);
}

TH1F hist::getRatio(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F hRatio(h1);
   hRatio.Divide(&h2);
   float N2,err;
   if (et==ONLY1){
      for (int i=0; i<=h1.GetNbinsX()+1;i++){
         N2 = h2.GetBinContent(i);
         err = (N2==0.0) ? 0.0 : h1.GetBinError(i)/N2;
         hRatio.SetBinError(i,err);
      }
   } else if (et==ONLY2){
      for (int i=0; i<=h1.GetNbinsX()+1;i++){
         N2 = h2.GetBinContent(i);
         err = (N2==0.0) ? 0.0 : h2.GetBinError(i)/N2;
         hRatio.SetBinError(i,err);
      }
   } else if (et==NOERR){
      for (int i=0; i<=h1.GetNbinsX()+1;i++){
         hRatio.SetBinError(i,0.);
      }
   } else assert(et==COMB);

   Double_t min = hRatio.GetBinContent(hRatio.GetMinimumBin());
   Double_t max = hRatio.GetBinContent(hRatio.GetMaximumBin());
   min=TMath::Min(min,0.5);
   max=TMath::Max(max,1.5);
   hRatio.GetYaxis()->SetRangeUser(min,max*1.1);
   hRatio.GetYaxis()->SetTitle(title);

   if (et==ONLY2 || et==COMB){
      hRatio.SetLineColor(kBlack);
      hRatio.SetLineWidth(2);
      hRatio.SetMarkerStyle(1);
      hRatio.SetFillStyle(1001);
      hRatio.SetFillColor(kGray);
   }

   return hRatio;
}

TH1F hist::getResidual(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   assert(et==COMB); // TODO implement COMB1/2
   TH1F hRes(h1);
   hRes.Add(&h2,-1);
   Double_t min = hRes.GetBinContent(hRes.GetMinimumBin());
   Double_t max = hRes.GetBinContent(hRes.GetMaximumBin());
   hRes.GetYaxis()->SetRangeUser(min,max*1.1);
   hRes.GetYaxis()->SetTitle(title);
   hRes.SetLineColor(kBlack);
   hRes.SetLineWidth(2);
   hRes.SetMarkerStyle(1);
   return hRes;
}

TH1F hist::getPull(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   assert(et==COMB); // TODO implement COMB1/2
   TH1F hPull=getResidual(h1,h2,title);
   for (int i=0; i<=h1.GetNbinsX()+1;i++){
      float const cont=hPull.GetBinContent(i);
      float const err=hPull.GetBinError(i);
      if (std::abs(cont)<1e-6 && err<1e-6){
         hPull.SetBinContent(i,0);
         hPull.SetBinError(i,0);
      } else {
         if (err<1e-6){
            debug_io<<"Zero error in non-emtpy bin!";
            throw;
         }
         hPull.SetBinContent(i,cont/err);
         hPull.SetBinError(i,1.);
      }
   }

   Double_t max = TMath::Max(
      std::abs(hPull.GetBinContent(hPull.GetMinimumBin())),
      std::abs(hPull.GetBinContent(hPull.GetMaximumBin())));
   max=TMath::Min(max,10.);
   max=(int)max+1.4;
   hPull.GetYaxis()->SetRangeUser(-max,max);
   hPull.SetFillStyle(1001);
   hPull.SetLineWidth(1);
   hPull.SetLineColor(kBlack);
   hPull.SetFillColor(kGray);
   return hPull;
}

TH1F hist::getRatio(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   return getRatio(h1,*(TH1F*)h2.GetStack()->Last(),title,et);
}

TH1F hist::getPull(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   return getPull(h1,*(TH1F*)h2.GetStack()->Last(),title,et);
}

TH1F hist::getRatio(THStack &h1,THStack &h2,TString title,ErrorType et)
{
   return getRatio(*(TH1F*)h1.GetStack()->Last(),*(TH1F*)h2.GetStack()->Last(),title,et);
}

TH1F hist::getPull(THStack &h1,THStack &h2,TString title,ErrorType et)
{
   return getPull(*(TH1F*)h1.GetStack()->Last(),*(TH1F*)h2.GetStack()->Last(),title,et);
}

TGraphErrors hist::getRatioGraph(TH1F const &h1,TH1F const &h2,TString title,ErrorType et)
{
   TH1F h = getRatio(h1,h2,title,et);
   return TGraphErrors(&h);
}

TGraphErrors hist::getRatioGraph(TH1F const &h1,THStack &h2,TString title,ErrorType et)
{
   TH1F h = getRatio(h1,h2,title,et);
   return TGraphErrors(&h);
}

TGraphErrors hist::getRatioGraph(THStack &h1,THStack &h2,TString title,ErrorType et)
{
   TH1F h = getRatio(h1,h2,title,et);
   return TGraphErrors(&h);
}

THStack hist::stackPrepend(THStack const& stOld, TH1F &h, Option_t *option){
   THStack st;
   st.SetTitle(stOld.GetTitle());
   st.Add(&h,option);
   for (int i=0; i<stOld.GetNhists(); i++){
      st.Add((TH1*)stOld.GetHists()->At(i),option);
   }
   return st;
}

std::pair<TH1F*,TH1F*> hist::getEnvelope(const TH1F* nominal, const std::vector<TH1F*> shifts){
   TH1F* hist_envelopeUP = (TH1F*)nominal->Clone();
   TH1F* hist_envelopeDOWN = (TH1F*)nominal->Clone();
   
   for (const TH1F* tempShift : shifts){
      for (int i=0; i<=tempShift->GetNbinsX(); i++){
         float content_shift = tempShift->GetBinContent(i);
         if (content_shift>hist_envelopeUP->GetBinContent(i)) hist_envelopeUP->SetBinContent(i,content_shift);
         else if (content_shift<hist_envelopeDOWN->GetBinContent(i)) hist_envelopeDOWN->SetBinContent(i,abs(content_shift));
      }
   }
   
   return std::make_pair(hist_envelopeDOWN,hist_envelopeUP);
}

std::pair<TH1F*,TH1F*> hist::getEnvelope(const TH1F* nominal, const std::vector<TH1F> &shifts){
   TH1F* hist_envelopeUP = (TH1F*)nominal->Clone();
   TH1F* hist_envelopeDOWN = (TH1F*)nominal->Clone();
   
   for (const TH1F tempShift : shifts){
      for (int i=0; i<=tempShift.GetNbinsX(); i++){
         float content_shift = tempShift.GetBinContent(i);
         if (content_shift>hist_envelopeUP->GetBinContent(i)) hist_envelopeUP->SetBinContent(i,content_shift);
         else if (content_shift<hist_envelopeDOWN->GetBinContent(i)) hist_envelopeDOWN->SetBinContent(i,abs(content_shift));
      }
   }
   
   return std::make_pair(hist_envelopeDOWN,hist_envelopeUP);
}

std::pair<TH2F*,TH2F*> hist::getEnvelope(const TH2F* nominal, const std::vector<TH2F*> shifts){
   TH2F* hist_envelopeUP = (TH2F*)nominal->Clone();
   TH2F* hist_envelopeDOWN = (TH2F*)nominal->Clone();
   
   for (const TH2F* tempShift : shifts){
      for (int i=0; i<=tempShift->GetNbinsX(); i++){
         for (int j=0; j<=tempShift->GetNbinsY(); j++){
            float content_shift = tempShift->GetBinContent(i,j);
            if (content_shift>hist_envelopeUP->GetBinContent(i,j)) hist_envelopeUP->SetBinContent(i,j,content_shift);
            else if (content_shift<hist_envelopeDOWN->GetBinContent(i,j)) hist_envelopeDOWN->SetBinContent(i,j,abs(content_shift));
         }
      }
   }
   
   return std::make_pair(hist_envelopeDOWN,hist_envelopeUP);
}

std::pair<TH2F*,TH2F*> hist::getEnvelope(const TH2F* nominal, const std::vector<TH2F> &shifts){
   TH2F* hist_envelopeUP = (TH2F*)nominal->Clone();
   TH2F* hist_envelopeDOWN = (TH2F*)nominal->Clone();
   
   for (const TH2F tempShift : shifts){
      for (int i=0; i<=tempShift.GetNbinsX(); i++){
         for (int j=0; j<=tempShift.GetNbinsY(); j++){
            float content_shift = tempShift.GetBinContent(i,j);
            if (content_shift>hist_envelopeUP->GetBinContent(i,j)) hist_envelopeUP->SetBinContent(i,j,content_shift);
            else if (content_shift<hist_envelopeDOWN->GetBinContent(i,j)) hist_envelopeDOWN->SetBinContent(i,j,abs(content_shift));
         }
      }
   }
   
   return std::make_pair(hist_envelopeDOWN,hist_envelopeUP);
}

// get graph with asym. errors from three histograms (shift=true if only shift and not shift+nominal is given)
TGraphAsymmErrors hist::getErrorGraph(TH1F* const &eDOWN, TH1F* const &eUP, TH1F* const &nominal, bool const shift, bool const eXzero){
   TGraphAsymmErrors asymmerrors(nominal);
   for (int i=0; i<=eUP->GetNbinsX(); i++){
      if (shift) {
         asymmerrors.SetPointEYhigh(i,eUP->GetBinContent(i+1));
         asymmerrors.SetPointEYlow(i,eDOWN->GetBinContent(i+1));
      }
      else{
         asymmerrors.SetPointEYhigh(i,abs(eUP->GetBinContent(i+1)-nominal->GetBinContent(i+1)));
         asymmerrors.SetPointEYlow(i,abs(eDOWN->GetBinContent(i+1)-nominal->GetBinContent(i+1)));
      }
      
      if (eXzero){
         asymmerrors.SetPointEXhigh(i,0.);
         asymmerrors.SetPointEXlow(i,0.);
      }
   }
   
   return asymmerrors;
}

TGraphAsymmErrors hist::getRatioAsymmGraph(TH1F const &down,TH1F const &up,TH1F const &nominal,TH1F const &denominator){

   TH1F* sysDown = (TH1F*)nominal.Clone();
   TH1F* sysUp = (TH1F*)nominal.Clone();
   sysDown->Add(&down,-1.);
   sysUp->Add(&up);
            
   TH1F ratio_nominal = getRatio(nominal,denominator,"data/MC",hist::ONLY1);
   TH1F ratio_systDown = getRatio(*sysDown,denominator,"data/MC",hist::ONLY1);
   TH1F ratio_systUp = getRatio(*sysUp,denominator,"data/MC",hist::ONLY1);
   
   return getErrorGraph(&ratio_systDown,&ratio_systUp,&ratio_nominal,false,true);
}
