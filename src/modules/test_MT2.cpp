//Script to test test the MT2 variable

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/MT2Functor.h"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TGraph.h>
#include <math.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   MT2Functor fctMT2_;
   double pa[3];
   double pb[3];
   double pmiss[3];

   pa[0]=0; pa[1]=TMath::Cos(M_PI/-4.); pa[2]=TMath::Sin(M_PI/-4.);
   pb[0]=0; pb[1]=TMath::Cos(M_PI/4.); pb[2]=TMath::Sin(M_PI/4.);
   
   TGraph curve;
   
   for (double phi=-1*M_PI; phi<M_PI; phi+=0.01){
      pmiss[0]=0; pmiss[1]=TMath::Cos(phi); pmiss[2]=TMath::Sin(phi);
      
      fctMT2_.set_mn(0.);
      fctMT2_.set_momenta(pa,pb,pmiss);
      
      curve.SetPoint(curve.GetN(),phi,fctMT2_.get_mt2());
   }
   curve.GetXaxis()->SetTitle("#phi(met)");
   curve.GetYaxis()->SetTitle("MT2");
   curve.SaveAs("test.root");
}
