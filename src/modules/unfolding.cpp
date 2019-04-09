#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>
#include <TUnfoldDensity.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   
   for (TString cat: {"ee","emu","mumu"}){
      TString sSelection="genParticles/"+cat+"/MetVSgenMet";
      for (TString sSample:{"TTbar"}){
         TH2F* migrMatrix=(TH2F*)histReader.read<TH2F>(sSelection+"/"+sSample);
         TUnfoldDensity unfold(migrMatrix,TUnfold::kHistMapOutputVert);
      }
   }
}
