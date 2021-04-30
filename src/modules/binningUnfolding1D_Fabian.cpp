//Raw script for 1D Binning Studies (Bachelor Thesis Fabian)

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TTreeReader.h>

#include <iomanip> 


Config const &cfg=Config::get();

extern "C"
void run()
{

   TString treeName="TTbar_diLepton";     //input sample name
   TFile file(TString::Format("/net/data_cms1b/user/dmeuser/top_analysis/%s/%s/minTrees/100.0/"+treeName+".root",cfg.year.Data(),cfg.treeVersion.Data()),"read");    //open input file
   TTreeReader reader("ttbar_res100.0/"+treeName, &file);      //open tree in root file
   
   //define variables taken from tree
   TTreeReaderValue<float> pTl1   (reader, "Lep1_pt");
   TTreeReaderValue<float> MET   (reader, "PuppiMET");
   TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
   TTreeReaderValue<float> MC_weight   (reader, "N");
   TTreeReaderValue<float> SF_weight   (reader, "SF");
   
   //define empty histograms
   hist::Histograms<TH1F> hs({treeName});
   hs.addHist("rec_selection/Lep1_pt"   ,";%pTl1;EventsBIN"           ,100,0,500);
   
   //define booleans for selection
   bool rec_selection;
   bool gen_selection;
   
   //set current sample for hist selection
   hs.setCurrentSample(treeName);
   
   //loop over events in tree
   int processEvents=cfg.processFraction*reader.GetEntries(true);
   int iEv=0;
   while (reader.Next()){
      iEv++;
      if (iEv>processEvents) break;
      if (iEv%(std::max(processEvents/10,1))==0){      //logging stuff
         io::log*".";
         io::log.flush();
      }
      
      //set weight for current event
      hs.setFillWeight((*SF_weight)*(*MC_weight));
     
      // set booleans for selection
      if(*PtNuNu>0) gen_selection=true;
      if(*MET>0) rec_selection=true;
     
      if(rec_selection){
         hs.fill("rec_selection/Lep1_pt",*pTl1);
      }
   }
   io::log<<"";
   
   //Merge overflow bins
   hs.mergeOverflow();
   
   file.Close();
   
   //Plot histograms
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"binningUnfolding1D_Fabian");
   TCanvas can;
   can.SetLogy();
   for(TString var : hs.getVariableNames()){
      TH1F* tempHist=hs.getHistogram(var,{treeName});
      
      tempHist->Draw("hist e");
      
      saver.save(can,var);
   }
}
