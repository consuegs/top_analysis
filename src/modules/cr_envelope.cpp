//Script to create envelopes for CR systematic and scale sample for mtop

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/systematics.hpp"

#include <TFile.h>
#include <TF1.h>
#include <memory>

Config const &cfg=Config::get();

void produce_cr_envelope()
{   
   io::RootFileReader histReader_nom(TString::Format("multiHists/%s/histograms_merged_%s.root","Nominal",cfg.treeVersion.Data()));
   io::RootFileReader histReader_CR1(TString::Format("multiHists/%s/histograms_merged_%s.root","CR1",cfg.treeVersion.Data()));
   io::RootFileReader histReader_CR2(TString::Format("multiHists/%s/histograms_merged_%s.root","CR2",cfg.treeVersion.Data()));
   io::RootFileReader histReader_ERDON(TString::Format("multiHists/%s/histograms_merged_%s.root","ERDON",cfg.treeVersion.Data()));
   
   io::RootFileSaver histSaver_ind_down(TString::Format("multiHists/%s/histograms_merged_%s.root","CR_ENVELOPE_IND_DOWN",cfg.treeVersion.Data()),"");
   io::RootFileSaver histSaver_ind_up(TString::Format("multiHists/%s/histograms_merged_%s.root","CR_ENVELOPE_IND_UP",cfg.treeVersion.Data()),"");
  
   // Store envelope per sample
   for (auto path : histReader_CR1.listPaths()){
      if(path.Contains("2d")){
         TH2F* hist_CR1 = histReader_CR1.read<TH2F>(path);
         TH2F* hist_CR2 = histReader_CR2.read<TH2F>(path.ReplaceAll("_CR1","_CR2"));
         TH2F* hist_ERDON = histReader_ERDON.read<TH2F>(path.ReplaceAll("_CR2","_ERDON"));
         TH2F* hist_nom = histReader_nom.read<TH2F>(path.ReplaceAll("_ERDON",""));
         
         std::pair<TH2F*,TH2F*> envelopes = hist::getEnvelope(hist_nom,{hist_CR1,hist_CR2,hist_ERDON});
         
         path.ReplaceAll("/distr","distr");
         histSaver_ind_down.save(*(envelopes.first),path+"_CR_ENVELOPE_IND_DOWN");
         histSaver_ind_up.save(*(envelopes.second),path+"_CR_ENVELOPE_IND_UP");
      }
      else {
         TH1F* hist_CR1 = histReader_CR1.read<TH1F>(path);
         TH1F* hist_CR2 = histReader_CR2.read<TH1F>(path.ReplaceAll("_CR1","_CR2"));
         TH1F* hist_ERDON = histReader_ERDON.read<TH1F>(path.ReplaceAll("_CR2","_ERDON"));
         TH1F* hist_nom = histReader_nom.read<TH1F>(path.ReplaceAll("_ERDON",""));
         
         std::pair<TH1F*,TH1F*> envelopes = hist::getEnvelope(hist_nom,{hist_CR1,hist_CR2,hist_ERDON});
         
         path.ReplaceAll("/distr","distr");
         histSaver_ind_down.save(*(envelopes.first),path+"_CR_ENVELOPE_IND_DOWN");
         histSaver_ind_up.save(*(envelopes.second),path+"_CR_ENVELOPE_IND_UP");
      }
   }
   
   histSaver_ind_down.closeFile();
   histSaver_ind_up.closeFile();
   
   // Store envelope for ttbar sum (will be saved in signal sample due to plotting reasons)
   io::RootFileSaver histSaver_down(TString::Format("multiHists/%s/histograms_merged_%s.root","CR_ENVELOPE_DOWN",cfg.treeVersion.Data()),"");
   io::RootFileSaver histSaver_up(TString::Format("multiHists/%s/histograms_merged_%s.root","CR_ENVELOPE_UP",cfg.treeVersion.Data()),"");
   
   for (auto path : histReader_CR1.listPaths(true)){
      if(path.Contains("2d")){
         TH2F* hist_CR1 = histReader_CR1.read<TH2F>(path+"/TTbar_diLepton_CR1");
         TH2F* hist_CR2 = histReader_CR2.read<TH2F>(path+"/TTbar_diLepton_CR2");
         TH2F* hist_ERDON = histReader_ERDON.read<TH2F>(path+"/TTbar_diLepton_ERDON");
         TH2F* hist_nom = histReader_nom.read<TH2F>(path+"/TTbar_diLepton");
         
         for(const TString sample : {"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"}){
            hist_CR1->Add(histReader_CR1.read<TH2F>(path+"/"+sample+"_CR1"));
            hist_CR2->Add(histReader_CR2.read<TH2F>(path+"/"+sample+"_CR2"));
            hist_ERDON->Add(histReader_ERDON.read<TH2F>(path+"/"+sample+"_ERDON"));
            hist_nom->Add(histReader_nom.read<TH2F>(path+"/"+sample));
         }
         
         std::pair<TH2F*,TH2F*> envelopes = hist::getEnvelope(hist_nom,{hist_CR1,hist_CR2,hist_ERDON});
         
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*(envelopes.first),path+"/TTbar_diLepton_CR_ENVELOPE_DOWN");
         histSaver_up.save(*(envelopes.second),path+"/TTbar_diLepton_CR_ENVELOPE_UP");
      }
      else {
         TH1F* hist_CR1 = histReader_CR1.read<TH1F>(path+"/TTbar_diLepton_CR1");
         TH1F* hist_CR2 = histReader_CR2.read<TH1F>(path+"/TTbar_diLepton_CR2");
         TH1F* hist_ERDON = histReader_ERDON.read<TH1F>(path+"/TTbar_diLepton_ERDON");
         TH1F* hist_nom = histReader_nom.read<TH1F>(path+"/TTbar_diLepton");
         
         for(const TString sample : {"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"}){
            hist_CR1->Add(histReader_CR1.read<TH1F>(path+"/"+sample+"_CR1"));
            hist_CR2->Add(histReader_CR2.read<TH1F>(path+"/"+sample+"_CR2"));
            hist_ERDON->Add(histReader_ERDON.read<TH1F>(path+"/"+sample+"_ERDON"));
            hist_nom->Add(histReader_nom.read<TH1F>(path+"/"+sample));
         }
         
         std::pair<TH1F*,TH1F*> envelopes = hist::getEnvelope(hist_nom,{hist_CR1,hist_CR2,hist_ERDON});
         
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*(envelopes.first),path+"/TTbar_diLepton_CR_ENVELOPE_DOWN");
         histSaver_up.save(*(envelopes.second),path+"/TTbar_diLepton_CR_ENVELOPE_UP");
      }
      
   }
   
   histReader_nom.closeFile();
   histReader_CR1.closeFile();
   histReader_CR2.closeFile();
   histReader_ERDON.closeFile();
   histSaver_down.closeFile();
   histSaver_up.closeFile();
  
}

void produce_topmass()
{
   io::RootFileReader histReader_nom_ind(TString::Format("multiHists/%s/histograms_merged_%s.root","Nominal",cfg.treeVersion.Data()));
   io::RootFileReader histReader_low_ind(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP169p5",cfg.treeVersion.Data()));
   io::RootFileReader histReader_high_ind(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP175p5",cfg.treeVersion.Data()));
   
   io::RootFileSaver histSaver_ind_down(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP_IND_DOWN",cfg.treeVersion.Data()),"");
   io::RootFileSaver histSaver_ind_up(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP_IND_UP",cfg.treeVersion.Data()),"");
   
   // Store shift per sample
   for (auto path : histReader_low_ind.listPaths()){
      if(path.Contains("2d")){
         TH2F* hist_low = (TH2F*)histReader_low_ind.read<TH2F>(path);
         TH2F* hist_high = (TH2F*)histReader_high_ind.read<TH2F>(path.ReplaceAll("_MTOP169p5","_MTOP175p5"));
         TH2F* hist_nom = (TH2F*)histReader_nom_ind.read<TH2F>(path.ReplaceAll("_MTOP175p5",""));
         
         hist_low->Add(hist_nom,-1.);
         hist_low->Scale(1/3.);
         hist_low->Add(hist_nom);
         
         hist_high->Add(hist_nom,-1.);
         hist_high->Scale(1/3.);
         hist_high->Add(hist_nom);
         
         path.ReplaceAll("/distr","distr");
         histSaver_ind_down.save(*hist_low,path+"_MTOP_IND_DOWN");
         histSaver_ind_up.save(*hist_high,path+"_MTOP_IND_UP");
      }
      else {
         TH1F* hist_low = (TH1F*)histReader_low_ind.read<TH1F>(path);
         TH1F* hist_high = (TH1F*)histReader_high_ind.read<TH1F>(path.ReplaceAll("_MTOP169p5","_MTOP175p5"));
         TH1F* hist_nom = (TH1F*)histReader_nom_ind.read<TH1F>(path.ReplaceAll("_MTOP175p5",""));
         
         hist_low->Add(hist_nom,-1.);
         hist_low->Scale(1/3.);
         hist_low->Add(hist_nom);
         
         hist_high->Add(hist_nom,-1.);
         hist_high->Scale(1/3.);
         hist_high->Add(hist_nom);
         
         path.ReplaceAll("/distr","distr");
         histSaver_ind_down.save(*hist_low,path+"_MTOP_IND_DOWN");
         histSaver_ind_up.save(*hist_high,path+"_MTOP_IND_UP");
      }
   }
   
   histReader_nom_ind.closeFile();
   histReader_low_ind.closeFile();
   histReader_high_ind.closeFile();
   histSaver_ind_down.closeFile();
   histSaver_ind_up.closeFile();
   
   io::RootFileReader histReader_nom(TString::Format("multiHists/%s/histograms_merged_%s.root","Nominal",cfg.treeVersion.Data()));
   io::RootFileReader histReader_low(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP169p5",cfg.treeVersion.Data()));
   io::RootFileReader histReader_high(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP175p5",cfg.treeVersion.Data()));
   
   io::RootFileSaver histSaver_down(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP_DOWN",cfg.treeVersion.Data()),"");
   io::RootFileSaver histSaver_up(TString::Format("multiHists/%s/histograms_merged_%s.root","MTOP_UP",cfg.treeVersion.Data()),"");
   
   
   // Store shift for ttbar sum (will be saved in signal sample due to plotting reasons)
   for (auto path : histReader_low.listPaths(true)){
      if(path.Contains("2d")){
         TH2F* hist_low = histReader_low.read<TH2F>(path+"/TTbar_diLepton_MTOP169p5");
         TH2F* hist_high = histReader_high.read<TH2F>(path+"/TTbar_diLepton_MTOP175p5");
         TH2F* hist_nom = histReader_nom.read<TH2F>(path+"/TTbar_diLepton");
               
         for(const TString sample : {"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"}){
            hist_low->Add(histReader_low.read<TH2F>(path+"/"+sample+"_MTOP169p5"));
            hist_high->Add(histReader_high.read<TH2F>(path+"/"+sample+"_MTOP175p5"));
            hist_nom->Add(histReader_nom.read<TH2F>(path+"/"+sample));
         }
         
         hist_low->Add(hist_nom,-1.);
         hist_low->Scale(1/3.);
         hist_low->Add(hist_nom);
         
         hist_high->Add(hist_nom,-1.);
         hist_high->Scale(1/3.);
         hist_high->Add(hist_nom);
               
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*hist_low,path+"/TTbar_diLepton_MTOP_DOWN");
         histSaver_up.save(*hist_high,path+"/TTbar_diLepton_MTOP_UP");
      }
      else {
         TH1F* hist_low = histReader_low.read<TH1F>(path+"/TTbar_diLepton_MTOP169p5");
         TH1F* hist_high = histReader_high.read<TH1F>(path+"/TTbar_diLepton_MTOP175p5");
         TH1F* hist_nom = histReader_nom.read<TH1F>(path+"/TTbar_diLepton");
               
         for(const TString sample : {"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"}){
            hist_low->Add(histReader_low.read<TH1F>(path+"/"+sample+"_MTOP169p5"));
            hist_high->Add(histReader_high.read<TH1F>(path+"/"+sample+"_MTOP175p5"));
            hist_nom->Add(histReader_nom.read<TH1F>(path+"/"+sample));
         }
         
         hist_low->Add(hist_nom,-1.);
         hist_low->Scale(1/3.);
         hist_low->Add(hist_nom);
         
         hist_high->Add(hist_nom,-1.);
         hist_high->Scale(1/3.);
         hist_high->Add(hist_nom);
               
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*hist_low,path+"/TTbar_diLepton_MTOP_DOWN");
         histSaver_up.save(*hist_high,path+"/TTbar_diLepton_MTOP_UP");
      }
   }
   
   histReader_nom.closeFile();
   histReader_low.closeFile();
   histReader_high.closeFile();
   histSaver_down.closeFile();
   histSaver_up.closeFile();
   
}

void produce_PDF_envelope(){  //currently only for 1D hists
      
   std::map<TString,std::vector<TH1F>> histMap;
   std::map<TString,std::vector<TH2F>> histMap2D;
   std::vector<TString> paths;
   
   for (int i=1; i<=50; i++){    // create reader for each shift
      std::cout<<"Reading PDF variation "+std::to_string(i)<<std::endl;
      io::RootFileReader histReader_up(TString::Format("multiHists/PDF_%i_UP/histograms_merged_%s.root",i,cfg.treeVersion.Data()));
      io::RootFileReader histReader_down(TString::Format("multiHists/PDF_%i_DOWN/histograms_merged_%s.root",i,cfg.treeVersion.Data()));
      for (auto path : histReader_up.listPaths(true)){
         if(path.Contains("2d")){
            histMap2D[path].push_back(*(histReader_up.read<TH2F>(path+"/TTbar_diLepton")));
            histMap2D[path].push_back(*(histReader_up.read<TH2F>(path+"/TTbar_diLepton")));
         }
         else{
            histMap[path].push_back(*(histReader_up.read<TH1F>(path+"/TTbar_diLepton")));
            histMap[path].push_back(*(histReader_down.read<TH1F>(path+"/TTbar_diLepton")));
         }
         if (i==1) paths.push_back(path);
      }
   }
   
   io::RootFileReader histReader_nom(TString::Format("multiHists/%s/histograms_merged_%s.root","Nominal",cfg.treeVersion.Data()));
   
   io::RootFileSaver histSaver_down(TString::Format("multiHists/%s/histograms_merged_%s.root","PDF_ENVELOPE_DOWN",cfg.treeVersion.Data()),"");
   io::RootFileSaver histSaver_up(TString::Format("multiHists/%s/histograms_merged_%s.root","PDF_ENVELOPE_UP",cfg.treeVersion.Data()),"");
   
   for (auto path : paths){   // loop over all histograms
      if(path.Contains("2d")){
         TH2F* hist_nom = histReader_nom.read<TH2F>(path+"/TTbar_diLepton");

         std::pair<TH2F*,TH2F*> envelopes = hist::getEnvelope(hist_nom,histMap2D[path]);
         
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*(envelopes.first),path+"/TTbar_diLepton");
         histSaver_up.save(*(envelopes.second),path+"/TTbar_diLepton");
      }
      else {
         TH1F* hist_nom = histReader_nom.read<TH1F>(path+"/TTbar_diLepton");
         std::pair<TH1F*,TH1F*> envelopes = hist::getEnvelope(hist_nom,histMap[path]);
         
         path.ReplaceAll("/distr","distr");
         histSaver_down.save(*(envelopes.first),path+"/TTbar_diLepton");
         histSaver_up.save(*(envelopes.second),path+"/TTbar_diLepton");
      }
      
   }
}
   
   

extern "C"
void run()
{   
   // ~produce_cr_envelope();
   produce_topmass();
   // ~produce_PDF_envelope();
  
}
