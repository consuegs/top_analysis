//Script to calulcare Rout/in

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>

Config const &cfg=Config::get();

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"Routin");
   
   TCanvas can;
   can.cd();
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()));
   
   std::vector<double> metBins = {0,20,40,60,80,100,120,140,160,195,230,400};
   std::vector<double> phiBins = {0,0.35,0.7,1.05,1.4,2.27,3.141};
   
   TH1F* ee_in = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/ee/met/DrellYan_HT");
   TH1F* mumu_in = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/mumu/met/DrellYan_HT");
   
   TH1F* ee_out = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/ee/met/DrellYan_HT");
   TH1F* mumu_out = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/mumu/met/DrellYan_HT");
   
   TH1F* ee_in_Puppi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/ee/PuppiMet/DrellYan_HT");
   TH1F* mumu_in_Puppi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/mumu/PuppiMet/DrellYan_HT");
   
   TH1F* ee_out_Puppi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/ee/PuppiMet/DrellYan_HT");
   TH1F* mumu_out_Puppi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/mumu/PuppiMet/DrellYan_HT");
   
   TH1F* ee_in_phi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/ee/dphi_metNearLep/DrellYan_HT");
   TH1F* mumu_in_phi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/mumu/dphi_metNearLep/DrellYan_HT");
   
   TH1F* ee_out_phi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/ee/dphi_metNearLep/DrellYan_HT");
   TH1F* mumu_out_phi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/mumu/dphi_metNearLep/DrellYan_HT");
   
   TH1F* ee_in_COMB = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/ee/met/DrellYanCOMB");
   TH1F* mumu_in_COMB = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/mumu/met/DrellYanCOMB");
   TH1F* mumu_in_COMB_tightRIso = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau_tightRIsoMu/mumu/met/DrellYanCOMB");
   
   TH1F* ee_out_COMB = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/ee/met/DrellYanCOMB");
   TH1F* mumu_out_COMB = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/mumu/met/DrellYanCOMB");
   TH1F* mumu_out_COMB_tightRIso = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau_tightRIsoMu/mumu/met/DrellYanCOMB");
   
   TH1F* ee_in_COMB_phi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/ee/dphi_metNearLep/DrellYanCOMB");
   TH1F* mumu_in_COMB_phi = histReader.read<TH1F>("distributions_Philipp100.0/Zpeak_noTau/mumu/dphi_metNearLep/DrellYanCOMB");
   
   TH1F* ee_out_COMB_phi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/ee/dphi_metNearLep/DrellYanCOMB");
   TH1F* mumu_out_COMB_phi = histReader.read<TH1F>("distributions_Philipp100.0/baseline_noTau/mumu/dphi_metNearLep/DrellYanCOMB");
   
   *ee_in = hist::rebinned(*ee_in,metBins);
   *ee_out = hist::rebinned(*ee_out,metBins);
   *mumu_in = hist::rebinned(*mumu_in,metBins);
   *mumu_out = hist::rebinned(*mumu_out,metBins);
   
   *ee_in_COMB = hist::rebinned(*ee_in_COMB,metBins);
   *ee_out_COMB = hist::rebinned(*ee_out_COMB,metBins);
   *mumu_in_COMB = hist::rebinned(*mumu_in_COMB,metBins);
   *mumu_out_COMB = hist::rebinned(*mumu_out_COMB,metBins);
   *mumu_in_COMB_tightRIso = hist::rebinned(*mumu_in_COMB_tightRIso,metBins);
   *mumu_out_COMB_tightRIso = hist::rebinned(*mumu_out_COMB_tightRIso,metBins);
   
   *ee_in_Puppi = hist::rebinned(*ee_in_Puppi,metBins);
   *ee_out_Puppi = hist::rebinned(*ee_out_Puppi,metBins);
   *mumu_in_Puppi = hist::rebinned(*mumu_in_Puppi,metBins);
   *mumu_out_Puppi = hist::rebinned(*mumu_out_Puppi,metBins);
   
   *ee_in_phi = hist::rebinned(*ee_in_phi,phiBins);
   *ee_out_phi = hist::rebinned(*ee_out_phi,phiBins);
   *mumu_in_phi = hist::rebinned(*mumu_in_phi,phiBins);
   *mumu_out_phi = hist::rebinned(*mumu_out_phi,phiBins);
   
   *ee_in_COMB_phi = hist::rebinned(*ee_in_COMB_phi,phiBins);
   *ee_out_COMB_phi = hist::rebinned(*ee_out_COMB_phi,phiBins);
   *mumu_in_COMB_phi = hist::rebinned(*mumu_in_COMB_phi,phiBins);
   *mumu_out_COMB_phi = hist::rebinned(*mumu_out_COMB_phi,phiBins);
   
   ee_out->Divide(ee_in);
   mumu_out->Divide(mumu_in);
   
   ee_out_COMB->Divide(ee_in_COMB);
   mumu_out_COMB->Divide(mumu_in_COMB);
   mumu_out_COMB_tightRIso->Divide(mumu_in_COMB_tightRIso);
   
   ee_out_Puppi->Divide(ee_in_Puppi);
   mumu_out_Puppi->Divide(mumu_in_Puppi);
   
   ee_out_phi->Divide(ee_in_phi);
   mumu_out_phi->Divide(mumu_in_phi);
   
   ee_out_COMB_phi->Divide(ee_in_COMB_phi);
   mumu_out_COMB_phi->Divide(mumu_in_COMB_phi);
   
   saver.save(*ee_out,"R_ee");
   saver.save(*mumu_out,"R_mumu");
   
   saver.save(*ee_out_COMB,"COMB/R_ee");
   saver.save(*mumu_out_COMB,"COMB/R_mumu");
   saver.save(*mumu_out_COMB_tightRIso,"COMB/R_mumu_tightRIso");
   
   saver.save(*ee_out_Puppi,"Puppi/R_ee");
   saver.save(*mumu_out_Puppi,"Puppi/R_mumu");
   
   saver.save(*ee_out_phi,"Phi/R_ee");
   saver.save(*mumu_out_phi,"Phi/R_mumu");
   
   saver.save(*ee_out_COMB_phi,"COMB/Phi/R_ee");
   saver.save(*mumu_out_COMB_phi,"COMB/Phi/R_mumu");
   
}
