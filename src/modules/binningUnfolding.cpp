//Script to test different binnings regarding stability and purity for unfolding based on minimal ttbarMC tree "~/top_analysis/framework/output/ttbar_res100.0.root"

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
#include <TFitResult.h>
#include <math.h>

#include <iomanip> 

void saveProfileFrom2D(TH2F* hist, TString name,io::RootFileSaver* saver){
   TProfile profile;
   profile=*(hist->ProfileX("SF_met"));
   saver->save(*hist,name);
   saver->save(profile,name+"_profile");
}


Config const &cfg=Config::get();

extern "C"
void run()
{
   double scaleFactor = 1.0;
   // ~double scaleFactor = 0.912;
   
   std::vector<float> met_bins1={0,100,200,300,400};
   std::vector<float> phi_bins1={0,0.8,1.6,2.4,3.2};
   
   std::vector<float> met_bins2={0,50,100,200,400};
   std::vector<float> phi_bins2={0,0.4,0.8,1.6,3.2};
   
   std::vector<float> met_bins3={0,70,140,250,400};
   std::vector<float> phi_bins3={0,0.4,0.8,1.2,3.14};
   
   // ~std::vector<float> met_bins4={0,40,120,230,400};
   // ~std::vector<float> met_bins4={0,40,80,120,230,400};
   std::vector<float> met_bins4={0,40,80,120,160,230,400};
   std::vector<float> phi_bins4={0,0.7,1.4,3.14};
   
   TH2F N_gen;
   TH2F N_rec;
   TH2F N_genrec;
   
   TH2F Evt_gen;  //Histograms with correct weights and not just MC N events
   TH2F Evt_rec;
   TH2F Evt_genrec;
   
   TH2F eff_gen;  //Histograms to calculate correct efficiency
   TH2F eff_gen_recSom;
   
   TH2F migration;   //Histogram to study the migration
   
   TProfile2D TopPt_profile;   //Profile to save mean lead top pT per bin
   
   int bin_gen;
   int bin_rec;
   
   int numberBinningSchemeMet=1;
   int numberBinningSchemePhi=1;
   
   float diffMET=0.;
   float diffPHI=0.;
   
   float truePtNuNu=0.;
   float trueDPhi_reco=0.;
   float MET_reco=0.;
   
   float diffMET_normed=0.;
   float diffPHI_normed=0.;
   
   TH1F dPhiMETnearJet_badReso("dphi_metnearJet"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETfarJet_badReso("dphi_metfarJet"   ,";|#Delta#phi|(p_{T}^{miss},farest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETleadJet_badReso("dphi_metleadJet"   ,";|#Delta#phi|(p_{T}^{miss},lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETlead2Jet_badReso("dphi_metlead2Jet"   ,";|#Delta#phi|(p_{T}^{miss},2nd lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETbJet_badReso("dphi_metbJet"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiLep1Lep2_badReso("dphi_Lep1Lep2"   ,";|#Delta#phi|(l_1,l_2);EventsBIN"           ,100,0,3.2);
   TH1F METsig_badReso("METsig"   ,";METsig;EventsBIN"           ,100,0,200);
   TH1F diffPuppi_badReso("diffPuppi"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,0.5);
   
   TH1F dPhiMETnearJet_goodReso("dphi_metnearJet_good"   ,";|#Delta#phi|(p_{T}^{miss},nearest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETfarJet_goodReso("dphi_metfarJet_good"   ,";|#Delta#phi|(p_{T}^{miss},farest jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETleadJet_goodReso("dphi_metleadJet_good"   ,";|#Delta#phi|(p_{T}^{miss},lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETlead2Jet_goodReso("dphi_metlead2Jet_good"   ,";|#Delta#phi|(p_{T}^{miss},2nd lead jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiMETbJet_goodReso("dphi_metbJet_good"   ,";|#Delta#phi|(p_{T}^{miss},b jet);EventsBIN"           ,100,0,3.2);
   TH1F dPhiLep1Lep2_goodReso("dphi_Lep1Lep2_good"   ,";|#Delta#phi|(l_1,l_2);EventsBIN"           ,100,0,3.2);
   TH1F METsig_goodReso("METsig_good"   ,";METsig;EventsBIN"           ,100,0,200);
   TH1F diffPuppi_goodReso("diffPuppi_good"   ,";|PFp_{T}^{miss}-PUPPIp_{T}^{miss}|/H_{T};EventsBIN"           ,100,0,0.5);
   
   //Histograms for 2D unfolding (current benchmark)
   std::vector<float>met_bins_unfold={0,20,40,60,80,100,120,140,160,195,230,400};
   std::vector<float>phi_bins_unfold={0,0.35,0.7,1.05,1.4,2.27,3.14};
   TH2F response("response",";binNumber;EventsBIN",66,-0.5,65.5,18,-0.5,17.5);
   TH2F response_sameBins("response_sameBins",";binNumber;EventsBIN",18,-0.5,17.5,18,-0.5,17.5);
   TH1F trueDistributions("trueDistributions",";binNumber;EventsBIN",18,-0.5,17.5);
   TH1F recoDistributions("recoDistributions",";binNumber;EventsBIN",66,-0.5,65.5);
   TH2F reco2D=hist::fromWidths_2d("reco2D",";p_{T}^{miss}(GeV);|#Delta#phi|(p_{T}^{miss},nearest l);",met_bins_unfold,hist::getWidths(met_bins_unfold),phi_bins_unfold,hist::getWidths(phi_bins_unfold));
   
   //Hist for METdiff vs. measured MET
   TH2F METdiff("METdiff_MET",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",100,0,400,600,-1,3);
   TH2F METdiff_PtNuNu("METdiff_PtNuNu",";p_{T}^{#nu#nu} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",100,0,400,600,-1,3);
   TH2F METdiff_GenMet("METdiff_GenMet",";genMet (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",100,0,400,600,-1,3);
   TH2F METdiffpuppi_PuppiMet("METdiffpuppi_PuppiMet",";PuppiMet (GeV);1+(p_{T}^{#nu#nu(+BSM)}-PuppiMet)/PuppiMet",100,0,400,600,-1,3);
   TH2F METdiffgen("METdiffgen",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-genMet)/genMet",100,0,400,600,-1,3);
   
   TH2F METdiff_phi("METdiff_MET_phi",";|#Delta#phi|(p_{T}^{miss},nearest l);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",100,0,3.14,600,-1,3);
   TH2F METdiffgen_phi("METdiffgen_phi",";|#Delta#phi|(p_{T}^{miss},nearest l);1+(p_{T}^{#nu#nu(+BSM)}-genMet)/genMet",100,0,3.14,600,-1,3);
   
   TH2F METdiff_phi1("METdiff_MET_phi1",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,400,600,-1,3);
   TH2F METdiff_phi2("METdiff_MET_phi2",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,400,600,-1,3);
   TH2F METdiff_phi3("METdiff_MET_phi3",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,400,600,-1,3);
   
   TH2F METdiffgen_phi1("METdiffgen_MET_phi1",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,400,600,-1,3);
   TH2F METdiffgen_phi2("METdiffgen_MET_phi2",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,400,600,-1,3);
   TH2F METdiffgen_phi3("METdiffgen_MET_phi3",";p_{T}^{miss} (GeV);1+(p_{T}^{#nu#nu(+BSM)}-genMet)/genMet}",50,0,400,600,-1,3);
   
   TH2F dPhiMETpTnunu_phi1("dPhiMETpTnunu_MET_phi1",";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},p_{T}^{#nu#nu})",50,0,400,600,-3,3);
   TH2F dPhiMETpTnunu_phi2("dPhiMETpTnunu_MET_phi2",";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},p_{T}^{#nu#nu})",50,0,400,600,-3,3);
   TH2F dPhiMETpTnunu_phi3("dPhiMETpTnunu_MET_phi3",";p_{T}^{miss} (GeV);|#Delta#phi|(p_{T}^{miss},p_{T}^{#nu#nu})",50,0,400,600,-3,3);
   
   TH2F TopPt_dPhi("TopPt_dPhi",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);p_{T}^{t1} (GeV)",100,0,3.14,1000,0,1000);
   TH2F TopPt_dPhi_met1("TopPt_dPhi_met1",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);p_{T}^{t1} (GeV)",50,0,3.14,1000,0,1000);
   TH2F TopPt_dPhi_met2("TopPt_dPhi_met2",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);p_{T}^{t1} (GeV)",50,0,3.14,1000,0,1000);
   TH2F TopPt_dPhi_met3("TopPt_dPhi_met3",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);p_{T}^{t1} (GeV)",50,0,3.14,1000,0,1000);
   TH2F TopPt_dPhi_met4("TopPt_dPhi_met4",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);p_{T}^{t1} (GeV)",50,0,3.14,1000,0,1000);
   
   TH2F METdiff_dPhiLep_phi1("METdiff_dPhiLep_phi1",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiLep_phi2("METdiff_dPhiLep_phi2",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiLep_phi3("METdiff_dPhiLep_phi3",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   
   TH2F METdiff_dPhiLep_met1("METdiff_dPhiLep_met1",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiLep_met2("METdiff_dPhiLep_met2",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiLep_met3("METdiff_dPhiLep_met3",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiLep_met4("METdiff_dPhiLep_met4",";|#Delta#phi|(l1,l2);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   
   TH2F METdiff_dPhiNuNu_phi1("METdiff_dPhiNuNu_phi1",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiNuNu_phi2("METdiff_dPhiNuNu_phi2",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiNuNu_phi3("METdiff_dPhiNuNu_phi3",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   
   TH2F METdiff_dPhiNuNu_met1("METdiff_dPhiNuNu_met1",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiNuNu_met2("METdiff_dPhiNuNu_met2",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiNuNu_met3("METdiff_dPhiNuNu_met3",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiNuNu_met4("METdiff_dPhiNuNu_met4",";|#Delta#phi|(#nu,#nu);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   
   TH2F METdiff_dPhiNuLep("METdiff_dPhiNuLep",";|#Delta#phi|(p_{T}^{#nu#nu},nearest l);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiff_dPhiMETLep("METdiff_dPhiMETLep",";|#Delta#phi|(p_{T}^{miss},nearest l);1+(p_{T}^{#nu#nu(+BSM)}-p_{T}^{miss})/p_{T}^{miss}",50,0,3.14,600,-1,3);
   TH2F METdiffgen_dPhiMETLep("METdiffgen_dPhiMETLep",";|#Delta#phi|(p_{T}^{miss},nearest l);p_{T}^{#nu#nu(+BSM)/genMet",50,0,3.14,600,-1,3);
   
   TH2F PtNuNuDiffGenMETRel_dPhiMETLep("PtNuNuDiffGenMETRel_dPhiMETLep",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffGenMETRel_dPhiMETLep_met1("PtNuNuDiffGenMETRel_dPhiMETLep_met1",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffGenMETRel_dPhiMETLep_met2("PtNuNuDiffGenMETRel_dPhiMETLep_met2",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffGenMETRel_dPhiMETLep_met3("PtNuNuDiffGenMETRel_dPhiMETLep_met3",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffGenMETRel_dPhiMETLep_met4("PtNuNuDiffGenMETRel_dPhiMETLep_met4",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep("GenMetDiffMETRel_dPhiMETLep",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met1("GenMetDiffMETRel_dPhiMETLep_met1",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met2("GenMetDiffMETRel_dPhiMETLep_met2",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met3("GenMetDiffMETRel_dPhiMETLep_met3",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met4("GenMetDiffMETRel_dPhiMETLep_met4",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{miss}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffMETRel_dPhiMETLep("PtNuNuDiffMETRel_dPhiMETLep",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffMETRel_dPhiMETLep_met1("PtNuNuDiffMETRel_dPhiMETLep_met1",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffMETRel_dPhiMETLep_met2("PtNuNuDiffMETRel_dPhiMETLep_met2",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffMETRel_dPhiMETLep_met3("PtNuNuDiffMETRel_dPhiMETLep_met3",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F PtNuNuDiffMETRel_dPhiMETLep_met4("PtNuNuDiffMETRel_dPhiMETLep_met4",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   
   TH2F GenMetDiffMETRel_dPhiMETLep_met120("PtNuNuDiffMETRel_dPhiMETLep_met120",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met120_ee("PtNuNuDiffMETRel_dPhiMETLep_met120_ee",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met120_emu("PtNuNuDiffMETRel_dPhiMETLep_met120_emu",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   TH2F GenMetDiffMETRel_dPhiMETLep_met120_mumu("PtNuNuDiffMETRel_dPhiMETLep_met120_mumu",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   
   TH2F GenMetDiffPuppiMETRel_dPhiMETLep_met120("GenMetDiffPuppiMETRel_dPhiMETLep_met120",";|#Delta#phi|(p_{T}^{miss},nearest l);|p_{T}^{#nu#nu(+BSM)}-GenMet|/genMET",50,0,3.14,6000,-5,5);
   
   TH2F METsig_dPhiMETLep_met120("METsig_dPhiMETLep_met120",";|#Delta#phi|(p_{T}^{miss},nearest l);metSig",50,0,3.14,6000,0,1000);
   
   // Plots for last bin
   // ~TH2F GenMETvsPTnunu_lastBin("GenMETvsPTnunu_lastBin",";p_{T}^{#nu#nu} (GeV);genMET (GeV);Events/Bin",100,0,400,100,230,400);
   TH2F GenMETvsPTnunu_lastBin("GenMETvsPTnunu_lastBin",";p_{T}^{#nu#nu} (GeV);genMET (GeV);Events/Bin",100,0,400,100,0,400);
   TH2F GenMETvsPTnunu_lastBin_nonPromptNu("GenMETvsPTnunu_lastBin_nonPromptNu",";p_{T}^{#nu#nu} (GeV);genMET (GeV);Events/Bin",100,0,400,100,0,400);
   TH2F GenMETvsPTnunu_lastBin_PromptNu("GenMETvsPTnunu_lastBin_PromptNu",";p_{T}^{#nu#nu} (GeV);genMET (GeV);Events/Bin",100,0,400,100,0,400);
   TH1F dPhiReco_lastBin_PromptNu("dPhiReco_lastBin_PromptNu",";|#Delta#phi|(p_{T}^{miss},nearest l);Events/Bin",30,0,3.14);
   TH1F dPhiReco_lastBin_nonPromptNu("dPhiReco_lastBin_nonPromptNu",";|#Delta#phi|(p_{T}^{miss},nearest l);Events/Bin",30,0,3.14);
   TH1F maxMuonF_bJET_lastBin_PromptNu("maxMuonF_bJET_lastBin_PromptNu",";max(MuonFraction_bJet);Events/Bin",100,0,1);
   TH1F maxMuonF_bJET_lastBin_nonPromptNu("maxMuonF_bJET_lastBin_nonPromptNu",";max(MuonFraction_bJet);Events/Bin",100,0,1);
   TH1F maxElectronF_bJET_lastBin_PromptNu("maxElectronF_bJET_lastBin_PromptNu",";max(ElectronFraction_bJet);Events/Bin",100,0,1);
   TH1F maxElectronF_bJET_lastBin_nonPromptNu("maxElectronF_bJET_lastBin_nonPromptNu",";max(ElectronFraction_bJet);Events/Bin",100,0,1);
   TH2F maxMuonF_bJET_vs_diff_lastBin("maxMuonF_bJET_vs_diff_lastBin",";max(ElectronFraction_bJet);|genMET-pTnunu|;Events/Bin",100,0,1,100,0,100);
   TH2F maxElectronF_bJET_vs_diff_lastBin("maxElectronF_bJET_vs_diff_lastBin",";max(ElectronFraction_bJet);|genMET-pTnunu|;Events/Bin",100,0,1,100,0,100);
   TH1F maxMuonEnergy_bJET_lastBin_PromptNu("maxMuonEnergy_bJET_lastBin_PromptNu",";max(MuonEnergy_bJet) (GeV);Events/Bin",100,0,100);
   TH1F maxMuonEnergy_bJET_lastBin_nonPromptNu("maxMuonEnergy_bJET_lastBin_nonPromptNu",";max(MuonEnergy_bJet) (GeV);Events/Bin",100,0,100);
   TH1F maxElectronEnergy_bJET_lastBin_PromptNu("maxElectronEnergy_bJET_lastBin_PromptNu",";max(ElectronEnergy_bJet) (GeV);Events/Bin",100,0,100);
   TH1F maxElectronEnergy_bJET_lastBin_nonPromptNu("maxElectronEnergy_bJET_lastBin_nonPromptNu",";max(ElectronEnergy_bJet) (GeV);Events/Bin",100,0,100);
   TH2F maxMuonEnergy_bJET_vs_diff_lastBin("maxMuonEnergy_bJET_vs_diff_lastBin",";max(ElectronEnergy_bJet);|genMET-pTnunu|;Events/Bin",100,0,100,100,0,100);
   TH2F maxElectronEnergy_bJET_vs_diff_lastBin("maxElectronEnergy_bJET_vs_diff_lastBin",";max(ElectronEnergy_bJet);|genMET-pTnunu|;Events/Bin",100,0,100,100,0,100);
   
   
   // ~for(std::vector<float> met_bins : {met_bins1,met_bins2,met_bins3,met_bins4}){    //Test every possible combination of the binning defined above
      // ~for(std::vector<float> phi_bins : {phi_bins1,phi_bins2,phi_bins3,phi_bins4}){
   for(std::vector<float> met_bins : {met_bins4}){    //Test every possible combination of the binning defined above
      for(std::vector<float> phi_bins : {phi_bins4}){
   
         std::vector<TH1F> met_res;    //Resolution for each bin
         std::vector<TH1F> phi_res;
         int nBinsMet=met_bins.size()-1;
         int nBinsPhi=phi_bins.size()-1;
         int nBins=nBinsMet*nBinsPhi;
         
         for(int i=0;i<nBins;i++){
            TH1F temp_met("","",200,-100,100);
            TH1F temp_phi("","",200,-3.2,3.2);
            met_res.push_back(temp_met);
            phi_res.push_back(temp_phi);
         }
         
         
         N_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         N_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         N_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         Evt_genrec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         Evt_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         Evt_rec=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         eff_gen=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         eff_gen_recSom=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         migration=hist::fromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         TopPt_profile=hist::ProfilefromWidths_2d("",";p_{T}^{#nu#nu}(GeV);|#Delta#phi|(p_{T}^{#nu#nu},nearest l);",met_bins,hist::getWidths(met_bins),phi_bins,hist::getWidths(phi_bins));
         
         // ~TString sampleName="";
         TString sampleName="dilepton";
         // ~TString sampleName="dilepton_CP5";
         // ~TString sampleName="MadGraph";
         // ~TString sampleName="T2tt_650_350";
         // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0.root","read");
         TFile file("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res100.0_new.root","read");
         TTreeReader reader((sampleName=="") ? "ttbar_res100.0/ttbar_res" : "ttbar_res100.0/ttbar_res_"+sampleName, &file);
         // ~TFile file("/net/data_cms1b/user/dmeuser/top_analysis/output/ttbar_res1.0_new.root","read");
         // ~TTreeReader reader((sampleName=="") ? "ttbar_res1.0/ttbar_res" : "ttbar_res1.0/ttbar_res_"+sampleName, &file);
         
         
         // ~TTreeReaderValue<float> MET   (reader, "MET");
         // ~TTreeReaderValue<float> MET   (reader, "DeepMET");
         TTreeReaderValue<float> MET   (reader, "PuppiMET");
         // ~TTreeReaderValue<float> MET   (reader, "genMET");
         // ~TTreeReaderValue<float> MET   (reader, "XYcorrMET");
         TTreeReaderValue<float> PtNuNu   (reader, "PtNuNu");
         // ~TTreeReaderValue<float> PtNuNu   (reader, "genMET");
         // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_rec");
         // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_recDeep");
         TTreeReaderValue<float> Phi_rec   (reader, "Phi_recPuppi");
         // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_gen");
         // ~TTreeReaderValue<float> Phi_rec   (reader, "Phi_recXYcorr");
         TTreeReaderValue<float> Phi_gen   (reader, "Phi_NuNu");
         // ~TTreeReaderValue<float> Phi_gen   (reader, "Phi_gen");
         TTreeReaderValue<float> dPhiMETnearJet   (reader, "dPhiMETnearJet");
         TTreeReaderValue<float> dPhiMETfarJet   (reader, "dPhiMETfarJet");
         TTreeReaderValue<float> dPhiMETleadJet   (reader, "dPhiMETleadJet");
         TTreeReaderValue<float> dPhiMETlead2Jet   (reader, "dPhiMETlead2Jet");
         TTreeReaderValue<float> dPhiMETbJet   (reader, "dPhiMETbJet");
         TTreeReaderValue<float> dPhiLep1Lep2   (reader, "dPhiLep1Lep2");
         TTreeReaderValue<float> METsig   (reader, "METsig");
         TTreeReaderValue<float> N   (reader, "N");
         TTreeReaderValue<UInt_t> runNo   (reader, "runNo");
         TTreeReaderValue<UInt_t> lumNo   (reader, "lumNo");
         TTreeReaderValue<ULong64_t> evtNo   (reader, "evtNo");
         TTreeReaderValue<UInt_t> genDecayMode   (reader, "genDecayMode");
         TTreeReaderValue<float> genMET   (reader, "genMET");
         TTreeReaderValue<float> PuppiMET   (reader, "PuppiMET");
         TTreeReaderValue<float> PFMET   (reader, "MET");
         TTreeReaderValue<float> Phi_recPuppi   (reader, "Phi_recPuppi");
         // ~TTreeReaderValue<float> PuppiMET   (reader, "MET");
         TTreeReaderValue<float> HT_tree   (reader, "HT");
         TTreeReaderValue<float> n_Interactions(reader, "n_Interactions");
         TTreeReaderValue<float> dPhiPtnunuMet(reader, "dPhiPtnunuMet");
         TTreeReaderValue<float> leadTop_pT(reader, "leadTop_pT");
         TTreeReaderValue<float> dPhiNuNu(reader, "dPhiNuNu");
         TTreeReaderValue<float> Lep1_flavor(reader, "Lep1_flavor");
         TTreeReaderValue<float> Lep2_flavor(reader, "Lep2_flavor");
         TTreeReaderValue<UInt_t> NnonPromptNeutrinos(reader, "NnonpromptNeutrinos");
         TTreeReaderValue<UInt_t> looseLeptonVeto(reader, "looseLeptonVeto");
         // ~TTreeReaderValue<std::vector<float>> v_bJet_muonFraction(reader, "bJet_muonFraction");
         // ~TTreeReaderValue<std::vector<float>> v_bJet_electronFraction(reader, "bJet_electronFraction");
         // ~TTreeReaderValue<std::vector<float>> v_Jet_muonFraction(reader, "Jet_muonFraction");
         // ~TTreeReaderValue<std::vector<float>> v_Jet_electronFraction(reader, "Jet_electronFraction");
         // ~TTreeReaderValue<std::vector<float>> v_bJet_muonEnergy(reader, "bJet_muonEnergy");
         // ~TTreeReaderValue<std::vector<float>> v_bJet_electronEnergy(reader, "bJet_electronEnergy");
         // ~TTreeReaderValue<std::vector<float>> v_Jet_muonEnergy(reader, "Jet_muonEnergy");
         // ~TTreeReaderValue<std::vector<float>> v_Jet_electronEnergy(reader, "Jet_electronEnergy");
         
          int migrated=0;
         
         
         
         while (reader.Next()){
            
            diffMET=*PtNuNu-*MET;//Save difference befor overflow handling
            diffPHI=*Phi_gen-*Phi_rec;
            diffMET_normed=diffMET/(*MET);
            diffPHI_normed=diffPHI/(*Phi_rec);
            truePtNuNu=*PtNuNu;
            trueDPhi_reco=*Phi_rec;
            MET_reco=*MET;
            
            if(*MET>=met_bins.back()) *MET=met_bins.back()-0.01;     //Handel overflow correctly
            if(*PtNuNu>=met_bins.back()) *PtNuNu=met_bins.back()-0.01;
            if(*Phi_rec>=phi_bins.back()) *Phi_rec=phi_bins.back()-0.01;
            if(*Phi_gen>=phi_bins.back()) *Phi_gen=phi_bins.back()-0.01;
            
            if (*PtNuNu>-1 && *Phi_gen>-1){
               eff_gen.Fill(*PtNuNu,*Phi_gen);
               if (*MET>-1 && *Phi_rec>-1) eff_gen_recSom.Fill(*PtNuNu,*Phi_gen);
            }
            
            if(*genDecayMode>3) continue;    //Remove tau events
            
            if(*genDecayMode!=3 && *PtNuNu<40) continue;   //Remove SF events if ptNuNu is smaler than 40GeV
            // ~if(*NnonPromptNeutrinos>0) continue;
            // ~if(*looseLeptonVeto==1) continue;
            
            // ~if (*MET<met_bins[0] || *PtNuNu<met_bins[0] || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
            if (*PFMET<met_bins[0] || *PtNuNu<met_bins[0] || *Phi_rec<0 || *Phi_gen<0) continue;    //Purity and stability based only on events which fullfill pseudo and reco selection
            
            /*
            std::vector<float>::iterator max_it_muon=std::max_element(v_bJet_muonFraction->begin(),v_bJet_muonFraction->end());
            std::vector<float>::iterator max_it_electron=std::max_element(v_bJet_electronFraction->begin(),v_bJet_electronFraction->end());
            float maxMuonF=(*v_bJet_muonFraction)[std::distance(v_bJet_muonFraction->begin(),max_it_muon)];
            float maxElectronF=(*v_bJet_electronFraction)[std::distance(v_bJet_electronFraction->begin(),max_it_electron)];
            // ~if(maxMuonF>0.1 || maxElectronF>0.1) continue;
            
            max_it_muon=std::max_element(v_bJet_muonEnergy->begin(),v_bJet_muonEnergy->end());
            max_it_electron=std::max_element(v_bJet_electronEnergy->begin(),v_bJet_electronEnergy->end());
            float maxMuonEnergy=(*v_bJet_muonEnergy)[std::distance(v_bJet_muonEnergy->begin(),max_it_muon)];
            float maxElectronEnergy=(*v_bJet_electronEnergy)[std::distance(v_bJet_electronEnergy->begin(),max_it_electron)];
            // ~if(maxMuonEnergy>20 || maxElectronEnergy>20) continue;
            // ~if(maxMuonEnergy>10) continue;
            */
            
            // ~if (*looseLeptonVeto==1) continue;
            // ~if(*MET<60) *MET=-1;    //Additional cuts, which might solve problem of poor dPhi resolution
            // ~if(*PtNuNu<60) *PtNuNu=-1;
            // ~if(*dPhiMETleadJet<1.0) continue;
            // ~if(*dPhiMETnearJet>0.5) continue;
            // ~if(*METsig<80.) continue;
            
            // ~if(*dPhiMETfarJet>2.64 || *dPhiMETnearJet<0.5) continue;  //Remove Events with met back to back/close with jet
            // ~if(*dPhiMETnearJet<0.5) continue;
            // ~if(*dPhiMETleadJet>2.5 || *dPhiMETlead2Jet>2.5) continue;
            
            // ~if(*dPhiMETlead2Jet<0.8) continue;
            // ~if(*dPhiMETnearJet<0.8) continue;
            // ~if(*dPhiLep1Lep2<2.1) continue;
            // ~if(*dPhiMETleadJet<1.6) continue;
            // ~if(*dPhiMETbJet>2.2) continue;
            
            float HT=*HT_tree;
            // ~if((abs(*MET-*PuppiMET)/HT)>0.05) continue;
            // ~if(abs(*MET-*PuppiMET)>30) continue;
            // ~if(*n_Interactions>35) continue;
            
            ///////////////////////////
            /////Ptmiss Scale Factor///
            ///////////////////////////
            *MET=scaleFactor * *MET;
            
            bin_gen=N_gen.Fill(*PtNuNu,*Phi_gen);
            bin_rec=N_rec.Fill(*MET,*Phi_rec);
            
            if(bin_gen==bin_rec){
               N_genrec.Fill(*PtNuNu,*Phi_gen);
               Evt_genrec.Fill(*PtNuNu,*Phi_gen,*N);
            }
            
            Evt_gen.Fill(*PtNuNu,*Phi_gen,*N);
            Evt_rec.Fill(*MET,*Phi_rec,*N);
            
            int realBin=bin_rec-((bin_rec/(nBinsMet+2)-1)*2+nBinsMet+3);   //Take into account over and underflow bin counting
            int realBin_gen=bin_gen-((bin_gen/(nBinsMet+2)-1)*2+nBinsMet+3);
            
            met_res[realBin].Fill(diffMET);
            phi_res[realBin].Fill(diffPHI);
            
            // ~METdiff.Fill(diffMET/diffMET_normed,1.+diffMET_normed);  //Same as Ptnunu/MET
            METdiff.Fill(diffMET/diffMET_normed,1.+diffMET_normed,*N);  //Same as Ptnunu/MET
            METdiff_PtNuNu.Fill(truePtNuNu,1.+diffMET_normed,*N);
            METdiff_GenMet.Fill(*genMET,1.+diffMET_normed,*N);
            METdiffpuppi_PuppiMet.Fill(*PuppiMET,truePtNuNu/(*PuppiMET),*N);
            METdiffgen.Fill(diffMET/diffMET_normed,truePtNuNu/(*genMET),*N);
            METdiff_dPhiNuLep.Fill(*Phi_gen,1.+diffMET_normed,*N);
            METdiff_dPhiMETLep.Fill(*Phi_rec,1.+diffMET_normed,*N);
            METdiffgen_dPhiMETLep.Fill(*Phi_rec,truePtNuNu/(*genMET),*N);
            
            METdiff_phi.Fill(trueDPhi_reco,1.+diffMET_normed,*N);
            METdiffgen_phi.Fill(trueDPhi_reco,truePtNuNu/(*genMET),*N);
            
            PtNuNuDiffGenMETRel_dPhiMETLep.Fill(*Phi_rec,(truePtNuNu-*genMET)/(*genMET),*N);
            GenMetDiffMETRel_dPhiMETLep.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
            PtNuNuDiffMETRel_dPhiMETLep.Fill(*Phi_rec,(truePtNuNu-MET_reco)/(*genMET),*N);
           
            TopPt_dPhi.Fill(*Phi_gen,*leadTop_pT,*N);
            
            //Save normalized METdiff for high MET as function of deltaPhi
            // ~if (MET_reco>120) {
            if (MET_reco>230) {
               GenMetDiffMETRel_dPhiMETLep_met120.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
               GenMetDiffPuppiMETRel_dPhiMETLep_met120.Fill(*Phi_recPuppi,(*genMET-*PuppiMET)/(*genMET),*N);
               METsig_dPhiMETLep_met120.Fill(*Phi_rec,*METsig,*N);
               
               switch(*genDecayMode) {
                  case(1):
                     GenMetDiffMETRel_dPhiMETLep_met120_ee.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
                     break;
                  case(2):
                     GenMetDiffMETRel_dPhiMETLep_met120_mumu.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
                     break;
                  case(3):
                     GenMetDiffMETRel_dPhiMETLep_met120_emu.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
                     break;
               }
            }
            
            //Save METdiff for different phi bins
            // ~if (realBin_gen<4) {
            if (realBin<4) {
               METdiff_phi1.Fill(diffMET/diffMET_normed,1.+diffMET_normed,*N);
               METdiffgen_phi1.Fill(diffMET/diffMET_normed,truePtNuNu/(*genMET),*N);
               dPhiMETpTnunu_phi1.Fill(diffMET/diffMET_normed,*dPhiPtnunuMet,*N);
               METdiff_dPhiLep_phi1.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_phi1.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
            }
            // ~else if (realBin_gen<8) {
            else if (realBin<8) {
               METdiff_phi2.Fill(diffMET/diffMET_normed,1.+diffMET_normed,*N);
               METdiffgen_phi2.Fill(diffMET/diffMET_normed,truePtNuNu/(*genMET),*N);
               dPhiMETpTnunu_phi2.Fill(diffMET/diffMET_normed,*dPhiPtnunuMet,*N);
               METdiff_dPhiLep_phi2.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_phi2.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
            }
            else {
               METdiff_phi3.Fill(diffMET/diffMET_normed,1.+diffMET_normed,*N);
               METdiffgen_phi3.Fill(diffMET/diffMET_normed,truePtNuNu/(*genMET),*N);
               dPhiMETpTnunu_phi3.Fill(diffMET/diffMET_normed,*dPhiPtnunuMet,*N);
               METdiff_dPhiLep_phi3.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_phi3.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
            }
            
            //Save ptTop for different met bins
            // ~if (realBin_gen%4==0) {
            if (realBin%4==0) {
               TopPt_dPhi_met1.Fill(*Phi_gen,*leadTop_pT,*N);
               METdiff_dPhiLep_met1.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_met1.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
               PtNuNuDiffGenMETRel_dPhiMETLep_met1.Fill(*Phi_rec,(truePtNuNu-*genMET)/(*genMET),*N);
               GenMetDiffMETRel_dPhiMETLep_met1.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
               PtNuNuDiffMETRel_dPhiMETLep_met1.Fill(*Phi_rec,(truePtNuNu-MET_reco)/(*genMET),*N);
            }
            // ~else if (realBin_gen%4==1) {
            else if (realBin%4==1) {
               TopPt_dPhi_met2.Fill(*Phi_gen,*leadTop_pT,*N);
               METdiff_dPhiLep_met2.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_met2.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
               PtNuNuDiffGenMETRel_dPhiMETLep_met2.Fill(*Phi_rec,(truePtNuNu-*genMET)/(*genMET),*N);
               GenMetDiffMETRel_dPhiMETLep_met2.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
               PtNuNuDiffMETRel_dPhiMETLep_met2.Fill(*Phi_rec,(truePtNuNu-MET_reco)/(*genMET),*N);
            }
            // ~else if (realBin_gen%4==2) {
            else if (realBin%4==2) {
               TopPt_dPhi_met3.Fill(*Phi_gen,*leadTop_pT,*N);
               METdiff_dPhiLep_met3.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_met3.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
               PtNuNuDiffGenMETRel_dPhiMETLep_met3.Fill(*Phi_rec,(truePtNuNu-*genMET)/(*genMET),*N);
               GenMetDiffMETRel_dPhiMETLep_met3.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
               PtNuNuDiffMETRel_dPhiMETLep_met3.Fill(*Phi_rec,(truePtNuNu-MET_reco)/(*genMET),*N);
            }
            else {
               TopPt_dPhi_met4.Fill(*Phi_gen,*leadTop_pT,*N);
               METdiff_dPhiLep_met4.Fill(*dPhiLep1Lep2,1.+diffMET_normed,*N);
               METdiff_dPhiNuNu_met4.Fill(*dPhiNuNu,1.+diffMET_normed,*N);
               PtNuNuDiffGenMETRel_dPhiMETLep_met4.Fill(*Phi_rec,(truePtNuNu-*genMET)/(*genMET),*N);
               GenMetDiffMETRel_dPhiMETLep_met4.Fill(*Phi_rec,(*genMET-MET_reco)/(*genMET),*N);
               PtNuNuDiffMETRel_dPhiMETLep_met4.Fill(*Phi_rec,(truePtNuNu-MET_reco)/(*genMET),*N);
            }
            
            // ~if (realBin==11 && realBin_gen!=11 && abs(*MET-*PtNuNu_true)>100) {
            if (realBin==11 && realBin_gen!=11) {
               migration.Fill(*PtNuNu,*Phi_gen);
               // ~std::cout<<"-------------------------------"<<std::endl;
               // ~std::cout<<*runNo<<":"<<*lumNo<<":"<<*evtNo<<std::endl;
               // ~std::cout<<"pT_NuNu="<<*PtNuNu<<"   "<<"MET="<<*MET<<"   "<<"genMET="<<*genMET<<"   "<<"puppiMET="<<*PuppiMET<<std::endl;
               // ~std::cout<<"phi_rec="<<*Phi_rec<<"   "<<"phi_gen="<<*Phi_gen<<std::endl;
               // ~std::cout<<*genDecayMode<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETleadJet)<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETlead2Jet)<<std::endl;
               // ~std::cout<<abs(3.14-*dPhiMETbJet)<<std::endl;
               // ~std::cout<<abs(*dPhiMETnearJet)<<std::endl;
               
               dPhiMETnearJet_badReso.Fill(*dPhiMETnearJet);
               dPhiMETfarJet_badReso.Fill(*dPhiMETfarJet);
               dPhiMETleadJet_badReso.Fill(*dPhiMETleadJet);
               dPhiMETlead2Jet_badReso.Fill(*dPhiMETlead2Jet);
               dPhiMETbJet_badReso.Fill(*dPhiMETbJet);
               dPhiLep1Lep2_badReso.Fill(*dPhiLep1Lep2);
               METsig_badReso.Fill(*METsig);
               diffPuppi_badReso.Fill(abs(*MET-*PuppiMET)/HT);
               migrated++;
            }
            
            else if (realBin==11 && realBin_gen==11) {
               dPhiMETnearJet_goodReso.Fill(*dPhiMETnearJet);
               dPhiMETfarJet_goodReso.Fill(*dPhiMETfarJet);
               dPhiMETleadJet_goodReso.Fill(*dPhiMETleadJet);
               dPhiMETlead2Jet_goodReso.Fill(*dPhiMETlead2Jet);
               dPhiMETbJet_goodReso.Fill(*dPhiMETbJet);
               dPhiLep1Lep2_goodReso.Fill(*dPhiLep1Lep2);
               METsig_goodReso.Fill(*METsig);
               diffPuppi_goodReso.Fill(abs(*MET-*PuppiMET)/HT);
            }
            
            //Fill histograms for unfolding
            int bin_rec_unfold=reco2D.Fill(*MET,*Phi_rec);
            int realBin_rec_unfold=bin_rec_unfold-((bin_rec_unfold/(2*nBinsMet+2)-1)*2+2*nBinsMet+3);
            trueDistributions.Fill(realBin_gen,*N);
            recoDistributions.Fill(realBin_rec_unfold,*N);
            response.Fill(realBin_rec_unfold,realBin_gen,*N);
            response_sameBins.Fill(realBin,realBin_gen);
            
            //Fill profile2d for lead top pt
            TopPt_profile.Fill(*PtNuNu,*Phi_gen,*leadTop_pT,*N);
            
            //Fill hist for last bin
            // ~if(realBin==17){
               // ~GenMETvsPTnunu_lastBin.Fill(*PtNuNu,*genMET);
               // ~maxMuonF_bJET_vs_diff_lastBin.Fill(maxMuonF,abs(*genMET-*PtNuNu));
               // ~maxElectronF_bJET_vs_diff_lastBin.Fill(maxElectronF,abs(*genMET-*PtNuNu));
               // ~maxMuonEnergy_bJET_vs_diff_lastBin.Fill(maxMuonEnergy,abs(*genMET-*PtNuNu));
               // ~maxElectronEnergy_bJET_vs_diff_lastBin.Fill(maxElectronEnergy,abs(*genMET-*PtNuNu));
               // ~if(*NnonPromptNeutrinos>0){
                  // ~GenMETvsPTnunu_lastBin_nonPromptNu.Fill(*PtNuNu,*genMET);
                  // ~dPhiReco_lastBin_nonPromptNu.Fill(*Phi_rec);
                  // ~maxMuonF_bJET_lastBin_nonPromptNu.Fill(maxMuonF);
                  // ~maxElectronF_bJET_lastBin_nonPromptNu.Fill(maxElectronF);
                  // ~maxMuonEnergy_bJET_lastBin_nonPromptNu.Fill(maxMuonEnergy);
                  // ~maxElectronEnergy_bJET_lastBin_nonPromptNu.Fill(maxElectronEnergy);
               // ~}
               // ~else {
                  // ~GenMETvsPTnunu_lastBin_PromptNu.Fill(*PtNuNu,*genMET);
                  // ~dPhiReco_lastBin_PromptNu.Fill(*Phi_rec);
                  // ~maxMuonF_bJET_lastBin_PromptNu.Fill(maxMuonF);
                  // ~maxElectronF_bJET_lastBin_PromptNu.Fill(maxElectronF);
                  // ~maxMuonEnergy_bJET_lastBin_PromptNu.Fill(maxMuonEnergy);
                  // ~maxElectronEnergy_bJET_lastBin_PromptNu.Fill(maxElectronEnergy);
               // ~}
            // ~}
            
            
         }
         file.Close();
         std::cout<<migrated<<std::endl;
         
         // ~io::RootFileSaver saver(TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100),"binningUnfolding");
         // ~io::RootFileSaver saver(TString::Format("binningUnfolding_dilepton%.1f.root",cfg.processFraction*100),"binningUnfolding");
         io::RootFileSaver saver((sampleName=="") ? TString::Format("binningUnfolding%.1f.root",cfg.processFraction*100) : TString::Format("binningUnfolding_"+sampleName+"%.1f.root",cfg.processFraction*100),"binningUnfolding");
         
         TH2F stability=Evt_genrec;
         TH2F purity=Evt_genrec;
         TH2F efficiency=eff_gen_recSom;
         TH2F efficiency_RG=N_rec;
         TH2F hist_res_phi=N_rec;
         TH2F hist_res_met=N_rec;
         
         stability.Divide(&Evt_gen);
         purity.Divide(&Evt_rec);
         efficiency.Divide(&eff_gen);
         efficiency_RG.Divide(&eff_gen);
         
         std::vector<TFitResultPtr> r_met;
         std::vector<TFitResultPtr> r_phi;
         
         
         TString Binning= (scaleFactor==1.0) ? "Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)
                                                : "Binning_Met"+std::to_string(numberBinningSchemeMet)+"_Phi"+std::to_string(numberBinningSchemePhi)+"_SF"+std::to_string(scaleFactor);
         for(int i=0;i<nBins;i++){
            r_met.push_back(met_res[i].Fit("gaus","SQ"));
            r_phi.push_back(phi_res[i].Fit("gaus","SQ"));
            saver.save(met_res[i],Binning+"/Resolution/met_bin"+std::to_string(i+1));
            saver.save(phi_res[i],Binning+"/Resolution/phi_bin"+std::to_string(i+1));
         }
         
         if (!r_met.empty() && !r_phi.empty()){
            int i=0;
            for(int y=1;y<=nBinsPhi;y++){
               for(int x=1;x<=nBinsMet;x++){
                  std::cout<<"-------------Bin Number "<<i+1<<"----Width_MET="<<hist::getWidths(met_bins)[i%nBinsMet]<<"----Width_PHI="<<hist::getWidths(phi_bins)[i%nBinsPhi]<<"-------------"<<std::endl;
                  if(r_met[i]!=-1) std::cout<<"Res_MET="<<r_met[i]->Parameter(2)<<"+-"<<r_met[i]->ParError(2)<<std::endl;
                  if(r_phi[i]!=-1) std::cout<<"Res_PHI="<<r_phi[i]->Parameter(2)<<"+-"<<r_phi[i]->ParError(2)<<std::endl;
                  std::cout<<"Purity="<<purity.GetBinContent(x,y)<<std::endl;
                  std::cout<<"Stability="<<stability.GetBinContent(x,y)<<std::endl;
                  std::cout<<"N_Events="<<met_res[i].GetEntries()<<std::endl;
                  std::cout<<"N_DataEvents="<<Evt_gen.GetBinContent(x,y)<<std::endl;
                  
                  if(r_met[i]!=-1) hist_res_phi.SetBinContent(x,y,r_phi[i]->Parameter(2));
                  if(r_phi[i]!=-1) hist_res_met.SetBinContent(x,y,r_met[i]->Parameter(2));
                  
                  i++;
               }
            }
         }
         
         std::vector<TString> z_axis={"stability","purity","efficiency","efficiency_RG","N_gen","N_rec","N_gen_rec","Evt_gen","Evt_rec","res_phi","res_met","migration"};
         int i=0;
         
         for(TH2F current_hist: {stability,purity,efficiency,efficiency_RG,N_gen,N_rec,N_genrec,Evt_gen,Evt_rec,hist_res_phi,hist_res_met,migration}){
            TCanvas can;
            can.cd();
            // ~can.SetLogz();
            
            gPad->SetRightMargin(0.2);
            gPad->SetLeftMargin(0.13);
            gPad->SetBottomMargin(0.11);
            
            current_hist.GetYaxis()->SetTitleOffset(1.3);
            current_hist.GetXaxis()->SetTitleOffset(0.9);
            // ~current_hist.GetZaxis()->SetTitleOffset(1.3);
            current_hist.GetYaxis()->SetTitleSize(0.05);
            current_hist.GetXaxis()->SetTitleSize(0.05);
            current_hist.GetZaxis()->SetTitleSize(0.05);
            current_hist.GetYaxis()->SetLabelSize(0.04);
            current_hist.GetXaxis()->SetLabelSize(0.04);
            current_hist.GetZaxis()->SetLabelSize(0.04);
            current_hist.GetZaxis()->SetLabelOffset(0.02);
            
            if (z_axis[i]=="Evt_rec") current_hist.GetZaxis()->SetTitle("expected events");
            else if (z_axis[i]=="res_phi") current_hist.GetZaxis()->SetTitle("#Delta#phi resolution");
            else if (z_axis[i]=="res_met") current_hist.GetZaxis()->SetTitle("p_{T}^{miss} resolution (GeV)");
            else current_hist.GetZaxis()->SetTitle(z_axis[i]);
            
            if (z_axis[i]=="purity" || z_axis[i]=="stability") current_hist.GetZaxis()->SetRangeUser(0.1,0.75);
            
            current_hist.SetStats(false);
            current_hist.SetMarkerColor(kRed);
            current_hist.SetMarkerSize(2.5);
            current_hist.Draw("colz text");
            
            TString plotLoc=Binning+"/"+z_axis[i];
            saver.save(can,plotLoc,true,true);
            can.Clear();
            i++;
         }
         
         std::vector<TH1F> goodHists={dPhiMETnearJet_goodReso,dPhiMETfarJet_goodReso,dPhiMETleadJet_goodReso,dPhiMETlead2Jet_goodReso,dPhiMETbJet_goodReso,dPhiLep1Lep2_goodReso,METsig_goodReso,diffPuppi_goodReso};
         std::vector<TH1F> badHists={dPhiMETnearJet_badReso,dPhiMETfarJet_badReso,dPhiMETleadJet_badReso,dPhiMETlead2Jet_badReso,dPhiMETbJet_badReso,dPhiLep1Lep2_badReso,METsig_goodReso,diffPuppi_badReso};
         std::vector<TString> names={"dPhiMETnearJet","dPhiMETfarJet","dPhiMETleadJet","dPhiMETlead2Jet","dPhiMETbJet","dPhiLep1Lep2","METsig","diffPuppi"};
         for (int i=0;i<goodHists.size();i++){
            TCanvas can;
            can.cd();
            
            // ~badHists[i].Rebin(4);
            // ~goodHists[i].Rebin(4);
            badHists[i].Rebin(5);
            goodHists[i].Rebin(5);
               
            goodHists[i].Scale(1.0/(goodHists[i].Integral()));   //Normalize the hist to the integral
            badHists[i].Scale(1.0/(badHists[i].Integral()));
            
            can.cd();
            can.SetLogz();
            
            badHists[i].GetYaxis()->SetTitle("normalized distribution");
            goodHists[i].SetLineColor(kBlue);
            badHists[i].SetLineColor(kRed);
            
            badHists[i].SetStats(0);
            goodHists[i].SetStats(0);
            
            badHists[i].Draw("hist");
            goodHists[i].Draw("same hist");
            
            badHists[i].SetMaximum(1.3*std::max(badHists[i].GetMaximum(),goodHists[i].GetMaximum()));
            
            gfx::LegendEntries legE;
            // ~legE.append(*hist_good,"good resolution (<0.3)","l");
            // ~legE.append(*hist_bad,"bad resolution (>0.3)","l");
            legE.append(goodHists[i],"true","l");
            legE.append(badHists[i],"migrated","l");
            TLegend leg=legE.buildLegend(.5,.65,0.75,.9,1);
            leg.SetTextSize(0.035);
            leg.Draw();
            
            TLatex label=gfx::cornerLabel("Last bin",1);
            // ~TLatex label=gfx::cornerLabel("MET>230 GeV, dPhi>1.4",1);
            // ~TLatex label2=gfx::cornerLabel("|#Delta#phi(met,near jet)|>0.5",2);
            label.Draw();
            // ~label2.Draw();
            
            can.RedrawAxis();
            // ~TString plotLoc=sPresel+"/"+sVar;
            TString plotLoc=Binning+"/METresolution/"+names[i];
            saver.save(can,plotLoc,true,true);
            can.Clear();
         }
         
         //Save Graph for METdiff vs met
         saveProfileFrom2D(&METdiff,"METdiff_vs_MET",&saver);
         saveProfileFrom2D(&METdiff_phi,"METdiff_vs_PHI",&saver);
         
         saveProfileFrom2D(&METdiff_PtNuNu,"METdiff_vs_PtNuNu",&saver);
         
         saveProfileFrom2D(&METdiff_GenMet,"METdiff_vs_GenMet",&saver);
         
         saveProfileFrom2D(&METdiffpuppi_PuppiMet,"METdiff_vs_PuppiMet",&saver);
         
         saveProfileFrom2D(&METdiffgen,"METdiffgen_vs_MET",&saver);
         saveProfileFrom2D(&METdiffgen_phi,"METdiffgen_vs_PHI",&saver);
         
         saveProfileFrom2D(&METdiff_phi1,"METdiff_phi1_vs_MET",&saver);
         saveProfileFrom2D(&METdiff_phi2,"METdiff_phi2_vs_MET",&saver);
         saveProfileFrom2D(&METdiff_phi3,"METdiff_phi3_vs_MET",&saver);
         
         saveProfileFrom2D(&METdiffgen_phi1,"METdiffgen_phi1_vs_MET",&saver);
         saveProfileFrom2D(&METdiffgen_phi2,"METdiffgen_phi2_vs_MET",&saver);
         saveProfileFrom2D(&METdiffgen_phi3,"METdiffgen_phi3_vs_MET",&saver);
         
         saveProfileFrom2D(&dPhiMETpTnunu_phi1,"dPhiMETpTnunu_phi1_vs_MET",&saver);
         saveProfileFrom2D(&dPhiMETpTnunu_phi2,"dPhiMETpTnunu_phi2_vs_MET",&saver);
         saveProfileFrom2D(&dPhiMETpTnunu_phi3,"dPhiMETpTnunu_phi3_vs_MET",&saver);
         
         saveProfileFrom2D(&TopPt_dPhi,"TopPt_vs_dPhi",&saver);
         saveProfileFrom2D(&TopPt_dPhi_met1,"TopPt_vs_dPhi_met1",&saver);
         saveProfileFrom2D(&TopPt_dPhi_met2,"TopPt_vs_dPhi_met2",&saver);
         saveProfileFrom2D(&TopPt_dPhi_met3,"TopPt_vs_dPhi_met3",&saver);
         saveProfileFrom2D(&TopPt_dPhi_met4,"TopPt_vs_dPhi_met4",&saver);
         
         saveProfileFrom2D(&METdiff_dPhiLep_phi1,"METdiff_dPhiLep_phi1",&saver);
         saveProfileFrom2D(&METdiff_dPhiLep_phi2,"METdiff_dPhiLep_phi2",&saver);
         saveProfileFrom2D(&METdiff_dPhiLep_phi3,"METdiff_dPhiLep_phi3",&saver);
         
         saveProfileFrom2D(&METdiff_dPhiLep_met1,"METdiff_dPhiLep_met1",&saver);
         saveProfileFrom2D(&METdiff_dPhiLep_met2,"METdiff_dPhiLep_met2",&saver);
         saveProfileFrom2D(&METdiff_dPhiLep_met3,"METdiff_dPhiLep_met3",&saver);
         saveProfileFrom2D(&METdiff_dPhiLep_met4,"METdiff_dPhiLep_met4",&saver);
         
         saveProfileFrom2D(&METdiff_dPhiNuNu_phi1,"METdiff_dPhiNuNu_phi1",&saver);
         saveProfileFrom2D(&METdiff_dPhiNuNu_phi2,"METdiff_dPhiNuNu_phi2",&saver);
         saveProfileFrom2D(&METdiff_dPhiNuNu_phi3,"METdiff_dPhiNuNu_phi3",&saver);
         
         saveProfileFrom2D(&METdiff_dPhiNuNu_met1,"METdiff_dPhiNuNu_met1",&saver);
         saveProfileFrom2D(&METdiff_dPhiNuNu_met2,"METdiff_dPhiNuNu_met2",&saver);
         saveProfileFrom2D(&METdiff_dPhiNuNu_met3,"METdiff_dPhiNuNu_met3",&saver);
         saveProfileFrom2D(&METdiff_dPhiNuNu_met4,"METdiff_dPhiNuNu_met4",&saver);
         
         saveProfileFrom2D(&METdiff_dPhiNuLep,"METdiff_dPhiNuLep",&saver);
         saveProfileFrom2D(&METdiff_dPhiMETLep,"METdiff_dPhiMETLep",&saver);
         saveProfileFrom2D(&METdiffgen_dPhiMETLep,"METdiffgen_dPhiMETLep",&saver);
         
         saveProfileFrom2D(&PtNuNuDiffGenMETRel_dPhiMETLep,"PtNuNuDiffGenMETRel_dPhiMETLep",&saver);
         saveProfileFrom2D(&PtNuNuDiffGenMETRel_dPhiMETLep_met1,"PtNuNuDiffGenMETRel_dPhiMETLep_met1",&saver);
         saveProfileFrom2D(&PtNuNuDiffGenMETRel_dPhiMETLep_met2,"PtNuNuDiffGenMETRel_dPhiMETLep_met2",&saver);
         saveProfileFrom2D(&PtNuNuDiffGenMETRel_dPhiMETLep_met3,"PtNuNuDiffGenMETRel_dPhiMETLep_met3",&saver);
         saveProfileFrom2D(&PtNuNuDiffGenMETRel_dPhiMETLep_met4,"PtNuNuDiffGenMETRel_dPhiMETLep_met4",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep,"GenMetDiffMETRel_dPhiMETLep",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met1,"GenMetDiffMETRel_dPhiMETLep_met1",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met2,"GenMetDiffMETRel_dPhiMETLep_met2",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met3,"GenMetDiffMETRel_dPhiMETLep_met3",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met4,"GenMetDiffMETRel_dPhiMETLep_met4",&saver);
         saveProfileFrom2D(&PtNuNuDiffMETRel_dPhiMETLep,"PtNuNuDiffMETRel_dPhiMETLep",&saver);
         saveProfileFrom2D(&PtNuNuDiffMETRel_dPhiMETLep_met1,"PtNuNuDiffMETRel_dPhiMETLep_met1",&saver);
         saveProfileFrom2D(&PtNuNuDiffMETRel_dPhiMETLep_met2,"PtNuNuDiffMETRel_dPhiMETLep_met2",&saver);
         saveProfileFrom2D(&PtNuNuDiffMETRel_dPhiMETLep_met3,"PtNuNuDiffMETRel_dPhiMETLep_met3",&saver);
         saveProfileFrom2D(&PtNuNuDiffMETRel_dPhiMETLep_met4,"PtNuNuDiffMETRel_dPhiMETLep_met4",&saver);
         
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met120,"GenMetDiffMETRel_dPhiMETLep_met120",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met120_ee,"GenMetDiffMETRel_dPhiMETLep_met120_ee",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met120_emu,"GenMetDiffMETRel_dPhiMETLep_met120_emu",&saver);
         saveProfileFrom2D(&GenMetDiffMETRel_dPhiMETLep_met120_mumu,"GenMetDiffMETRel_dPhiMETLep_met120_mumu",&saver);
         
         saveProfileFrom2D(&GenMetDiffPuppiMETRel_dPhiMETLep_met120,"GenMetDiffPuppiMETRel_dPhiMETLep_met120",&saver);
         
         saveProfileFrom2D(&METsig_dPhiMETLep_met120,"METsig_dPhiMETLep_met120",&saver);
         
         saver.save(TopPt_profile,"TopPt_profile");
         
         //Save histograms for unfolding
         for (int i=0; i<trueDistributions.GetNbinsX(); i++){
            trueDistributions.SetBinError(i,sqrtf(trueDistributions.GetBinContent(i)));
         }
         saver.save(trueDistributions,Binning+"/Unfolding/trueDistribution");
         saver.save(recoDistributions,Binning+"/Unfolding/recoDistribution");
         saver.save(response,Binning+"/Unfolding/response");
         saver.save(response_sameBins,Binning+"/Unfolding/response_sameBins");
         saver.save(reco2D,Binning+"/Unfolding/reco2D");
         
         //Save hists for last bin
         saver.save(GenMETvsPTnunu_lastBin,"lastBin/GenMETvsPTnunu");
         saver.save(GenMETvsPTnunu_lastBin_nonPromptNu,"lastBin/GenMETvsPTnunu_nonPromptNu");
         saver.save(GenMETvsPTnunu_lastBin_PromptNu,"lastBin/GenMETvsPTnunu_PromptNu");
         dPhiReco_lastBin_PromptNu.Scale(1./dPhiReco_lastBin_PromptNu.Integral());
         dPhiReco_lastBin_nonPromptNu.Scale(1./dPhiReco_lastBin_nonPromptNu.Integral());
         saver.save(dPhiReco_lastBin_PromptNu,"lastBin/dPhiReco_lastBin_PromptNu");
         saver.save(dPhiReco_lastBin_nonPromptNu,"lastBin/dPhiReco_lastBin_nonPromptNu");
         
         maxMuonF_bJET_lastBin_nonPromptNu.Scale(1./maxMuonF_bJET_lastBin_nonPromptNu.Integral());
         maxMuonF_bJET_lastBin_PromptNu.Scale(1./maxMuonF_bJET_lastBin_PromptNu.Integral());
         saver.save(maxMuonF_bJET_lastBin_nonPromptNu,"lastBin/maxMuonF_bJET_lastBin_nonPromptNu");
         saver.save(maxMuonF_bJET_lastBin_PromptNu,"lastBin/maxMuonF_bJET_lastBin_PromptNu");
         maxElectronF_bJET_lastBin_nonPromptNu.Scale(1./maxElectronF_bJET_lastBin_nonPromptNu.Integral());
         maxElectronF_bJET_lastBin_PromptNu.Scale(1./maxElectronF_bJET_lastBin_PromptNu.Integral());
         saver.save(maxElectronF_bJET_lastBin_nonPromptNu,"lastBin/maxElectronF_bJET_lastBin_nonPromptNu");
         saver.save(maxElectronF_bJET_lastBin_PromptNu,"lastBin/maxElectronF_bJET_lastBin_PromptNu");
         saver.save(maxElectronF_bJET_vs_diff_lastBin,"lastBin/maxElectronF_bJET_vs_diff_lastBin");
         saver.save(maxMuonF_bJET_vs_diff_lastBin,"lastBin/maxMuonF_bJET_vs_diff_lastBin");
         
         // ~maxMuonEnergy_bJET_lastBin_nonPromptNu.Scale(1./maxMuonEnergy_bJET_lastBin_nonPromptNu.Integral());
         // ~maxMuonEnergy_bJET_lastBin_PromptNu.Scale(1./maxMuonEnergy_bJET_lastBin_PromptNu.Integral());
         saver.save(maxMuonEnergy_bJET_lastBin_nonPromptNu,"lastBin/maxMuonEnergy_bJET_lastBin_nonPromptNu");
         saver.save(maxMuonEnergy_bJET_lastBin_PromptNu,"lastBin/maxMuonEnergy_bJET_lastBin_PromptNu");
         // ~maxElectronEnergy_bJET_lastBin_nonPromptNu.Scale(1./maxElectronEnergy_bJET_lastBin_nonPromptNu.Integral());
         // ~maxElectronEnergy_bJET_lastBin_PromptNu.Scale(1./maxElectronEnergy_bJET_lastBin_PromptNu.Integral());
         saver.save(maxElectronEnergy_bJET_lastBin_nonPromptNu,"lastBin/maxElectronEnergy_bJET_lastBin_nonPromptNu");
         saver.save(maxElectronEnergy_bJET_lastBin_PromptNu,"lastBin/maxElectronEnergy_bJET_lastBin_PromptNu");
         saver.save(maxElectronEnergy_bJET_vs_diff_lastBin,"lastBin/maxElectronEnergy_bJET_vs_diff_lastBin");
         saver.save(maxMuonEnergy_bJET_vs_diff_lastBin,"lastBin/maxMuonEnergy_bJET_vs_diff_lastBin");
         
         numberBinningSchemePhi++;
      }
      
      numberBinningSchemeMet++;
      numberBinningSchemePhi=1;
   }
   
}
