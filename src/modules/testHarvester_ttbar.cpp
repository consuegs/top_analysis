#include <string>
#include <map>
#include <set>
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib>
#include <boost/filesystem.hpp>

#include <CombineHarvester/CombineTools/interface/CombineHarvester.h>
#include <CombineHarvester/CombineTools/interface/Observation.h>
#include <CombineHarvester/CombineTools/interface/Process.h>
#include <CombineHarvester/CombineTools/interface/Utilities.h>
#include <CombineHarvester/CombineTools/interface/Systematics.h>
#include <CombineHarvester/CombineTools/interface/BinByBin.h>

#include "Config.hpp"

std::vector<double> linspace(double a, double b, std::size_t N){
        double h = (b - a) / static_cast<double>(N-1);
        std::vector<double> xs(N);
        std::vector<double>::iterator x;
        double val;
        for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
            *x = val;
        }
        return xs;
}

using namespace std;
using ch::syst::SystMap;
using ch::syst::SystMapAsymm;
using ch::syst::era;
using ch::syst::bin_id;
using ch::syst::process;
Config const &cfg=Config::get();

extern "C"
void run(){
  
  // Location of the ROOT histogram files
  string histLoc = "/net/data_cms1b/user/dmeuser/top_analysis/" + string("2018") + "/v06/output_framework/multiHists/";

  // ~vector<string> test_variables = {"nJets", "Lep1_pt"};
  // ~vector<string> test_variables = {"Lep1_pt"};
  // ~vector<string> test_variables = {"n_Interactions"};
  // ~vector<string> test_variables = {"Lep1_phi"};
  vector<string> test_variables = {"Lep1_pt*cos(Lep1_phi)"};
  
  for (auto test_variable : test_variables) {
    // Create an empty CombineHarvester instance that will hold all of the
    // datacard configuration and histograms etc.
    ch::CombineHarvester cb;
    // Uncomment this next line to see a *lot* of debug information
    // cb.SetVerbosity(3);

    // Define semileptonic channels as categories
    ch::Categories cats = {
        {1, "mumu"},
        {2, "ee"},
        {3, "emu"}
      };
      
    
    // Mass declaration needed by the program; Not used in any calculation
    vector<string> masses = {"120"};
    
    cb.AddObservations({"*"}, {"ttbar"}, {"2018"}, {""}, cats);
    
    //Add backgrounds
    // ~vector<string> bkg_procs = {"DrellYan_comb", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic", "SingleTop", "WJetsToLNu", "WW", "WZ", "ZZ", "ttZ_QQ", "ttZ_2L", "ttW"};
    vector<string> bkg_procs = {"DrellYan_comb", "TTbar_other", "SingleTop", "otherBKG"};
    cb.AddProcesses({"*"}, {"ttbar"}, {"2018"}, {""}, bkg_procs, cats, false);
    
    //Add signals
    // ~vector<string> sig_procs = {"TTbar_diLepton"};
    vector<string> sig_procs = {"TTbar_diLepton"};
    cb.AddProcesses({"*"}, {"ttbar"}, {"2018"}, {""}, sig_procs, cats, true);

    // Read in Uncertainties:
    // Lumi:
    float lumiUnc = cfg.systUncFactor.at("LUMI").first;
    cb.cp().process(ch::JoinStr({sig_procs, bkg_procs})).AddSyst(cb, "lumi", "lnN", SystMap<>::init(1+lumiUnc));
    
    // Add xsec unceratinties
    cb.cp().process({"DrellYan_comb"}).AddSyst(cb, "DY_xsec", "lnN", SystMap<>::init(1+cfg.systUncFactor.at("XSEC_DY").first));
    cb.cp().process({"TTbar_other"}).AddSyst(cb, "TTother_xsec", "lnN", SystMap<>::init(1+cfg.systUncFactor.at("XSEC_TTOTHER").first));
    cb.cp().process({"SingleTop"}).AddSyst(cb, "ST_xsec", "lnN", SystMap<>::init(1+cfg.systUncFactor.at("XSEC_ST").first));
    cb.cp().process({"otherBKG"}).AddSyst(cb, "other_xsec", "lnN", SystMap<>::init(1+cfg.systUncFactor.at("XSEC_OTHER").first));
    
    // Uncertainties, that are applied to all processes:
    // ~vector<string> shapeUncAllMC = {"BTAGBC", "BTAGL", "ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "JER", "JESTotal", "MEFACSCALE", "MERENSCALE", "MUON_ID", "MUON_ISO", "MUON_SCALE", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "TRIG", "UNCLUSTERED"}; 
    vector<string> shapeUncAllMC = {"BTAGBC", "BTAGL", "ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "JER", "JESTotal", "MUON_ID", "MUON_ISO", "MUON_SCALE", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "TRIG", "UNCLUSTERED"}; 
    for (auto systNameAll : shapeUncAllMC) {
      cb.cp().process(ch::JoinStr({sig_procs, bkg_procs})).AddSyst(cb, systNameAll, "shape", SystMap<>::init(1.00));
    }
    
    // Uncertainties, that are not applied to DY processes:
    vector<string> shapeUncNoDY = {"MEFACSCALE", "MERENSCALE"};
    for (auto systNameTTbar : shapeUncNoDY) {
      cb.cp().process({"TTbar_diLepton", "TTbar_other", "SingleTop", "otherBKG"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
      // ~cb.cp().process({"TTbar_diLepton", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
    }
    
    // Uncertainties, that are applied only to ttbar processes:
    vector<string> shapeUncTTbarOnly = {"BFRAG", "BSEMILEP", "CR1", "CR2", "ERDON", "MTOP"};
    for (auto systNameTTbar : shapeUncTTbarOnly) {
      cb.cp().process({"TTbar_diLepton", "TTbar_other"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
      // ~cb.cp().process({"TTbar_diLepton", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
    }
    
    // "MATCH" and "UETUNE" have counts of 0.0 in ttbar_hadronic for nJets in the ee channel; prevents combine from functioning
    cb.cp().bin({"emu", "mumu"}).process({"TTbar_diLepton", "TTbar_other"}).AddSyst(cb, "MATCH", "shape", SystMap<>::init(1.00));
    cb.cp().bin({"emu", "mumu"}).process({"TTbar_diLepton", "TTbar_other"}).AddSyst(cb, "UETUNE", "shape", SystMap<>::init(1.00));
    // ~cb.cp().bin({"emu", "mumu"}).process({"TTbar_diLepton", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic"}).AddSyst(cb, "MATCH", "shape", SystMap<>::init(1.00));
    // ~cb.cp().bin({"emu", "mumu"}).process({"TTbar_diLepton", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic"}).AddSyst(cb, "UETUNE", "shape", SystMap<>::init(1.00));
    
    
    // Extract shapes and create datacards for all variables in test_variables
    cb.cp().ExtractShapes(
        histLoc + "combine/combineInput_v06.root",
        // ~histLoc + "combine/combineInput_v06_old.root",
        "distributions100.0/baseline/$BIN/"+string(test_variable)+"/$PROCESS",
        "distributions100.0/baseline/$BIN/"+string(test_variable)+"/$PROCESS_$SYSTEMATIC"
        );
    
    std::cout << "rebin: " << endl;
    // ~cb.VariableRebin(linspace(-3.2, 3.2, 51));
    // ~cb.VariableRebin(linspace(-3.2, 3.2, 26));
    // ~cb.VariableRebin(linspace(0, 80, 41));
    // ~cb.VariableRebin(linspace(0, 420, int(420/30+1))); //pt
    cb.VariableRebin(linspace(-250, 250, 26)); //pt
    
    ch::SetStandardBinNames(cb, "$BIN");

    set<string> bins = cb.bin_set();
    
    cb.AddDatacardLineAtEnd("* autoMCStats 10 1 1");
    
    if (test_variable=="Lep1_pt*cos(Lep1_phi)") test_variable ="Lep1_pX";
    
    boost::filesystem::create_directories((cfg.outputDirectory+"/datacards/").Data());
    TFile output(string(cfg.outputDirectory+"/datacards/ttbar_"+test_variable+"_1d.input.root").c_str(), "RECREATE");
    for (auto m : masses) {
      cout << ">> Writing datacard for variable: " << test_variable << "\n";
      cb.cp().mass({m, "*"}).WriteDatacard((cfg.outputDirectory+"/datacards/ttbar_" + test_variable + "_1d.txt").Data(), output);
    }
  }
}

