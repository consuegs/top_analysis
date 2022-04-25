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

#include <boost/filesystem.hpp>

std::vector<double> linspace(double a, double b, std::size_t N, double f = -20000, double e = -20000){
    double h = (b - a) / static_cast<double>(N-1);
    std::vector<double> xs(N);
    std::vector<double>::iterator x;
    double val;
    for (x = xs.begin(), val = a; x != xs.end(); ++x, val += h) {
        *x = val;
    }
    if (f > -20000){
        xs.front() = f;
    }
    if (e > -20000){
        xs.back() = e;
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
    // ~string histLoc = "/net/data_cms1b/user/dmeuser/top_analysis/" + string("2018") + "/v06/output_framework/multiHists/";
    string histLoc = "/net/data_cms1b/user/dmeuser/top_analysis/" + string("2018") + "/v07/output_framework/multiHists/";
    
    vector<tuple<string,string>> DNN_2Dhists_vec;
    
    DNN_2Dhists_vec.push_back(tuple("PuppiMET_xy*cos(PuppiMET_xy_phi)", "PuppiMET_xy_X"));
    DNN_2Dhists_vec.push_back(tuple("PuppiMET_xy*sin(PuppiMET_xy_phi)", "PuppiMET_xy_Y"));
    DNN_2Dhists_vec.push_back(tuple("MET_xy*cos(MET_xy_phi)", "MET_xy_X"));
    DNN_2Dhists_vec.push_back(tuple("MET_xy*sin(MET_xy_phi)", "MET_xy_Y"));
    DNN_2Dhists_vec.push_back(tuple("vecsum_pT_allJet*cos(HT_phi)", "vecsum_pT_allJet_X"));
    DNN_2Dhists_vec.push_back(tuple("vecsum_pT_allJet*sin(HT_phi)", "vecsum_pT_allJet_Y"));
    DNN_2Dhists_vec.push_back(tuple("mass_l1l2_allJet", "mass_l1l2_allJet"));
    DNN_2Dhists_vec.push_back(tuple("Jet1_pt*sin(Jet1_phi)", "Jet1_pY"));
    DNN_2Dhists_vec.push_back(tuple("MHT", "MHT"));
    DNN_2Dhists_vec.push_back(tuple("Lep1_pt*cos(Lep1_phi)", "Lep1_pX"));
    DNN_2Dhists_vec.push_back(tuple("Lep1_pt*sin(Lep1_phi)", "Lep1_pY"));
    DNN_2Dhists_vec.push_back(tuple("Jet1_pt*cos(Jet1_phi)", "Jet1_pX"));
    DNN_2Dhists_vec.push_back(tuple("CaloMET", "CaloMET"));
    DNN_2Dhists_vec.push_back(tuple("MT2", "MT2"));
    DNN_2Dhists_vec.push_back(tuple("mjj", "mjj"));
    DNN_2Dhists_vec.push_back(tuple("nJets", "nJets"));
    DNN_2Dhists_vec.push_back(tuple("Jet1_E", "Jet1_E"));
    DNN_2Dhists_vec.push_back(tuple("HT", "HT"));
    DNN_2Dhists_vec.push_back(tuple("Jet2_pt*cos(Jet2_phi)", "Jet2_pX"));
    DNN_2Dhists_vec.push_back(tuple("Jet2_pt*sin(Jet2_phi)", "Jet2_pY"));
    
    
    for(int i=0; i<DNN_2Dhists_vec.size(); ++i){
        for(int j=i+1; j<DNN_2Dhists_vec.size(); ++j){ 
            string short_n1 = get<1>(DNN_2Dhists_vec.at(i));
            string short_n2 = get<1>(DNN_2Dhists_vec.at(j));
            string long_n1 = get<0>(DNN_2Dhists_vec.at(i));
            string long_n2 = get<0>(DNN_2Dhists_vec.at(j));
            // ~cout << n1 << "   " << n2 << endl;
            string testVarName = short_n1+"_VS_"+short_n2;
            string test_variable = long_n1+"_VS_"+long_n2;
            
            
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
            vector<string> masses = {"125"};
            
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
            // ~vector<string> shapeUncAllMC = {"BTAGBC", "BTAGL", "ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "JER", "JESTotal", "MUON_ID", "MUON_ISO", "MUON_SCALE", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "TRIG"}; 
            vector<string> shapeUncAllMC = {"BTAGBC", "BTAGL", "ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "JER", "JESAbsoluteMPFBias", "JESAbsoluteScale", "JESAbsoluteStat", "JESFlavorPureCharm", "JESFlavorPureBottom", "JESFlavorPureGluon", "JESFlavorPureQuark", "JESFlavorQCD", "JESFlavorRealistic", "JESFragmentation", "JESPileUpDataMC", "JESPileUpPtBB", "JESPileUpPtEC1", "JESPileUpPtRef", "JESRelativeBal", "JESRelativeFSR", "JESRelativeJEREC1", "JESRelativePtBB", "JESRelativePtEC1", "JESRelativeSample", "JESRelativeStatEC", "JESRelativeStatFSR", "JESSinglePionECAL", "JESSinglePionHCAL", "JESTimePtEta", "MUON_ID", "MUON_ISO", "MUON_SCALE", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "TRIG", "UNCLUSTERED", "L1PREFIRING"};                     
            // ~vector<string> shapeUncAllMC = {"BTAGBC", "BTAGL", "ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "JER", "JESTotal", "MUON_ID", "MUON_ISO", "MUON_SCALE", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "TRIG", "UNCLUSTERED", "L1PREFIRING"}; 
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
            // ~vector<string> shapeUncTTbarOnly = {"BFRAG", "BSEMILEP", "CR1", "CR2", "ERDON", "MTOP"};
            vector<string> shapeUncTTbarOnly = {"BFRAG", "BSEMILEP", "CR1", "CR2", "ERDON", "MTOP", "MATCH", "UETUNE"};
            for (auto systNameTTbar : shapeUncTTbarOnly) {
                cb.cp().process({"TTbar_diLepton", "TTbar_other"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
                // ~cb.cp().process({"TTbar_diLepton", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
            }
            
            for (int i=1; i<51; i++) {
                string systName = "PDF_"+std::to_string(i);
                cb.cp().process({"TTbar_diLepton"}).AddSyst(cb, systName, "shape", SystMap<>::init(1.00));
            }
            
            
            // Extract shapes and create datacards for all variables in test_variables
            cout << "Not printing \"Warning: process shape has negative bins\" in CombineHarvester/CombineTools/src/CombineHarvester.cc ll.406 and ll.551" << endl;
            cb.cp().ExtractShapes(
                histLoc + "combine/combineInput_v07.root",
                // ~histLoc + "combine/combineInput_v06_old.root",
                "distributions100.0/baseline_GOF2D/$BIN/"+test_variable+"/$PROCESS",
                "distributions100.0/baseline_GOF2D/$BIN/"+test_variable+"/$PROCESS_$SYSTEMATIC"
                );
            
            // ~std::cout << "\nrebin " << testVarName << " (beginning/end/binwidth): " << to_string(rebinVec[0]) << "/" << to_string(rebinVec.back()) << "/" << to_string(rebinVec[1]-rebinVec[0]) << endl;
            // ~cb.VariableRebin(rebinVec); //pt

            
            ch::SetStandardBinNames(cb, "$BIN");

            set<string> bins = cb.bin_set();
            
            cb.AddDatacardLineAtEnd("* autoMCStats 10 1 1");
            
            boost::filesystem::create_directories((cfg.outputDirectory+"/datacards/"+testVarName+"/").Data());
            TFile output(string(cfg.outputDirectory+"/datacards/"+testVarName+"/ttbar_"+testVarName+".input.root").c_str(), "RECREATE");
            for (auto m : masses) {
                cout << ">> Writing datacard for variable: " << testVarName << "\n";
                cb.cp().mass({m, "*"}).WriteDatacard((cfg.outputDirectory+"/datacards/"+testVarName+"/ttbar_" + testVarName + ".txt").Data(), output);
            }
        }
    }
}

