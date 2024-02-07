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
    string year = string(cfg.year);
    // Location of the ROOT histogram files
    // ~string histLoc = "/net/data_cms1b/user/dmeuser/top_analysis/" + string("2018") + "/v06/output_framework/multiHists/";
    string histLoc = "";
    if (year=="2016_preVFP" or year=="2016_postVFP"){
        // ~histLoc+="/net/data_cms1b/user/dmeuser/top_analysis/" + year + "/"+string(cfg.treeVersion.Data())+"/output_framework/multiHists/combine_new/";
        histLoc+="/net/data_cms1b/user/dmeuser/top_analysis/" + year + "/"+string(cfg.treeVersion.Data())+"/output_framework/multiHists/combine/";
    } else{
        histLoc+="/net/data_cms1b/user/dmeuser/top_analysis/" + year + "/"+string(cfg.treeVersion.Data())+"/output_framework/multiHists/combine/";
    }
    
    map<string,tuple<vector<double>,string>> test_variables;
    
    // variables not used in DNN:
    // ~test_variables["Lep1_pt*cos(Lep1_phi)"] = tuple(linspace(-250, 250, 26), "Lep1_pX");
    // ~test_variables["Lep1_phi"] = tuple(linspace(-3.2, 3.2, 26), "Lep1_phi");
    // ~test_variables["Lep1_pt"] = tuple(linspace(0, 420, int(420/30+1)), "Lep1_pt");
    // ~test_variables["vecsum_pT_allJet"] = tuple(linspace(0, 500, 26), "vecsum_pT_allJet");
    // ~test_variables["METunc_Puppi"] = tuple(linspace(0, 40, 21, 0., 100.), "METunc_Puppi");
    // ~test_variables["CaloMET*cos(CaloMET)"] = tuple(linspace(-150, 150, 31, -250., 250.), "CaloMET_X");
    // ~test_variables["CaloMET*sind(CaloMET)"] = tuple(linspace(-150, 150, 31, -250., 250.), "CaloMET_Y");
    // ~test_variables["n_Interactions"] = tuple(linspace(0, 80, 41, 0., 100.), "n_Interactions");
    // ~test_variables["dPhiMETleadJet_Puppi"] = tuple(linspace(0, 3.2, 33), "dPhiMETleadJet_Puppi");
    // ~test_variables["Lep2_pt*cos(Lep2_phi)"] = tuple(linspace(-250, 250, 26), "Lep2_pX");
    // ~test_variables["Lep2_pt*sin(Lep2_phi)"] = tuple(linspace(-250, 250, 26), "Lep2_pY");
    // ~test_variables["PuppiMET*cos(PuppiMET_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "PuppiMET_X");
    // ~test_variables["PuppiMET*sin(PuppiMET_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "PuppiMET_Y");
    // ~test_variables["MET*cos(PFMET_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "MET_X");
    // ~test_variables["MET*sin(PFMET_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "MET_Y");
    
    // DNN input variables:
    
    // ~test_variables["PuppiMET_xy*cos(PuppiMET_xy_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "PuppiMET_xy_X");
    // ~test_variables["PuppiMET_xy*sin(PuppiMET_xy_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "PuppiMET_xy_Y");
    // ~test_variables["MET_xy*cos(MET_xy_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "MET_xy_X");
    // ~test_variables["MET_xy*sin(MET_xy_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "MET_xy_Y");
    // ~test_variables["vecsum_pT_allJet*cos(HT_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "vecsum_pT_allJet_X");
    // ~test_variables["vecsum_pT_allJet*sin(HT_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "vecsum_pT_allJet_Y");
    // ~test_variables["mass_l1l2_allJet"] = tuple(linspace(120, 1800, 29, 0.), "mass_l1l2_allJet");
    // ~test_variables["Jet1_pt*sin(Jet1_phi)"] = tuple(linspace(-400, 400, 41), "Jet1_pY");
    // ~test_variables["MHT"] = tuple(linspace(60, 1800, 30, 0., 2000.), "MHT");
    // ~test_variables["Lep1_pt*cos(Lep1_phi)"] = tuple(linspace(-250, 250, 26), "Lep1_pX");
    // ~test_variables["Lep1_pt*sin(Lep1_phi)"] = tuple(linspace(-250, 250, 26), "Lep1_pY");
    // ~test_variables["Jet1_pt*cos(Jet1_phi)"] = tuple(linspace(-400, 400, 41), "Jet1_pX");
    // ~test_variables["CaloMET"] = tuple(linspace(0, 200, 41, 0., 500.), "CaloMET");
    // ~test_variables["MT2"] = tuple(linspace(0, 160, 41, 0., 200.), "MT2");
    // ~test_variables["mjj"] = tuple(linspace(60, 1800, 30, 0., 2000.), "mjj");
    // ~test_variables["nJets"] = tuple(linspace(1.5, 8.5, 8, 1.5, 12.5), "nJets");
    // ~test_variables["Jet1_E"] = tuple(linspace(40, 500, 24, 0.), "Jet1_E");
    // ~test_variables["HT"] = tuple(linspace(50, 1200, 24, 0., 2500.), "HT");
    // ~test_variables["Jet2_pt*cos(Jet2_phi)"] = tuple(linspace(-250, 250, 26), "Jet2_pX");
    // ~test_variables["Jet2_pt*sin(Jet2_phi)"] = tuple(linspace(-250, 250, 26), "Jet2_pY");
    // ~test_variables["DeepMET_reso*cos(DeepMET_reso_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "DeepMET_reso_X");
    // ~test_variables["DeepMET_reso*sin(DeepMET_reso_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "DeepMET_reso_Y");
    // ~test_variables["DeepMET_resp*cos(DeepMET_resp_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "DeepMET_resp_X");
    // ~test_variables["DeepMET_resp*sin(DeepMET_resp_phi)"] = tuple(linspace(-150, 150, 31, -250., 250.), "DeepMET_resp_Y");
    
    std::vector<double> metBinning = {0,20,40,54,68,84,100,120,140,168,196,228,260,296,332,371,410,455,500};
    test_variables["DNN_MET_pT"] = tuple(metBinning, "DNN_MET_pT");
    // ~test_variables["DNN_MET_pT"] = tuple(linspace(0., 500., 25 , 0., 500.), "DNN_MET_pT");


    map<string,tuple<vector<double>,string>>::iterator it;
    
    for(it=test_variables.begin(); it!=test_variables.end(); ++it){
        // Create an empty CombineHarvester instance that will hold all of the
        // datacard configuration and histograms etc.
        ch::CombineHarvester cb;
        // Uncomment this next line to see a *lot* of debug information
        // cb.SetVerbosity(3);
        
        string test_variable = it->first;
        std::vector<double> rebinVec = get<0>(it->second);
        string testVarName = get<1>(it->second);

        // Define semileptonic channels as categories
        ch::Categories cats = {
            {1, "mumu"},
            {2, "ee"},
            {3, "emu"}
        };
        
        // Mass declaration needed by the program; Not used in any calculation
        vector<string> masses = {"125"};
        
        cb.AddObservations({"*"}, {"ttbar"}, {year}, {""}, cats);
        
        //Add backgrounds
        // ~vector<string> bkg_procs = {"DrellYan_comb", "TTbar_diLepton_tau", "TTbar_singleLepton", "TTbar_hadronic", "SingleTop", "WJetsToLNu", "WW", "WZ", "ZZ", "ttZ_QQ", "ttZ_2L", "ttW"};
        vector<string> bkg_procs = {"DrellYan_comb", "TTbar_other", "SingleTop", "otherBKG"};
        cb.AddProcesses({"*"}, {"ttbar"}, {year}, {""}, bkg_procs, cats, false);
        
        //Add signals
        // ~vector<string> sig_procs = {"TTbar_diLepton"};
        vector<string> sig_procs = {"TTbar_diLepton"};
        cb.AddProcesses({"*"}, {"ttbar"}, {year}, {""}, sig_procs, cats, true);

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
        vector<string> shapeUncAllMC;
        shapeUncAllMC = {"BTAGBC_CORR", "BTAGBC_UNCORR", "BTAGL_CORR", "BTAGL_UNCORR", "JEREta0", "JEREta1", "JESAbsolute", "JESAbsoluteYear", "JESBBEC1", "JESBBEC1Year", "JESFlavorRealistic", "JESRelativeBalreg", "JESRelativeSampleYear", "JETPILEUPID", "PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "TOP_PT", "L1PREFIRING", "MESCALE", "MEFACSCALE", "MERENSCALE"};
        // ~shapeUncAllMC = {"BTAGBC", "BTAGL", "JER", "JESAbsolute", "JESAbsoluteYear", "JESBBEC1", "JESBBEC1Year", "JESFlavorRealistic", "JESRelativeBalreg", "JESRelativeSampleYear", "JETPILEUPID","PDF_ALPHAS", "PSFSRSCALE", "PSISRSCALE", "PU", "L1PREFIRING", "MESCALE", "MEFACSCALE", "MERENSCALE"};
        for (auto systNameAll : shapeUncAllMC) {
            cb.cp().process(ch::JoinStr({sig_procs, bkg_procs})).AddSyst(cb, systNameAll, "shape", SystMap<>::init(1.00));
        }
        
        for (string eeUncAllMC : {"ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "UNCLUSTERED"}){
            cb.cp().bin({"ee"}).AddSyst(cb, eeUncAllMC, "shape", SystMap<>::init(1.00));
        }
        for (string emuUncAllMC : {"ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "MUON_ID_STAT", "MUON_ID_SYST", "MUON_ISO_STAT", "MUON_ISO_SYST", "MUON_SCALE"}){
        // ~for (string emuUncAllMC : {"ELECTRON_ID", "ELECTRON_RECO", "ELECTRON_SCALESMEARING", "MUON_ID", "MUON_ISO", "MUON_SCALE"}){
            cb.cp().bin({"emu"}).AddSyst(cb, emuUncAllMC, "shape", SystMap<>::init(1.00));
        }
        for (string mumuUncAllMC : {"MUON_ID_STAT", "MUON_ID_SYST", "MUON_ISO_STAT", "MUON_ISO_SYST", "MUON_SCALE", "UNCLUSTERED"}){
        // ~for (string mumuUncAllMC : {"MUON_ID", "MUON_ISO", "MUON_SCALE", "UNCLUSTERED"}){
            cb.cp().bin({"mumu"}).AddSyst(cb, mumuUncAllMC, "shape", SystMap<>::init(1.00));
        }
        for (string bin_name : {"ee", "mumu", "emu"}){
            cb.cp().bin({bin_name}).AddSyst(cb, "TRIG_$BIN", "shape", SystMap<>::init(1.00));
        }
        
        
        // Uncertainties, that are not applied to DY processes:
        // ~vector<string> shapeUncNoDY = {"MEFACSCALE", "MERENSCALE", "MESCALE"};
        // ~for (auto systNameTTbar : shapeUncNoDY) {
            // ~cb.cp().process({"TTbar_diLepton", "TTbar_other", "SingleTop", "otherBKG"}).AddSyst(cb, systNameTTbar, "shape", SystMap<>::init(1.00));
        // ~}
        
        
        // Uncertainties, that are applied only to ttbar processes:
        vector<string> shapeUncTTbarOnly = {"BSEMILEP", "CR1", "CR2", "ERDON", "MTOP", "MATCH", "TOP_PT", "UETUNE"};
        
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
            histLoc + "combineInput_"+string(cfg.treeVersion.Data())+".root",
            // ~histLoc + "combine/combineInput_v06_old.root",
            "distributions100.0/baseline/$BIN/"+string(test_variable)+"/$PROCESS",
            "distributions100.0/baseline/$BIN/"+string(test_variable)+"/$PROCESS_$SYSTEMATIC"
            );
        
        std::cout << "\nrebin " << testVarName << " (beginning/end/binwidth): " << to_string(rebinVec[0]) << "/" << to_string(rebinVec.back()) << "/" << to_string(rebinVec[1]-rebinVec[0]) << endl;
        cb.VariableRebin(rebinVec); //pt

        
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

