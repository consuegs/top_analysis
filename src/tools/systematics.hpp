//Mainly taken from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/common/include/sampleHelpers.h
#ifndef SYTEMATICS_HPP__
#define SYTEMATICS_HPP__

#include <vector>
#include <string>

#include <TString.h>

/// Namespace to treat systematics as enum types
namespace Systematic{

    /// All systematic types as needed in any part of the framework
    enum Type{
        nominal,            // nominal, i.e. no systematic variation applied
        mH110,              // Higgs mass of 110 GeV
        mH115,              // Higgs mass of 115 GeV
        mH120,              // Higgs mass of 120 GeV
        mH1225,             // Higgs mass of 122.5 GeV
        mH1275,             // Higgs mass of 127.5 GeV
        mH130,              // Higgs mass of 130 GeV
        mH135,              // Higgs mass of 135 GeV
        mH140,              // Higgs mass of 140 GeV
        mTop169p5,          // top quark mass of 169.5 GeV
        mTop175p5,          // top quark mass of 175.5 GeV
        mtop,               // top quark mass scaled by +/- 1 GeV using samples with 169.5 and 175.5 GeV combining all ttbar processes
        mtop_ind,           // top quark mass scaled by +/- 1 GeV using samples with 169.5 and 175.5 GeV individually for ttbar processes
        lept,               // scale lepton ID/ISO data-to-MC scale factors
        ele,                // scale electron ID/ISO data-to-MC scale factors
        eleID,              // scale electron ID/ISO data-to-MC scale factors
        eleIDStat,          // scale electron ID/ISO data-to-MC scale factors
        eleIDSyst,          // scale electron ID/ISO data-to-MC scale factors
        eleReco,            // scale electron ID/ISO data-to-MC scale factors
        eleRecoStat,        // scale electron ID/ISO data-to-MC scale factors
        eleRecoSyst,        // scale electron ID/ISO data-to-MC scale factors
        muon,               // scale muon ID/ISO data-to-MC scale factors
        muonID,             // scale muon ID/ISO data-to-MC scale factors
        muonIDStat,         // scale muon ID/ISO data-to-MC scale factors
        muonIDSyst,         // scale muon ID/ISO data-to-MC scale factors
        muonIso,            // scale muon ID/ISO data-to-MC scale factors
        muonIsoStat,        // scale muon ID/ISO data-to-MC scale factors
        muonIsoSyst,        // scale muon ID/ISO data-to-MC scale factors
        eleScale,           // electron scale correction
        eleSmearingPhi,     // electron smearing correction
        eleSmearingRho,     // electron smearing correction
        eleScaleSmearing,   // electron scale+smearing enevlope
        muonScaleStat,      // muon scale correction
        muonScaleZpt,       // muon scale correction
        muonScaleEwk,       // muon scale correction
        muonScaleDeltaM,    // muon scale correction
        muonScaleEwk2,      // muon scale correction
        muonScale,          // muon scale correction envelope
        trig,               // scale trigger data-to-MC scale factors
        trigEta,            // scale trigger data-to-MC scale factors wrt eta (barrel-or-endcap) in antagonistic way
        pu,                 // scale pileup data-to-MC scale factors
        dy,                 // uncertainty on the Drell-Yan same-flavour background
        bg,                 // general background uncertainty
        dynorm,             // uncertainty on the Drell-Yan background estimation normalization
        kin,                // scale kinematic reconstruction scale factors
        jetPileupID,        // jet pileup-ID Data/MC scale factors
        
        btagBC,
        btagL,
        btagBCcorr,
        btagBCuncorr,
        btagLcorr,
        btagLuncorr,
        
        jerEta0,    // scale jet energy resolution scale factors (pt/eta bin0)
        jerEta1,    // scale jet energy resolution scale factors (pt/eta bin1)
        jerEta2Pt0, // scale jet energy resolution scale factors (pt/eta bin2)
        jerEta2Pt1, // scale jet energy resolution scale factors (pt/eta bin3)
        jerEta3Pt0, // scale jet energy resolution scale factors (pt/eta bin4)
        jerEta3Pt1, // scale jet energy resolution scale factors (pt/eta bin5)
        jer,        // scale jet energy resolution scale factors
        jesTotal,        // scale jet energy scale scale factors
        jesAbsoluteStat, //0
        jesAbsoluteScale,//1
        jesAbsoluteFlavMap,//2
        jesAbsoluteMPFBias,//3
        jesFragmentation,   // "HighPtExtra,    //4 -----
        jesSinglePionECAL, //5
        jesSinglePionHCAL, //6
        jesFlavorQCD,        //7
        jesTimePtEta,               // "Time" 8 -----
        jesRelativeJEREC1,    //9
        jesRelativeJEREC2,    //10
        jesRelativeJERHF,   //11
        jesRelativePtBB,       //12
        jesRelativePtEC1,    //13
        jesRelativePtEC2,    //14
        jesRelativePtHF,    //15
        jesRelativeBal, //new
        jesRelativeFSR,       //16
        jesRelativeSample,
        jesRelativeStatFSR,   //new
        jesRelativeStatEC,    //2,    //17 -----
        jesRelativeStatHF,   //18
        jesPileUpDataMC,       //19
        jesPileUpPtRef,
        jesPileUpPtBB,        //20
        jesPileUpPtEC1,        //21-----PileUpPtRef]
        jesPileUpPtEC2,
        jesPileUpPtHF,       //22
            //"PileUpBias,       //23-----
        jesPileUpMuZero,
        jesPileUpEnvelope,
        jesSubTotalPileUp,   //24
        jesSubTotalRelative,   //25
        jesSubTotalPt,       //26
        jesSubTotalScale,       //26
        jesSubTotalMC,       //27
        jesSubTotalAbsolute,       //26
        jesTotalNoFlavor,   //29
            //"TotalNoTime"       // ignoring for the moment these subtotals
            //"TotalNoFlavorNoTime"
            //"Time A-D
        jesFlavorZJet,       //30
        jesFlavorPhotonJet,   //31
        jesFlavorPureGluon,   //32
        jesFlavorPureQuark,   //33
        jesFlavorPureCharm,   //34
        jesFlavorPureBottom,//35
        jesCorrelationGroupMPFInSitu,//36
        jesCorrelationGroupIntercalibration,//37
        jesCorrelationGroupbJES,//38
        jesCorrelationGroupFlavor,//39
        jesCorrelationGroupUncorrelated,//40
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
        jesFlavorRealistic,

        frac_tthf,          // correction factor for the fraction of tt+HF events from the template fit
        frac_ttother,       // correction factor for the fraction of tt+Other events from the template fit
        lumi,               // luminosity uncertainty
        xsec_ttother,       // cross-section uncertainty of tt bkg. process
        xsec_dy,            // cross-section uncertainty of Drell-Yan process
        xsec_st,            // cross-section uncertainty of SingleTop process
        xsec_other,         // cross-section uncertainty of other bkg. process
        topPtTheory,        // scale top pt as predicted in theoretical ttbar differential cross-section calculations
        topPt,              // scale top pt as estimated in ttbar differential cross-section measurements, uncertainty via switching off and on
        mass,               // variations of masses used in process generation (here top quark mass)
        match,              // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter
        match_ttbb,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_tt2b,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttb,          // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttcc,         // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        match_ttother,      // matching uncertainty in process generation, associated to powheg-pythia hdamp parameter per ttbar+XX
        erdon,
        CR1,
        CR2,
        CR_envelope,         // envelope for color reconnection uncertainty if enevelope is taken for sum of ttbar samples
        CR_envelope_ind,    // envelope for color reconnection uncertainty if enevelope is taken for each ttbar sample individually
        tw_ds,
        meScale,            // Q2 scale uncertainty in process generation on Matrix Element only
        meScale_ttbb,       // Q2 scale uncertainty in process generation on Matrix Element only (ttbb process)
        meScale_ttb,        // Q2 scale uncertainty in process generation on Matrix Element only (ttb process)
        meScale_tt2b,       // Q2 scale uncertainty in process generation on Matrix Element only(tt2b process)
        meScale_ttcc,       // Q2 scale uncertainty in process generation on Matrix Element only (ttcc process)
        meScale_ttother,    // Q2 scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meScale_z,          // Q2 scale uncertainty in process generation on Matrix Element only (z+jets  jets process)
        meScale_st,         // Q2 scale uncertainty in process generation on Matrix Element only (single top jets process)
        meFacScale,            // Q2 factorization scale uncertainty in process generation on Matrix Element only
        meFacScale_ttbb,       // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbb process)
        meFacScale_ttb,        // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttb process)
        meFacScale_tt2b,       // Q2 factorization scale uncertainty in process generation on Matrix Element only(tt2b process)
        meFacScale_ttcc,       // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttcc process)
        meFacScale_ttother,    // Q2 factorization scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meFacScale_tt,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar process)
        meFacScale_z,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (z+jets process)
        meFacScale_w,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (w+jets process)
        meFacScale_st,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (single top process)
        meFacScale_ttv,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+X process)
        meFacScale_ttz,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+Z process)
        meFacScale_ttw,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+W process)
        meFacScale_ttg,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (ttbar+Gamma process)
        meFacScale_vv,     // Q2 factorization scale uncertainty in process generation on Matrix Element only (diboson process)
        meFacScale_ttdm,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (tt+DM process)
        meFacScale_htott_res,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (H/A->tt resonant process)
        meFacScale_htott_int,   // Q2 factorization scale uncertainty in process generation on Matrix Element only (H/A->tt interference process)
        meRenScale,            // Q2 renormalization scale uncertainty in process generation on Matrix Element only
        meRenScale_ttbb,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbb process)
        meRenScale_ttb,        // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttb process)
        meRenScale_tt2b,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only(tt2b process)
        meRenScale_ttcc,       // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttcc process)
        meRenScale_ttother,    // Q2 renormalization scale uncertainty in process generation on Matrix Element only (tt+light jets process)
        meRenScale_tt,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar process)
        meRenScale_z,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (z+jets process)
        meRenScale_w,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (w+jets process)
        meRenScale_st,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (single top process)
        meRenScale_ttv,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+X process)
        meRenScale_ttz,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+Z process)
        meRenScale_ttw,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+W process)
        meRenScale_ttg,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (ttbar+Gamma process)
        meRenScale_vv,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (diboson process)
        meRenScale_ttdm,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (tt+DM process)
        meRenScale_htott_res,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (H/A->tt resonant process)
        meRenScale_htott_int,     // Q2 renormalization scale uncertainty in process generation on Matrix Element only (H/A->tt interference process)

        psScale,            // Q2 scale uncertainty in process generation on Parton Shower only
        psScale_ttbb,       // Q2 scale uncertainty in process generation on Parton Shower only (ttbb process)
        psScale_ttb,        // Q2 scale uncertainty in process generation on Parton Shower only (ttb process)
        psScale_tt2b,       // Q2 scale uncertainty in process generation on Parton Shower only (tt2b process)
        psScale_ttcc,       // Q2 scale uncertainty in process generation on Parton Shower only (ttcc process)
        psScale_ttother,    // Q2 scale uncertainty in process generation on Parton Shower only (tt+light jets process)
        psISRScale,         // alpha_s^ISR scale uncertainty in process generation on Parton Shower only
        psISRScale_ttbb,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_tt2b,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttb,     // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttcc,    // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psISRScale_ttother, // alpha_s^ISR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale,         // alpha_s^FSR scale uncertainty in process generation on Parton Shower only
        psFSRScale_ttbb,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_tt2b,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttb,     // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttcc,    // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psFSRScale_ttother, // alpha_s^FSR scale uncertainty in process generation on Parton Shower only per ttbar+XX process
        psScaleWeight,      // alpha_s scale uncertainties in process generation on Parton Shower only
        psScaleWeight_st,   // alpha_s scale uncertainties in process generation on Parton Shower only per single top process
        scale,              // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower
        scale_ttbb,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttbb process)
        scale_ttb,          // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttb process)
        scale_tt2b,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (tt2b process)
        scale_ttcc,         // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (ttcc process)
        scale_ttother,      // Q2 scale uncertainty in process generation on Matrix Element and Parton Shower (tt+light jets process)
        bFrag,              // b quark fragmentation function (up and down)
        bFrag_central,      // b quark fragmentation function (so called central: it is not nominal)
        bFrag_Peterson,     // b quark fragmentation function (alternative Peterson fragmentation function)
        bSemilep,           // b hadron semileptonic decay branching ratios
        ueTune,             // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune
        ueTune_ttbb,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_tt2b,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttb,         // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttcc,        // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        ueTune_ttother,     // underlying event uncertainty stemming from orthogonal variations of the tuned parameters in the UE tune per tt+XX process
        powhegv2,           // POWHEGV2 event generator matched to PYTHIA8 shower
        powheg,             // POWHEG event generator matched to PYTHIA shower
        powhegv2Herwig,     // POWHEGV2 event generator matched to HERWIG++ shower
        powhegHerwig,       // POWHEG event generator matched to HERWIG shower
        powhegHelac,        // POWHEG ME generator with Helac as one loop provider (usually interfaced with Pythia8)
        powhegOpenloops,     // POWHEG ME generator with Openloops as one loop provider
        amcatnlofxfx,       // aMC@NLO event generator matched via FxFx to PYTHIA8 shower
        mcatnlo,            // MC@NLO event generator
        madgraphmlm,        // Madgraph event generator matched via MLM to PYTHIA8 shower
        cp5,                // cp5 tune for 2016 94X samples, used for Btag reweighting
        perugia11,          // Perugia11 parton shower tune
        perugia11NoCR,      // Perugia11 parton shower tune, no colour-reconnection
        alphasPdf,          // Variation of strong coupling in nominal PDF
        l1prefiring,          // Variation of L1 prefiring correction for 2016/2017
        normPdfGg,          // Normalization due to pdf in gg production
        normPdfGq,          // Normalization due to pdf in gqbar production
        normPdfQq,          // Normalization due to pdf in qqbar production
        normPdfTth,         // Normalization due to pdf gg prodution applied to ttH only (decoupled sys. from tt+XX)
        pdf_pca_1,          // PDF variation from PCA decomposition
        pdf_pca_2,          // PDF variation from PCA decomposition
        pdf,                // PDF variations
        pdf_envelope,       // Envelope of PDF variations
        unclustered,        // unclustered energy variation in the MET
        uncorrelatedType,   // Variations that aren't correlated between the FullRunII datasets
        closure,            // Closure test
        allAvailable,       // All systematics which are available
        all,                // All allowed systematics
        undefinedType,       // No systematic defined (also not nominal)
        
        jetPileupIDapplied,             // Check impact of jetPileupID
        jetLooseCleaningApplied,        // Check impact of loose cleaning
        met40Cut                        // Check impact of MET>40GeV in emu channel
        
    };



    /// Convert a type from string to enum
    Type convertType(const TString& type, bool const &quiet=false);

    /// Convert a type from enum to string
    TString convertType(const Type& type);
    
    /// Convert a type from string to string (removes variation from string)
    TString convertTypeString(const TString& type);

    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<TString>& types);

    /// Convert a vector of types from string to enum
    std::vector<Type> convertType(const std::vector<std::string>& types);

    /// Convert a vector of types from string to enum
    std::vector<TString> convertType(const std::vector<Type>& types);

    /// All variations as needed in any part of the framework
    enum Variation{up, down, central, undefinedVariation};

    /// Convert a variation from string to enum
    Variation convertVariation(const TString& variation);

    /// Convert a variation from enum to string
    TString convertVariation(const Variation& variation);

    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<TString>& variations);

    /// Convert a vector of variations from string to enum
    std::vector<Variation> convertVariation(const std::vector<std::string>& variations);

    /// Convert a vector of variations from string to enum
    std::vector<TString> convertVariation(const std::vector<Variation>& variations);







    /// Define for which systematics up/down variations are allowed
    const std::vector<Type> upDownTypes{
        lept, trig, trigEta, pu,
        ele, eleID, eleIDSyst, eleIDStat, eleReco, eleRecoSyst, eleRecoStat,
        muon, muonID, muonIDSyst, muonIDStat, muonIso, muonIsoSyst, muonIsoStat,
        eleScale, eleSmearingPhi, eleSmearingRho,
        muonScaleEwk, muonScaleStat, muonScaleZpt, muonScaleDeltaM, muonScaleEwk2,
        eleScaleSmearing, muonScale,
        dy, bg, kin,
        dynorm,
        jetPileupID,
        btagBC, btagL, btagBCcorr, btagBCuncorr, btagLcorr, btagLuncorr,
        jer, jerEta0, jerEta1, jerEta2Pt0, jerEta2Pt1, jerEta3Pt0, jerEta3Pt1,
        jesTotal, jesAbsoluteStat, jesAbsoluteScale, jesAbsoluteFlavMap, jesAbsoluteMPFBias, jesFragmentation, jesSinglePionECAL,
        jesSinglePionHCAL, jesFlavorQCD, jesTimePtEta, jesRelativeJEREC1, jesRelativeJEREC2, jesRelativeJERHF, jesRelativePtBB, jesRelativePtEC1,
        jesRelativePtEC2, jesRelativePtHF, jesRelativeBal, jesRelativeSample, jesRelativeFSR, jesRelativeStatFSR, jesRelativeStatEC, jesRelativeStatHF, jesPileUpDataMC,
        jesPileUpPtRef, jesPileUpPtEC1, jesPileUpPtEC2, jesPileUpPtHF, jesPileUpPtBB, jesPileUpMuZero, jesPileUpEnvelope, jesSubTotalPileUp,
        jesSubTotalRelative, jesSubTotalPt, jesSubTotalScale, jesSubTotalMC, jesSubTotalAbsolute, jesTotalNoFlavor,
        //"TotalNoTime"       // ignoring for the moment these subtotals
        //"TotalNoFlavorNoTime"
        //"Time A-D
        jesFlavorZJet, jesFlavorPhotonJet, jesFlavorPureGluon, jesFlavorPureQuark, jesFlavorPureCharm, jesFlavorPureBottom, jesFlavorRealistic,
        jesCorrelationGroupMPFInSitu, jesCorrelationGroupIntercalibration, jesCorrelationGroupbJES, jesCorrelationGroupFlavor, jesCorrelationGroupUncorrelated,
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
        frac_tthf, frac_ttother,
        lumi,
        xsec_ttother,xsec_dy,xsec_st,xsec_other,
        mass,
        match,
        match_ttbb, match_ttb, match_tt2b, match_ttcc, match_ttother,
        meScale, meScale_ttbb, meScale_ttb, meScale_tt2b, meScale_ttcc, meScale_ttother, meScale_z, meScale_st,
        meFacScale, meFacScale_ttbb, meFacScale_ttb, meFacScale_tt2b, meFacScale_ttcc, meFacScale_ttother,
        meRenScale, meRenScale_ttbb, meRenScale_ttb, meRenScale_tt2b, meRenScale_ttcc, meRenScale_ttother,
        meFacScale_tt, meFacScale_z, meFacScale_w, meFacScale_st, meFacScale_vv, meFacScale_ttv,
        meFacScale_ttz, meFacScale_ttw, meFacScale_ttg, meFacScale_ttdm, meFacScale_htott_res, meFacScale_htott_int,
        meRenScale_tt, meRenScale_z, meRenScale_w, meRenScale_st, meRenScale_vv, meRenScale_ttv,
        meRenScale_ttz, meRenScale_ttw, meRenScale_ttg, meRenScale_ttdm, meRenScale_htott_res, meRenScale_htott_int,
        psScale,
        psScale_ttbb, psScale_ttb, psScale_tt2b, psScale_ttcc, psScale_ttother,
        psISRScale,
        psISRScale_ttbb, psISRScale_ttb, psISRScale_tt2b, psISRScale_ttcc, psISRScale_ttother,
        psFSRScale,
        psFSRScale_ttbb, psFSRScale_ttb, psFSRScale_tt2b, psFSRScale_ttcc, psFSRScale_ttother,
        scale,
        //scale_ttbb, scale_ttb, scale_tt2b, scale_ttcc, scale_ttother,
        bFrag, bSemilep,bFrag_Peterson,
        erdon, CR1, CR2, tw_ds, CR_envelope, CR_envelope_ind,
        mtop,mtop_ind,
        ueTune,
        ueTune_ttbb, ueTune_ttb, ueTune_tt2b, ueTune_ttbb, ueTune_ttcc, ueTune_ttother,
        alphasPdf,l1prefiring,
        normPdfGg, normPdfGq, normPdfQq, normPdfTth,
        pdf_pca_1, pdf_pca_2, pdf_envelope, pdf, psScaleWeight,
        unclustered
    };

    /// Define for which systematics central variations are allowed
    /// This is also used to identify for which systematics variation numbers can be assigned
    const std::vector<Type> centralTypes{
        pdf, psScaleWeight, psScaleWeight_st
    };



    /// Check the validity of a variation for a given type
    void isValid(const Type& type, const Variation& variation, const int variationNumber =-1);

    /// Define jes systematics
    const std::vector<Type> jesTypes{
        jesTotal, jesAbsoluteStat, jesAbsoluteScale, jesAbsoluteFlavMap, jesAbsoluteMPFBias, jesFragmentation, jesSinglePionECAL,
        jesSinglePionHCAL, jesFlavorQCD, jesTimePtEta, jesRelativeJEREC1, jesRelativeJEREC2, jesRelativeJERHF, jesRelativePtBB, jesRelativePtEC1,
        jesRelativePtEC2, jesRelativePtHF, jesRelativeBal, jesRelativeSample, jesRelativeFSR, jesRelativeStatFSR, jesRelativeStatEC, jesRelativeStatHF, jesPileUpDataMC,
        jesPileUpPtRef, jesPileUpPtEC1, jesPileUpPtEC2, jesPileUpPtHF, jesPileUpPtBB, jesPileUpMuZero, jesPileUpEnvelope, jesSubTotalPileUp,
        jesSubTotalRelative, jesSubTotalPt, jesSubTotalScale, jesSubTotalMC, jesSubTotalAbsolute, jesTotalNoFlavor,
        //"TotalNoTime"       // ignoring for the moment these subtotals
        //"TotalNoFlavorNoTime"
        //"Time A-D
        jesFlavorZJet, jesFlavorPhotonJet, jesFlavorPureGluon, jesFlavorPureQuark, jesFlavorPureCharm, jesFlavorPureBottom, jesFlavorRealistic,
        jesCorrelationGroupMPFInSitu, jesCorrelationGroupIntercalibration, jesCorrelationGroupbJES, jesCorrelationGroupFlavor, jesCorrelationGroupUncorrelated,
        jesHttEta0to5, jesHttEta0to3, jesHttEta3to5, jesHttEC2,
        jesAbsolute, jesAbsoluteYear, jesBBEC1, jesBBEC1Year, jesEC2, jesEC2Year, jesHF, jesHFYear,
        jesUserDefinedHEM1516,
    };
    
    /// Define jes systematics for pure flavor
    const std::vector<Type> jesTypes_pureFlavor{
        jesFlavorPureGluon, jesFlavorPureQuark, jesFlavorPureCharm, jesFlavorPureBottom
    };
    
    /// Define jer systematics
    const std::vector<Type> jerTypes{
        jer, jerEta0, jerEta1, jerEta2Pt0, jerEta2Pt1, jerEta3Pt0, jerEta3Pt1
    };
    
    /// Define lepton scale and resolution systematics
    const std::vector<Type> leptonScaleResTypes{
        eleScale,eleSmearingPhi,eleSmearingRho,eleScaleSmearing,
        muonScaleStat,muonScaleZpt,muonScaleEwk,muonScaleDeltaM,muonScaleEwk2,muonScale
    };


    /// Define b-tag systematics, valid for all b-tag corrections
    const std::vector<Type> btagTypes{
        btagBC,btagL,
        btagBCcorr,btagBCuncorr,
        btagLcorr,btagLuncorr,
    };
    
    /// Define lepton SF systematics
    const std::vector<Type> leptonsfTypes{
        eleID,eleReco,
        muonID,muonIDStat,muonIDSyst,muonIso,muonIsoStat,muonIsoSyst,
    };

    /// Define ttbar systematics, i.e. variations of the ttbar sample (e.g. mass or scale variations)
    const std::vector<Type> ttbarTypes{
        topPtTheory, topPt,
        mass,
        match,
        //match_ttbb, match_ttb, match_tt2b, match_ttcc, match_ttother,
        erdon, CR1, CR2, CR_envelope, CR_envelope_ind,
        meScale, meScale_ttbb, meScale_ttb, meScale_tt2b, meScale_ttcc, meScale_ttother,
        meFacScale, meFacScale_ttbb, meFacScale_ttb, meFacScale_tt2b, meFacScale_ttcc, meFacScale_ttother,
        meRenScale, meRenScale_ttbb, meRenScale_ttb, meRenScale_tt2b, meRenScale_ttcc, meRenScale_ttother,
        psScale, psScale_ttbb, psScale_ttb, psScale_tt2b, psScale_ttcc, psScale_ttother,
        psISRScale,
        //psISRScale_ttbb, psISRScale_ttb, psISRScale_tt2b, psISRScale_ttcc, psISRScale_ttother,
        psFSRScale,
        //psFSRScale_ttbb, psFSRScale_ttb, psFSRScale_tt2b, psFSRScale_ttcc, psFSRScale_ttother,
        scale, scale_ttbb, scale_ttb, scale_tt2b, scale_ttcc, scale_ttother,
        bFrag, bFrag_central, bFrag_Peterson, bSemilep,
        ueTune,
        ueTune_ttbb, ueTune_ttb, ueTune_tt2b, ueTune_ttcc, ueTune_ttother,
        powhegv2, powheg, powhegv2Herwig, powhegHerwig, powhegHelac, powhegOpenloops, amcatnlofxfx, mcatnlo, madgraphmlm, cp5, perugia11, perugia11NoCR,
        alphasPdf, pdf, pdf_envelope, psScaleWeight,
        closure,
    };

    /// Define cross-section uncertainty systematics, which use nominal samples, and change only the scaling
    const std::vector<Type> crossSectionTypes{
        xsec_ttother,xsec_dy,xsec_st,xsec_other,
    };

    /// Define uncertainties due to tt+HF fraction scale factor from the fit, which use nominal samples, and change only the scaling
    const std::vector<Type> tthfFractionTypes{
        frac_tthf, frac_ttother
    };

    /// Define systematics that do not require dedicated root files
    const std::vector<Type> fileIndependentTypes{
        xsec_ttother,xsec_dy,xsec_st,xsec_other,
        dynorm,
        frac_tthf, frac_ttother,
        lumi,
        normPdfGg, normPdfGq, normPdfQq, normPdfTth
    };
    
    /// Define systematics that do not require dedicated b tagging efficiencies
    const std::vector<Type> noBtagEffTypes{
        nominal,lumi,
        btagBC,btagL,
        btagBCcorr,btagBCuncorr,
        btagLcorr,btagLuncorr,
        eleID,eleReco,
        muonID,muonIDStat,muonIDSyst,muonIso,muonIsoStat,muonIsoSyst,
        jetPileupIDapplied,jetLooseCleaningApplied,met40Cut,
        trig,
        pu,
                
        eleScale,eleSmearingPhi,eleSmearingRho,eleScaleSmearing,
        muonScaleStat,muonScaleZpt,muonScaleEwk,muonScaleDeltaM,muonScaleEwk2,muonScale,
        unclustered,
        
        pdf,
        xsec_ttother,xsec_dy,xsec_st,xsec_other
    };
    
    ///Define systematics that are applied by varying the nominal MC weight
    const std::vector<Type> mcWeightTypes{
        meScale,meFacScale,meRenScale,
        psISRScale,psFSRScale,
        bFrag,bSemilep,
        alphasPdf,pdf
    };
    
    ///Define systematics that require rescaling of lumi weight
    const std::vector<Type> lumiRescalingTypes{
        meScale,meFacScale,meRenScale,
        psISRScale,psFSRScale,
        bFrag,bSemilep,
        alphasPdf,
        pu

    };
    
    ///Define systematics that are applied by using an alternative sample
    const std::vector<Type> altSampleTypes{
        erdon,CR1,CR2,CR_envelope,CR_envelope_ind,ueTune,match,mTop169p5,mTop175p5,mtop,mtop_ind

    };
    
    ///Define systematic that only change the event weight, not the kinematics (same events as in nominal are selected!)
    //!!!!!!!!!!!!!!!If order is changed, order of weights in minTree production (distributions.cpp) has to be changes as well!!!!!!!!!!!!!!!!!!!!!!!!!!!
    const std::vector<Type> weightTypes{
        eleID,eleReco,
        muonID,muonIDStat,muonIDSyst,muonIso,muonIsoStat,muonIsoSyst,
        
        btagBC,btagL,
        // ~btagBCcorr,btagBCuncorr,
        // ~btagLcorr,btagLuncorr,
        
        trig,
        
        meScale,meFacScale,meRenScale,
        psISRScale,psFSRScale,
        bFrag,bSemilep,
        alphasPdf,pdf,
        pu,
        topPt,
        lumi
    };

    const std::vector<Type> uncorrelatedTypes{
        jer, jerEta0, jerEta1, jerEta2Pt0, jerEta2Pt1, jerEta3Pt0, jerEta3Pt1,
        jesAbsoluteStat, jesRelativeStatEC, jesRelativeStatFSR,
        jesRelativeJEREC1, jesRelativeJEREC2,
        jesRelativePtEC1, jesRelativePtEC2, jesTimePtEta,
        jesAbsoluteYear, jesBBEC1Year, jesEC2Year, jesHFYear, jesRelativeSample

    };

    /// Class for proper handling of systematic
    class Systematic{

     public:

        Systematic();

        Systematic(const Type& type, const Variation& variation, const int variationNumber =-1);

        Systematic(const TString& systematicName);

        ~Systematic(){}

        bool operator<(const Systematic& rhs)const{return this->name() < rhs.name();}

        TString name()const;

        Type type()const{return type_;}

        Variation variation()const{return variation_;}

        int variationNumber()const{return variationNumber_;}

        std::string type_str() const { return std::string(convertType(type_).Data()); }

    private:

        Type type_;

        Variation variation_;

        int variationNumber_;
    };

    /// Set all systematics from a list of allowed types, using the defined variations
    std::vector<Systematic> allowedSystematicsAnalysis(const std::vector<Type>& allowedTypes);

    /// Set up systematics from vector of systematicNames
    std::vector<Systematic> setSystematics(const std::vector<std::string>& systematicNames);

    /// Set up systematic for nominal (i.e. no systematic variation)
    Systematic nominalSystematic();

    /// Set up undefined systematic
    Systematic undefinedSystematic();
    
    /// Get correct pileup weight name
    TString puWeightName(const Systematic & systematic);
    
    /// Get correct met name (to derive unc. due to unclustered Energy)
    TString metNameAddition(const Systematic & systematic);
    
    /// Check if alternative sample and selected systematic match
    void checkAlternativeSample(const Systematic & systematic, const TString &dsSyst, const TString &dsName);
    
    ///Check if top pT reweighting should be applied based on current systematic;
    bool checkTopPTreweighting(const Systematic & systematic);
    
    ///Check the number of weights types systematics taking into accout that up down types require two weights
    int numberOfWeightTypes();
    
}

#endif
