//Mainly taken from https://gitlab.cern.ch/cms-desy-top/TopAnalysis/-/blob/master/Configuration/analysis/common/src/sampleHelpers.cc
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <algorithm>
#include <sstream>

#include <TString.h>
#include <TSystem.h>
#include <TObjArray.h>
#include <Rtypes.h>

#include "systematics.hpp"

// --------------------- Functions defined in namespace Systematic for Type -------------------------

Systematic::Type Systematic::convertType(const TString& type, bool const &quiet)
{
    // Attention: the order here is important, since the first line where the BeginsWith is true is returned
    if(type.BeginsWith("Nominal")) return nominal;
    if(type.BeginsWith("mH110")) return mH110;
    if(type.BeginsWith("mH115")) return mH115;
    if(type.BeginsWith("mH120")) return mH120;
    if(type.BeginsWith("mH1225")) return mH1225;
    if(type.BeginsWith("mH1275")) return mH1275;
    if(type.BeginsWith("mH130")) return mH130;
    if(type.BeginsWith("mH135")) return mH135;
    if(type.BeginsWith("mH140")) return mH140;
    if(type.BeginsWith("MTOP169p5")) return mTop169p5;
    if(type.BeginsWith("MTOP175p5")) return mTop175p5;
    if(type.BeginsWith("MTOP_IND")) return mtop_ind;
    if(type.BeginsWith("MTOP")) return mtop;
    if(type.BeginsWith("LEPT")) return lept;
    if(type.BeginsWith("ELECTRON_SCALESMEARING")) return eleScaleSmearing;
    if(type.BeginsWith("ELECTRON_SCALE")) return eleScale;
    if(type.BeginsWith("ELECTRON_SMEARING_PHI")) return eleSmearingPhi;
    if(type.BeginsWith("ELECTRON_SMEARING_RHO")) return eleSmearingRho;
    if(type.BeginsWith("ELECTRON_ID_STAT")) return eleIDStat;
    if(type.BeginsWith("ELECTRON_ID_SYST")) return eleIDSyst;
    if(type.BeginsWith("ELECTRON_ID")) return eleID;
    if(type.BeginsWith("ELECTRON_RECO_STAT")) return eleRecoStat;
    if(type.BeginsWith("ELECTRON_RECO_SYST")) return eleRecoSyst;
    if(type.BeginsWith("ELECTRON_RECO")) return eleReco;
    if(type.BeginsWith("ELECTRON")) return ele;
    if(type.BeginsWith("MUON_SCALE_EWK2")) return muonScaleEwk2;
    if(type.BeginsWith("MUON_SCALE_STAT")) return muonScaleStat;
    if(type.BeginsWith("MUON_SCALE_ZPT")) return muonScaleZpt;
    if(type.BeginsWith("MUON_SCALE_DELTAM")) return muonScaleDeltaM;
    if(type.BeginsWith("MUON_SCALE_EWK")) return muonScaleEwk;
    if(type.BeginsWith("MUON_SCALE")) return muonScale;
    if(type.BeginsWith("MUON_ISO_STAT")) return muonIsoStat;
    if(type.BeginsWith("MUON_ISO_SYST")) return muonIsoSyst;
    if(type.BeginsWith("MUON_ISO")) return muonIso;
    if(type.BeginsWith("MUON_ID_STAT")) return muonIDStat;
    if(type.BeginsWith("MUON_ID_SYST")) return muonIDSyst;
    if(type.BeginsWith("MUON_ID")) return muonID;
    if(type.BeginsWith("MUON")) return muon;
    if(type.BeginsWith("TRIG_ETA")) return trigEta;
    if(type.BeginsWith("TRIG")) return trig;
    if(type.BeginsWith("PU")) return pu;
    if(type.BeginsWith("DYNORM")) return dynorm;
    if(type.BeginsWith("DY")) return dy;
    if(type.BeginsWith("BG")) return bg;
    if(type.BeginsWith("KIN")) return kin;
    if(type.BeginsWith("JetPileupID")) return jetPileupID;
    if(type.BeginsWith("BTAGBC_CORR")) return btagBCcorr;
    if(type.BeginsWith("BTAGBC_UNCORR")) return btagBCuncorr;
    if(type.BeginsWith("BTAGL_CORR")) return btagLcorr;
    if(type.BeginsWith("BTAGL_UNCORR")) return btagLuncorr;
    if(type.BeginsWith("BTAGBC")) return btagBC;
    if(type.BeginsWith("BTAGL")) return btagL;
    if(type.BeginsWith("JEREta0")) return jerEta0;
    if(type.BeginsWith("JEREta1")) return jerEta1;
    if(type.BeginsWith("JEREta2Pt0")) return jerEta2Pt0;
    if(type.BeginsWith("JEREta2Pt1")) return jerEta2Pt1;
    if(type.BeginsWith("JEREta3Pt0")) return jerEta3Pt0;
    if(type.BeginsWith("JEREta3Pt1")) return jerEta3Pt1;
    if(type.BeginsWith("JER")) return jer;
    if(type.BeginsWith("JESRelativeBalreg")) return jesRelativeBal_reg;
    if(type.BeginsWith("JESFlavorQCDreg")) return jesFlavorQCD_reg;
    if(type.BeginsWith("JESFlavorRealisticreg")) return jesFlavorRealistic_reg;
    if(type.BeginsWith("JESRelativeSampleYear")) return jesRelativeSampleYear;
    if(type.BeginsWith("JESAbsoluteStat")) return jesAbsoluteStat;
    if(type.BeginsWith("JESAbsoluteScale")) return jesAbsoluteScale;
    if(type.BeginsWith("JESAbsoluteFlavMap")) return jesAbsoluteFlavMap;
    if(type.BeginsWith("JESAbsoluteMPFBias")) return jesAbsoluteMPFBias;
    if(type.BeginsWith("JESFragmentation")) return jesFragmentation;
    if(type.BeginsWith("JESSinglePionECAL")) return jesSinglePionECAL;
    if(type.BeginsWith("JESSinglePionHCAL")) return jesSinglePionHCAL;
    if(type.BeginsWith("JESFlavorQCD")) return jesFlavorQCD;
    if(type.BeginsWith("JESTimePtEta")) return jesTimePtEta;
//    if(type.BeginsWith("JESTimePt")) return jesTimePt;
    if(type.BeginsWith("JESRelativeJEREC1")) return jesRelativeJEREC1;
    if(type.BeginsWith("JESRelativeJEREC2")) return jesRelativeJEREC2;
    if(type.BeginsWith("JESRelativeJERHF")) return jesRelativeJERHF;
    if(type.BeginsWith("JESRelativePtBB")) return jesRelativePtBB;
    if(type.BeginsWith("JESRelativePtEC1")) return jesRelativePtEC1;
    if(type.BeginsWith("JESRelativePtEC2")) return jesRelativePtEC2;
    if(type.BeginsWith("JESRelativePtHF")) return jesRelativePtHF;
    if(type.BeginsWith("JESRelativeBal")) return jesRelativeBal;
    if(type.BeginsWith("JESRelativeSample")) return jesRelativeSample;
    if(type.BeginsWith("JESRelativeFSR")) return jesRelativeFSR;
    if(type.BeginsWith("JESRelativeStatFSR")) return jesRelativeStatFSR;
    if(type.BeginsWith("JESRelativeStatEC")) return jesRelativeStatEC;
    if(type.BeginsWith("JESRelativeStatHF")) return jesRelativeStatHF;
    if(type.BeginsWith("JESPileUpDataMC")) return jesPileUpDataMC;
    if(type.BeginsWith("JESPileUpPtRef")) return jesPileUpPtRef;
    if(type.BeginsWith("JESPileUpPtBB")) return jesPileUpPtBB;
    if(type.BeginsWith("JESPileUpPtEC1")) return jesPileUpPtEC1;
    if(type.BeginsWith("JESPileUpPtEC2")) return jesPileUpPtEC2;
    if(type.BeginsWith("JESPileUpPtHF")) return jesPileUpPtHF;
            //"PileUpBias,       //23-----
    if(type.BeginsWith("JESPileUpMuZero")) return jesPileUpMuZero;
    if(type.BeginsWith("JESPileUpEnvelope")) return jesPileUpEnvelope;
    if(type.BeginsWith("JESSubTotalPileUp")) return jesSubTotalPileUp;
    if(type.BeginsWith("JESSubTotalRelative")) return jesSubTotalRelative;
    if(type.BeginsWith("JESSubTotalPt")) return jesSubTotalPt;
    if(type.BeginsWith("JESSubTotalScale")) return jesSubTotalScale;
    if(type.BeginsWith("JESSubTotalMC")) return jesSubTotalMC;
    if(type.BeginsWith("JESSubTotalAbsolute")) return jesSubTotalAbsolute;
    if(type.BeginsWith("JESTotalNoFlavor")) return jesTotalNoFlavor;
            //"TotalNoTime"       // ignoring for the moment these subtotals
            //"TotalNoFlavorNoTime"
            //"Time A-D
    if(type.BeginsWith("JESFlavorZJet")) return jesFlavorZJet;
    if(type.BeginsWith("JESFlavorPhotonJet")) return jesFlavorPhotonJet;
    if(type.BeginsWith("JESFlavorPureGluon")) return jesFlavorPureGluon;
    if(type.BeginsWith("JESFlavorPureQuark")) return jesFlavorPureQuark;
    if(type.BeginsWith("JESFlavorPureCharm")) return jesFlavorPureCharm;
    if(type.BeginsWith("JESFlavorPureBottom")) return jesFlavorPureBottom;
    if(type.BeginsWith("JESFlavorRealistic")) return jesFlavorRealistic;
    if(type.BeginsWith("JESCorrelationGroupMPFInSitu")) return jesCorrelationGroupMPFInSitu;
    if(type.BeginsWith("JESCorrelationGroupIntercalibration")) return jesCorrelationGroupIntercalibration;
    if(type.BeginsWith("JESCorrelationGroupbJES")) return jesCorrelationGroupbJES;
    if(type.BeginsWith("JESCorrelationGroupFlavor")) return jesCorrelationGroupFlavor;
    if(type.BeginsWith("JESCorrelationGroupUncorrelated")) return jesCorrelationGroupUncorrelated;
    if(type.BeginsWith("JESHttEta0to5")) return jesHttEta0to5;
    if(type.BeginsWith("JESHttEta0to3")) return jesHttEta0to3;
    if(type.BeginsWith("JESHttEta3to5")) return jesHttEta3to5;
    if(type.BeginsWith("JESHttEC2")) return jesHttEC2;
    if(type.BeginsWith("JESAbsoluteYear")) return jesAbsoluteYear;
    if(type.BeginsWith("JESAbsolute")) return jesAbsolute;
    if(type.BeginsWith("JESBBEC1Year")) return jesBBEC1Year;
    if(type.BeginsWith("JESBBEC1")) return jesBBEC1;
    if(type.BeginsWith("JESEC2Year")) return jesEC2Year;
    if(type.BeginsWith("JESEC2")) return jesEC2;
    if(type.BeginsWith("JESHFYear")) return jesHFYear;
    if(type.BeginsWith("JESHF")) return jesHF; 
    if(type.BeginsWith("JESUserDefinedHEM1516")) return jesUserDefinedHEM1516;
    if(type.BeginsWith("JESTotal")) return jesTotal;

    if(type.BeginsWith("FRAC_TTHF")) return frac_tthf;
    if(type.BeginsWith("FRAC_TTOTHER")) return frac_ttother;
    if(type.BeginsWith("LUMI")) return lumi;
    if(type.BeginsWith("XSEC_TTOTHER")) return xsec_ttother;
    if(type.BeginsWith("XSEC_DY")) return xsec_dy;
    if(type.BeginsWith("XSEC_ST")) return xsec_st;
    if(type.BeginsWith("XSEC_OTHER")) return xsec_other;
    if(type.BeginsWith("TOP_PT_THEORY")) return topPtTheory;
    if(type.BeginsWith("TOP_PT")) return topPt;
    if(type.BeginsWith("MASS")) return mass;
    if(type.BeginsWith("MATCH_TTBB"))    return match_ttbb;
    if(type.BeginsWith("MATCH_TT2B"))    return match_tt2b;
    if(type.BeginsWith("MATCH_TTB"))     return match_ttb;
    if(type.BeginsWith("MATCH_TTCC"))    return match_ttcc;
    if(type.BeginsWith("MATCH_TTOTHER")) return match_ttother;
    if(type.BeginsWith("MATCH")) return match;
    if(type.BeginsWith("CR1")) return CR1;
    if(type.BeginsWith("CR2")) return CR2;
    if(type.BeginsWith("CR_ENVELOPE_IND")) return CR_envelope_ind;
    if(type.BeginsWith("CR_ENVELOPE")) return CR_envelope;
    if(type.BeginsWith("ERDON")) return erdon;
    if(type.BeginsWith("DS")) return tw_ds;
    if(type.BeginsWith("MESCALE_ENVELOPE_IND")) return meScale_envelope_ind;
    if(type.BeginsWith("MESCALE_ENVELOPE")) return meScale_envelope;
    if(type.BeginsWith("MESCALE_TTBB")) return meScale_ttbb;
    if(type.BeginsWith("MESCALE_TTB")) return meScale_ttb;
    if(type.BeginsWith("MESCALE_TT2B")) return meScale_tt2b;
    if(type.BeginsWith("MESCALE_TTCC")) return meScale_ttcc;
    if(type.BeginsWith("MESCALE_TTOTHER")) return meScale_ttother;
    if(type.BeginsWith("MESCALE_Z")) return meScale_z;
    if(type.BeginsWith("MESCALE_ST")) return meScale_st;
    if(type.BeginsWith("MESCALE")) return meScale;
    if(type.BeginsWith("MEFACSCALE_TTBB")) return meFacScale_ttbb;
    if(type.BeginsWith("MEFACSCALE_TTB")) return meFacScale_ttb;
    if(type.BeginsWith("MEFACSCALE_TT2B")) return meFacScale_tt2b;
    if(type.BeginsWith("MEFACSCALE_TTCC")) return meFacScale_ttcc;
    if(type.BeginsWith("MEFACSCALE_TTOTHER")) return meFacScale_ttother;
    if(type.BeginsWith("MEFACSCALE_TTV")) return meFacScale_ttv;
    if(type.BeginsWith("MEFACSCALE_TTZ")) return meFacScale_ttz;
    if(type.BeginsWith("MEFACSCALE_TTW")) return meFacScale_ttw;
    if(type.BeginsWith("MEFACSCALE_TTG")) return meFacScale_ttg;
    if(type.BeginsWith("MEFACSCALE_TTDM")) return meFacScale_ttdm;
    if(type.BeginsWith("MEFACSCALE_HTOTT_RES")) return meFacScale_htott_res;
    if(type.BeginsWith("MEFACSCALE_HTOTT_INT")) return meFacScale_htott_int;
    if(type.BeginsWith("MEFACSCALE_TT")) return meFacScale_tt;
    if(type.BeginsWith("MEFACSCALE_Z")) return meFacScale_z;
    if(type.BeginsWith("MEFACSCALE_W")) return meFacScale_w;
    if(type.BeginsWith("MEFACSCALE_ST")) return meFacScale_st;
    if(type.BeginsWith("MEFACSCALE_VV")) return meFacScale_vv;
    if(type.BeginsWith("MEFACSCALE")) return meFacScale;
    if(type.BeginsWith("MERENSCALE_TTBB")) return meRenScale_ttbb;
    if(type.BeginsWith("MERENSCALE_TTB")) return meRenScale_ttb;
    if(type.BeginsWith("MERENSCALE_TT2B")) return meRenScale_tt2b;
    if(type.BeginsWith("MERENSCALE_TTCC")) return meRenScale_ttcc;
    if(type.BeginsWith("MERENSCALE_TTOTHER")) return meRenScale_ttother;
    if(type.BeginsWith("MERENSCALE_TTV")) return meRenScale_ttv;
    if(type.BeginsWith("MERENSCALE_TTZ")) return meRenScale_ttz;
    if(type.BeginsWith("MERENSCALE_TTW")) return meRenScale_ttw;
    if(type.BeginsWith("MERENSCALE_TTG")) return meRenScale_ttg;
    if(type.BeginsWith("MERENSCALE_TTDM")) return meRenScale_ttdm;
    if(type.BeginsWith("MERENSCALE_HTOTT_RES")) return meRenScale_htott_res;
    if(type.BeginsWith("MERENSCALE_HTOTT_INT")) return meRenScale_htott_int;
    if(type.BeginsWith("MERENSCALE_TT")) return meRenScale_tt;
    if(type.BeginsWith("MERENSCALE_Z")) return meRenScale_z;
    if(type.BeginsWith("MERENSCALE_W")) return meRenScale_w;
    if(type.BeginsWith("MERENSCALE_ST")) return meRenScale_st;
    if(type.BeginsWith("MERENSCALE_VV")) return meRenScale_vv;
    if(type.BeginsWith("MERENSCALE")) return meRenScale;
    if(type.BeginsWith("PSSCALE_TTBB")) return psScale_ttbb;
    if(type.BeginsWith("PSSCALE_TTB")) return psScale_ttb;
    if(type.BeginsWith("PSSCALE_TT2B")) return psScale_tt2b;
    if(type.BeginsWith("PSSCALE_TTCC")) return psScale_ttcc;
    if(type.BeginsWith("PSSCALE_TTOTHER")) return psScale_ttother;
    if(type.BeginsWith("PSSCALE_WEIGHT_ST")) return psScaleWeight_st;
    if(type.BeginsWith("PSSCALE_WEIGHT")) return psScaleWeight;
    if(type.BeginsWith("PSSCALE")) return psScale;
    if(type.BeginsWith("PSISRSCALE_TTBB"))    return psISRScale_ttbb;
    if(type.BeginsWith("PSISRSCALE_TT2B"))    return psISRScale_tt2b;
    if(type.BeginsWith("PSISRSCALE_TTB"))     return psISRScale_ttb;
    if(type.BeginsWith("PSISRSCALE_TTCC"))    return psISRScale_ttcc;
    if(type.BeginsWith("PSISRSCALE_TTOTHER")) return psISRScale_ttother;
    if(type.BeginsWith("PSISRSCALE")) return psISRScale;
    if(type.BeginsWith("PSFSRSCALE_TTBB"))    return psFSRScale_ttbb;
    if(type.BeginsWith("PSFSRSCALE_TT2B"))    return psFSRScale_tt2b;
    if(type.BeginsWith("PSFSRSCALE_TTB"))     return psFSRScale_ttb;
    if(type.BeginsWith("PSFSRSCALE_TTCC"))    return psFSRScale_ttcc;
    if(type.BeginsWith("PSFSRSCALE_TTOTHER")) return psFSRScale_ttother;
    if(type.BeginsWith("PSFSRSCALE")) return psFSRScale;
    if(type.BeginsWith("SCALE_TTBB")) return scale_ttbb;
    if(type.BeginsWith("SCALE_TTB")) return scale_ttb;
    if(type.BeginsWith("SCALE_TT2B")) return scale_tt2b;
    if(type.BeginsWith("SCALE_TTCC")) return scale_ttcc;
    if(type.BeginsWith("SCALE_TTOTHER")) return scale_ttother;
    if(type.BeginsWith("SCALE")) return scale;
    if(type.BeginsWith("BFRAG_CENTRAL")) return bFrag_central;
    if(type.BeginsWith("BFRAG_PETERSON")) return bFrag_Peterson;
    if(type.BeginsWith("BFRAG")) return bFrag;
    if(type.BeginsWith("BSEMILEP")) return bSemilep;
    if(type.BeginsWith("UETUNE_TTBB"))    return ueTune_ttbb;
    if(type.BeginsWith("UETUNE_TT2B"))    return ueTune_tt2b;
    if(type.BeginsWith("UETUNE_TTB"))     return ueTune_ttb;
    if(type.BeginsWith("UETUNE_TTCC"))    return ueTune_ttcc;
    if(type.BeginsWith("UETUNE_TTOTHER")) return ueTune_ttother;
    if(type.BeginsWith("UETUNE")) return ueTune;
    if(type.BeginsWith("POWHEGV2HERWIG")) return powhegv2Herwig;
    if(type.BeginsWith("POWHEGHERWIG")) return powhegHerwig;
    if(type.BeginsWith("POWHEGHELAC")) return powhegHelac;
    if(type.BeginsWith("POWHEGOPENLOOPS")) return powhegOpenloops;
    if(type.BeginsWith("POWHEGV2")) return powhegv2;
    if(type.BeginsWith("POWHEG")) return powheg;
    if(type.BeginsWith("AMCATNLOFXFX")) return amcatnlofxfx;
    if(type.BeginsWith("MCATNLO")) return mcatnlo;
    if(type.BeginsWith("MADGRAPHMLM")) return madgraphmlm;
    if(type.BeginsWith("CP5")) return cp5;
    if(type.BeginsWith("PERUGIA11NoCR")) return perugia11NoCR;
    if(type.BeginsWith("PERUGIA11")) return perugia11;
    if(type.BeginsWith("PDF_ALPHAS")) return alphasPdf;
    if(type.BeginsWith("L1PREFIRING")) return l1prefiring;
    if(type.BeginsWith("NORMPDFGG")) return normPdfGg;
    if(type.BeginsWith("NORMPDFGQ")) return normPdfGq;
    if(type.BeginsWith("NORMPDFQQ")) return normPdfQq;
    if(type.BeginsWith("NORMPDFTTH")) return normPdfTth;
    if(type.BeginsWith("PDF_PCA_1")) return pdf_pca_1;
    if(type.BeginsWith("PDF_PCA_2")) return pdf_pca_2;
    if(type.BeginsWith("PDF_ENVELOPE")) return pdf_envelope;
    if(type.BeginsWith("PDF")) return pdf;
    if(type.BeginsWith("UNCLUSTERED")) return unclustered;
    if(type.BeginsWith("closure")) return closure;
    if(type.BeginsWith("allAvailable")) return allAvailable;
    if(type.BeginsWith("all")) return all;
    if(type.BeginsWith("jetPileupIDapplied")) return jetPileupIDapplied;
    if(type.BeginsWith("jetLooseCleaningApplied")) return jetLooseCleaningApplied;
    if(type.BeginsWith("met40Cut")) return met40Cut;
    if(!quiet) std::cout<<"Warning in Systematic::convertType()! Following conversion is not implemented: "
             <<type<<std::endl<<std::endl;
    return undefinedType;
}



TString Systematic::convertType(const Type& type)
{
    if(type == nominal) return "Nominal";
    if(type == mH110) return "mH110";
    if(type == mH115) return "mH115";
    if(type == mH120) return "mH120";
    if(type == mH1225) return "mH1225";
    if(type == mH1275) return "mH1275";
    if(type == mH130) return "mH130";
    if(type == mH135) return "mH135";
    if(type == mH140) return "mH140";
    if(type == mTop169p5) return "MTOP169p5";
    if(type == mTop175p5) return "MTOP175p5";
    if(type == mtop_ind) return "MTOP_IND";
    if(type == mtop) return "MTOP";
    if(type == lept) return "LEPT";
    if(type == ele) return "ELE";
    if(type == eleID) return "ELECTRON_ID";
    if(type == eleIDSyst) return "ELECTRON_ID_SYST";
    if(type == eleIDStat) return "ELECTRON_ID_STAT";
    if(type == eleReco) return "ELECTRON_RECO";
    if(type == eleRecoSyst) return "ELECTRON_RECO_SYST";
    if(type == eleRecoStat) return "ELECTRON_RECO_STAT";
    if(type == eleScale) return "ELECTRON_SCALE";
    if(type == eleScaleSmearing) return "ELECTRON_SCALESMEARING";
    if(type == eleSmearingRho) return "ELECTRON_SMEARING_RHO";
    if(type == eleSmearingPhi) return "ELECTRON_SMEARING_PHI";
    if(type == muon) return "MUON";
    if(type == muonID) return "MUON_ID";
    if(type == muonIDSyst) return "MUON_ID_SYST";
    if(type == muonIDStat) return "MUON_ID_STAT";
    if(type == muonIso) return "MUON_ISO";
    if(type == muonIsoSyst) return "MUON_ISO_SYST";
    if(type == muonIsoStat) return "MUON_ISO_STAT";
    if(type == muonScale) return "MUON_SCALE";
    if(type == muonScaleZpt) return "MUON_SCALE_ZPT";
    if(type == muonScaleEwk) return "MUON_SCALE_EWK";
    if(type == muonScaleEwk2) return "MUON_SCALE_EWK2";
    if(type == muonScaleStat) return "MUON_SCALE_STAT";
    if(type == muonScaleDeltaM) return "MUON_SCALE_DELTAM";
    if(type == trigEta) return "TRIG_ETA";
    if(type == trig) return "TRIG";
    if(type == pu) return "PU";
    if(type == dy) return "DY";
    if(type == bg) return "BG";
    if(type == dynorm) return "DYNORM";
    if(type == kin) return "KIN";
    if(type == jetPileupID) return "JetPileupID";
    if(type == btagBC) return "BTAGBC";
    if(type == btagL) return "BTAGL";
    if(type == btagBCcorr) return "BTAGBC_CORR";
    if(type == btagBCuncorr) return "BTAGBC_UNCORR";
    if(type == btagLcorr) return "BTAGL_CORR";
    if(type == btagLuncorr) return "BTAGL_UNCORR";
    if(type == jerEta0) return "JEREta0";
    if(type == jerEta1) return "JEREta1";
    if(type == jerEta2Pt0) return "JEREta2Pt0";
    if(type == jerEta2Pt1) return "JEREta2Pt1";
    if(type == jerEta3Pt0) return "JEREta3Pt0";
    if(type == jerEta3Pt1) return "JEREta3Pt1";
    if(type == jer) return "JER";
    if(type == jesTotal) return "JESTotal";
    if(type == jesAbsoluteStat) return "JESAbsoluteStat";
    if(type == jesAbsoluteScale) return "JESAbsoluteScale";
    if(type == jesAbsoluteFlavMap) return "JESAbsoluteFlavMap";
    if(type == jesAbsoluteMPFBias) return "JESAbsoluteMPFBias";
    if(type == jesFragmentation) return "JESFragmentation";
    if(type == jesSinglePionECAL) return "JESSinglePionECAL";
    if(type == jesSinglePionHCAL) return "JESSinglePionHCAL";
    if(type == jesFlavorQCD) return "JESFlavorQCD";
    if(type == jesTimePtEta) return "JESTimePtEta";
    if(type == jesRelativeJEREC1) return "JESRelativeJEREC1";
    if(type == jesRelativeJEREC2) return "JESRelativeJEREC2";
    if(type == jesRelativeJERHF) return "JESRelativeJERHF";
    if(type == jesRelativePtBB) return "JESRelativePtBB";
    if(type == jesRelativePtEC1) return "JESRelativePtEC1";
    if(type == jesRelativePtEC2) return "JESRelativePtEC2";
    if(type == jesRelativePtHF) return "JESRelativePtHF";
    if(type == jesRelativeBal) return "JESRelativeBal";
    if(type == jesRelativeSample) return "JESRelativeSample";
    if(type == jesRelativeFSR) return "JESRelativeFSR";
    if(type == jesRelativeStatFSR) return "JESRelativeStatFSR";
    if(type == jesRelativeStatEC) return "JESRelativeStatEC";
    if(type == jesRelativeStatHF) return "JESRelativeStatHF";
    if(type == jesPileUpDataMC) return "JESPileUpDataMC";
    if(type == jesPileUpPtRef) return "JESPileUpPtRef";
    if(type == jesPileUpPtBB) return "JESPileUpPtBB";
    if(type == jesPileUpPtEC1) return "JESPileUpPtEC1";
    if(type == jesPileUpPtEC2) return "JESPileUpPtEC2";
    if(type == jesPileUpPtHF) return "JESPileUpPtHF";
            //"PileUpBias,       //23-----
    if(type == jesPileUpMuZero) return "JESPileUpMuZero";
    if(type == jesPileUpEnvelope) return "JESPileUpEnvelope";
    if(type == jesSubTotalPileUp) return "JESSubTotalPileUp";
    if(type == jesSubTotalRelative) return "JESSubTotalRelative";
    if(type == jesSubTotalPt) return "JESSubTotalPt";
    if(type == jesSubTotalScale) return "JESSubTotalScale";
    if(type == jesSubTotalMC) return "JESSubTotalMC";
    if(type == jesSubTotalAbsolute) return "JESSubTotalAbsolute";
    if(type == jesTotalNoFlavor) return "JESTotalNoFlavor";
            //"TotalNoTime"       // ignoring for the moment these subtotals
            //"TotalNoFlavorNoTime"
            //"Time A-D
    if(type == jesFlavorZJet) return "JESFlavorZJet";
    if(type == jesFlavorPhotonJet) return "JESFlavorPhotonJet";
    if(type == jesFlavorPureGluon) return "JESFlavorPureGluon";
    if(type == jesFlavorPureQuark) return "JESFlavorPureQuark";
    if(type == jesFlavorPureCharm) return "JESFlavorPureCharm";
    if(type == jesFlavorPureBottom) return "JESFlavorPureBottom";
    if(type == jesFlavorRealistic) return "JESFlavorRealistic";
    if(type == jesCorrelationGroupMPFInSitu) return "JESCorrelationGroupMPFInSitu";
    if(type == jesCorrelationGroupIntercalibration) return "JESCorrelationGroupIntercalibration";
    if(type == jesCorrelationGroupbJES) return "JESCorrelationGroupbJES";
    if(type == jesCorrelationGroupFlavor) return "JESCorrelationGroupFlavor";
    if(type == jesCorrelationGroupUncorrelated) return "JESCorrelationGroupUncorrelated";
    if(type == jesHttEta0to5) return "JESHttEta0to5";
    if(type == jesHttEta0to3) return "JESHttEta0to3";
    if(type == jesHttEta3to5) return "JESHttEta3to5";
    if(type == jesHttEC2) return "JESHttEC2";
    if(type == jesAbsoluteYear) return "JESAbsoluteYear";
    if(type == jesAbsolute) return "JESAbsolute";
    if(type == jesBBEC1Year) return "JESBBEC1Year";
    if(type == jesBBEC1) return "JESBBEC1";
    if(type == jesEC2Year) return "JESEC2Year";
    if(type == jesEC2) return "JESEC2";
    if(type == jesHFYear) return "JESHFYear";
    if(type == jesHF) return "JESHF";
    if(type == jesRelativeBal_reg) return "JESRelativeBalreg";
    if(type == jesFlavorQCD_reg) return "JESFlavorQCDreg";
    if(type == jesFlavorRealistic_reg) return "JESFlavorRealisticreg";
    if(type == jesRelativeSampleYear) return "JESRelativeSampleYear";
    if(type == jesUserDefinedHEM1516) return "JESUserDefinedHEM1516";

    if(type == frac_tthf) return "FRAC_TTHF";
    if(type == frac_ttother) return "FRAC_TTOTHER";
    if(type == lumi) return "LUMI";
    if(type == xsec_ttother) return "XSEC_TTOTHER";
    if(type == xsec_dy) return "XSEC_DY";
    if(type == xsec_st) return "XSEC_ST";
    if(type == xsec_other) return "XSEC_OTHER";
    if(type == topPtTheory) return "TOP_PT_THEORY";
    if(type == topPt) return "TOP_PT";
    if(type == mass) return "MASS";
    if(type == match_ttbb)    return "MATCH_TTBB";
    if(type == match_tt2b)    return "MATCH_TT2B";
    if(type == match_ttb)     return "MATCH_TTB";
    if(type == match_ttcc)    return "MATCH_TTCC";
    if(type == match_ttother) return "MATCH_TTOTHER";
    if(type == match) return "MATCH";
    if(type == CR1) return "CR1";
    if(type == CR2) return "CR2";
    if(type == CR_envelope_ind) return "CR_ENVELOPE_IND";
    if(type == CR_envelope) return "CR_ENVELOPE";
    if(type == erdon) return "ERDON";
    if(type == tw_ds) return "DS";
    if(type == meScale_envelope_ind) return "MESCALE_ENVELOPE_IND";
    if(type == meScale_envelope) return "MESCALE_ENVELOPE";
    if(type == meScale_ttbb) return "MESCALE_TTBB";
    if(type == meScale_ttb) return "MESCALE_TTB";
    if(type == meScale_tt2b) return "MESCALE_TT2B";
    if(type == meScale_ttcc) return "MESCALE_TTCC";
    if(type == meScale_ttother) return "MESCALE_TTOTHER";
    if(type == meScale_z) return "MESCALE_Z";
    if(type == meScale_st) return "MESCALE_ST";
    if(type == meScale) return "MESCALE";
    if(type == meFacScale_ttbb) return "MEFACSCALE_TTBB";
    if(type == meFacScale_ttb) return "MEFACSCALE_TTB";
    if(type == meFacScale_tt2b) return "MEFACSCALE_TT2B";
    if(type == meFacScale_ttcc) return "MEFACSCALE_TTCC";
    if(type == meFacScale_ttother) return "MEFACSCALE_TTOTHER";
    if(type == meFacScale_tt) return "MEFACSCALE_TT";
    if(type == meFacScale_z) return "MEFACSCALE_Z";
    if(type == meFacScale_w) return "MEFACSCALE_W";
    if(type == meFacScale_st) return "MEFACSCALE_ST";
    if(type == meFacScale_vv) return "MEFACSCALE_VV";
    if(type == meFacScale_ttv) return "MEFACSCALE_TTV";
    if(type == meFacScale_ttz) return "MEFACSCALE_TTZ";
    if(type == meFacScale_ttw) return "MEFACSCALE_TTW";
    if(type == meFacScale_ttg) return "MEFACSCALE_TTG";
    if(type == meFacScale_ttdm) return "MEFACSCALE_TTDM";
    if(type == meFacScale_htott_res) return "MEFACSCALE_HTOTT_RES";
    if(type == meFacScale_htott_int) return "MEFACSCALE_HTOTT_INT";
    if(type == meFacScale) return "MEFACSCALE";
    if(type == meRenScale_ttbb) return "MERENSCALE_TTBB";
    if(type == meRenScale_ttb) return "MERENSCALE_TTB";
    if(type == meRenScale_tt2b) return "MERENSCALE_TT2B";
    if(type == meRenScale_ttcc) return "MERENSCALE_TTCC";
    if(type == meRenScale_ttother) return "MERENSCALE_TTOTHER";
    if(type == meRenScale_tt) return "MERENSCALE_TT";
    if(type == meRenScale_z) return "MERENSCALE_Z";
    if(type == meRenScale_w) return "MERENSCALE_W";
    if(type == meRenScale_st) return "MERENSCALE_ST";
    if(type == meRenScale_vv) return "MERENSCALE_VV";
    if(type == meRenScale_ttv) return "MERENSCALE_TTV";
    if(type == meRenScale_ttz) return "MERENSCALE_TTZ";
    if(type == meRenScale_ttw) return "MERENSCALE_TTW";
    if(type == meRenScale_ttg) return "MERENSCALE_TTG";
    if(type == meRenScale_ttdm) return "MERENSCALE_TTDM";
    if(type == meRenScale_htott_res) return "MERENSCALE_HTOTT_RES";
    if(type == meRenScale_htott_int) return "MERENSCALE_HTOTT_INT";
    if(type == meRenScale) return "MERENSCALE";
    if(type == psScale_ttbb) return "PSSCALE_TTBB";
    if(type == psScale_ttb) return "PSSCALE_TTB";
    if(type == psScale_tt2b) return "PSSCALE_TT2B";
    if(type == psScale_ttcc) return "PSSCALE_TTCC";
    if(type == psScale_ttother) return "PSSCALE_TTOTHER";
    if(type == psScale) return "PSSCALE";
    if(type == psISRScale_ttbb)    return "PSISRSCALE_TTBB";
    if(type == psISRScale_tt2b)    return "PSISRSCALE_TT2B";
    if(type == psISRScale_ttb)     return "PSISRSCALE_TTB";
    if(type == psISRScale_ttcc)    return "PSISRSCALE_TTCC";
    if(type == psISRScale_ttother) return "PSISRSCALE_TTOTHER";
    if(type == psISRScale) return "PSISRSCALE";
    if(type == psFSRScale_ttbb)    return "PSFSRSCALE_TTBB";
    if(type == psFSRScale_tt2b)    return "PSFSRSCALE_TT2B";
    if(type == psFSRScale_ttb)     return "PSFSRSCALE_TTB";
    if(type == psFSRScale_ttcc)    return "PSFSRSCALE_TTCC";
    if(type == psFSRScale_ttother) return "PSFSRSCALE_TTOTHER";
    if(type == psScaleWeight_st) return "PSSCALE_WEIGHT_ST";
    if(type == psScaleWeight) return "PSSCALE_WEIGHT";
    if(type == psFSRScale) return "PSFSRSCALE";
    if(type == scale_ttbb) return "SCALE_TTBB";
    if(type == scale_ttb) return "SCALE_TTB";
    if(type == scale_tt2b) return "SCALE_TT2B";
    if(type == scale_ttcc) return "SCALE_TTCC";
    if(type == scale_ttother) return "SCALE_TTOTHER";
    if(type == scale) return "SCALE";
    if(type == bFrag) return "BFRAG";
    if(type == bFrag_central) return "BFRAG_CENTRAL";
    if(type == bFrag_Peterson) return "BFRAG_PETERSON";
    if(type == bSemilep) return "BSEMILEP";
    if(type == ueTune_ttbb)    return "UETUNE_TTBB";
    if(type == ueTune_tt2b)    return "UETUNE_TT2B";
    if(type == ueTune_ttb)     return "UETUNE_TTB";
    if(type == ueTune_ttcc)    return "UETUNE_TTCC";
    if(type == ueTune_ttother) return "UETUNE_TTOTHER";
    if(type == ueTune) return "UETUNE";
    if(type == powhegv2Herwig) return "POWHEGV2HERWIG";
    if(type == powhegHerwig) return "POWHEGHERWIG";
    if(type == powhegHelac) return "POWHEGHELAC";
    if(type == powhegOpenloops) return "POWHEGOPENLOOPS";
    if(type == powhegv2) return "POWHEGV2";
    if(type == powheg) return "POWHEG";
    if(type == amcatnlofxfx) return "AMCATNLOFXFX";
    if(type == mcatnlo) return "MCATNLO";
    if(type == madgraphmlm) return "MADGRAPHMLM";
    if(type == cp5) return "CP5";
    if(type == perugia11NoCR) return "PERUGIA11NoCR";
    if(type == perugia11) return "PERUGIA11";
    if(type == alphasPdf) return "PDF_ALPHAS";
    if(type == l1prefiring) return "L1PREFIRING";
    if(type == normPdfGg) return "NORMPDFGG";
    if(type == normPdfGq) return "NORMPDFGQ";
    if(type == normPdfQq) return "NORMPDFQQ";
    if(type == normPdfTth) return "NORMPDFTTH";
    if(type == pdf_pca_1) return "PDF_PCA_1";
    if(type == pdf_pca_2) return "PDF_PCA_2";
    if(type == pdf_envelope) return "PDF_ENVELOPE";
    if(type == pdf) return "PDF";
    if(type == unclustered) return "UNCLUSTERED";
    if(type == closure) return "closure";
    if(type == allAvailable) return "allAvailable";
    if(type == all) return "all";
    if(type == undefinedType) return "";
    
    if(type == jetPileupIDapplied) return "jetPileupIDapplied";
    if(type == jetLooseCleaningApplied) return "jetLooseCleaningApplied";
    if(type == met40Cut) return "met40Cut";
    
    std::cerr<<"Error in Systematic::convertType()! Conversion is not implemented\n...break\n"<<std::endl;
    exit(99);
}

TString Systematic::convertTypeString(const TString& type)
{
    return convertType(convertType(type));
}



std::vector<Systematic::Type> Systematic::convertType(const std::vector<TString>& types)
{
    std::vector<Type> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}



std::vector<Systematic::Type> Systematic::convertType(const std::vector<std::string>& types)
{
    std::vector<Type> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}



std::vector<TString> Systematic::convertType(const std::vector<Type>& types)
{
    std::vector<TString> v_type;
    for(const auto& type : types) v_type.push_back(convertType(type));
    return v_type;
}

TString Systematic::getPrintName(const TString& type)
{
    if(type == "BSEMILEP") return "B semi-leptonic BR";
    else if(type == "BTAG") return "b tagging";
    else if(type == "CR_ENVELOPE") return "Color reconnection";
    else if(type == "JES") return "Jet energy scale";
    else if(type == "JER") return "Jet energy resolution";
    else if(type == "L1PREFIRING") return "L1 prefiring";
    else if(type == "LEPTON") return "Lepton reconstruction";
    else if(type == "LUMI") return "Luminosity";
    else if(type == "MATCH") return "ME-PS matching";
    else if(type == "MESCALE_ENVELOPE") return "ME scale";
    else if(type == "MTOP") return "Top mass";
    else if(type == "PDF_ALPHAS") return "PDF #alpha_{s}";
    else if(type == "PDF_ENVELOPE") return "PDF replica";
    else if(type == "PS") return "Parton shower";
    else if(type == "PU") return "Pileup";
    else if(type == "TOP_PT") return "Top p_{T}";
    else if(type == "TRIG") return "Trigger";
    else if(type == "UETUNE") return "Underlying event";
    else if(type == "UNCLUSTERED") return "Unclustered energy";
    else if(type == "XSEC BKG") return "Background cross sections";
    else return type;       // return argument if print name not definied
}








// --------------------- Functions defined in namespace Systematic for Variation -------------------------





Systematic::Variation Systematic::convertVariation(const TString& variation)
{
    if(variation.EndsWith("_UP")) return up;
    if(variation.EndsWith("_DOWN")) return down;
    if(variation.EndsWith("_CENTRAL")) return central;
    //std::cout<<"Warning in Systematic::convertVariation()! Following conversion is not implemented: "
    //         <<variation<<std::endl<<std::endl;
    return undefinedVariation;
}



TString Systematic::convertVariation(const Variation& variation)
{
    if(variation == up) return "_UP";
    if(variation == down) return "_DOWN";
    if(variation == central) return "_CENTRAL";
    if(variation == undefinedVariation) return "";
    std::cerr<<"Error in Systematic::convertVariation()! Conversion is not implemented\n...break\n"<<std::endl;
    exit(99);
}



std::vector<Systematic::Variation> Systematic::convertVariation(const std::vector<TString>& variations)
{
    std::vector<Variation> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}



std::vector<Systematic::Variation> Systematic::convertVariation(const std::vector<std::string>& variations)
{
    std::vector<Variation> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}



std::vector<TString> Systematic::convertVariation(const std::vector<Variation>& variations)
{
    std::vector<TString> v_variation;
    for(const auto& variation : variations) v_variation.push_back(convertVariation(variation));
    return v_variation;
}










// --------------------- Further functions defined in namespace Systematic -------------------------



void Systematic::isValid(const Type& type, const Variation& variation, const int variationNumber)
{

    // Check validity of variationNumber
    if(variationNumber >= 0){
        if( (std::find(centralTypes.begin(), centralTypes.end(), type) == centralTypes.end()) && (std::find(uncorrelatedTypes.begin(), uncorrelatedTypes.end(), type) == uncorrelatedTypes.end())){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation numbers (type, variationNumber): "
                     <<convertType(type)<<" , "<<variationNumber<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    
    // Check validity of variation
    if(variation == undefinedVariation) return;
    else if(variation==up || variation==down){
        if(std::find(upDownTypes.begin(), upDownTypes.end(), type) == upDownTypes.end()){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation (type, variation): "
                     <<convertType(type)<<" , "<<convertVariation(variation)<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    else if(variation == central){
        if(std::find(centralTypes.begin(), centralTypes.end(), type) == centralTypes.end()){
            std::cerr<<"ERROR in Systematic::isValid()! Given type does not allow variation (type, variation): "
                     <<convertType(type)<<" , "<<convertVariation(variation)<<"\n...break\n"<<std::endl;
            exit(7);
        }
    }
    else{
        std::cerr<<"ERROR in Systematic::isValid()! Variation is not defined for validity check: "
                 <<convertVariation(variation)<<"\n...break\n"<<std::endl;
        exit(7);
    }
}




std::vector<Systematic::Systematic> Systematic::allowedSystematicsAnalysis(const std::vector<Type>& allowedTypes)
{
    std::vector<Systematic> result;

    for(const Type& type : allowedTypes){
        // Exclude non-real types
        if(type==all || type==allAvailable) continue;

        if(std::find(centralTypes.begin(), centralTypes.end(), type) != centralTypes.end()){
            // Central types need specific treatment using variation numbers, e.g. PDF variations
            // They require detailed specifications at the place where they are used
            // FIXME: This implementation of 26 hardcoded PDF variations should be made more generic
            if(type == pdf){
                result.push_back(Systematic(type, central, 0));
                for(int id = 1; id <= 26; ++id){
                    result.push_back(Systematic(type, up, id));
                    result.push_back(Systematic(type, down, id));
                }
            }
            else
                result.push_back(Systematic(type, undefinedVariation));
        }
        else if(std::find(upDownTypes.begin(), upDownTypes.end(), type) != upDownTypes.end()){
            // Up/down types need the two variations
            result.push_back(Systematic(type, up));
            result.push_back(Systematic(type, down));
        }
        else{
            // All others have no variations
            result.push_back(Systematic(type, undefinedVariation));
        }
    }

    return result;
}



std::vector<Systematic::Systematic> Systematic::setSystematics(const std::vector<std::string>& systematicNames)
{
    std::vector<Systematic> result;

    for(const auto& name : systematicNames) result.push_back(Systematic(name));

    return result;
}


Systematic::Systematic Systematic::nominalSystematic()
{
    return Systematic(nominal, undefinedVariation);
}



Systematic::Systematic Systematic::undefinedSystematic()
{
    return Systematic(undefinedType, undefinedVariation);
}

TString Systematic::puWeightName(const Systematic & systematic){
    std::cout<<"--- Beginning preparation of pileup weight\n";
    const Type type = systematic.type();
    TString weightName = "pu_weight";
    if(type==pu){
        std::cout<<"Use systematic of type: "<<convertType(type)<<"\n";
        if(systematic.variation() == up) {
            std::cout<<"Apply systematic variation: up\n";
            weightName = "pu_weight_up";
        }
        else if(systematic.variation() == down) {
            std::cout<<"Apply systematic variation: down\n";
            weightName = "pu_weight_down";
        }
        else {
            std::cerr << "ERROR in constructor of pileup weights! Systematic variation is invalid: "
                  << convertVariation(systematic.variation()) << "\n...break\n\n";
            exit(98);
      }
   }
   else std::cout<<"Do not apply systematic variation\n";
   
   return weightName;
}

TString Systematic::prefiringWeightName(const Systematic & systematic){
    std::cout<<"--- Beginning preparation of prefiring weight\n";
    const Type type = systematic.type();
    TString weightName = "prefiring_weight";
    if(type==l1prefiring){
        std::cout<<"Use systematic of type: "<<convertType(type)<<"\n";
        if(systematic.variation() == up) {
            std::cout<<"Apply systematic variation: up\n";
            weightName = "prefiring_weight_up";
        }
        else if(systematic.variation() == down) {
            std::cout<<"Apply systematic variation: down\n";
            weightName = "prefiring_weight_down";
        }
        else {
            std::cerr << "ERROR in constructor of pileup weights! Systematic variation is invalid: "
                  << convertVariation(systematic.variation()) << "\n...break\n\n";
            exit(98);
      }
   }
   else std::cout<<"Do not apply systematic variation\n";
   
   return weightName;
}

TString Systematic::metNameAddition(const Systematic & systematic){
    std::cout<<"--- Beginning preparation of unclustered energy shift\n";
    const Type type = systematic.type();
    TString nameAddition = "";
    if(type==unclustered){
        std::cout<<"Use systematic of type: "<<convertType(type)<<"\n";
        if(systematic.variation() == up) {
            std::cout<<"Apply systematic variation: up\n";
            nameAddition = "_UnclE_up";
        }
        else if(systematic.variation() == down) {
            std::cout<<"Apply systematic variation: down\n";
            nameAddition = "_UnclE_down";
        }
        else {
            std::cerr << "ERROR in constructor of unclustered energy shift! Systematic variation is invalid: "
                  << convertVariation(systematic.variation()) << "\n...break\n\n";
            exit(98);
      }
   }
   else std::cout<<"Do not apply systematic variation\n";
   
   return nameAddition;
}

void Systematic::checkAlternativeSample(const Systematic & systematic, const TString &dsSyst, const TString &dsName){
    std::cout<<"--- Beginning check of alternative samples\n";
    const Type type = systematic.type();
    if(std::find(altSampleTypes.begin(), altSampleTypes.end(), type) != altSampleTypes.end()){
        if (systematic.name() == dsSyst){
            std::cout<<"Use systematic: "<<dsSyst<<"\n";
        }
        else {
            std::cerr << "ERROR in check of alternativ samples! Systematic variation is not valid for sample: "
                    << dsName <<"\n...break\n\n";
            exit(150);
        }
    }
    else if(dsSyst != "Nominal"){
        std::cerr << "ERROR in check of alternativ samples! Alternative sample does not match selected systematic: "
            << dsName <<"\n...break\n\n";
        exit(151);
    }
}

bool Systematic::checkTopPTreweighting(const Systematic & systematic){
    std::cout<<"--- Beginning preparation of top PT reweighting\n";
    const Type type = systematic.type();
    bool apply = true;
    if(type==topPt){
        std::cout<<"Use systematic of type: "<<convertType(type)<<"\n";
        std::cout<<"Top PT reweighting will not be applied\n";
        apply = false;
    }
    else std::cout<<"Do not apply systematic variation\n";
    
    return apply;
}

int Systematic::numberOfWeightTypes(){
    int count = 0;
    for (const Type type : weightTypes){
        if (type == pdf){
            count +=100;
        }
        else{
            if (std::find(upDownTypes.begin(), upDownTypes.end(), type) != upDownTypes.end()){
                count +=2;
            }
            else count++;
        }
    }
    return count;
}

bool Systematic::isCorrelated(const TString& type){
    if(std::find(correlatedTypes.begin(), correlatedTypes.end(), convertType(type)) != correlatedTypes.end()){
        return true;
    }
    else{
        return false;
    }
}
    
// --------------------- Methods of class Systematic in namespace Systematic -------------------------



Systematic::Systematic::Systematic():
type_(undefinedType),
variation_(undefinedVariation),
variationNumber_(-1)
{}



Systematic::Systematic::Systematic(const Type& type, const Variation& variation, const int variationNumber):
type_(type),
variation_(variation),
variationNumber_(variationNumber)
{}



Systematic::Systematic::Systematic(const TString& systematicName):
type_(undefinedType),
variation_(undefinedVariation),
variationNumber_(-1)
{

    TString fragment(systematicName);
    type_ = convertType(systematicName);
    fragment.ReplaceAll(convertType(type_), "");
    variation_ = convertVariation(fragment);
    fragment.ReplaceAll(convertVariation(variation_), "");
    if(fragment != ""){
        fragment.ReplaceAll("_", "");
        int variationNumber = -1;
        std::stringstream stream(fragment.Data());
        if(!(stream>>variationNumber)){
            std::cerr<<"ERROR in constructor of Systematic! Could not fragment systematic name (name --- type, variation, variationNumber): "
                     <<systematicName<<" --- "<<convertType(type_)<<" , "<<convertVariation(variation_)<<" , "<<variationNumber_<<"\n...break\n"<<std::endl;
            exit(8);
        }
        variationNumber_ = variationNumber;
        if(variationNumber_ < 0){
            std::cerr<<"ERROR in constructor of Systematic! Variation numbers must be >=0, but is (systematicName, extracted variationNumber): "
                     <<systematicName<<" , "<<variationNumber_<<"\n...break\n"<<std::endl;
            exit(8);
        }
    }
    isValid(type_, variation_, variationNumber_);
}



TString Systematic::Systematic::name()const
{
    TString result = convertType(type_);
    if(variationNumber_ >= 0){
        std::stringstream stream;
        stream<<"_"<<variationNumber_;
        result.Append(stream.str());
    }
    result.Append(convertVariation(variation_));
    return result;
}
