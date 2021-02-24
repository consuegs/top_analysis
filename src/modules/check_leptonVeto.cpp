//Script to check for adapted lepton veto based on separating events with nonPrompt neutrinos
#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>

Config const &cfg=Config::get();

extern "C"
void run()
{
   // ~std::vector<TString> vsDatasubsets({"TT_TuneCUETP8M2T4_13TeV-powheg-pythia8_merged"});
   std::vector<TString> vsDatasubsets({"TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8"});
   
   hist::Histograms<TH1F> hs(vsDatasubsets);    //Define histograms in the following
   hist::Histograms<TH2F> hs2D(vsDatasubsets);    //Define histograms in the following
   for(TString base:{"baseline/","baseline_met120/"}){
      for(TString channel:{"ee/","mumu/","emu/","all/"}){
         for (TString selection:{"NonPromptNeutrinos/","Clean/","Combined/"}){
            hs.addHist(base+channel+selection+"N_addMuon", ";N_{add #mu};EventsBIN" ,3,-0.5,2.5);
            hs.addHist(base+channel+selection+"N_addElectron", ";N_{add e};EventsBIN" ,3,-0.5,2.5);
            hs.addHist(base+channel+selection+"N_addLepton", ";N_{add lepton};EventsBIN" ,3,-0.5,2.5);
            
            hs.addHist(base+channel+selection+"dPhi_met_bJet", ";#|Delta #phi(met,bJet)|;EventsBIN" ,100,0,4);
            
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_lepVeto", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_lepVetoPt40", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAllignedBJetMet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAllignedBJetMet_lepVeto", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAllignedBJetMet_addLeptonInBJet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_lepVetoIfaddLeptonInBJet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_lepVetoIfaddLeptonInAnyBJet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyBJetInMETdirection", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyJetInMETdirection", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyBJetInMETdirection_addLepton", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyJetInMETdirection_addLepton", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyBJetInMETdirection_addLeptonInJet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            hs.addHist(base+channel+selection+"diff(pTnunu-genMET)_VetoAnyJetInMETdirection_addLeptonInJet", ";p_{T}^{#nu#nu}-GenMET (GeV);EventsBIN" ,150,-30,30);
            
            hs2D.addHist(base+channel+selection+"2D_Puppi_vs_genMET", ";p_{T}^{miss} (GeV);GenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
            hs2D.addHist(base+channel+selection+"2D_pTnunu_vs_genMET", ";p_{T}^{#nu#nu} (GeV);GenMET (GeV);EventsBIN" ,100,0,600,100,0,600);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_pT_bJet1", ";|p_{T}^{#nu#nu}-GenMET| (GeV);p_{T}^{bJet1} (GeV);EventsBIN" ,50,0,100,50,0,300);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_pT_bJet1_withAddLep", ";|p_{T}^{#nu#nu}-GenMET| (GeV);p_{T}^{bJet1} (GeV);EventsBIN" ,50,0,100,50,0,300);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_pT_bJet1+AddLep", ";|p_{T}^{#nu#nu}-GenMET| (GeV);p_{T}^{bJet1+addLepton} (GeV);EventsBIN" ,50,0,100,50,0,300);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_ratio_pT_bJet1_AddLep", ";|p_{T}^{#nu#nu}-GenMET| (GeV);p_{T}^{addLepton}/p_{T}^{bJet1};EventsBIN" ,50,0,100,50,0,1);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_dPhi_met_bJet1", ";|p_{T}^{#nu#nu}-GenMET| (GeV);#|Delta #phi(p_{T}^{miss},bJet);EventsBIN" ,50,0,100,50,0,4);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_dPhi_met_bJet1_withAddLep", ";|p_{T}^{#nu#nu}-GenMET| (GeV);#|Delta #phi(p_{T}^{miss},bJet);EventsBIN" ,50,0,100,50,0,4);
            hs2D.addHist(base+channel+selection+"2D_diff(pTnunu-genMET)_vs_METsig", ";|p_{T}^{#nu#nu}-GenMET| (GeV);MET_sig;EventsBIN" ,50,0,100,50,0,100);
            
            for (TString addLep:{"","addMuon/","addElectron/"}){
               hs.addHist(base+channel+selection+addLep+"pT_addLepton1", ";p_{T}^{addLepton1} (GeV);EventsBIN" ,50,0,150);
               hs.addHist(base+channel+selection+addLep+"eta_addLepton1", ";|#eta^{addLepton1}|;EventsBIN" ,50,0,2.5);
               hs.addHist(base+channel+selection+addLep+"dR_addLepton1_bJet", ";#Delta R(addLepton1,bJet);EventsBIN" ,100,0,4);
               hs.addHist(base+channel+selection+addLep+"dPhi_addLepton1_bJet", ";#|Delta #phi(addLepton1,bJet)|;EventsBIN" ,100,0,4);
               
               hs2D.addHist(base+channel+selection+addLep+"2D_diff(pTnunu-genMET)_vs_pT_addLepton1", ";|p_{T}^{#nu#nu}-GenMET| (GeV);p_{T}^{addLepton1} (GeV);EventsBIN" ,50,0,100,50,0,150);
               hs2D.addHist(base+channel+selection+addLep+"2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet", ";|p_{T}^{#nu#nu}-GenMET| (GeV);#Delta R(addLepton1,bJet);EventsBIN" ,50,0,100,50,0,4);
            }
         }
      }
   }
   
   TProfile2D Profile_clean=TProfile2D("",";p_{T}^{addLepton1} (GeV);#Delta R(addLepton1,bJet);mean |p_{T}^{#nu#nu}-GenMET| (GeV)" ,30,0,150,30,0,4);
   TProfile2D Profile_nonPromptNu=TProfile2D("",";p_{T}^{addLepton1} (GeV);#Delta R(addLepton1,bJet);mean |p_{T}^{#nu#nu}-GenMET| (GeV)" ,30,0,150,30,0,4);
   TProfile2D Profile_combined=TProfile2D("",";p_{T}^{addLepton1} (GeV);#Delta R(addLepton1,bJet);mean |p_{T}^{#nu#nu}-GenMET| (GeV)" ,30,0,150,30,0,4);

   
   auto const dss = cfg.datasets.getDatasubset(vsDatasubsets[0]);
   TFile file(dss.getPath(),"read");
   if (file.IsZombie()) {
      return;
   }
   io::log * ("Processing '"+dss.datasetName+"' ");
;
   
   hs.setCurrentSample(dss.name);
   hs2D.setCurrentSample(dss.name);
   
   //Lumi weight for current sample
   float lumi_weight=dss.xsec/float(dss.Ngen)*cfg.lumi;

   TTreeReader reader(cfg.treeName, &file);
   TTreeReaderValue<float> w_pu(reader, "pu_weight");
   TTreeReaderValue<UInt_t> runNo(reader, "runNo");
   TTreeReaderValue<UInt_t> lumNo(reader, "lumNo");
   TTreeReaderValue<ULong64_t> evtNo(reader, "evtNo");
   TTreeReaderValue<Char_t> w_mc(reader, "mc_weight");
   TTreeReaderValue<std::vector<float>> w_pdf(reader, "pdf_weights");
   TTreeReaderValue<std::vector<tree::Muon>>     muons    (reader, "muons");
   TTreeReaderValue<std::vector<tree::Electron>> electrons_add(reader, "electrons_add");
   TTreeReaderValue<std::vector<tree::Muon>>     muons_add    (reader, "muons_add");
   TTreeReaderValue<std::vector<tree::Electron>> electrons(reader, "electrons");
   TTreeReaderValue<std::vector<tree::Jet>>      jets     (reader, "jets");
   TTreeReaderValue<std::vector<tree::GenParticle>> genParticles(reader, "genParticles");
   TTreeReaderValue<std::vector<tree::IntermediateGenParticle>> intermediateGenParticles(reader, "intermediateGenParticles");     
   TTreeReaderValue<std::vector<tree::Particle>> triggerObjects(reader, "hltEG165HE10Filter");
   TTreeReaderValue<tree::MET> MET(reader, "met");
   TTreeReaderValue<tree::MET> GENMET(reader, "met_gen");
   TTreeReaderValue<tree::MET> MET_JESu(reader, "met_JESu");
   TTreeReaderValue<tree::MET> MET_JESd(reader, "met_JESd");
   TTreeReaderValue<tree::MET> MET_Puppi(reader, "metPuppi");
   TTreeReaderValue<tree::MET> MET_NoHF(reader, "metNoHF");
   TTreeReaderValue<tree::MET> MET_Calo(reader, "metCalo");
   TTreeReaderValue<tree::MET> MET_Raw(reader, "met_raw");
   TTreeReaderValue<int> n_Interactions(reader, "true_nPV");
   TTreeReaderValue<float> HTgen(reader, "genHt");
   TTreeReaderValue<bool> is_ee   (reader, "ee");
   TTreeReaderValue<bool> is_emu   (reader, "emu");
   TTreeReaderValue<bool> is_mumu   (reader, "mumu");
   TTreeReaderValue<float> mll   (reader, "mll");
   TTreeReaderValue<float> mt2   (reader, "MT2");
   TTreeReaderValue<float> genMT2   (reader, "genMT2");
   TTreeReaderValue<float> genMT2neutrino   (reader, "genMT2neutrino");
   TTreeReaderValue<float> sf_lep1(reader, "lepton1SF");
   TTreeReaderValue<float> sf_lep2(reader, "lepton2SF");
   TTreeReaderValue<TLorentzVector> genTop(reader, "pseudoTop");
   TTreeReaderValue<TLorentzVector> genAntiTop(reader, "pseudoAntiTop");
   TTreeReaderValue<TLorentzVector> genLepton(reader, "pseudoLepton");
   TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "pseudoAntiLepton");
   TTreeReaderValue<TLorentzVector> genTau(reader, "pseudoTau");
   TTreeReaderValue<TLorentzVector> genAntiTau(reader, "pseudoAntiTau");
   TTreeReaderValue<int> genLeptonPdgId(reader, "pseudoLeptonPdgId");
   TTreeReaderValue<int> genAntiLeptonPdgId(reader, "pseudoAntiLeptonPdgId");
   TTreeReaderValue<TLorentzVector> genB(reader, "pseudoBJet");
   TTreeReaderValue<TLorentzVector> genAntiB(reader, "pseudoAntiBJet");
   TTreeReaderValue<TLorentzVector> genNeutrino(reader, "pseudoNeutrino");
   TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "pseudoAntiNeutrino");
   TTreeReaderValue<TLorentzVector> genWMinus(reader, "pseudoWMinus");
   TTreeReaderValue<TLorentzVector> genWPlus(reader, "pseudoWPlus");
   TTreeReaderValue<int> genDecayMode_pseudo(reader, "ttbarPseudoDecayMode");
   // ~TTreeReaderValue<TLorentzVector> genTop(reader, "genTop");    //Alternative genParticles, which are not based on RIVET
   // ~TTreeReaderValue<TLorentzVector> genAntiTop(reader, "genAntiTop");
   // ~TTreeReaderValue<TLorentzVector> genLepton(reader, "genLepton");
   // ~TTreeReaderValue<TLorentzVector> genAntiLepton(reader, "genAntiLepton");
   // ~TTreeReaderValue<TLorentzVector> genTau(reader, "genTau");
   // ~TTreeReaderValue<TLorentzVector> genAntiTau(reader, "genAntiTau");
   // ~TTreeReaderValue<int> genLeptonPdgId(reader, "genLeptonPdgId");
   // ~TTreeReaderValue<int> genAntiLeptonPdgId(reader, "genAntiLeptonPdgId");
   // ~TTreeReaderValue<TLorentzVector> genB(reader, "genB");
   // ~TTreeReaderValue<TLorentzVector> genAntiB(reader, "genAntiB");
   // ~TTreeReaderValue<TLorentzVector> genNeutrino(reader, "genNeutrino");
   // ~TTreeReaderValue<TLorentzVector> genAntiNeutrino(reader, "genAntiNeutrino");
   // ~TTreeReaderValue<TLorentzVector> genWMinus(reader, "genWMinus");
   // ~TTreeReaderValue<TLorentzVector> genWPlus(reader, "genWPlus");
   TTreeReaderValue<int> genDecayMode(reader, "ttbarDecayMode");
   

   int iEv=0;
   int processEvents=cfg.processFraction*dss.entries;
   while (reader.Next()){
      iEv++;
      if (iEv>processEvents) break;
      if (iEv%(std::max(processEvents/10,1))==0){
         io::log*".";
         io::log.flush();
      }
      
      bool leptonVeto=false;
      bool lepVetoPt40=false;
      bool VetoAllignedBJetMet=false;
      bool VetoAllignedBJetMet_lepVeto=false;
      bool VetoAllignedBJetMet_addLeptonInBJet=false;
      bool lepVetoIfaddLeptonInBJet=false;
      bool lepVetoIfaddLeptonInAnyBJet=false;
      bool VetoAnyBJetInMETdirection=false;
      bool VetoAnyJetInMETdirection=false;
      bool VetoAnyBJetInMETdirection_addLepton=false;
      bool VetoAnyJetInMETdirection_addLepton=false;
      bool VetoAnyBJetInMETdirection_addLeptonInJet=false;
      bool VetoAnyJetInMETdirection_addLeptonInJet=false;
      
      //Do not use tau events in signal sample
      if (*genDecayMode>3) continue;
      
      float fEventWeight=*w_pu * *w_mc * *sf_lep1 * *sf_lep2;     //Set event weight also taking lepton scale factors into account
      hs.setFillWeight(fEventWeight);
      
      float const met=MET->p.Pt();
      float const met_puppi=MET_Puppi->p.Pt();
      float const genMet=GENMET->p.Pt();
      
      bool rec_selection=true;
      bool pseudo_selection=true;
      
      //Baseline selection (separation into ee, emu, mumu already done at TreeWriter)
      std::vector<bool> channel={*is_ee,*is_mumu,*is_emu};
      TLorentzVector p_l1;
      TLorentzVector p_l2;
      int flavor_l1=0;  //1 for electron and 2 for muon
      int flavor_l2=0;
      bool muonLead=true; //Boolean for emu channel
      TString cat="";
      
      rec_selection=selection::diLeptonSelection(*electrons,*muons,channel,p_l1,p_l2,flavor_l1,flavor_l2,cat,muonLead);
            
      std::vector<tree::Jet> cjets;
      std::vector<tree::Jet> BJets;
      std::vector<bool> ttbarSelection=selection::ttbarSelection(p_l1,p_l2,met_puppi,channel,*jets,cjets,BJets);
      if(!std::all_of(ttbarSelection.begin(), ttbarSelection.end(), [](bool v) { return v; })) rec_selection=false;
      
      if (*genDecayMode_pseudo==0) pseudo_selection=false; //pseudo baseline selection
            
      if(!rec_selection || !pseudo_selection) continue;  //fill the following histograms only with events selected by the reco baseline selection and genLevel selection
      
      // Get pT of Neutrino Pair
      TLorentzVector neutrinoPair(0,0,0,0);
      neutrinoPair=(*genNeutrino)+(*genAntiNeutrino);
      
      //Get number of (non)prompt neutrinos in event
      int NpromptNeutrinos=0;
      int NnonpromptNeutrinos=0;
      for (auto const &genParticle : *genParticles){
         if(abs(genParticle.pdgId)==12 || abs(genParticle.pdgId)==14 || abs(genParticle.pdgId)==16){
            if(genParticle.isPrompt) NpromptNeutrinos++;
            else {
               NnonpromptNeutrinos++;
            }
         }
      }
                
      TString selection="Clean";
      if (NnonpromptNeutrinos>0) selection="NonPromptNeutrinos";
      
      TString path_cat="ee";
      if (*is_emu) path_cat="emu";
      else if (*is_mumu) path_cat="mumu";
      
      bool baseline=true;
      for(TString base:{"baseline/","baseline_met120/"}){
         
         if (!baseline && met_puppi<120) continue;
      
         hs2D.fill(base+path_cat+"/"+selection+"/2D_Puppi_vs_genMET",met_puppi,genMet);
         hs2D.fill(base+path_cat+"/"+selection+"/2D_pTnunu_vs_genMET",neutrinoPair.Pt(),genMet);
         hs.fill(base+path_cat+"/"+selection+"/N_addMuon",muons_add->size());
         hs.fill(base+path_cat+"/"+selection+"/N_addElectron",electrons_add->size());
         hs.fill(base+path_cat+"/"+selection+"/N_addLepton",muons_add->size()+electrons_add->size());
         hs.fill(base+path_cat+"/"+selection+"/diff(pTnunu-genMET)",(neutrinoPair.Pt()-genMet));
         hs.fill(base+path_cat+"/"+selection+"/dPhi_met_bJet",abs(BJets[0].p.DeltaPhi(MET_Puppi->p)));

         hs2D.fill(base+"all/"+selection+"/2D_Puppi_vs_genMET",met_puppi,genMet);
         hs2D.fill(base+"all/"+selection+"/2D_pTnunu_vs_genMET",neutrinoPair.Pt(),genMet);
         hs.fill(base+"all/"+selection+"/N_addMuon",muons_add->size());
         hs.fill(base+"all/"+selection+"/N_addElectron",electrons_add->size());
         hs.fill(base+"all/"+selection+"/N_addLepton",muons_add->size()+electrons_add->size());
         hs.fill(base+"all/"+selection+"/diff(pTnunu-genMET)",(neutrinoPair.Pt()-genMet));
         hs.fill(base+"all/"+selection+"/dPhi_met_bJet",abs(BJets[0].p.DeltaPhi(MET_Puppi->p)));
         
         hs2D.fill(base+"all/Combined/2D_Puppi_vs_genMET",met_puppi,genMet);
         hs2D.fill(base+"all/Combined/2D_pTnunu_vs_genMET",neutrinoPair.Pt(),genMet);
         hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_pT_bJet1",abs(neutrinoPair.Pt()-genMet),BJets[0].p.Pt());
         hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_dPhi_met_bJet1",abs(neutrinoPair.Pt()-genMet),abs(BJets[0].p.DeltaPhi(MET_Puppi->p)));
         hs.fill(base+"all/Combined/N_addMuon",muons_add->size());
         hs.fill(base+"all/Combined/N_addElectron",electrons_add->size());
         hs.fill(base+"all/Combined/N_addLepton",muons_add->size()+electrons_add->size());
         hs.fill(base+"all/Combined/diff(pTnunu-genMET)",(neutrinoPair.Pt()-genMet));
         hs.fill(base+"all/Combined/dPhi_met_bJet",abs(BJets[0].p.DeltaPhi(MET_Puppi->p)));
         
         TLorentzVector leadAddLepton(0,0,0,0);
         TLorentzVector leadAddMuon(0,0,0,0);
         TLorentzVector leadAddElectron(0,0,0,0);
         if(electrons_add->size()>0 && muons_add->size()>0){
            if((*electrons_add)[0].p.Pt()>(*muons_add)[0].p.Pt()) {
               leadAddLepton=(*electrons_add)[0].p;
               leadAddElectron=(*electrons_add)[0].p;
            }
            else{
               leadAddLepton=(*muons_add)[0].p;
               leadAddMuon=(*muons_add)[0].p;
            }
         }
         else if(electrons_add->size()>0) {
            leadAddLepton=(*electrons_add)[0].p;
            leadAddElectron=(*electrons_add)[0].p;
         }
         else if(muons_add->size()>0){
            leadAddLepton=(*muons_add)[0].p;
            leadAddMuon=(*muons_add)[0].p;
         }
         
         TLorentzVector METalignedBJet(0,0,0,0);
         for(auto bJet: BJets){
            if(abs(bJet.p.DeltaPhi(MET_Puppi->p))<0.8){
               METalignedBJet=bJet.p;
               VetoAnyBJetInMETdirection=true;
               break;
            }
         }
         
         TLorentzVector METalignedJet(0,0,0,0);
         for(auto jet: cjets){
            if(abs(jet.p.DeltaPhi(MET_Puppi->p))<0.8){
               METalignedJet=jet.p;
               VetoAnyJetInMETdirection=true;
               break;
            }
         }
         
         if(abs(BJets[0].p.DeltaPhi(MET_Puppi->p))<0.4 || (3.14-abs(BJets[0].p.DeltaPhi(MET_Puppi->p)))<0.4) VetoAllignedBJetMet=true;
                  
         if(electrons_add->size()>0 || muons_add->size()>0){
            leptonVeto=true;
            hs.fill(base+"all/"+selection+"/pT_addLepton1",leadAddLepton.Pt());
            hs.fill(base+"all/"+selection+"/eta_addLepton1",leadAddLepton.Eta());
            hs.fill(base+"all/"+selection+"/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
            hs.fill(base+"all/"+selection+"/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
            hs2D.fill(base+"all/"+selection+"/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
            hs2D.fill(base+"all/"+selection+"/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
            
            hs.fill(base+"all/Combined/pT_addLepton1",leadAddLepton.Pt());
            hs.fill(base+"all/Combined/eta_addLepton1",leadAddLepton.Eta());
            hs.fill(base+"all/Combined/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
            hs.fill(base+"all/Combined/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
            hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
            hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
            hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_METsig",abs(neutrinoPair.Pt()-genMet),MET_Puppi->sig);
            
            if(NnonpromptNeutrinos>0) Profile_nonPromptNu.Fill(leadAddLepton.Pt(),leadAddLepton.DeltaR(BJets[0].p),abs(neutrinoPair.Pt()-genMet));
            else Profile_clean.Fill(leadAddLepton.Pt(),leadAddLepton.DeltaR(BJets[0].p),abs(neutrinoPair.Pt()-genMet));
            Profile_combined.Fill(leadAddLepton.Pt(),leadAddLepton.DeltaR(BJets[0].p),abs(neutrinoPair.Pt()-genMet));
            
            if(leadAddLepton.DeltaR(BJets[0].p)<0.4){
               hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_pT_bJet1_withAddLep",abs(neutrinoPair.Pt()-genMet),BJets[0].p.Pt());
               hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_pT_bJet1+AddLep",abs(neutrinoPair.Pt()-genMet),(BJets[0].p+leadAddLepton).Pt());
               hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_ratio_pT_bJet1_AddLep",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt()/BJets[0].p.Pt());
               hs2D.fill(base+"all/Combined/2D_diff(pTnunu-genMET)_vs_dPhi_met_bJet1_withAddLep",abs(neutrinoPair.Pt()-genMet),abs(BJets[0].p.DeltaPhi(MET_Puppi->p)));
               
               if(abs(BJets[0].p.DeltaPhi(MET_Puppi->p))<0.4 || (3.14-abs(BJets[0].p.DeltaPhi(MET_Puppi->p)))<0.4) VetoAllignedBJetMet_addLeptonInBJet=true;
            }
            
            if(leadAddLepton.Pt()>30) lepVetoPt40=true;
            
            for(auto bJet: BJets){
               if(leadAddLepton.DeltaR(bJet.p)<0.4) lepVetoIfaddLeptonInAnyBJet=true;
            }
            
            if(leadAddLepton.DeltaR(BJets[0].p)<0.4) lepVetoIfaddLeptonInBJet=true;
            
            
            if(VetoAnyBJetInMETdirection && METalignedBJet.DeltaR(leadAddLepton)<0.4) VetoAnyBJetInMETdirection_addLeptonInJet=true;
            if(VetoAnyJetInMETdirection && METalignedJet.DeltaR(leadAddLepton)<0.4) VetoAnyJetInMETdirection_addLeptonInJet=true;
            
            if(electrons_add->size()>0){
               hs.fill(base+"all/"+selection+"/addElectron/pT_addLepton1",leadAddLepton.Pt());
               hs.fill(base+"all/"+selection+"/addElectron/eta_addLepton1",leadAddLepton.Eta());
               hs.fill(base+"all/"+selection+"/addElectron/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
               hs.fill(base+"all/"+selection+"/addElectron/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
               hs2D.fill(base+"all/"+selection+"/addElectron/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
               hs2D.fill(base+"all/"+selection+"/addElectron/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
               
               hs.fill(base+"all/Combined/addElectron/pT_addLepton1",leadAddLepton.Pt());
               hs.fill(base+"all/Combined/addElectron/eta_addLepton1",leadAddLepton.Eta());
               hs.fill(base+"all/Combined/addElectron/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
               hs.fill(base+"all/Combined/addElectron/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
               hs2D.fill(base+"all/Combined/addElectron/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
               hs2D.fill(base+"all/Combined/addElectron/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
            }
            if(muons_add->size()>0){
               hs.fill(base+"all/"+selection+"/addMuon/pT_addLepton1",leadAddLepton.Pt());
               hs.fill(base+"all/"+selection+"/addMuon/eta_addLepton1",leadAddLepton.Eta());
               hs.fill(base+"all/"+selection+"/addMuon/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
               hs.fill(base+"all/"+selection+"/addMuon/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
               hs2D.fill(base+"all/"+selection+"/addMuon/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
               hs2D.fill(base+"all/"+selection+"/addMuon/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
               
               hs.fill(base+"all/Combined/addMuon/pT_addLepton1",leadAddLepton.Pt());
               hs.fill(base+"all/Combined/addMuon/eta_addLepton1",leadAddLepton.Eta());
               hs.fill(base+"all/Combined/addMuon/dR_addLepton1_bJet",leadAddLepton.DeltaR(BJets[0].p));
               hs.fill(base+"all/Combined/addMuon/dPhi_addLepton1_bJet",abs(leadAddLepton.DeltaPhi(BJets[0].p)));
               hs2D.fill(base+"all/Combined/addMuon/2D_diff(pTnunu-genMET)_vs_pT_addLepton1",abs(neutrinoPair.Pt()-genMet),leadAddLepton.Pt());
               hs2D.fill(base+"all/Combined/addMuon/2D_diff(pTnunu-genMET)_vs_dR_addLepton1_bJet",abs(neutrinoPair.Pt()-genMet),leadAddLepton.DeltaR(BJets[0].p));
            }
            
         }
         
         VetoAllignedBJetMet_lepVeto=leptonVeto || VetoAllignedBJetMet;
         VetoAnyBJetInMETdirection_addLepton=VetoAnyBJetInMETdirection && leptonVeto;
         VetoAnyJetInMETdirection_addLepton=VetoAnyJetInMETdirection && leptonVeto;
         
         if(!leptonVeto) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_lepVeto",(neutrinoPair.Pt()-genMet));
         if(!lepVetoPt40) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_lepVetoPt40",(neutrinoPair.Pt()-genMet));
         if(!VetoAllignedBJetMet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAllignedBJetMet",(neutrinoPair.Pt()-genMet));
         if(!VetoAllignedBJetMet_lepVeto) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAllignedBJetMet_lepVeto",(neutrinoPair.Pt()-genMet));
         if(!VetoAllignedBJetMet_addLeptonInBJet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAllignedBJetMet_addLeptonInBJet",(neutrinoPair.Pt()-genMet));
         if(!lepVetoIfaddLeptonInBJet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_lepVetoIfaddLeptonInBJet",(neutrinoPair.Pt()-genMet));
         if(!lepVetoIfaddLeptonInAnyBJet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_lepVetoIfaddLeptonInAnyBJet",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyBJetInMETdirection) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyBJetInMETdirection",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyJetInMETdirection) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyJetInMETdirection",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyBJetInMETdirection_addLepton) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyBJetInMETdirection_addLepton",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyJetInMETdirection_addLepton) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyJetInMETdirection_addLepton",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyBJetInMETdirection_addLeptonInJet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyBJetInMETdirection_addLeptonInJet",(neutrinoPair.Pt()-genMet));
         if(!VetoAnyJetInMETdirection_addLeptonInJet) hs.fill(base+"all/Combined/diff(pTnunu-genMET)_VetoAnyJetInMETdirection_addLeptonInJet",(neutrinoPair.Pt()-genMet));
         
         baseline=false;
      }
      
            
   }// evt loop
   io::log<<"";
   
   Profile_nonPromptNu.SaveAs("Profile_nonPromptNu.root");
   Profile_clean.SaveAs("Profile_clean.root");
   Profile_combined.SaveAs("Profile_combined.root");
      
   hs.scaleLumi();
   hs.mergeOverflow();
   hs.normHists();
   hs2D.scaleLumi();
   hs2D.mergeOverflow();
   file.Close();
   
   std::vector<TString> samplesToCombine={"TTbar_diLepton"};
   hs.combineFromSubsamples(samplesToCombine);
   hs2D.combineFromSubsamples(samplesToCombine);
   
   // Save histograms
   io::RootFileSaver saver_hist(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("check_leptonVeto%.1f",cfg.processFraction*100),false);
   hs.saveHistograms(saver_hist,samplesToCombine);
   hs2D.saveHistograms(saver_hist,samplesToCombine);
   
}
