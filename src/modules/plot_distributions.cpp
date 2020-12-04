//Script to (re-)plot distributions from distributions.cpp 

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
#include <iomanip>

Config const &cfg=Config::get();

TH1F HistTrafo_2D(TH2F* const &hist2D, std::vector<float> binedges_x, std::vector<float> binedges_y){   
   *hist2D=hist::rebinned(*hist2D, binedges_x, binedges_y);
   
   int numBins_x=binedges_x.size()-1;
   int numBins_y=binedges_y.size()-1;
   int numBins=numBins_x*numBins_y;
   double binedges_1d[numBins+1];
	binedges_1d[0]=0;
   int phi_bin = 0;
   for (int i=0; i<(numBins); i++)   {
      binedges_1d[i+1] = binedges_x[i%numBins_x+1]+phi_bin*binedges_x[numBins_x];
      if (i%numBins_x==numBins_x-1) phi_bin++;
   }
   
   
   TH1F* tempHist= new TH1F("", "", numBins_x*numBins_y, binedges_1d);
   int binNew = 1;
   for (int j=1; j<=numBins_y; j++){
      for (int i=1; i<=numBins_x; i++){
         tempHist->SetBinContent(binNew, hist2D->GetBinContent(i,j));
         tempHist->SetBinError(binNew, hist2D->GetBinError(i,j));
         binNew++;
      }
   }
   
	return *tempHist;
}

void add_Categories(TString const path, io::RootFileReader const &reader_hist, TH2F &out_hist) {   //Function to add the three different categories
   TH2F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH2F*) reader_hist.read<TH2F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH2F) *(reader_hist.read<TH2F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   // ~io::RootFileReader histReader(TString::Format("../output_v19/histograms_v19.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   
   // ~std::vector<TString> samplesToPlot={"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW",
   std::vector<TString> samplesToPlot={"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan_NLO","WW","WZ","ZZ","ttZ","ttW",
      "T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron","TTbar_diLepton_CUETP8M2"};
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
      {"baseline/ee/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metLeadJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/emu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metLeadJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      {"baseline/mumu/",{"met","met_puppi","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metLeadJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"baseline_Met200/ee/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"baseline_Met200/emu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"baseline_Met200/mumu/",{"met","met1000","mll","pTlep1","pTlep2","pTsumlep","sumpTlep","pTJet1","pTJet2","pTbJet","dphi_metJet","dphi_metBJet","dphi_bJetLep1","dR_bJetLep1","dphi_bJetLep2","dphi_bJetnearLep","dphi_b1b2","dR_b1b2","dphi_metLep1","dphi_metLep2","dphi_Lep1Lep2","nBjets","nJets","mt2","dR_Lep1Lep2","ST","HT","sum_STHT","mt_MetLep1","mt_MetLep2","mt_MetNextLep","sum_mlb","conMt_Lep1Lep2","dphi_metNearLep","COSdphi_metNearLep","SINdphi_metNearLep","dphi_metLepsum"}},
      // ~{"genParticles/ee/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      // ~{"genParticles/emu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      // ~{"genParticles/mumu/",{"dphi_NeutrinoLep","dR_NeutrinoLep","pTtop1","pTtop2"}},
      {"cutflow/",{"ee","emu","mumu"}}
   };
   
   hist::Histograms<TH1F> hs(samplesToPlot);
   
   //Import 1D hists
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      
      for (auto const &sPresel_vVars:msPresel_vVars){
         TString const &sPresel=sPresel_vVars.first;
         for (TString sVar:sPresel_vVars.second){
            TString loc;
            loc=sPresel+sVar;
            TH1F* tempHist=histReader.read<TH1F>(loc+"/"+sSample);
            // ~if (sVar!="nBjets" and sVar!="nJets") tempHist->Rebin(5);
            // ~if (sVar=="dphi_metNearLep") tempHist->Rebin(4);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
      }
   }
   
   //Import 2D hist (trafo to 1D)
   std::vector<float> binedges_x = {0, 20, 40, 60, 80, 100, 120, 140, 160, 195, 230, 400};
   std::vector<float> binedges_y = {0, 0.35, 0.7, 1.05, 1.4, 2.27, 3.14};
   
   for (TString sSample : samplesToPlot){
      hs.setCurrentSample(sSample);
      TH1F sumHist;
      for (TString channel:{"ee","mumu","emu"}){
         TString loc="baseline/"+channel+"/2d_MetVSdPhiMetNearLep_Puppi/";
         TH2F* tempHist=histReader.read<TH2F>(loc+sSample);
         TH1F tempHist1D=HistTrafo_2D(tempHist,binedges_x,binedges_y);
         
         for(int i=1; i<=tempHist1D.GetNbinsX(); i++){
            if(tempHist1D.GetBinContent(i)<0) tempHist1D.SetBinContent(i,0);
         }
         
         hs.addFilledHist(loc,sSample,tempHist1D);
         if(channel=="ee") sumHist=tempHist1D;
         else sumHist.Add(&tempHist1D);
      }
      hs.addFilledHist("baseline/all/2d_MetVSdPhiMetNearLep_Puppi/",sSample,sumHist);
   }
   
   hs.combineSamples("Diboson",{"WW","WZ","ZZ"});
   hs.combineSamples("ttW/Z",{"ttW","ttZ"});
   hs.combineSamples("tt other",{"TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic"});
   hs.combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
   hs.combineSamples("MC",{"TTbar_diLepton","tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW"});
   hs.combineSamples("MC_withCUETP8M2",{"TTbar_diLepton_CUETP8M2","TTbar_hadronic","TTbar_singleLepton","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW"});
   hs.combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW"});
   std::map<const TString,Color_t> colormap = {{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   
   TCanvas can;
   can.SetLogy();
   gfx::SplitCan sp_can;
   sp_can.pU_.SetLogy();
   for (auto const &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      for (TString sVar:sPresel_vVars.second){
         sp_can.pU_.cd();
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs.getStack(loc,{"ttW/Z","WJetsToLNu","Diboson","DrellYan_NLO","SingleTop","tt other","TTbar_diLepton"},colormap);
         gfx::LegendEntries le=hs.getLegendEntries();
         TString cat;
         if (sPresel.Contains("ee/")) cat="ee";
         else if (sPresel.Contains("emu/")) cat="e#mu";
         else if (sPresel.Contains("mumu/")) cat="#mu#mu";
         if (sPresel.Contains("Met200")) cat+="  p_{T}^{miss}>200 GeV";
         TLatex label=gfx::cornerLabel(cat,1);
         if (sVar.Contains("phi")){
            st_mc.SetMinimum(1);
            st_mc.SetMaximum(1e6);
         }
         st_mc.Draw();
         
         auto hist_data = hs.getHistogram(loc,{"data"});
         if(sPresel.Contains("cutflow")) hist_data->Draw("same");
         
         auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
         for (auto const &h: hists) {
            h->Draw("same hist");
            h->SetLineWidth(4);
         }
         le+=hs.getLegendEntries();
         auto BSM_legend=hs.getLegendEntries();
         TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         TH1F axis=*(hists[0]);
         axis.SetStats(0);
         axis.Draw("same axis");
         leg.Draw();
         label.Draw();
         
         sp_can.pL_.cd();
         TH1F ratio=hist::getRatio(*hist_data,st_mc,"data/MC",hist::ONLY1);
         TH1F ratio_mc=hist::getRatio(st_mc,st_mc,"data/MC",hist::ONLY1);
         if(sPresel.Contains("cutflow")){
            ratio_mc.GetXaxis()->SetBinLabel(1,"Ntuple");
            ratio_mc.GetXaxis()->SetBinLabel(2,"Lepton eta,ID,ISO");
            ratio_mc.GetXaxis()->SetBinLabel(3,"Lepton pT");
            ratio_mc.GetXaxis()->SetBinLabel(4,"mll");
            ratio_mc.GetXaxis()->SetBinLabel(5,"jets");
            ratio_mc.GetXaxis()->SetBinLabel(6,"met");
            ratio_mc.GetXaxis()->SetBinLabel(7,"btag");
            ratio_mc.GetXaxis()->SetBinLabel(8,"triggerSF");
            ratio_mc.GetXaxis()->SetBinLabel(9,"(addLepton veto)");
         }
         ratio_mc.GetYaxis()->SetTitleOffset(0.3);
         ratio_mc.SetStats(0);
         ratio.SetLineColor(kBlack);
         ratio_mc.SetMaximum(1.04);
         ratio_mc.SetMinimum(0.85);
         ratio_mc.SetFillColor(kGray);
         ratio_mc.SetMarkerSize(0);
         ratio_mc.Draw("e2");
         ratio.SetMarkerSize(4);
         if(sPresel.Contains("cutflow"))ratio.Draw("pe1 same text");
         
         saver.save(sp_can,loc);
         // ~saver.save(can,loc);
         
         
         //normalized distributions
         can.cd();
         auto ttbar_hist=hs.getHistogram(loc,"TTbar_diLepton");
         ttbar_hist->Scale(1.0/(ttbar_hist->Integral()));
         ttbar_hist->SetFillColor(ttbar_hist->GetLineColor());
         ttbar_hist->SetFillStyle(1001);
         auto SMbkg_hist=hs.getHistogram(loc,"SM bkg.");
         SMbkg_hist->Scale(1.0/(SMbkg_hist->Integral()));
         SMbkg_hist->SetFillColor(kGray);
         SMbkg_hist->SetFillStyle(1001);
         THStack st_norm;
         st_norm.Add(SMbkg_hist);
         st_norm.Add(ttbar_hist);
         axis.SetStats(0);
         axis.SetMaximum(1.0);
         axis.SetMinimum(1e-3);
         axis.GetYaxis()->SetTitle("normalized distribution");
         axis.Draw("axis");
         st_norm.Draw("hist same");
         for (auto const &h: hists) {
            h->Scale(1.0/(h->Integral()));
            h->Draw("same hist");
         }
         axis.Draw("axis same");
         le.clear();
         le.prepend(*ttbar_hist,"tt ll","f");
         le.append(*SMbkg_hist,"SM bkg.","f");
         le+=BSM_legend;
         TLegend leg2=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         leg2.Draw();
         label.Draw();
         saver.save(can,"normalized/"+loc);
         
      }
   }
   
   for (TString cat:{"ee","emu","mumu"}){    //Get the number of events per category
      TH1F* mc_total=hs.getHistogram("baseline/"+cat+"/met","MC");
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      for (TString sample:{"TTbar_diLepton","tt other","Diboson","SingleTop","WJetsToLNu","DrellYan_NLO","ttZ","ttW","data","MC","MC_withCUETP8M2"}){
         TH1F* temp_hist=hs.getHistogram("baseline/"+cat+"/met",sample);
         // ~std::cout<<std::fixed<<sample<<"&"<<temp_hist->Integral()<<"&"<<std::setprecision(1)<<temp_hist->Integral()/mc_total->Integral()*100<<"\\\\"<<std::endl;
         std::cout<<sample<<"   "<<temp_hist->Integral()<<"   "<<temp_hist->Integral()/mc_total->Integral()*100<<std::endl;
      }
   }
   
   //Plot 2D SignalRegion Plot (reco level)
   for (TString channel:{"ee","mumu","emu","all"}){
      TString loc="baseline/"+channel+"/2d_MetVSdPhiMetNearLep_Puppi/";
      THStack st_mc=hs.getStack(loc,{"ttW/Z","WJetsToLNu","Diboson","DrellYan_NLO","SingleTop","tt other","TTbar_diLepton"},colormap);
      st_mc.SetMaximum(st_mc.GetMaximum()*3);
      st_mc.SetTitle(";p_{#scale[.8]{T}}^{#scale[.8]{miss}}(GeV);Events/Bin");
      st_mc.Draw();
      
      gfx::LegendEntries le=hs.getLegendEntries();
      st_mc.GetXaxis()->SetNdivisions(24);
      st_mc.GetXaxis()->ChangeLabel(25,-1,-1,-1,-1,-1," ");
      for (int i=0; i<=5; i++){
         st_mc.GetXaxis()->ChangeLabel(i*4+1,-1,-1,-1,-1,-1,"  ");
         st_mc.GetXaxis()->ChangeLabel(i*4+2,-1,-1,-1,-1,-1,"100");
         st_mc.GetXaxis()->ChangeLabel(i*4+3,-1,-1,-1,-1,-1," ");
         st_mc.GetXaxis()->ChangeLabel(i*4+4,-1,-1,-1,-1,-1,"300");
      }
      
      auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
      for (auto const &h: hists) {
         h->SetLineWidth(4);
         h->Draw("same hist");
      }
      le+=hs.getLegendEntries();
      
      TLine * aline = new TLine();
      TLatex * atext = new TLatex();
      atext->SetTextSize(0.03);
      aline->SetLineWidth(2);
      float upperEnd=st_mc.GetHistogram()->GetMaximum();
      float lowerEnd=st_mc.GetHistogram()->GetMinimum();
      for (int i=1; i<=5; i++){
         aline->DrawLine(i*400,lowerEnd,i*400,upperEnd);
      }
      atext->DrawLatex(100,0.2*upperEnd,"|#Delta#phi|<0.35");
      atext->DrawLatex(475,0.2*upperEnd,"0.35<|#Delta#phi|<0.7");
      atext->DrawLatex(875,0.2*upperEnd,"0.7<|#Delta#phi|<1.05");
      atext->DrawLatex(1275,0.2*upperEnd,"1.05<|#Delta#phi|<1.4");
      atext->DrawLatex(1675,0.2*upperEnd,"1.4<|#Delta#phi|<2.27");
      atext->DrawLatex(2075,0.2*upperEnd,"2.27<|#Delta#phi|");
      TString cat=channel;
      if(channel=="emu") cat="e#mu";
      else if(channel=="mmu") cat="#mu#mu";
      atext->DrawLatex(50,0.3*upperEnd,cat);
      
      TLegend leg=le.buildLegend(gPad->GetLeftMargin(),.9,1-gPad->GetRightMargin(),1-gPad->GetTopMargin(),10);
      leg.SetTextSize(0.025);
      leg.SetBorderSize(2);
      leg.SetShadowColor(0);
      leg.SetFillColor(kWhite);
      leg.SetFillStyle(1001);
      leg.Draw();
   
      saver.save(can,loc);
   }
}
