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

void add_Categories(TString const path, io::RootFileReader const &reader_hist, TH1F &out_hist) {   //Function to add the three different categories
   TH1F *hist;
   for (TString cat:{"ee","emu","mumu"}){
      hist = (TH1F*) reader_hist.read<TH1F>("baseline/"+cat+"/"+path);
      if (cat=="ee") out_hist=(TH1F) *(reader_hist.read<TH1F>("baseline/"+cat+"/"+path));
      else out_hist.Add(hist);
   }
}

extern "C"
void run()
{
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_distributions");
   // ~io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));
   io::RootFileReader histReader(TString::Format("multiHists/histograms_%s.root",cfg.treeVersion.Data()),TString::Format("distributions%.1f",cfg.processFraction*100));

   
   // ~std::vector<TString> samplesToPlot={"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW","T1tttt_1200_800","T1tttt_1500_100","T2tt_650_350","T2tt_850_100","DM_pseudo_50_50","DM_scalar_10_10","DM_scalar_1_200","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
   
   std::vector<TString> samplesToPlot={};
   switch(cfg.year_int){
      case(3): //2018
      samplesToPlot={"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttZ","ttW","DoubleMuon","EGamma","MuonEG","SingleMuon"};
      break;
      case(2): //2017
      samplesToPlot={"TTbar_diLepton","TTbar_diLepton_tau","TTbar_singleLepton","TTbar_hadronic","SingleTop","WJetsToLNu","DrellYan","WW","WZ","ZZ","ttW","ttZ","DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"};
      break;
   }    
   
   std::map<TString,std::vector<TString>> msPresel_vVars={
   {"cutflow/",{"ee","emu","mumu"}}
   };
   
   // ~for(TString selection:{"baseline","baseline_Met200"}){ //Reco 1D Histograms
   for(TString selection:{"baseline"}){ //Reco 1D Histograms
      for(TString channel:{"/ee/","/mumu/","/emu/"}){
         msPresel_vVars.insert(std::pair<TString,std::vector<TString>>(selection+channel,
         {"MET"
         ,"PuppiMET"
         ,"DNN_MET_pT"
         ,"DNN_MET_dPhi_nextLep"
         ,"met1000"
         ,"mLL"
         ,"Lep1_pt"
         ,"Lep2_pt"
         ,"pTsumlep"
         ,"sumpTlep"
         ,"pTbJet"
         ,"Jet1_pt"
         ,"Jet2_pt"
         ,"dPhiMETnearJet"
         ,"dPhiMETleadJet"
         ,"dPhiMETlead2Jet"
         ,"dphi_metNearLep"
         ,"dphi_metNearLep_puppi"
         ,"COSdphi_metNearLep"
         ,"SINdphi_metNearLep"
         ,"dPhiMETbJet"
         ,"dPhiLep1bJet"
         ,"dR_bJetLep1"
         ,"dphi_bJetLep2"
         ,"dphi_bJetnearLep"
         ,"dphi_b1b2"
         ,"dR_b1b2"
         ,"dphi_metLep1"
         ,"dphi_metLep2"
         ,"dphi_metLepsum"
         ,"dPhiLep1Lep2"
         ,"dR_Lep1Lep2"
         ,"nJets"
         ,"nBjets"
         ,"MT2"
         ,"MT"
         ,"mt_MetLep2"
         ,"mt_MetNextLep"
         ,"conMt_Lep1Lep2"
         ,"ST"
         ,"HT"
         ,"sum_STHT"
         ,"sum_mlb"
         ,"METunc_Puppi"
         ,"n_Interactions"
         ,"Lep1_flavor"
         ,"Lep2_flavor"
         ,"Lep1_phi"
         ,"Lep2_phi"
         ,"Lep1_eta"
         ,"Lep2_eta"
         ,"Lep1_E"
         ,"Lep2_E"
         ,"Jet1_phi"
         ,"Jet2_phi"
         ,"Jet1_eta"
         ,"Jet2_eta"
         ,"Jet1_E"
         ,"Jet2_E"
         ,"dPhiMETfarJet"
         ,"dPhiJet1Jet2"
         ,"METsig"
         ,"MHT"
         ,"looseLeptonVeto"
         ,"dPhiMETnearJet_Puppi"
         ,"dPhiMETfarJet_Puppi"
         ,"dPhiMETleadJet_Puppi"
         ,"dPhiMETlead2Jet_Puppi"
         ,"dPhiMETbJet_Puppi"
         ,"dPhiLep1Jet1"
         ,"PFMET_phi"
         ,"PuppiMET_phi"
         ,"CaloMET"
         ,"CaloMET_phi"
         ,"vecsum_pT_allJet"
         ,"vecsum_pT_l1l2_allJet"
         ,"mass_l1l2_allJet"
         ,"ratio_vecsumpTlep_vecsumpTjet"
         ,"mjj"}));
      }
   }
   
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
            if (tempHist->GetNbinsX()>25) tempHist->Rebin(2);
            if (sVar=="dphi_metNearLep") tempHist->Rebin(4);
            hs.addFilledHist(loc,sSample,*(tempHist));
         }
         if(sPresel.Contains("cutflow")){
            TH1F* tempHist=histReader.read<TH1F>(sPresel+"ee/"+sSample);
            tempHist->Add(histReader.read<TH1F>(sPresel+"emu/"+sSample));
            tempHist->Add(histReader.read<TH1F>(sPresel+"mumu/"+sSample));
            hs.addFilledHist(sPresel+"all",sSample,*(tempHist));
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
   // ~hs.combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
   switch(cfg.year_int){
      case(3): //2018
      hs.combineSamples("data",{"DoubleMuon","EGamma","MuonEG","SingleMuon"});
      break;
      case(2): //2017
      hs.combineSamples("data",{"DoubleMuon","DoubleEG","MuonEG","SingleMuon","SingleElectron"});
      break;
   }
   hs.combineSamples("MC",{"TTbar_diLepton","tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttZ","ttW"});
   // ~hs.combineSamples("MC_withCUETP8M2",{"TTbar_diLepton_CUETP8M2","TTbar_hadronic","TTbar_singleLepton","Diboson","SingleTop","WJetsToLNu","DrellYan","ttZ","ttW"});
   hs.combineSamples("SM bkg.",{"tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttZ","ttW"});
   std::map<const TString,Color_t> colormap = {{"Diboson",kCyan-8},{"ttW/Z",kGreen-7},{"tt other",kRed-9}};
   
   
   TCanvas can;
   can.SetLogy();
   gfx::SplitCan sp_can;
   sp_can.pU_.SetLogy();
   for (auto &sPresel_vVars:msPresel_vVars){
   TString const &sPresel=sPresel_vVars.first;
      if(sPresel.Contains("cutflow")) sPresel_vVars.second.push_back("all");
      for (TString sVar:sPresel_vVars.second){
         sp_can.pU_.cd();
         TString loc;
         loc=sPresel+sVar;
         THStack st_mc=hs.getStack(loc,{"ttW/Z","WJetsToLNu","Diboson","DrellYan","SingleTop","tt other","TTbar_diLepton"},colormap);
         gfx::LegendEntries le=hs.getLegendEntries();
         TString cat;
         if (loc.Contains("ee")) cat="ee";
         else if (loc.Contains("emu")) cat="e#mu";
         else if (loc.Contains("mumu")) cat="#mu#mu";
         else if (loc.Contains("all")) cat="all";
         if (sPresel.Contains("Met200")) cat+="  p_{T}^{miss}>200 GeV";
         TLatex label=gfx::cornerLabel(cat,1);
         if (sVar.Contains("phi")){
            st_mc.SetMinimum(1);
            st_mc.SetMaximum(1e6);
         }
         st_mc.SetMinimum(1);
         st_mc.SetMaximum(1e3*st_mc.GetMaximum());
         st_mc.Draw();
         if(sPresel.Contains("cutflow")) st_mc.GetXaxis()->SetRangeUser(0.5,6.5);
         
         auto hist_data = hs.getHistogram(loc,{"data"});
         hist_data->SetLineColor(kBlack);
         hist_data->SetMarkerSize(0.5);
         if(!(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep"))) {
            hist_data->Draw("same");
            le.append(*hist_data,"data","lep");
         }
         
         /*
         auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
         for (auto const &h: hists) {
            h->Draw("same hist");
            // ~h->SetLineWidth(4);
         }
         le+=hs.getLegendEntries();
         auto BSM_legend=hs.getLegendEntries();
         */
         auto hists_SM=hs.getHistograms(loc,{"TTbar_diLepton"});
         TH1F axis=*(hists_SM[0]);
         axis.SetStats(0);
         axis.Draw("same axis");
         TLegend leg=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         leg.Draw();
         label.Draw();
         
         sp_can.pL_.cd();
         TH1F ratio=hist::getRatio(*hist_data,st_mc,"data/MC",hist::ONLY1);
         TH1F ratio_mc=hist::getRatio(st_mc,st_mc,"data/MC",hist::ONLY1);
         if(sPresel.Contains("cutflow")){
            ratio_mc.GetXaxis()->SetBinLabel(1,"DiLepton");
            ratio_mc.GetXaxis()->SetBinLabel(2,"mll");
            ratio_mc.GetXaxis()->SetBinLabel(3,"jets");
            ratio_mc.GetXaxis()->SetBinLabel(4,"met");
            ratio_mc.GetXaxis()->SetBinLabel(5,"btag");
            ratio_mc.GetXaxis()->SetBinLabel(6,"ScaleFactors");
            // ~ratio_mc.GetXaxis()->SetBinLabel(7,"(addLepton veto)");
            ratio_mc.GetXaxis()->SetRangeUser(0.5,6.5);
         }
         ratio_mc.GetYaxis()->SetTitleOffset(0.45);
         ratio_mc.SetStats(0);
         ratio.SetLineColor(kBlack);
         // ~ratio_mc.SetMaximum(1.04);
         // ~ratio_mc.SetMinimum(0.9);
         ratio_mc.SetMaximum(1.25);
         ratio_mc.SetMinimum(0.75);
         ratio_mc.SetFillColor(kGray);
         ratio_mc.SetMarkerSize(0);
         ratio_mc.Draw("e2");
         // ~ratio.SetMarkerSize(4);
         if(!(sVar.Contains("MET") || sVar.Contains("met1000")|| sVar.Contains("dphi_metNearLep")))ratio.Draw("pe1 same");
         
         gfx::LegendEntries le_low;
         le_low.append(ratio_mc,"stat.","f");
         TLegend leg_low=le_low.buildLegend(.2,.8,0.5,0.95,2);
         leg_low.Draw();
         
         saver.save(sp_can,loc,false,true);
         // ~saver.save(can,loc);
         
         
         //normalized distributions
         can.cd();
         auto ttbar_hist=hs.getHistogram(loc,"TTbar_diLepton");
         ttbar_hist->Scale(1.0/(ttbar_hist->Integral()));
         ttbar_hist->SetFillColor(ttbar_hist->GetLineColor());
         // ~ttbar_hist->SetFillStyle(1001);
         auto SMbkg_hist=hs.getHistogram(loc,"SM bkg.");
         SMbkg_hist->Scale(1.0/(SMbkg_hist->Integral()));
         SMbkg_hist->SetFillColor(kGray);
         // ~SMbkg_hist->SetFillStyle(1001);
         THStack st_norm;
         st_norm.Add(SMbkg_hist);
         st_norm.Add(ttbar_hist);
         axis.SetStats(0);
         axis.SetMaximum(1.0);
         axis.SetMinimum(1e-3);
         axis.GetYaxis()->SetTitle("normalized distribution");
         axis.Draw("axis");
         // ~st_norm.Draw("hist same");
         ttbar_hist->SetMarkerSize(0);
         SMbkg_hist->SetMarkerSize(0);
         ttbar_hist->Draw("hist e same");
         SMbkg_hist->Draw("hist e same");
         /*
         for (auto const &h: hists) {
            h->Scale(1.0/(h->Integral()));
            h->SetMarkerSize(0);
            h->Draw("same hist e");
         }
         */
         axis.Draw("axis same");
         le.clear();
         le.prepend(*ttbar_hist,"tt ll","l");
         le.append(*SMbkg_hist,"SM bkg.","l");
         // ~le+=BSM_legend;
         TLegend leg2=le.buildLegend(.42,.7,1-(gPad->GetRightMargin()+0.02),-1,2);
         leg2.Draw();
         label.Draw();
         saver.save(can,"normalized/"+loc,true,false,true);
         
      }
   }
   
   for (TString cat:{"ee","emu","mumu"}){    //Get the number of events per category
      TH1F* mc_total=hs.getHistogram("baseline/"+cat+"/MET","MC");
      std::cout<<"----------------"<<cat<<"-----------------------"<<std::endl;
      for (TString sample:{"TTbar_diLepton","tt other","Diboson","SingleTop","WJetsToLNu","DrellYan","ttW/Z","MC","data"}){
         TH1F* temp_hist=hs.getHistogram("baseline/"+cat+"/MET",sample);
         std::cout<<std::fixed<<sample<<"&"<<temp_hist->Integral()<<"&"<<std::setprecision(1)<<temp_hist->Integral()/mc_total->Integral()*100<<"\\\\"<<std::endl;
         // ~std::cout<<sample<<"   "<<temp_hist->Integral()<<"   "<<temp_hist->Integral()/mc_total->Integral()*100<<std::endl;
      }
   }
   
   //Plot 2D SignalRegion Plot (reco level)
   for (TString channel:{"ee","mumu","emu","all"}){
      TString loc="baseline/"+channel+"/2d_MetVSdPhiMetNearLep_Puppi/";
      THStack st_mc=hs.getStack(loc,{"ttW/Z","WJetsToLNu","Diboson","DrellYan","SingleTop","tt other","TTbar_diLepton"},colormap);
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
      
      /*
      auto hists=hs.getHistograms(loc,{"T1tttt_1200_800","T2tt_650_350","DM_scalar_1_200"});
      for (auto const &h: hists) {
         // ~h->SetLineWidth(4);
         h->Draw("same hist");
      }
      le+=hs.getLegendEntries();
      */
      
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
   
      saver.save(can,loc,true,false,true);
   }
}
