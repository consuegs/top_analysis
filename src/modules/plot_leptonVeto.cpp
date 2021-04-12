//Script to plot genMET<->Ptnunu resolution for different leptonVetos

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
   io::RootFileSaver saver(TString::Format("plots%.1f.root",cfg.processFraction*100),"plot_leptonVeto");
   io::RootFileReader histReader(TString::Format("histograms_%s.root",cfg.treeVersion.Data()),TString::Format("check_leptonVeto%.1f",cfg.processFraction*100));
   TCanvas can;
   can.cd(); 
            
   TH1F *hist;
   TH1F *hist_v1;
   TH1F *hist_v2;
   TH1F *hist_v3;
   TH1F *hist_v4;
   hist=(TH1F*) histReader.read<TH2F>("baseline/all/Combined/diff(pTnunu-genMET)/TTbar_diLepton");
   hist_v1=(TH1F*) histReader.read<TH2F>("baseline/all/Combined/diff(pTnunu-genMET)_lepVeto/TTbar_diLepton");
   hist_v2=(TH1F*) histReader.read<TH2F>("baseline/all/Combined/diff(pTnunu-genMET)_lepVetoIfaddLeptonInAnyBJet/TTbar_diLepton");
   hist_v3=(TH1F*) histReader.read<TH2F>("baseline/all/Combined/diff(pTnunu-genMET)_VetoAnyBJetInMETdirection/TTbar_diLepton");
   hist_v4=(TH1F*) histReader.read<TH2F>("baseline/all/Combined/diff(pTnunu-genMET)_VetoAnyBJetInMETdirection_addLeptonInJet/TTbar_diLepton");
   
   hist_v1->SetLineColor(kBlue+2);
   hist_v2->SetLineColor(kMagenta+2);
   hist_v3->SetLineColor(kOrange+1);
   hist_v4->SetLineColor(kCyan+2);
   hist->SetStats(0);
   hist->Draw("hist");
   hist_v4->Draw("hist same");
   hist_v2->Draw("hist same");
   hist_v1->Draw("hist same");
   hist_v3->Draw("hist same");
   
   gfx::LegendEntries le;
   le.append(*hist,TString::Format("Nominal(#mu=%.1f #sigma=%.1f evt_frac=%.0f%)",hist->GetMean(),hist->GetRMS(),hist->GetEntries()*100/hist->GetEntries()),"l");
   le.append(*hist_v1,TString::Format("V1(#mu=%.1f #sigma=%.1f evt_frac=%.0f%)",hist_v1->GetMean(),hist_v1->GetRMS(),hist_v1->GetEntries()*100/hist->GetEntries()),"l");
   le.append(*hist_v2,TString::Format("V2(#mu=%.1f #sigma=%.1f evt_frac=%.0f%)",hist_v2->GetMean(),hist_v2->GetRMS(),hist_v2->GetEntries()*100/hist->GetEntries()),"l");
   le.append(*hist_v3,TString::Format("V3(#mu=%.1f #sigma=%.1f evt_frac=%.0f%)",hist_v3->GetMean(),hist_v3->GetRMS(),hist_v3->GetEntries()*100/hist->GetEntries()),"l");
   le.append(*hist_v4,TString::Format("V4(#mu=%.1f #sigma=%.1f evt_frac=%.0f%)",hist_v4->GetMean(),hist_v4->GetRMS(),hist_v4->GetEntries()*100/hist->GetEntries()),"l");
   
   TLegend leg=le.buildLegend(.55,.75,1-1.5*gPad->GetRightMargin(),-1,1);
   leg.Draw();
      
   TString cat_label="all";
   TLatex label=gfx::cornerLabel(cat_label,1);
   label.Draw();
   can.RedrawAxis();
   saver.save(can,"diff(pTnunu-genMET)_compareVetoes",true,true);
}
