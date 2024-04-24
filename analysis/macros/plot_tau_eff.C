//ROOT includes

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TEfficiency.h"

//file streaming includes
#include <fstream>
#include <iostream>

void plot_tau_eff(){

  TH1 *_h_all = new TH1F("All", ";;", 30,0,270);
  TH1 *_h_good = new TH1F("Pass", ";;", 30,0,270);
  TH1 *_h_all_theta = new TH1F("AllTheta", ";;", 30,0,3.14);
  TH1 *_h_good_theta = new TH1F("PassTheta", ";;", 30,0,3.14);

  bool comp=true;


  TH1*_h_allc = new TH1F("AllComp", ";;", 30,0,270);
  TH1*_h_passc = new TH1F("PassComp", ";;", 30,0,270);
  TH1*_h_allcth = new TH1F("AllCTh", ";;", 30,0,3.14);
  TH1*_h_passcth = new TH1F("PassCTh", ";;", 30,0,3.14);
  

  double pt;
  double eff;
  double theta;

  std::fstream inFile;
  inFile.open("tauEff_August.txt");
  while(inFile >> pt >> eff >> theta){
    _h_all->Fill(pt);
    _h_all_theta->Fill(theta);
    if(eff){
      _h_good->Fill(pt);
      _h_good_theta->Fill(theta);
    }
  }
  inFile.close();

  if(comp){
    std::fstream cfile;
    cfile.open("tauEff_August_COMP.txt");
    while(cfile >> pt >> eff >> theta){
      _h_allc->Fill(pt);
      _h_allcth->Fill(theta);
      if(eff){
	_h_passc->Fill(pt);
	_h_passcth->Fill(theta);
      }
    }
    cfile.close();
  }





  TEfficiency *Eff = new TEfficiency(*_h_good, *_h_all);
  Eff->SetTitle("Tau PFO Efficiency; pT [GeV]; #epsilon");
  TEfficiency *compEff = new TEfficiency(*_h_passc, *_h_allc);

  Int_t n=1;
  Double_t x[n],y[n];
  x[0]=0;
  y[0]=0;

  auto canvas = new TCanvas;
  auto leg1 = new TLegend(0.1,0.1,0.4,0.2);
  auto gr1 = new TGraph(n,x,y);
  auto axis = gr1->GetXaxis();
  axis->SetLimits(0.,290.);
  gr1->GetHistogram()->SetMaximum(1.0);
  gr1->GetHistogram()->SetMinimum(0.0);
  gr1->SetTitle("Tau PFO Efficiency; pT [GeV]; #epsilon");
  gr1->Draw("AEP");
  Eff->SetMarkerColor(4);
  Eff->SetLineColor(2);
  Eff->SetMarkerStyle(20);
  leg1->AddEntry(Eff, "Raised Energy Cut");
  Eff->Draw("EZP SAME");
  compEff->SetMarkerColor(kAzure+7);
  compEff->SetLineColor(kViolet-5);
  compEff->SetMarkerStyle(20);
  if(comp) {
    leg1->AddEntry(compEff, "Flat #theta sample");
    leg1->Draw();
    compEff->Draw("EZP SAME");
  }
  canvas->Print("TauEff_compplots.png");

  TEfficiency *EffTh = new TEfficiency(*_h_good_theta, *_h_all_theta);
  EffTh->SetTitle("Tau PFO Efficiency; #theta [rad]; #epsilon");
  TEfficiency *thetaCompEff = new TEfficiency(*_h_passcth, *_h_allcth);


  auto thetacan = new TCanvas;
  auto leg2 = new TLegend(0.1,0.1,0.4,0.2);
  EffTh->SetMarkerColor(kCyan);
  EffTh->SetLineColor(kMagenta);
  EffTh->SetMarkerStyle(20);
  leg2->AddEntry(EffTh, "Raised Energy Cut");
  EffTh->Draw("AEZP");
  thetaCompEff->SetMarkerColor(kGreen+2);
  thetaCompEff->SetLineColor(kViolet+8);
  thetaCompEff->SetMarkerStyle(20);
  if(comp){
    leg2->AddEntry(thetaCompEff, "Flat #theta sample");
    leg2->Draw();
    thetaCompEff->Draw("EZP SAME");
  }
  thetacan->Print("TauEff_compthetaplots.png");

}

int main(){
  plot_tau_eff();
  return(0);
}
