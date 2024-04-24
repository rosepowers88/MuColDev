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

void plot_trkeff(){

  //declare hists
  TH1 *_h_all = new TH1F("All", ";;", 30,0,120);
  TH1 *_h_pass = new TH1F("Pass", ";;", 30,0,120);

  //now for theta
  TH1 *_h_all_th = new TH1F("All_th", ";;", 30,0,3);
  TH1 *_h_pass_th = new TH1F("Pass_th", ";;", 30,0,3);

  //read in the file and fill the hists

  double pt;
  double eff;
  double pdg;
  double theta;

  std::fstream inFile;
  inFile.open("trkEff_new.txt", std::ios::in);
  while(inFile >> pt >> eff >> pdg >> theta){
    if(abs(pdg)==211){
      _h_all->Fill(pt);
      _h_all_th->Fill(theta);
      if(eff){
	_h_pass->Fill(pt);
	_h_pass_th->Fill(theta);
      }
    }
  }

  //make the TEfficiency
  TEfficiency *Eff = new TEfficiency(*_h_pass, *_h_all);
  Eff->SetTitle("Tracking Efficiency (pions); pT [GeV]; #epsilon");
  auto legend = new TLegend(0.1,0.8,0.2,0.9);
  legend->AddEntry(Eff, "Trk Eff");

  TEfficiency *Eff_th = new TEfficiency(*_h_pass_th, *_h_all_th);
  Eff_th->SetTitle("Tracking Efficiency (pions); theta [rad]; #epsilon");
  auto legend2 = new TLegend(0.1,0.8,0.2,0.9);
  legend2->AddEntry(Eff_th, "Trk Eff");

  //style
  auto canvas = new TCanvas;
  Eff->SetMarkerColor(2);
  Eff->SetLineColor(4);
  Eff->SetMarkerStyle(20);

  //draw and save
  Eff->Draw("AEPZ");
  legend->Draw();
  canvas->Print("piTrkEff_pt.png");

  //style
  auto canvas2 = new TCanvas;
  Eff_th->SetMarkerColor(kCyan);
  Eff_th->SetLineColor(kMagenta);
  Eff_th->SetMarkerStyle(20);

  //draw and save
  Eff_th->Draw("AEPZ");
  legend2->Draw();
  canvas2->Print("piTrkEff_theta.png");

}

int main(){
  plot_trkeff();
  return(0);
}
