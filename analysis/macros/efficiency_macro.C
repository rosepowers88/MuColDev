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
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"

//file streaming includes
#include <fstream>
#include <iostream>



void efficiency_macro(){
  //get the histograms from the root file we pass in
  TFile* effFile= TFile::Open("taugun_pieff.root") ;
  //TList* hlist=(TList*)effFile->Get(Form("CaloAna"));
  TDirectory* dir = (TDirectory*)effFile->Get("eff");

  TH1* _h_pass_pt = (TH1F*)dir->Get<TH1F>("pass, pT");
  TH1* _h_all_pt = (TH1F*)dir->Get<TH1F>("all, pT");

  TH1* _h_pass_theta = (TH1F*)dir->Get<TH1F>("pass, theta");
  TH1* _h_all_theta = (TH1F*)dir->Get<TH1F>("all, theta");

 


  TEfficiency *trks_1 = new TEfficiency(*_h_pass_pt, *_h_all_pt);
  trks_1->SetTitle("Tracking Efficiency, pi+/pi- (tau gun); pT [GeV]; #epsilon");
   
  TCanvas *c_trk = new TCanvas;
  auto leg_trk = new TLegend(0.7,0.1,0.9,0.2);
  //leg_trk->AddEntry(trks_1, "");
  trks_1->SetMarkerColor(4);
  trks_1->SetLineColor(2);
  trks_1->SetMarkerStyle(20);
  trks_1->Draw("PEAZ");
  //leg_trk->Draw();
  c_trk->Print("trkeff_taugun_pt.png");

  TEfficiency *theta = new TEfficiency(*_h_pass_theta,*_h_all_theta);
  theta->SetTitle("Tracking Efficiency, pi+/pi- (tau gun); #theta; #epsilon");
  
  TCanvas *c_theta = new TCanvas;
  auto leg_theta = new TLegend(0.7,0.1,0.9,0.2);
  theta->SetMarkerColor(kMagenta);
  theta->SetLineColor(kCyan);
  theta->SetMarkerStyle(20);
  theta->Draw("PEAZ");
  c_theta->Print("trkeff_taugun_theta.png");

  effFile->Close();


}
int main(){
  efficiency_macro();
  return(0);
}
