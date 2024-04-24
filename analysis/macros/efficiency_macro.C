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
  TFile* effFile= TFile::Open("/work/rosep8/TauRecoDev/analysis/pftest.root") ;
  //TList* hlist=(TList*)effFile->Get(Form("CaloAna"));
  TDirectory* dir = (TDirectory*)effFile->Get("eff");

  TH1* _h_pass_pt = (TH1F*)dir->Get<TH1F>("pass, pT");
  TH1* _h_all_pt = (TH1F*)dir->Get<TH1F>("all, pT");

  TH1* _h_pass_theta = (TH1F*)dir->Get<TH1F>("pass, theta");
  TH1* _h_all_theta = (TH1F*)dir->Get<TH1F>("all, theta");

  //TH1* _h_pass_phi = (TH1F*)dir->Get<TH1F>("pass, phi");
  //TH1* _h_all_phi = (TH1F*)dir->Get<TH1F>("all, phi");

  // TH1* _h_pass_z = (TH1F*)dir->Get<TH1F>("pass, z");
  //TH1* _h_all_z = (TH1F*)dir->Get<TH1F>("all, z");

 


  TEfficiency *trks_1 = new TEfficiency(*_h_pass_pt, *_h_all_pt);
  //trks_1->SetTitle("Cluster Efficiency, pions (tau gun); Truth pT [GeV]; #epsilon");

  TGraph *point = new TGraph();
  point->SetPoint(0,50,0);
  point->SetPoint(1,50,1);
  point->SetPoint(2,300,0);
  point->SetTitle("Cluster Efficiency, pions (tau gun) (dR < 1); Truth E [GeV]; #epsilon");


  TCanvas *c_trk = new TCanvas;
  auto leg_trk = new TLegend(0.7,0.1,0.9,0.2);
  //leg_trk->AddEntry(trks_1, "");
  trks_1->SetMarkerColor(4);
  trks_1->SetLineColor(2);
  trks_1->SetMarkerStyle(20);
  point->Draw("PEAZ");
  trks_1->Draw("PEZ SAME");
  //point->Draw("PEZ SAME");
  //leg_trk->Draw();
  c_trk->Print("cleff_taugun_E_100.png");

  TEfficiency *theta = new TEfficiency(*_h_pass_theta,*_h_all_theta);
  //theta->SetTitle("PFO Efficiency, pions (tau gun); #theta; #epsilon");
  
  TGraph *point2 = new TGraph();
  point2->SetPoint(0,-6,0);
  point2->SetPoint(1,0,1);
  point2->SetPoint(2,6,0);
  point2->SetTitle("Cluster Efficiency, pions (tau gun), dR < 1; #eta; #epsilon");


  TCanvas *c_theta = new TCanvas;
  auto leg_theta = new TLegend(0.7,0.1,0.9,0.2);
  theta->SetMarkerColor(4);
  theta->SetLineColor(2);
  theta->SetMarkerStyle(20);
  point2->Draw("PEAZ");
  theta->Draw("PEZ SAME");
  c_theta->Print("cleff_taugun_theta_100.png");
  /*
  TEfficiency *phi = new TEfficiency(*_h_pass_phi, *_h_all_phi);
  phi->SetTitle("PFO Efficiency, pi+/pi- (tau gun); #phi; #epsilon");

  TCanvas *c_phi = new TCanvas;
  auto leg_phi = new TLegend(0.7, 0.1, 0.9, 0.2);
  phi->SetMarkerColor(kViolet);
  phi->SetLineColor(kGreen);
  phi->SetMarkerStyle(20);
  phi->Draw("PEAZ");
  c_phi->Print("trkeff_taugun_phi.png");

  TEfficiency *z = new TEfficiency(*_h_pass_z, *_h_all_z);
  z->SetTitle("PFO Efficiency, pi+/pi- (tau gun); vertex z[mm]; #epsilon");

  TCanvas *c_z = new TCanvas;
  auto leg_z = new TLegend(0.7, 0.1, 0.9, 0.2);
  z->SetMarkerColor(kViolet);
  z->SetLineColor(kGreen);
  z->SetMarkerStyle(20);
  z->Draw("PEAZ");
  c_z->Print("trkeff_taugun_z.png");
  */
  effFile->Close();


}
int main(){
  efficiency_macro();
  return(0);
}
