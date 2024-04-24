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


/*void plot_efficiency() {

  double pt;
  double eff;

  std::vector<double> pT;
  std::vector<double> pTNew;
  std::vector<double> Eff;
  std::vector<double> EffNew;
  std::vector<double> Xerrs;
  std::vector<double> Yerrs;
  std::vector<double> newXerrs;
  std::vector<double> newYerrs;

  std::vector<double> YSTDs;
  std::vector<double> newYSTDs;

  std::vector<double> pTAll;
  std::vector<double> pTAllNew;
  std::vector<double> effAll;
  std::vector<double> effAllNew;

  std::fstream inFile1;
  inFile1.open("piPFOs.txt",std::ios::in);
  while ( inFile1 >> pt >> eff){
  pTAll.push_back(pt);
  effAll.push_back(eff);
  }
  inFile1.close();

  std::fstream inFile2;
  inFile2.open("piPFOs_new.txt",std::ios::in);
  while (inFile2 >> pt >> eff){
  pTAllNew.push_back(pt);
  effAllNew.push_back(eff);
  }
  inFile2.close();
  
  int bins = 20;
  double max = *max_element(pTAll.begin(), pTAll.end());
  double min = *min_element(pTAll.begin(), pTAll.end());
  double bin_width = (max-min)/bins;
  std::vector<double> edges;
  for(int i=0; i<bins; i++){
  edges.push_back(bin_width*i);
  }
  for(int j=0; j<edges.size()-1; j++){
  std::vector<double> oldSelectedpT; //for calculating standard deviation
  std::vector<double> oldSelectedEff; //for calculating stdev
  std::vector<double> newSelectedpT;
  std::vector<double> newSelectedEff;
  double pTsum=0;
  double effsum=0;
  double newpTsum=0;
  double newEffsum=0;

  int onumel=0;
  int nnumel=0;

  for(int i=0; i<pTAll.size(); i++){
  if(edges[j]<=pTAll[i] && pTAll[i]<edges[j+1]){
  oldSelectedpT.push_back(pTAll[i]);
  oldSelectedEff.push_back(effAll[i]);
  pTsum+=pTAll[i];
  effsum+=effAll[i];
  onumel++;
  }
  }
  pT.push_back(pTsum/double(onumel));
  Eff.push_back(effsum/double(onumel));
  //now get the errors -- use standard deviation
  double xerrsum=0;
  double yerrsum=0;
  for(int t=0; t<oldSelectedpT.size(); t++){
  xerrsum+=std::pow((oldSelectedpT[t]-pT.back()),2);
  yerrsum+=std::pow((oldSelectedEff[t]-Eff.back()),2);
  }
  YSTDs.push_back(std::sqrt(yerrsum/onumel));
  Yerrs.push_back(std::sqrt(1/double(onumel)));
  Xerrs.push_back(std::sqrt(xerrsum/onumel));
  

  
    
  for(int i=0; i<pTAllNew.size(); i++){
  if(edges[j]<=pTAllNew[i] && pTAllNew[i]<edges[j+1]){
  newSelectedpT.push_back(pTAllNew[i]);
  newSelectedEff.push_back(effAllNew[i]);
  newpTsum+=pTAllNew[i];
  newEffsum+=effAllNew[i];
  nnumel++;
  }
  }
  pTNew.push_back(newpTsum/double(nnumel));
  EffNew.push_back(newEffsum/double(nnumel));
  //now get the errors -- use standard deviation
  xerrsum=0;
  yerrsum=0;
  for(int t=0; t<newSelectedpT.size(); t++){
  xerrsum+=std::pow((newSelectedpT[t]-pTNew.back()),2);
  yerrsum+=std::pow((newSelectedEff[t]-EffNew.back()),2);
  }
  newYSTDs.push_back(std::sqrt(yerrsum/nnumel));
  newYerrs.push_back(std::sqrt(1/double(nnumel)));
  newXerrs.push_back(std::sqrt(xerrsum/nnumel));
    
  }
  //now convert the vectors to arrays so tgraph will accept them
  Double_t pT_arr[bins];
  std::copy(pT.begin(), pT.end(), pT_arr);
  Double_t Eff_arr[bins];
  std::copy(Eff.begin(), Eff.end(), Eff_arr);
  Double_t Xerrs_arr[bins];
  std::copy(Xerrs.begin(), Xerrs.end(), Xerrs_arr);
  Double_t Yerrs_arr[bins];
  std::copy(Yerrs.begin(), Yerrs.end(), Yerrs_arr);
  Double_t YSTDs_arr[bins];
  std::copy(YSTDs.begin(), YSTDs.end(), YSTDs_arr);

  
  Double_t npT_arr[bins];
  std::copy(pTNew.begin(), pTNew.end(), npT_arr);
  Double_t nEff_arr[bins];
  std::copy(EffNew.begin(), EffNew.end(), nEff_arr);
  Double_t nXerrs_arr[bins];
  std::copy(newXerrs.begin(), newXerrs.end(), nXerrs_arr);
  Double_t nYerrs_arr[bins];
  std::copy(newYerrs.begin(), newYerrs.end(), nYerrs_arr);
  Double_t nYSTDs_arr[bins];
  std::copy(newYSTDs.begin(), newYSTDs.end(), nYSTDs_arr);

  //TCanvas *c1 = new TCanvas("c1", "Efficiency", 200,100,700,500);


  // c1->SetGrid();
  //c1->GetFrame()->SetBorderSize(12);
  //int *xstd;
  // xstd=(int*)malloc(bins * sizeof(int));
  // for(int i=0; i<bins; i++){
  //  &xstd[i]=0;
  // }

  const double* xstd=0;

  TGraphErrors *std = new TGraphErrors(bins, pT_arr, Eff_arr,xstd, YSTDs_arr);
  std->SetMarkerColor(4);
  std->SetMarkerStyle(1);
  std->SetLineColorAlpha(kYellow-10,0.999999);
  std->SetLineWidth(5);
  std->GetHistogram()->SetMaximum(1.);
  std->GetHistogram()->SetMinimum(0.);
  std->GetYaxis()->SetTitle("#epsilon");
  std->GetXaxis()->SetTitle("pT [GeV]");
  std->GetXaxis()->SetRangeUser(1.,145);
  std->SetTitle("Pion PFO Efficiency");

  TGraphErrors *nstd = new TGraphErrors(bins,npT_arr, nEff_arr, xstd, nYSTDs_arr);
  nstd->SetMarkerColor(4);
  nstd->SetMarkerStyle(1);
  nstd->SetLineColorAlpha(kYellow-9,0.999999);
  nstd->SetLineWidth(5);
  

  TGraphErrors *gr = new TGraphErrors(bins,pT_arr,Eff_arr,Xerrs_arr,Yerrs_arr);
  gr->SetMarkerColor(4);
  gr->SetLineColor(2);
  gr->SetMarkerStyle(20);
  gr->SetTitle("Pion PFO Efficiency");
  gr->GetHistogram()->SetMaximum(1.);
  gr->GetHistogram()->SetMinimum(0.);
  gr->GetYaxis()->SetTitle("#epsilon");
  gr->GetXaxis()->SetTitle("pT [GeV]");
  gr->GetXaxis()->SetRangeUser(1.,145);

  TGraphErrors *ngr = new TGraphErrors(bins, npT_arr, nEff_arr, nXerrs_arr, nYerrs_arr);
  ngr->SetMarkerColor(2);
  ngr->SetLineColor(kGreen+3);
  ngr->SetMarkerStyle(20);

  auto canvas = new TCanvas();

  std->Draw("AEP");
  nstd->Draw("PE");
  gr->Draw("PE");
  ngr->Draw("PE");
  

  auto legend = new TLegend(0.1,0.8,0.4,0.9);
  legend->AddEntry(gr, "Pion Efficiency (PFOs,Original DeDupe)");
  legend->AddEntry(std, "Standard Deviation");
  legend->AddEntry(ngr, "Modified DeDuper");
  legend->Draw();

  canvas->Print("pfoPiEff_comp.png");


  }*/



//NEW technique with TEfficiency to get more accurate error bars
void plot_efficiency(){

  //declare histograms
  TH1 *_h_all = new TH1F("All", ";;", 40,0,270);
  TH1 *_h_good = new TH1F("Pass", ";;",40,0,270);
  TH1 *_h_ID = new TH1F("ID'd",  ";;", 40,0,270);
  TH1 *_h_trks = new TH1F("HasTracks", ";;", 40,0,270);
  TH1 *_h_charged = new TH1F("Good", ";;", 40,0,270);
  TH1 *_h_ID_trk = new TH1F("GoodID_trkd", ";;", 40, 0, 270);

  TH1 *_h_all1 = new TH1F("All1", ";;", 40,0,270);
  TH1 *_h_good1 = new TH1F("Pass1", ";;", 40,0,270);
  TH1 *_h_ID1 = new TH1F("ID'd1",  ";;", 40,0,270);
  TH1 *_h_trks1 = new TH1F("HasTracks1", ";;", 40,0,270);
  TH1 *_h_charged1 = new TH1F("Good1", ";;", 40,0,270);
  TH1 *_h_ID_trk1 = new TH1F("GoodID_trkd1", ";;", 40, 0, 270);

  TH1 *_h_all2 = new TH1F("All2", ";;", 30,0,3);
  TH1 *_h_good2 = new TH1F("Pass2", ";;", 30,0,3);
  TH1 *_h_ID2 = new TH1F("ID'd2",  ";;", 30,0,3);
  TH1 *_h_trks2 = new TH1F("HasTracks2", ";;", 30,0,3);
  TH1 *_h_charged2 = new TH1F("Good2", ";;", 30,0,3);
  TH1 *_h_ID_trk2 = new TH1F("GoodID_trkd2", ";;", 30, 0, 3);

  TH1 *_h_all3 = new TH1F("All3", ";;", 30,0,3);
  TH1 *_h_good3 = new TH1F("Pass3", ";;", 30,0,3);
  TH1 *_h_ID3 = new TH1F("ID'd3",  ";;", 30,0,3);
  TH1 *_h_trks3 = new TH1F("HasTracks3", ";;", 30,0,3);
  TH1 *_h_charged3 = new TH1F("Good3", ";;", 30,0,3);
  TH1 *_h_ID_trk3 = new TH1F("GoodID_trkd3", ";;", 30, 0, 3);
 
  //read in the file, fill the histograms

  double pt;
  double eff;
  double pdg_e;
  double trks;
  double q;
  double theta;
  double pdg;
  double cls;

  std::vector<int> pdgs;
  int nTrks(0);
  int ncharged(0);
  std::fstream inFile;
  inFile.open("PFO_CONTROL.txt", std::ios::in);
  while(inFile >> pt >> eff >> pdg_e >> trks >> q >> theta >> pdg >> cls){
    //std::cout<<pdg<<std::endl;
    if(abs(pdg)!=211) continue;
    pdgs.push_back(pdg);
    _h_all->Fill(pt);
    _h_all3->Fill(theta);
    if(eff){
      _h_good->Fill(pt);
      _h_good3->Fill(theta);
      if(pdg_e){
	_h_ID->Fill(pt);
	_h_ID3->Fill(theta);
      }
      if(abs(q)==1){//if the truth charge is actually a charged particle
	ncharged++;
	_h_charged->Fill(pt); //the baseline "truth" charged particle histogram (all)
	_h_charged3->Fill(theta);
	if(trks&&cls){
	  _h_trks->Fill(pt); //of the truth charged particles, how many PFO have tracks AND cl associated?
	  _h_trks3->Fill(theta);
	  if(pdg_e){
	    _h_ID_trk->Fill(pt);
	    _h_ID_trk3->Fill(theta);
	  }
	}
      }
    }
  }

  inFile.close();


  std::fstream inFile1;
  inFile1.open("PFO_Z0.txt", std::ios::in);
  while(inFile1 >> pt >> eff >> pdg_e >> trks >> q >> theta >> pdg){
    if(abs(pdg)!=211) continue;
    pdgs.push_back(pdg_e);
    _h_all1->Fill(pt);
    _h_all2->Fill(theta);
    if(eff){
      _h_good1->Fill(pt);
      _h_good2->Fill(theta);
      if(pdg_e){
	_h_ID1->Fill(pt);
	_h_ID2->Fill(theta);
      }
      if(abs(q)==1){ //if the truth charge is actually a charged particle
	_h_charged1->Fill(pt); //the baseline "truth" charged particle histogram (all)
	_h_charged2->Fill(theta);
	if(trks){
	  _h_trks1->Fill(pt); //of the truth charged particles, how many PFO have tracks associated?
	  _h_trks2->Fill(theta);
	  if(pdg_e){
	    _h_ID_trk1->Fill(pt);
	    _h_ID_trk2->Fill(theta);
	  }
	}
      }
    }
  }
    
  inFile1.close();

  //now make the TEfficiency

  TEfficiency *Eff = new TEfficiency(*_h_good, *_h_all);
  Eff->SetTitle("PFO Efficiency; pT [GeV]; #epsilon");
  TEfficiency *compEff = new TEfficiency(*_h_good1, *_h_all1);
  compEff->SetTitle("PFO Efficiency; pT [GeV]; #epsilon");

  TEfficiency *eff_th = new TEfficiency(*_h_good2, *_h_all2);
  eff_th->SetTitle("PFO Efficiency; #theta [rad]; #epsilon");
  TEfficiency *ceff_th = new TEfficiency(*_h_good3, *_h_all3);
  ceff_th->SetTitle("PFO Efficiency; #theta [rad]; #epsilon");

  TEfficiency *pdg_eff = new TEfficiency(*_h_ID, *_h_good);
  pdg_eff->SetTitle("PFO ID Efficiency; pT [GeV]; #epsilon");
  TEfficiency *compPdg_eff = new TEfficiency(*_h_ID1, *_h_good1);
  compPdg_eff->SetTitle("Identification Efficiency; pT [GeV]; #epsilon");

  TEfficiency *pdg_eff_th = new TEfficiency(*_h_ID2, *_h_good2);
  pdg_eff_th->SetTitle("PFO ID Efficiency; #theta [rad]; #epsilon");
  TEfficiency *cpdg_eff_th = new TEfficiency(*_h_ID3, *_h_good3);
  cpdg_eff_th->SetTitle("PFO ID Efficiency; #theta [rad]; #epsilon");

  TEfficiency *trk_eff = new TEfficiency(*_h_trks, *_h_charged);
  trk_eff->SetTitle("PFO Track-Association Efficiency; pT [GeV]; #epsilon");
  TEfficiency *compTrk_eff = new TEfficiency(*_h_trks1, *_h_charged1);
  compTrk_eff->SetTitle("PFO Track-Association Efficiency; pT [GeV]; #epsilon");

  TEfficiency *trk_eff_th = new TEfficiency(*_h_trks2, *_h_charged2);
  trk_eff_th->SetTitle("PFO Track-Association Efficiency; #theta [rad]; #epsilon");
  TEfficiency *ctrk_eff_th = new TEfficiency(*_h_trks3, *_h_charged3);
  ctrk_eff_th->SetTitle("PFO Track-Association Efficiency; #theta [rad]; #epsilon");

  TEfficiency *trk_pdg_eff = new TEfficiency(*_h_ID_trk, *_h_trks);
  trk_pdg_eff->SetTitle("ID Efficiency for PFO with Associated Tracks; pT [GeV]; #epsilon"); 
  TEfficiency *compTrk_pdg_eff = new TEfficiency(*_h_ID_trk1, *_h_trks1);
  compTrk_pdg_eff->SetTitle("ID Efficiency for PFO with Associated Tracks; pT [GeV]; #epsilon");

  TEfficiency *trk_pdg_eff_th = new TEfficiency(*_h_ID_trk2, *_h_trks2);
  trk_pdg_eff_th->SetTitle("ID Efficiency for PFO with Associated Tracks; #theta [rad]; #epsilon");
  TEfficiency *ctrk_pdg_eff_th = new TEfficiency(*_h_ID_trk3, *_h_trks3);
  ctrk_pdg_eff_th->SetTitle("ID Efficiency for PFO with Associated Tracks; #theta [rad]; #epsilon");

  auto legend = new TLegend(0.9,0.8,0.6,0.9);
  legend->AddEntry(Eff, "PFO Creation Efficiency");
  legend->AddEntry(compEff, "PFO Efficiency (flat #theta)");
  //legend->AddEntry(pdg_eff, "Particle ID Efficiency");
  //legend->AddEntry(trk_eff, "Track Association Efficiency");
  
  auto canvas = new TCanvas;

  
  Eff->SetMarkerColor(4);
  Eff->SetLineColor(2);
  Eff->SetMarkerStyle(20);
  ceff_th->SetMarkerColor(4);
  ceff_th->SetLineColor(2);
  ceff_th->SetMarkerStyle(20);
  compEff->SetMarkerColor(kBlue-9);
  compEff->SetLineColor(kRed-2);
  compEff->SetMarkerStyle(20);
  eff_th->SetMarkerColor(kBlue-9);
  eff_th->SetLineColor(kRed-2);
  eff_th->SetMarkerStyle(20);

  pdg_eff->SetMarkerColor(2);
  pdg_eff->SetLineColor(kGreen+3);
  pdg_eff->SetMarkerStyle(20);
  cpdg_eff_th->SetMarkerColor(2);
  cpdg_eff_th->SetLineColor(kGreen+3);
  cpdg_eff_th->SetMarkerStyle(20);
  compPdg_eff->SetMarkerColor(kRed+3);
  compPdg_eff->SetLineColor(kCyan-6);
  compPdg_eff->SetMarkerStyle(20);
  pdg_eff_th->SetMarkerColor(kRed+3);
  pdg_eff_th->SetLineColor(kCyan-6);
  pdg_eff_th->SetMarkerStyle(20);

  trk_eff->SetMarkerColor(kMagenta+3);
  trk_eff->SetMarkerStyle(20);
  trk_eff->SetLineColor(kCyan+2);
  ctrk_eff_th->SetMarkerColor(kMagenta+3);
  ctrk_eff_th->SetMarkerStyle(20);
  ctrk_eff_th->SetLineColor(kCyan+2);
  compTrk_eff->SetMarkerColor(kMagenta-6);
  compTrk_eff->SetMarkerStyle(20);
  compTrk_eff->SetLineColor(kGreen+2);
  trk_eff_th->SetMarkerColor(kMagenta-6);
  trk_eff_th->SetMarkerStyle(20);
  trk_eff_th->SetLineColor(kGreen+2);
  
  trk_pdg_eff->SetMarkerColor(kBlue+2);
  trk_pdg_eff->SetMarkerStyle(20);
  trk_pdg_eff->SetLineColor(kRed-7);
  ctrk_pdg_eff_th->SetMarkerColor(kBlue+2);
  ctrk_pdg_eff_th->SetMarkerStyle(20);
  ctrk_pdg_eff_th->SetLineColor(kRed-7);
  compTrk_pdg_eff->SetMarkerColor(kBlue-8);
  compTrk_pdg_eff->SetMarkerStyle(20);
  compTrk_pdg_eff->SetLineColor(kOrange+10);
  trk_pdg_eff_th->SetMarkerColor(kBlue-8);
  trk_pdg_eff_th->SetMarkerStyle(20);
  trk_pdg_eff_th->SetLineColor(kOrange+10);


  //Eff->SetStatisticOption(TEfficiency::kFFC);
  //Eff->Draw("AEPZ");
  compEff->Draw("AEPZ");
  //pdg_eff->Draw("AEPZ");
  trk_eff->Draw("EPZ SAME");
  //Eff->Draw("EPZ SAME");
  
  //legend->Draw();
  //canvas->Print("eff_control.png");
  canvas->Print("piPFOEff_t.png");

  auto canvas2 = new TCanvas;

  auto legend2 = new TLegend(0.1,0.8,0.4,0.9);
  legend2->AddEntry(pdg_eff, "PFO ID Efficiency");
  legend2->AddEntry(compPdg_eff, "PFO ID Eff (flat #theta)");
  pdg_eff->Draw("ZAEP");
  compPdg_eff->Draw("EPZ SAME");
  legend2->Draw();
  canvas2->Print("pdgEff_new_pi_t.png");

  auto canvas3 = new TCanvas;
  auto legend3 = new TLegend(0.1, 0.8, 0.4, 0.9);
  legend3->AddEntry(trk_eff, "PFO Track Efficiency");
  legend3->AddEntry(compTrk_eff, "PFO Trk Eff (flat #theta)");
  trk_eff->Draw("ZAEP");
  compTrk_eff->Draw("EZP SAME");
  legend3->Draw();
  canvas3->Print("trkEff_new_pi_t.png");

  auto canvas4 = new TCanvas;
  auto legend4 = new TLegend(0.1,0.1,0.4,0.2);
  legend4->AddEntry(trk_pdg_eff, "ID Efficiency for PFO w Trk");
  //legend4->AddEntry(compTrk_pdg_eff, "ID Eff for PFO w Trk (flat #theta)");
  trk_pdg_eff->Draw("ZAEP");
  //compTrk_pdg_eff->Draw("AEZP");
  legend4->Draw();
  canvas4->Print("pdgtrkEff_new_pi.png");

  auto canvas5 = new TCanvas;
  auto legend5 = new TLegend(0.1,0.1,0.4,0.2);
  legend5->AddEntry(eff_th, "PFO Creation Efficiency (flat #theta)");
  //legend5->AddEntry(ceff_th, "PFO Creation Efficiency");
  eff_th->Draw("ZAEP");
  //ceff_th->Draw("PEZ SAME");
  legend5->Draw();
  canvas5->Print("eff_theta_pi_t.png");

  auto canvas6 = new TCanvas;
  auto legend6 = new TLegend(0.1,0.1,0.4,0.2);
  legend6->AddEntry(pdg_eff_th, "PID Efficiency (flat #theta)");
  //legend6->AddEntry(cpdg_eff_th, "PID Efficiency");
  pdg_eff_th->Draw("ZAEP");
  //cpdg_eff_th->Draw("ZEP SAME");
  legend6->Draw();
  canvas6->Print("pdgeff_theta_pi_t.png");

  auto canvas7 = new TCanvas;
  auto legend7 = new TLegend(0.1,0.1,0.4,0.2);
  legend7->AddEntry(trk_eff_th, "Track Association Efficiency (flat #theta)");
  //legend7->AddEntry(ctrk_eff_th, "Track Association Efficiency");
  trk_eff_th->Draw("ZAEP");
  //ctrk_eff_th->Draw("ZEP SAME");
  legend7->Draw();
  canvas7->Print("trkeff_theta_pi_t.png");

  auto canvas8 = new TCanvas;
  auto legend8 = new TLegend(0.1,0.1,0.4,0.2);
  legend8->AddEntry(trk_pdg_eff_th, "ID Eff for PFO w Trk (flat #theta)");
  //legend8->AddEntry(ctrk_pdg_eff_th, "ID Eff for PFO w Trk");
  trk_pdg_eff_th->Draw("ZAEP");
  //ctrk_pdg_eff_th->Draw("PEZ SAME");
  legend8->Draw();
  canvas8->Print("eff_theta_pi_t.png");
  




}

void plot_efficiency(){
  TH1 *_h_all = new TH1F("All", ";;", 20,0,270);
  TH1 *_h_good = new TH1F("Pass", ";;",20,0,270);

  int pdg;
  double pt;
  int eff;

  std::fstream clFile; 
  clFile.open("clAna.txt", std::ios::in);
  while(clFile >> eff >> pt >> pdg){
    if(abs(pdg)!=211) continue;
    _h_all->Fill(pt);
    if(eff){
      _h_good->Fill(pt);
    }
  }

  TEfficiency * clEff = new TEfficiency(*_h_good, *_h_all);
  clEff->SetTitle("Cluster Efficiency; pT [GeV]; #epsilon");
  clEff->SetMarkerColor(4);
  clEff->SetLineColor(2);
  clEff->SetMarkerStyle(20);
  auto canvas = new TCanvas;
  clEff->Draw("ZAEP");
  canvas->Print("cl_pi.txt");
  

}

int main(){
  plot_efficiency();

  return(0);
}
