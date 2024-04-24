//ROOT includes

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TMath.h"

//file streaming includes
#include <fstream>
#include <iostream>


void plot_resolution() {

  double pt;
  double res;

  std::vector<double> pT;
  std::vector<double> pTNew;
  std::vector<double> Res;
  std::vector<double> ResNew;
  std::vector<double> Xerrs;
  std::vector<double> Yerrs;
  std::vector<double> newXerrs;
  std::vector<double> newYerrs;

  std::vector<double> pTAll;
  std::vector<double> pTAllNew;
  std::vector<double> resAll;
  std::vector<double> resAllNew;

  std::vector<double> numEl;

  std::fstream inFile1;
  inFile1.open("PFO_resolution.txt",std::ios::in);
  while ( inFile1 >> pt >> res){
    pTAll.push_back(pt);
    resAll.push_back(res/pt);
    if(res/pt<-1){
      std::cout<<"pt: "<<pt<<", res: "<<res<<std::endl;
    }
  }
  inFile1.close();

  //pop outliers because they definitely came from a coding issue in the performance file
  for(int i=0; i<pTAll.size(); i++){
    if(resAll[i]>1){
      resAll.erase(resAll.begin()+i);
      pTAll.erase(pTAll.begin()+i);
    }
  }

  /* std::fstream inFile2;
  inFile2.open("tauEff_new.txt",std::ios::in);
  while (inFile2 >> pt >> eff){
    pTAllNew.push_back(pt);
    effAllNew.push_back(eff);
  }
  inFile2.close();*/
  
  int bins = 10000;
  double max = *max_element(resAll.begin(), resAll.end());
  std::cout<<max<<std::endl;
  double min = *min_element(resAll.begin(), resAll.end());
  double bin_width = (max-min)/bins;
  std::vector<double> edges;
  for(int m=0; m<bins; m++){
    edges.push_back(min+bin_width*m);
  }
  for(int j=0; j<edges.size()-1; j++){
    std::vector<double> oldSelectedpT; //for calculating standard deviation
    std::vector<double> oldSelectedRes; //for calculating stdev
    std::vector<double> newSelectedpT;
    std::vector<double> newSelectedRes;
    double pTsum=0;
    double ressum=0;
    double newpTsum=0;
    double newRessum=0;

    int onumel=0;
    int nnumel=0;

    for(int i=0; i<resAll.size(); i++){
      if(edges[j]<=resAll[i] && resAll[i]<edges[j+1]){
	oldSelectedRes.push_back(resAll[i]);
        pTsum+=pTAll[i];
	ressum+=resAll[i];
	onumel++;
      }
    }
    if(onumel!=0){
      Res.push_back(ressum/double(onumel)); 
    }
    else{
      Res.push_back(edges[j]+bin_width/2);
    }// this is the x axis now
    numEl.push_back(onumel); //we're basically manually histogramming
    //now get the errors -- use standard deviation
    double xerrsum=0;
    double yerrsum=0;
    for(int t=0; t<oldSelectedpT.size(); t++){
      xerrsum+=std::pow((oldSelectedRes[t]-Res.back()),2);
    }
    if(onumel!=0){
      Yerrs.push_back(std::sqrt(1/double(onumel))); 
      Xerrs.push_back(std::sqrt(xerrsum/onumel));
    }
    else{
      Xerrs.push_back(bin_width);
      Yerrs.push_back(0);
    }
  

  
    /*
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
  newYerrs.push_back(std::sqrt(yerrsum/nnumel)); 
  newXerrs.push_back(std::sqrt(xerrsum/nnumel));
    */
  }
  //now convert the vectors to arrays so tgraph will accept them
  Double_t pT_arr[bins];
  std::copy(pT.begin(), pT.end(), pT_arr);
  Double_t Res_arr[bins];
  std::copy(Res.begin(), Res.end(), Res_arr);
  Double_t Xerrs_arr[bins];
  std::copy(Xerrs.begin(), Xerrs.end(), Xerrs_arr);
  Double_t Yerrs_arr[bins];
  std::copy(Yerrs.begin(), Yerrs.end(), Yerrs_arr);

  Double_t numel_arr[bins];
  std::copy(numEl.begin(), numEl.end(), numel_arr);

  double par[9];

  //make a function for two overlaid Gaussians
  TF1 *f1 = new TF1("f1", "gaus", -0.8, -0.1);
  TF1 *f2 = new TF1("f2", "gaus", -0.05, 0.05);
  TF1 *f3 = new TF1("f3", "gaus", 0.1, 1);
  //TF1 *f4 = new TF1("f4", "gaus", min, -0.8);
  TF1 *fitAll = new TF1("fitAll", "gaus(0)+gaus(3)+gaus(6)",-0.8,0.8);
  fitAll->SetLineColor(kGreen-3);

  TH1F *h1 = new TH1F("#delta pT/pT", "pT Resolution; #delta pT/pT; ", 500,-1,1);
  for(int z=0; z<resAll.size(); z++){
    h1->Fill(resAll[z]);
  }
  TGraphErrors *gr = new TGraphErrors(bins,Res_arr,numel_arr,Xerrs_arr,Yerrs_arr);
  gr->SetMarkerColor(kRed+1);
  gr->SetLineColor(kBlue-4);
  gr->SetMarkerStyle(22);
  gr->SetTitle("All PFO pT Resolution");
  //gr->GetHistogram()->SetMaximum(1.);
  //gr->GetHistogram()->SetMinimum(0.);
  gr->GetXaxis()->SetTitle("#Delta pT/pT");
  gr->GetYaxis()->SetTitle("N");
  gr->GetXaxis()->SetRangeUser(-5,5);

  //h1->Fit("gaus","V","E1",-5,5);
  h1->Fit(fitAll,"V","E1",-0.8,0.8);
  h1->Fit(f1, "R");
  h1->Fit(f2, "R+");
  h1->Fit(f3, "R+");
  //h1->Fit(f4, "R+");

  f1->GetParameters(&par[0]);
  f2->GetParameters(&par[3]);
  f3->GetParameters(&par[6]);

  fitAll->SetParameters(par);
  h1->Fit(fitAll, "R+");
  


  //TF1 *fit = h1->GetFunction("fitAll");

  //TGraphErrors *ngr = new TGraphErrors(bins, npT_arr, nEff_arr, nXerrs_arr, nYerrs_arr);
  //ngr->SetMarkerColor(2);
  //ngr->SetLineColor(46);
  //ngr->SetMarkerStyle(20);


  auto canvas = new TCanvas();
  gStyle->SetOptFit(1110);



  gr->Draw("APE");
  h1->Draw("EP");
  //ngr->Draw("EP");

  auto legend = new TLegend(0.1,0.8,0.4,0.9);
  //legend->AddEntry(gr, "Pion pT Resolution (PFOs)");
  //legend->AddEntry(ngr, "Modified DeDuper");
  //legend->AddEntry(fit);
  //legend->Draw();

  canvas->Print("pfoRes.png");

  std::cout<<min<<std::endl;



}

int main(){
  plot_resolution();

  return(0);
}
