//ROOT includes

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"

//file streaming includes
#include <fstream>
#include <iostream>

void plot_dupes(){


  double hits;
  double dupes;

  std::vector<float> oldHitAxis;
  std::vector<float> oldDupeAxis;
  std::vector<float> oldXerrs;
  std::vector<float> oldYerrs;


  std::vector<float> oldHitsAll;
  std::vector<float> oldDupesAll;

  std::vector<float> newHitAxis;
  std::vector<float> newDupeAxis;
  std::vector<float> newXerrs;
  std::vector<float> newYerrs;


  std::vector<float> newHitsAll;
  std::vector<float> newDupesAll;

  std::fstream inFile;
  inFile.open("oldDupeTrks.txt",std::ios::in);
  while ( inFile >> dupes >> hits){
    oldHitsAll.push_back(hits);
    oldDupesAll.push_back(dupes);
  }
  inFile.close();

  std::fstream inFile2;
  inFile2.open("newDupeTrks.txt", std::ios::in);
  while( inFile2 >>dupes >> hits){
    newHitsAll.push_back(hits);
    newDupesAll.push_back(dupes);
  }
  inFile2.close();
  

  int bins = 20;
  double max = *max_element(oldHitsAll.begin(), oldHitsAll.end());
  double min = *min_element(oldHitsAll.begin(), oldHitsAll.end());
  double bin_width = (max-min)/bins;
  std::vector<double> edges;
  for(int i=0; i<bins; i++){
    edges.push_back(bin_width*i+2);
  }

  for(int j=0; j<edges.size()-1; j++){
    std::vector<double> oldSelectedHits; //for calculating standard deviation
    std::vector<double> oldSelectedDupes; //for calculating stdev
    std::vector<double> newSelectedHits; //for calculating standard deviation
    std::vector<double> newSelectedDupes; //for calculating stdev
    double oldhitsum=0;
    double newhitsum=0;
    double olddupesum=0;
    double newdupesum=0;
    int onumel=0;
    int nnumel=0;
    for(int i=0; i<oldHitsAll.size(); i++){
      if(edges[j]<=oldHitsAll[i] && oldHitsAll[i]<edges[j+1]){
	oldSelectedHits.push_back(oldHitsAll[i]);
	oldSelectedDupes.push_back(oldDupesAll[i]);
	oldhitsum+=oldHitsAll[i];
	olddupesum+=oldDupesAll[i];
	onumel++;
      }
    }
    for(int i=0; i<newHitsAll.size(); i++){
      if(edges[j]<=newHitsAll[i] && newHitsAll[i]<edges[j+1]){
	newSelectedHits.push_back(newHitsAll[i]);
	newSelectedDupes.push_back(newDupesAll[i]);
	newhitsum+=newHitsAll[i];
	newdupesum+=newDupesAll[i];
	nnumel++;
      }
    }
    if(onumel!=0){
      oldHitAxis.push_back(oldhitsum/double(onumel));
      oldDupeAxis.push_back(olddupesum/double(onumel));
    }
    else{
      oldHitAxis.push_back(edges[j]+bin_width);
      oldDupeAxis.push_back(0);
    }
    if(nnumel!=0){
      newHitAxis.push_back(newhitsum/double(nnumel));
      newDupeAxis.push_back(newdupesum/double(nnumel));
    }
    else{
      newHitAxis.push_back(edges[j]+bin_width);
      newDupeAxis.push_back(0);
    }


    //now get the errors -- use standard deviation
    double oxerrsum=0;
    double oyerrsum=0;
    double nxerrsum=0;
    double nyerrsum=0;
    for(int t=0; t<oldSelectedHits.size(); t++){
      oxerrsum+=std::pow((oldSelectedHits[t]-oldHitAxis.back()),2);
      oyerrsum+=std::pow((oldSelectedDupes[t]-oldDupeAxis.back()),2);
    }
    if(onumel!=0){
      oldYerrs.push_back(std::sqrt(oyerrsum/onumel)); 
      oldXerrs.push_back(std::sqrt(oxerrsum/onumel));
    }
    else{
      oldYerrs.push_back(0);
      oldXerrs.push_back(bin_width);
    }
    for(int t=0; t<newSelectedHits.size(); t++){
      nxerrsum+=std::pow((newSelectedHits[t]-newHitAxis.back()),2);
      nyerrsum+=std::pow((newSelectedDupes[t]-newDupeAxis.back()),2);
    }
    if(nnumel!=0){
      newYerrs.push_back(std::sqrt(nyerrsum/nnumel)); 
      newXerrs.push_back(std::sqrt(nxerrsum/nnumel));
    }
    else{
      newYerrs.push_back(0);
      newXerrs.push_back(bin_width);
    }

  }

  Double_t ohits_arr[bins];
  std::copy(oldHitAxis.begin(), oldHitAxis.end(), ohits_arr);
  Double_t odupe_arr[bins];
  std::copy(oldDupeAxis.begin(), oldDupeAxis.end(), odupe_arr);
  Double_t oXerrs_arr[bins];
  std::copy(oldXerrs.begin(), oldXerrs.end(), oXerrs_arr);
  Double_t oYerrs_arr[bins];
  std::copy(oldYerrs.begin(), oldYerrs.end(), oYerrs_arr);

  for(int i=0; i<bins-1; i++){
    //std::cout<<ohits_arr[i]<<std::endl;
  }
  
  Double_t nhits_arr[bins];
  std::copy(newHitAxis.begin(), newHitAxis.end(), nhits_arr);
  Double_t ndupe_arr[bins];
  std::copy(newDupeAxis.begin(), newDupeAxis.end(), ndupe_arr);
  Double_t nXerrs_arr[bins];
  std::copy(newXerrs.begin(), newXerrs.end(), nXerrs_arr);
  Double_t nYerrs_arr[bins];
  std::copy(newYerrs.begin(), newYerrs.end(), nYerrs_arr);

  TGraphErrors *ogr = new TGraphErrors(bins,ohits_arr,odupe_arr,oXerrs_arr,oYerrs_arr);
  ogr->SetMarkerColor(4);
  ogr->SetLineColor(38);
  ogr->SetMarkerStyle(20);
  ogr->SetTitle("nDuplicates as a fn of nHitsPerTrk");
  ogr->GetXaxis()->SetRangeUser(0.5,max+0.5);
  ogr->GetHistogram()->SetMinimum(0.1);
  ogr->GetYaxis()->SetTitle("nDupes");
  ogr->GetXaxis()->SetTitle("nHits");

  TGraphErrors *ngr = new TGraphErrors(bins,nhits_arr,ndupe_arr,nXerrs_arr,nYerrs_arr);
  ngr->SetMarkerColor(2);
  ngr->SetLineColor(46);
  ngr->SetMarkerStyle(20);
  //ngr->SetTitle("New DeDuper");



  auto canvas = new TCanvas();
  ogr->Draw("AEP");
  ngr->Draw("EP");

  auto legend = new TLegend(0.6,0.8,0.9,0.9); // xleft,ybottom , xright, ytop
  legend->AddEntry(ogr, "Old DeDuper");
  legend->AddEntry(ngr, "Modified DeDuper");
  legend->Draw();
  
  canvas->Print("DupesPerHit_100k.png");


}



int main(){
  plot_dupes();

  return(0);
}
