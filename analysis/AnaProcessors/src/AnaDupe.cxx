#include "AnaProcessors/AnaDupe.hxx"

#include <fstream>
#include <iostream>

#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <EVENT/Track.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHit.h>

#include <marlin/AIDAProcessor.h>
#include <marlin/Statusmonitor.h>

#include <cmath>
#include <algorithm>

//ROOT includes

#include "TCanvas.h"
#include "TROOT.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TFile.h"

AnaDupe aAnaDupe ;

AnaDupe::AnaDupe()
  : Processor("AnaDupe") {
  _description = "AnaDupe analyzes duplicate tracks" ;
  /******************************************
   * Rose Powers (2023)
   * Use this processor for:
   * - Identifying duplicate tracks
   * - Plotting them as a function of nTrkHits
   * 
   *
   *
   *
   *
   *
   *
   **************************************************/

  // register steering parameters: name, description, class-variable, default value


  registerProcessorParameter("MinPt",
			     "Minimum transverse momentum",
			     _minPt,
			     _minPt);
 
  
}

void AnaDupe::init() {
  // Print the initial parameters
  printParameters() ;

  //Basic histograms
  marlin::AIDAProcessor::histogramFactory(this);
  _h_effPFO  = new TH1F("PFO Efficiency", ";nPFO/nMCP", 100,0,1);

}

void AnaDupe::plot(std::vector<int> dupes, std::vector<int> hits, const char* title)
{
  std::vector<float> hitAxis;
  std::vector<float> dupeAxis;
  std::vector<float> Xerrs;
  std::vector<float> Yerrs;
  int bins = 20;
  double max = *max_element(hits.begin(), hits.end());
  double min = *min_element(hits.begin(), hits.end());
  double bin_width = (max-min)/bins;
  std::vector<double> edges;
  for(int i=0; i<bins; i++){
    edges.push_back(bin_width*i);
  }

  for(int j=0; j<edges.size()-1; j++){
    std::vector<double> selectedHits; //for calculating standard deviation
    std::vector<double> selectedDupes; //for calculating stdev
    double hitsum=0;
    double dupesum=0;
    int numel=0;
    for(int i=0; i<hits.size(); i++){
      if(edges[j]<=hits[i] && hits[i]<edges[j+1]){
	selectedHits.push_back(hits[i]);
	selectedDupes.push_back(dupes[i]);
	hitsum+=hits[i];
	dupesum+=dupes[i];
	numel++;
      }
    }

    hitAxis.push_back(hitsum/double(numel));
    dupeAxis.push_back(dupesum/double(numel));
    //now get the errors -- use standard deviation
    double xerrsum=0;
    double yerrsum=0;
    for(int t=0; t<selectedHits.size(); t++){
      xerrsum+=std::pow((selectedHits[t]-hitAxis.back()),2);
      yerrsum+=std::pow((selectedDupes[t]-dupeAxis.back()),2);
    }
    Yerrs.push_back(std::sqrt(yerrsum/numel)); 
    Xerrs.push_back(std::sqrt(xerrsum/numel));
  }

  Double_t hits_arr[bins];
  std::copy(hitAxis.begin(), hitAxis.end(), hits_arr);
  Double_t dupe_arr[bins];
  std::copy(dupeAxis.begin(), dupeAxis.end(), dupe_arr);
  Double_t Xerrs_arr[bins];
  std::copy(Xerrs.begin(), Xerrs.end(), Xerrs_arr);
  Double_t Yerrs_arr[bins];
  std::copy(Yerrs.begin(), Yerrs.end(), Yerrs_arr);

  TFile *fout = new TFile(title, "RECREATE");

  TGraphErrors *gr = new TGraphErrors(bins,hits_arr,dupe_arr,Xerrs_arr,Yerrs_arr);
  gr->SetMarkerColor(2);
  gr->SetMarkerStyle(20);
  gr->SetTitle("Duplicates as fn of TrkHits");
  gr->GetYaxis()->SetTitle("nDupes");
  gr->GetXaxis()->SetTitle("nHits");



  auto canvas = new TCanvas();
  gr->Draw("APE");
  canvas->Print(title);
  canvas->Update();


}


std::vector<std::vector<int>> AnaDupe::findDuplicates( LCCollection *trkCol ) {
  std::vector<std::vector<int>> Dupes;
  // takes a collection of tracks and identifies duplicates
  // identical pT is a tipoff, so is identical Z0 and phi, so check those first
  // then check nSharedHits
  // loop over the collection, starting with the first element, and check for duplicate candidacy
  int ntrk = trkCol->getNumberOfElements();
  for(uint32_t i=0; i<ntrk; i++){
    const EVENT::Track *trk=static_cast<const EVENT::Track*>(trkCol->getElementAt(i));
    float phi = trk->getPhi();
    float Z0 = trk->getZ0();
    const EVENT::TrackerHitVec trkhits = trk->getTrackerHits();
    int nhits = trkhits.size();
    //get the IDs of individual hits so we can cross-check for shared hits
    std::vector<int> hitIDs;
    for(uint32_t l=0; l<nhits; l++){
      const EVENT::TrackerHit* th = static_cast<const EVENT::TrackerHit*>(trkhits[l]);
      int ID = th->getCellID0();
      hitIDs.push_back(ID);
    }
    int nDupes = 0;
    std::vector<int> nHitsSmallest;
    if(ntrk>1){
      for(uint32_t j=0; j<ntrk; j++){
	if(j>i){
	  const EVENT::Track *dupeTrk = static_cast<const EVENT::Track*>(trkCol->getElementAt(j));
	  float dupePhi = dupeTrk->getPhi();
	  float dupeZ0 = dupeTrk->getZ0();
	  const EVENT::TrackerHitVec dupetrkhits = dupeTrk->getTrackerHits();
	  int nDupeHits = dupetrkhits.size();

	  //now check for shared hits
	  int nshared = 0;
	  for(uint32_t z=0; z<nDupeHits; z++){
	    const EVENT::TrackerHit* candidate = static_cast<const EVENT::TrackerHit*>(dupetrkhits[z]);
	    int candID = candidate->getCellID0();
	    if(std::find(hitIDs.begin(), hitIDs.end(), candID) != hitIDs.end()) nshared++;
	  }
      
	  //find track with most hits
	  int bigHits=0; int smallHits=0;
	  if( nDupeHits <= nhits ){
	    bigHits = nhits;
	    smallHits = nDupeHits;
	  }
	  else{
	    bigHits = nDupeHits;
	    smallHits = nhits;
	  }
      
      
	  if(nshared>=0.25*smallHits || dupePhi==phi || dupeZ0 == Z0){
	    nDupes++;
	    nHitsSmallest.push_back(smallHits);
	  }
	}
      }
      //what we want to push back is the smallest ntrkhits, so we minimize the nHitsSmallest vector
      if(nHitsSmallest.size()>0){
	int smallest = *std::min_element(nHitsSmallest.begin(), nHitsSmallest.end());
	std::vector<int> row;
	row.push_back(nDupes);
	row.push_back(smallest);
	Dupes.push_back(row);
      }
    }
  }
  return(Dupes);
}


void AnaDupe::processRunHeader( LCRunHeader* /*run*/) {
}

void AnaDupe::processEvent( LCEvent * evt ) {
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.


  LCCollection* oldTrkCol = evt->getCollection("SiTracks");
  LCCollection* newTrkCol = evt->getCollection("MySiTracks");
  
  std::vector<std::vector<int>> oldDupes = findDuplicates(oldTrkCol);
  std::vector<std::vector<int>> newDupes = findDuplicates(newTrkCol);
  int evtnum=evt->getEventNumber();

  //stream out the information for further analysis -- if you'd rather plot it yourself or do something with it
  //but also aggregate to the big vecs we'll use in the end function to plot with root if you want that too
  std::fstream dupeFile0;
  dupeFile0.open("oldDupeTrks.txt", std::ios::app);
  for(int i=0; i<oldDupes.size(); i++){
    dupeFile0<<oldDupes[i][0]<<' '<<oldDupes[i][1]<<'\n';
    _oldDupeTotal.push_back(oldDupes[i][0]);
    _oldHitsTotal.push_back(oldDupes[i][1]);
  }

  dupeFile0.close();
  std::fstream dupeFile1;
  dupeFile1.open("newDupeTrks.txt", std::ios::app);
  for(int i=0; i<newDupes.size(); i++){
    dupeFile1<<newDupes[i][0]<<' '<<newDupes[i][1]<<'\n';
    _newDupeTotal.push_back(newDupes[i][0]);
    _newHitsTotal.push_back(newDupes[i][1]);
  }
  dupeFile1.close();

  
}
void AnaDupe::check( LCEvent * /*evt*/ )
{ }

void AnaDupe::end()
{ 

  // now we plot (need to debug this further)

  //plot(_oldDupeTotal, _oldHitsTotal, "oldDupesperHit.root");
  //plot(_newDupeTotal, _newHitsTotal, "newDupesperHit.root");


}
