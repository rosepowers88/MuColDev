#pragma once

#include <TH1F.h>
#include <TGraph.h>

#include <marlin/Processor.h>


class RecoPerformance : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new RecoPerformance ; }

  RecoPerformance(const RecoPerformance &) = delete ;
  RecoPerformance& operator =(const RecoPerformance &) = delete ;
  RecoPerformance() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */

  virtual std::vector<std::vector<int>> getSimTrackerHitIDs(LCCollection* trackhitCol);
  virtual double anaPiEff(const LCObject *inputTrk, LCCollection* trkHitCol);
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */

  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;  

private:
  //! Input Collection
  std::string _inputCollectionNameR {};
  std::string _inputCollectionNameC {};
  std::string _inputCollectionNameT {};
  std::string _inputCollectionNameMCP {};
  std::string _inputCollectionSimHits {};

  //! Output Collection
  //std::string _outputCollectionName {};

  //! PDG selection
  int _PDG {};

  float _minPt {};
  float _minTheta {};

  //! Output histogram
  // declare all histograms
  TH1 *_h_effPFO = nullptr;
  TH1 *_h_effCl = nullptr;
  TH1 *_h_effTrk = nullptr;

  TH1 *_h_effPFO_cp = nullptr;
  TH1 *_h_effPFO_np = nullptr;

  TH1 *_h_effPFO_n = nullptr;
  TH1 *_h_effPFO_gam = nullptr;



};
