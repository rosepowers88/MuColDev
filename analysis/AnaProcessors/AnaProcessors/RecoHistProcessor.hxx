#pragma once

#include <TH1F.h>

#include <marlin/Processor.h>

//! A template processor that should be modified to do what you want.
/**
 * Implements a loop over the MCParticle collection and creates an output
 * colleciton passing certain criteria. This provides rovides examples of:
 *  - using parameters to configure a processor
 *  - opening collections
 *  - outputing histograms
 */
class RecoHistProcessor : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new RecoHistProcessor ; }

  RecoHistProcessor(const RecoHistProcessor &) = delete ;
  RecoHistProcessor& operator =(const RecoHistProcessor &) = delete ;
  RecoHistProcessor() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void book_pfo_histograms();
  virtual void book_cluster_histograms();
  virtual void book_track_histograms();

  virtual void fill_pfo_histograms(LCCollection* inputCol);
  virtual void fill_cluster_histograms(LCCollection* inputCol);
  virtual void fill_track_histograms(LCCollection* inputCol);
  
  virtual void init() ;

  /** Called for every run.
   */

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

  

  //! Output Collection
  //std::string _outputCollectionName {};

  //! Minimum pT for particle filter
  float _minPt = 1;

  //! Output histogram
  // declare all histograms

  //pfo histograms
  TH1 *_h_ptpf = nullptr;
  TH1 *_h_ppf = nullptr;
  TH1 *_h_Npf = nullptr;
  TH1 *_h_pdgpf = nullptr;
  TH1 *_h_Epf = nullptr;

  //cluster histograms
  TH1 *_h_Ecl = nullptr;
  TH1 *_h_Xcl = nullptr;
  TH1 *_h_Ycl = nullptr;
  TH1 *_h_Zcl = nullptr;
  TH1 *_h_iPhicl = nullptr;
  TH1 *_h_iThetacl = nullptr;
  TH1 *_h_PIDcl = nullptr;
  TH1 *_h_Ncl   = nullptr;

  //track histograms
  TH1 *_h_D0trk = nullptr;
  TH1 *_h_phitrk = nullptr;
  TH1 *_h_omegatrk = nullptr;
  TH1 *_h_Z0trk = nullptr;
  TH1 *_h_lamtrk = nullptr;
  TH1 *_h_dEdxtrk = nullptr;
  TH1 *_h_Ntrk     = nullptr;

  

};
