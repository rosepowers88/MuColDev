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
  std::string _inputCollectionName {};

  //! Output Collection
  std::string _outputCollectionName {};

  //! Minimum pT for particle filter
  float _minPt = 1;

  //! Output histogram
  // declare all histograms
  TH1 *_h_pt = nullptr;
  TH1 *_h_p = nullptr;
  TH1 *_h_N = nullptr;
  TH1 *_h_pdg = nullptr;
  TH1 *_h_E = nullptr;

  TH1 *_h_phi = nullptr;
  TH1 *_h_tanLam = nullptr;
  //pis only
  TH1 *_h_p_pi= nullptr;
  TH1 *_h_pt_pi= nullptr;
  TH1 *_h_N_pi = nullptr;
  TH1 *_h_pdg_pi = nullptr;
  TH1 *_h_theta = nullptr;
  TH1 *_h_E_pi  = nullptr;
  TH1 *_h_x_pi = nullptr;
  TH1 *_h_y_pi = nullptr;
  TH1 *_h_z_pi = nullptr;
  TH1 *_h_r_pi = nullptr;
  

};
