#pragma once

#include <TH1F.h>
#include <TH2F.h>
#include <TGraph.h>

#include <marlin/Processor.h>


class EfficiencyProcessor : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new EfficiencyProcessor ; }

  EfficiencyProcessor(const EfficiencyProcessor &) = delete ;
  EfficiencyProcessor& operator =(const EfficiencyProcessor &) = delete ;
  EfficiencyProcessor() ;

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
 
  virtual const double eta( const double theta);

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;  

private:
  //! Input Collection
  std::string _inputCollectionNameR {};
  std::string _inputCollectionNameMCP {};

  //! Output Collection
  //std::string _outputCollectionName {};

  //! PDG selection
  float _minPt {};
  float _minTheta {};

  //! Output histogram
  // declare all histograms
  TH1 *_h_pass_pT = nullptr;
  TH1 *_h_all_pT = nullptr;
  TH1 *_h_pass_theta = nullptr;
  TH1 *_h_all_theta = nullptr;
  TH1 *_h_pass_phi = nullptr;
  TH1 *_h_all_phi = nullptr;
  TH1 *_h_pass_z = nullptr;
  TH1 *_h_all_z = nullptr;

  TH1 *_h_notrk_pt = nullptr;
  TH1 *_h_notrk_theta = nullptr;
  TH1 *_h_notrk_phi = nullptr;
  TH1 *_h_notrk_pdg = nullptr;

  TH1 *_h_matched_d0 = nullptr;
  TH1 *_h_unmatched_d0 = nullptr;
  TH1 *_h_matched_z0 = nullptr;
  TH1 *_h_unmatched_z0 = nullptr;

  TH1 *_h_trkclusterdist = nullptr;
  TH1 *_h_paralleldist = nullptr;
  TH1 *_h_clEnergyRes = nullptr;
  TH1 *_h_noCL_taupt = nullptr;
  TH1 *_h_noCL_theta = nullptr;


  TH2 *_h_pt_vs_theta = nullptr;

};
