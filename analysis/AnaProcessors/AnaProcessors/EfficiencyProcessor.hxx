#pragma once

#include <TH1F.h>
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



};
