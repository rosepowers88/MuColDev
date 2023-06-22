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
class BranchingProcessor : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new BranchingProcessor ; }

  BranchingProcessor(const BranchingProcessor &) = delete ;
  BranchingProcessor& operator =(const BranchingProcessor &) = delete ;
  BranchingProcessor() ;

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

  //! PDG selection
  int _PDG {};

  //! Output histogram
  // declare all histograms
  TH1 *_h_pdg = nullptr;
  TH1 *_h_MID = nullptr;
  TH1 *_h_q = nullptr;

};
