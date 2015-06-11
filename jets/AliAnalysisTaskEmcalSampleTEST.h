#ifndef ALIANALYSISTASKEMCALSAMPLETEST_H
#define ALIANALYSISTASKEMCALSAMPLETEST_H

// $Id$

class TH1;
class TH2;
class TH3;
class AliParticleContainer;
class AliClusterContainer;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalSampleTEST : public AliAnalysisTaskEmcal {
 public:

  AliAnalysisTaskEmcalSampleTEST();
  AliAnalysisTaskEmcalSampleTEST(const char *name);
  virtual ~AliAnalysisTaskEmcalSampleTEST();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms()   ;
  Bool_t                      Run()              ;
  void                        CheckClusTrackMatching();

  // General histograms
  TH1                       **fHistTracksPt;            //!Track pt spectrum
  TH1                       **fHistClustersPt;          //!Cluster pt spectrum
  TH3                        *fHistPtDEtaDPhiTrackClus; //!track pt, delta eta, delta phi to matched cluster
  TH3                        *fHistPtDEtaDPhiClusTrack; //!cluster pt, delta eta, delta phi to matched track

  AliParticleContainer       *fTracksCont;                 //!Tracks
  AliClusterContainer        *fCaloClustersCont;           //!Clusters  

 private:
  AliAnalysisTaskEmcalSampleTEST(const AliAnalysisTaskEmcalSampleTEST&);            // not implemented
  AliAnalysisTaskEmcalSampleTEST &operator=(const AliAnalysisTaskEmcalSampleTEST&); // not implemented

  ClassDef(AliAnalysisTaskEmcalSampleTEST, 1) // emcal sample analysis task
};
#endif
