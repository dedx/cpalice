#ifndef AliAnalysisTaskPtMC_cxx
#define AliAnalysisTaskPtMC_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPtMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPtMC() : AliAnalysisTaskSE(), fOutputList(0), fHistPt(0) {}
  AliAnalysisTaskPtMC(const char *name);
  virtual ~AliAnalysisTaskPtMC() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  TList       *fOutputList; //! Output list
  TH1F        *fHistPt; //!Pt spectrum
   
  AliAnalysisTaskPtMC(const AliAnalysisTaskPtMC&); // not implemented
  AliAnalysisTaskPtMC& operator=(const AliAnalysisTaskPtMC&); // not implemented

  ClassDef(AliAnalysisTaskPtMC, 1); // example of analysis
};

#endif
