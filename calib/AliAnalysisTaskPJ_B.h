#ifndef AliAnalysisTaskPJ_B_cxx
#define AliAnalysisTaskPJ_B_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPJ_B : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPJ_B(); //default constructor
  AliAnalysisTaskPJ_B(const char *name, TChain *pass);
  virtual ~AliAnalysisTaskPJ_B() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliVEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistPt; //! Track Pt spectrum
  TH1F        *fHistTC; //!
  TH1F        *fHistTOFMatch; //!
  TChain      *pass1tree;
  AliAnalysisTaskPJ_B(const AliAnalysisTaskPJ_B&); // not implemented
  AliAnalysisTaskPJ_B& operator=(const AliAnalysisTaskPJ_B&); // not implemented
  
  ClassDef(AliAnalysisTaskPJ_B, 1); // example of analysis
};

#endif
