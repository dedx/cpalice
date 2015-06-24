#ifndef AliAnalysisTaskPJ_cxx
#define AliAnalysisTaskPJ_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TH2F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskPJ : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPJ(); //default constructor
  AliAnalysisTaskPJ(const char *name);
  virtual ~AliAnalysisTaskPJ() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual Double_t CalcMass(Double_t MassHyp,Double_t time, Double_t eta);
  virtual void   Terminate(Option_t *);
  
 private:
  AliVEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistPt; //! Track Pt spectrum
  TH2F        *fHistTOF; //Track TOF Spectrum
  TH1F        *fHistNumTOFTOT;
  TH1F        *fHistNumTOFTDC;
  TH1F        *fHistUnmatchedClusTOF;
  TH1F        *fHistUnmatchedClusEM;
  TH1F        *fHistUnmatchedTOFEnergy;
  TH1F        *fHistUnmatchedEMEnergy;
  TH2F        *fHistTOFEMEnergyMatch;
  TH1F        *fHistTotalClusTOF;
  TH1F        *fHistTotalClusALLTOF;
  TH1F        *fHistTotalClusEM;
  TH1F        *fHistUnmatchedTOF;
  TH2F        *fHistDeltaE;
  TH2F        *fHistDeltaADC;
  TH2F        *fHistEtaPhiCC[100];
  TH2F        *fHistEtaPhiTC[100]; // my plot
  TH1F        *fHistNumCC;
  TH1F        *fHistNumTC;
  Int_t       fEvtNum; // keep track of event number
  Int_t       fHistNum; // store every Nth event's info in histogram

  AliAnalysisTaskPJ(const AliAnalysisTaskPJ&); // not implemented
  AliAnalysisTaskPJ& operator=(const AliAnalysisTaskPJ&); // not implemented
  
  ClassDef(AliAnalysisTaskPJ, 1); // example of analysis
};

#endif
