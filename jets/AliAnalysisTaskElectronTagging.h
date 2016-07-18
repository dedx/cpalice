/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
/* AliAnalysisTaskEx01.h
 *
 * Template task producing a P_t spectrum and pseudorapidity distribution.
 * Includes explanations of physics and primary track selections
 *
 * Based on tutorial example from offline pages
 * Edited by Arvinder Palaha
 */
#ifndef ALIANALYSISTASKELECTRONTAGGING_H
#define ALIANALYSISTASKELECTRONTAGGING_H

class TH1F;
class TList;
class AliESDtrackCuts;

#ifndef ALIANALYSISTASKSE_H
#include <AliAnalysisTaskSE.h>
#endif

class AliAnalysisTaskElectronTagging : public AliAnalysisTaskSE {
 public:
    AliAnalysisTaskElectronTagging();
    AliAnalysisTaskElectronTagging(const char *name);
    virtual ~AliAnalysisTaskElectronTagging();
    
    virtual void     UserCreateOutputObjects();
    virtual void     UserExec(Option_t *option);
    virtual void     Terminate(Option_t *);
    
 private:
    TList           *fOutput;        //! 
    AliESDtrackCuts *fTrackCuts;     //
    TH2F            *fHistDeDx;        //! 
    TH1F            *fHistRej;
    TH2F            *fHistTOF;
    TH2F            *fHistTRD;
    // NEW HISTO to be declared here
    
    AliAnalysisTaskElectronTagging(const AliAnalysisTaskElectronTagging&); // not implemented
    AliAnalysisTaskElectronTagging& operator=(const AliAnalysisTaskElectronTagging&); // not implemented
    
    ClassDef(AliAnalysisTaskElectronTagging, 1); // example of analysis
};

#endif

