/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

/* AliAnalysisTaskElectronTagging.cxx
 *
 * Macro for Electron PID for heavy flavour electrons
 * -Patrick Steffanic
 */
#include "AliAnalysisTaskElectronTagging.h"

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TList.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"
#include "AliStack.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDEvent.h"
#include "AliPID.h"

ClassImp(AliAnalysisTaskElectronTagging)

//________________________________________________________________________
AliAnalysisTaskElectronTagging::AliAnalysisTaskElectronTagging() // All data members should be initialised here
   :AliAnalysisTaskSE(),
    fOutput(0),
    fTrackCuts(0),
    fHistDeDx(0),
    fHistRej(0),
    fHistTOF(0),
    fHistTRD(0)// The last in the above list should not have a comma after it
{
    // Dummy constructor ALWAYS needed for I/O.
}

//________________________________________________________________________
AliAnalysisTaskElectronTagging::AliAnalysisTaskElectronTagging(const char *name) // All data members should be initialised here
   :AliAnalysisTaskSE(name),
    fOutput(0),
    fTrackCuts(0),
    fHistDeDx(0),
    fHistRej(0), 
    fHistTOF(0),
    fHistTRD(0)// The last in the above list should not have a comma after it
{
    // Constructor
    // Define input and output slots here (never in the dummy constructor)
    // Input slot #0 works with a TChain - it is connected to the default input container
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());                                            // for output list
}

//________________________________________________________________________
AliAnalysisTaskElectronTagging::~AliAnalysisTaskElectronTagging()
{
    // Destructor. Clean-up the output list, but not the histograms that are put inside
    // (the list is owner and will clean-up these histograms). Protect in PROOF case.
    if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
        delete fOutput;
    }
    delete fTrackCuts;
}

//________________________________________________________________________
void AliAnalysisTaskElectronTagging::UserCreateOutputObjects()
{
    // Create histograms
    // Called once (on the worker node)
        
    fOutput = new TList();
    OpenFile(1);
    fOutput->SetOwner();  // IMPORTANT!
    
    //Strong cuts for heavy flavour, maybe I will explain them in the future
    
    fTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(kTRUE);
    fTrackCuts->SetMinNClustersTPC(120);
    fTrackCuts->SetMinNClustersITS(4);
    fTrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kFirst);
    fTrackCuts->SetMaxChi2PerClusterTPC(2);
    fTrackCuts->SetMaxDCAToVertexXY(1);
    fTrackCuts->SetMaxDCAToVertexZ(2);
    fTrackCuts->SetEtaRange(-.6, .6);
    fTrackCuts->SetAcceptKinkDaughters(kFALSE);

    // Create histograms
    // TPC DE/Dx
    fHistDeDx = new TH2F("fHistDeDx", "De/Dx v. Pt", 100, 0, 10, 300, 0, 150);
    fHistDeDx->GetXaxis()->SetTitle("Pt");
    fHistDeDx->GetYaxis()->SetTitle("TPC De/Dx");
    // TOF time of flight
    fHistTOF = new TH2F("fHistTOF", "TOF signal v. Pt", 100, 0, 10, 15, -5, 10);
    fHistTOF->GetXaxis()->SetTitle("Pt");
    fHistTOF->GetYaxis()->SetTitle("TOF Signal");
    // TRD PID Probability
    fHistTRD = new TH2F("fHistTRD", "TRD PID Likelihood v. Pt", 100, 0, 10, 100, 0, 1);
    fHistTRD->GetXaxis()->SetTitle("Pt");
    fHistTRD->GetYaxis()->SetTitle("TRD electron likelihood");
    // Number of rejected events (change to a category histo)
    fHistRej = new TH1F("fHistRej", "Number of rejected tracks", 1,0,1);
    fHistRej->GetXaxis()->SetTitle("Rejected");
    fHistRej->GetYaxis()->SetTitle("Num");
        
    // NEW HISTO should be defined here, with a sensible name,
        
    fOutput->Add(fHistDeDx);
    fOutput->Add(fHistRej);
    fOutput->Add(fHistTOF);
    fOutput->Add(fHistTRD);
    // NEW HISTO added to fOutput here
    PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}

//________________________________________________________________________
void AliAnalysisTaskElectronTagging::UserExec(Option_t *) 
{
    // Main loop
    // Called for each event
        
        
    // Create pointer to reconstructed event
    AliVEvent *event = InputEvent();
    if (!event) { Printf("ERROR: Could not retrieve event"); return; }
        
    // create pointer to event
    AliESDEvent* esd = dynamic_cast<AliESDEvent*>(event);
    if (!esd) {
        AliError("Cannot get the ESD event");
        return;
    }  
//    AliESDHeader* esdheader = (AliESDHeader*)esd->GetHeader();
        
    // === Physics Selection Task ===
    // 
    // To perform a physics selection here, a bitwise operation is used against
    // the UInt_t mask which is extracted in the following way:
    //
    //  UInt_t mask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();   
    //
    // This can be tested to produce the following
    //
    //  Bool_t bMinBias = (mask == AliVEvent::kMB) ? 1 : 0; // check if minimum bias trigger class fired
    //  Bool_t bHighMult = (mask == AliVEvent::kHighMult) ? 1 : 0; // check if high multiplicity trigger class fired
    //
    // For more complicated trigger selections, one can directly test both
    // trigger classes and fired trigger inputs for a particular event, for e.g.
    //
    //  Bool_t bCSH1 = (esd->IsTriggerClassFired("CSH1-B-NOPF-ALLNOTRD")) ? 1 : 0;
    //  Bool_t b0SH1 = (esdheader->IsTriggerInputFired("0SH1")) ? 1 : 0;
    //
    // These booleans can then be used to fill different histograms for specific
    // conditions, or summed to make one cut for all events that fill histos.


        
    // Track loop for reconstructed event
    Int_t ntracks = esd->GetNumberOfTracks();
    for(Int_t i = 0; i < ntracks; i++) {
        AliESDtrack* esdtrack = esd->GetTrack(i); // pointer to reconstructed to track          
        if(!esdtrack) { 
            AliError(Form("ERROR: Could not retrieve esdtrack %d",i)); 
            continue; 
        }
                
        if(!fTrackCuts->AcceptTrack(esdtrack)){
            fHistRej->Fill(0);
            continue;
        }
                
        fHistDeDx->Fill(esdtrack->Pt(),esdtrack->GetTPCsignal());
        fHistTOF->Fill(esdtrack->Pt(), esdtrack->GetTOFsignal());
        fHistTRD->Fill(esdtrack->Pt(), esdtrack->GetTRDpid(AliPID::kElectron));
    }
    // NEW HISTO should be filled before this point, as PostData puts the
    // information for this iteration of the UserExec in the container
    PostData(1, fOutput);
}


//________________________________________________________________________
void AliAnalysisTaskElectronTagging::Terminate(Option_t *) 
{
    // Draw result to screen, or perform fitting, normalizations
    // Called once at the end of the query
    fOutput = dynamic_cast<TList*> (GetOutputData(1));
    if(!fOutput) { Printf("ERROR: could not retrieve TList fOutput"); return; }
        
    fHistDeDx = dynamic_cast<TH2F*> (fOutput->FindObject("fHistDeDx"));
    if (!fHistDeDx) { Printf("ERROR: could not retrieve fHistDeDx"); return;}
        
    //Get the physics selection histograms with the selection statistics
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliESDInputHandler *inputH = dynamic_cast<AliESDInputHandler*>(mgr->GetInputEventHandler());
    TH2F *histStat = (TH2F*)inputH->GetStatistics();
   
   
    // NEW HISTO should be retrieved from the TList container in the above way,
    // so it is available to draw on a canvas such as below

    TCanvas *c = new TCanvas("AliAnalysisTaskElectronTagging","P_{T} & #eta",10,10,1020,510);
    fHistDeDx->DrawCopy("E");
}
