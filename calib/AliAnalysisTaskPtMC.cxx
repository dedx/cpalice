#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTaskPtMC.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPtMC)

//________________________________________________________________________
AliAnalysisTaskPtMC::AliAnalysisTaskPtMC(const char *name) 
  : AliAnalysisTaskSE(name), fOutputList(0), fHistPt(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskPtMC::UserCreateOutputObjects() 
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 15, 0.1, 3.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  fOutputList->Add(fHistPt);
  PostData(1,fOutputList); // important for merging
}

//________________________________________________________________________
void AliAnalysisTaskPtMC::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  // Process MC truth

  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent) {
     Printf("ERROR: Could not retrieve MC event");
     return;
  }

  Printf("MC particles: %d", mcEvent->GetNumberOfTracks());

  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) {
    AliVParticle* track = mcEvent->GetTrack(iTracks);
    if (!track) {
      Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
      continue;
    }
      
    fHistPt->Fill(track->Pt());
  } //track loop 

  // Post output data.
  PostData(1, fOutputList);
}      

//________________________________________________________________________
void AliAnalysisTaskPtMC::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: Output list not available");
    return;
  }

  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {
    Printf("ERROR: fHistPt not available");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysisTaskPtMC","Pt MC",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
}
