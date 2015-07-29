
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDTZERO.h"
#include "AliESDInputHandlerRP.h"
#include "AliTOFcluster.h"
#include "AliCluster.h"
#include "AliTOFGeometry.h"

#include "AliAnalysisTaskPJ_B.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPJ_B)

//________________________________________________________________________
AliAnalysisTaskPJ_B::AliAnalysisTaskPJ_B(const char *name) 
: AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fHistPt(0), fHistTC(0)


{
  // Constructor
  

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskPJ_B::AliAnalysisTaskPJ_B() 
  : AliAnalysisTaskSE(), fESD(0), fOutputList(0), fHistPt(0), fHistTC(0)


{
  // Default Constructor
  
  
}

//________________________________________________________________________
void AliAnalysisTaskPJ_B::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging
  
  fHistPt = new TH1F("fHistPt", "EMCAL Cluster-Matching distribution", 20, 0, 1);
  fHistPt->GetXaxis()->SetTitle("Percentage of Clusters Matched per Event");
  fHistPt->GetYaxis()->SetTitle("Counts");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  fHistTC = new TH1F("fHistTC", "TOF Cluster-Matching distribution", 200, 0, 100);
  fHistTC->GetXaxis()->SetTitle("Percentage of Clusters Matched per Event");
  fHistTC->GetYaxis()->SetTitle("Counts");
  fHistTC->SetMarkerStyle(kFullCircle);
   
  fHistTOFMatch = new TH1F("fHistTOFMatch", "# of TOF clusters matched to tracks per event", 10, 0, 10);
  fHistTOFMatch->GetXaxis()->SetTitle("# of TOF clusters");
  fHistTOFMatch->GetYaxis()->SetTitle("Counts");

  

  fOutputList->Add(fHistPt);
  fOutputList->Add(fHistTC);
  fOutputList->Add(fHistTOFMatch);
  PostData(1,fOutputList); // important for merging
}

//________________________________________________________________________
void AliAnalysisTaskPJ_B::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fInputEvent) {
    printf("ERROR: Input event not available\n");
    return;
  }
  // Get rec points
  AliESDInputHandlerRP *handler = 
    dynamic_cast<AliESDInputHandlerRP*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) { printf("No RP handler\n"); return; }
  AliESDEvent *esd  = handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  
  //get reconstructed vertex position
  Double_t vertex_position[3];
  esd->GetVertex()->GetXYZ(vertex_position);
  
  TTree* tofClusterTree = handler->GetTreeR("TOF");
  if (!tofClusterTree) {printf("<WARN>No TOF clusters!\n");return;}
  
  TBranch* tofClusterBranch = tofClusterTree->GetBranch("TOF");
  if (!tofClusterBranch) {printf("Can't get the branch with the TOF digits !\n"); return;}
  
  TClonesArray dummy("AliTOFcluster",1000), * tofClusters = &dummy;
  tofClusterBranch->SetAddress(&tofClusters);
  
  tofClusterTree->GetEvent(0);//  this may be the error
  
  Double_t nClusters = tofClusters->GetEntriesFast();
  printf("<INFO>NTOFclusters=%d\n",nClusters);
  cout<<"TOF Clusters: " <<nClusters<<endl;
  
  
  printf("<INFO>Ntracks=%d\n",fInputEvent->GetNumberOfTracks());
  TRefArray* caloClusters = new TRefArray();
  fInputEvent->GetEMCALClusters(caloClusters);
  Double_t nclus = caloClusters->GetEntries(); 
  Double_t EMcount=0;
  Double_t TOFcount=0;

  TClonesArray* clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));

  for(Int_t ic=0;ic<clusArr->GetEntries();ic++){
    AliESDCaloCluster* clus = (AliESDCaloCluster*)clusArr->At(ic);
    if(!clus){break;}
    fHistTC->Fill(clus->GetTrackMatchedIndex());
  }
  /*
  // Track loop to fill a pT spectrum
    if(nClusters>0 && fInputEvent->GetNumberOfTracks()>0 && nclus>0){
  for (Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++) {
    AliESDtrack* track = (AliESDtrack*)fInputEvent->GetTrack(iTracks);
    if (!track) { printf("ERROR: Could not receive track %d\n", iTracks);continue;}
    //printf("Track %d has been matched to EMCAL %d", track->GetID(), track->GetEMCALcluster());
    //AliTOFcluster* tof = (AliTOFcluster*); 
     fHistTOFMatch->Fill(track->GetTOFcluster());
      for(Int_t iclus=0;iclus<nclus;iclus++){
          AliESDCaloCluster *clus = (AliESDCaloCluster*)caloClusters->At(iclus);
          if(clus->GetID()==track->GetEMCALcluster()){cout<<"NICE!"<<EMcount/nclus; EMcount++;}
      }
      for(Int_t iclus=0;iclus<nClusters;iclus++){
          AliTOFcluster *tclus = (AliTOFcluster*)tofClusters->UncheckedAt(iclus);
          cout<<tclus->GetESDID();
          if(track->GetTOFcluster()!=-1 && tclus->GetIndex()!=-1){
          if(tclus->GetIndex()==track->GetTOFcluster()){cout<<TOFcount/nClusters<<"\n"; TOFcount++;}
          }
      }
    //EmCal = new AliESDCaloCluster();
    //EmCal->SetID(track->GetEMCALcluster());
    //cout<<"EMCCal Energy "<<EmCal->E();
    //Int_t tofIndex = track->GetTOFcluster();
    // printf("\t<I>Track %d has tofIndex %d\n",iTracks,tofIndex);
    
  }
    
    fHistPt->Fill(EMcount/nclus);
    fHistTC->Fill(TOFcount/nClusters);
    }
  
  for (Int_t icl = 0; icl < nclus; icl++) {    
    AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
   
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 vpos(pos[0],pos[1],pos[2]);
    
    Double_t cphi = vpos.Phi();
    if(cphi < 0) cphi +=TMath::TwoPi();
    Double_t ceta = vpos.Eta();
    
    // Time of Flight (TOF)
    Double_t EmCaltof = clus->GetTOF();
    Double_t EmCalEnergy = clus->E();
    //Print basic cluster information
    cout << "Cluster: " << icl+1 << "/" << nclus << "Phi: " << 
       cphi*TMath::RadToDeg() << "; Eta: " << ceta << endl;

  }
  */

  //clean up to avoid mem leaks
  delete tofClusterTree;
  PostData(1, fOutputList);
  
}
//________________________________________________________________________
void AliAnalysisTaskPJ_B::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {printf("ERROR: Output list not available\n");return;}
  
  fHistTC = dynamic_cast<TH1F*> (fOutputList->At(1));
  if (!fHistTC) {printf("ERROR: fHistPt not available\n");return;}  
}
