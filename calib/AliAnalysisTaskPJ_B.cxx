
#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TColor.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESD.h"
#include "AliESDInputHandler.h"
#include "AliESDTZERO.h"
#include "AliESDInputHandlerRP.h"
#include "AliTOFcluster.h"
#include "AliCluster.h"
#include "AliTOFGeometry.h"
#include "AliTrackerBase.h"
#include "AliCDBManager.h"
#include "AliGeomManager.h"

#include "AliAnalysisTaskPJ_B.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPJ_B)

//________________________________________________________________________
AliAnalysisTaskPJ_B::AliAnalysisTaskPJ_B(const char* name, TChain* pass1) 
: AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fHistPt(0), fHistTC(0)


{
  // Constructor
  

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
  pass1tree=pass1;
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
  
  fHistPt = new TH1F("fHistPt", "EMCAL Cluster-Matching distribution", 20, 0, 2);
  fHistPt->GetXaxis()->SetTitle("Percentage of Clusters Matched per Event");
  fHistPt->GetYaxis()->SetTitle("Counts");
  fHistPt->SetMarkerStyle(kFullCircle);
  
  fHistTC = new TH1F("fHistTC", "TOF Cluster-Matching distribution", 10000, -500, 500);
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
  // Get Event
  AliESDInputHandlerRP *handler = 
    dynamic_cast<AliESDInputHandlerRP*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) { printf("No RP handler\n"); return; }
  AliESDEvent *esd  = handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  
  //get reconstructed vertex position
  Double_t vertex_position[3];
  esd->GetVertex()->GetXYZ(vertex_position);
    
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  cdbManager->SetRun(esd->GetRunNumber());  
  AliGeomManager::LoadGeometry();
    
  /*TTree* tree = (TTree*)pass1tree->Get("esdTree");
    
  AliESDEvent *RecoESD = new AliESDEvent();
  RecoESD->ReadFromTree(tree);
  RecoESD->Print();
  cout<<RecoESD->GetNumberOfTracks()<<"\n";
    for(Int_t itrk=0;itrk<RecoESD->GetNumberOfTracks();itrk++){
  fHistTC->Fill(RecoESD->GetTrack(itrk)->GetOuterParam()->GetZ());
    }*/
    
  //Get all the TOF Clusters from the TOFRecPoints File
    TTree* tofClusterTree = handler->GetTreeR("TOF");
  if (!tofClusterTree) {printf("<WARN>No TOF clusters!\n");return;}
  
  TBranch* tofClusterBranch = tofClusterTree->GetBranch("TOF");
  if (!tofClusterBranch) {printf("Can't get the branch with the TOF digits !\n"); return;}
  
  TClonesArray dummy("AliTOFcluster",1000), * tofClusters = &dummy;
  tofClusterBranch->SetAddress(&tofClusters);
  
  tofClusterTree->GetEvent(0);//  this may be the error
  
  Int_t nClusters = tofClusters->GetEntriesFast();
  printf("<INFO>NTOFclusters=%d\n",nClusters);
  cout<<"TOF Clusters: " <<nClusters<<endl;
  
  AliESDEvent *RecoESD = new AliESDEvent();
  RecoESD->ReadFromTree((TTree*)pass1tree);

  /*pass1tree->SetMakeClass(1);
  pass1tree->SetBranchStatus("*", 0);
  pass1tree->SetBranchStatus("AliESDHeader.fTimeStamp", 1);
  pass1tree->SetBranchStatus("Tracks*",1);
    
  if(!pass1tree){cout<<"No Chain";}
  
  UInt_t recoESD = 0;
  pass1tree->SetBranchAddress("AliESDHeader.fTimeStamp",&recoESD);

  //TClonesArray tdummy("AliESDtrack");
  TClonesArray* recoTrack = new TClonesArray("AliESDtrack");
  pass1tree->SetBranchAddress("Tracks", recoTrack);*/
  Bool_t foundOne=kFALSE;
  Int_t mEv = -999;
  Int_t numev = pass1tree->GetEntries();
    for(Int_t iev=0;iev<numev;iev++){
        pass1tree->GetEntry(iev);
        if(RecoESD->GetTimeStamp()==esd->GetTimeStamp()){
            printf("Found one @ %d", iev);
            fHistPt->Fill(1);
            foundOne=kTRUE;
            mEv = iev;
            break;
        }
    }
    
    
    if(foundOne){
    pass1tree->GetEntry(mEv);
    cout<<RecoESD->GetNumberOfTracks();
    for(Int_t icl = 0; icl<RecoESD->GetNumberOfTracks();icl++){
        AliESDtrack* trk = (AliESDtrack*)RecoESD->GetTrack(icl);
        cout<<icl;
        if(!trk){break;}
        AliExternalTrackParam* trkpar = (AliExternalTrackParam*)trk->GetOuterParam();
        if(!trkpar){break;}
        if(!AliTrackerBase::PropagateTrackTo(trkpar, 370, trk->M(), 1, kTRUE, -1, trkpar->GetSign(), kFALSE, kTRUE)){break;}
        //if(!trkpar->PropagateTo(370, RecoESD->GetMagneticField())){break;}
        Double_t trkVec[3];
        Float_t tofVec[3];
        trkpar->GetXYZ(trkVec);
        
        for(Int_t itof = 0; itof<tofClusters->GetEntriesFast(); itof++){
            AliCluster* tofClus = (AliCluster*)tofClusters->At(itof);
            tofClus->GetGlobalXYZ(tofVec);
            //if(sqrt(pow(trkVec[0]-tofVec[0],2)+pow(trkVec[1]-tofVec[1],2)+pow(trkVec[2]-tofVec[2],2))<10){
            fHistTC->Fill(trkVec[0]-tofVec[0]);
    }
    }
    }
    
  
  /*
  
  TRefArray* caloClusters = new TRefArray();
  esd->GetEMCALClusters(caloClusters);
  //Double32_t nclus = esd->GetT0amplitude(); 
  printf("<INFO>Nclusters=%f\n",vertex_position[0]);
  Double_t EMcount=0;
  Double_t TOFcount=0;
  TClonesArray* clusArr = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("EmcCaloClusters"));

  for(Int_t ic=0;ic<clusArr->GetEntries();ic++){
    AliESDCaloCluster* clus = (AliESDCaloCluster*)clusArr->At(ic);
    if(!clus){break;}
    fHistTC->Fill(clus->GetTrackMatchedIndex());
  }
  
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
    //delete tofClusterTree;
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
