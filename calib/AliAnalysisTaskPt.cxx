#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDInputHandlerRP.h"
#include "AliTOFcluster.h"
#include "AliCluster.h"

#include "AliAnalysisTaskPt.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPt)

//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt(const char *name) 
: AliAnalysisTaskSE(name), fOutputList(0), fHistPt(0), fHistNumCC(0),fHistNumTC(0),fEvtNum(0),fHistNum(0)
{
  // Constructor
  for (int i = 0; i < 100; i++) {
    fHistEtaPhiCC[i] = 0;
    fHistEtaPhiTC[i] = 0;
    fHistTOFCC[i] = 0;
    fHistTOFTC[i] = 0;
  }

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskPt::AliAnalysisTaskPt() 
  : AliAnalysisTaskSE(), fESD(0), fOutputList(0), fHistPt(0), fHistNumCC(0),fHistNumTC(0),fEvtNum(0),fHistNum(0)
{
  // Default Constructor
  for (int i = 0; i < 100; i++) {
    fHistEtaPhiCC[i] = 0;
    fHistEtaPhiTC[i] = 0;
    fHistTOFCC[i] = 0;
    fHistTOFTC[i] = 0;
    
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskPt::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 100, 0.1, 10.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);  

  for (int i = 0; i < 100; i++) {
    char name[100];
    sprintf(name,"fHistEtaPhiCC%d",i);
    fHistEtaPhiCC[i] = new TH2F(name,name,100,1.35,3.15,100,-1.,1);
    fHistEtaPhiCC[i]->GetXaxis()->SetTitle("#phi");
    fHistEtaPhiCC[i]->GetYaxis()->SetTitle("#eta");

    sprintf(name,"fHistEtaPhiTC%d",i);
    fHistEtaPhiTC[i] = new TH2F(name,name,100,1.35,3.15,100,-1.,1);
    fHistEtaPhiTC[i]->GetXaxis()->SetTitle("#phi");
    fHistEtaPhiTC[i]->GetYaxis()->SetTitle("#eta");

    sprintf(name,"fHistTOFCC%d",i);
    fHistTOFCC[i] = new TH1F(name,name,100,0.,100.);
    fHistTOFCC[i]->GetXaxis()->SetTitle("Calocluster time (ns)");

    sprintf(name,"fHistTOFTC%d",i);
    fHistTOFTC[i] = new TH1F(name,name,100,0.,100.);
    fHistTOFTC[i]->GetXaxis()->SetTitle("TOF cluster time (ns)");
  }

  fHistNumCC = new TH1F("fHistNumCC", "# of Calo clusters per event", 100,0,100);
  fHistNumCC->GetXaxis()->SetTitle("# Calo Clusters");
  fHistNumCC->GetYaxis()->SetTitle("Frequency");

  fHistNumTC = new TH1F("fHistNumTC", "# of TOF clusters per event", 100,0,100);
  fHistNumTC->GetXaxis()->SetTitle("# TOF Clusters");
  fHistNumTC->GetYaxis()->SetTitle("Frequency");

  fOutputList->Add(fHistPt);
  for (int i = 0; i < 100; i++) {
    fOutputList->Add(fHistEtaPhiCC[i]);
    fOutputList->Add(fHistEtaPhiTC[i]);
    fOutputList->Add(fHistTOFCC[i]);
    fOutputList->Add(fHistTOFTC[i]);
  }
  fOutputList->Add(fHistNumCC);
  fOutputList->Add(fHistNumTC);

  PostData(1,fOutputList); // important for merging
}

//________________________________________________________________________
void AliAnalysisTaskPt::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event

  if (!fInputEvent) {
    printf("ERROR: Input event not available\n");
    return;
  }
  
  // Get rec points
  AliESDInputHandlerRP *handler = dynamic_cast<AliESDInputHandlerRP*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) { printf("No RP handler\n"); return; }
  AliESDEvent *esd  = handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }


  //get reconstructed vertex position
  Double_t vertex_position[3];
  esd->GetVertex()->GetXYZ(vertex_position);

  TTree* tofClusterTree = handler->GetTreeR("TOF");
  if (!tofClusterTree) {
    printf("<WARN>No TOF clusters!\n");
    return;
  }

  TBranch* tofClusterBranch = tofClusterTree->GetBranch("TOF");
  if (!tofClusterBranch) {
    printf("Can't get the branch with the TOF digits !\n");
    return;
  }

  TClonesArray dummy("AliTOFcluster",10000), * tofClusters = &dummy;
  tofClusterBranch->SetAddress(&tofClusters);

  tofClusterTree->GetEvent(0); // this may be the error

  Int_t nClusters = tofClusters->GetEntriesFast();
  printf("<INFO>NTOFclusters=%d\n",nClusters);

  int countinEMCal = 0;
  // loop over clusters
  while (nClusters--) {
    AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(nClusters);
    if(cluster) {

      if (cluster->GetPhi() > 1.35 && cluster->GetPhi() < 3.15) {
	countinEMCal++;

	//TOF x,y,z
        Float_t x = cluster->GetR()*TMath::Cos(cluster->GetPhi());
	Float_t y = cluster->GetR()*TMath::Sin(cluster->GetPhi());
	Float_t z = cluster->GetZ();

	//JLK NEED TO FIGURE OUT HOW TO CONVERT TDC to TOF!!!!!
	Float_t tof = cluster->GetTDC(); //NOT time of flight in ns
	//JLK

	printf("CLUSTER=%d, r=%.2f,phi=%.2f,(x,y,z)=(%.2f,%.2f,%.2f)\n",cluster->GetIndex(),cluster->GetR(),cluster->GetPhi(),x,y,z);

	TVector3 vpos(x,y,z);
	Double_t cphi = vpos.Phi();
	if(cphi < 0) cphi +=TMath::TwoPi();
	Double_t ceta = vpos.Eta();
	printf("cphi,ceta = %.3f,%.3f\n",cphi,ceta);
	if (fEvtNum %2000 == 0 && fHistNum < 100) {
	  fHistEtaPhiTC[fHistNum]->Fill(cphi,ceta);
	  fHistTOFTC[fHistNum]->Fill(tof);
	}

      }
    }
  }

  fHistNumTC->Fill(countinEMCal);


  printf("<INFO>Ntracks=%d\n",fInputEvent->GetNumberOfTracks());

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++) {
    //AliVParticle* track = fInputEvent->GetTrack(iTracks);
    AliESDtrack* track = (AliESDtrack*)fInputEvent->GetTrack(iTracks);
    if (!track) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    
    Int_t tofIndex = track->GetTOFcluster();
    printf("\t<I>Track %d has tofIndex %d\n",iTracks,tofIndex);
    fHistPt->Fill(track->Pt());
  } //track loop


  //Init geometry and array that will contain the clusters
  TRefArray* caloClusters = new TRefArray();
  //Pass the geometry transformation matrix from ESDs to geometry
  
  // cell array
  //  AliVCaloCells &cells= *(fInputEvent->GetEMCALCells());
  
  //Get CaloClusters Array
  // GetEMCALClusters will provide only the EMCAL clusters where GetNCaloClusters
  // returns the total number EMCAL + PHOS
  // Bool_t IsEMCAL();
  fInputEvent->GetEMCALClusters(caloClusters);
  
  
  Int_t nclus = caloClusters->GetEntries();
    
  //  Int_t nCells = clus->GetNCells();
  fHistNumCC->Fill(nclus);
  
  for (Int_t icl = 0; icl < nclus; icl++) {
    
    AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
    //Float_t energy = clus->E();
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 vpos(pos[0],pos[1],pos[2]);
    
    //JLK Want cluster time in ns - not sure of the units here
    Double_t ctime = clus->GetTOF();
    //JLK

    Double_t cphi = vpos.Phi();
    if(cphi < 0) cphi +=TMath::TwoPi();
    Double_t ceta = vpos.Eta();
    //    cout<<"Calo Cluster Phi Value is:  "<<cphi<<endl;
    //cout<<"Calo Cluster Eta Value is: "<<ceta<<endl;

    if (fEvtNum%2000 == 0 && fHistNum < 100) {
      fHistEtaPhiCC[fHistNum]->Fill(cphi, ceta);
      fHistTOFCC[fHistNum]->Fill(ctime);
    }
    //ADD a nested loop over TOF clusters to compute Delta-R and
    //Delta-t between TOF and EMCAL





    //Print basic cluster information
    //    cout << "Cluster: " << icl+1 << "/" << nclus << "; Phi: " << cphi*TMath::RadToDeg() << "; Eta: " << ceta << "; NCells: " << nCells << endl;
  } //calocluster loop 
  
  
  //clean up to avoid mem leaks
  delete tofClusterTree;
  //Keep every 2000th event's info in histogram
  if (fEvtNum%2000 == 0) fHistNum++;
  //increment event counter
  fEvtNum++;

  PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskPt::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {
    printf("ERROR: fHistPt not available\n");
    return;
  }

  TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");
}
