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
#include "AliESDInputHandlerRP.h"
#include "AliTOFcluster.h"
#include "AliCluster.h"
#include "AliTOFGeometry.h"

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

  fHistTOF = new TH1F("fHistTOF", "TOF distribution", 100, 0, 1E-6);
  fHistTOF->GetXaxis()->SetTitle("TOF (Sec)");
  fHistTOF->GetYaxis()->SetTitle("Counts");

  fHistDeltaR = new TH1F("fHistDeltaR", "2D DeltaE-DeltaR distribution", 10000, 0, 5);
  fHistDeltaR->GetXaxis()->SetTitle("Delta R ()");
  fHistDeltaR->GetYaxis()->SetTitle("cts");

  fHistDeltaE = new TH2F("fHistDeltaE", "2D EMCal Energy-Delta R", 10000, 0, 5, 10000, .1, .5);
  fHistDeltaE->GetXaxis()->SetTitle("Delta R");
  fHistDeltaE->GetYaxis()->SetTitle("EMCal Energy");

  fHistDeltaT = new TH2F("fHistDeltaT", "2D DeltaT-DeltaR", 10000, 0, 5, 10000, 0, 500);
  fHistDeltaT->GetXaxis()->SetTitle("Delta R");
  fHistDeltaT->GetYaxis()->SetTitle("Delta T");

  fHistDeltaADC = new TH2F("fHistDeltaADC", "2D ADC-DeltaR", 10000, 0, 5, 10000, 100, 300);
  fHistDeltaADC->GetXaxis()->SetTitle("Delta R");
  fHistDeltaADC->GetYaxis()->SetTitle("TOF ADC");

  fHistEADC = new TH2F("fHistEADC", "Ratio of EMCAL Energy to ADC at small Delta R", 10000, 0, 500, 10000, 0, .05);
  fHistEADC->GetXaxis()->SetTitle("TOF Cluster");
  fHistEADC->GetYaxis()->SetTitle("E/ADC");

  for (int i = 0; i < 100; i++) {
    char name[100];sprintf(name,"fHistEtaPhiCC%d",i);
    fHistEtaPhiCC[i] = new TH2F(name,name,100,1.35,3.15,100,-1.,1);
    fHistEtaPhiCC[i]->GetXaxis()->SetTitle("#phi");
    fHistEtaPhiCC[i]->GetYaxis()->SetTitle("#eta"); }


  for (int i = 0; i < 100; i++) {
    char name[100]; sprintf(name,"fHistEtaPhiTC%d",i);
    fHistEtaPhiTC[i] = new TH2F(name,name,100,1.35,3.15,100,-1.,1);
    fHistEtaPhiTC[i]->GetXaxis()->SetTitle("#phi");
    fHistEtaPhiTC[i]->GetYaxis()->SetTitle("#eta");
  }

  fHistNumCC = new TH1F("fHistNumCC", "# of Calo clusters per event", 100,0,100);
  fHistNumCC->GetXaxis()->SetTitle("# Calo Clusters");
  fHistNumCC->GetYaxis()->SetTitle("Frequency");

  fHistNumTC = new TH1F("fHistNumTC", "# of EMCAL TOF clusters per event", 100,0,100);
  fHistNumTC->GetXaxis()->SetTitle("# EMCAL TOF Clusters");
  fHistNumTC->GetYaxis()->SetTitle("Frequency");

  fHistNumTOFTOT =  new TH1F("fHistNumTC", "# of TOF TOT clusters per event", 100,0,100);
  fHistNumTOFTOT -> GetYaxis()->SetTitle("# TOF TOT Clusters");
  fHistNumTOFTOT -> GetXaxis()->SetTitle("TOF TOT (ns)");

  fHistNumTOFTDC =  new TH1F("fHistNumTC", "# of TOF TDC clusters per event", 100,0,100);
  fHistNumTOFTDC -> GetYaxis()->SetTitle("# TOF TDC Clusters");
  fHistNumTOFTDC -> GetXaxis()->SetTitle("TOF TDC (ns)");

  fOutputList->Add(fHistPt);
  for (int i = 0; i < 100; i++) {fOutputList->Add(fHistEtaPhiCC[i]); fOutputList->Add(fHistEtaPhiTC[i]);}
  fOutputList->Add(fHistNumCC);
  fOutputList->Add(fHistNumTC);
  fOutputList->Add(fHistNumTOFTOT);
  fOutputList->Add(fHistNumTOFTDC);
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
  
  TClonesArray dummy("AliTOFcluster",1000), * tofClusters = &dummy; tofClusterBranch->SetAddress(&tofClusters);
  
  tofClusterTree->GetEvent(0);//  this may be the error
  
  Int_t nClusters = tofClusters->GetEntriesFast(); printf("<INFO>NTOFclusters=%d\n",nClusters);
  int countinEMCal = 0;

  cout<<"TOF Clusters: " <<tofClusters->GetEntriesFast()<<endl;
  // loop over clusters
  //  while (nClusters--) {
  for (Int_t iClusters = 0; iClusters < nClusters; iClusters++) { 
  AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iClusters);
    if(cluster) {
      if (cluster->GetPhi() > 1.35 && cluster->GetPhi() < 3.15) {
	countinEMCal++;
	
	//TOF x,y,z
        Float_t x = cluster->GetR()*TMath::Cos(cluster->GetPhi());
	Float_t y = cluster->GetR()*TMath::Sin(cluster->GetPhi());
	Float_t z = cluster->GetZ();
	
	//	printf("CLUSTER=%d, r=%.2f,phi=%.2f,(x,y,z)=(%.2f,%.2f,%.2f)\n",cluster->GetIndex(),cluster->GetR(),cluster->GetPhi(),x,y,z);
	
	TVector3 vpos(x,y,z);
	Double_t cphi = vpos.Phi();
	if(cphi < 0) cphi +=TMath::TwoPi();
	Double_t ceta = vpos.Eta();
	//	Float_t time =(AliTOFGeometry::TdcBinWidth()*cluster->GetTDC())*1E-3; // in ns
        //Float_t tot = (AliTOFGeometry::TdcBinWidth()*cluster->GetToT())*1E-3;//in ns
	//printf("cphi,ceta = %.3f,%.3f\n",cphi,ceta);
	if (fEvtNum %2000 == 0 && fHistNum < 100) {fHistEtaPhiTC[fHistNum]->Fill(cphi,ceta);}
	//fHistNumTOFTOT-> Fill(tot);
	//fHistNumTOFTDC->Fill(time);
      }
    }
  }
  fHistNumTC->Fill(countinEMCal);
  
  
  printf("<INFO>Ntracks=%d\n",fInputEvent->GetNumberOfTracks());
  
  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++) {
    //AliVParticle* track = fInputEvent->GetTrack(iTracks);
    AliESDtrack* track = (AliESDtrack*)fInputEvent->GetTrack(iTracks);
    if (!track) { printf("ERROR: Could not receive track %d\n", iTracks);continue;}
    //AliTOFcluster* tof = (AliTOFcluster*);
    //tof.fIdx = track->GetTOFcluster();
    //EmCal = new AliESDCaloCluster();
    //EmCal->SetID(track->GetEMCALcluster());
    //cout<<"EMCCal Energy "<<EmCal->E();
    //Int_t tofIndex = track->GetTOFcluster();
    // printf("\t<I>Track %d has tofIndex %d\n",iTracks,tofIndex);
    fHistPt->Fill(track->Pt());
    //fHistDeltaE->Fill(tof->GetADC(), EmCal->E());
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
  //Int_t nCells = clus->GetNCells();
  fHistNumCC->Fill(nclus);
  for (Int_t icl = 0; icl < nclus; icl++) {    
    AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
    Int_t nCells = clus->GetNCells();
    //Float_t energy = clus->E();
    Float_t pos[3];
    clus->GetPosition(pos);
    TVector3 vpos(pos[0],pos[1],pos[2]);
    //cout<< "Pos: " <<"X: "<<pos[0]<<","<< "Y: "<<pos[1]<<"Z: "<<pos[2] <<endl;
    //cout<<"Pos R: "<< sqrt((pos[0]*pos[0])+(pos[1]*pos[1])+(pos[2]*pos[2]))<<endl;
    Double_t cphi = vpos.Phi();
    if(cphi < 0) cphi +=TMath::TwoPi();
    Double_t ceta = vpos.Eta();
    // cout<<"Calo Cluster Phi Value is:  "<<cphi<<endl; cout<<"Calo Cluster Eta Value is: "<<ceta<<endl;
    
    if (fEvtNum%2000 == 0 && fHistNum < 100) {fHistEtaPhiCC[fHistNum]->Fill(cphi, ceta);}
    
    // Time of Flight (TOF)
    Double_t EmCaltof = clus->GetTOF();
    Double_t EmCalEnergy = clus->E();
    fHistTOF->Fill(EmCaltof);
    //Print basic cluster information
    cout << "Cluster: " << icl+1 << "/" << nclus << "Phi: " << 
       cphi*TMath::RadToDeg() << "; Eta: " << ceta << "; NCells: " << nCells << endl;
    cout<< "Number of TofClusters: "<< tofClusters->GetEntriesFast()<<endl<<endl;   /* } calocluster loop */

    //Finding Delta R for Emcal and Tof points
    for (Int_t iToFTrack = 0; iToFTrack <  tofClusters->GetEntriesFast(); iToFTrack++) 
      {//cout<<"iToFTrack = "<<iToFTrack<<endl;
	AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iToFTrack);
	Float_t time =(AliTOFGeometry::TdcBinWidth()*cluster->GetTDC())*1E-3; // in ns
	Float_t tot = (AliTOFGeometry::TdcBinWidth()*cluster->GetToT())*1E-3;//in ns
	fHistNumTOFTOT-> Fill(tot);
	fHistNumTOFTDC->Fill(time);
	//TOF x,y,z
	Float_t TOFx = cluster->GetR()*TMath::Cos(cluster->GetPhi());
	Float_t TOFy = cluster->GetR()*TMath::Sin(cluster->GetPhi());
	Float_t TOFz = cluster->GetZ();
	
	//Double_t p[3];
	//cluster->GetPxPyPz(p);
	//cout<<"Px= "<<p[0]<<"Py="<<p[1]<<"Pz= "<<p[2];
	//Double_t momentum = TMath::Sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
	
	TVector3 TOFvpos(TOFx,TOFy,TOFz);
	Double_t TOFphi = cluster->GetPhi();
	if(TOFphi < 0) {TOFphi +=TMath::TwoPi();}
	//if((momentum != TMath::Abs(p[2]))&&(momentum != 0)){TOFeta = 0.5*TMath::Log((momentum + p[2])/(momentum - p[2]));}
	Double_t TOFeta = TOFvpos.Eta();
	Double_t TOFADC = cluster->GetADC();
	//cout<<"ceta: "<<ceta<<endl;
	//cout<<"TOFeta: "<<TOFeta<<endl;
	//cout<<"cphi: "<<cphi<<endl;
	//cout<<"TOFphi: "<<TOFphi<<endl;
	if (TOFphi < 3.15 && TOFphi > 1.35 && TOFeta < .7 && TOFeta > -.7)
	  {
	    Double_t R1 = ceta - TOFeta;
	//cout<<"R1: "<<R1<<endl;
	    Double_t R2 = cphi-TOFphi;
	//cout<<"R2: "<<R2<<endl;
	    Double_t DeltaR = sqrt(pow(R1, 2) + pow(R2,2));
	    //cout<<"DeltaR: "<<DeltaR<<endl;
	    fHistDeltaR->Fill(DeltaR);
	    fHistDeltaT->Fill(DeltaR, EmCaltof-time);
	    if (DeltaR < .05 && DeltaR > .01){fHistDeltaE->Fill(DeltaR, EmCalEnergy);fHistDeltaADC->Fill(DeltaR, TOFADC);
	      if (TOFADC!=0){fHistEADC->Fill(iToFTrack, EmCalEnergy/TOFADC);}}
	  }
	  }}
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
  if (!fOutputList) {printf("ERROR: Output list not available\n");return;}
  
  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {printf("ERROR: fHistPt not available\n");return;}
 
  //TCanvas *c1 = new TCanvas("AliAnalysisTaskPt","Pt",10,10,510,510);
  //c1->cd(1)->SetLogy();fHistPt->DrawCopy("E");
  
  //TCanvas *c2 = new TCanvas("histo","TOF",10,10,510,510);
  //c2->cd(); fHistTOF->Draw();
  
  //TCanvas *c3 = new TCanvas("histoTDC","TOF TDC",10,10,510,510);
  //c3->cd(); fHistNumTOFTDC->Draw();
  
  //TCanvas *c4 = new TCanvas("histoTOT","TOF TOT",10,10,510,510);
  //c4->cd(); fHistNumTOFTOT->Draw();
  
  TCanvas *c5 = new TCanvas("histoDeltaR", "Tof DeltaR", 10,10,510,510);
  c5->cd(); fHistDeltaR->Draw();

  TCanvas *c9 = new TCanvas("histoDeltaT", "DeltaR-DeltaT", 10,10,510,510);
  c9->cd(); fHistDeltaT->Draw();

  TCanvas *c6 = new TCanvas("histoDeltaE", "TOF-EMACAL Energy", 10,10,510,510);
  c6->cd();fHistDeltaE->Draw();

  TCanvas *c7 = new TCanvas("histoDeltaADC", "TOF ADC-Delta R", 10, 10, 510, 510);
  c7->cd();fHistDeltaADC->Draw();

  TCanvas *c8 = new TCanvas("histoEADC", "E/ADC", 10,10,510,510);
  c8->cd();fHistEADC->Draw();}
