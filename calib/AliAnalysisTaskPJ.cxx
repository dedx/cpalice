
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

#include "AliAnalysisTaskPJ.h"

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPJ)

//________________________________________________________________________
AliAnalysisTaskPJ::AliAnalysisTaskPJ(const char *name) 
: AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fHistPt(0), fHistTOF(0),fHistNumTOFTOT(0),fHistNumTOFTDC(0),
  fHistTotalClusEM(0),fHistTotalClusTOF(0),fHistTotalClusALLTOF(0),fHistUnmatchedClusEM(0),fHistDeltaE(0),fHistUnmatchedClusTOF(0),
    fHistDeltaADC(0),fHistNumCC(0),fHistNumTC(0),fHistUnmatchedEMEnergy(0),fHistUnmatchedTOF(0),
  fHistUnmatchedTOFEnergy(0),fEvtNum(0),fHistNum(0),fHistTOFEMEnergyMatch(0), fHistUnmatchedClusPair(0), fHistVel(0), fHistEMClusDist(0), fHistVelT(0), fHistVelD(0), fHistEMClusMatchedDist(0),fHistMatchedClusNum(0),fHistUnMatchedClusNum(0)


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
AliAnalysisTaskPJ::AliAnalysisTaskPJ() 
  : AliAnalysisTaskSE(), fESD(0), fOutputList(0), fHistPt(0), fHistTOF(0),fHistNumTOFTOT(0),fHistNumTOFTDC(0),
    fHistTotalClusEM(0),fHistTotalClusTOF(0),fHistTotalClusALLTOF(0),fHistUnmatchedClusEM(0),fHistDeltaE(0),fHistUnmatchedClusTOF(0),
    fHistDeltaADC(0),fHistNumCC(0),fHistNumTC(0),fHistUnmatchedEMEnergy(0),fHistUnmatchedTOF(0),
    fHistUnmatchedTOFEnergy(0),fEvtNum(0),fHistNum(0),fHistTOFEMEnergyMatch(0), fHistUnmatchedClusPair(0), fHistVel(0), fHistEMClusDist(0), fHistVelT(0), fHistVelD(0), fHistEMClusMatchedDist(0), fHistMatchedClusNum(0),fHistUnMatchedClusNum(0)


{
  // Default Constructor
  for (int i = 0; i < 100; i++) {
    fHistEtaPhiCC[i] = 0;
    fHistEtaPhiTC[i] = 0;
  }
  
}

//________________________________________________________________________
void AliAnalysisTaskPJ::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  fOutputList = new TList();
  fOutputList->SetOwner(); // otherwise it produces leaks in merging
  
  fHistPt = new TH1F("fHistPt", "P_{T} distribution", 100, 0.1, 10.1);
  fHistPt->GetXaxis()->SetTitle("P_{T} (GeV/c)");
  fHistPt->GetYaxis()->SetTitle("dN/dP_{T} (c/GeV)");
  fHistPt->SetMarkerStyle(kFullCircle);  

  fHistTOF = new TH2F("fHistTOF", "Unmatched TOF Eta and Phi", 100, -.7, .7, 100, 0, .75*TMath::Pi());
  fHistTOF->GetXaxis()->SetTitle("Eta");
  fHistTOF->GetYaxis()->SetTitle("Phi");

  fHistTotalClusEM = new TH1F("fHistTotalClusEM", "Total # of Clusters per Event in EMCAL", 200, 0, 200);
  fHistTotalClusEM->GetXaxis()->SetTitle("# of Clusters per Event");
  fHistTotalClusEM->GetYaxis()->SetTitle("Counts");

  fHistTotalClusTOF = new TH1F("fHistTotalClusTOF", "Total # of Clusters per Event in TOF Restricted to EMCAL ROI", 200, 0, 200);
  fHistTotalClusTOF->GetXaxis()->SetTitle("# of Clusters per Event");
  fHistTotalClusTOF->GetYaxis()->SetTitle("Counts");

  fHistTotalClusALLTOF = new TH1F("fHistTotalClusTOF", "Total # of Clusters per Ev\
ent in TOF", 200, 0, 200);
  fHistTotalClusALLTOF->GetXaxis()->SetTitle("# of Clusters per Event");
  fHistTotalClusALLTOF->GetYaxis()->SetTitle("Counts");


  fHistUnmatchedClusEM = new TH1F("fHistUnmatchedClusEM", "Total # of Unmatched Clusters per Event in EMCAL", 50, 0,50);
  fHistUnmatchedClusEM->GetXaxis()->SetTitle("# of Unmatched Clusters per Event");
  fHistUnmatchedClusEM->GetYaxis()->SetTitle("Counts");

  fHistUnmatchedClusPair = new TH1F("fHistUnmatchedClusPair", "Total # of Unmatched Clusters per Event With a Close TOF Neighbor", 50, 0,50);
  fHistUnmatchedClusPair->GetXaxis()->SetTitle("# of Unmatched Clusters per Event");
  fHistUnmatchedClusPair->GetYaxis()->SetTitle("Counts");


  fHistTOFEMEnergyMatch = new TH2F("fHistTOFEMEnergyMatch", "Delta R vs. Delta E for Unmatched Clusters", 100, 0, .005, 100, 0, 3);
  fHistTOFEMEnergyMatch->GetXaxis()->SetTitle("Delta E(GeV)");
  fHistTOFEMEnergyMatch->GetYaxis()->SetTitle("Delta R");

  fHistUnmatchedClusTOF = new TH1F("fHistUnmatchedClusTOF", "Total # of Unmatched Clusters per Event in TOF Restricted to EMCAL ROI", 50, 0, 50);
  fHistUnmatchedClusTOF->GetXaxis()->SetTitle("# of Unmatched Clusters per Event");
  fHistUnmatchedClusTOF->GetYaxis()->SetTitle("Counts");

  fHistUnmatchedTOFEnergy = new TH1F("fHistUnmatchedTOFEnergy", "TOF Calculated Energy Distribution for Unmatched Clusters with an $e^-$ Mass Assumption", 1000, .0005, 50);
  fHistUnmatchedTOFEnergy->GetXaxis()->SetTitle("Calculated Energy(GeV)");
  fHistUnmatchedTOFEnergy->GetYaxis()->SetTitle("Counts");

  fHistDeltaE = new TH2F("fHistDeltaE", "Unmatched EmCal Eta and Phi", 100, -.7, .7, 100, 0, .75*TMath::Pi());
  fHistDeltaE->GetXaxis()->SetTitle("Eta");
  fHistDeltaE->GetYaxis()->SetTitle("Phi");

  fHistUnmatchedEMEnergy = new TH1F("fHistUnmatchedEMEnergy", "Unmatched EMCAL Cluster Energy Distribution", 1000, .0005, 50);
  fHistUnmatchedEMEnergy->GetXaxis()->SetTitle("Energy(GeV)");
  fHistUnmatchedEMEnergy->GetYaxis()->SetTitle("Counts");

  fHistDeltaADC = new TH2F("fHistDeltaADC", "2D ADC-DeltaR", 10000, 0, 5, 10000, 100, 300);
  fHistDeltaADC->GetXaxis()->SetTitle("Delta R");
  fHistDeltaADC->GetYaxis()->SetTitle("TOF ADC");

  fHistUnmatchedTOF = new TH1F("fHistUnmatchedTOF", "Time of Flight to TOF for Unmatched Clusters", 250000, 0, 500);
  fHistUnmatchedTOF->GetXaxis()->SetTitle("Time(ns)");
  fHistUnmatchedTOF->GetYaxis()->SetTitle("Counts");

  fHistVel = new TH1F("fHistVel", "Velocity Distribution of Unmatched TOF Clusters", 1000, 0, 3);
  fHistVel->GetXaxis()->SetTitle("Velocity (in units of c)");
  fHistVel->GetYaxis()->SetTitle("Counts");
  
  fHistVelT = new TH1F("fHistVelT", "Velocity Distribution of Unmatched TOF Clusters After Time Cut", 1000, 0, 3);
  fHistVelT->GetXaxis()->SetTitle("Velocity (in units of c)");
  fHistVelT->GetYaxis()->SetTitle("Counts");

  fHistVelD = new TH1F("fHistVel", "Velocity Distribution of Unmatched TOF Clusters After Time Cut and Distance Cut", 1000, 0, 3);
  fHistVelD->GetXaxis()->SetTitle("Velocity (in units of c)");
  fHistVelD->GetYaxis()->SetTitle("Counts");

  fHistEMClusMatchedDist = new TH2F("fHistEMClusMatchedDist", "The Eta Phi Distribution of Matched EMCal Clusters", 50,  -.7, .7, 50, 80*TMath::Pi()/180, TMath::Pi());
  fHistEMClusMatchedDist->GetXaxis()->SetTitle("Eta");
  fHistEMClusMatchedDist->GetYaxis()->SetTitle("Phi");


  fHistEMClusDist = new TH2F("fHistEMClusDist", "THe Eta Phi Distribution of Unmatched EMCal Clusters", 50, -.7, .7, 50, 80*TMath::Pi()/180, TMath::Pi());
  fHistEMClusDist->GetXaxis()->SetTitle("Eta");
  fHistEMClusDist->GetYaxis()->SetTitle("Phi");

  fHistMatchedClusNum = new TH1F("fHistMatchedClusNum", "The Number of Towers per Cluster in Matched Clusters", 10, 0, 10);
  fHistMatchedClusNum->GetXaxis()->SetTitle("# of Towers per Cluster");
  fHistMatchedClusNum->GetYaxis()->SetTitle("Counts");

  fHistUnMatchedClusNum = new TH1F("fHistMatchedClusNum", "The Number of Towers per Cluster in Unmatched Clusters", 10, 0, 10);
  fHistUnMatchedClusNum->GetXaxis()->SetTitle("# of Towers per Cluster");
  fHistUnMatchedClusNum->GetYaxis()->SetTitle("Counts");


  for (int i = 0; i < 100; i++) {
    char name[100];sprintf(name,"fHistEtaPhiCC%d",i);
    fHistEtaPhiCC[i] = new TH2F(name,name,100,-.9,.9,100,100,190);
    fHistEtaPhiCC[i]->GetXaxis()->SetTitle("#eta");
    fHistEtaPhiCC[i]->GetYaxis()->SetTitle("#phi"); }


  for (int i = 0; i < 100; i++) {
    char name[100]; sprintf(name,"fHistEtaPhiTC%d",i);
    fHistEtaPhiTC[i] = new TH2F(name,name,100,-.9,.9,100,100,190);
    fHistEtaPhiTC[i]->GetXaxis()->SetTitle("#eta");
    fHistEtaPhiTC[i]->GetYaxis()->SetTitle("#phi");
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
  fOutputList->Add(fHistTOF);
  fOutputList->Add(fHistTotalClusEM);
  fOutputList->Add(fHistTotalClusTOF);
  fOutputList->Add(fHistTotalClusALLTOF);
  fOutputList->Add(fHistUnmatchedClusEM);
  fOutputList->Add(fHistUnmatchedClusPair);
  fOutputList->Add(fHistUnmatchedClusTOF);
  fOutputList->Add(fHistTOFEMEnergyMatch);
  fOutputList->Add(fHistUnmatchedEMEnergy);
  fOutputList->Add(fHistUnmatchedTOFEnergy);
  fOutputList->Add(fHistVel);
  fOutputList->Add(fHistVelT);
  fOutputList->Add(fHistVelD);
  fOutputList->Add(fHistEMClusDist);
  fOutputList->Add(fHistEMClusMatchedDist);
  fOutputList->Add(fHistUnMatchedClusNum);
  fOutputList->Add(fHistMatchedClusNum);
  PostData(1,fOutputList); // important for merging
}

//________________________________________________________________________
void AliAnalysisTaskPJ::UserExec(Option_t *) 
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
  
  Double32_t T0 = esd->GetT0();
  //get reconstructed vertex position
  Double_t vertex_position[3];
  esd->GetVertex()->GetXYZ(vertex_position);
  
  TTree* tofClusterTree = handler->GetTreeR("TOF");
  TTree* Tree = handler->GetTreeR("TOF");
  if (!tofClusterTree) {printf("<WARN>No TOF clusters!\n");return;}
  
  TBranch* tofClusterBranch = tofClusterTree->GetBranch("TOF");
  if (!tofClusterBranch) {printf("Can't get the branch with the TOF digits !\n"); return;}
  
  TClonesArray dummy("AliTOFcluster",1000), * tofClusters = &dummy; tofClusterBranch->SetAddress(&tofClusters);
  
  tofClusterTree->GetEvent(0);//  this may be the error
  
  Int_t nClusters = tofClusters->GetEntriesFast(); printf("<INFO>NTOFclusters=%d\n",nClusters);
  int countinEMCal = 0;
  fHistTotalClusALLTOF->Fill(nClusters);
  cout<<"TOF Clusters: " <<tofClusters->GetEntriesFast()<<endl;
  // loop over clusters
  //  while (nClusters--) {
  for (Int_t iClusters = 0; iClusters < nClusters; iClusters++) { 
  AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iClusters);
    if(cluster) {
      Float_t x = cluster->GetR()*TMath::Cos(cluster->GetPhi());
      Float_t y = cluster->GetR()*TMath::Sin(cluster->GetPhi());
      Float_t z = cluster->GetZ();
      TVector3 vpos(x,y,z);
      Double_t cphi = vpos.Phi();
      if(cphi < 0) cphi +=TMath::TwoPi();
      Double_t ceta = vpos.Eta();

      if (cluster->GetPhi() > 1.35 && cluster->GetPhi() < 3.15 && ceta < 0.7 && ceta > -0.7) {
	countinEMCal++; }
    }
  }
  fHistTotalClusTOF->Fill(countinEMCal);
  
  
  printf("<INFO>Ntracks=%d\n",fInputEvent->GetNumberOfTracks());

  // Track loop to fill a pT spectrum
  for (Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++) {
    //AliVParticle* track = fInputEvent->GetTrack(iTracks);
    AliESDtrack* track = (AliESDtrack*)fInputEvent->GetTrack(iTracks);
    if (!track) { printf("ERROR: Could not receive track %d\n", iTracks);continue;}
    printf("Track %d has been matched to EMCAL %d", track->GetID(), track->GetEMCALcluster());
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
  fHistTotalClusEM->Fill(nclus);
  Bool_t matchedTOF[tofClusters->GetEntriesFast()];
  Bool_t unMatchedTOF[tofClusters->GetEntriesFast()];
  Bool_t unMatchedEmCal[nclus];
  Bool_t matchedEmCal;
  Int_t UnmatchedTOF=0;
  Int_t UnmatchedEM=0;
  Double_t DeltaR=0;
  cout<<"A new loop should happen here"<<"\n";
  matchedEmCal = false;
  Double_t Rthresh=0.1;
  Int_t lastMatched=-999;
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
    
   // if (fEvtNum%2000 == 0 && fHistNum < 100) {fHistEtaPhiCC[fHistNum]->Fill(cphi, ceta);}
    
    // Time of Flight (TOF)
    Double_t EmCaltof = clus->GetTOF();
    Double_t EmCalEnergy = clus->E();
    //fHistTOF->Fill(EmCaltof);
    //Print basic cluster information
    cout << "Cluster: " << icl+1 << "/" << nclus << "Phi: " << 
       cphi*TMath::RadToDeg() << "; Eta: " << ceta << "; NCells: " << nCells << endl;
    cout<< "Number of TofClusters: "<< tofClusters->GetEntriesFast()<<endl<<endl;   /* } calocluster loop */

    //Finding Delta R for Emcal and Tof points
    for (Int_t iToFTrack = 0; iToFTrack <  tofClusters->GetEntriesFast(); iToFTrack++) 
      {//cout<<"iToFTrack = "<<iToFTrack<<endl;
	AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iToFTrack);
	Float_t time =(AliTOFGeometry::TdcBinWidth()*cluster->GetTDC()); // in ps
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
	if (TOFphi < 3.15 && TOFphi > 1.35 && TOFeta < .7 && TOFeta > -.7 && time > 12.33333)
	  {
	    Double_t R1 = TOFeta;
	//cout<<"R1: "<<R1<<endl;
	    Double_t R2 = TOFphi;
	//cout<<"R2: "<<R2<<endl;
	    Double_t Rtof = sqrt(pow(R1, 2) + pow(R2,2));
	    //cout<<"DeltaR: "<<DeltaR<<endl;
        DeltaR = sqrt(pow(ceta-R1, 2) + pow(cphi-R2, 2));
        //These two if statements are for deciding which data you want to see.
	//Putting the fill in the first one will show all matched cluster and the second, all unmatched clusters
	if (abs(DeltaR)<Rthresh && matchedTOF[iToFTrack] == false)
        {
	  if (abs(EmCalEnergy-CalcMass(.511, time, TOFeta, false))<1||abs(EmCalEnergy-CalcMass(139.570, time, TOFeta, false))<1||abs(EmCalEnergy-CalcMass(938.272, time, TOFeta, false))<1||abs(EmCalEnergy-CalcMass(493.667, time, TOFeta, false))<1){
	  if(lastMatched != -999){
	  matchedTOF[lastMatched] = false;
	  }
	  matchedTOF[iToFTrack] = true;
        matchedEmCal = true;
	fHistMatchedClusNum->Fill(clus->GetNCells());
	if(clus->GetNCells()<4){
	fHistEMClusMatchedDist->Fill(ceta, cphi);
	}
	Rthresh = DeltaR;
	lastMatched = iToFTrack;
	  }
        }
        if ( icl+1==nclus && matchedTOF[iToFTrack] == false)
	{
	  fHistUnmatchedTOFEnergy->Fill(CalcMass(.511, time, TOFeta, true));
	  fHistTOFEMEnergyMatch->Fill(CalcMass(.511, time, TOFeta, false)-EmCalEnergy, DeltaR);
	  UnmatchedTOF++;
		fHistUnmatchedTOF->Fill(time);
		if (fEvtNum%2000 == 0 && fHistNum < 100){fHistEtaPhiTC[fHistNum]->Fill(TOFeta, TOFphi*180.0/TMath::Pi());}
	  unMatchedTOF[iToFTrack] = true;
	}

	  }
      }
    if(matchedEmCal == false){
      if (fEvtNum%2000 == 0 && fHistNum < 100){
        fHistEtaPhiCC[fHistNum]->Fill(ceta, cphi*180.0/TMath::Pi());}
      UnmatchedEM++;
      cout<<"HEY LISTEN"<<UnmatchedEM<<"\n";
      unMatchedEmCal[icl]=true;
      fHistEMClusDist->Fill(ceta, cphi);
      fHistUnMatchedClusNum->Fill(clus->GetNCells());
      fHistUnmatchedEMEnergy->Fill(EmCalEnergy);}

  }
  Int_t ClosePair=0;
  for(Int_t iToFTrack = 0; iToFTrack<tofClusters->GetEntriesFast(); iToFTrack++)
    {
      AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iToFTrack);
      Float_t TOFx = cluster->GetR()*TMath::Cos(cluster->GetPhi());
      Float_t TOFy = cluster->GetR()*TMath::Sin(cluster->GetPhi());
      Float_t TOFz = cluster->GetZ();

      TVector3 TOFvpos(TOFx,TOFy,TOFz);
      Double_t TOFeta = TOFvpos.Eta();
      Double32_t time =(AliTOFGeometry::TdcBinWidth()*cluster->GetTDC\
			())*1E-3; // in ns

      if(unMatchedTOF[iToFTrack]==true)
	{
	  for(Int_t iToFTrack2 = iToFTrack+1; iToFTrack2<tofClusters->GetEntriesFast(); iToFTrack2++)
	    {
	      AliTOFcluster *cluster2=(AliTOFcluster*)tofClusters->UncheckedAt(iToFTrack2);
	      Float_t TOFx2 = cluster2->GetR()*TMath::Cos(cluster2->GetPhi());
	      Float_t TOFy2 = cluster2->GetR()*TMath::Sin(cluster2->GetPhi());
	      Float_t TOFz2 = cluster2->GetZ();

	      TVector3 TOFvpos2(TOFx2,TOFy2,TOFz2);
	      Double_t TOFeta2 = TOFvpos2.Eta();
	      Double32_t time2 =(AliTOFGeometry::TdcBinWidth()*cluster2->GetTDC\
				())*1E-3; // in ns
	      Double_t DeltaRpair = sqrt(pow(TOFeta2-TOFeta,2)+pow(cluster2->GetPhi()-cluster->GetPhi(),2));
	      Double_t DeltaT = abs(time2-time);
	      if(unMatchedTOF[iToFTrack2]==true && DeltaRpair<0.001 && DeltaT<10)
		{
		  ClosePair++;
		}
	    }
	}
    }
  fHistUnmatchedClusPair->Fill(ClosePair);

//This is all TOF e- Mass Assumotion Stuff
  for(Int_t icl=0; icl<nclus; icl++)
    {
      AliVCluster* clus = (AliVCluster*)caloClusters->At(icl);
      if(unMatchedEmCal[icl]==true)
	{
	  Double_t energy = clus->E();
	  for(Int_t iToFTrack = 0; iToFTrack<tofClusters->GetEntriesFast(); iToFTrack++)
	      {
		AliTOFcluster *cluster=(AliTOFcluster*)tofClusters->UncheckedAt(iToFTrack);
		Float_t TOFx = cluster->GetR()*TMath::Cos(cluster->GetPhi());
		Float_t TOFy = cluster->GetR()*TMath::Sin(cluster->GetPhi());
		Float_t TOFz = cluster->GetZ();

		TVector3 TOFvpos(TOFx,TOFy,TOFz);
		Double_t TOFeta = TOFvpos.Eta();
		Double32_t time =(AliTOFGeometry::TdcBinWidth()*cluster->GetTDC()); // in ps
		Float_t tot = (AliTOFGeometry::TdcBinWidth()*cluster->GetToT())*1E-3;//in ns
		if(unMatchedTOF[iToFTrack]==true)
		  {
		    Double_t elecenergycalc = CalcMass(.511, time, TOFeta, false);
		    /*
		    if (time!=0){
		      Double_t c = (3.00*1E8);
		      Double_t veloc = ((3.70/(TMath::Sin(2*TMath::ATan(TMath::Exp(-TOFeta)))))/((time)*1E-9))/c;
		     Double_t elecmass = (.511);
		     if((pow((veloc/c),2))<1){
		       Double_t elecenergycalc = sqrt(pow((elecmass),2)+pow((elecmass*veloc/sqrt(1-pow((veloc/c),2))),2));*/
		    //fHistUnmatchedTOFEnergy->Fill(elecenergycalc);
		    // fHistTOFEMEnergyMatch->Fill(elecenergycalc-energy, DeltaR);
		     
		  
		    }
	      }
	      }
    }
  cout<<"Unmatched Emcal"<<UnmatchedEM<<"Unmatched TOT"<<UnmatchedTOF;
  fHistUnmatchedClusTOF->Fill(UnmatchedTOF);
  fHistUnmatchedClusEM->Fill(UnmatchedEM);
  //clean up to avoid mem leaks
  delete tofClusterTree;
  //Keep every 2000th event's info in histogram
  if (fEvtNum%2000 == 0) fHistNum++;
  //increment event counter
  fEvtNum++;
  PostData(1, fOutputList);
  
}

//------------------------------------------------------------------------
Double_t AliAnalysisTaskPJ::CalcMass(Double_t MassHyp, Double_t time, Double_t eta, Bool_t Histo)
{
  if(Histo){
  Double_t c = 0.000299792458;
  Double_t dist = (3.70/(TMath::Sin(2*TMath::ATan(TMath::Exp(-eta)))));
  Double_t veloc = (dist/time)/c;
  fHistVel->Fill(veloc);
  if (time>15491){
    Double_t c = 0.000299792458;
    Double_t dist = (3.70/(TMath::Sin(2*TMath::ATan(TMath::Exp(-eta)))));
    Double_t veloc = (dist/time)/c;
    fHistVelT->Fill(veloc);
  }
  if (time>15491){
    Double_t c = 0.000299792458;
    Double_t dist = (3.70/(TMath::Sin(2*TMath::ATan(TMath::Exp(-eta)))));
    if(abs(dist)<4.644){
      Double_t veloc = (dist/time)/c;
      fHistVelD->Fill(veloc);}}
  PostData(1, fOutputList);}
  if (time>15491){
    Double_t c = 0.000299792458;
    Double_t dist = (3.70/(TMath::Sin(2*TMath::ATan(TMath::Exp(-eta)))));
    Double_t veloc = (dist/time)/c;
    if(abs(veloc)<1){
    cout<<"Velocity: "<<veloc<<"Distance: "<<dist;
      Double_t energycalc = sqrt(pow(MassHyp*1E-3,2)+pow((MassHyp*1E-3*veloc/sqrt(1-pow((veloc),2))),2));
      
    return energycalc;
    }}
  cout<<"The velocity is too high or time is zero";
  return -999;
}
//________________________________________________________________________
void AliAnalysisTaskPJ::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  
  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {printf("ERROR: Output list not available\n");return;}
  
  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {printf("ERROR: fHistPt not available\n");return;}
 
  /*  TCanvas *h = new TCanvas("histoEtaPhi", "EtaPhi", 10,10,510,510);
  h->cd();
  char str[32];
  for(int i=0;i<101;i++){
  sprintf(str, "./EtaPhiHistos/TofEtaPhi%d.pdf", i);
  h->Clear();
  fHistEtaPhiTC[i]->Draw();
  h->SaveAs(str);
  sprintf(str, "./EtaPhiHistos/EMCalEtaPhi%d.pdf", i);
  h->Clear(); 
  fHistEtaPhiCC[i]->Draw();
  h->SaveAs(str);
  }

  */  
  /*TCanvas *c1 = new TCanvas("AliAnalysisTaskPJ","Pt",10,10,510,510);
  c1->cd();fHistTotalClusALLTOF->Draw();
  
  TCanvas *c8 = new TCanvas("histo","TOF",10,10,510,510);
  c8->cd(); fHistUnmatchedClusPair->Draw();
  
  TCanvas *c3 = new TCanvas("histototalclusem", "Total # of Clusters per Event in Emcal", 10,10,510,510);
  c3->cd(); fHistTotalClusEM->Draw();
  
  TCanvas *c2 = new TCanvas("histoTOFunmatchedclus","Total # of Unmatched Clusters per Event in TOF Restricted to EMCAL ROI",10,10,510,510);
  c2->cd(); fHistUnmatchedClusTOF->Draw();
  
  TCanvas *c4 = new TCanvas("histoEMEnergy","Energy Distribution for Unmatched Clusters in Emcal",10,10,510,510);
  c4->cd(); fHistUnmatchedEMEnergy->Draw();
  
  TCanvas *c5 = new TCanvas("histoTOFtotalclus", "Total # of Clusters per Event in TOF Restricted to EMCAL ROI", 10,10,510,510);
  c5->cd(); fHistTotalClusTOF->Draw();
  
  TCanvas *c9 = new TCanvas("histoEMUnmatchedclus", "Total # of Unmatched Clusters per event in Emcal", 10,10,510,510);
  c9->cd(); fHistUnmatchedClusEM->Draw();
  
  TCanvas *c6 = new TCanvas("histoTOFEnergy", "Calculated Energy Distribution for Unmatched Clusters in TOF based on an $e^-$ Mass Assumption", 10,10,510,510);
  c6->cd();fHistUnmatchedTOFEnergy->Draw();
  
  TCanvas *c7 = new TCanvas("histoUnmatchedDeltaRDeltaE", "Delta R vs. Delta E for Unmatched Clusters", 10, 10, 510, 510);
  c7->cd();fHistTOFEMEnergyMatch->Draw("COLZ");
*/  
}
