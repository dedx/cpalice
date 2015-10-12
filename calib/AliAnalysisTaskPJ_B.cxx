
#include <Rtypes.h>
#include <TROOT.h>

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TGeoManager.h>
#include <TTree.h>

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
#include "AliTOFRecoParam.h"
#include "AliTOFtrack.h"
#include "AliESDTOFCluster.h"

#include "AliAnalysisTaskPJ_B.h"
    
extern TGeoManager *gGeoManager;

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing
// Reviewed: A.Gheata (19/02/10)

ClassImp(AliAnalysisTaskPJ_B)

//________________________________________________________________________
AliAnalysisTaskPJ_B::AliAnalysisTaskPJ_B(const char* name, TChain* pass1) 
: AliAnalysisTaskSE(name), fESD(0), fOutputList(0), fHistPt(0), fHistTC(0), fWrittenInPos(0x0)


{
  // Constructor
  
  fWrittenInPos = new Int_t[77777];  
    
    
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
    
  
  // Get ESD from SPC Calo run
  AliESDInputHandlerRP *handler = 
    dynamic_cast<AliESDInputHandlerRP*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if (!handler) { printf("No RP handler\n"); return; }
  AliESDEvent *esd  = handler->GetEvent();
  if (!esd) { printf("No AliESDEvent\n"); return; }
  
  //get reconstructed vertex position
  Double_t vertex_position[3];
  esd->GetVertex()->GetXYZ(vertex_position);
    
  //Instantiate CDB and Geometry
  AliCDBManager *cdbManager = AliCDBManager::Instance();
  cdbManager->SetDefaultStorage("alien://folder=/alice/data/2011/OCDB");
  cdbManager->SetRun(esd->GetRunNumber());  
  AliGeomManager::LoadGeometry();
    
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
  
  //Get ESD from reconstructed pass  
  AliESDEvent *RecoESD = new AliESDEvent();
  RecoESD->ReadFromTree((TTree*)pass1tree);
  
  //Match Timestamps, TODO::::::Match other things from email!!
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
    
    //Match TOF clusters and tracks
    if(foundOne){
    pass1tree->GetEntry(mEv);
    TClonesArray* fTracks = RecoESD->fTracks;
        
    //-------------------------------------------------------------------------
    cout<<"Matching Tracks";
    Double_t* fTimesAr[9];// = new Double_t[9];
    Float_t* fTrackPos[4];// = new Float_t[4];
    Int_t kMaxCluster = tofClusters->GetEntries();
    TClonesArray* tmpESDclus;
    Int_t fnunmatch, fnmatch;
    
    for (Int_t ii=0; ii<kMaxCluster; ii++){
        fWrittenInPos[ii] = -1;
    }

    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++)
        fTimesAr[isp] = NULL;

    for (Int_t ii=0; ii<4; ii++)
        fTrackPos[ii] = NULL;
  //Instantiate reconstruction parameters
    static Float_t detDepth=18.;
    static Float_t padDepth=0.5;
    
    const Float_t kSpeedOfLight = 2.99792458e-2;
    
    Float_t dY=AliTOFGeometry::XPad();
    Float_t dZ=AliTOFGeometry::ZPad();
    
    AliTOFRecoParam* TOFrp = new AliTOFRecoParam();
    
    Float_t sensRadius = TOFrp->GetSensRadius();
    Float_t stepSize = TOFrp->GetStepSize();
    Float_t scaleFact = TOFrp->GetWindowScaleFact();
    Float_t dyMax = TOFrp->GetWindowSizeMaxY();
    Float_t dzMax = TOFrp->GetWindowSizeMaxZ();
    Float_t dCut = 10.;
    Double_t maxChi2 = TOFrp->GetMaxChi2TRD();
    Bool_t timeWalkCorr = TOFrp->GetTimeWalkCorr();
    
    Int_t fNsteps = 0;

    
    if(Int_t(detDepth/stepSize) > fNsteps){ 
        fNsteps =(Int_t)(detDepth/stepSize);

        for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
            if(fTimesAr[isp]) delete[] fTimesAr[isp];
        }

        for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
            fTimesAr[isp] = new Double_t[fNsteps];
        }

        for (Int_t ii=0; ii<4; ii++) if(fTrackPos[ii]) delete [] fTrackPos[ii];

        for (Int_t ii=0; ii<4; ii++) fTrackPos[ii] = new Float_t[fNsteps];
    }

    const Int_t kNclusterMax = 1000;
    TGeoHMatrix global[kNclusterMax];
    Int_t clind[kNclusterMax];
    Bool_t isClusterMatchable[kNclusterMax];
  //-------------------------------------------------------------------------
    
    AliTOFtrack trackTOFin;
    
    //The matching loop
    for (Int_t iseed=0; iseed<fTracks->GetEntriesFast(); iseed++) {
        AliESDtrack *t =(AliESDtrack*)fTracks->At(iseed); // ciao replace with loop on ESD + kTOFin
        if( (t->GetStatus()&AliESDtrack::kTOFin) == 0 ) continue;

        trackTOFin = *t;

        for (Int_t ii=0; ii<4; ii++)
          for (Int_t jj=0; jj<fNsteps; jj++) fTrackPos[ii][jj]=0.;

        for (Int_t ii=0; ii<kNclusterMax; ii++) clind[ii]=-1;
        for (Int_t ii=0; ii<kNclusterMax; ii++) global[ii] = 0x0;
        for (Int_t ii=0; ii<kNclusterMax; ii++) isClusterMatchable[ii] = kFALSE;	  	

        Double_t timesOr[AliPID::kSPECIESC]; t->GetIntegratedTimes(timesOr,AliPID::kSPECIESC); // in ps
        
        // Determine a window around the track
        Double_t x,par[5]; 
        trackTOFin.GetExternalParameters(x,par);
        Double_t cov[15]; 
        trackTOFin.GetExternalCovariance(cov);

        if (cov[0]<0. || cov[2]<0.) {
          AliWarning(Form("Very strange track (%d)! At least one of its covariance matrix diagonal elements is negative!",iseed));
          continue;
        }

        Double_t dphi=
          scaleFact*
          ((5*TMath::Sqrt(TMath::Abs(cov[0])) + 0.5*dY + 2.5*TMath::Abs(par[2]))/sensRadius); 
        Double_t dz=
           scaleFact*
           (5*TMath::Sqrt(TMath::Abs(cov[2])) + 0.5*dZ + 2.5*TMath::Abs(par[3]));

        Double_t phi=TMath::ATan2(par[0],x) + trackTOFin.GetAlpha();
        if (phi<-TMath::Pi())phi+=2*TMath::Pi();
        if (phi>=TMath::Pi())phi-=2*TMath::Pi();
        Double_t z=par[1];
        
    //upper limit on window's size.
    if (dz> dzMax) dz=dzMax;
    if (dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;


    // find the clusters in the window of the track
    Int_t nc=0;
    Int_t k;
     
    Int_t n = tofClusters->GetEntriesFast();
    if (tofClusters){
        if (n==0) k=0;
        if (tofClusters->UncheckedAt(0)){
            if (z <= (dynamic_cast<AliCluster *>(tofClusters->UncheckedAt(0)))->GetZ()) k=0;}
        if (tofClusters->UncheckedAt(n-1)){
            if (z > (dynamic_cast<AliCluster *>(tofClusters->UncheckedAt(n-1)))->GetZ()) k=n;}
        Int_t b=0, e=n-1, m=(b+e)/2;
        for (; b<e; m=(b+e)/2) {
            if (tofClusters->UncheckedAt(m)){
                if (z > (dynamic_cast<AliCluster *>(tofClusters->UncheckedAt(m)))->GetZ()) b=m+1;
                else e=m; }
  }
  k=m;
    }
        
    for (; k<tofClusters->GetEntriesFast(); k++) {

      if (nc>=kNclusterMax) {AliWarning("No more matchable clusters can be stored! Please, increase the corresponding vectors size.");
      break;
      }

      AliESDTOFCluster *c=(AliESDTOFCluster *) tofClusters->UncheckedAt(k);
        if(c){
      if (c->GetZ() > z+dz) break;
      if (!c->GetStatus()) {
	AliDebug(1,"Cluster in channel declared bad!");
	continue; // skip bad channels as declared in OCDB
      }


      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;

      Double_t yc=(c->GetPhi() - trackTOFin.GetAlpha())*c->GetR();
      Double_t p[2]={yc, c->GetZ()};
      Double_t cov2[3]= {dY*dY/12., 0., dZ*dZ/12.};
      if (trackTOFin.AliExternalTrackParam::GetPredictedChi2(p,cov2) > maxChi2)continue;

      clind[nc] = k;      
      Char_t path[200];
      Int_t ind[5]; AliTOFGeometry::GetVolumeIndices(c->GetTOFchannel(),ind);
      AliTOFGeometry::GetVolumePath(ind,path);
      gGeoManager->cd(path);
      global[nc] = *gGeoManager->GetCurrentMatrix();
      nc++;
        }
    }


    if (nc == 0 ) {
      AliDebug(1,Form("No available clusters for the track number %d",iseed));
      fnunmatch++;
      continue;
    }
        
    //start fine propagation 

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<fNsteps; istep++){ 
      
      // First of all, propagate the track...
      Float_t xs = AliTOFGeometry::RinTOF()+istep*stepSize;
      if (!(trackTOFin.PropagateTo(xs))) break;

      //  ...and then, if necessary, rotate the track
      Double_t ymax = xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
      Double_t ysect = trackTOFin.GetY();
      if (ysect > ymax) {
	if (!(trackTOFin.Rotate(AliTOFGeometry::GetAlpha()))) break;
      } else if (ysect <-ymax) {
	if (!(trackTOFin.Rotate(-AliTOFGeometry::GetAlpha()))) break;
      }

      Double_t mom = trackTOFin.P();

      if(istep == 0){
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = stepSize/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }
      else{
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = fTimesAr[isp][nStepsDone-1] + (trackTOFin.GetIntegratedLength()-fTrackPos[3][nStepsDone-1])/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }

      // store the running point (Globalrf) - fine propagation     

      Double_t r[3]; trackTOFin.GetXYZ(r);
      fTrackPos[0][nStepsDone]= (Float_t) r[0];
      fTrackPos[1][nStepsDone]= (Float_t) r[1];
      fTrackPos[2][nStepsDone]= (Float_t) r[2];   
      fTrackPos[3][nStepsDone]= trackTOFin.GetIntegratedLength();

      nStepsDone++;
      AliDebug(3,Form(" current step %d (%d) - nStepsDone=%d",istep,fNsteps,nStepsDone));
    }

    if ( nStepsDone == 0 ) {
      AliDebug(1,Form(" No track points for track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(3,Form(" Number of steps done for the track number %d: %d",iseed,nStepsDone));

    if(nc){
      for (Int_t i=0; i<nc; i++) isClusterMatchable[i] = kFALSE;	  	
    }
        
    Int_t nfound = 0;
    Bool_t accept = kFALSE;
    for (Int_t istep=0; istep<nStepsDone; istep++) {
      Float_t ctrackPos[3];     
      ctrackPos[0] = fTrackPos[0][istep];
      ctrackPos[1] = fTrackPos[1][istep];
      ctrackPos[2] = fTrackPos[2][istep];

      //now see whether the track matches any of the TOF clusters            

      Float_t dist3d[3]={0.,0.,0.};
      accept = kFALSE;

      for (Int_t i=0; i<nc; i++) {

        AliTOFGeometry::IsInsideThePad((TGeoHMatrix*)(&global[i]),ctrackPos,dist3d);

        // check multiple hit cases
        AliESDTOFCluster *cmatched=(AliESDTOFCluster *) tofClusters->UncheckedAt(clind[i]);

        if(cmatched->GetNTOFhits() > 1){ // correct residual for mean position of the clusters (w.r.t. the first pad/hit)
          Float_t zmain = cmatched->GetTOFchannel(0)/48;
          Float_t xmain = cmatched->GetTOFchannel(0)%48;
          for(Int_t ihit=1;ihit < cmatched->GetNTOFhits();ihit++){
            Float_t deltaz = (cmatched->GetTOFchannel(ihit)/48 - zmain) * 3.5;
            Float_t deltax = (cmatched->GetTOFchannel(ihit)%48 - xmain) * 2.5;
            dist3d[0] -= deltax / cmatched->GetNTOFhits();
            dist3d[2] -= deltaz / cmatched->GetNTOFhits();
          }
        }
          
        // ***** NEW *****
        /* if track is inside this cluster set flags which will then
         * inhibit to add track points for the other clusters */
    Float_t yLoc = dist3d[1];
    Float_t rLoc = TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[2]*dist3d[2]);
	accept = (TMath::Abs(yLoc)<padDepth*0.5 && rLoc<dCut);

	//***** NEW *****
	/* add point everytime that:
	 * - the tracks is within dCut from the cluster
	 */
        if (accept) {

	  Double_t timesCurrent[AliPID::kSPECIESC];
	  AliDebug(3,Form(" Momentum for track %d -> %f", iseed,t->P()));
	  for (Int_t j=0;j<AliPID::kSPECIESC;j++) {
	    timesCurrent[j] = timesOr[j] + fTimesAr[j][istep];
	  }


	  if (TMath::Abs(dist3d[1])<stepSize && !isClusterMatchable[i]) {
	    isClusterMatchable[i] = kTRUE;
	    
	    Int_t currentpos = tmpESDclus->GetEntriesFast(); // position of cluster in ESD
	    if(fWrittenInPos[clind[i]] != -1){
	      currentpos = fWrittenInPos[clind[i]];
	      cmatched = (AliESDTOFCluster *) tmpESDclus->At(currentpos); // update the new one in the ESDEvent
	    }
	    else{ // add as a new cluster in the ESD TClonesArray
	      AliESDTOFCluster *clnew =  new( (*tmpESDclus)[currentpos] ) AliESDTOFCluster(*cmatched);
	      clnew->SetESDID(currentpos);

	     /* // remap also TOF hit in the filtered array
	      for(Int_t ii=0;ii < cmatched->GetNTOFhits();ii++){
            Int_t index = cmatched->GetHitIndex(ii);
            AliESDTOFHit *hitOld = (AliESDTOFHit *) esdTOFHitArr->At(index);
            Int_t index_new = fHitsESD->GetEntriesFast();
            AliESDTOFHit *hitNew = new( (*fHitsESD)[index_new] ) AliESDTOFHit(*hitOld);
            hitNew->SetESDTOFClusterIndex(currentpos);
            clnew->SetHitIndex(ii,index_new);
	      }*/

	      fWrittenInPos[clind[i]] = currentpos;
	      cmatched = clnew; // update the new one added to the ESDEvent
	    }

	    if(cmatched->GetNMatchableTracks() < AliESDTOFCluster::kMaxMatches){
	      cmatched->Update(t->GetID(),dist3d[0],dist3d[1],dist3d[2],fTrackPos[3][istep],timesCurrent);//x,y,z -> tracking RF
	      t->AddTOFcluster(currentpos);
	      t->SetStatus(AliESDtrack::kTOFout);
	    }
	  }
          AliDebug(2,Form(" dist3dLoc[0] = %f, dist3dLoc[1] = %f, dist3dLoc[2] = %f ",dist3d[0],dist3d[1],dist3d[2]));

          nfound++;
       // ***** NEW *****
        }//end if accept
        
      } //end for on the clusters
    } //end for on the steps     


    if (nfound == 0 ) {
      AliDebug(1,Form(" No matchable track points for the track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" Number of track points for the track number %d: %d",iseed,nfound));

    Int_t nMatchedClusters = t->GetNTOFclusters();
 
    if (nMatchedClusters==0) {
      AliDebug(1,Form("Reconstructed track %d doesn't match any TOF cluster", iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" %d - matched (%d)",iseed,nMatchedClusters));

    fnmatch++;

    /*
    AliTOFcluster cTOF = AliTOFcluster(volIdClus,
    (Float_t)posClus[0],(Float_t)posClus[1],(Float_t)posClus[2],
    (Float_t)covClus[0],(Float_t)covClus[1],(Float_t)covClus[2],
    (Float_t)covClus[3],(Float_t)covClus[4],(Float_t)covClus[5],
    tofLabels,volIndices,parClu,kTRUE,index[i]);

    // Fill the track residual histograms.
    FillResiduals(trackTOFin,c,kFALSE);
    */
  } // loop on fSeeds
        
        
  //---------------------------------------------------------------------
    
    for(Int_t itof = 0; itof<tofClusters->GetEntriesFast(); itof++){
            AliTOFcluster* tofClus = (AliTOFcluster*)tofClusters->At(itof);
        
            fHistTC->Fill(tofClus->GetESDID());
            }
     /*
    for(Int_t icl = 0; icl<RecoESD->GetNumberOfTracks();icl++){
        AliESDtrack* trk = (AliESDtrack*)RecoESD->GetTrack(icl);
        cout<<icl;
        if(!trk){continue;}
        AliExternalTrackParam* trkpar = (AliExternalTrackParam*)trk->GetOuterParam();
        if(!trkpar){continue;}
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
    }*/
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

void AliAnalysisTaskPJ_B::MatchTracks(TClonesArray* ESDtrk, TClonesArray* ESDclus)
{
    //-------------------------------------------------------------------------
    cout<<"Matching Tracks";
    Double_t* fTimesAr[9];// = new Double_t[9];
    Float_t* fTrackPos[4];// = new Float_t[4];
    Int_t kMaxCluster = ESDclus->GetEntries();
    TClonesArray* tmpESDclus;
    Int_t fnunmatch, fnmatch;
    
    for (Int_t ii=0; ii<kMaxCluster; ii++){
        fWrittenInPos[ii] = -1;
    }

    for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++)
        fTimesAr[isp] = NULL;

    for (Int_t ii=0; ii<4; ii++)
        fTrackPos[ii] = NULL;
  //Instantiate reconstruction parameters
    static Float_t detDepth=18.;
    static Float_t padDepth=0.5;
    
    const Float_t kSpeedOfLight = 2.99792458e-2;
    
    Float_t dY=AliTOFGeometry::XPad();
    Float_t dZ=AliTOFGeometry::ZPad();
    
    AliTOFRecoParam* TOFrp = new AliTOFRecoParam();
    
    Float_t sensRadius = TOFrp->GetSensRadius();
    Float_t stepSize = TOFrp->GetStepSize();
    Float_t scaleFact = TOFrp->GetWindowScaleFact();
    Float_t dyMax = TOFrp->GetWindowSizeMaxY();
    Float_t dzMax = TOFrp->GetWindowSizeMaxZ();
    Float_t dCut = 10.;
    Double_t maxChi2 = TOFrp->GetMaxChi2TRD();
    Bool_t timeWalkCorr = TOFrp->GetTimeWalkCorr();
    
    Int_t fNsteps = 0;

    
    if(Int_t(detDepth/stepSize) > fNsteps){ 
        fNsteps =(Int_t)(detDepth/stepSize);

        for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
            if(fTimesAr[isp]) delete[] fTimesAr[isp];
        }

        for(Int_t isp=0;isp < AliPID::kSPECIESC;isp++){
            fTimesAr[isp] = new Double_t[fNsteps];
        }

        for (Int_t ii=0; ii<4; ii++) if(fTrackPos[ii]) delete [] fTrackPos[ii];

        for (Int_t ii=0; ii<4; ii++) fTrackPos[ii] = new Float_t[fNsteps];
    }

    const Int_t kNclusterMax = 1000;
    TGeoHMatrix global[kNclusterMax];
    Int_t clind[kNclusterMax];
    Bool_t isClusterMatchable[kNclusterMax];
  //-------------------------------------------------------------------------
    
    AliTOFtrack trackTOFin;
    
    //The matching loop
    for (Int_t iseed=0; iseed<ESDtrk->GetEntriesFast(); iseed++) {
        AliESDtrack *t =(AliESDtrack*)ESDtrk->At(iseed); // ciao replace with loop on ESD + kTOFin
        if( (t->GetStatus()&AliESDtrack::kTOFin) == 0 ) continue;

        trackTOFin = *t;

        for (Int_t ii=0; ii<4; ii++)
          for (Int_t jj=0; jj<fNsteps; jj++) fTrackPos[ii][jj]=0.;

        for (Int_t ii=0; ii<kNclusterMax; ii++) clind[ii]=-1;
        for (Int_t ii=0; ii<kNclusterMax; ii++) global[ii] = 0x0;
        for (Int_t ii=0; ii<kNclusterMax; ii++) isClusterMatchable[ii] = kFALSE;	  	

        Double_t timesOr[AliPID::kSPECIESC]; t->GetIntegratedTimes(timesOr,AliPID::kSPECIESC); // in ps
        
        // Determine a window around the track
        Double_t x,par[5]; 
        trackTOFin.GetExternalParameters(x,par);
        Double_t cov[15]; 
        trackTOFin.GetExternalCovariance(cov);

        if (cov[0]<0. || cov[2]<0.) {
          AliWarning(Form("Very strange track (%d)! At least one of its covariance matrix diagonal elements is negative!",iseed));
          continue;
        }

        Double_t dphi=
          scaleFact*
          ((5*TMath::Sqrt(TMath::Abs(cov[0])) + 0.5*dY + 2.5*TMath::Abs(par[2]))/sensRadius); 
        Double_t dz=
           scaleFact*
           (5*TMath::Sqrt(TMath::Abs(cov[2])) + 0.5*dZ + 2.5*TMath::Abs(par[3]));

        Double_t phi=TMath::ATan2(par[0],x) + trackTOFin.GetAlpha();
        if (phi<-TMath::Pi())phi+=2*TMath::Pi();
        if (phi>=TMath::Pi())phi-=2*TMath::Pi();
        Double_t z=par[1];
        
    //upper limit on window's size.
    if (dz> dzMax) dz=dzMax;
    if (dphi*sensRadius> dyMax) dphi=dyMax/sensRadius;


    // find the clusters in the window of the track
    Int_t nc=0;
    Int_t k;
    for (Int_t k=FindClusterIndex(z-dz, ESDclus); k<ESDclus->GetEntriesFast(); k++) {

      if (nc>=kNclusterMax) {AliWarning("No more matchable clusters can be stored! Please, increase the corresponding vectors size.");
      break;
      }

      AliESDTOFCluster *c=(AliESDTOFCluster *) ESDclus->UncheckedAt(k);
        if(c){
      if (c->GetZ() > z+dz) break;
      if (!c->GetStatus()) {
	AliDebug(1,"Cluster in channel declared bad!");
	continue; // skip bad channels as declared in OCDB
      }


      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2.*TMath::Pi();
      if (TMath::Abs(dph)>dphi) continue;

      Double_t yc=(c->GetPhi() - trackTOFin.GetAlpha())*c->GetR();
      Double_t p[2]={yc, c->GetZ()};
      Double_t cov2[3]= {dY*dY/12., 0., dZ*dZ/12.};
      if (trackTOFin.AliExternalTrackParam::GetPredictedChi2(p,cov2) > maxChi2)continue;

      clind[nc] = k;      
      Char_t path[200];
      Int_t ind[5]; AliTOFGeometry::GetVolumeIndices(c->GetTOFchannel(),ind);
      AliTOFGeometry::GetVolumePath(ind,path);
      gGeoManager->cd(path);
      global[nc] = *gGeoManager->GetCurrentMatrix();
      nc++;
        }
    }


    if (nc == 0 ) {
      AliDebug(1,Form("No available clusters for the track number %d",iseed));
      fnunmatch++;
      continue;
    }
        
    //start fine propagation 

    Int_t nStepsDone = 0;
    for( Int_t istep=0; istep<fNsteps; istep++){ 
      
      // First of all, propagate the track...
      Float_t xs = AliTOFGeometry::RinTOF()+istep*stepSize;
      if (!(trackTOFin.PropagateTo(xs))) break;

      //  ...and then, if necessary, rotate the track
      Double_t ymax = xs*TMath::Tan(0.5*AliTOFGeometry::GetAlpha());
      Double_t ysect = trackTOFin.GetY();
      if (ysect > ymax) {
	if (!(trackTOFin.Rotate(AliTOFGeometry::GetAlpha()))) break;
      } else if (ysect <-ymax) {
	if (!(trackTOFin.Rotate(-AliTOFGeometry::GetAlpha()))) break;
      }

      Double_t mom = trackTOFin.P();

      if(istep == 0){
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = stepSize/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }
      else{
	for(Int_t isp=0;isp<AliPID::kSPECIESC;isp++){
	  Double_t mass=AliPID::ParticleMass(isp);
	  Double_t momz = mom*AliPID::ParticleCharge(isp);
	  fTimesAr[isp][nStepsDone] = fTimesAr[isp][nStepsDone-1] + (trackTOFin.GetIntegratedLength()-fTrackPos[3][nStepsDone-1])/kSpeedOfLight*TMath::Sqrt(momz*momz+mass*mass)/momz;
	}
      }

      // store the running point (Globalrf) - fine propagation     

      Double_t r[3]; trackTOFin.GetXYZ(r);
      fTrackPos[0][nStepsDone]= (Float_t) r[0];
      fTrackPos[1][nStepsDone]= (Float_t) r[1];
      fTrackPos[2][nStepsDone]= (Float_t) r[2];   
      fTrackPos[3][nStepsDone]= trackTOFin.GetIntegratedLength();

      nStepsDone++;
      AliDebug(3,Form(" current step %d (%d) - nStepsDone=%d",istep,fNsteps,nStepsDone));
    }

    if ( nStepsDone == 0 ) {
      AliDebug(1,Form(" No track points for track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(3,Form(" Number of steps done for the track number %d: %d",iseed,nStepsDone));

    if(nc){
      for (Int_t i=0; i<nc; i++) isClusterMatchable[i] = kFALSE;	  	
    }
        
    Int_t nfound = 0;
    Bool_t accept = kFALSE;
    for (Int_t istep=0; istep<nStepsDone; istep++) {
      Float_t ctrackPos[3];     
      ctrackPos[0] = fTrackPos[0][istep];
      ctrackPos[1] = fTrackPos[1][istep];
      ctrackPos[2] = fTrackPos[2][istep];

      //now see whether the track matches any of the TOF clusters            

      Float_t dist3d[3]={0.,0.,0.};
      accept = kFALSE;

      for (Int_t i=0; i<nc; i++) {

        AliTOFGeometry::IsInsideThePad((TGeoHMatrix*)(&global[i]),ctrackPos,dist3d);

        // check multiple hit cases
        AliESDTOFCluster *cmatched=(AliESDTOFCluster *) ESDclus->UncheckedAt(clind[i]);

        if(cmatched->GetNTOFhits() > 1){ // correct residual for mean position of the clusters (w.r.t. the first pad/hit)
          Float_t zmain = cmatched->GetTOFchannel(0)/48;
          Float_t xmain = cmatched->GetTOFchannel(0)%48;
          for(Int_t ihit=1;ihit < cmatched->GetNTOFhits();ihit++){
            Float_t deltaz = (cmatched->GetTOFchannel(ihit)/48 - zmain) * 3.5;
            Float_t deltax = (cmatched->GetTOFchannel(ihit)%48 - xmain) * 2.5;
            dist3d[0] -= deltax / cmatched->GetNTOFhits();
            dist3d[2] -= deltaz / cmatched->GetNTOFhits();
          }
        }
          
        // ***** NEW *****
        /* if track is inside this cluster set flags which will then
         * inhibit to add track points for the other clusters */
    Float_t yLoc = dist3d[1];
    Float_t rLoc = TMath::Sqrt(dist3d[0]*dist3d[0]+dist3d[2]*dist3d[2]);
	accept = (TMath::Abs(yLoc)<padDepth*0.5 && rLoc<dCut);

	//***** NEW *****
	/* add point everytime that:
	 * - the tracks is within dCut from the cluster
	 */
        if (accept) {

	  Double_t timesCurrent[AliPID::kSPECIESC];
	  AliDebug(3,Form(" Momentum for track %d -> %f", iseed,t->P()));
	  for (Int_t j=0;j<AliPID::kSPECIESC;j++) {
	    timesCurrent[j] = timesOr[j] + fTimesAr[j][istep];
	  }


	  if (TMath::Abs(dist3d[1])<stepSize && !isClusterMatchable[i]) {
	    isClusterMatchable[i] = kTRUE;
	    
	    Int_t currentpos = tmpESDclus->GetEntriesFast(); // position of cluster in ESD
	    if(fWrittenInPos[clind[i]] != -1){
	      currentpos = fWrittenInPos[clind[i]];
	      cmatched = (AliESDTOFCluster *) tmpESDclus->At(currentpos); // update the new one in the ESDEvent
	    }
	    else{ // add as a new cluster in the ESD TClonesArray
	      AliESDTOFCluster *clnew =  new( (*tmpESDclus)[currentpos] ) AliESDTOFCluster(*cmatched);
	      clnew->SetESDID(currentpos);

	     /* // remap also TOF hit in the filtered array
	      for(Int_t ii=0;ii < cmatched->GetNTOFhits();ii++){
            Int_t index = cmatched->GetHitIndex(ii);
            AliESDTOFHit *hitOld = (AliESDTOFHit *) esdTOFHitArr->At(index);
            Int_t index_new = fHitsESD->GetEntriesFast();
            AliESDTOFHit *hitNew = new( (*fHitsESD)[index_new] ) AliESDTOFHit(*hitOld);
            hitNew->SetESDTOFClusterIndex(currentpos);
            clnew->SetHitIndex(ii,index_new);
	      }*/

	      fWrittenInPos[clind[i]] = currentpos;
	      cmatched = clnew; // update the new one added to the ESDEvent
	    }

	    if(cmatched->GetNMatchableTracks() < AliESDTOFCluster::kMaxMatches){
	      cmatched->Update(t->GetID(),dist3d[0],dist3d[1],dist3d[2],fTrackPos[3][istep],timesCurrent);//x,y,z -> tracking RF
	      t->AddTOFcluster(currentpos);
	      t->SetStatus(AliESDtrack::kTOFout);
	    }
	  }
          AliDebug(2,Form(" dist3dLoc[0] = %f, dist3dLoc[1] = %f, dist3dLoc[2] = %f ",dist3d[0],dist3d[1],dist3d[2]));

          nfound++;
       // ***** NEW *****
        }//end if accept
        
      } //end for on the clusters
    } //end for on the steps     


    if (nfound == 0 ) {
      AliDebug(1,Form(" No matchable track points for the track number %d",iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" Number of track points for the track number %d: %d",iseed,nfound));

    Int_t nMatchedClusters = t->GetNTOFclusters();
 
    if (nMatchedClusters==0) {
      AliDebug(1,Form("Reconstructed track %d doesn't match any TOF cluster", iseed));
      fnunmatch++;
      continue;
    }

    AliDebug(1,Form(" %d - matched (%d)",iseed,nMatchedClusters));

    fnmatch++;

    /*
    AliTOFcluster cTOF = AliTOFcluster(volIdClus,
    (Float_t)posClus[0],(Float_t)posClus[1],(Float_t)posClus[2],
    (Float_t)covClus[0],(Float_t)covClus[1],(Float_t)covClus[2],
    (Float_t)covClus[3],(Float_t)covClus[4],(Float_t)covClus[5],
    tofLabels,volIndices,parClu,kTRUE,index[i]);

    // Fill the track residual histograms.
    FillResiduals(trackTOFin,c,kFALSE);
    */
  } // loop on fSeeds
}
//_________________________________________________________________________
Int_t AliAnalysisTaskPJ_B::FindClusterIndex(Double_t z, TClonesArray* ESDclus) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  TClonesArray* TOFClArr = ESDclus; // use temporary array
  Int_t n = TOFClArr->GetEntriesFast();
    if (TOFClArr){
  if (n==0) return 0;
  if (z <= ((AliTOFcluster *) TOFClArr->UncheckedAt(0))->GetZ()) return 0;
  if (z > ((AliTOFcluster *) TOFClArr->UncheckedAt(n-1))->GetZ()) return n;
  Int_t b=0, e=n-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > ((AliTOFcluster *) TOFClArr->UncheckedAt(m))->GetZ()) b=m+1;

    else e=m; 
  }
  return m;
    }
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
