#include "TFile.h"
#include "iostream"
#include "fstream"
#include "algorithm"
#include "vector"
    
void evaluation()
{

    THStack * hs = new THStack("hs", "Distance Traveled");
    
    gStyle->SetOptStat(0);    
    //TH3F * hits3d = 
    Int_t n = 5;
    char * vrsnames[] = {"v4", "v5", "v6", "v7", "v8"};
    for(Int_t v=0; v<n; v++){
        TString vrs = vrsnames[v];
        cout<<endl<<string(37,'-')<<" " + vrs + " Geometry Evaluation "<<string(37, '-')<<endl;
        cout<<" "<<string(96,'_')<<" "<<endl;
        cout<<"|"<<string(36,' ')<<" gAlice Header Messages "<<string(36,' ')<<"|"<<endl;
        
        
        TH1F * hpdgcode = (TH1F*) new TH1F(vrs + "hpdgcode",vrs,6001, -3000,3000); 
        hpdgcode->SetXTitle("Pdg Code");
        
        TH1F * hdist = (TH1F*) new TH1F(vrs + "hdist",vrs,1000, 80,100); 
        hdist->SetXTitle("#sqrt{x^{2} + y^{2} + z^{2}} [cm]");
        
        TH2F * hdist2D = (TH2F*) new TH2F(vrs + "hdist2D", vrs, 1000, -20, 20, 1000, -20, 20);
        hdist2D->SetXTitle("x [cm]");
        hdist2D->SetYTitle("y [cm]");
        
        TH2F * htof = (TH2F*) new TH2F(vrs + "htof","tof" + vrs, 54, 0,53, 1000, 0, 30000);
        htof->SetXTitle("MCP");
        htof->SetYTitle("Time");

        TH3F * hits3d = (TH3F*) new TH3F(vrs + "hits3d","XYZ of hits for " + vrs,100, -110, -60, 100, -25, 25, 100, -25, 25);
        hits3d->SetXTitle("  Z(cm)");
        hits3d->SetYTitle("  X(cm)"); //Note, X and Z are flipped to get the beam perspective right
        hits3d->SetZTitle("  Y(cm)");
        
        TH1F * ehist = (TH1F*) new TH1F(vrs + "ehist", "Eloss summed " + vrs + ", E = 0 Excluded", 1000, 0, 0.12);
        ehist->SetXTitle("Energy (GeV)");
        ehist->SetYTitle("Number of Tracks");
        
        TH1F * ehist2 = (TH1F*) new TH1F(vrs + "ehist2", "Etot summed", 1000, 0, 10);
        ehist2->SetXTitle("Energy (GeV)");
        ehist2->SetYTitle("Number of Particles");
        
                
        TH1F * ehist3 = (TH1F*) new TH1F(vrs + "ehist3", "Etot", 1000, 0, 0.000000012);
        ehist3->SetXTitle("Energy (GeV)");
        ehist3->SetYTitle("Number of Particles");
        
        TH1F * ehist4 = (TH1F*) new TH1F(vrs + "ehist4", "Eloss", 1000, 0, 0.05);
        ehist4->SetXTitle("Energy (GeV)");
        ehist4->SetYTitle("Number of Particles");
        
        TH2F * PtEta = (TH2F*) new TH2F(vrs + "PtEta", "Hit Eta vs. Transverse Momentum for " + vrs, 1000, -10, 10, 1000, 0, 10);
        PtEta->SetXTitle("#eta");
        PtEta->SetYTitle("P_{t} #left[#frac{GeV}{c}#right]");
        
        TH2F * pPtEta = (TH2F*) new TH2F(vrs + "pPtEta", "Primary Hit Eta vs. Transverse Momentum for " + vrs, 1000, -10, 10, 1000, 0, 10);
        pPtEta->SetXTitle("#eta");
        pPtEta->SetYTitle("P_{t} #left[#frac{GeV}{c}#right]");
        
        
        TH1F * nPrim = (TH1F*) new TH1F(vrs + "nPrim", vrs + "nPrimaries per event", 500,0,1000);
        nPrim->SetXTitle("Number of Primaries per event");
        nPrim->SetYTitle("Instances");
        
        TH1F * TracksH = (TH1F*) new TH1F(vrs + "TracksH", vrs + "TracksH per event", 500,0,1000);
        TracksH->SetXTitle("Number of TracksH per event");
        TracksH->SetYTitle("Instances");
        
        TH2F * kPtEta = (TH2F*) new TH2F(vrs + "kPtEta", "Primary Eta vs. Transverse Momentum for " + vrs, 1000, -10, 10, 1000, 0, 10);
        kPtEta->SetXTitle("#eta");
        kPtEta->SetYTitle("P_{t} #left[#frac{GeV}{c}#right]");
        
        TH2F * effichist = (TH2F*) new TH2F(vrs + "effichist", "& Efficiency vs. number of Primaries for " + vrs, 1000,0,1000, 101, 0, 1);
        effichist->SetXTitle("Number of Primaries");
        effichist->SetYTitle("Efficiency");
        
        TH1F * nptrue = (TH1F*) new TH1F(vrs + "nptrue", "Number of true primaries for " + vrs, 1000, 0, 1000);
        nptrue->SetXTitle("Number of Primaries");
        nptrue->SetYTitle("Number of Instances");
        nptrue->Sumw2();    
        
        TH1F * npmeas = (TH1F*) new TH1F(vrs + "npmeas", "Number of & measured primaries for " + vrs, 1000, 0, 1000);
        npmeas->SetXTitle("Number of Primaries");
        npmeas->SetYTitle("Number of Instances");
        npmeas->Sumw2();   
        
        TH1F * effich = (TH1F*) new TH1F(vrs + "effich", "& Efficiency for " + vrs, 1000, 0, 1000);
        effich->SetXTitle("Number of Primaries");
        effich->SetYTitle("Number of Instances");
        effich->Sumw2(); 
        
        TH1F * nptrue2 = (TH1F*) new TH1F(vrs + "nptrue2", "Number of true primaries for " + vrs, 1000, 0, 1000);
        nptrue2->SetXTitle("Number of Primaries");
        nptrue2->SetYTitle("Number of Instances");
        nptrue2->Sumw2();    
        
        TH1F * npmeas2 = (TH1F*) new TH1F(vrs + "npmeas2", "Number of || measured primaries for " + vrs, 1000, 0, 1000);
        npmeas2->SetXTitle("Number of Primaries");
        npmeas2->SetYTitle("Number of Instances");
        npmeas2->Sumw2();   
        
        TH1F * effich2 = (TH1F*) new TH1F(vrs + "effich2", "|| Efficiency for " + vrs, 1000, 0, 1000);
        effich2->SetXTitle("Number of Primaries");
        effich2->SetYTitle("Number of Instances");
        effich2->Sumw2();
        
        THStack * hs2 = new THStack("hs2", vrs + " nPrimaries vs. TracksH");
        
        THStack * hs3 = new THStack("hs3", vrs + " Measured vs. Expected Primaries");
            
        
        TString galice = "/galice.root";
        TString gpath = vrs + galice;
        AliRunLoader* rl = AliRunLoader::Open(gpath,vrs + AliConfig::GetDefaultEventFolderName(),
        "read");
        if (rl == 0x0)
        {
            cerr<<"Can not open session for file galice.root\n";
            return;
        }
        TParticle *particle;
        rl->LoadgAlice();
        gAlice = rl->GetAliRun();
        AliFIT* fT0  = (AliFIT*)gAlice->GetDetector("FIT");
        rl->LoadHeader();
        rl->LoadKinematics("READ");
        TTree *TK = rl->TreeK(); 
        TTree *te = rl->TreeE();
        /*
        //TK->Scan("fPx:fPy:fPz");
        vector<Double_t>* brPx;
        vector<Double_t>* brPy;
        vector<Double_t>* brPz;
        
        TK->SetBranchAddress("fPx", &brPx);
        TK->SetBranchAddress("fPy", &brPy);
        TK->SetBranchAddress("fPz", &brPz);
        
        Int_t npart = TK->GetEntries();

        for(Int_t i=0; i<npart; i++){

            vector<Double_t>brPx->at(i);
            cout<<kPx<<endl;
        
            
            Float_t kPx= fPx;
            
            Float_t kPy=par->GetPy();
            Float_t kPz=par->GetPz();
            Float_t kP=TMath::Sqrt(kPx*kPx+kPy*kPy+kPz*kPz);
            Float_t kPt=TMath::Sqrt(kPx*kPx+kPy*kPy);
            Float_t kEta=TMath::ATanH(kPz/kP);
            cout<<"Transverse Momentum : "<<kPt<<endl;
            cout<<"Eta : "<<kEta<<endl;
            
        }
        */
        //Initialize Run Loader
        Int_t retval;
        AliLoader* lstart = rl->GetLoader("FITLoader");
        Int_t iNevents= rl->GetNumberOfEvents();
        
        //brKin->Print();
        
        Int_t primaries = 0;
        Int_t tracks = 0;
        Int_t htracks = 0;
        Int_t stracks = 0;
        Int_t sprimaries = 0;
        Int_t exprimc = 0;
        Int_t hatcside = 0;
        Int_t hstotal = 0;
        Int_t exprima = 0;
        Int_t hataside = 0;
        Float_t effic = 0;
        Float_t effic2 = 0;
        Float_t expdet  = 0;
        Float_t expdet2  = 0;
        cout<<"|"<<string(96,'_')<<"|"<<endl<<endl;
        for (Int_t j=0; j<iNevents; j++){
            //Progress Bar
            Int_t plng = 68;                            // Set progress bar length
            Float_t per = (Float_t) (j+1)/iNevents;     // Define percentage
            Int_t prog = (Int_t) (per * plng);          // Declare progress as percent of progress bar length
            Int_t prol = plng - prog;                   // Remaining space after progress in the bar
            cout<<"\r"<<string(4,' ');                  // Set cursor in front of the line
            printf("Event: %5d [",j+1);                 // Printing current event number
            cout<<string(prog, '|')<<string(prol, '-'); // Printing Progress bar
            printf("] %5.1f%%", per*100);               // Printing percent progress
            cout.flush();                               // Flush buffer
            
            
            //Intializing header for event
            //cout<<"-----------------------------------------------"<<endl;
            Int_t phitsc = 0;
            Int_t phitsa = 0;
            Int_t pexpc = 0;
            Int_t pexpa = 0;
            rl->GetEvent(j);
            AliHeader *header=rl->GetHeader();
            //header->Print();
            Int_t nprim = header->GetNprimary();
            
            primaries+=nprim;
            //hNtracks->Fill(nprim);
            //if(nprim>0) nev++;

            lstart->LoadHits("READ");
            TClonesArray *fHits = fT0->Hits ();
            TTree *th = lstart -> TreeH();
            //th->Print();
            brHits = th->GetBranch("FIT");
            const TObjArray *Particles = gAlice->GetMCApp()->Particles();
            //Particles->Print();
            AliStack* stack = rl->Stack();
            stracks+=stack->GetNtrack();
            sprimaries+=stack->GetNprimary();                     //Total Primaries in AliStack
            Int_t curprim = 0;
            for(Int_t i = 0; i<stack->GetNtrack(); i++){
                //cout<<endl<<"-----------------------------------"<<endl;
                //cout<<"Track : "<<i<<endl;
                
                TParticle * par= stack->Particle(i);
                if(par->IsPrimary() && curprim != stack->GetPrimary(i)){
                    //cout<<"par->IsPrimary == True"<<endl;
                    curprim = stack->GetPrimary(i);
                    Float_t kPx =  par->Px();
                    Float_t kPy =  par->Py();
                    Float_t kPz =  par->Pz();
                    Float_t kP  =  TMath::Sqrt(kPx*kPx+kPy*kPy+kPz*kPz);
                    Float_t kPt =  TMath::Sqrt(kPx*kPx+kPy*kPy);
                    Float_t kEta = par->Eta();//TMath::ATanH(kPz/kP);
                    
                    if(kPt > 0 && TMath::Abs(kPz) > 0){
                        
                        if(vrs  == "v4")
                        {
                            
                            if(/*C-side*/ kEta>=-3.345 && kEta<=-2.087)
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprimc++;           //Primaries in v4 AliStack whose Eta falls within the C-Side's range
                                pexpc++;
                            }

                            if(/*A-side*/ kEta>=3.468 && kEta<=5.000 )
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprima++;           //Primaries in v4 AliStack whose Eta falls within the A-Side's range
                                pexpa++;
                            }
                        }
                        if(vrs  == "v8")
                        {
                            
                            if(/*C-side*/ kEta>=-3.367 && kEta<=-2.075 )
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprimc++;           //Primaries in v8 AliStack whose Eta falls within the C-Side's range
                                pexpc++;
                            }
                            
                            if(/*A-side*/ kEta>=3.468 && kEta<=5.000 )
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprima++;
                                pexpa++;
                            }  
                        }
                        else 
                        {
                            
                            if(/*C-side*/ kEta>=-3.317 && kEta<=-2.025 )
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprimc++;      //Primaries in v5, v6, and v7 AliStack whose Eta falls within the C-Side's range
                                pexpc++;
                            }
                            
                            if(/*A-side*/ kEta>=3.468 && kEta<=5.000 )
                            {
                                kPtEta->Fill(kEta, kPt);
                                exprima++;
                                pexpa++;
                            }
                        }
                    } 
                    
                }
                //cout<<"Px : "<<kPx<<endl;
            }
            //cout<<"nprim : "<<nprim<<endl;

                
                
            Int_t ntracks  = (Int_t) th->GetEntries();
            htracks +=ntracks;
            //cout<<"TracksH : "<<ntracks<<endl;
            //cout<<"Event "<<j<<" : tracks "<<ntracks<<endl;
            //cout<<"------------------------------------------------------------------------"<<endl;
            nPrim->Fill(nprim);
            TracksH->Fill(ntracks);
            
            
            
            for (Int_t track=0; track<ntracks;track++) {
                brHits->GetEntry(track);
                Int_t nhits = fHits->GetEntriesFast();

                if(nhits > 0){
                  //
                  //cout<<"Track "<<track<<" : "<<" hits "<<nhits<<endl;
                }
                
                Float_t etrack = 0;
                Double_t epart = 0;
               
                //cout<<"Track : "<<brHits->GetEntry(track)<<endl;
                //cout<<"--------------------------------"<<endl;
                
                Int_t curtrack = 0;
                for (Int_t hit=0;hit<nhits;hit++) {
                    startHit   = (AliFITHits*) fHits->UncheckedAt(hit);
                    if (!startHit) {
                      ::Error("Exec","The unchecked hit doesn't exist");
                      break;
                    }
                    Int_t itrack=startHit->GetTrack();
                    
                    
                    Int_t part = startHit->Particle();
                    
                    Int_t volume = startHit->Volume();

                    Float_t x=startHit->X();
                    Float_t y=startHit->Y();
                    Float_t z=startHit->Z();
                    if(z < 0){hdist2D->Fill(x, y);}  
                    hits3d->Fill(z,x,y);  
                    Float_t r = TMath::Sqrt(x*x + y*y + z*z);
                    Float_t time=startHit->Time();
                    Float_t mcp = startHit->MCP();
                    
                    Float_t ehit = startHit->Eloss();
                    Double_t etot = startHit->Etot();
                    etrack+=ehit;
                    if(z<0){epart+=etot;}
                    htof->Fill(mcp, time);   
                    hdist->Fill(r);
                    //cout<<"Etot : "<<etot<<endl;
                    if(etot>0){ehist3->Fill(etot);} //Cherenkov Photon Energy
                    if(ehit>0){ehist4->Fill(ehit);}

                    //cout<<"iTrack : "<<itrack<<endl;
                    
                    /* 
                    cout<<"-----------------------------------------------"<<endl;
                    cout<<"Part : "<<part<<endl;
                    cout<<"Volume : "<<volume<<endl;
                    cout<<"Hit "<<hit<<endl;
                    cout<<"X : "<<x<<" | Y : "<<y<<" | Z : "<<z<<endl;
                    */
                    if (itrack != curtrack){     //Check Track Uniqueness
                        curtrack = itrack;
                        TParticle* particle = stack->Particle(itrack);

                        if(particle) {
                            Int_t pdg=particle->GetPdgCode();
                            hpdgcode->Fill(pdg);
                            Float_t Vx=particle->Vx();
                            Float_t Vy=particle->Vy();
                            Float_t Vz=particle->Vz();
                            Float_t Vr=TMath::Sqrt(Vx*Vx+Vy*Vy+Vz*Vz);
                            Float_t Px=particle->Px();
                            Float_t Py=particle->Py();
                            Float_t Pz=particle->Pz();
                            Float_t P=TMath::Sqrt(Px*Px+Py*Py+Pz*Pz);
                            Float_t Pt=TMath::Sqrt(Px*Px+Py*Py);
                            Float_t Eta= particle->Eta();//TMath::ATanH(Pz/P);
                            
                            PtEta->Fill(Eta,Pt);
                            if(particle->IsPrimary()){
                                if(vrs != "v8" && r<88){
                                    hstotal++;
                                }
                                else {hstotal++;}
                                
                                if(vrs  == "v4")
                                {
                                    if(/*C-side*/ Eta>=-3.345 && Eta<=-2.087 && (z > -85 && z < -80) && (mcp>23 && mcp<52))
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hatcside++;         //Primaries in v4 AliStack whose Eta falls within the C-Side's range
                                        phitsc++;
                                    }
                                    
                                    if(/*A-side*/ Eta>=3.468 && Eta<=5.000  && (z > 320 && z < 350) && (mcp>-1 && mcp<24))
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hataside++;         //Primaries in v4 AliStack whose Eta falls within the A-Side's range
                                        phitsa++;
                                    }
                                }
                                if(vrs  == "v8")
                                {
                                   
                                    if(/*C-side*/ Eta>=-3.367 && Eta<=-2.07 )
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hatcside++;         //Primaries in v8 AliStack whose Eta falls within the C-Side's range
                                        phitsc++;
                                    }
                                    
                                    if(/*A-side*/ Eta>=3.468 && Eta<=5.000 )
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hataside++;         //Primaries in v8 AliStack whose Eta falls within the A-Side's range
                                        phitsa++;
                                    }  
                                }
                                else 
                                {
                                    
                                    if(/*C-side*/ Eta>=-3.317 && Eta<=-2.025 )
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hatcside++;//Primaries in v5 v6 or v7 AliStack whose Eta falls within the C-Side's range
                                        phitsc++;
                                    }
                                    
                                    if(/*A-side*/  Eta>=3.468 &&  Eta<=5.000 )
                                    {
                                        pPtEta->Fill(Eta, Pt);
                                        hataside++;//Primaries in v5 v6 or v7 AliStack whose Eta falls within the C-Side's range
                                        phitsa++;
                                    }
                                } 
                            }
                            
                            

                            //Float_t E = particle->E();
                            //epart += E;

                            //Float_t E=particle->E();


                            //cout<<"Vx : "<<Vx<<" | Vy : "<<Vy<<" | Vz : "<<Vz<<endl;

                        }
                    }
                    /*
                    cout<<"\r"<<hatcside<<" C-Side Hits Identified";
                    cout.flush();*/
                    //cout<<"-----------------------------------------------"<<endl;
                    
                }
                if(epart>0){tracks++;}
                if(etrack>0){ehist->Fill(etrack);}
                if(epart>0){ehist2->Fill(epart);}
                
                //ehist2->Fill(epart);
                /*
                if(isV > 0){
                    cout<<"has velocity for "<<isV<<" hits"<<endl;
                }
                if(nhits > 0){
                    cout<<"------------------------"<<endl;
                }
                */
                
            }
            Int_t NPrim = stack->GetNprimary();
                
            if(phitsc>0 && phitsa>0){
                npmeas->Fill(NPrim);
                effich->Fill(NPrim);
            }
            
            if(phitsc>0 || phitsa>0){
                npmeas2->Fill(NPrim);
                effich2->Fill(NPrim);
            }
            
            if(pexpc>0 && pexpa>0){
                nptrue->Fill(NPrim);
                effichist->Fill(NPrim, effic);
                expdet += 1;
                if(phitsc>0 && phitsa>0){
                    effic += 1;
                }
            }
            
            if(pexpc>0 || pexpa>0){
                nptrue2->Fill(NPrim);
                expdet2+=1;
                if(phitsc>0 || phitsa>0){
                    effic2 += 1;
                }
            }
            
            
            
           //cout<<"------------------------------------------------------------------------"<<endl<<endl;
        }
        effic  /= (Float_t) expdet;
        effic2 /= (Float_t) expdet2;
        
        effich->Divide(nptrue);
        effich2->Divide(nptrue2);
        
        cout<<endl<<endl<<string(4,' ')<<string(33,'-')<<string(5,' ')<<"End Statistics"<<string(5,' ')<<string(33,'-')<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"header->GetNprimary()                                            :";
        printf("%8d ",primaries) ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"if(epart>0){tracks++;}                                           :";
        printf("%8d ",tracks)    ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"th->GetEntries()                                                 :";
        printf("%8d ",htracks)   ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"stack->GetNtrack()                                               :";
        printf("%8d ",stracks)   ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"stack->GetNprimary()                                             :";
        printf("%8d ",sprimaries);
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Expected Primaries C-Side                                        :";
        printf("%8d ",exprimc)    ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Expected Primaries A-Side                                        :";
        printf("%8d ",exprima)    ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Primary Hits at C-Side                                           :";
        printf("%8d ",hatcside)  ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Primary Hits at A-Side                                           :";
        printf("%8d ",hataside)  ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Total Unique Primary Hits                                        :";
        printf("%8d ",hstotal)   ;
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Percentage of Primary Hits within Eta range                      :";
        printf("%8.3f\%",(hatcside+hataside)/(Float_t) hstotal * 100);
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"Percentage of Primary Hits out of expected primaries in AliStack :";
        printf("%8.3f\%",hstotal/ (Float_t) (exprima+exprimc) * 100);
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"&& Efficiency?                                                   :";
        printf("%8.3f\%",effic * 100);
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(6,' ')<<"|| Efficiency?                                                   :";
        printf("%8.3f\%",effic2 * 100);
        cout<<string(7, ' ')<<'|'<<endl;
        cout<<string(4,' ')<<'|'<<string(88,'_')<<'|'<<endl<<endl;
        
        
        TString outfile = "/someting.root";
        TString opath = vrs+outfile;
        TFile *rfile = TFile::Open(opath, "RECREATE");
        if(!rfile) {return;}
        hpdgcode->Write();
        htof->Write();
        hdist2D->Write();
        ehist->Write();
        ehist2->Write();
        ehist3->Write();
        ehist4->Write();
        hdist->Write();
        hdist->SetLineColor(v+2);
        nPrim->Write();
        effichist->Write();
        nptrue->Write();
        nptrue2->Write();
        npmeas->Write();
        npmeas->SetLineColor(2);
        npmeas2->Write();
        npmeas2->SetLineColor(2);
        effich->Write();
        effich2->Write();
        TracksH->Write();
        TracksH->SetLineColor(2);
        hs->Add(hdist);

        hs2->Add(nPrim, "nostack");
        hs2->Add(TracksH, "nostack");
        hs2->Write();
        hs3->Add(nptrue);
        hs3->Add(npmeas, "nostack");
        hs3->Write();
        
        hits3d->Write();
        PtEta->Write();
        pPtEta->Write();
        kPtEta->Write();
        rfile->Close();
        
        
        ofstream tfile (vrs + "/effic.txt");
        if(tfile.is_open())
        {
            tfile<<"#|| Efficiency"<<endl;
            tfile<<("%9.5f", effic2 * 100)<<endl;
            tfile<<"#&& Efficiency"<<endl;
            tfile<<("%9.5f", effic* 100)<<endl;
            tfile.close();
        }
        else
        {
            cout<<"Unable to open effic.txt"<<endl;
        }
    }


    //auto c1 = new TCanvas();
/*
    for(Int_t v=n-1; v>=0; v--){
	TString vrs = vrsnames[v];
	cout<<"version : "<<vrs<<endl;
	TString outfile = "/someting.root";
        TString opath = vrs+outfile;
        TFile *file = TFile::Open(opath);
        if(!file) {return;}
	cout<<"file : "<<opath<<endl;

	TString vdist = vrs + "hdist";
	cout<<"TH1F : "<< vdist<<endl;
	TH1F * dist = file->Get(vdist);
	dist->SetLineColor(v+1);
	hs->Add(dist);
	file->Close();
    }
*/
    //hs->Print();
    //hs->Draw("nostack");
    //c1->SetGrid(5,5);
    //c1->BuildLegend(0.68,0.775,0.98, 0.95);
    
    //c1->Draw();
}


