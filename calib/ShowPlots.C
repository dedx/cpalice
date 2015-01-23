void ShowPlots() {

  TFile *f = new TFile("Pt.ESD.1.root"); 
  TList *myList = (TList*) f->Get("chistpt");

  TH1F* fHistEtaPhiCC[100];
  TH1F* fHistEtaPhiTC[100];

  char nameTC[100];
  char nameCC[100];
  for (int i = 0; i < 100; i++) {
    sprintf(nameCC,"fHistEtaPhiCC%d",i);
    sprintf(nameTC,"fHistEtaPhiTC%d",i);
    
    fHistEtaPhiCC[i] = (TH1F*)myList->FindObject(nameCC);
    fHistEtaPhiCC[i]->SetMarkerStyle(20);
    fHistEtaPhiTC[i] = (TH1F*)myList->FindObject(nameTC);
    fHistEtaPhiTC[i]->SetMarkerStyle(24);
  }

  for (int i = 0; i < 100; i++){
    if(fHistEtaPhiCC[i]->GetEntries() > 0 && fHistEtaPhiTC[i]->GetEntries() > 0) {
      TCanvas *c1 = new TCanvas();
      c1->cd();
      gStyle->SetOptStat(0);
      fHistEtaPhiCC[i]->Draw();
      fHistEtaPhiTC[i]->Draw("same");

      Double_t xl1=0.7, yl1=0.85, xl2=xl1+0.3, yl2=yl1+0.125;
      TLegend *leg = new TLegend(xl1,yl1,xl2,yl2);
      leg->AddEntry(fHistEtaPhiCC[i],"CaloClusters","p");
      leg->AddEntry(fHistEtaPhiTC[i],"TOFClusters","p");
      leg->Draw();

    }
  }

}
