void plotit(){

  TFile* f=new TFile("AnalysisResults.root","R");
  f->ls();
  TList* ell1=(TList*) f->Get("EmcalClusterMaker_tmpCaloClusters_EmcCaloClusters");
  ell1->Print();
  TH1* zvert=(TH1*) ell1->FindObject("fHistZVertex");
  zvert->Draw();























}
