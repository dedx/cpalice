void PhaseSpace() {
// example of use of TGenPhaseSpace  
//Author: Valerio Filippini

   if (!gROOT->GetClass("TGenPhaseSpace")) gSystem.Load("libPhysics");

   TLorentzVector target(0.0, 0.0, 0.0, 0.938);
   TLorentzVector beam(0.0, 0.0, .65, .65);
   TLorentzVector W = beam + target;

   //(Momentum, Energy units are Gev/C, GeV)
   Double_t masses[2] = { 0.511*1E-3, 0.511*1E-3} ;

   TGenPhaseSpace event;
   event.SetDecay(W, 2, masses);

   TH2F *h2 = new TH2F("h2","h2", 100, -5, 5, 100, -5, 5);

   for (Int_t n=0;n<100000;n++) {
      Double_t weight = event.Generate();

      TLorentzVector *plepton = event.GetDecay(0);

      TLorentzVector *plepton2    = event.GetDecay(1);
      h2->Fill(plepton.PseudoRapidity(), plepton.PseudoRapidity());
   }
   h2->Draw();
}
