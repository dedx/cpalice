
void jetana()
{
/*
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libANALYSIS.so");
*/
  
  // Load libs (more needed when running with root instead of aliroot)


  // AliAnalyTaskJets derives from AnalysisTaskSE:
  gSystem->Load("libANALYSISalice.so");

  // AliAnalyTaskJets is in libJETAN
  // Other user tasks by gROOT->LoadMacro("AliAnalysisTaskXYZ.cxx+g");
  gSystem->Load("libJETAN.so");
  
  
  // Manual chaining 
  // TChain *chain = new TChain("esdTree");
  // chain->Add("/afs/cern.ch/user/k/kleinb/public/tutorial/local/data/AliESDs.root");

  gROOT->LoadMacro("$ALICE_ROOT/PWGUD/macros/CreateESDChain.C");
//  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("esd.txt", 2);

  
  // Create the analysis manager
  AliAnalysisManager *mgr  = new AliAnalysisManager("My Manager", "My Manager");
  
  // Input
  AliESDInputHandler* inpHandler = new AliESDInputHandler();
  mgr->SetInputEventHandler  (inpHandler);
  
  // output 
  AliAODHandler* aodHandler   = new AliAODHandler();
  aodHandler->SetOutputFileName("aod.root");
  mgr->SetOutputEventHandler(aodHandler);
  
  // Jet analysis     
  AliAnalysisTaskJets *jetana = new AliAnalysisTaskJets("JetAnalysis");
  mgr->AddTask(jetana);
  
  //
  // Create containers for input/output
  // Common ESD input for all tasks
  AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
  // Common AOD output
  AliAnalysisDataContainer *coutput1 = mgr->GetCommonOutputContainer();
  // Custom container for output histograms
  AliAnalysisDataContainer *coutput2 = 
     mgr->CreateContainer("histos", TList::Class(),
                           AliAnalysisManager::kOutputContainer, "histos.root");
  // Connect all input/output slots of the task to containers
  mgr->ConnectInput  (jetana,  0, cinput1  );
  mgr->ConnectOutput (jetana,  0, coutput1 );
  mgr->ConnectOutput (jetana,  1, coutput2 );
  
  //
  // Run the analysis
  //    
  mgr->InitAnalysis();
  mgr->PrintStatus();
  mgr->StartAnalysis("local",chain);
}
