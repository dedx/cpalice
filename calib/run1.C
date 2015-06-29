AliAnalysisTask* run1()
{
  // If running at root prompt, manually load these libraries
  /*
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libMinuit.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so"); 
  gSystem->Load("libTOFrec.so");
  // load analysis framework
  gSystem->Load("libANALYSIS");
  */
  // if running in aliroot, load only this library
  //  gSystem->Load("libANALYSISalice");

  //  gROOT->LoadMacro("$ROOTSYS/macros/CreateESDChain.C");
//  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
//  TChain* chain = CreateChain("esdTree","run_149881NoB.txt",

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  



  // Add ESD input handler
  AliVEventHandler* rphandler = new AliESDInputHandlerRP();
  
  // Register input handler to manager
  mgr->SetInputEventHandler(rphandler);

  // Create task
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPhysicsSelection.C");
  //AliPhysicsSelectionTask *physSelTask = AddTaskEmcalPhysicsSelection(kTRUE, kFALSE, AliVEvent::kAny, 5, 5, 10, kTRUE, -1, -1, -1,-1);


  //gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  //AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kTRUE);

  gROOT->LoadMacro("AliAnalysisTaskPJ.cxx++g");
  AliAnalysisTask *taskpj = new AliAnalysisTaskPJ("OurTask");
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalSetup.C");
  //gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPreparation.C");
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMatchingChain.C");

  //AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
  //setupTask->SetGeoPath("$ALICE_PHYSICS/OADB/EMCAL");
  //AliAnalysisTaskSE *clusm = AddTaskEmcalPreparation("LHC11h","run_149881NoB.txt");
  AddTaskMatchingChain("LHC11h", AliVEvent::kAny,"EmcCaloClusters",1.,kTRUE,0.1,kFALSE,kFALSE);
  gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
  AliCentralitySelectionTask *centralityTask = AddTaskCentrality(kTRUE);


  
  // Add task(s)
  mgr->AddTask(taskpj);

  gROOT->LoadMacro("./CreateESDChain.C");
  //  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateChain("esdTree","run_149881NoB.txt", 10000);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputpj = mgr->CreateContainer("chistPJ", TList::Class(),   
      AliAnalysisManager::kOutputContainer, "PJ.Calib.root");

  // Connect input/output
  mgr->ConnectInput(taskpj, 0, cinput);

  // No need to connect to a common AOD output container if the task does not
  // fill AOD info.
  //mgr->ConnectOutput(task, 0, coutput0);
  mgr->ConnectOutput(taskpj, 1, coutputpj);


  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, 10000000);
}
