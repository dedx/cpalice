void run1()
{
  // If running at root prompt, manually load these libraries

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

  // if running in aliroot, load only this library
    gSystem->Load("libANALYSISalice");

  gROOT->LoadMacro("macros/CreateESDChain.C");
//  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateChain("esdTree","test.txt", 99);

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

  gROOT->LoadMacro("AliAnalysisTaskPt.cxx++g");
  AliAnalysisTask *taskpt = new AliAnalysisTaskPt("TaskPt");

  // Add task(s)
  mgr->AddTask(taskpt);


  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputpt = mgr->CreateContainer("chistpt", TList::Class(),   
      AliAnalysisManager::kOutputContainer, "Pt.ESD.1.root");

  // Connect input/output
  mgr->ConnectInput(taskpt, 0, cinput);

  // No need to connect to a common AOD output container if the task does not
  // fill AOD info.
  //  mgr->ConnectOutput(task, 0, coutput0);
  mgr->ConnectOutput(taskpt, 1, coutputpt);


  // Enable debug printouts
  mgr->SetDebugLevel(2);

  if (!mgr->InitAnalysis()) return;
  mgr->PrintStatus();
  mgr->StartAnalysis("local", chain, 100000);

}
