AliAnalysisTask* run1()
{
  gSystem->Load("libTree");
  gSystem->Load("libVMC");
  gSystem->Load("libGeom");
  gSystem->Load("libGui");
  gSystem->Load("libXMLParser");
  gSystem->Load("libMinuit");
  gSystem->Load("libMinuit2");
  gSystem->Load("libProof");
  gSystem->Load("libPhysics");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libCDB");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libSTEER");
  gSystem->Load("libEVGEN");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libCORRFW");
  gSystem->Load("libTOFbase");
  gSystem->Load("libTOFrec");
  gSystem->Load("libRAWDatabase");
  gSystem->Load("libRAWDatarec");
  gSystem->Load("libTPCbase");
  gSystem->Load("libTPCrec");
  gSystem->Load("libITSbase");
  gSystem->Load("libITSrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libTender");
  gSystem->Load("libSTAT");
  gSystem->Load("libTRDrec");
  gSystem->Load("libHMPIDbase");
  gSystem->Load("libPWGPP");
  gSystem->Load("libPWGHFbase");
  gSystem->Load("libPWGDQdielectron");
  gSystem->Load("libPWGHFhfe");
  gSystem->Load("libEMCALUtils");
  gSystem->Load("libPHOSUtils");
  gSystem->Load("libPWGCaloTrackCorrBase");
  gSystem->Load("libEMCALraw");
  gSystem->Load("libEMCALbase");
  gSystem->Load("libEMCALrec");
  gSystem->Load("libTRDbase");
  gSystem->Load("libVZERObase");
  gSystem->Load("libVZEROrec");
  gSystem->Load("libTender");
  gSystem->Load("libTenderSupplies");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libESDfilter");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGCFCorrelationsBase");
  gSystem->Load("libPWGCFCorrelationsDPhi");
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  



  // Add ESD input handler
  AliVEventHandler* rphandler = new AliESDInputHandlerRP();
  
  // Register input handler to manager
  mgr->SetInputEventHandler(rphandler);

  gROOT->LoadMacro("AliAnalysisTaskPJ_B.cxx++g");
  AliAnalysisTask *taskpj = new AliAnalysisTaskPJ_B("OurTask");
    
  gROOT->LoadMacro("./AddTaskEmcalSetup.C");
  AliEmcalSetupTask *setupTask = AddTaskEmcalSetup();
  setupTask->SetGeoPath("$ALICE_PHYSICS/OADB/EMCAL");
  setupTask->SetOcdbPath("");  
    
  //gROOT->LoadMacro("./AddTaskEMCALTenderUsingDatasetDef.C");
  //AliAnalysisTaskSE *emcTender = AddTaskEMCALTenderUsingDatasetDef("LHC12b");
  
  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskEmcalPreparation.C");
  AliAnalysisTaskSE *clusm = AddTaskEmcalPreparation();

  gROOT->LoadMacro("$ALICE_PHYSICS/PWG/EMCAL/macros/AddTaskMatchingChain.C");
  AddTaskMatchingChain("LHC11b",AliVEvent::kAny,"EmcCaloClusters",1.,kTRUE,0.1,kTRUE,kTRUE);

  // Add task(s)
  mgr->AddTask(taskpj);

  gROOT->LoadMacro("./CreateESDChain.C");
  //  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateChain("esdTree","run_152367BF.txt", 100);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutputpj = mgr->CreateContainer("PJ_Bhist", TList::Class(),   
      AliAnalysisManager::kOutputContainer, "PJB.Calib.root");

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
  mgr->StartAnalysis("local", chain, 10000);
}
