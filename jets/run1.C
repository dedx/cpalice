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

  // Load common libraries (better too many than too few)                       
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
  //gSystem->Load("libTOFrec");                                                 
  gSystem->Load("libRAWDatabase.so");
  gSystem->Load("libRAWDatarec.so");
  gSystem->Load("libTPCbase.so");
  gSystem->Load("libTPCrec.so");
  gSystem->Load("libITSbase.so");
  gSystem->Load("libITSrec.so");
  gSystem->Load("libTRDbase.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libSTAT.so");
  gSystem->Load("libTRDrec.so");
  gSystem->Load("libHMPIDbase.so");
  gSystem->Load("libPWGPP.so");
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
  gSystem->Load("libTENDER");
  gSystem->Load("libTENDERSupplies");
  gSystem->Load("libPWGEMCAL");
  gSystem->Load("libPWGGAEMCALTasks");
  gSystem->Load("libPWGTools");
  gSystem->Load("libPWGCFCorrelationsBase");
  gSystem->Load("libPWGCFCorrelationsDPhi");

  //load CGAL, Fastjet and SISCone                                              
  gSystem->Load("libCGAL");
  gSystem->Load("libfastjet");
  gSystem->Load("libSISConePlugin");
  gSystem->Load("libJETAN");
  gSystem->Load("libFASTJETAN");
  gSystem->Load("libPWGJEEMCALJetTasks");


  // load analysis framework
  gSystem->Load("libANALYSIS");

  // if running in aliroot, load only this library
  gSystem->Load("libANALYSISalice");

  gROOT->LoadMacro("./CreateESDChain.C");
//  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateESDChain("pass1", 13);

  // for includes use either global setting in $HOME/.rootrc
  // ACLiC.IncludePaths: -I$(ALICE_ROOT)/include
  // or in each macro
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TOF");
  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT -I$ALICE_ROOT/include -I$ALICE_ROOT/EM\
CAL");
  gSystem->AddIncludePath("-I$ALICE_ROOT/PWGDQ/dielectron -I$ALICE_ROOT/PWGHF/h\
fe");
  gSystem->AddIncludePath("-I$ALICE_ROOT/JETAN -I$ALICE_ROOT/JETAN/fastjet");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Add ESD input handler
  AliVEventHandler* rphandler = new AliESDInputHandlerRP();

  // Register input handler to manager
  mgr->SetInputEventHandler(rphandler);

  // Create task

  gROOT->LoadMacro("AliAnalysisTaskFullppJetkdt.cxx+g");
  AliAnalysisTask *taskpt = new AliAnalysisTaskFullppJetkdt("TaskPt");

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
