#include <ctime>
#include "TGrid.h"

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir,const char* gridMode, const char* runNumbers,const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker,Int_t workerTTL, Bool_t isMC);


void runGRIDtest(const char* gridMode = "full",const char* pattern = "*ESDs/spc_calo/", const char* gridDir = "/alice/data/2011/LHC11c", const char* runNumbers = "154808 154796", const char* runPeriod = "LHC11c", const char* uniqueName = "PJB_EM_ana",  const char* additionalCXXs = "AliAnalysisTaskPJ_B.cxx",const char* additionalHs  = "AliAnalysisTaskPJ_B.h",Int_t maxFilesPerWorker = 4, Int_t workerTTL = 7200, Bool_t isMC = kFALSE)
{

  gSystem->SetFPEMask();
  gSystem->Setenv("ETRAIN_ROOT", ".");
  gSystem->Setenv("ETRAIN_PERIOD", runPeriod);

  gSystem->AddIncludePath("-Wno-deprecated");
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");

  LoadLib();


// Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");
  

  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/train/AddESDHandler.C");
  AliESDInputHandler* rphandler = AddESDHandler();


  // Add ESD input handler
  //AliVEventHandler* rphandler = new AliESDInputHandlerRP();
  
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
  

  /*----------------------------------------
  gROOT->LoadMacro("./CreateESDChain.C");
  //  gROOT->LoadMacro("$ALICE_ROOT/PWG0/CreateESDChain.C");
  TChain* chain = CreateChain("esdTree","run_152367BF.txt", 100);
  ------------------------------------------*/
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
  AliAnalysisGrid *plugin = CreateAlienHandler(uniqueName, gridDir, gridMode,runNumbers, pattern, additionalCXXs, additionalHs, maxFilesPerWorker, workerTTL,isMC);
  mgr->SetGridHandler(plugin);

  // start analysis
  cout << "Starting GRID Analysis...";
  mgr->SetDebugLevel(2);
  mgr->StartAnalysis("grid",10);

}
 void LoadLib(){
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
}

AliAnalysisGrid* CreateAlienHandler(const char* uniqueName, const char* gridDir,const char* gridMode, const char* runNumbers, const char* pattern, TString additionalCode, TString additionalHeaders, Int_t maxFilesPerWorker,Int_t workerTTL, Bool_t isMC)
 {
  TDatime currentTime;
  TString tmpName(uniqueName);

  // Only add current date and time when not in terminate mode! In this case the
  // exact name has to be supplied by the user
   if(strcmp(gridMode, "terminate"))
     {
  tmpName += "_";
  tmpName += currentTime.GetDate();
  tmpName += "_";
  tmpName += currentTime.GetTime();
}

 TString tmpAdditionalLibs("");
 tmpAdditionalLibs = Form("libTree.so libVMC.so libGeom.so libGui.so libXMLParser.so libMinuit.so libMinuit2.so libProof.so libPhysics.so libSTEERBase.so libESD.so libAOD.so libOADB.so libANALYSIS.so libCDB.so libRAWDatabase.so libSTEER.so libANALYSISalice.so libCORRFW.so libTOFbase.so libRAWDatabase.so libRAWDatarec.so libTPCbase.so libTPCrec.so libITSbase.so libITSrec.so libTRDbase.so libTender.so libSTAT.so libTRDrec.so libHMPIDbase.so libPWGPP.so libPWGHFbase.so libPWGDQdielectron.so libPWGHFhfe.so libEMCALUtils.so libPHOSUtils.so libPWGCaloTrackCorrBase.so libEMCALraw.so libEMCALbase.so libEMCALrec.so libTRDbase.so libVZERObase.so libVZEROrec.so libTender.so libTenderSupplies.so libESDfilter.so libPWGTools.so libPWGEMCAL.so libPWGGAEMCALTasks.so libPWGCFCorrelationsBase.so libPWGCFCorrelationsDPhi.so %s %s",additionalCode.Data(),additionalHeaders.Data());


 TString macroName("");
 TString execName("");
 TString jdlName("");
 macroName = Form("%s.C", tmpName.Data());
 execName = Form("%s.sh", tmpName.Data());
 jdlName = Form("%s.jdl", tmpName.Data());

 AliAnalysisAlien *plugin = new AliAnalysisAlien();
 plugin->SetOverwriteMode();
 plugin->SetRunMode(gridMode);

 // Here you can set the (Ali)ROOT version you want to use
 plugin->SetAPIVersion("V1.1x");
 plugin->SetROOTVersion("v5-34-26");
 plugin->SetAliROOTVersion("v5-06-16");
 plugin->SetAliPhysicsVersion("vAN-20150506");


 plugin->SetGridDataDir(gridDir); // e.g. "/alice/sim/LHC10a6"
 plugin->SetDataPattern(pattern); //dir structure in run directory
 if (!isMC)
   plugin->SetRunPrefix("000");

 plugin->AddRunList(runNumbers);

 plugin->SetGridWorkingDir(Form("work/%s",tmpName.Data()));
 plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output

 plugin->SetAnalysisSource(additionalCode.Data());
 plugin->SetAdditionalLibs(tmpAdditionalLibs.Data());

 plugin->SetDefaultOutputs(kTRUE);
 //plugin->SetMergeExcludes("");
 plugin->SetAnalysisMacro(macroName.Data());
 plugin->SetSplitMaxInputFileNumber(maxFilesPerWorker);
 plugin->SetExecutable(execName.Data());
 plugin->SetTTL(workerTTL);
 plugin->SetInputFormat("xml-single");
 plugin->SetJDLName(jdlName.Data());
 plugin->SetPrice(1);
 plugin->SetSplitMode("se");
 plugin->SetNtestFiles(1);

 return plugin;
}
