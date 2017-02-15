// run.C
//
// Template run macro for AliBasicTask.cxx/.h with example layout of
// physics selections and options, in macro and task.
//
// Author: Arvinder Palaha
//
class AliAnalysisGrid;
#include <TSystem.h>
#include <TROOT.h>
#include <TGrid.h>
/*#include <AliAnalysisGrid.h>
#include <AliAnalysisTaskSE.h>
#include <AliAnalysisManager.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisTaskElectronTagging.h>
#include <AliAnalysisAlien.h>*/


AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset);

//______________________________________________________________________________
void runETag(
         const char* runtype = "grid", // local, proof or grid
         const char *gridmode = "terminate", // Set the run mode (can be "full", "test", "offline", "submit" or "terminate"). Full & Test work for proof
         const bool bMCtruth = 0, // 1 = MCEvent handler is on (MC truth), 0 = MCEvent handler is off (MC reconstructed/real data)
         const bool bMCphyssel = 0, // 1 = looking at MC truth or reconstructed, 0 = looking at real data
         const Long64_t nentries = 2000, // for local and proof mode, ignored in grid mode. Set to 1234567890 for all events.
         const Long64_t firstentry = 0, // for local and proof mode, ignored in grid mode
         const char *proofdataset = "/alice/data/LHC10c_000120821_p1", // path to dataset on proof cluster, for proof analysis
         const char *proofcluster = "alice-caf.cern.ch", // which proof cluster to use in proof mode
         const char *taskname = "example_task" // sets name of grid generated macros
         )
{
    printf("Made it here");
    // check run type
    if(runtype != "local" && runtype != "proof" && runtype != "grid"){
        printf("\n\tIncorrect run option, check first argument of run macro");
        printf("\tint runtype = local, proof or grid\n");
        return;
    }
    printf("%s analysis chosen",runtype);
  
  
  gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
    // load libraries
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
    
    if (!TGeoGlobalMagField::Instance()->GetField()) printf("Field map is not set.");
    TGeoGlobalMagField::Instance()->SetField(ConfigCalibTrain());
    if (!TGeoGlobalMagField::Instance()->GetField()) printf("Field map is not set after.");
    // add aliroot indlude path
    //gROOT->ProcessLine(Form(".include %s/include",gSystem->ExpandPathName("$ALICE_ROOT")));
    //gROOT->SetStyle("Plain");
        
    // analysis manager
    AliAnalysisManager* mgr = new AliAnalysisManager(taskname);
    
    // create the alien handler and attach it to the manager
    AliAnalysisGrid *plugin = CreateAlienHandler(taskname, gridmode, proofcluster, proofdataset); 
    mgr->SetGridHandler(plugin);
    
    AliESDInputHandler* iH = new AliESDInputHandler();
    iH->SetNeedField(kTRUE);
    //AliAODInputHandler* iH = new AliAODInputHandler();
//    iH->SetInactiveBranches("tracks. vertices. v0s. cascades. jets. caloClusters. fmdClusters. pmdClusters. dimuons. AliAODZDC");
//    iH->SetInactiveBranches("*");
//    iH->SetCheckStatistics(kTRUE);
    mgr->SetInputEventHandler(iH);
        
    // === Physics Selection Task ===
    //
    // In SelectCollisionCandidate(), default is kMB, so the task UserExec() 
    // function is only called for these events.
    // Options are:
    //    kMB             Minimum Bias trigger
    //    kMBNoTRD        Minimum bias trigger where the TRD is not read out
    //    kMUON           Muon trigger
    //    kHighMult       High-Multiplicity Trigger
    //    kUserDefined    For manually defined trigger selection
    //
    // Multiple options possible with the standard AND/OR operators && and ||
    // These all have the usual offline SPD or V0 selections performed.
    //
    // With a pointer to the physics selection object using physSelTask->GetPhysicsSelection(),
    // one can manually set the selected and background classes using:
    //    AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL")
    //    AddBGTriggerClass("+CINT1A-ABCE-NOPF-ALL");
    //
    // One can also specify multiple classes at once, or require a class to NOT
    // trigger, for e.g.
    //    AddBGTriggerClass("+CSMBA-ABCE-NOPF-ALL -CSMBB-ABCE-NOPF-ALL");
    //
    // NOTE that manually setting the physics selection overrides the standard
    // selection, so it must be done in completeness.
    //
    // ALTERNATIVELY, one can make the physics selection inside the task
    // UserExec().
    // For this case, comment out the task->SelectCol.... line, 
    // and see AliBasicTask.cxx UserExec() function for details on this.

//    gROOT->LoadMacro("$ALICE_ROOT/OADB/macros/AddTaskPhysicsSelection.C");
//    AliPhysicsSelectionTask *physSelTask = AddTaskPhysicsSelection(bMCphyssel);
//    if(!physSelTask) { Printf("no physSelTask"); return; }
    //AliPhysicsSelection *physSel = physSelTask->GetPhysicsSelection();
    //physSel->AddCollisionTriggerClass("+CINT1B-ABCE-NOPF-ALL");// #3119 #769");
               
    // create task
    gROOT->ProcessLine(".L AliAnalysisTaskElectronTagging.cxx+g");
    gROOT->LoadMacro("AliAnalysisTaskElectronTagging.cxx");
    AliAnalysisTaskSE* task = new AliAnalysisTaskElectronTagging(taskname);
    task->SelectCollisionCandidates(AliVEvent::kMB); // if physics selection performed in UserExec(), this line should be commented
    mgr->AddTask(task);
    
    // set output root file name for different analysis
    TString outfilename = "phistos.root";
  
    // create containers for input/output
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("phistos", TList::Class(), AliAnalysisManager::kOutputContainer, outfilename);
        
    // connect input/output
    mgr->ConnectInput(task, 0, cinput);
    mgr->ConnectOutput(task, 1, coutput1);
        
    // enable debug printouts
    mgr->SetDebugLevel(2);
    mgr->SetNSysInfo(100);
    if (!mgr->InitAnalysis()) return;
    mgr->PrintStatus();
  
    // start analysis
    Printf("Starting Analysis....");
    mgr->StartAnalysis(runtype,nentries,firstentry);
}

//______________________________________________________________________________
AliAnalysisGrid* CreateAlienHandler(const char *taskname, const char *gridmode, const char *proofcluster, const char *proofdataset)
{
    AliAnalysisAlien *plugin = new AliAnalysisAlien();
    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    plugin->SetRunMode(gridmode);

    // Set versions of used packages
    plugin->SetAPIVersion("V1.1x");
    plugin->SetROOTVersion("v5-34-30-alice-12");
    plugin->SetAliROOTVersion("v5-08-01-1");

    // Declare input data to be processed.
//    plugin->SetCheckCopy(kFALSE);

    // Method 1: Create automatically XML collections using alien 'find' command.
    // Define production directory LFN
    plugin->SetGridDataDir("/alice/data/2012/LHC12d");
    // On real reconstructed data:
    // plugin->SetGridDataDir("/alice/data/2009/LHC09d");
    // Set data search pattern
    //plugin->SetDataPattern("*ESDs.root"); // THIS CHOOSES ALL PASSES
    // Data pattern for reconstructed data
    plugin->SetDataPattern("pass2/*ESDs.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    //    plugin->SetDataPattern("ESDs/pass2/AOD038/*AliAOD.root"); // CHECK LATEST PASS OF DATA SET IN ALIENSH
    plugin->SetRunPrefix("000");   // real data
    // ...then add run numbers to be considered
    Int_t runlist[75]={184127, 184132, 184135, 184137, 184138, 184188, 184208, 184209, 184215, 184371, 184673, 184678, 184687, 184784, 184786, 184928, 184933, 184938, 184964, 184967, 184968, 184987, 184988, 184990, 185029, 185031, 185126, 185127, 185134, 185164, 185189, 185196, 185203, 185206, 185208, 185282, 185288, 185289, 185291, 185292, 185299, 185300, 185302, 185303, 185349, 185350, 185351, 185356, 185465, 185574, 185575, 185578, 185580, 185581, 185582, 185583, 185588, 185589, 185695, 185697, 185698, 185699, 185701, 185738, 185764, 185765, 185775, 185776, 185784, 186163, 186164, 186165, 186205, 186319, 186320};  
    for (Int_t ind=0; ind<1; ind++) {
//     plugin->AddRunNumber(183913);
      plugin->AddRunNumber(runlist[ind]);
    }
    //plugin->SetRunRange(114917,115322);
    plugin->SetNrunsPerMaster(10); // 1
    //plugin->SetOutputToRunNo(0);
    plugin->SetOverwriteMode(kTRUE);
    // comment out the next line when using the "terminate" option, unless
    // you want separate merged files for each run
    plugin->SetMergeViaJDL();

    // Method 2: Declare existing data files (raw collections, xml collections, root file)
    // If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
    // XML collections added via this method can be combined with the first method if
    // the content is compatible (using or not tags)
    //   plugin->AddDataFile("tag.xml");
    //   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");

    // Define alien work directory where all files will be copied. Relative to alien $HOME.
    plugin->SetGridWorkingDir(taskname);

    // Declare alien output directory. Relative to working directory.
    plugin->SetGridOutputDir("out"); // In this case will be $HOME/taskname/out

    // Declare the analysis source files names separated by blancs. To be compiled runtime
    // using ACLiC on the worker nodes.
    plugin->SetAnalysisSource("AliAnalysisTaskElectronTagging.cxx");

    // Declare all libraries (other than the default ones for the framework. These will be
    // loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
    plugin->SetAdditionalLibs("AliAnalysisTaskElectronTagging.h AliAnalysisTaskElectronTagging.cxx");

    // Declare the output file names separated by blancs.
    // (can be like: file.root or file.root@ALICE::Niham::File)
    // To only save certain files, use SetDefaultOutputs(kFALSE), and then
    // SetOutputFiles("list.root other.filename") to choose which files to save
    //plugin->SetDefaultOutputs(kTRUE);
    //plugin->SetOutputFiles("phistos.root");

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    plugin->SetAnalysisMacro(Form("%s.C",taskname));

    // Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
    plugin->SetSplitMaxInputFileNumber(100);

    // Optionally modify the executable name (default analysis.sh)
    plugin->SetExecutable(Form("%s.sh",taskname));

    // set number of test files to use in "test" mode
    plugin->SetNtestFiles(10);

    // Optionally resubmit threshold.
    plugin->SetMasterResubmitThreshold(90);

    // Optionally set time to live (default 30000 sec)
    plugin->SetTTL(30000);

    // Optionally set input format (default xml-single)
    plugin->SetInputFormat("xml-single");

    // Optionally modify the name of the generated JDL (default analysis.jdl)
    plugin->SetJDLName(Form("%s.jdl",taskname));

    // Optionally modify job price (default 1)
    plugin->SetPrice(1);      

    // Optionally modify split mode (default 'se')    
    plugin->SetSplitMode("se");
    
    //----------------------------------------------------------
    //---      PROOF MODE SPECIFIC SETTINGS         ------------
    //---------------------------------------------------------- 
    // Proof cluster
    plugin->SetProofCluster(proofcluster);
    // Dataset to be used   
    plugin->SetProofDataSet(proofdataset);
    // May need to reset proof. Supported modes: 0-no reset, 1-soft, 2-hard
    plugin->SetProofReset(0);
    // May limit number of workers
    plugin->SetNproofWorkers(0);
    // May limit the number of workers per slave
    plugin->SetNproofWorkersPerSlave(1);   
    // May use a specific version of root installed in proof
    plugin->SetRootVersionForProof("current");
    // May set the aliroot mode. Check http://aaf.cern.ch/node/83 
    plugin->SetAliRootMode("default"); // Loads AF libs by default
    // May request ClearPackages (individual ClearPackage not supported)
    plugin->SetClearPackages(kFALSE);
    // Plugin test mode works only providing a file containing test file locations, used in "local" mode also
    plugin->SetFileForTestMode("files.txt"); // file should contain path name to a local directory containg *ESDs.root etc
    // Request connection to alien upon connection to grid
    plugin->SetProofConnectGrid(kFALSE);
    // Other PROOF specific parameters
    plugin->SetProofParameter("PROOF_UseMergers","-1");
    printf("Using: PROOF_UseMergers   : %s\n", plugin->GetProofParameter("PROOF_UseMergers"));
    return plugin;
}

void LoadLibs()
{
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

  // include paths
//gSystem->AddIncludePath("-Wno-deprecated");
  //gSystem->AddIncludePath("-I$ALICE_PHYSICS -I$ALICE_PHYSICS/include -I$ALICE_PHYSICS/EMCAL");
  //gSystem->AddIncludePath("-I$ALICE_PHYSICS/PWGDQ/dielectron -I$ALICE_PHYSICS/PWGHF/hfe");
}

AliMagF* ConfigCalibTrain(const char *ocdb="raw://"){

  // OCDB
 
  printf("setting run to");
  //AliCDBManager::Instance()->SetDefaultStorage("raw://");

  AliCDBManager::Instance()->SetDefaultStorage("alien://folder=/alice/data/2012/OCDB/");
  AliCDBManager::Instance()->SetSpecificStorage("GRP/GRP/Data", "alien://folder=/alice/data/2012/OCDB/");
  AliCDBManager::Instance()->SetRun(183913); 

  // geometry
  AliGeomManager::LoadGeometry();
  AliGeomManager::ApplyAlignObjsFromCDB("GRP");



  AliCDBEntry *entry = AliCDBManager::Instance()->Get("GRP/GRP/Data");
  AliGRPObject* grpData = (AliGRPObject*)entry->GetObject();

  Bool_t ok=kTRUE;
  Float_t l3Current = grpData->GetL3Current((AliGRPObject::Stats)0);
  if (l3Current == AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the L3 current !");
    ok = kFALSE;
  }
  
  Char_t l3Polarity = grpData->GetL3Polarity();
  if (l3Polarity == AliGRPObject::GetInvalidChar()) {
    printf("GRP/GRP/Data entry:  missing value for the L3 polarity !");
    ok = kFALSE;
  }
  
  // Dipole
  Float_t diCurrent = grpData->GetDipoleCurrent((AliGRPObject::Stats)0);
  if (diCurrent == AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the dipole current !");
    ok = kFALSE;
  }
  
  Char_t diPolarity = grpData->GetDipolePolarity();
  if (diPolarity == AliGRPObject::GetInvalidChar()) {
    printf("GRP/GRP/Data entry:  missing value for the dipole polarity !");
    ok = kFALSE;
  }

  TString beamType = grpData->GetBeamType();
  if (beamType==AliGRPObject::GetInvalidString()) {
    printf("GRP/GRP/Data entry:  missing value for the beam type ! Using UNKNOWN");
    beamType = "UNKNOWN";
  }

  Float_t beamEnergy = grpData->GetBeamEnergy();
  if (beamEnergy==AliGRPObject::GetInvalidFloat()) {
    printf("GRP/GRP/Data entry:  missing value for the beam energy ! Using 0");
    beamEnergy = 0;
  }

  // read special bits for the polarity convention and map type
  //Int_t  polConvention = grpData->IsPolarityConventionLHC() ? AliMagF::kConvLHC : AliMagF::kConvDCS2008;
  Int_t  polConvention = grpData->IsPolarityConventionLHC() ? 0 : 1;
  Bool_t uniformB = grpData->IsUniformBMap();
  
  if (ok) {
    AliMagF* fld = AliMagF::CreateFieldMap(TMath::Abs(l3Current) * (l3Polarity ? -1:1),
					   TMath::Abs(diCurrent) * (diPolarity ? -1:1),
					   polConvention,uniformB,beamEnergy, beamType.Data());
    if (fld) {
        //TGeoGlobalMagField::Instance()->SetField( fld );
      //TGeoGlobalMagField::Instance()->Lock();
      printf("Running with the B field constructed out of GRP !");
    }
  }
    else{
  printf("Problem with magnetic field setup\n");
    }
    return fld;
}