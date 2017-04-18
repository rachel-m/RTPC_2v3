// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Main program
// Startup the simulation
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif
#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "SteppingVerbose.hh"

//#include "LHEP_BIC.hh"


int main(int argc,char** argv) {

  G4int IsInteractive = 0;
  if(argc == 1) IsInteractive = 1; //no parameters so interactive
  if(argc > 1){
    if( !strcmp(argv[1],"NULL") ) IsInteractive = 1;
  }
  // choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  //my Verbose output class
  G4VSteppingVerbose::SetInstance(new SteppingVerbose);
  // Construct the default run manager
  G4RunManager * runManager = new G4RunManager;
  // set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction( );
  detector->SetIsInteractive(IsInteractive);
  runManager->SetUserInitialization(detector);
  //use below insted if cannot install physics_list
  runManager->SetUserInitialization(new PhysicsList);
  
#ifdef G4VIS_USE
  // visualization manager
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif
    
  // set user action classes
  PrimaryGeneratorAction* pga=new PrimaryGeneratorAction(detector);
  runManager->SetUserAction(pga);

  RunAction* runaction = new RunAction;  
  runManager->SetUserAction(runaction);
  EventAction* eventaction = new EventAction(runaction);
  eventaction->SetIsInteractive(IsInteractive);
  runManager->SetUserAction(eventaction);
  runManager->SetUserAction(new SteppingAction(detector, eventaction));
  eventaction->SetDet(detector);
  
  G4UImanager * UImanager = G4UImanager::GetUIpointer();  
  if( argc == 1 ){           // interactive
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute macros/vis.mac");     
#endif
    ui->SessionStart();
    delete ui;
#endif
  }
  else           // Batch mode
    { 
#ifdef G4VIS_USE
      visManager->SetVerboseLevel("quiet");
#endif
      //Run set up macro
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UImanager->ApplyCommand(command+fileName);
      //do ROOT run
      if(pga->GetMode()==EPGA_ROOT)
	runManager->BeamOn(pga->GetNEvents());
    }

  // job termination
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}
