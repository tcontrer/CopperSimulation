//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1RunAction.cc
/// \brief Implementation of the B1RunAction class

#include "B1RunAction.hh"
#include "B1PrimaryGeneratorAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1AnalysisManager.hh"
//#include "B1Analysis.hh"
// #include "B1Run.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include <string.h>
#include <iostream>
using namespace std;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::B1RunAction(B1AnalysisManager* ana)
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  fEdepe(0.),
  fEdepc(0.),
  fEdepw(0.),
  //fargc(argc),
  //fargv(argv),
  fAnalysisManager(ana)
{ 
  // add new units for dose
  // 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray); 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2);
  accumulableManager->RegisterAccumulable(fEdepe);
  accumulableManager->RegisterAccumulable(fEdepc);
  accumulableManager->RegisterAccumulable(fEdepw);

  // Create analysis manager                                                                    
  // The choice of analysis technology is done via selectin of a namespace                      
  // in B4Analysis.hh                                                                           
  //auto analysisManager = G4AnalysisManager::Instance();
  //G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories                                                                         
  //analysisManager->SetHistoDirectoryName("histograms");                                       
  //analysisManager->SetNtupleDirectoryName("ntuple");
  /*
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output                                 
  
  // Book histograms, ntuple                                                                    
  //                                                                                            

  // Creating histograms                                                                        
  analysisManager->CreateH1("Edep","Edep in Xe", 200, 0., 200*MeV);
  analysisManager->CreateH1("xinit", "xinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("yinit", "yinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("zinit", "zinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("xfin", "xfin", 150., -150., 150.*cm);
  analysisManager->CreateH1("yfin", "yfin", 150., -150., 150.*cm);
  analysisManager->CreateH1("zfin", "zfin", 150., -150., 150.*cm);
  // Creating ntuple                                                                            
  //
  analysisManager->CreateNtuple("B1", "Info");
  analysisManager->CreateNtupleDColumn("Edep");
  analysisManager->CreateNtupleDColumn("xinit");
  analysisManager->CreateNtupleDColumn("yinit");
  analysisManager->CreateNtupleDColumn("xinit");
  analysisManager->CreateNtupleDColumn("xfin");
  analysisManager->CreateNtupleDColumn("yfin");
  analysisManager->CreateNtupleDColumn("zfin");
  analysisManager->FinishNtuple();
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1RunAction::~B1RunAction()
{
  //delete G4AnalysisManager::Instance();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // Get analysis manager                                                                       
  //auto analysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->Book();

  
  // reset accumulables to their initial values
  //G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  //accumulableManager->Reset();

  // Open an output file                                                                        
  //
  /*
  G4String fileName = "B1";
  
  if (fargc>=2){
   G4String sub = "_";  
   fileName.append(sub);
    
    for (int i=2; i<fargc; i++){
      fileName.append(fargv[i]);
    }

  }

  analysisManager->OpenFile(fileName);   
 */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::EndOfRunAction(const G4Run* run)
{

  //auto analysisManager = G4AnalysisManager::Instance();
  
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  
  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const B1DetectorConstruction* detectorConstruction
   = static_cast<const B1DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  G4double mass = detectorConstruction->GetScoringVolume()->GetMass();
  G4double dose = edep/mass;
  G4double rmsDose = rms/mass;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const B1PrimaryGeneratorAction* generatorAction
   = static_cast<const B1PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  if (generatorAction)
  {
    const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
    runCondition += particleGun->GetParticleDefinition()->GetParticleName();
    runCondition += " of ";
    G4double particleEnergy = particleGun->GetParticleEnergy();
    runCondition += G4BestUnit(particleEnergy,"Energy");
  }
        
  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     << " Cumulated dose per run, in scoring volume : " 
     << G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  //fAnalysisManager->FillNtupleEdep();
  fAnalysisManager->Save();
  //analysisManager->Write();
  //analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void B1RunAction::AddEdepe(G4double edepe)
{
  fEdepe += edepe;
}

void B1RunAction::AddEdepc(G4double edepc){
  fEdepc += edepc;
}

void B1RunAction::AddEdepw(G4double edepw){
  fEdepw += edepw;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
