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
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"
#include "B1AnalysisManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4TrackStatus.hh"
#include "G4Track.hh"
#include "G4ThreeVector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction, B1AnalysisManager* ana)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0),
  fCopper(0),
  fEnv(0),
  fnumf(0),
  fAnalysisManager(ana)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) { 
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();   

    fCopper = detectorConstruction->GetCopper();
    fEnv = detectorConstruction->GetEnv();
  }

  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // Get position starting position of step
  //auto analysisManager = G4AnalysisManager::Instance();
  G4Track* track = step->GetTrack();
  //const G4ThreeVector& vtx = track->GetVertexPosition();
  //analysisManager->FillH1(1, vtx.getX());
  
  G4int id = track->GetTrackID();

  // Get final position
  G4TrackStatus status = track->GetTrackStatus();
  if (status != fAlive){
    fnumf += 1;
    G4cout <<"DEAD! "<<fnumf<<"\n";
    
    const G4ThreeVector& pos = track->GetPosition();
    const G4ParticleDefinition* pid = track->GetParticleDefinition();

   
    //fAnalysisManager->FillHisto(4, pos.getX());
    G4double x = pos.getX();
    G4double y = pos.getY();
    G4double z = pos.getZ();
      
    fAnalysisManager->FillHisto(4, x);
    fAnalysisManager->FillHisto(5, y);
    fAnalysisManager->FillHisto(6, z);
    fAnalysisManager->Setxfin(x);
    fAnalysisManager->Setyfin(y);
    fAnalysisManager->Setzfin(z);
    fAnalysisManager->Settrackid(id);
    fAnalysisManager->Setpid(pid->GetPDGEncoding());
    fAnalysisManager->FillNtuple1();
  }
  
  
  // check if we are in scoring volume
  if (volume != fScoringVolume && volume != fCopper && volume != fEnv) return;

  // collect energy deposited Xe in this step
  if (volume == fScoringVolume){
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdep(edepStep);  

    if (id == 11){
      fEventAction->AddEdepe(edepStep);
    }
  }

  if (volume == fCopper){
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdepc(edepStep);
  }

  if (volume == fEnv){
    G4double edepStep = step->GetTotalEnergyDeposit();
    fEventAction->AddEdepw(edepStep);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

