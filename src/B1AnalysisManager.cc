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
/// \file analysis/AnaEx01/src/TestSimAnalysisManager.cc
/// \brief Implementation of the TestSimAnalysisManager class
//
//
// $Id: TestSimAnalysisManager.cc 105494 2018-08-23 09:02:56Z $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "B1AnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1AnalysisManager::B1AnalysisManager(int argc, char** argv)
  :fFactoryOn(false)
{
  fEdep = 0.;
  fEdepe = 0.;
  fEdepc = 0.;
  fEdepw = 0.;
  fxinit = -9999.;
  fyinit = -9999.;
  fzinit = -9999.;
  fxfin = -9999.;
  fyfin = -9999.;
  fzfin = -9999.;
  fid = -1;
  fpid = -1;
  fargc = argc;
  fargv = argv;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1AnalysisManager::Book()
{
  // Create or get analysis manager
  // The choice of analysis technology is done via selection of a namespace
  // in TestSimAnalysisManager.hh
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
      
  // Create directories 
  analysisManager->SetHistoDirectoryName("histo");
  //analysisManager->SetNtupleDirectoryName("
  //                                     
  G4String fileName = "B1";

  if (fargc>=2){
   G4String sub = "_";
   fileName.append(sub);

    for (int i=2; i<fargc; i++){
      fileName.append(fargv[i]);
    }
  }

  analysisManager->OpenFile(fileName);
  
  // Create histograms
  analysisManager->CreateH1("Edep","Edep in Xe", 200, 0., 400*MeV);
  analysisManager->CreateH1("Edepe", "Edepe in Xe", 200, 0, 400*MeV);
  analysisManager->CreateH1("Edepc", "Edep in Copper", 200., 0, 400*MeV);
  analysisManager->CreateH1("Edepw", "Edep in World", 200., 0, 400*MeV);
  analysisManager->CreateH1("xinit", "xinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("yinit", "yinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("zinit", "zinit", 150., -150., 150.*cm);
  analysisManager->CreateH1("xfin", "xfin", 150., -150., 150.*cm);
  analysisManager->CreateH1("yfin", "yfin", 150., -150., 150.*cm);
  analysisManager->CreateH1("zfin", "zfin", 150., -150., 150.*cm);

  analysisManager->CreateNtuple("Ntuple0", "init");  // id = 0
  analysisManager->CreateNtupleDColumn("Edep");  // column id = 0
  analysisManager->CreateNtupleDColumn("Edepe"); // column id = 1
  analysisManager->CreateNtupleDColumn("Edepc"); // column id = 2
  analysisManager->CreateNtupleDColumn("Edepw"); // column id = 3
  analysisManager->CreateNtupleDColumn("xinit");  // column id = 4
  analysisManager->CreateNtupleDColumn("yinit");  // column id = 5
  analysisManager->CreateNtupleDColumn("zinit");  // column id = 6
  analysisManager->FinishNtuple();
  
  analysisManager->CreateNtuple("Ntuple1", "fin");  // id = 1
  analysisManager->CreateNtupleDColumn("xfin");  // column id = 0
  analysisManager->CreateNtupleDColumn("yfin");  // column id = 1
  analysisManager->CreateNtupleDColumn("zfin");  // column id = 2
  analysisManager->CreateNtupleDColumn("trackid");
  analysisManager->CreateNtupleDColumn("pid");
  analysisManager->FinishNtuple();
  
  fFactoryOn = true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1AnalysisManager::Save()
{
  if (! fFactoryOn) return;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
   
  G4cout << "\n----> Ntuples are saved\n" << G4endl;
      
  delete G4AnalysisManager::Instance();
  fFactoryOn = false;
}

void B1AnalysisManager::FillHisto(G4int id, G4double weight)
{
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->FillH1(id, weight);
}

void B1AnalysisManager::FillNtuple0()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->FillNtupleDColumn(0, 0, fEdep);
  analysisManager->FillNtupleDColumn(0, 1, fEdepe);
  analysisManager->FillNtupleDColumn(0, 2, fEdepc);
  analysisManager->FillNtupleDColumn(0, 3, fEdepw);
  analysisManager->FillNtupleDColumn(0, 4, fxinit);
  analysisManager->FillNtupleDColumn(0, 5, fyinit);
  analysisManager->FillNtupleDColumn(0, 6, fzinit);
  
  analysisManager->AddNtupleRow(0);
}

void B1AnalysisManager::FillNtuple1()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->FillNtupleDColumn(1, 0, fxfin);
  analysisManager->FillNtupleDColumn(1, 1, fyfin);
  analysisManager->FillNtupleDColumn(1, 2, fzfin);
  analysisManager->FillNtupleDColumn(1, 3, fid);
  analysisManager->FillNtupleDColumn(1, 4, fpid);

  analysisManager->AddNtupleRow(1);
}

void B1AnalysisManager::Reset()
{
  fEdep = 0.;
  fEdepe = 0.;
  fEdepc = 0.;
  fEdepw = 0.;
  fxinit = -9999.;
  fyinit = -9999.;
  fzinit = -9999.;
  fxfin = -9999.;
  fyfin = -9999.;
  fzfin = -9999.;
  fid = -1;
  fpid = -1;
}
