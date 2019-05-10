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
/// \file B4Analysis.hh                                                                         
/// \brief Selection of the analysis technology                                                 

#ifndef B1AnalysisManager_h
#define B1AnalysisManager_h 1

#include "g4root.hh"
//#include "g4cvs.hh"                                                                           
//#include "g4xml.hh"                                                                           
class B1AnalysisManager
{
public:
  B1AnalysisManager(int argc, char** argv);
  ~B1AnalysisManager();

  void Book();
  void Save();
  void Reset();
  void FillHisto(G4int id, G4double weight);
  void FillNtuple0();
  void FillNtuple1();
  
  void SetEdep(G4double edep){
    fEdep = edep;}
  void SetEdepe(G4double edepe){
    fEdepe = edepe;}
  void SetEdepc(G4double edepc){
    fEdepc = edepc;}
  void SetEdepw(G4double edepw){
    fEdepw = edepw;}
  void Setxinit(G4double xinit){
    fxinit = xinit;}
  void Setyinit(G4double yinit){
    fyinit = yinit;}
  void Setzinit(G4double zinit){
    fzinit = zinit;}
  void Setxfin(G4double xfin){
    fxfin = xfin;}
  void Setyfin(G4double yfin){
    fyfin = yfin;}
  void Setzfin(G4double zfin){
    fzfin = zfin;}
  void Settrackid(G4int id){
    fid = id;}
  void Setpid(G4int pid){
    fpid = pid;}
  
private:
  G4bool fFactoryOn;

  G4double fEdep;
  G4double fEdepe; // energy from electrons only
  G4double fEdepc; // energy depositied in copper
  G4double fEdepw; // energy in world

  G4double fxinit;
  G4double fyinit;
  G4double fzinit;

  G4double fxfin;
  G4double fyfin;
  G4double fzfin;

  G4int fid;
  G4int fpid;
  
  int fargc;
  char** fargv;
};

  
#endif




