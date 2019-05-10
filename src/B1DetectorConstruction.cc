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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4UnionSolid.hh"
#include "G4VSolid.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fCopper(0),
  fWorld(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 150.*cm, env_sizeZ = 300.*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_AIR");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
  /*
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 
  //     
  // Shape 1: Xenon Chamber
  //
  // Making Xe 136 isotope
  G4double atomic_weight = 131.293*g/mole;
  G4double density = 0.08856*g/cm3;
  G4Isotope* isoXe136 = new G4Isotope("Xe136",
				      54, // Z
				      131, // A
				      atomic_weight);
  
  G4Element* enrichedXe = new G4Element("enrichedXe", "Xe", 1);
  enrichedXe->AddIsotope(isoXe136, 100*perCent);

  // Now make chamber
  G4Material* shape1_mat = new G4Material("shape1_mat", density, 1);
  shape1_mat->AddElement(enrichedXe, 1);
  
  //G4Material* shape1_mat = nist->FindOrBuildMaterial("G4_Xe");

  G4ThreeVector pos1 = G4ThreeVector(0, 0, 0); // at center of envelope

  G4double shape1_pRmin = 0.*cm, shape1_pRmax = 68.*cm;
  G4double shape1_pSPhi = 0.*deg, shape1_pDPhi = 360.*deg;
  G4double shape1_pDz = 52.5*cm;
  G4Tubs* solidShape1 =
    new G4Tubs("Shape1",
  		 shape1_pRmin, shape1_pRmax, shape1_pDz, shape1_pSPhi,
 		 shape1_pDPhi);

                      
  G4LogicalVolume* logicShape1 =                         
    new G4LogicalVolume(solidShape1,         //its solid
                        shape1_mat,          //its material
                        "Shape1");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos1,                    //at position
                    logicShape1,             //its logical volume
                    "Shape1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  */
  //     
  // Shape 3: forward end cap
  //


  // first make base
  G4double base_pRmin = 0., base_pRmax = 68.*cm;//shape1_pRmax;
  G4double base_pSPhi = 0*deg, base_pDPhi = 360*deg;
  G4double base_pDz = 3.*cm;
  
  G4ThreeVector posbase = G4ThreeVector(0, 0, 0);//-(shape1_pDz+base_pDz*2)); 
     
  G4Tubs* solidbase = 
    new G4Tubs("Base",                      //its name
	       base_pRmin, base_pRmax,
	       base_pDz, base_pSPhi,
	       base_pDPhi);   
  
  // then make hats to subtract to make base with holes
  
  G4double hole_pRmin = 0., hole_pRmax = 17.62*cm;
  G4double hole_pSPhi = 0*deg, hole_pDPhi = 360*deg;
  G4double hole_pDz = 6.*cm;

  G4ThreeVector poshole = G4ThreeVector(0,0,0);
  
  G4Tubs* solidhole = 
    new G4Tubs("Hole",                      //its name
	       hole_pRmin, hole_pRmax,
	       hole_pDz, hole_pSPhi,
	       hole_pDPhi);

  G4SubtractionSolid* subtract = new G4SubtractionSolid("Base-Hole", solidbase, solidhole);
  
  // then make cans to cover holes
  /*G4double can_pRmin = 0., can_pRmax = shape1_pRmax;
  G4double can_pSPhi = 0*deg, can_pDPhi = 360*deg;
  G4double can_pDz = 3.*cm;
  
  G4Material* can_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4ThreeVector poscan = G4ThreeVector(0, 0, -(shape1_pDz+base_pDz*2)); // FIXME
     
  G4Tubs* solidcan = 
    new G4Tubs("Can",                      //its name
	       base_pRmin, base_pRmax,
	       base_pDz, base_pSPhi,
	       base_pDPhi);  
  */

  // then combine
  G4double endcap_pDz = 6.*cm;
  G4ThreeVector posEndCap = G4ThreeVector(0, 0, 0); //-(shape1_pDz+endcap_pDz));
  G4Material* endcap_mat = nist->FindOrBuildMaterial("G4_Cu");
  
  //G4VSolid* CopperEndCap = new G4UnionSolid("CopperEndCap", solidbase, solidcan);

  G4LogicalVolume* logicCopperEndCap =
    new G4LogicalVolume(subtract,
			endcap_mat,
			"logicCopperEndCap");
  new G4PVPlacement(0,
		    posEndCap,
		    logicCopperEndCap,
		    "ShapeCopperEndCap",
		    logicWorld,
		    false,
		    0,
		    checkOverlaps);
  /*
  //     
  // Shape 4: (not forward?) endcap
  //
  G4double shape4_pRmin = 0., shape4_pRmax = shape1_pRmax;
  G4double shape4_pSPhi = 0*deg, shape4_pDPhi = 360*deg;
  G4double shape4_pDz = 6.*cm;
  
  G4Material* shape4_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4ThreeVector pos4 = G4ThreeVector(0, 0, shape1_pDz+shape4_pDz);
     
  G4Tubs* solidShape4 = 
    new G4Tubs("Shape4",                      //its name
	       shape4_pRmin, shape4_pRmax,
	       shape4_pDz, shape4_pSPhi,
	       shape4_pDPhi);
                
  G4LogicalVolume* logicShape4 =                         
    new G4LogicalVolume(solidShape4,         //its solid
                        shape4_mat,          //its material
                        "Shape4");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicShape4,             //its logical volume
                    "Shape4",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




  //
  // Plastic barrel
  //
  G4double shape6_pRmin = shape1_pRmax, shape6_pRmax = shape1_pRmax+2.5*cm;
  G4double shape6_pDz = shape1_pDz;
  G4double shape6_pSPhi = 0*deg, shape6_pDPhi = 360*deg;
  
  G4Material* shape6_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos6 = G4ThreeVector(0,0,0);

  G4Tubs* solidShape6 = new G4Tubs("Shape6", shape6_pRmin, shape6_pRmax, shape6_pDz, shape6_pSPhi, shape6_pDPhi);

  G4LogicalVolume* logicShape6 =                         
    new G4LogicalVolume(solidShape6,         //its solid
                        shape6_mat,          //its material
                        "Shape6");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    logicShape6,             //its logical volume
                    "Shape6",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Copper barrel
  //
  G4double shape5_pRmin = shape6_pRmax, shape5_pRmax = shape6_pRmax+2.5*cm;
  G4double shape5_pDz = shape1_pDz;
  G4double shape5_pSPhi = 0*deg, shape5_pDPhi = 360*deg;
  
  G4Material* shape5_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4ThreeVector pos5 = G4ThreeVector(0,0,0);

  G4Tubs* solidShape5 = new G4Tubs("Shape5", shape5_pRmin, shape5_pRmax, shape5_pDz, shape5_pSPhi, shape5_pDPhi);

  G4LogicalVolume* logicShape5 =                         
    new G4LogicalVolume(solidShape5,         //its solid
                        shape5_mat,          //its material
                        "Shape5");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logicShape5,             //its logical volume
                    "Shape5",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  
  
  // Color shapes for easier viewing
  //G4VisAttributes *Red = new G4VisAttributes( G4Colour(255/255., 0/255., 0/255.));
  G4VisAttributes *Yellow = new G4VisAttributes(G4Colour(255/255., 255/255., 0/255., .98));
  G4VisAttributes *LightBlue = new G4VisAttributes(G4Colour(0/255.,204/255.,204/255., .1));
  G4VisAttributes *Grey = new G4VisAttributes(G4Colour(153/255., 153./255, 153/255., .1));
  
  logicShape1->SetVisAttributes(LightBlue);
  //logicShape2->SetVisAttributes(Yellow);
  //logicShape3->SetVisAttributes(Yellow);
  logicCopperEndCap->SetVisAttributes(Yellow);
  logicShape4->SetVisAttributes(Grey);
  logicShape5->SetVisAttributes(Grey);
  logicShape6->SetVisAttributes(Grey);
  */
  // Set Shape1 as scoring volume
  //
  fScoringVolume = logicCopperEndCap; //logicShape1;
  fCopper = logicCopperEndCap;
  fWorld = logicWorld;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
