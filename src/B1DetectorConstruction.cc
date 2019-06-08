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
#include <math.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0),
  fCopper(0),
  fEnv(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
     
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 150.*cm;
  G4double world_sizeZ  = 150.*cm;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
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

  //                                                                  
  // Envelope                                                                       
  //                                                                                  
  G4Box* solidEnv =
    new G4Box("Envelope",                    //its name                      
              0.95*world_sizeXY, 0.95*world_sizeXY, 0.95*world_sizeZ); //its size   
  
  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,            //its solid                                        
                        world_mat,             //its material                                 
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
  // Xenon Chamber
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
  G4Material* XeChamber_mat = new G4Material("XeChamber_mat", density, 1);
  XeChamber_mat->AddElement(enrichedXe, 1);
  
  //G4Material* XeChamber_mat = nist->FindOrBuildMaterial("G4_Xe");

  G4ThreeVector posXeChamber = G4ThreeVector(0, 0, 0); // at center of world

  G4double XeChamber_pRmin = 0.*cm, XeChamber_pRmax = 51.5*cm;
  G4double XeChamber_pSPhi = 0.*deg, XeChamber_pDPhi = 360.*deg;
  G4double XeChamber_pDz = 68*cm;
  G4Tubs* solidXeChamber =
    new G4Tubs("XeChamber",
  		 XeChamber_pRmin, XeChamber_pRmax, XeChamber_pDz, XeChamber_pSPhi,
 		 XeChamber_pDPhi);

                      
  G4LogicalVolume* logicXeChamber =                         
    new G4LogicalVolume(solidXeChamber,         //its solid
                        XeChamber_mat,          //its material
                        "XeChamber");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    logicXeChamber,             //its logical volume
                    "XeChamber",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking





  //
  // Plastic barrel
  //
  G4double PlasticBarrel_pRmin = XeChamber_pRmax, PlasticBarrel_pRmax = XeChamber_pRmax+2*cm;
  G4double PlasticBarrel_pDz = XeChamber_pDz;
  G4double PlasticBarrel_pSPhi = 0*deg, PlasticBarrel_pDPhi = 295*deg;
  
  G4Material* PlasticBarrel_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
  G4ThreeVector pos6 = G4ThreeVector(0,0,0);

  G4Tubs* solidPlasticBarrel = new G4Tubs("PlasticBarrel", PlasticBarrel_pRmin, PlasticBarrel_pRmax, PlasticBarrel_pDz, PlasticBarrel_pSPhi, PlasticBarrel_pDPhi);

  G4LogicalVolume* logicPlasticBarrel =                         
    new G4LogicalVolume(solidPlasticBarrel,         //its solid
                        PlasticBarrel_mat,          //its material
                        "PlasticBarrel");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos6,                    //at position
                    logicPlasticBarrel,             //its logical volume
                    "PlasticBarrel",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  //
  // Copper barrel
  //
  G4double CopperBarrel_pRmin = PlasticBarrel_pRmax, CopperBarrel_pRmax = PlasticBarrel_pRmax+12.*cm;
  G4double CopperBarrel_pDz = XeChamber_pDz;
  G4double CopperBarrel_pSPhi = 0*deg, CopperBarrel_pDPhi = 240*deg;

  G4Material* shielding_mat = nist->FindOrBuildMaterial("G4_Cu");
  G4ThreeVector pos5 = G4ThreeVector(0,0,0);

  G4Tubs* solidCopperBarrel = new G4Tubs("CopperBarrel", CopperBarrel_pRmin, CopperBarrel_pRmax, CopperBarrel_pDz, CopperBarrel_pSPhi, CopperBarrel_pDPhi);

  G4LogicalVolume* logicCopperBarrel =                         
    new G4LogicalVolume(solidCopperBarrel,         //its solid
                        shielding_mat,          //its material
                        "CopperBarrel");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos5,                    //at position
                    logicCopperBarrel,             //its logical volume
                    "CopperBarrel",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  
  //     
  // forward end cap
  //
 
  G4double base_pRmax = CopperBarrel_pRmax;
  G4double base_pDz = 6.*cm; //52.5*cm;
  G4double hole_pRmax = (7.9/2.)*cm;

  // Make base to put holes into (holes are for pmts)
  G4VSolid * solidfinal;
  G4Tubs* solidbase = new G4Tubs("CopperBase", 0., base_pRmax, base_pDz, 0., 360.*deg);
  G4Tubs* solidhole = new G4Tubs("Hole", 0., hole_pRmax, base_pDz*2., 0., 360.*deg);
  solidfinal = solidbase->Clone();
  delete solidbase;
  G4VisAttributes* CanColor = new G4VisAttributes(G4Colour(1.0,0.8,0,1));
  CanColor->SetForceSolid(true); // force the object to be visualized with a surface
  CanColor->SetForceAuxEdgeVisible(true); // force auxiliary edges to be shown
  
  // Loop to make hexagonal rings
  // r is defined as length from center of hexagon to a corner
  double nmax = 1.; // number of cans/holes on each side
  double nrings = 4.;
  double rmin = hole_pRmax*2. + 1.*cm;
  double rmax = base_pRmax - hole_pRmax*2.5;
  double rspacing = rmax/nrings; //hole_pRmax*3.5;
  double r = rmin;
  for (int ring=1; ring<=nrings;ring++){

    // Make can that will cover holes
    G4double can_pRmax = hole_pRmax;
    G4double can_pDz = 6.*cm;
    // Cans are either all flat or get shorter (if input argv[4] = "y")
    //if (fargc > 3 && strcmp(fargv[5],"y")==0){
    //  // cans shrink in height as you go out in rings     
    //if (ring == 3){
    //  can_pDz = 5*cm;
    //}
    //if (ring == 4){
    //  can_pDz = 3*cm;
    //}
    //}
    G4Tubs* solidcan = new G4Tubs("solidcan",0, can_pRmax, can_pDz, 0, 360.*deg);
    double dx=0,dy=0,x=0,y=0; 
    
    // loop over each side of the hexagon, each have diff equations to place cans 
    for (int side=1;side<7;side++){
      
      // starting at corner at y=0, then moving counter clockwise around hexagon
      // with x,y initiated as position of corner
      if (side==1){
	x = r;
	y = 0;
	dx = -r/(2*nmax);
	dy = -sqrt(3)*dx;
      }
      if (side==2){
	x = r/2.;
	y = sqrt(3)*r/2.;
	dx = -r/nmax;
	dy = 0;
      }
      if (side==3){
	x = -r/2.;
	y = sqrt(3)*r/2.;
	dx = -r/(2.*nmax);
	dy = sqrt(3)*dx;
      }
      if (side==4){
	x = -r;
	y = 0;
	dx = r/(2.*nmax);
	dy = -sqrt(3)*dx;
      }
      if (side==5){
	x = -r/2.;
	y = -sqrt(3)*r/2.;
	dx = r/nmax;
	dy = 0;
      }
      if (side==6){
	x = r/2.;
	y = -sqrt(3)*r/2.;
	dx = r/(2.*nmax);
	dy = sqrt(3)*dx;
      }

      // Loop over each hole/can to place on this side of this hexagon
      for (int n=1;n<=nmax;n++){
	
	// Position holes in copper	
	G4ThreeVector poshole = G4ThreeVector(x, y, 0);
	solidfinal = new G4SubtractionSolid("solidfinal",solidfinal, solidhole, 0, poshole);

	// Place cans ontop of holes
	// make new can for each position
	G4ThreeVector poscan = G4ThreeVector(x, y, -(XeChamber_pDz+base_pDz*2+can_pDz));
	G4LogicalVolume* logiccan = new G4LogicalVolume(solidcan, shielding_mat, "logiccan");
	logiccan->SetVisAttributes(CanColor);
	new G4PVPlacement(0, poscan, logiccan, "Can", logicEnv, false, 0, true);
	
	// or union with base to make one solid
	//G4ThreeVector poscan = G4ThreeVector(x, y, -(base_pDz+can_pDz));
	//solidfinal = new G4UnionSolid("solidfinal", solidfinal, solidcan, 0, poscan);

	x += dx;
	y += dy;
      }
    }
    nmax += 1.; // each side has one more can/hole when moving out in radius
    r += rspacing;
  }
  
  // Place Holey Base
  G4LogicalVolume* logicCopperEndcap = new G4LogicalVolume(solidfinal, shielding_mat, "logicCopperEndcap");  
  new G4PVPlacement(0, G4ThreeVector(0, 0, -(XeChamber_pDz+base_pDz)), logicCopperEndcap, "CopperEndcap", logicEnv, false, 0, true);
      

  
  //     
  // Back endcap
  //
  G4double BackEndcap_pRmin = 0., BackEndcap_pRmax = CopperBarrel_pRmax;
  G4double BackEndcap_pSPhi = 0*deg, BackEndcap_pDPhi = 360*deg;
  G4double BackEndcap_pDz = 6.*cm;
  G4ThreeVector pos4 = G4ThreeVector(0, 0, XeChamber_pDz+BackEndcap_pDz);
     
  G4Tubs* solidBackEndcap = 
    new G4Tubs("BackEndcap",                      //its name
	       BackEndcap_pRmin, BackEndcap_pRmax,
	       BackEndcap_pDz, BackEndcap_pSPhi,
	       BackEndcap_pDPhi);
                
  G4LogicalVolume* logicBackEndcap =                         
    new G4LogicalVolume(solidBackEndcap,         //its solid
                        shielding_mat,          //its material
                        "BackEndcap");           //its name
               
  new G4PVPlacement(0,                       //no rotation
                    pos4,                    //at position
                    logicBackEndcap,             //its logical volume
                    "BackEndcap",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  
  ///
  // Set vis attributes like color
  ///
  
  // Make everything else invisible
  //logicEnv->SetVisAttributes(G4VisAttributes::Invisible);

  // Make colors
  G4VisAttributes* Blue = new G4VisAttributes(G4Colour(0,0.5,1.,1));
  G4VisAttributes* Yellow = new G4VisAttributes(G4Colour(1.0,1.0,0,1));
  G4VisAttributes* Grey = new G4VisAttributes(G4Colour(.5,.5,.5,1));  
  G4VisAttributes* White = new G4VisAttributes(G4Colour(1.0,1.0,1.0,0));

  
  // set colors
  logicEnv->SetVisAttributes(White);
  logicCopperEndcap->SetVisAttributes(Yellow);
  logicBackEndcap->SetVisAttributes(Yellow);
  logicCopperBarrel->SetVisAttributes(Yellow);
  logicXeChamber->SetVisAttributes(Blue);
  logicPlasticBarrel->SetVisAttributes(Grey);
  //logiccan->SetVisAttributes(tmpVisAtt);

  
  //
  // Set scoring volumes
  //
  fScoringVolume = logicXeChamber;
  fCopper = logicCopperEndcap;
  fEnv = logicEnv;
  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
