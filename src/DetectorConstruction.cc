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
/// \file electromagnetic/TestEm4/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4NistManager.hh"
#include <TVector3.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
:G4VUserDetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //
  // define a material from its elements.   case 1: chemical molecule
  // 
  G4double a, z;
  G4double density;  
  G4int ncomponents, natoms;
 
  G4Element* C = new G4Element("Carbon"  ,"C" , z= 6., a= 12.01*g/mole);
  G4Element* F = new G4Element("Fluorine","N" , z= 9., a= 18.99*g/mole);
 
  G4Material* C6F6 = 
  new G4Material("FluorCarbonate", density= 1.61*g/cm3, ncomponents=2);
  C6F6->AddElement(C, natoms=6);
  C6F6->AddElement(F, natoms=6);
  
  G4cout << C6F6 << G4endl;

  G4Material* world =
  G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
  
  //     
  // Container
  //
  
  G4double world_x = 0.5*5*m;
  G4double world_y = 0.5*5*m;
  G4double world_z  = 0.5*10*m;
  G4Box* solidWorld =
  new G4Box("world", world_x, world_y, world_z); //its size

  //G4double Rmin=0., Rmax=5*m, deltaZ= 10*m, Phimin=0., deltaPhi=360*degree;

  //G4Tubs*  
  //solidWorld = new G4Tubs("world",                        //its name
  //                 Rmin,Rmax,deltaZ,Phimin,deltaPhi);        //its size

  G4LogicalVolume*                         
  logicWorld = new G4LogicalVolume(solidWorld,                //its solid
                                   world,                //its material
                                   "world");                //its name
  G4VPhysicalVolume*                                   
  physiWorld = new G4PVPlacement(0,                        //no rotation
                                   G4ThreeVector(),        //at (0,0,0)
                                 logicWorld,                //its logical volume
                                 "world",                //its name
                                 0,                        //its mother  volume
                                 false,                        //no boolean operation
                                 0);                        //copy number

  //
  // Detector
  //
  G4Material* det_mat = G4NistManager::Instance()->FindOrBuildMaterial("G4_Ge");

  //G4ThreeVector det_pos = G4ThreeVector(0*cm, 0*cm, 10*cm);

  //TVector3 det_posr(0, 0, 10);
  //det_posr.RotateY(45*TMath::DegToRad());
  //G4ThreeVector det_pos = G4ThreeVector(det_posr.X()*cm, det_posr.Y()*cm, det_posr.Z()*cm);

  // Cubelike Detector
  //
  //G4double det_x = 0.5*5*cm;
  //G4double det_y = 0.5*5*cm;
  //G4double det_z  = 0.5*10*cm;
  //G4Box* det =
  //new G4Box("det", det_x, det_y, det_z); //its size

  G4double det_Rmin=0., det_Rmax=4.5*cm, det_z= 0.5*10*cm, det_Phimin=0., det_deltaPhi=360*degree;

  G4Tubs* det = new G4Tubs("det",                        //its name
                   det_Rmin,det_Rmax,det_z,det_Phimin,det_deltaPhi);        //its size

  G4RotationMatrix rotm  = G4RotationMatrix();
  G4double theta = 45*deg;
  rotm.rotateY(theta);
  G4ThreeVector uz = G4ThreeVector(std::sin(theta),  0., std::cos(theta));
  G4ThreeVector position = (10*cm+det_z)*uz;
  G4Transform3D transform = G4Transform3D(rotm,position);

  G4LogicalVolume* detlogicShape =
    new G4LogicalVolume(det,         //its solid
                        det_mat,          //its material
                        "det");           //its name

  new G4PVPlacement(transform,                    //at rotation, position
                    detlogicShape,             //its logical volume
                    "det",                //its name
                    logicWorld,                //its mother  volume
                    false,                   //no boolean operation
                    0);                       //copy number          //overlaps checking

  // make scoring volume
  fScoringVolume = detlogicShape;

  //
  //always return the physical World
  //  
  return physiWorld;


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
