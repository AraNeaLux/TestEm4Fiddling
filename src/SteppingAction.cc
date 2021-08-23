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
/// \file electromagnetic/TestEm4/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//
//
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "G4SteppingManager.hh"
#include "G4RunManager.hh"

#include "DetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "Analysis.hh"

#include <math.h>
#include <TVector3.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(EventAction* EvAct)
:G4UserSteppingAction(),fEventAction(EvAct),fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
 G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

 if (!fScoringVolume) {
    const DetectorConstruction* detectorConstruction
      = static_cast<const DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

 // get volume of the current step
 G4LogicalVolume* volume
   = aStep->GetPreStepPoint()->GetTouchableHandle()
     ->GetVolume()->GetLogicalVolume();

 // check if we are in scoring volume
 //if (volume != fScoringVolume) return;

 G4int score;
 if (volume == fScoringVolume){
    score = 1;
  }
  else {
    score = 0;
  }

 analysisManager->FillH1(1, score);

 // add energies of steps inside detector
 G4double EdepStep = aStep->GetTotalEnergyDeposit();
 if (EdepStep > 0.) fEventAction->AddEdep(EdepStep);

 // get presteppoint and pos
 G4StepPoint* prestep = aStep->GetPreStepPoint();
 G4ThreeVector presteppos = prestep->GetPosition(); // position

 // get track and particle
 G4Track* track = aStep->GetTrack();
 G4String partname = track->GetParticleDefinition()->GetParticleName();
 G4String typelim = "gamma";
 if (partname != typelim) return;

 // convert to polar coords 
 TVector3 v1(presteppos.x(),presteppos.y(),presteppos.z());
 TVector3 v2(0,0,0);
 TVector3 v3= v1-v2;

 double theta = v3.Theta()*TMath::RadToDeg();
 double phi = v3.Phi()*TMath::RadToDeg();

 
 //#define PI 3.14159265
 //double r = sqrt(pow(presteppos.x()/CLHEP::cm,2)+pow(presteppos.y()/CLHEP::cm,2)+pow(presteppos.z()/CLHEP::cm,2));
 //double theta = acos((presteppos.z()/CLHEP::cm)/r)*180/PI;
 //double phi = atan2((presteppos.y()/CLHEP::cm),(presteppos.x()/CLHEP::cm))*180/PI;

 double x_flat = presteppos.x()/CLHEP::cm;
 double y_flat = presteppos.y()/CLHEP::cm;
 double z_flat = presteppos.z()/CLHEP::cm;

 analysisManager->FillH2(3, x_flat, y_flat);
 analysisManager->FillH2(4, x_flat, z_flat);
 analysisManager->FillH2(5, y_flat, z_flat);

 analysisManager->FillH2(0, theta, phi);

 if (score==1) {
  analysisManager->FillH2(1, theta, phi);
 }
 else if (score==0) {
  analysisManager->FillH2(2, theta, phi);
 }

 //example of saving random number seed of this event, under condition
 //// if (condition) G4RunManager::GetRunManager()->rndmSaveThisEvent();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

