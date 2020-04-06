/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 *
 * BBC Event Plane Builder
 * \author Ding Chen
 *
 * Updated to Calculate v2
 * \date Dec 23, 2019
 *
 * Tweak to analyze BES-II FXT 3GeV, 7.2GeV and more
 * \author Ding Chen
 * \date Feb 19, 2020
 */

// This is needed for calling standalone classes (not needed on RACF)
//#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <stdio.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
//EPD header files
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"


// Define global constants
// const Int_t daynumber     = 6;
// const Int_t Ncentralities = 7;
// const Int_t order         = 20;
// const Int_t twoorder      = 2 * order;

//////////////////////////////// Main Function /////////////////////////////////
void PicoAnalyzer(const Char_t *inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_EpdEP"//,
                      //Int_t   inputp1 = 1
                    )
{

  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("EpdHit",1);

  if( !picoReader->chain() ) {
      std::cout << "No chain has been found." << std::endl;
  }

  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Number of events to read: " << events2read << std::endl;

  outFile.Append(".picoDst.result.root");
  TFile *outputFile = new TFile(outFile,"recreate");
  TH2D *hEpdRawHitsAll = new TH2D("hEpdRawHits","Tile Center of EPD hits",200,-100.0,100.0,200,-100.0,100.0);
  hEpdRawHitsAll->GetXaxis()->SetTitle("X [cm]");
  hEpdRawHitsAll->GetYaxis()->SetTitle("Y [cm]");
  TH2D *hEpdRawHits[16];
  for(int i=0; i<16; i++)
  {
    hEpdRawHits[i] = new TH2D(Form("hEpdRawHitsRow%d",i+1),Form("Tile Center of EPD hits Row %d",i+1),200,-100.0,100.0,200,-100.0,100.0);
    hEpdRawHits[i]->GetXaxis()->SetTitle("X [cm]");
    hEpdRawHits[i]->GetYaxis()->SetTitle("Y [cm]");
  }


  // Event loop
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
  {
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1)
    << "/" << events2read << "]" << std::endl;

    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..."
        << std::endl;
        break;
    }

    StPicoDst *dst = picoReader->picoDst();
    int nEpdHits = dst->numberOfEpdHits();
    // std::cout<< "# of Epd Hits = " << nEpdHits << std::endl;
    StEpdGeom *mEpdGeom = new StEpdGeom();

    // Loop over EPD hits
    for (int iEpdHit = 0; iEpdHit < nEpdHits; iEpdHit++){
      StPicoEpdHit* epdHit =dst->epdHit(iEpdHit);

      int tileId,ring,TT,PP,EW,ADC;
      float nMip;

      nMip = epdHit->nMIP();
      tileId = epdHit->id();
      EW = (tileId<0)?0:1;
      ring = epdHit->row();
      if(nMip<0.3) continue;// Threshold
      if(EW!=0) continue;//Epd East
      TVector3 tileCenter = mEpdGeom->TileCenter(tileId);

      hEpdRawHitsAll->Fill(tileCenter.X(),tileCenter.Y());
      for(int row=1;row<17;row++)
      {
        if(ring == row) hEpdRawHits[row-1]->Fill(tileCenter.X(),tileCenter.Y());
      }
      // see: https://www.star.bnl.gov/webdata/dox/html/classStPicoEpdHit.html
      } // loop over EPD hits
  }  // Event Loop
  outputFile->Write();
}
