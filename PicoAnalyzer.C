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
#include "TMath.h"

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
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"

// Define global constants
// const Int_t daynumber     = 6;
// const Int_t Ncentralities = 7;
// const Int_t order         = 20;
// const Int_t twoorder      = 2 * order;

//////////////////////////////// Main Function /////////////////////////////////
void PicoAnalyzer(const Char_t *inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_EpdEP",
                      Int_t   inputp1 = 1
                    )
{

  //Int_t EpOrder = inputp1; // EpOrder = 1, 2, 3

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
  // EPD EP finder to get EPD event plane
  TString EpdEpOutputName = "EpdEpCorrectionHistograms_OUTPUT_";
  EpdEpOutputName += outFile;
  EpdEpOutputName += ".root";
  StEpdEpFinder *mEpFinder = new StEpdEpFinder(1,EpdEpOutputName,"/star/u/dchen/GitHub/EpdAna/EpdEpCorrectionHistograms_INPUT.root");
  int format = 2;
  mEpFinder->SetEpdHitFormat(format);    // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
  mEpFinder->SetnMipThreshold(0.3);    // recommended by EPD group
  mEpFinder->SetMaxTileWeight(3.0);     // recommended by EPD group 3.0
  TClonesArray * mEpdHits = new TClonesArray("StPicoEpdHit");
  unsigned int found;
  // Retrieve picoDst TChain*
  TChain *mPicoDst = picoReader->chain();
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  std::cout << "EpdHit Branch returned found= " << found << std::endl;
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);

  // output root files
  TFile *outputFile = new TFile(outFile,"recreate");
  TH2D *hEastRingRawQDotProduct = new TH2D("hEastRingRawQDotProduct","EPD east Raw Q dot product between different rings",120,0.,120.,400,-2.0,2.0);
  hEastRingRawQDotProduct->GetXaxis()->SetTitle("(Ring_{a},Ring_{b})");
  hEastRingRawQDotProduct->GetYaxis()->SetTitle("Dot Product");
  TH2D *hEastRingWeightedQDotProduct = new TH2D("hEastRingWeightedQDotProduct","EPD east Weighted Q dot product between different rings",120,0.,120.,400,-2.0,2.0);
  hEastRingWeightedQDotProduct->GetXaxis()->SetTitle("(Ring_{a},Ring_{b})");
  hEastRingWeightedQDotProduct->GetYaxis()->SetTitle("Dot Product");
  int iBin = 0;
  for(int i = 0;i<15;i++){
    for(int j = i+1; j<16;j++){
      iBin++;
      hEastRingRawQDotProduct->GetXaxis()->SetBinLabel(iBin,Form("(%d, %d)",i+1,j+1));
      hEastRingWeightedQDotProduct->GetXaxis()->SetBinLabel(iBin,Form("(%d, %d)",i+1,j+1));
    }
  }
  TProfile *tpEastRingRawQDotProduct = nullptr;
  TProfile *tpEastRingWeightedQDotProduct = nullptr;

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
    // Retrieve event information
    StPicoEvent *event = dst->event();
    if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
        break;
    }

    TVector3 pVtx     = event->primaryVertex();
    StEpdEpInfo result = mEpFinder->Results(mEpdHits,pVtx,0);  // and now you have all the EP info you could ever want :-)

    TVector2 eastRingRawQ[16];
    TVector2 eastRingWeightedQ[16];
    for(int i=0; i<16; i++)
    {
      eastRingRawQ[i] = (TVector2) result.EastRingRawQ(1,i+1);
      eastRingWeightedQ[i] = (TVector2) result.EastRingPhiWeightedQ(1,i+1);
    }

    int ibin = 0;
    for(int i = 0;i<15;i++){
      TVector2 RawQa = eastRingRawQ[i];
      TVector2 WeightedQa = eastRingWeightedQ[i];
      for(int j = i+1; j<16;j++){
        TVector2 RawQb = eastRingRawQ[j];
        TVector2 WeightedQb = eastRingWeightedQ[j];
        hEastRingRawQDotProduct->Fill(ibin,RawQa.X()*RawQb.X()+RawQa.Y()*RawQb.Y());
        hEastRingWeightedQDotProduct->Fill(ibin,WeightedQa.X()*WeightedQb.X()+WeightedQa.Y()*WeightedQb.Y());
        ibin++;
        if(ibin == 1) std::cout << "first bin dot product = " << (RawQa.X()*RawQb.X()+RawQa.Y()*RawQb.Y()) << std::endl;
      }
    }

  }  // Event Loop
  tpEastRingRawQDotProduct = hEastRingRawQDotProduct->ProfileX();
  tpEastRingWeightedQDotProduct = hEastRingWeightedQDotProduct->ProfileX();
  outputFile->Write();
  mEpFinder->Finish();
}
