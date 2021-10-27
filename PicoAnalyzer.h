#ifndef PicoAnalyzer_h
#define PicoAnalyzer_h

// C++ headers
#include <fstream>
#include <iostream>
#include <map>
#include <iterator>
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
#include "TVector2.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom.h"
#include "TRandom3.h"
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
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"

// Define global constants
// const Int_t daynumber     = 6;
const Int_t _Ncentralities = 9; // 9 centrality bins
const Int_t _EpTermsMaxIni = 20; // Shift Order
const Int_t _nEventTypeBins = 5; // 5 etaRange
const Int_t _nEventTypeBins_tpc = 6; // 5 etaRange for TPC
const Double_t _massPion     = 0.13957039;
const Double_t _massKaon     = 0.493677;
const Double_t _massProton   = 0.938272081;
const Double_t _massPhi = 1.019461;
const Double_t _y_mid = -2.02; // mid rapidity
//--------------------------------------------------------------------
// Histogram
// pt bin
//                                       0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17 ,18 ,19 ,10 ,21 ,22
const Float_t pt_low_phi[23] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2};
const Float_t pt_up_phi[23]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.4,5.8,6.2,6.6};
// Centrality bin
Int_t cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString Centrality_01[4] = {"0080","0010","1040","4080"};
TString Centrality_23[4] = {"0060","0010","1040","4060"};

// phi-Psi bin
Double_t phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
Double_t phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};

Double_t Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
Double_t Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};

Int_t pt_total_phi = 23;

Int_t Centrality_total = 4;    // shaowei
Int_t Centrality_start = 0;
Int_t Centrality_stop  = 4;    // shaowei

Int_t Phi_Psi_total = 7;
// flow analysis
// 0 = pt bin
// 1 = centrality: 0 = 0-80%(0-70%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-70%)
// 2 = SE, RT
// 3 = phi - Psi
TH1F *h_mMass2_EP[23][4][2][7]; //reweight/Res2

TH3F *h_pt_y_mass_se;
TH3F *h_pt_y_mass_rt;
// raw pt spectra
// 0 = pt bin
// 1 = centrality: 0 = 0-80%(0-60%), 1 = 0-10%, 2 = 10-40%, 3 = 40-80%(40-60%)
// 2 = SE, RT
TH1F *h_mMass_Spec[23][4][2]; //reweight

// event plane resolution correction
// 0 = centrality
// 1 = SE, RT
TH1F *h_mMass_Yields[9][2]; //reweight
