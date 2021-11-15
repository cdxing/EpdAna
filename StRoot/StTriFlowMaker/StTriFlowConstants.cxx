#include "StTriFlowConstants.h"
#include "TMath.h"
#include "StarClassLibrary/SystemOfUnits.h"

ClassImp(TriFlow)

// event cut
std::map<float,float> TriFlow::mVzMaxMap = TriFlow::createVzMaxMap();
Float_t TriFlow::mVrMax = 2.0;
Float_t TriFlow::mVzMax = 70.0;
Int_t TriFlow::mMatchedToFMin = 2;

// track cut
std::map<float,float> TriFlow::mSigScaleMap = TriFlow::createSigScaleMap();
Float_t TriFlow::mDcaEPMax[5] = {3.0,1.0,1.0,1.0,1.0}; // for event plane reconstruction: 3.0 for 200GeV, 1.0 for BES
Float_t TriFlow::mDcaTrMax = 1.0; // for pion, kaon, proton mDcaTrMax = 1.0 for flow;
Float_t TriFlow::mDcaTrMax_phi = 2.0; // for phi meson mDcaTrMax = 2.0 to fill a tree and apply an additional cut
Int_t TriFlow::mHitsDedxMin = 5;
Int_t TriFlow::mHitsFitTPCMin = 15;
Int_t TriFlow::mHitsMaxTPCMin = 0;
Float_t TriFlow::mHitsRatioTPCMin = 0.51;
Float_t TriFlow::mEtaMax = 1.0;
Float_t TriFlow::mPrimPtMin[5] = {0.15,0.2,0.2,0.2,0.2}; // for event plane reconstruction and for pion, kaon, proton: 0.15 for 200 GeV, 0.2 for BES
Float_t TriFlow::mGlobPtMin = 0.1; // for phi, Lambda, K0s
Float_t TriFlow::mPrimPtMax = 2.0;
Float_t TriFlow::mPrimPtWeight = 2.0;
Float_t TriFlow::mPrimMomMax = 10.0; // also use for gMom
Float_t TriFlow::mMass2Min = -10.0;
Double_t TriFlow::MAGFIELDFACTOR = kilogauss;
Int_t TriFlow::mTrackMin = 2;
Int_t TriFlow::mTrackMin_Full = 4;
Float_t TriFlow::mToFYLocalMax = 1.8;
Float_t TriFlow::mToFZLocalMax = 1.8;
Float_t TriFlow::mNSigmaElectronMax = 2.5;
Float_t TriFlow::mNSigmaPionMax = 2.5;
Float_t TriFlow::mNSigmaKaonMax = 2.5;
Float_t TriFlow::mNSigmaProtonMax = 2.5;
Float_t TriFlow::mMassPion = 0.13957;
Float_t TriFlow::mMassKaon = 0.49368;
Float_t TriFlow::mMassProton = 0.93827;
Float_t TriFlow::mMassPhi = 1.019461;
Float_t TriFlow::mVzVpdDiffMax = 3.0;
Float_t TriFlow::mEta_Gap[4] = {0.05,0.10,0.20,0.50};

// used constant
Float_t TriFlow::mShiftOrder2[5] = {2.0, 4.0, 6.0, 8.0, 10.0};
Float_t TriFlow::mShiftOrder3[5] = {3.0, 6.0, 9.0, 12.0, 15.0};
std::map<float,float> TriFlow::mExScaleMap = TriFlow::createExScaleMap();

//   change pT bin from 11 to 14  line 52,53,58,59,61,62,81 changed
// pt bin
//                              0 , 1 , 2 , 3 , 4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14
Float_t TriFlow::pt_low[14] = {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8};
Float_t TriFlow::pt_up[14]  = {0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0};
Float_t TriFlow::pt_low_spectra[32] = {0.15,0.25,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.8};
Float_t TriFlow::pt_up_spectra[32]  = {0.25,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.2,3.4,3.8,4.2};
Float_t TriFlow::pt_low_phi[3] = {0.2,1.0,1.4};
Float_t TriFlow::pt_up_phi[3]  = {1.0,1.4,3.4};
// x and y range
//                              0 ,  1 ,  2 ,  3 ,  4 ,  5 ,  6 ,  7 ,  8 ,  9 ,  10
Float_t TriFlow::x_low[14] = {-0.4,-0.6,-0.6,-0.6,-0.6,-0.6,-0.6,-0.8,-0.8,-1.0,-1.0,-1.0,-1.0,-1.0};
Float_t TriFlow::x_up[14]  = { 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.6, 1.6, 2.0, 2.0, 2.0, 2.0, 2.4};

Float_t TriFlow::y_low[14] = {-0.15,-0.2,-0.3,-0.45,-0.6,-0.8,-1.0,-1.3,-1.5,-1.7,-1.8,-1.8,-1.8,-1.8};
Float_t TriFlow::y_up[14]  = { 0.15, 0.2, 0.2, 0.20, 0.3, 0.4, 0.6, 1.0, 1.0, 1.0, 1.4, 1.4, 1.4, 1.4};

// Centrality bin
Int_t TriFlow::cent_low[4] = {0,7,4,0}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
Int_t TriFlow::cent_up[4]  = {8,8,6,3}; // 0 = 0-80%, 1 = 0-10%, 2 = 10-40%, 3 = 40-80%
TString TriFlow::Centrality_01[4] = {"0080","0010","1040","4080"};
TString TriFlow::Centrality_23[4] = {"0070","0010","1040","4070"};
/*
Int_t TriFlow::cent_low[3] = {0,5,0}; // 0 = 0-80%, 1 = 0-30%, 2 = 30-80%
Int_t TriFlow::cent_up[3]  = {8,8,4}; // 0 = 0-80%, 1 = 0-30%, 2 = 30-80%
TString TriFlow::Centrality_01[3] = {"0080","0030","3080"};
TString TriFlow::Centrality_23[3] = {"0010","1040","4070"};
*/


// phi-Psi bin
Double_t TriFlow::phi_Psi2_low[7] = {0.0,TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0};
Double_t TriFlow::phi_Psi2_up[7]  = {TMath::Pi()/14.0,2.0*TMath::Pi()/14.0,3.0*TMath::Pi()/14.0,4.0*TMath::Pi()/14.0,5.0*TMath::Pi()/14.0,6.0*TMath::Pi()/14.0,7.0*TMath::Pi()/14.0};
Double_t TriFlow::phi_Psi3_low[7] = {0.0,TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0};
Double_t TriFlow::phi_Psi3_up[7]  = {TMath::Pi()/21.0,2.0*TMath::Pi()/21.0,3.0*TMath::Pi()/21.0,4.0*TMath::Pi()/21.0,5.0*TMath::Pi()/21.0,6.0*TMath::Pi()/21.0,7.0*TMath::Pi()/21.0};

Double_t TriFlow::Psi2_low[3] = {-3.0*TMath::Pi()/2.0,-1.0*TMath::Pi()/2.0,1.0*TMath::Pi()/2.0};
Double_t TriFlow::Psi2_up[3]  = {-1.0*TMath::Pi()/2.0, 1.0*TMath::Pi()/2.0,3.0*TMath::Pi()/2.0};
Double_t TriFlow::Psi3_low[5] = {-4.0*TMath::Pi()/3.0,-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0,1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0};
Double_t TriFlow::Psi3_up[5]  = {-3.0*TMath::Pi()/3.0,-1.0*TMath::Pi()/3.0, 1.0*TMath::Pi()/3.0,3.0*TMath::Pi()/3.0,4.0*TMath::Pi()/3.0};

Int_t TriFlow::pt_total = 14;
Int_t TriFlow::pt_total_spectra = 32;
Int_t TriFlow::Centrality_total = 3;
Int_t TriFlow::Centrality_start = 0;
Int_t TriFlow::Centrality_stop  = 3;

Int_t TriFlow::EtaGap_total = 4;
Int_t TriFlow::EtaGap_start = 0;
Int_t TriFlow::EtaGap_stop  = 1;

Int_t TriFlow::Charge_total = 2;  // shaowei
Int_t TriFlow::Charge_start = 0;
Int_t TriFlow::Charge_stop  = 2;  // shaowei


TString TriFlow::Energy[5] = {"200GeV","54GeV","39GeV","7GeV","15GeV"};

// mix event
Int_t TriFlow::Bin_Centrality = 9;
Int_t TriFlow::Bin_Centrality_01 = 4;
Int_t TriFlow::Bin_pT = 3;
Int_t TriFlow::Bin_VertexZ = 10;
Int_t TriFlow::Bin_Phi_Psi = 5;
Int_t TriFlow::Buffer_depth = 5;
Int_t TriFlow::Flag_ME = 0; // 0 for same event, 1 for mixed event
TString TriFlow::MixEvent[2] = {"SE","ME"};
