#ifndef StTriFolwConstants_h
#define StTriFolwConstants_h
#include "TString.h"
#include "Rtypes.h"
#include <map>

class TriFlow
{
  public:

  // event cut
  static std::map<float,float> mVzMaxMap;
  static Float_t mVrMax; // 2.0
  static Float_t mVzMax; // only for 54GeV   shaowei  add
  static Int_t mMatchedToFMin; // 2

  // track cut
  static Float_t mDcaEPMax[5]; // 0: 200 GeV, 1: 39 GeV, 2: 27 GeV, 3: 14.5 GeV, 4: 62.4 GeV
  static Float_t mDcaTrMax; // 1.0 for pion, kaon, proton flow calculation
  static Float_t mDcaTrMax_phi; // 2.0 for phi flow calculation
  static Int_t mHitsDedxMin; // 15
  static Int_t mHitsFitTPCMin; // 15
  static Int_t mHitsMaxTPCMin; // 0
  static Float_t mHitsRatioTPCMin; // 0.51
  static Float_t mEtaMax; // 1.0
  static Float_t mPrimPtMin[5]; // 0: 0.15(200 GeV), 1-4: 0.2 (BES)
  static Float_t mGlobPtMin; // 0.1
  static Float_t mPrimPtMax; // 2.0
  static Float_t mPrimPtWeight; // 2.0
  static Float_t mPrimMomMax; // 10.0
  static Float_t mMass2Min; // -10.0
  static Double_t MAGFIELDFACTOR; // kilogauss
  static Int_t mTrackMin; // 2
  static Int_t mTrackMin_Full; // 4
  static Float_t mToFYLocalMax; // 1.8
  static Float_t mToFZLocalMax; // 1.8
  static Float_t mNSigmaElectronMax; // 2.5
  static Float_t mNSigmaPionMax; // 2.5
  static Float_t mNSigmaKaonMax; // 2.5
  static Float_t mNSigmaProtonMax; // 2.5
  static Float_t mMassPion; // 0.13957
  static Float_t mMassKaon; // 0.49368
  static Float_t mMassProton; // 0.93827
  static Float_t mMassPhi; // 1.019461
  static Float_t mVzVpdDiffMax; // 3
  static std::map<float,float> mSigScaleMap;

  // used constant
  static Float_t mShiftOrder2[5];
  static Float_t mShiftOrder3[5];
  static std::map<float,float> mExScaleMap;

  static Float_t mEta_Gap[4]; // 0 = 0.05, 1 = 0.1, 2 = 0.2, 3 = 0.5

  static Float_t pt_low[14];
  static Float_t pt_up[14];
  static Float_t pt_low_spectra[32];
  static Float_t pt_up_spectra[32];
  static Float_t pt_low_phi[3];
  static Float_t pt_up_phi[3];

  static Float_t x_low[14];
  static Float_t x_up[14];
  static Float_t y_low[14];
  static Float_t y_up[14];

  static Int_t cent_low[4]; // shaowei
  static Int_t cent_up[4];
  static TString Centrality_01[4];
  static TString Centrality_23[4]; // shaowei

  static Double_t phi_Psi2_low[7];
  static Double_t phi_Psi2_up[7];
  static Double_t phi_Psi3_low[7];
  static Double_t phi_Psi3_up[7];
  static Double_t Psi2_low[3];
  static Double_t Psi2_up[3];
  static Double_t Psi3_low[5];
  static Double_t Psi3_up[5];

  static Int_t pt_total;
  static Int_t pt_total_spectra;
  static Int_t Centrality_total;
  static Int_t Centrality_start;
  static Int_t Centrality_stop;

  static Int_t EtaGap_total;
  static Int_t EtaGap_start;
  static Int_t EtaGap_stop;

  static Int_t Charge_total;
  static Int_t Charge_start;
  static Int_t Charge_stop;

  static TString Energy[5];

  // mix event
  static Int_t Bin_Centrality; // 9
  static Int_t Bin_Centrality_01; // 9
  static Int_t Bin_pT ;
  static Int_t Bin_VertexZ; // 10
  static Int_t Bin_Phi_Psi; // 5
  static Int_t Buffer_depth; // 3
  static Int_t Flag_ME; // 0 for same event, 1 for mix event
  static TString MixEvent[2];


  private:

  static std::map<float,float> createVzMaxMap()
  {
    std::map<float,float> m;

    m[7.7] = 70.0;
    m[11.5] = 50.0;
    m[19.6] = 70.0;
    m[27.0] = 70.0;
    m[39.0] = 40.0;
    m[62.4] = 40.0;
    m[200.0] = 30.0;
    m[14.5] = 70.0;
    m[54] = 40.0;

    return m;
  }

  static std::map<float,float> createSigScaleMap()
  {
    std::map<float,float> m;

    m[7.7] = 1.0;
    m[11.5] = 1.0;
    m[19.6] = 1.0;
    m[27.0] = 1.9;
    m[39.0] = 1.0;
    m[62.4] = 1.0;
    m[200.0] = 1.0;
    m[14.5] = 1.0;

    return m;
  }

  static std::map<float,float> createExScaleMap()
  {
    std::map<float,float> m;

    m[7.7] = 1.0;
    m[11.5] = 1.0;
    m[19.6] = 1.0;
    m[27.0] = 1.0;
    m[39.0] = 1.0;
    m[62.4] = 1.0;
    m[200.0] = 0.8273;
    m[14.5] = 1.0;

    return m;
  }

  ClassDef(TriFlow, 1)
};

#endif
