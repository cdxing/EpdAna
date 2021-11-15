#ifndef StTriFlowV0_h
#define StTriFlowV0_h

#include "StMessMgr.h"
#include "StPhysicalHelixD.hh"
#include "StTriFlowMEKey.h"
#include <vector>
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TProfile.h"

class StPicoDst;
class StAlexPhiMesonEvent;
class StAlexPhiMesonTrack;
class StV0TofCorrection;
class StV0Event;
class StV0Track;
class StTriFlowCut;
class TH1F;
class TH2F;
class TTree;
class TVector2;
class TVector3;

class StTriFlowV0
{
  public:
    StTriFlowV0(Int_t energy);
    ~StTriFlowV0();

    void InitPhi();

    void doPhi(Int_t,Int_t,Int_t,Int_t,Float_t,Float_t,Float_t);
    void MixEvent_Phi(Int_t,StPicoDst*,Int_t,Float_t,Float_t,Float_t,Float_t,Float_t);
    void clear_phi(Int_t,Int_t,Int_t);
    void size_phi(Int_t,Int_t,Int_t);

    void WritePhiMass2();

    void clearEvent();
    void passEvent(Int_t,Int_t,Int_t); // N_prim,N_non_prim,N_Tof_match
    void passEventPlane2East(TVector2,TVector2,TVector2,TVector2); // qVector ater re-center: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50
    void passEventPlane2West(TVector2,TVector2,TVector2,TVector2); // qVector ater re-center: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50
    void passEventPlane3East(TVector2,TVector2,TVector2,TVector2); // qVector ater re-center: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50
    void passEventPlane3West(TVector2,TVector2,TVector2,TVector2); // qVector ater re-center: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50
    void passNumTrackEast(Int_t,Int_t,Int_t,Int_t); // Number of East Track: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50
    void passNumTrackWest(Int_t,Int_t,Int_t,Int_t); // Number of West Track: eta_gap = 0.05, eta_gap = 0.10, eta_gap = 0.20, eta_gap = 0.50

    bool passPhiEtaEast(TLorentzVector,Int_t,Int_t);
    bool passPhiEtaWest(TLorentzVector,Int_t,Int_t);
  private:
    StTriFlowCut *mTriFlowCut;
    StV0TofCorrection *mTofCorr;
    TH2F *h_Mass2;
    TH2F *h_Mass2_rot;
    TH2F *h2_azimu_eta;
    TH2F *h2_pT_eta;
    TH2F *h2_pT_y;
    TH1F *h_mT;
    TH1F *h_dip_angle;
    TH1F *h_psi2_east;
    TH1F *h_psi2_west;
    TH1F *mHist_SE_InvM_ptSetA_centSetA[3][4];
    TH1F *mHist_rotation_InvM_ptSetA_centSetA[3][4];
    TProfile *mProfile_v2_reso_ptSetA_centSetA[3][4];

    Int_t mEventCounter2[9][10][5]; // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin

    // 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin, 3 = mixed event bin, 4 = charge bin(0 for pos, 1 for neg) || push_back->track
    vectorHelixMap mHelix_Kaon;
    vectorFloatMap mMomentum;
    vectorFloatMap mMass2;
    vectorFloatMap mDca;
    vectorFloatMap mNHitsFit;
    vectorFloatMap mNSigmaKaon;

    TTree *mTree_Phi;
    StAlexPhiMesonEvent *mXuPhiMesonEvent;
    StAlexPhiMesonTrack *mXuPhiMesonTrack;

    StV0Event *mV0Event;
    StV0Track *mV0Track;

    // event information | 0 = centrality bin, 1 = vertexZ bin, 2 = EP bin, 3 = eta_gap || push_back->event
    std::vector<TVector3> mPrimaryvertex[9][10][5];
    std::vector<Int_t> mRefMult[9][10][5];
    std::vector<Float_t> mReweight[9][10][5];
    std::vector<Int_t> mCentrality[9][10][5];
    std::vector<Int_t> mRunId[9][10][5];
    std::vector<Int_t> mN_prim[9][10][5];
    std::vector<Int_t> mN_non_prim[9][10][5];
    std::vector<Int_t> mN_Tof_match[9][10][5];
    std::vector<Float_t> mZDCx[9][10][5];
    std::vector<Float_t> mBBCx[9][10][5];
    std::vector<Float_t> mVzVpd[9][10][5];
    std::vector<Float_t> mField[9][10][5];
    std::vector<UShort_t> mNumTracks[9][10][5];
    std::vector<TVector2> mQ2East[9][10][5][4];
    std::vector<TVector2> mQ2West[9][10][5][4];
    std::vector<TVector2> mQ3East[9][10][5][4];
    std::vector<TVector2> mQ3West[9][10][5][4];
    std::vector<Int_t> mNumTrackEast[9][10][5][4];
    std::vector<Int_t> mNumTrackWest[9][10][5][4];

    // passing variable
    Int_t mNumber_prim, mNumber_non_prim, mNumber_Tof_match;
    TVector2 mQVector2East[4], mQVector2West[4], mQVector3East[4], mQVector3West[4];
    Int_t mTrackEtaEast[4], mTrackEtaWest[4];
    Int_t mEnergy;

  ClassDef(StTriFlowV0,1)
};
#endif
