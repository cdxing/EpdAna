#include "StTriFlowV0.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StRoot/StPicoEvent/StPicoDst.h"   //shaowei
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"  // shaowei
#include "StRoot/StAlexPhiMesonEvent/StAlexPhiMesonEvent.h"
//#include "StRoot/StV0Event/StV0Event.h"
//#include "StRoot/StV0TofCorrection/StV0TofCorrection.h"
#include <vector>
#include "TLorentzVector.h"
#include "StPicoEvent/StPicoPhysicalHelix.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#include "StarClassLibrary/SystemOfUnits.h"
#include "TMath.h"
#include "TObject.h"
#include "StLorentzVectorD.hh"

ClassImp(StTriFlowV0)

    //------------------------------------------------------------------------------------------------------------------
StTriFlowV0::StTriFlowV0(Int_t energy)
{
    mEnergy = energy;
}

StTriFlowV0::~StTriFlowV0()
{
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::InitPhi()
{
    mTriFlowCut = new StTriFlowCut(mEnergy);
    TString HistName = "Mass2_pt", HistName_rot = "Mass2_rot_pt";
    TString HistName_mt = "mT", HistName_dip = "dip_angle";
    TString HistName_azimuEta = "azimuth_eta", HistName_azimuthY = "azimuth_y";
    TString HistName_ptEta = "pT_eta", HistName_ptY = "pT_y";
    TString HistName_psi2East = "Psi2_East", HistName_psi2West = "Psi2_West";
    h_Mass2 = new TH2F(HistName.Data(),HistName.Data(),20,0.2,5.0,200,0.98,1.08);
    h_Mass2_rot = new TH2F(HistName_rot.Data(),HistName_rot.Data(),20,0.2,5.0,200,0.98,1.08);
    h2_azimu_eta = new TH2F(HistName_azimuEta.Data(),HistName_azimuEta.Data(),200,-2.0,2.0,200,-2.0*TMath::Pi(),2.0*TMath::Pi());
    h2_pT_eta = new TH2F(HistName_ptEta.Data(),HistName_ptEta.Data(),200,-2.0,2.0,20,0.2,5.0);
    h2_pT_y = new TH2F(HistName_ptY.Data(),HistName_ptY.Data(),200,-2.0,2.0,20,0.2,5.0);
    h_mT = new TH1F(HistName_mt.Data(),HistName_mt.Data(),200,0.0,10);
    h_dip_angle = new TH1F(HistName_dip.Data(),HistName_dip.Data(),200,-(1.0/4.0)*TMath::Pi(),1.0*TMath::Pi());
    h_psi2_east = new TH1F(HistName_psi2East.Data(),HistName_psi2East.Data(),200,-1.0*TMath::Pi(),1.0*TMath::Pi());
    h_psi2_west = new TH1F(HistName_psi2West.Data(),HistName_psi2West.Data(),200,-1.0*TMath::Pi(),1.0*TMath::Pi());
    for(Int_t cent = 0; cent < TriFlow::Bin_Centrality_01; cent++)
    {
        for(Int_t pt_bin = 0; pt_bin < TriFlow::Bin_pT; pt_bin++)
        {
          TString hist_name_SE = Form("InvMass_SE_ptbin%d_cent%s",pt_bin+1,TriFlow::Centrality_01[cent].Data());
          mHist_SE_InvM_ptSetA_centSetA[pt_bin][cent] = new TH1F(hist_name_SE.Data() ,
          hist_name_SE.Data() ,
          200,0.98,1.08);
          mHist_SE_InvM_ptSetA_centSetA[pt_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
          TString hist_name_rot = Form("InvMass_rot_ptbin%d_cent%s",pt_bin+1,TriFlow::Centrality_01[cent].Data());
          mHist_rotation_InvM_ptSetA_centSetA[pt_bin][cent] = new TH1F(hist_name_rot.Data() ,
          hist_name_rot.Data() ,
          200,0.98,1.08);
          mHist_rotation_InvM_ptSetA_centSetA[pt_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
          TString hist_name_profile = Form("flow_InvMass_ptbin%d_cent%s",pt_bin+1,TriFlow::Centrality_01[cent].Data());
          mProfile_v2_reso_ptSetA_centSetA[pt_bin][cent] = new TProfile(hist_name_profile.Data(),
          hist_name_profile.Data(),
          100,0.98,1.08,
          0,0,"");
          mProfile_v2_reso_ptSetA_centSetA[pt_bin][cent]->GetXaxis()->SetTitle("m_{inv} [GeV/c^{2}]");
          mProfile_v2_reso_ptSetA_centSetA[pt_bin][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
    }

    for(Int_t cent = 0; cent < TriFlow::Bin_Centrality; cent++)
    {
        for(Int_t vz = 0; vz < TriFlow::Bin_VertexZ; vz++)
        {
            for(Int_t phi_psi = 0; phi_psi < TriFlow::Bin_Phi_Psi; phi_psi++)
            {
                mEventCounter2[cent][vz][phi_psi] = 0;
                clear_phi(cent,vz,phi_psi);
            }
        }
    }

    mXuPhiMesonEvent = new StAlexPhiMesonEvent();
    mTree_Phi = new TTree("XuPhiMesonEvent","XuPhiMesonEvent");
    mTree_Phi->Branch("phi_flow_branch","StAlexPhiMesonEvent",&mXuPhiMesonEvent);
    mTree_Phi->SetAutoSave(5000000);
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::WritePhiMass2()
{
    h_Mass2->Write();
    h_Mass2_rot->Write();
    h2_azimu_eta->Write();
    h2_pT_eta->Write();
    h2_pT_y->Write();
    h_mT->Write();
    h_dip_angle->Write();
    h_psi2_east->Write();
    h_psi2_west->Write();
    for(Int_t cent = 0; cent < TriFlow::Bin_Centrality_01; cent++)
    {
        for(Int_t pt_bin = 0; pt_bin < TriFlow::Bin_pT; pt_bin++)
        {
          mHist_SE_InvM_ptSetA_centSetA[pt_bin][cent]->Write();
          mHist_rotation_InvM_ptSetA_centSetA[pt_bin][cent]->Write();
          mProfile_v2_reso_ptSetA_centSetA[pt_bin][cent]->Write();
        }
    }
    //mTree_Phi->Write("",TObject::kOverwrite);
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::clear_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
    mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].clear();
    mRefMult[cent9][Bin_vz][Bin_Psi2].clear();
    mReweight[cent9][Bin_vz][Bin_Psi2].clear();
    mCentrality[cent9][Bin_vz][Bin_Psi2].clear();
    mRunId[cent9][Bin_vz][Bin_Psi2].clear();
    mN_prim[cent9][Bin_vz][Bin_Psi2].clear();
    mN_non_prim[cent9][Bin_vz][Bin_Psi2].clear();
    mN_Tof_match[cent9][Bin_vz][Bin_Psi2].clear();
    mZDCx[cent9][Bin_vz][Bin_Psi2].clear();
    mBBCx[cent9][Bin_vz][Bin_Psi2].clear();
    mVzVpd[cent9][Bin_vz][Bin_Psi2].clear();

    for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
    {
        mQ2East[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ2West[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ3East[cent9][Bin_vz][Bin_Psi2][j].clear();
        mQ3West[cent9][Bin_vz][Bin_Psi2][j].clear();
    }

    for(Int_t Bin_Event = 0; Bin_Event < TriFlow::Buffer_depth; Bin_Event++)
    {
        for(Int_t charge = 0; charge < 2; charge++)
        {
            MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
            mHelix_Kaon[key].clear();
            mMomentum[key].clear();
            mMass2[key].clear();
            mDca[key].clear();
            mNHitsFit[key].clear();
            mNSigmaKaon[key].clear();
        }
    }
    mEventCounter2[cent9][Bin_vz][Bin_Psi2] = 0;
}

void StTriFlowV0::size_phi(Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2)
{
    LOG_INFO << "Event Buffer: Centrality = " << cent9 << ", VertexZ = " << Bin_vz << ", psi2_east = " << Bin_Psi2 << endm;
    LOG_INFO << "Buffer_depth = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endm;

    LOG_INFO << "Size of primaryVertex = " << mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    LOG_INFO << "Size of refMult       = " << mRefMult[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    LOG_INFO << "Size of reweight       = " << mReweight[cent9][Bin_vz][Bin_Psi2].size() << endm;;
    LOG_INFO << "---------------------------------------------------------------------------" << endm;

    for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
    {
        MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
        LOG_INFO << "Event Number " << Bin_Event << ":" << endm;
        LOG_INFO << "Positive Particle:" << endm;
        LOG_INFO << "  Size of Helix_Kplus  = " << mHelix_Kaon[key].size() << endm;;
        LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
        LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
        LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
        LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
        LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;

        key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);
        LOG_INFO << "Negative Particle:" << endm;
        LOG_INFO << "  Size of Helix_Kminus = " << mHelix_Kaon[key].size() << endm;
        LOG_INFO << "  Size of Momentum     = " << mMomentum[key].size() << endm;
        LOG_INFO << "  Size of Mass2        = " << mMass2[key].size() << endm;
        LOG_INFO << "  Size of dca          = " << mDca[key].size() << endm;
        LOG_INFO << "  Size of nHitsFit     = " << mNHitsFit[key].size() << endm;
        LOG_INFO << "  Size of nSigmaKaon   = " << mNSigmaKaon[key].size() << endm;
        LOG_INFO << "---------------------------------------------------------------------------" << endm;
    }
}

//------------------------------------------------------------------------------------------------------------------

void StTriFlowV0::doPhi(Int_t Flag_ME, Int_t cent9, Int_t Bin_vz, Int_t Bin_Psi2, Float_t Psi2_East, Float_t Psi2_West, Float_t Res_EP) // 0: Same Event, 1: Mix Event
{
    if(Flag_ME == 0) // same event
    {
            //cout << "Bin_Event = " << mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endl;
        for(Int_t Bin_Event = 0; Bin_Event < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event++)
        {
            // event header
            mXuPhiMesonEvent->clearTrackList();
            mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setReweight(mReweight[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

            for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
            {
                // QVector
                mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                // Number of Tracks
                mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
            }

            mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
            mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

            // start to select phi candidate in a event
            MEKey key_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,0);
            MEKey key_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,1);

            TLorentzVector ltrackA, ltrackB;
            //cout <<"K+ number = " <<  mHelix_Kaon[key_plus].size() << endl;
            for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_plus].size(); n_kplus++) // first track loop over K+ candidates
            {
                TVector3 p_vecA = mHelix_Kaon[key_plus][n_kplus].cat(mHelix_Kaon[key_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
                p_vecA *= mMomentum[key_plus][n_kplus];
                ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),TriFlow::mMassKaon);
                Double_t d_ptA = ltrackA.Perp(), d_pzA = ltrackA.Pz(), d_momA = ltrackA.P();

                for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_minus].size(); n_kminus++) // second track loop over K- candidates
                {
                    TVector3 p_vecB = mHelix_Kaon[key_minus][n_kminus].cat(mHelix_Kaon[key_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]));  // primary momentum
                    p_vecB *= mMomentum[key_minus][n_kminus];
                    ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon);
                    Double_t d_ptB = ltrackB.Perp(), d_pzB = ltrackB.Pz(), d_momB = ltrackB.P();

                    Double_t d_dip_angle = TMath::ACos((d_ptA*d_ptB+d_pzA*d_pzB) / (d_momA*d_momB) );
                    TLorentzVector trackAB      = ltrackA+ltrackB;
                    Double_t InvMassAB          = trackAB.M();
                    Double_t pt = trackAB.Perp();
                    Double_t eta = trackAB.Eta();
                    Double_t rap = trackAB.Rapidity();
                    Double_t d_mT_phi = sqrt(pt*pt + TriFlow::mMassPhi*TriFlow::mMassPhi );
                    // do the rotation for track B
                    Double_t randomNumber = gRandom->Uniform(1);
                    // std::cout << "randomNumber = " << randomNumber << std::endl;
                    Double_t d_randAngle = TMath::Pi()*randomNumber;

                    TLorentzVector ltrackB_rot = ltrackB;
                    ltrackB_rot.RotateZ(d_randAngle);
                    TLorentzVector trackAB_rot      = ltrackA+ltrackB_rot;
                    Double_t InvMassAB_rot          = trackAB_rot.M();
                    Double_t pt_rot = trackAB_rot.Perp();
              // fill phi candidate into mTree_Phi
                    if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05)
                    {
                        mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                        mXuPhiMesonTrack->setMass2A(mMass2[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setMass2B(mMass2[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setDcaA(mDca[key_plus][n_kplus]); // K+
                        mXuPhiMesonTrack->setDcaB(mDca[key_minus][n_kminus]); // K-
                        mXuPhiMesonTrack->setTrackA(ltrackA); // K+
                        mXuPhiMesonTrack->setTrackB(ltrackB); // K-
                        mXuPhiMesonTrack->setFlagA(Bin_Event); // K+
                        mXuPhiMesonTrack->setFlagB(Bin_Event); // K-
                    }
                    Float_t phi_Psi2 = -999.0, flow2 = -999.0;
                    if(passPhiEtaEast(trackAB,0,1))
                    {
                        Float_t phi_East = trackAB.Phi();
                        phi_Psi2 = phi_East - Psi2_West;
                        flow2 = TMath::Cos( 2. * phi_Psi2 )/(Float_t)Res_EP;
                        // Float_t phi_Psi3 = phi_East - Psi3_West;
                        // std::cout << "phi = " << phi_East << std::endl;
                        // std::cout << "psi = " << Psi2_West << std::endl;
                        // std::cout << "flow2 = " << flow2 << std::endl;
                    }
                    if(passPhiEtaWest(trackAB,0,1))
                    {
                        Float_t phi_West = trackAB.Phi();
                        phi_Psi2 = phi_West - Psi2_East;
                        flow2 = TMath::Cos( 2. * phi_Psi2 )/(Float_t)Res_EP;
                        // Float_t phi_Psi3 = phi_West - Psi3_East;
                        // std::cout << "phi = " << phi_West << std::endl;
                        // std::cout << "psi = " << Psi2_East << std::endl;
                        // std::cout << "flow2 = " << flow2 << std::endl;
                    }
                    // Fill histogram with InvMassAB information
                    h_Mass2->Fill(pt,InvMassAB);
                    h_Mass2_rot->Fill(pt_rot,InvMassAB_rot);
                    // Fill histograms with for QA cuts: dip angle, etc.
                    h2_pT_y->Fill(rap,pt);
                    h2_azimu_eta->Fill(eta,trackAB.Phi());
		    h2_pT_eta->Fill(eta,pt);
                    h_mT->Fill(d_mT_phi);
                    h_dip_angle->Fill(d_dip_angle);
                    h_psi2_east->Fill(Psi2_East);
                    h_psi2_west->Fill(Psi2_West);
                    // Fill hisograms for invM fit Method
                    if(mTriFlowCut->passDipAngle(d_dip_angle))
                    {
                      for(Int_t cent = 0; cent < TriFlow::Bin_Centrality_01; cent++)
                      {
                          for(Int_t pt_bin = 0; pt_bin < TriFlow::Bin_pT; pt_bin++)
                          {
                            if(TriFlow::cent_low[cent]<= cent9 && cent9 <= TriFlow::cent_up[cent] &&
                               TriFlow::pt_low_phi[pt_bin] <= pt && pt < TriFlow::pt_up_phi[pt_bin])
                               {
                                 mHist_SE_InvM_ptSetA_centSetA[pt_bin][cent]->Fill(InvMassAB);
                                 // std::cout << "invM = " << InvMassAB << std::endl;
                                 if(!(phi_Psi2 == -999.0 || Res_EP == -999.0 || flow2 == -999.0))
                                 {
                                   // if(pt_bin==0)std::cout << "phi_Psi2  = " << phi_Psi2 << std::endl;
                                   // if(pt_bin==0)std::cout << "Res_EP  = " << Res_EP << std::endl;
                                   // if(pt_bin==0)std::cout << "ptbin1 flow2 = " << flow2 << std::endl;
                                   mProfile_v2_reso_ptSetA_centSetA[pt_bin][cent]->Fill(InvMassAB,flow2);
                                 }
                               }
                            if(TriFlow::cent_low[cent]<= cent9 && cent9 <= TriFlow::cent_up[cent] &&
                               TriFlow::pt_low_phi[pt_bin] <= pt_rot && pt_rot <= TriFlow::pt_up_phi[pt_bin])
                               {
                                 // std::cout << "invM rot = " << InvMassAB_rot << std::endl;
                                 mHist_rotation_InvM_ptSetA_centSetA[pt_bin][cent]->Fill(InvMassAB_rot);
                               }

                          }
                      }
                    }
                }
            }
        }
        //mTree_Phi->Fill();
    }

    if(Flag_ME == 1) // mixed event
    {
        for(Int_t Bin_Event_A = 0; Bin_Event_A < mEventCounter2[cent9][Bin_vz][Bin_Psi2]-1; Bin_Event_A++)
        {
            MEKey key_A_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,0);
            MEKey key_A_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_A,1);
            for(Int_t Bin_Event_B = Bin_Event_A+1; Bin_Event_B < mEventCounter2[cent9][Bin_vz][Bin_Psi2]; Bin_Event_B++)
            {
                MEKey key_B_plus  = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,0);
                MEKey key_B_minus = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event_B,1);

                if(Bin_Event_A == 0 && Bin_Event_B == 1)
                {
                    Int_t Bin_Event = Bin_Event_A;
                    // event header
                    mXuPhiMesonEvent->clearTrackList();
                    mXuPhiMesonEvent->setPrimaryVertex(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setRunId(mRunId[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setRefMult(mRefMult[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setReweight(mReweight[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setCentrality(mCentrality[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setN_prim(mN_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setN_non_prim(mN_non_prim[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setN_Tof_match(mN_Tof_match[cent9][Bin_vz][Bin_Psi2][Bin_Event]);

                    for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
                    {
                        // QVector
                        mXuPhiMesonEvent->setQ2East(mQ2East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ2West(mQ2West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ3East(mQ3East[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setQ3West(mQ3West[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);

                        // Number of Tracks
                        mXuPhiMesonEvent->setNumTrackEast(mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                        mXuPhiMesonEvent->setNumTrackWest(mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j][Bin_Event],j);
                    }

                    mXuPhiMesonEvent->setZDCx(mZDCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setBBCx(mBBCx[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                    mXuPhiMesonEvent->setVzVpd(mVzVpd[cent9][Bin_vz][Bin_Psi2][Bin_Event]);
                }

                TLorentzVector ltrackA, ltrackB;

                // start to mix events
                // mix K+ candidates from A event with K- candidates from B event
                for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_A_plus].size(); n_kplus++) // first track loop over K+ candidates from event A
                {
                    TVector3 p_vecA(mHelix_Kaon[key_A_plus][n_kplus].cat(mHelix_Kaon[key_A_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
                    p_vecA *= mMomentum[key_A_plus][n_kplus];
                    ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),TriFlow::mMassKaon); // K+
                    Double_t d_ptA = ltrackA.Perp(), d_pzA = ltrackA.Pz(), d_momA = ltrackA.P();

                    for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_B_minus].size(); n_kminus++) // second track loop over K- candidates from event B
                    {
                        TVector3 p_vecB(mHelix_Kaon[key_B_minus][n_kminus].cat(mHelix_Kaon[key_B_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
                        p_vecB *= mMomentum[key_B_minus][n_kminus];
                        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon); // K-
                        Double_t d_ptB = ltrackB.Perp(), d_pzB = ltrackB.Pz(), d_momB = ltrackB.P();

                        Double_t d_dip_angle = TMath::ACos((d_ptA*d_ptB+d_pzA*d_pzB) / (d_momA*d_momB) );
                        TLorentzVector trackAB      = ltrackA+ltrackB;
                        Double_t InvMassAB          = trackAB.M();
                        Double_t pt = trackAB.Perp();

                        // fill phi candidate background into mTree_Phi
                        if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05)
                        {
                            mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                            mXuPhiMesonTrack->setMass2A(mMass2[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setMass2B(mMass2[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setDcaA(mDca[key_A_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setDcaB(mDca[key_B_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setTrackA(ltrackA); // K+
                            mXuPhiMesonTrack->setTrackB(ltrackB); // K-
                            mXuPhiMesonTrack->setFlagA(Bin_Event_A); // K+
                            mXuPhiMesonTrack->setFlagB(Bin_Event_B); // K-
                        }

                        // Fill histogram with InvMassAB information
                        h_Mass2->Fill(pt,InvMassAB);
                        // Fill histograms with for QA cuts: dip angle, etc.
                        h_dip_angle->Fill(d_dip_angle);
                    }
                }

                // mix K- candidates from A event with K+ candidates from B event
                for(Int_t n_kminus = 0; n_kminus < mHelix_Kaon[key_A_minus].size(); n_kminus++) // first track loop over K- candidates from event A
                {
                    TVector3 p_vecA(mHelix_Kaon[key_A_minus][n_kminus].cat(mHelix_Kaon[key_A_minus][n_kminus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_A])));  // primary momentum
                    p_vecA *= mMomentum[key_A_minus][n_kminus];
                    ltrackA.SetXYZM(p_vecA.X(),p_vecA.Y(),p_vecA.Z(),TriFlow::mMassKaon); // K-
                    Double_t d_ptA = ltrackA.Perp(), d_pzA = ltrackA.Pz(), d_momA = ltrackA.P();

                    for(Int_t n_kplus = 0; n_kplus < mHelix_Kaon[key_B_plus].size(); n_kplus++) // second track loop over K+ candidates from event B
                    {
                        TVector3 p_vecB(mHelix_Kaon[key_B_plus][n_kplus].cat(mHelix_Kaon[key_B_plus][n_kplus].pathLength(mPrimaryvertex[cent9][Bin_vz][Bin_Psi2][Bin_Event_B])));  // primary momentum
                        p_vecB *= mMomentum[key_B_plus][n_kplus];
                        ltrackB.SetXYZM(p_vecB.x(),p_vecB.y(),p_vecB.z(),TriFlow::mMassKaon); // K+
                        Double_t d_ptB = ltrackB.Perp(), d_pzB = ltrackB.Pz(), d_momB = ltrackB.P();

                        Double_t d_dip_angle = TMath::ACos((d_ptA*d_ptB+d_pzA*d_pzB) / (d_momA*d_momB) );
                        TLorentzVector trackAB      = ltrackA+ltrackB;
                        Double_t InvMassAB          = trackAB.M();
                        Double_t pt = trackAB.Perp();

                        // fill phi candidate background into mTree_Phi
                        if(InvMassAB > TriFlow::mMassKaon*2 && InvMassAB < 1.05)
                        {
                            mXuPhiMesonTrack = mXuPhiMesonEvent->createTrack();
                            mXuPhiMesonTrack->setMass2A(mMass2[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setMass2B(mMass2[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setNSigKaonA(mNSigmaKaon[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setNSigKaonB(mNSigmaKaon[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setDcaA(mDca[key_B_plus][n_kplus]); // K+
                            mXuPhiMesonTrack->setDcaB(mDca[key_A_minus][n_kminus]); // K-
                            mXuPhiMesonTrack->setTrackA(ltrackB); // K+
                            mXuPhiMesonTrack->setTrackB(ltrackA); // K-
                            mXuPhiMesonTrack->setFlagA(Bin_Event_B); // K+
                            mXuPhiMesonTrack->setFlagB(Bin_Event_A); // K-
                        }

                        // Fill histogram with InvMassAB information
                        h_Mass2->Fill(pt,InvMassAB);
                        // Fill histograms with for QA cuts: dip angle, etc.
                        h_dip_angle->Fill(d_dip_angle);
                    }
                }
            }
        }
        //mTree_Phi->Fill();
    }
}

//------------------------------------------------------------------------------------------------------------------


void StTriFlowV0::MixEvent_Phi(Int_t Flag_ME, StPicoDst *pico, Int_t cent9, Float_t vz, Float_t psi2_east, Float_t psi2_west,Float_t reweight,Float_t resolution2)
{
    StPicoEvent *event = (StPicoEvent*)pico->event();

    Int_t Bin_vz, Bin_Psi2;

    //Float_t vz_start = TriFlow::mVzMaxMap[event->energy()];
    Float_t vz_start = 198.0; //   shaowei   same with 39 62GeV
    Float_t vz_bin = 4./TriFlow::Bin_VertexZ;

    Float_t psi2_start = TMath::Pi()/2.0;
    Float_t psi2_bin = 2*psi2_start/TriFlow::Bin_Phi_Psi;

    for(Int_t i = 0; i < TriFlow::Bin_VertexZ; i++)
    {
        if((vz > vz_start+i*vz_bin) && (vz <= vz_start+(i+1)*vz_bin))
        {
            Bin_vz = i;
        }
    }
    for(Int_t i = 0; i < TriFlow::Bin_Phi_Psi; i++)
    {
        if((psi2_east > -1.0*psi2_start+i*psi2_bin) && (psi2_east <= -1.0*psi2_start+(i+1)*psi2_bin))
        {
            Bin_Psi2 = i;
        }
    }

    //mEventCounter2[cent9][Bin_vz][Bin_Psi2] =0;
    Int_t Bin_Event = mEventCounter2[cent9][Bin_vz][Bin_Psi2];

    const Double_t MAGFIELDFACTOR = kilogauss;
    const Int_t nTracks = pico->numberOfTracks();

    // store Enent Information
    mPrimaryvertex[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<TVector3>(event->primaryVertex()));
    mRefMult[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->refMult()));
    mReweight[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(reweight));
    mCentrality[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(cent9));
    mRunId[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(event->runId()));
    mN_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_prim));
    mN_non_prim[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_non_prim));
    mN_Tof_match[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Int_t>(mNumber_Tof_match));
    mZDCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->ZDCx()));
    mBBCx[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->BBCx()));
    mVzVpd[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<Float_t>(event->vzVpd()));
    mNumTracks[cent9][Bin_vz][Bin_Psi2].push_back(static_cast<UShort_t>(pico->numberOfTracks()));
    for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
    {
        mQ2East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2East[j]));
        mQ2West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector2West[j]));
        mQ3East[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3East[j]));
        mQ3West[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<TVector2>(mQVector3West[j]));
        mNumTrackEast[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaEast[j]));
        mNumTrackWest[cent9][Bin_vz][Bin_Psi2][j].push_back(static_cast<Int_t>(mTrackEtaWest[j]));
    }

    // store Track Information
    TVector3 mVertexPos = event->primaryVertex();
    float mField = event->bField();
    // cout << "Res2_EP = " << resolution2 << endl;
    //cout << "total track = " << nTracks << endl;
    //cout << "total tof track =  " << pico->numberOfBTofPidTraits() << endl;
    for(Int_t i = 0; i < nTracks; i ++) // loop over all particles in event
    {
        StPicoTrack *track = pico->track(i);
        if(!track)
        {
            cout << "No PicoTrack Skip! " << endl;
            return;
        }
        StPicoPhysicalHelix helix = track->helix(mField);
        //Float_t dca = helix.geometricSignedDistance(mVertexPos);
        Float_t dca=track->gDCA(mVertexPos).Mag();

        //cout << "mBTofPidTraitsIndex = " << track->bTofPidTraitsIndex() << endl;

        if(mTriFlowCut->passTrackPhi(track,dca))     //shaowei
        {
            Float_t Mass2 = mTriFlowCut->getMass2(track, pico); // shaowei
            //Float_t scale_nSigma_factor = TriFlow::mSigScaleMap[event->energy()];
            Float_t scale_nSigma_factor = 1.0;    // shaowei
            Float_t Polarity = static_cast<Float_t>(track->charge());
            Float_t momentum = track->pMom().Mag();
            Float_t Mass2_low;
            Float_t Mass2_up;
            if(momentum < 0.5)
            {
                Mass2_low = 0.4*0.4;
                Mass2_up = 0.6*0.6;
            }
            if(momentum >= 0.5)
            {
                Mass2_low = 0.277205 - 0.0812931*momentum;
                Mass2_up = 0.215517 + 0.076801*momentum;
            }

            Int_t charge = 0; // k+
            if(Polarity < 0) charge = 1; // k-


            if(mTriFlowCut->passSigKaonCut(track,scale_nSigma_factor))
            {
               if(
                        (momentum < 0.65 && ((Mass2 > Mass2_low && Mass2 < Mass2_up) || Mass2 < -10.0)) // dE/dx + ToF
                        || (momentum >= 0.65 && (Mass2 > Mass2_low && Mass2 < Mass2_up)) // dE/dx + ToF(always)
                  )
                {
                    MEKey key = MEKey(cent9,Bin_vz,Bin_Psi2,Bin_Event,charge);
                    //mMass2[key].push_back(static_cast<Float_t>(mTriFlowCut->getMass2(track, tof))); // mass2
                    mMass2[key].push_back(static_cast<Float_t>(Mass2)); // mass2
                    mDca[key].push_back(static_cast<Float_t>(dca*track->charge())); // dca*charge //shaowei
                    mNHitsFit[key].push_back(static_cast<Float_t>(track->nHitsFit())); // nHitsFit
                    mNSigmaKaon[key].push_back(static_cast<Float_t>((track->nSigmaKaon())*scale_nSigma_factor)); // nSigmaKaon
                    //mHelix_Kaon[key].push_back(static_cast<StPhysicalHelixD>(StPhysicalHelixD(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom
                    mHelix_Kaon[key].push_back(static_cast<StPicoPhysicalHelix>(StPicoPhysicalHelix(track->pMom(),event->primaryVertex(),event->bField()*MAGFIELDFACTOR,track->charge())));// get helix from the pMom
                    mMomentum[key].push_back(static_cast<Float_t>(track->pMom().Mag()));// get helix from the pMom
                }
            }
        }
    }

    mEventCounter2[cent9][Bin_vz][Bin_Psi2]++;
    //cout << "mEventCounters = "<<mEventCounter2[cent9][Bin_vz][Bin_Psi2] << endl;

    if(Flag_ME == 0) // same event
    {
        doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2,psi2_east,psi2_west,resolution2);
        clear_phi(cent9,Bin_vz,Bin_Psi2);
    }

    if(Flag_ME == 1) // mix event
    {
        if(mEventCounter2[cent9][Bin_vz][Bin_Psi2] == TriFlow::Buffer_depth)
        {
            doPhi(Flag_ME,cent9,Bin_vz,Bin_Psi2,psi2_east,psi2_west,resolution2);
            clear_phi(cent9,Bin_vz,Bin_Psi2);
        }
    }
}

//------------------------------------------------------------------------------------------------------------------

// pass event information from Maker
void StTriFlowV0::clearEvent()
{
    mNumber_prim = 0;
    mNumber_non_prim = 0;
    mNumber_Tof_match = 0;

    for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
    {
        mQVector2East[j].Set(-999.9,-999.9);
        mQVector2West[j].Set(-999.9,-999.9);
        mQVector3East[j].Set(-999.9,-999.9);
        mQVector3West[j].Set(-999.9,-999.9);
    }
}

void StTriFlowV0::passEvent(Int_t N_prim, Int_t N_non_prim, Int_t N_Tof_match)
{
    mNumber_prim = N_prim;
    mNumber_non_prim = N_non_prim;
    mNumber_Tof_match = N_Tof_match;
}

void StTriFlowV0::passEventPlane2East(TVector2 Q2East_0, TVector2 Q2East_1, TVector2 Q2East_2, TVector2 Q2East_3)
{
    mQVector2East[0] = Q2East_0;
    mQVector2East[1] = Q2East_1;
    mQVector2East[2] = Q2East_2;
    mQVector2East[3] = Q2East_3;
}

void StTriFlowV0::passEventPlane2West(TVector2 Q2West_0, TVector2 Q2West_1, TVector2 Q2West_2, TVector2 Q2West_3)
{
    mQVector2West[0] = Q2West_0;
    mQVector2West[1] = Q2West_1;
    mQVector2West[2] = Q2West_2;
    mQVector2West[3] = Q2West_3;
}

void StTriFlowV0::passEventPlane3East(TVector2 Q3East_0, TVector2 Q3East_1, TVector2 Q3East_2, TVector2 Q3East_3)
{
    mQVector3East[0] = Q3East_0;
    mQVector3East[1] = Q3East_1;
    mQVector3East[2] = Q3East_2;
    mQVector3East[3] = Q3East_3;
}

void StTriFlowV0::passEventPlane3West(TVector2 Q3West_0, TVector2 Q3West_1, TVector2 Q3West_2, TVector2 Q3West_3)
{
    mQVector3West[0] = Q3West_0;
    mQVector3West[1] = Q3West_1;
    mQVector3West[2] = Q3West_2;
    mQVector3West[3] = Q3West_3;
}

void StTriFlowV0::passNumTrackEast(Int_t NumTrackEast_0, Int_t NumTrackEast_1, Int_t NumTrackEast_2, Int_t NumTrackEast_3)
{
    mTrackEtaEast[0] = NumTrackEast_0;
    mTrackEtaEast[1] = NumTrackEast_1;
    mTrackEtaEast[2] = NumTrackEast_2;
    mTrackEtaEast[3] = NumTrackEast_3;
}

void StTriFlowV0::passNumTrackWest(Int_t NumTrackWest_0, Int_t NumTrackWest_1, Int_t NumTrackWest_2, Int_t NumTrackWest_3)
{
    mTrackEtaWest[0] = NumTrackWest_0;
    mTrackEtaWest[1] = NumTrackWest_1;
    mTrackEtaWest[2] = NumTrackWest_2;
    mTrackEtaWest[3] = NumTrackWest_3;
}
//---------------------------------------------------------------------------------

bool StTriFlowV0::passPhiEtaEast(TLorentzVector track, Int_t i, Int_t Mode) // neg || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
    Float_t eta = track.Eta();

    if(Mode == 0) // Event Plane Mode
    {
        // eta cut
        // eta_gap between two sub event plane is 2*mEta_Gap[i]
        if(!(eta > -1.0*TriFlow::mEtaMax && eta < -1.0*TriFlow::mEta_Gap[i]))
        {
            return kFALSE;
        }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > -1.0*TriFlow::mEtaMax && eta < 0.0))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTriFlowV0::passPhiEtaWest(TLorentzVector track, Int_t i, Int_t Mode) // pos || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track.Eta();

  if(Mode == 0) // Event Plane Mode
  {
    // eta cut
    // eta_gap between two sub event plane is 2*mEta_Gap[i]
    if(!(eta > TriFlow::mEta_Gap[i] && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }
  if(Mode == 1) // Flow Mode
  {
    // eta cut
    // eta_gap between two sub event plane is mEta_Gap[i]
    if(!(eta > 0.0 && eta < TriFlow::mEtaMax))
    {
      return kFALSE;
    }

    return kTRUE;
  }

  return kTRUE;
}

//------------------------------------------------------------------------------------------------------------------
