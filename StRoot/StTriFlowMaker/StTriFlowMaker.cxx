#include "StTriFlowMaker.h"
#include "StTriFlowConstants.h"
#include "StTriFlowCut.h"
#include "StTriFlowCorrection.h"
#include "StTriFlowV0.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
//#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "StRoot/run/run.h"
//#include "StRoot/run/badrun.h"
#include "StThreeVectorF.hh"
#include "TH1F.h"
#include "TProfile.h"
#include "TH2F.h"
#include "TFile.h"
#include "TVector3.h"
#include "TMath.h"
#include "StMessMgr.h"
#include <algorithm>
#include <array>
ClassImp(StTriFlowMaker)

    //StRefMultCorr* StTriFlowMaker::mRefMultCorr = NULL;
    //-----------------------------------------------------------------------------
    StTriFlowMaker::StTriFlowMaker(const char* name, StPicoDstMaker *picoMaker, char* jobid)
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
    mEnergy = 3;

    mOutPut_Phi = Form("%s_.root",jobid);

}

//-----------------------------------------------------------------------------
StTriFlowMaker::~StTriFlowMaker()
{ /*  */ }

//-----------------------------------------------------------------------------
Int_t StTriFlowMaker::Init()
{

    //if(!mRefMultCorr)
    //{
    //    mRefMultCorr = CentralityMaker::instance()->getRefMultCorr();
    //}


    mTriFlowCut = new StTriFlowCut(mEnergy);
    mTriFlowCorrection = new StTriFlowCorrection(mEnergy);

    mTriFlowV0 = new StTriFlowV0(mEnergy);
    mTriFlowV0->InitPhi();
    mTriFlowCorrection->InitReCenterCorrection(mEnergy);
    mTriFlowCorrection->InitShiftCorrection(mEnergy);
    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StTriFlowMaker::Finish()
{
    mFile_Phi = new TFile(mOutPut_Phi.Data(),"RECREATE");
    //mFile_Phi->cd();
    if(mOutPut_Phi != "")
    {
        mFile_Phi->cd();
        mTriFlowV0->WritePhiMass2();
        mFile_Phi->Close();
    }
    return kStOK;
}

//-----------------------------------------------------------------------------
void StTriFlowMaker::Clear(Option_t *opt) {
}

//-----------------------------------------------------------------------------
Int_t StTriFlowMaker::Make()
{
    if(!mPicoDstMaker)
    {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst)
    {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent)
    {
        LOG_WARN << " No PicoEvent! Skip! " << endm;
        return kStWarn;
    }

    // RefMult
    Int_t runId = mPicoEvent->runId();
    /*
    for(int ii=0; ii<286; ii++)
    {
        //if(runId == badrun[ii]) return kStOK;
    }
    */
    //cout << "runID = " << runId << endl;
    Int_t refMult = mPicoEvent->refMult();
    Float_t vz = mPicoEvent->primaryVertex().Z();
    Float_t vx = mPicoEvent->primaryVertex().X();
    Float_t vy = mPicoEvent->primaryVertex().Y();

    Float_t vzvpd = mPicoEvent->vzVpd();
    Int_t TOF_Mul = mPicoEvent->btofTrayMultiplicity();
    Float_t zdcX = mPicoEvent->ZDCx();

    // vz sign
    Int_t vz_sign;
    if(vz > 0.0)
    {
        vz_sign = 0;
    }
    else
    {
        vz_sign = 1;
    }

    // runIndex
    const int runIndex = GetRunIndex(runId);
    // Event Cut
    if(mTriFlowCut->passEventCut(mPicoDst)) // event cut
    {
        //mRefMultCorr->init(runId);
        //mRefMultCorr->initEvent(refMult, vz);
        //const Int_t cent9 = mRefMultCorr->getCentralityBin9();
        //const Double_t reweight = mRefMultCorr->getWeight();
        const Int_t cent9 = Centrality(refMult);
        const Double_t reweight = 1.0;
        if(cent9 >  8 || cent9 < 0) return 0;

        const Int_t nTracks = mPicoDst->numberOfTracks();
        float nMatchedToF=0, nN_prim=0, nN_non_prim=0;

        TVector3 mVertexPos = mPicoDst->event()->primaryVertex();
        float mField = mPicoEvent->bField();

        nMatchedToF = mPicoEvent->nBTOFMatch();
       // cout << "nTracks = " << nTracks << endl;
        for(Int_t i = 0; i < nTracks; i++) // track loop
        {
            StPicoTrack *track = (StPicoTrack*)mPicoDst->track(i);
            //StPicoPhysicalHelix helix = track->helix(mField);
            //Float_t dca = helix.geometricSignedDistance(mVertexPos);
            Float_t dca=track->gDCA(mVertexPos).Mag();
            if(dca > 3){// non primary track
                nN_non_prim++;
            }
            else{
                nN_prim++;
            }
            if(mTriFlowCut->passTrackEP(track,dca)) // track cut   //shaowei
            {
	//	cout << "pass track EP cut" << endl;
                // calculate Q Vector after recentering for full event and eta sub event
                if(mTriFlowCorrection->passTrackFull(track))
                {
                    mTriFlowCorrection->addTrack_Full(track,cent9,runIndex,vz_sign);
                }
                for(Int_t j = 0; j < 4; j++) // eta_gap loop
                {
                    if(mTriFlowCorrection->passTrackEtaEast(track,j,0)) // neg eta sub
                    {
                        mTriFlowCorrection->addTrack_East(track,cent9,runIndex,vz_sign,j);
                    }
                    if(mTriFlowCorrection->passTrackEtaWest(track,j,0)) // pos eta sub
                    {
                        mTriFlowCorrection->addTrack_West(track,cent9,runIndex,vz_sign,j);
                    }
                }
            }
        }
        Float_t Psi2_East, Psi2_West, Res2_EP;
        TVector2 Q2East[TriFlow::EtaGap_total], Q2West[TriFlow::EtaGap_total];
        TVector2 Q3East[TriFlow::EtaGap_total], Q3West[TriFlow::EtaGap_total];
        Int_t NumTrackEast[TriFlow::EtaGap_total], NumTrackWest[TriFlow::EtaGap_total];
        for(Int_t j = 0; j < TriFlow::EtaGap_total; j++)
        {
            Q2East[j].Set(-999.9,-999.9); // initialize Q Vector to unreasonable value
            Q2West[j].Set(-999.9,-999.9);
            Q3East[j].Set(-999.9,-999.9);
            Q3West[j].Set(-999.9,-999.9);
            NumTrackEast[j] = 0;
            NumTrackWest[j] = 0;
            if(mTriFlowCorrection->passTrackEtaNumCut(j))
            {
                // get QVector of sub event
                Q2East[j] = mTriFlowCorrection->getQVector(j,0,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q2West[j] = mTriFlowCorrection->getQVector(j,0,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q3East[j] = mTriFlowCorrection->getQVector(j,1,0); // 0 = eta_gap, 1 = flow type, 2 = east/west
                Q3West[j] = mTriFlowCorrection->getQVector(j,1,1); // 0 = eta_gap, 1 = flow type, 2 = east/west
                NumTrackEast[j] = mTriFlowCorrection->getNumTrack(j,0); // 0 = eta_gap, 1 = east/west
                NumTrackWest[j] = mTriFlowCorrection->getNumTrack(j,1); // 0 = eta_gap, 1 = east/west
            }
        }

        //if(mTriFlowCorrection->passTrackEtaNumCut(0))
        {
		//cout << "pass track eta Numb cut" << endl; 
            // Event Plane method
            Psi2_East = mTriFlowCorrection->calShiftAngle2East_EP(runIndex,cent9,vz_sign,0);
            Psi2_West = mTriFlowCorrection->calShiftAngle2West_EP(runIndex,cent9,vz_sign,0);
            Res2_EP = mTriFlowCorrection->getResolution2_EP(cent9,0);
            
	    // get N_prim, N_non_prim, N_Tof_match
            Int_t N_prim = nN_prim;
            Int_t N_non_prim = nN_non_prim;
            Int_t N_Tof_match = nMatchedToF;
            
	    // pass the event information to StTriFlowV0
            mTriFlowV0->clearEvent();
            mTriFlowV0->passEvent(N_prim, N_non_prim, N_Tof_match);

            // 2nd sub event plane
            mTriFlowV0->passEventPlane2East(Q2East[0],Q2East[1],Q2East[2],Q2East[3]);
            mTriFlowV0->passEventPlane2West(Q2West[0],Q2West[1],Q2West[2],Q2West[3]);

            // 3rd sub event plane
            mTriFlowV0->passEventPlane3East(Q3East[0],Q3East[1],Q3East[2],Q3East[3]);
            mTriFlowV0->passEventPlane3West(Q3West[0],Q3West[1],Q3West[2],Q3West[3]);
            // Number of Track in East and West part of TPC
            mTriFlowV0->passNumTrackEast(NumTrackEast[0],NumTrackEast[1],NumTrackEast[2],NumTrackEast[3]);
            mTriFlowV0->passNumTrackWest(NumTrackWest[0],NumTrackWest[1],NumTrackWest[2],NumTrackWest[3]);
	/*
	cout << "test 4" << endl;
	cout << "ME flag: " << TriFlow::Flag_ME << endl;
	cout << "cent9: "<< cent9 <<endl;
	cout << "vz: " << vz << endl;
	cout << "psi_east: " << Psi2_East << endl;
	cout << "Psi2_West: " << Psi2_West << endl;
	cout << "reweight: " << reweight << endl;
	cout << "Res2_EP: " << Res2_EP << endl;
    	*/    
    	mTriFlowV0->MixEvent_Phi(TriFlow::Flag_ME,mPicoDst,cent9,vz,Psi2_East,Psi2_West,reweight,Res2_EP);
            // cout<<"finish phi meson reconstruction"<<endl;
        }
        mTriFlowCorrection->clear();
    }
    return kStOK;
}

int StTriFlowMaker::GetRunIndex(int runID)
{
    int runIndex=-999;
    for(int i=0; i<139; i++)
    {
        if(runID==numbers[i])
        {
            runIndex=i;
        }
    }
    if(runIndex == -999) cout << "Run numbers are not found!!!" << endl;
    return runIndex;
}
int StTriFlowMaker::Centrality(int gRefMult )
{
    int centrality;
    int centFull[9]={6,12,22,39,64,100,154,191,241};
    if      (gRefMult>=centFull[8]) centrality=8;
    else if (gRefMult>=centFull[7]) centrality=7;
    else if (gRefMult>=centFull[6]) centrality=6;
    else if (gRefMult>=centFull[5]) centrality=5;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=3;
    else if (gRefMult>=centFull[2]) centrality=2;
    else if (gRefMult>=centFull[1]) centrality=1;
    else if (gRefMult>=centFull[0]) centrality=0;
    else centrality = 9;

    return centrality;
}
