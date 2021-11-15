#include "StTriFlowCut.h"
#include "StTriFlowConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"  // shaowei  18c
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StMessMgr.h"

ClassImp(StTriFlowCut)

    //StRefMultCorr* StTriFlowCut::mRefMultCorr = NULL; // shaowei
    //---------------------------------------------------------------------------------

StTriFlowCut::StTriFlowCut(Int_t energy)
{
    mEnergy = energy;
}

//---------------------------------------------------------------------------------

StTriFlowCut::~StTriFlowCut()
{
}

//---------------------------------------------------------------------------------

bool StTriFlowCut::isGoodTrigger(StPicoDst *pico)
{
    array<unsigned int, 1> const triggers = {
        630052
    };
    StPicoEvent *event = pico->event();
    for(unsigned i=0; i<event->triggerIds().size(); i++)
    {
      cout << "trigger ID = " << event->triggerIds()[i] << endl;
    }
    for(auto trg: triggers)
    {
        if(event->isTrigger(trg)) return kTRUE;
    }
}
bool StTriFlowCut::passEventCut(StPicoDst *pico)
{
    // initialize mMatchedToF
    //if(!isGoodTrigger(pico)) return kFALSE;
    StPicoEvent *event = pico->event();
    if(!event)
    {
        return kFALSE;
    }

    // initialize StRefMultCorr
    const Int_t runId = event->runId();
    const Int_t refMult = event->refMult();
    const Int_t tofMult = event->btofTrayMultiplicity();
    const Int_t tofmatch= event->nBTOFMatch();

    //if(tofMult < (2.8  * refMult - 150)) return 0;
    //if(tofMult > (3.45 * refMult + 300)) return 0;

    //if(tofmatch < (0.88 * refMult - 40)) return 0;
    //if(tofmatch > (1.25 * refMult + 50)) return 0;

    const Float_t vx = event->primaryVertex().X();
    const Float_t vy = event->primaryVertex().Y();
    const Float_t vz = event->primaryVertex().Z();
    /*
    if( (vx < 1.e-5 && vx > -1.e-5) &&
            (vy < 1.e-5 && vy > -1.e-5) &&
            (vz < 1.e-5 && vz > -1.e-5)  )
        return kFALSE;
    */
    // vz cut
    if(fabs(vz-200.0) > 2.0)
    {
        return kFALSE;
    }
    if(event->btofTrayMultiplicity()<2)return kFALSE;
    // vr cut
    if(sqrt(vx*vx+(vy+2.0)*(vy+2.0)) > 2.0)
    {
        return kFALSE;
    }
    return kTRUE;
}

//---------------------------------------------------------------------------------

bool StTriFlowCut::passTrackBasic(StPicoTrack *track)
{
    // nHitsFit cut
    if(track->nHitsFit() < TriFlow::mHitsFitTPCMin)
    {
        return kFALSE;
    }

    // nHitsRatio cut
    if(track->nHitsMax() <= TriFlow::mHitsMaxTPCMin)
    {
        return kFALSE;
    }
    if((Float_t)track->nHitsFit()/(Float_t)track->nHitsMax() < TriFlow::mHitsRatioTPCMin)
    {
        return kFALSE;
    }

    // eta cut
    Float_t eta = track->pMom().PseudoRapidity();
    /*
    if(fabs(eta) > TriFlow::mEtaMax)
    {
        return kFALSE;
    }
    */
    return kTRUE;
}

bool StTriFlowCut::passTrackEP(StPicoTrack *track, float dca)
{
    if(!track) return kFALSE;

    if(!passTrackBasic(track)) return kFALSE;

    // dca cut
    // dca cut for event plane reconstruction: 200GeV = 3.0, BES = 1.0
    //if(dca > TriFlow::mDcaEPMax[mEnergy])   // change by shaowei
    //cout << "dcaaaaaaaa = "<< dca << endl;
    if(dca > 3.0)   // change by shaowei
    {
        return kFALSE;
    }
    //cout << "dca = "<< dca << endl;

    // pt cut for event plane
    // pt cut 0.2 - 2.0 GeV/c
    Float_t pt = track->pMom().Perp();
    Float_t p  = track->pMom().Mag();
    if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax && p < TriFlow::mPrimMomMax))
    {
        return kFALSE;
    }

    return kTRUE;
}
bool StTriFlowCut::passTrackPhi(StPicoTrack *track, float dca)
{
    if(!track) return kFALSE;

    if(!passTrackBasic(track)) return kFALSE;

    // dca cut for flow analysis: 2.0
    if(dca > TriFlow::mDcaTrMax_phi)
    {
        return kFALSE;
    }

    // primary pt and momentum cut: PtMin = 0.1
    if(!(track->pMom().Perp() > TriFlow::mGlobPtMin && track->pMom().Mag() < TriFlow::mPrimMomMax))
    {
        return kFALSE;
    }

    return kTRUE;
}

// Dip angle cut
bool StTriFlowCut::passDipAngle(Double_t dipAngle)
{
  if(dipAngle <= 0.04)
  {
    return kFALSE;
  }

  return kTRUE;
}
Float_t StTriFlowCut::getMass2(StPicoTrack *track, StPicoDst *pico)
{
    Float_t Mass2 = -999.0;
    int tofindex=track->bTofPidTraitsIndex();
    if(tofindex >= 0)
    {
        StPicoBTofPidTraits *tof= pico->btofPidTraits(tofindex);
        Float_t Beta = tof->btofBeta();
        //cout<<"Beta = "<<Beta<<endl;
        Float_t Momentum = track->pMom().Mag(); // primary momentum

        if(tof->btofMatchFlag() > 0 && tof->btof() != 0 && Beta != 0)
        {
            Mass2 = Momentum*Momentum*(1.0/(Beta*Beta) - 1.0);
        }
    }
    return Mass2;
}

bool StTriFlowCut::passSigKaonCut(StPicoTrack* track, Float_t scale_nSigma_factor)
{
    Float_t nSigmaKaon = track->nSigmaKaon();
    if(fabs(nSigmaKaon*scale_nSigma_factor) > TriFlow::mNSigmaKaonMax)
    {
        return kFALSE;
    }
    return kTRUE;
}
