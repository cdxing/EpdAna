#include "StTriFlowCorrection.h"
#include "StTriFlowConstants.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StMessMgr.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TNtuple.h"
#include "TRandom3.h"
#include "TF1.h"

Double_t Resolution_Full(Double_t *x_val, Double_t *par)
{
  Double_t y;
  Double_t chi = x_val[0];
  Double_t arg = chi*chi/4.0;
  Double_t norm = TMath::Sqrt(TMath::Pi()/2.0)/2.0;

  y = norm*chi*TMath::Exp(-1.0*arg)*(TMath::BesselI0(arg)+TMath::BesselI1(arg));

  return y;
}

ClassImp(StTriFlowCorrection)

TString StTriFlowCorrection::mVStr[2] = {"pos","neg"};
TString StTriFlowCorrection::mOrder[2] = {"2nd","3rd"};
TString StTriFlowCorrection::mMethod[2] = {"EP","SP"};
//---------------------------------------------------------------------------------

StTriFlowCorrection::StTriFlowCorrection(Int_t energy)
{
  mEnergy = energy;
}

//---------------------------------------------------------------------------------

StTriFlowCorrection::~StTriFlowCorrection()
{
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::InitReCenterCorrection(Int_t mEnergy)
{
    TString InPutFile = Form("/star/u/slan/pwg/fastoffline/7p7gev/recenterpar/file_%s_ReCenterPar.root",TriFlow::Energy[mEnergy].Data());

  mInPutFile = TFile::Open(InPutFile.Data());

  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::InitShiftCorrection(Int_t mEnergy)
{
    TString InPutFile_Shift = Form("/star/u/slan/pwg/fastoffline/7p7gev/recenter/production/out_Shift/shift.root");
    mInPutFile_Shift = TFile::Open(InPutFile_Shift.Data());

    //TString InPutFile_Res = Form("/star/u/slan/pwg/27GeVflow/resolution/out_Resolution/file_%s_Resolution.root",TriFlow::Energy[mEnergy].Data());
    TString InPutFile_Res = Form("/star/u/slan/pwg/fastoffline/7p7gev/resolution/out_Resolution/Resolution.root");
    mInPutFile_Res = TFile::Open(InPutFile_Res.Data());
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackEtaEast(StPicoTrack *track, Int_t i, Int_t Mode) // neg || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
    Float_t eta = track->pMom().PseudoRapidity();

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

bool StTriFlowCorrection::passTrackEtaWest(StPicoTrack *track, Int_t i, Int_t Mode) // pos || i = different eta_gap || Mode = 0 Event Plane Mode, Mode = 1 Flow Mode
{
  Float_t eta = track->pMom().PseudoRapidity();

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
//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackFull(StPicoTrack *track) // Full Event Plane pt Cut
{
  // pt cut 0.2 - 2.0 GeV/c
  Float_t pt = track->pMom().Perp();
  if(!(pt > TriFlow::mPrimPtMin[mEnergy] && pt < TriFlow::mPrimPtMax))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method and Scalar Product method
TVector2 StTriFlowCorrection::calq2Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q2Vector(0.0,0.0);

  const Float_t q2x = TMath::Cos(2.0*phi);
  const Float_t q2y = TMath::Sin(2.0*phi);
  q2Vector.Set(q2x,q2y);

  return q2Vector;
}

TVector2 StTriFlowCorrection::calq3Vector(StPicoTrack *track)
{
  const Float_t phi = track->pMom().Phi();
  TVector2 q3Vector(0.0,0.0);

  const Float_t q3x = TMath::Cos(3.0*phi);
  const Float_t q3y = TMath::Sin(3.0*phi);
  q3Vector.Set(q3x,q3y);

  return q3Vector;
}

Float_t StTriFlowCorrection::getWeight(StPicoTrack *track)
{
  Float_t pt = track->pMom().Perp();
  Float_t w;
  if(pt <= TriFlow::mPrimPtWeight)
  {
    w = pt;
  }
  if(pt > TriFlow::mPrimPtWeight)
  {
    w = TriFlow::mPrimPtWeight;
  }

  return w;
}
//------------------------------------------------------------------------------


TVector2 StTriFlowCorrection::getReCenterPar_East(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_East_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTriFlowCorrection::getReCenterPar_West(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t eta_gap, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_EtaGap_%d_West_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),eta_gap,mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

TVector2 StTriFlowCorrection::getReCenterPar_Full(Int_t order, Int_t Cent9, Int_t RunIndex, Int_t vz_sign, Int_t method)
{
  Float_t mean_qx, mean_qy;
  TVector2 qVector(0.0,0.0);

  TString ProName_x = Form("qx_%s_Vertex_%s_Full_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),mMethod[method].Data());
  TString ProName_y = Form("qy_%s_Vertex_%s_Full_%s",mOrder[order].Data(),mVStr[vz_sign].Data(),mMethod[method].Data());

  TProfile2D *p_x = (TProfile2D*)mInPutFile->Get(ProName_x.Data());
  TProfile2D *p_y = (TProfile2D*)mInPutFile->Get(ProName_y.Data());

  mean_qx = p_x->GetBinContent(p_x->FindBin((Double_t)RunIndex,(Double_t)Cent9));
  mean_qy = p_y->GetBinContent(p_y->FindBin((Double_t)RunIndex,(Double_t)Cent9));

  qVector.Set(mean_qx,mean_qy);

  return qVector;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_East(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_East_EP[j] += w*(calq2Vector(track) - getReCenterPar_East(0,Cent9,RunIndex,i,j,0));
  mQ3Vector_East_EP[j] += w*(calq3Vector(track) - getReCenterPar_East(1,Cent9,RunIndex,i,j,0));

  mQCounter_East[j]++;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_West(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i, Int_t j)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_West_EP[j] += w*(calq2Vector(track) - getReCenterPar_West(0,Cent9,RunIndex,i,j,0));
  mQ3Vector_West_EP[j] += w*(calq3Vector(track) - getReCenterPar_West(1,Cent9,RunIndex,i,j,0));

  mQCounter_West[j]++;
}

//---------------------------------------------------------------------------------

void StTriFlowCorrection::addTrack_Full(StPicoTrack *track, Int_t Cent9, Int_t RunIndex, Int_t i)
{
  // Event Plane method
  Float_t w = getWeight(track);
  mQ2Vector_Full_EP += w*(calq2Vector(track) - getReCenterPar_Full(0,Cent9,RunIndex,i,0));
  mQ3Vector_Full_EP += w*(calq3Vector(track) - getReCenterPar_Full(1,Cent9,RunIndex,i,0));

  mQCounter_Full++;

  Float_t eta = track->pMom().PseudoRapidity();
  if(eta >= 0.0)
  {
    mQCounter_Full_West++;
  }
  if(eta < 0.0)
  {
    mQCounter_Full_East++;
  }
}

Float_t StTriFlowCorrection::AngleShift(Float_t Psi_raw, Float_t order)
{
    Float_t Psi_Corr = Psi_raw;
    if(Psi_raw > TMath::Pi()/order)
    {
        Psi_Corr = Psi_raw - 2.0*TMath::Pi()/order;
    }
    if(Psi_raw < -1.0*TMath::Pi()/order)
    {
        Psi_Corr = Psi_raw + 2.0*TMath::Pi()/order;
    }

    return Psi_Corr;
}


Float_t StTriFlowCorrection::calShiftAngle2East_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
    Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_East_EP[eta_gap].Y(),mQ2Vector_East_EP[eta_gap].X())/2.0;
    Float_t mean_sin[5], mean_cos[5];
    Float_t delta_Psi = 0.0;
    Float_t Psi_Shift;

    for(Int_t k = 0; k < 5; k++) // Shift Order loop
    {
        TString ProName_cos, ProName_sin;
        TProfile2D *p_cos, *p_sin;

        ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
        p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
        mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

        ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_East_EP",mVStr[vz_sign].Data(),eta_gap,k);
        p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
        mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

        delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
    }

    Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
    Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

    return Psi_Shift;
}

Float_t StTriFlowCorrection::calShiftAngle2West_EP(Int_t runIndex, Int_t Cent9, Int_t vz_sign, Int_t eta_gap)
{
    Float_t Psi_ReCenter = TMath::ATan2(mQ2Vector_West_EP[eta_gap].Y(),mQ2Vector_West_EP[eta_gap].X())/2.0;
    Float_t mean_sin[5], mean_cos[5];
    Float_t delta_Psi = 0.0;
    Float_t Psi_Shift;

    for(Int_t k = 0; k < 5; k++) // Shift Order loop
    {
        TString ProName_cos, ProName_sin;
        TProfile2D *p_cos, *p_sin;

        ProName_cos = Form("CosPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
        p_cos = (TProfile2D*)mInPutFile_Shift->Get(ProName_cos.Data());
        mean_cos[k] = p_cos->GetBinContent(p_cos->FindBin((Double_t)runIndex,(Double_t)Cent9));

        ProName_sin = Form("SinPsi2_Vertex_%s_EtaGap_%d_Order_%d_West_EP",mVStr[vz_sign].Data(),eta_gap,k);
        p_sin = (TProfile2D*)mInPutFile_Shift->Get(ProName_sin.Data());
        mean_sin[k] = p_sin->GetBinContent(p_sin->FindBin((Double_t)runIndex,(Double_t)Cent9));

        delta_Psi += (1.0/2.0)*(2.0/(Float_t)(k+1))*(-1.0*mean_sin[k]*TMath::Cos(TriFlow::mShiftOrder2[k]*Psi_ReCenter)+mean_cos[k]*TMath::Sin(TriFlow::mShiftOrder2[k]*Psi_ReCenter));
    }

    Float_t Psi_Shift_raw = Psi_ReCenter + delta_Psi;
    Psi_Shift = AngleShift(Psi_Shift_raw,2.0);

    return Psi_Shift;
}
//Event Plane method
Float_t StTriFlowCorrection::getResolution2_EP(Int_t Cent9, Int_t eta_gap)
{
    TString ProName = Form("Res2_EtaGap_%d_EP",eta_gap);
    TProfile *p_res2 = (TProfile*)mInPutFile_Res->Get(ProName.Data());
    Float_t Res_raw = p_res2->GetBinContent(p_res2->FindBin(Cent9));
    if(Res_raw <= 0)
    {
        return -999.9;
    }
    else
    {
        Float_t Res = TMath::Sqrt(Res_raw);
        return Res;
    }
}
//---------------------------------------------------------------------------------

void StTriFlowCorrection::clear()
{
  for(Int_t i = 0; i < 4; i++)
  {
    mQ2Vector_East_EP[i].Set(0.0,0.0);
    mQ3Vector_East_EP[i].Set(0.0,0.0);
    mQCounter_East[i] = 0;

    mQ2Vector_West_EP[i].Set(0.0,0.0);
    mQ3Vector_West_EP[i].Set(0.0,0.0);
    mQCounter_West[i] = 0;
  }
}

//---------------------------------------------------------------------------------

bool StTriFlowCorrection::passTrackEtaNumCut(Int_t j)
{
  if(!(mQCounter_East[j] > TriFlow::mTrackMin && mQCounter_West[j] > TriFlow::mTrackMin))
  {
    return kFALSE;
  }

  return kTRUE;
}

//---------------------------------------------------------------------------------
// Event Plane method
// 2nd
//---------------------------------------------------------------------------------
TVector2 StTriFlowCorrection::getQVector(Int_t j, Int_t k, Int_t l) // 0 = eta_gap, 1 = flow type, 2 = east/west
{
  if(k == 0 && l == 0) return mQ2Vector_East_EP[j];
  if(k == 0 && l == 1) return mQ2Vector_West_EP[j];
  if(k == 1 && l == 0) return mQ3Vector_East_EP[j];
  if(k == 1 && l == 1) return mQ3Vector_West_EP[j];
}

Int_t StTriFlowCorrection::getNumTrack(Int_t j, Int_t l) // 0 = eta_gap, 1 = east/west
{
  if(l == 0) return mQCounter_East[j];
  if(l == 1) return mQCounter_West[j];
}
//---------------------------------------------------------------------------------
