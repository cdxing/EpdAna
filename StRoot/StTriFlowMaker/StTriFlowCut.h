#ifndef StTriFlowCut_h
#define StTriFlowCut_h

#include "TObject.h"
#include "TString.h"
//#include <array>

class StPicoDst;
class StPicoTrack;
class StPicoBTofPidTraits;

using namespace std;


class StTriFlowCut : public TObject
{
  public:
    StTriFlowCut(Int_t energy);
    ~StTriFlowCut();

    bool isGoodTrigger(StPicoDst* );
    bool passEventCut(StPicoDst*);
    bool passTrackBasic(StPicoTrack*);
    bool passTrackEP(StPicoTrack*, float);
    bool passTrackPhi(StPicoTrack*, float);
    bool passDipAngle(Double_t);
    bool passSigKaonCut(StPicoTrack*, Float_t);
    Float_t getMass2(StPicoTrack*, StPicoDst*);

  private:
    Int_t mEnergy;


    ClassDef(StTriFlowCut,1)
};
#endif
