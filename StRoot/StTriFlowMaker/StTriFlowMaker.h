#ifndef StTriFlowMaker_h
#define StTriFlowMaker_h

#include "StMaker.h"
#include "TString.h"
class StPicoDst;
class StPicoDstMaker;
class StPicoEvent;
class StTriFlowCut;
class StTriFlowProManger;
class StTriFlowCorrection;
class StTriFlowHistoManger;
class StTriFlowV0;
class StRefMultCorr;

// run QA
class TProfile;
class TH1F;
class TH2F;

class StTriFlowMaker : public StMaker {
    public:
        StTriFlowMaker(const char *name, StPicoDstMaker *picoMaker, char* jobid);
        virtual ~StTriFlowMaker();

        virtual Int_t Init();
        virtual Int_t Make();
        virtual void  Clear(Option_t *opt="");
        virtual Int_t Finish();

        int Centrality(int gRefMult);
        int GetRunIndex(int runId);
        bool isGoodTrigger(StPicoEvent const*) const;

    private:
        StPicoDstMaker *mPicoDstMaker;
        StPicoDst      *mPicoDst;
        StPicoEvent *mPicoEvent;
        static StRefMultCorr *mRefMultCorr;
        StTriFlowCut *mTriFlowCut;
        StTriFlowProManger *mTriFlowProManger;
        StTriFlowCorrection *mTriFlowCorrection;
        StTriFlowHistoManger *mTriFlowHistoManger;
        StTriFlowV0 *mTriFlowV0;


        Int_t mEnergy;
        TString mOutPut_Phi;
        TString mOutPut_Lambda;
        TString mOutPut_AntiLambda;
        TFile *mFile_Phi;

        ClassDef(StTriFlowMaker, 1)
};

#endif
