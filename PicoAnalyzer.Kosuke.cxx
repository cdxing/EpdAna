#include <TFile.h>
#include <TTree.h>
#include <StMessMgr.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TROOT.h>
#include <TRandom.h>
#include <StThreeVectorF.hh>
#include <StHelix.hh>
#include <TLorentzVector.h>

//#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoEvent/StPicoEvent.h"
#include "StPicoEvent/StPicoDst.h"
#include "StPicoEvent/StPicoTrack.h"
#include "StPicoEvent/StPicoBTofPidTraits.h"
#include "StPicoEvent/StPicoEpdHit.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StarClassLibrary/StPhysicalHelixD.hh"
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StEpdUtil/StEpdGeom.h"
#include "StEpdUtil/StEpdEpFinder.h"
#include "StThreeVectorF.hh"

#include "../RunNumber.h"
#include "PicoAnalyzer.h"

const float C_C_LIGHT = 299792458;//(m/s)

ClassImp(PicoAnalyzer)

		//__________________________________________________________________________________
		PicoAnalyzer::PicoAnalyzer( const Char_t *name, StPicoDstMaker *picoMaker, const Char_t *outName ) : StMaker(name) { 
				mPicoDstMaker = picoMaker;
				mPicoDst = 0;
				mEventCounter = 0; 

				mQvectorFileName        = Form("Qvector_%s.root",outName);
				mCalibrationFileName    = Form("Calibration_%s.root",outName);
				mLambdaAnalisysFileName = Form("LambdaAnalisys_%s.root",outName);
		}

//__________________________________________________________________________________
Int_t PicoAnalyzer::Init() {
		cout << "Init" << endl;
		cout << "runnumber = " << runnumber << endl;
		mEpdGeom = new StEpdGeom();

		//===============================                          
		//  Define Histograms                           
		//===============================                                            

		hVz       = new TH1F("Vztpc","Vztpc",120,130,250);
		hVzAfter  = new TH1F("VztpcAfter","Vztpc after event cut",120,130,250);
		hTPCmultiplicity                              = new TH1F("TPCmultiplicity","TPCmultiplicity;TPCmultiplicity",450,0,450);
		hTPCmultiplicityAfterPileup                   = new TH1F("TPCmultiplicityAfterPileup","TPCmultiplicity after pile up;TPCmultiplicity",450,0,450);
		hMultiplicityvsTofmult                        = new TH2F("MultiplicityvsTofmult","MultiplicityvsTofmult;Multiplicity;Tofmult",100,0,350,100,0,600);
		hMultiplicityvsTofMatch                       = new TH2F("MultiplicityvsTofmatch","MultiplicityvsTofmatch;Multiplicity;Tofmatch",100,0,350,100,0,200);
		hNumOfPionProtonvsMultiplicity                = new TH2F("NumOfPionProtonvsMultiplicity","Multiplicity vs # of Pion&Proton;# of Pion&Proton;TPC multiplicity",450,0,450,200,0,200);
		hNumOfPrimaryvsNumOfGlobal                    = new TH2F("NumOfPrimaryvsNumOfGlobal","# of primary vs # of global",100,0,350,100,0,1000);
		hMultiplicityvsTofmultAfterPileup             = new TH2F("MultiplicityvsTofmultAfterPileup","MultiplicityvsTofmult after pile up;Multiplicity;Tofmult",100,0,350,100,0,600);
		hMultiplicityvsTofMatchAfterPileup            = new TH2F("MultiplicityvsTofmatchAfterPileup","MultiplicityvsTofmatch after pile up;Multiplicity;Tofmatch",100,0,350,100,0,200);
		hNumOfPrimaryvsNumOfGlobalAfterPileup         = new TH2F("NumOfPrimaryvsNumOfGlobalAfterPileup","# of primary vs # of global",100,0,350,100,0,1000);
		hNumOfPionProtonvsMultiplicityAfterPileup     = new TH2F("NumOfPionProtonvsMultiplicityAfterPileup","Multiplicity vs # of Pion&Proton;# of Pion&Proton;TPC multiplicity",450,0,450,200,0,200);
		hCountCentrality                              = new TH1F("CountCentrality","CountCentrality",9,-1,8);

		//for re-centering && flattening
		for(int rapidity=0;rapidity<3;rapidity++){
				pTpcQv[rapidity][0]      = new TProfile(Form("TpcQx_rapidity%d",rapidity),Form("TpcQx_rapidity%d",rapidity),8,0,8,"s");
				pTpcQv[rapidity][1]      = new TProfile(Form("TpcQy_rapidity%d",rapidity),Form("TpcQy_rapidity%d",rapidity),8,0,8,"s");
				pTPC_shift_sin[rapidity] = new TProfile2D(Form("TPC_shiftsin_rapidity%d",rapidity),Form("TPC_shiftsin_rapidity%d",rapidity),8,-0.5,7.5,8,-0.5,7.5);
				pTPC_shift_cos[rapidity] = new TProfile2D(Form("TPC_shiftcos_rapidity%d",rapidity),Form("TPC_shiftcos_rapidity%d",rapidity),8,-0.5,7.5,8,-0.5,7.5);
		}
		for(int nmip=0;nmip<4;nmip++){
				for(int iring=0;iring<4;iring++){
						pEpdQv[nmip][iring][0]      = new TProfile(Form("EpdQx_nmipmax%d_ring%d",nmip+2,iring),Form("EpdQx_nmipmax%d_ring%d",nmip+2,iring),8,0,8,"s");
						pEpdQv[nmip][iring][1]      = new TProfile(Form("EpdQy_nmipmax%d_ring%d",nmip+2,iring),Form("EpdQy_nmipmax%d_ring%d",nmip+2,iring),8,0,8,"s");
						pEPD_shift_sin[nmip][iring] = new TProfile2D(Form("EPD_shiftsin_nmipmax%d_iring%d",nmip+2,iring),Form("EPD_shiftsin_nmipmax%d_iring%d",nmip+2,iring),8,-0.5,7.5,8,-0.5,7.5);
						pEPD_shift_cos[nmip][iring] = new TProfile2D(Form("EPD_shiftcos_nmipmax%d_iring%d",nmip+2,iring),Form("EPD_shiftcos_nmipmax%d_iring%d",nmip+2,iring),8,-0.5,7.5,8,-0.5,7.5);
				}
		}

		//Calculate v1
		for(int psi=0;psi<7;psi++){
				for(int cent=0;cent<8;cent++){
						for(int pt=0;pt<7;pt++){
								for(int particle=0;particle<2;particle++){
										CalcV1PionTPC[psi][cent][pt][particle]   = new TProfile(Form("CalcV1TPC_Psi%dCent%dPt%dPion%d",psi,cent,pt,particle),Form("CalcV1TPC_Psi%dCent%dPt%dPion%d",psi,cent,pt,particle),6,-1.5,0);
										CalcV1ProtonTPC[psi][cent][pt][particle] = new TProfile(Form("CalcV1TPC_Psi%dCent%dPt%dProton%d",psi,cent,pt,particle),Form("CalcV1TPC_Psi%dCent%dPt%dProton%d",psi,cent,pt,particle),6,-1.5,0);
								}
						}
						CalcV1EPD[psi][cent] = new TProfile(Form("CalcV1EPD_Psi%dCent%d",psi,cent),Form("CalcV1EPD_Psi%dCent%d",psi,cent),16,-5.53,-2.52);
				}
		}



		const int EPcorr_epdA[6]   = {0,0,0,1,1,2};
		const int EPcorr_epdB[6]   = {1,2,3,2,3,3};
		const int EPcorr_tpcA[3]   = {0,0,1};
		const int EPcorr_tpcB[3]   = {1,2,2};
		const int EPcorr_evstA[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
		const int EPcorr_evstB[12] = {0,1,2,0,1,2,0,1,2,0,1,2};
		//for calculate event plane
		for(int cent=0;cent<8;cent++){
				for(int etaId=0;etaId<3;etaId++){
						hQvTpcRaw[cent][etaId]   = new TH2F(Form("QxvsQytpcRaw_cent%d_etaID%d",cent,etaId),Form("QxvsQytpcRaw_cent%d_etaID%d",cent,etaId),50,-10,10,50,-10,10);
						hQvTpcRece[cent][etaId]  = new TH2F(Form("QxvsQytpcRece_cent%d_etaID%d",cent,etaId),Form("QxvsQytpcRece_cent%d_etaID%d",cent,etaId),50,-10,10,50,-10,10);
						hPsiTpcRaw[cent][etaId]  = new TH1F(Form("PsiTpcRaw_cent%d_etaID%d",cent,etaId),Form("PsiTpcRaw_cent%d_etaId%d",cent,etaId),300,-TMath::Pi(),TMath::Pi());
						hPsiTpcRece[cent][etaId] = new TH1F(Form("PsiTpcRece_cent%d_etaID%d",cent,etaId),Form("PsiTpcRece_cent%d_etaId%d",cent,etaId),300,-TMath::Pi(),TMath::Pi());
						hPsiTpcFlat[cent][etaId] = new TH1F(Form("PsiTpcFlat_cent%d_etaID%d",cent,etaId),Form("PsiTpcFlat_cent%d_etaId%d",cent,etaId),300,-TMath::Pi(),TMath::Pi());
				}

				for(int nmip=0;nmip<4;nmip++){
						for(int ringId=0;ringId<4;ringId++){
								hQvEpdRaw[cent][nmip][ringId]   = new TH2F(Form("QxvsQyepdRaw_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),Form("QxvsQyepdRaw_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),50,-10,10,50,-10,10);
								hQvEpdRece[cent][nmip][ringId]  = new TH2F(Form("QxvsQyepdRece_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),Form("QxvsQyepdRece_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),50,-10,10,50,-10,10);
								hPsiEpdRaw[cent][nmip][ringId]  = new TH1F(Form("PsiEpdRaw_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),Form("PsiEpdRaw_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),300,-TMath::Pi(),TMath::Pi());
								hPsiEpdRece[cent][nmip][ringId] = new TH1F(Form("PsiEpdRece_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),Form("PsiEpdRece_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),300,-TMath::Pi(),TMath::Pi());
								hPsiEpdFlat[cent][nmip][ringId] = new TH1F(Form("PsiEpdFlat_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),Form("PsiEpdFlat_cent%d_nMipMax%d_ringID%d",cent,nmip,ringId),300,-TMath::Pi(),TMath::Pi());
						}
				}
		}
		for(int epcorr=0;epcorr<12;epcorr++){
				if(epcorr<3){
						pEPcorr_TPC[epcorr] = new TProfile(Form("EPcorr_TPC%dvsTPC%d",EPcorr_tpcA[epcorr],EPcorr_tpcB[epcorr]),Form("EPcorr_TPC%dvsTPC%d",EPcorr_tpcA[epcorr],EPcorr_tpcB[epcorr]),8,0,8);
						pEPcorrsin_TPC[epcorr] = new TProfile(Form("EPcorrSin_TPC%dvsTPC%d",EPcorr_tpcA[epcorr],EPcorr_tpcB[epcorr]),Form("EPcorrSin_TPC%dvsTPC%d",EPcorr_tpcA[epcorr],EPcorr_tpcB[epcorr]),8,0,8);
				}
				for(int nmip=0;nmip<4;nmip++){
						if(epcorr<6){
								pEPcorr_EPD[epcorr][nmip] = new TProfile(Form("EPcorr_EPD%dvsEPD%d_nMipMax%d",EPcorr_epdA[epcorr],EPcorr_epdB[epcorr],nmip+2),Form("EPcorr_EPD%dvsEPD%d_nMipMax%d",EPcorr_epdA[epcorr],EPcorr_epdB[epcorr],nmip+2),8,0,8);
								pEPcorrsin_EPD[epcorr][nmip] = new TProfile(Form("EPcorrSin_EPD%dvsEPD%d_nMipMax%d",EPcorr_epdA[epcorr],EPcorr_epdB[epcorr],nmip+2),Form("EPcorrSin_EPD%dvsEPD%d_nMipMax%d",EPcorr_epdA[epcorr],EPcorr_epdB[epcorr],nmip+2),8,0,8);
						}
						pEPcorr_EPDvsTPC[epcorr][nmip] = new TProfile(Form("EPcorr_EPD%dvsTPC%d_nMipMax%d",EPcorr_evstA[epcorr],EPcorr_evstB[epcorr],nmip+2),Form("EPcorr_EPD%dvsTPC%d_nMipMax%d",EPcorr_evstA[epcorr],EPcorr_evstB[epcorr],nmip+2),8,0,8);
						pEPcorrsin_EPDvsTPC[epcorr][nmip] = new TProfile(Form("EPcorrSin_EPD%dvsTPC%d_nMipMax%d",EPcorr_evstA[epcorr],EPcorr_evstB[epcorr],nmip+2),Form("EPcorrSin_EPD%dvsTPC%d_nMipMax%d",EPcorr_evstA[epcorr],EPcorr_evstB[epcorr],nmip+2),8,0,8);
				}
		}
		pV1RawFullRegion = new TProfile2D("V1RawFullRegion","V1FullRegion",70,-6.5,0.5,8,0,80);


		hTPCphivsTOFphi             = new TH2F("TPCphivsTOFphi","TPCphivsTOFphi;#phi^{TPC};#phi^{TOF}",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
		hTPCphivsTOFphiDproton      = new TH2F("TPCphivsTOFphiDproton","TPCphivsTOFphiDproton;#phi^{TPC};#phi^{TOF}",100,-TMath::Pi(),TMath::Pi(),100,-TMath::Pi(),TMath::Pi());
		hDeltaphiTpcTof             = new TH1F("DeltaPhiTpcTof","DeltaPhiTpcTof;#phi^{TOF}-#phi^{TPC}",300,-2*TMath::Pi(),2*TMath::Pi());
		hTofXvsY                    = new TH2F("TofXvsY","TofXvsY;TofX;TofY",100,-250,250,100,-250,250);
		hTofXvsZ                    = new TH2F("TofXvsZ","TofXvsZ;TofX;TofZ",100,-250,250,100,-250,250);
		hTofYvsZ                    = new TH2F("TofYvsZ","TofYvsZ;TofY;TofZ",100,-250,250,100,-250,250);
		hDeltaPhi2D                 = new TH2F("DeltaPhi2D","DeltaPhi2D;#phi_{#Lambda}-#phi*_{p};#phi^{TOF}-#phi^{TPC}",300,-2*TMath::Pi(),2*TMath::Pi(),300,-2*TMath::Pi(),2*TMath::Pi());
		for(int dphi=0;dphi<2;dphi++){
				hDeltaPhiPionvsProton[dphi] = new TH2F(Form("DeltaPhiPionvsProton_dPhi%d",dphi),Form("DeltaPhiPionvsProton_dPhi%d",dphi),300,-2*TMath::Pi(),2*TMath::Pi(),300,-2*TMath::Pi(),2*TMath::Pi());
		}
		hDeltaLambdaProtonPhi_LaboV0RF = new TH2F("DeltaLambdaProtonPhi_LaboV0RF","DeltaLambdaProtonPhi_LaboV0RF",300,-TMath::Pi(),TMath::Pi(),300,-TMath::Pi(),TMath::Pi());
		hDeltaLambdaProtonPhivsPionPhi = new TH2F("DeltaLambdaProtonPhivsPionPhi","DeltaLambdaProtonPhivsPionPhi",300,-TMath::Pi(),TMath::Pi(),300,-TMath::Pi(),TMath::Pi());
		for(int charge=0;charge<2;charge++){
				hDeltaPhivsPt[charge] = new TH2F(Form("DeltaPhivsPt_Charge%d",charge),Form("DeltaPhivsPt_Charge%d",charge),300,-TMath::Pi(),TMath::Pi(),300,0,5);
				hDeltaEtavsPt[charge] = new TH2F(Form("DeltaEtavsPt_Charge%d",charge),Form("DeltaEtavsPt_Charge%d",charge),300,-2,2,300,0,5);
		}

		//PID plot
		hDedx                       = new TH2F("Dedx","p/q vs de/dx;p/q;de/dx",100,-3,3,150,0,15);
		hDedx_pion                  = new TH2F("Dedx_pion","p/q vs de/dx pion;p/q;de/dx of pion",100,-3,3,150,0,15);
		hDedx_proton                = new TH2F("Dedx_proton","p/q vs de/dx proton;p/q;de/dx of proton",100,-3,3,150,0,15);
		hMM                         = new TH2F("MM","p/q vs m^{2};p/q;m^{2}",100,-3,3,100,-0.5,2.5);
		hMM_pion                    = new TH2F("MM_pion","p/q vs m^{2} pion;p/q;m^{2} of pion",100,-3,3,100,-0.5,2.5);
		hMM_proton                  = new TH2F("MM_proton","p/q vs m^{2} proton;p/q;m^{2} of proton",100,-3,3,100,-0.5,2.5);

		hprimaryDedx_pion           = new TH2F("PrimaryDedx_pion","p/q vs de/dx Primary pion;p/q;de/dx of Primary pion",100,-3,3,150,0,15);
		hprimaryDedx_proton         = new TH2F("PrimaryDedx_proton","p/q vs de/dx Primary proton;p/q;de/dx of Primary proton",100,-3,3,150,0,15);
		hprimaryMM_pion             = new TH2F("PrimaryMM_pion","p/q vs m^{2} Primary pion;p/q;m^{2} of Primary pion",100,-3,3,100,-0.5,2.5);
		hprimaryMM_proton           = new TH2F("PrimaryMM_proton","p/q vs m^{2} Primary proton;p/q;m^{2} of Primary proton",100,-3,3,100,-0.5,2.5);

		hLambdaEta[0]                  = new TH1F("LambdaEta","LambadEta",250,-4,1);
		hLambdaEta[1]                  = new TH1F("antiLambdaEta","antiLambadEta",250,-4,1);
		hLambdaPt[0]                   = new TH1F("LambdaPt","LambdaPt",240,0,6);
		hLambdaPt[1]                   = new TH1F("antiLambdaPt","antiLambdaPt",240,0,6);
		hLambdaDeltaPhi[0]             = new TH1F("LambdaDeltaPhi","LambdaDeltaPhi",300,-TMath::Pi(),TMath::Pi());
		hLambdaDeltaPhi[1]             = new TH1F("antiLambdaDeltaPhi","antiLambdaDeltaPhi",300,-TMath::Pi(),TMath::Pi());
		hLambdaPhi[0]                  = new TH1F("LambdaPhi","LambdaPhi",300,-TMath::Pi(),TMath::Pi());
		hLambdaPhi[1]                  = new TH1F("antiLambdaPhi","antiLambdaPhi",300,-TMath::Pi(),TMath::Pi());
		hLambdaRapidity[0]             = new TH1F("LambdaRapidity","LambadRapidity",250,-4,1);
		hLambdaRapidity[1]             = new TH1F("antiLambdaRapdity","antiLambadRapidity",250,-4,1);
		hLambdaRapidityvsPt[0]         = new TH2F("LambdaRapidityvsPt","LambadRapidityvsPt",300,-2.5,0.5,300,0,3.0);
		hLambdaRapidityvsPt[1]         = new TH2F("antiLambdaRapidityvsPt","antiLambadRapidityvsPt",300,-2.5,0.5,300,0,3.0);
		hinvMLambda_topolo[0][0]       = new TH1F("invMLambda_Raw","invM Lambda Raw",340,1.05,1.18);
		hinvMLambda_topolo[0][1]       = new TH1F("invMLambda_ppiDCA","invM Lambda ppiDCA",340,1.05,1.18);
		hinvMLambda_topolo[0][2]       = new TH1F("invMLambda_DaughterDCA","invM Lambda DaughterDCA",340,1.05,1.18);
		hinvMLambda_topolo[0][3]       = new TH1F("invMLambda_LambdaDCA","invM Lambda LambdaDCA",340,1.05,1.18);
		hinvMLambda_topolo[0][4]       = new TH1F("invMLambda_DecayLength","invM Lambda DecayLength",340,1.05,1.18);
		hinvMLambda_topolo[1][0]       = new TH1F("invMantiLambda_Raw","invM antiLambda Raw",340,1.05,1.18);
		hinvMLambda_topolo[1][1]       = new TH1F("invMantiLambda_ppiDCA","invM antiLambda ppiDCA",340,1.05,1.18);
		hinvMLambda_topolo[1][2]       = new TH1F("invMantiLambda_DaughterDCA","invM antiLambda DaughterDCA",340,1.05,1.18);
		hinvMLambda_topolo[1][3]       = new TH1F("invMantiLambda_LambdaDCA","invM antiLambda LambdaDCA",340,1.05,1.18);
		hinvMLambda_topolo[1][4]       = new TH1F("invMantiLambda_DecayLength","invM antiLambda DecayLength",340,1.05,1.18);

		//Change Topological cut Now using EPD0 event plane
		for(int Lambda=0;Lambda<2;Lambda++){
				hDaughterPionRapidityvsPt[Lambda] = new TH2F(Form("DaughterPionRapidityvsPt%d",Lambda),Form("DaughterPionRapidityvsPt%d;y;p_{T}",Lambda),200,-3,3,200,0,5);
				hDaughterProtonRapidityvsPt[Lambda] = new TH2F(Form("DaughterProtonRapidityvsPt%d",Lambda),Form("DaughterProtonRapidityvsPt%d;y;p_{T}",Lambda),200,-3,3,200,0,5);
				hDaughterPionRapidityvsPtRF[Lambda] = new TH2F(Form("DaughterPionRapidityvsPtRF%d",Lambda),Form("DaughterPionRapidityvsPtRF%d;y;p_{T}",Lambda),200,-3,3,200,0,5);
				hDaughterProtonRapidityvsPtRF[Lambda] = new TH2F(Form("DaughterProtonRapidityvsPtRF%d",Lambda),Form("DaughterProtonRapidityvsPtRF%d;y;p_{T}",Lambda),200,-3,3,200,0,5);
				hDaughterPionDeltaPhi[Lambda] = new TH1F(Form("DaughterPionDeltaPhi%d",Lambda),Form("DaughterPionDeltaPhi%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				hDaughterProtonDeltaPhi[Lambda] = new TH1F(Form("DaughterProtonDeltaPhi%d",Lambda),Form("DaughterProtonDeltaPhi%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				hDaughterPionPhi[Lambda] = new TH1F(Form("DaughterPionPhi%d",Lambda),Form("DaughterPionPhi%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				hDaughterProtonPhi[Lambda] = new TH1F(Form("DaughterProtonPhi%d",Lambda),Form("DaughterProtonPhi%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				hDaughterPionPhiRF[Lambda] = new TH1F(Form("DaughterPionPhiRF%d",Lambda),Form("DaughterPionPhiRF%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				hDaughterProtonPhiRF[Lambda] = new TH1F(Form("DaughterProtonPhiRF%d",Lambda),Form("DaughterProtonPhiRF%d",Lambda),300,-TMath::Pi(),TMath::Pi());
				for(int mass=0;mass<7;mass++){
						RapidityvsDeltaphi[Lambda][mass]             = new TH2F(Form("RapidityvsDeltaphi_Lambda%dMass%d",Lambda,mass),Form("RapidityvsDeltaphi_Lambda%dMass%d;y_{#Lambda};#phi-#Psi_{1}",Lambda,mass),15,-1.5,0,16,-TMath::Pi(),TMath::Pi());
						if(Lambda==0)RapidityvsDeltaphiAfterCorrect[mass] = new TH2F(Form("RapidityvsDeltaphiAfterCorrect_Lambda%dMass%d",Lambda,mass),Form("RapidityvsDeltaphi_Lambda%dMass%d;y_{#Lambda};#phi-#Psi_{1}",Lambda,mass),15,-1.5,0,16,-TMath::Pi(),TMath::Pi());
				}
				for(int cent=0;cent<8;cent++){
						for(int topolo=0;topolo<10;topolo++){
								hinvMLambdaTopolo[Lambda][cent][topolo]    = new TH1F(Form("invMLambda%dCent%d_pDCA%d",Lambda,cent,topolo),Form("invMLambda%dCent%d_pDCA%d",Lambda,cent,topolo),340,1.05,1.18);
								hinvMLambdaTopolo[Lambda][cent][topolo+10] = new TH1F(Form("invMLambda%dCent%d_piDCA%d",Lambda,cent,topolo),Form("invMLambda%dCent%d_piDCA%d",Lambda,cent,topolo),340,1.05,1.18);
								hinvMLambdaTopolo[Lambda][cent][topolo+20] = new TH1F(Form("invMLambda%dCent%d_ppiDCA%d",Lambda,cent,topolo),Form("invMLambda%dCent%d_ppiDCA%d",Lambda,cent,topolo),340,1.05,1.18);
								hinvMLambdaTopolo[Lambda][cent][topolo+30] = new TH1F(Form("invMLambda%dCent%d_LambdaDCA%d",Lambda,cent,topolo),Form("invMLambda%dCent%d_LambdaDCA%d",Lambda,cent,topolo),340,1.05,1.18);
								hinvMLambdaTopolo[Lambda][cent][topolo+40] = new TH1F(Form("invMLambda%dCent%d_DecayL%d",Lambda,cent,topolo),Form("invMLambda%dCent%d_DecayL%d",Lambda,cent,topolo),340,1.05,1.18);
						}
				}
		}

		//InvM vs Average quantities
		for(int Lambda=0;Lambda<2;Lambda++){
				pInvMvsDaughterPionDeltaPhi[Lambda]         = new TProfile(Form("invMvsDaughterPionDeltaPhi%d",Lambda),Form("invMvsDaughterPionDeltaPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterPionPhi[Lambda]              = new TProfile(Form("invMvsDaughterPionPhi%d",Lambda),Form("invMvsDaughterPionPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterPionPhiFlat[Lambda]          = new TProfile(Form("invMvsDaughterPionPhiFlat%d",Lambda),Form("invMvsDaughterPionPhiFlat%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterPionMom[Lambda]              = new TProfile(Form("invMvsDaughterPionMom%d",Lambda),Form("invMvsDaughterPionMom%d;invM[GeV^{2}/c^{4}];p",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonDeltaPhi[Lambda]       = new TProfile(Form("invMvsDaughterProtonDeltaPhi%d",Lambda),Form("invMvsDaughterProtonDeltaPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonPhi[Lambda]            = new TProfile(Form("invMvsDaughterProtonPhi%d",Lambda),Form("invMvsDaughterProtonPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonPhiFlat[Lambda]        = new TProfile(Form("invMvsDaughterProtonPhiFlat%d",Lambda),Form("invMvsDaughterProtonPhiFlat%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonMom[Lambda]            = new TProfile(Form("invMvsDaughterProtonMom%d",Lambda),Form("invMvsDaughterProtonMom%d;invM[GeV^{2}/c^{4}];p",Lambda),30,1.10,1.13);
				pInvMvsLambdaDeltaPhi[Lambda]               = new TProfile(Form("invMvsLambdaDeltaPhi%d",Lambda),Form("invMvsLambdaDeltaPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsLambdaPhi[Lambda]                    = new TProfile(Form("invMvsLambdaPhi%d",Lambda),Form("invMvsLambdaPhi%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsLambdaPhiFlat[Lambda]                = new TProfile(Form("invMvsLambdaPhiFlat%d",Lambda),Form("invMvsLambdaPhiFlat%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsLambdaMom[Lambda]                    = new TProfile(Form("invMvsLambdaMom%d",Lambda),Form("invMvsLambdaMom%d;invM[GeV^{2}/c^{4}];p",Lambda),30,1.10,1.13);
				pInvMvsDaughterPionPhiRF[Lambda]            = new TProfile(Form("invMvsDaughterPionPhiRF%d",Lambda),Form("invMvsDaughterPionPhiRF%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterPionMomRF[Lambda]            = new TProfile(Form("invMvsDaughterPionMomRF%d",Lambda),Form("invMvsDaughterPionMomRF%d;invM[GeV^{2}/c^{4}];p",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonPhiRF[Lambda]          = new TProfile(Form("invMvsDaughterProtonPhiRF%d",Lambda),Form("invMvsDaughterProtonPhiRF%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonDeltaPhiRF[Lambda]     = new TProfile(Form("invMvsDaughterProtonDeltaPhiRF%d",Lambda),Form("invMvsDaughterProtonDeltaPhiRF%d;invM[GeV^{2}/c^{4}];<#Psi_{1}-#phi>",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonMomRF[Lambda]          = new TProfile(Form("invMvsDaughterProtonMomRF%d",Lambda),Form("invMvsDaughterProtonMomRF%d;invM[GeV^{2}/c^{4}];p",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonPhiRFFlat[Lambda]      = new TProfile(Form("invMvsDaughterProtonPhiRFFlat%d",Lambda),Form("invMvsDaughterProtonPhiRFFlat%d;invM[GeV^{2}/c^{4}];#phi",Lambda),30,1.10,1.13);
				pInvMvsDaughterProtonDeltaPhiRFFlat[Lambda] = new TProfile(Form("invMvsDaughterProtonDeltaPhiRFFlat%d",Lambda),Form("invMvsDaughterProtonDeltaPhiRFFlat%d;invM[GeV^{2}/c^{4}];<#Psi_{1}-#phi>",Lambda),30,1.10,1.13);
				for(int y=0;y<6;y++){
						pInvMvsDaughterPionPt[Lambda][y]         = new TProfile(Form("invMvsDaughterPionRapidity%dPt%d",y,Lambda),Form("invMvsDaughterPionRapidity%dPt%d;invM[GeV^{2}/c^{4};pT[GeV/c^{2}]",y,Lambda),30,1.10,1.13);
						pInvMvsDaughterProtonPt[Lambda][y]       = new TProfile(Form("invMvsDaughterProtonRapidity%dPt%d",y,Lambda),Form("invMvsDaughterProtonRapidity%dPt%d;invM[GeV^{2}/c^{4};pT[GeV/c^{2}]",y,Lambda),30,1.10,1.13);
						pInvMvsLambdaPt[Lambda][y]               = new TProfile(Form("invMvsLambdaRapidity%dPt%d",y,Lambda),Form("invMvsLambdaRapidity%dPt%d;invM[GeV^{2}/c^{4};pT[GeV/c^{2}]",y,Lambda),30,1.10,1.13);
						pInvMvsDaughterPionPtRF[Lambda][y]         = new TProfile(Form("invMvsDaughterPionRFRapidity%dPt%d",y,Lambda),Form("invMvsDaughterPionRFRapidity%dPt%d;invM[GeV^{2}/c^{4};pT[GeV/c^{2}]",y,Lambda),30,1.10,1.13);
						pInvMvsDaughterProtonPtRF[Lambda][y]       = new TProfile(Form("invMvsDaughterProtonRFRapidity%dPt%d",y,Lambda),Form("invMvsDaughterProtonRFRapidity%dPt%d;invM[GeV^{2}/c^{4};pT[GeV/c^{2}]",y,Lambda),30,1.10,1.13);
						for(int pt=0;pt<5;pt++){
								pInvMvsDaughterPionPtDivide[Lambda][y][pt] = new TProfile(Form("invMvsDaughterPion%dRapidity%dPt%d",Lambda,y,pt),Form("invMvsDaughterPion%dRapidity%dPt%d;invM;p_{T}",Lambda,y,pt),30,1.10,1.13);
								pInvMvsDaughterProtonPtDivide[Lambda][y][pt] = new TProfile(Form("invMvsDaughterProton%dRapidity%dPt%d",Lambda,y,pt),Form("invMvsDaughterProton%dRapidity%dPt%d;invM;p_{T}",Lambda,y,pt),30,1.10,1.13);
								pInvMvsLambdaPtDivide[Lambda][y][pt] = new TProfile(Form("invMvsLambda%dRapidity%dPt%d",Lambda,y,pt),Form("invMvsLambda%dRapidity%dPt%d;invM;p_{T}",Lambda,y,pt),30,1.10,1.13);
						}
				}
				for(int pt=0;pt<12;pt++){
						pInvMvsDaughterPionRapidity[Lambda][pt]   = new TProfile(Form("invMvsDaughterPionPt%dRapidity%d",pt,Lambda),Form("invMvsDaughterPionPt%dRapidity%d;invM[GeV^{2}/c^{4}];y",pt,Lambda),30,1.10,1.13);
						pInvMvsDaughterProtonRapidity[Lambda][pt] = new TProfile(Form("invMvsDaughterProtonPt%dRapidity%d",pt,Lambda),Form("invMvsDaughterProtonPt%dRapidity%d;invM[GeV^{2}/c^{4}];y",pt,Lambda),30,1.10,1.13);
						pInvMvsLambdaRapidity[Lambda][pt]         = new TProfile(Form("invMvsLambdaPt%dRapidity%d",pt,Lambda),Form("invMvsLambdaPt%dRapidity%d;invM[GeV^{2}/c^{4}];y",pt,Lambda),30,1.10,1.13);
						pInvMvsDaughterPionRapidityRF[Lambda][pt]   = new TProfile(Form("invMvsDaughterPionRFPt%dRapidity%d",pt,Lambda),Form("invMvsDaughterPionRFPt%dRapidity%d;invM[GeV^{2}/c^{4}];y",pt,Lambda),30,1.10,1.13);
						pInvMvsDaughterProtonRapidityRF[Lambda][pt] = new TProfile(Form("invMvsDaughterProtonRFPt%dRapidity%d",pt,Lambda),Form("invMvsDaughterProtonRFPt%dRapidity%d;invM[GeV^{2}/c^{4}];y",pt,Lambda),30,1.10,1.13);
						for(int y=0;y<6;y++){
								pInvMvsDaughterPionRapidityDivide[Lambda][pt][y]   = new TProfile(Form("invMvsDaughterPion%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsDaughterPion%dPt%dRapidity%d;invM;y",Lambda,pt,y),30,1.10,1.13);
								pInvMvsDaughterProtonRapidityDivide[Lambda][pt][y] = new TProfile(Form("invMvsDaughterProton%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsDaughterProton%dPt%dRapidity%d;invM;y",Lambda,pt,y),30,1.10,1.13);
								pInvMvsLambdaRapidityDivide[Lambda][pt][y] = new TProfile(Form("invMvsLambda%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsLambda%dPt%dRapidity%d;invM;y",Lambda,pt,y),30,1.10,1.13);
								pInvMvsDaughterPionDeltaPhiDivide[Lambda][pt][y]   = new TProfile(Form("invMvsDaughterPionDeltaPhi%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsDaughterPionDeltaPhi%dPt%dRapidity%d;invM;<#phi-#Psi_{1}>",Lambda,pt,y),30,1.10,1.13);
								pInvMvsDaughterProtonDeltaPhiDivide[Lambda][pt][y] = new TProfile(Form("invMvsDaughterProtonDeltaPhi%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsDaughterProtonDeltaPhi%dPt%dRapidity%d;invM;<#phi-#Psi_{1}>",Lambda,pt,y),30,1.10,1.13);
								pInvMvsLambdaDeltaPhiDivide[Lambda][pt][y] = new TProfile(Form("invMvsLambdaDeltaPhi%dPt%dRapidity%d",Lambda,pt,y),Form("invMvsLambdadeltaPhi%dPt%dRapidity%d;invM;<#phi-#Psi_{1}>",Lambda,pt,y),30,1.10,1.13);
						}
				}
		}

		for(int cent=0;cent<8;cent++){
				for(int rapidity=0;rapidity<6;rapidity++){
						for(int pt=0;pt<12;pt++){
								for(int psi=0;psi<7;psi++){
										//	hinvMLambda[0][cent][rapidity][pt][psi] = new TH1F(Form("invMLambda_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invMLambda_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),340,1.05,1.18);
										//	hinvMLambda[1][cent][rapidity][pt][psi] = new TH1F(Form("invMAntiLambda_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invMAntiLambda_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),340,1.05,1.18);
										//	pInvMassvsSinpol[0][cent][rapidity][pt][psi] = new TProfile(Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dPsi",cent,rapidity,pt,psi),30,1.10,1.13);
										//	pInvMassvsSinpol[1][cent][rapidity][pt][psi] = new TProfile(Form("invMantiLambdavsSinpol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invantiMLambdavsSinpol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),30,1.10,1.13);
										pInvMassvsCospol[0][cent][rapidity][pt][psi] = new TProfile(Form("invMLambdavsCospol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invMLambdavsCospol_Cent%dRapidity%dPt%dPsi",cent,rapidity,pt,psi),50,1.09,1.14);
										pInvMassvsCospol[1][cent][rapidity][pt][psi] = new TProfile(Form("invMantiLambdavsCospol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invantiMLambdavsCospol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),50,1.09,1.14);
										pInvMassvsCos2pol[0][cent][rapidity][pt][psi] = new TProfile(Form("invMLambdavsCos2pol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invMLambdavsCos2pol_Cent%dRapidity%dPt%dPsi",cent,rapidity,pt,psi),50,1.09,1.14);
										pInvMassvsCos2pol[1][cent][rapidity][pt][psi] = new TProfile(Form("invMantiLambdavsCos2pol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),Form("invantiMLambdavsCos2pol_Cent%dRapidity%dPt%dPsi%d",cent,rapidity,pt,psi),50,1.09,1.14);
								}
								for(int dphi=0;dphi<20;dphi++){
										hinvMLambda[0][cent][rapidity][pt][dphi] = new TH1F(Form("invMLambda_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),Form("invMLambda_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),340,1.05,1.18);
										hinvMLambda[1][cent][rapidity][pt][dphi] = new TH1F(Form("invMAntiLambda_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),Form("invMAntiLambda_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),340,1.05,1.18);
										pInvMassvsSinpol[0][cent][rapidity][pt][dphi] = new TProfile(Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dDphi",cent,rapidity,pt,dphi),70,1.08,1.15);
										pInvMassvsSinpol[1][cent][rapidity][pt][dphi] = new TProfile(Form("invMantiLambdavsSinpol_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),Form("invantiMLambdavsSinpol_Cent%dRapidity%dPt%dDphi%d",cent,rapidity,pt,dphi),70,1.08,1.15);
								}
								for(int dPsiPhi=0;dPsiPhi<8;dPsiPhi++){
										hinvMLambdaDpsiphi[0][cent][rapidity][pt][dPsiPhi] = new TH1F(Form("invMLambda_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),Form("invMLambda_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),340,1.05,1.18);
										hinvMLambdaDpsiphi[1][cent][rapidity][pt][dPsiPhi] = new TH1F(Form("invMAntiLambda_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),Form("invMAntiLambda_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),340,1.05,1.18);
										pInvMassvsSinpolDpsiphi[0][cent][rapidity][pt][dPsiPhi] = new TProfile(Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),Form("invMLambdavsSinpol_Cent%dRapidity%dPt%dDpsiphi",cent,rapidity,pt,dPsiPhi),30,1.10,1.13);
										pInvMassvsSinpolDpsiphi[1][cent][rapidity][pt][dPsiPhi] = new TProfile(Form("invMantiLambdavsSinpol_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),Form("invantiMLambdavsSinpol_Cent%dRapidity%dPt%dDpsiphi%d",cent,rapidity,pt,dPsiPhi),30,1.10,1.13);
								}
								for(int Lambda=0;Lambda<2;Lambda++){
										pInvMassvsSinpolPhiCorrect[Lambda][cent][rapidity][pt] = new TProfile(Form("invMLambda%dvsSinpolPhiCorrect_Cent%dRapidity%dPt%d",Lambda,cent,rapidity,pt),Form("invMLambdavsSinpolPhiCorrect_Cent%dRapidity%dPt%d",cent,rapidity,pt),30,1.10,1.13);
								}
								hinvMLambdaMassWidthCorrect[0][cent][rapidity][pt] = new TH1F(Form("invMLambdaMassWidthCorrect_Cent%dRapidity%dPt%d",cent,rapidity,pt),Form("invMLambdaMassWidthCorrect_Cent%dRapidity%dPt%d",cent,rapidity,pt),340,1.05,1.18);
								pInvMassvsSinpolMassWidthCorrect[0][cent][rapidity][pt] = new TProfile(Form("invMLambdavsSinpolMassWidthCorrect_Cent%dRapidity%dPt%d",cent,rapidity,pt),Form("invMLambdavsSinpolMassWidthCorrect_Cent%dRapidity%dPt%d",cent,rapidity,pt),30,1.10,1.13);

								hinvMLambdadEtaCut[cent][rapidity][pt] = new TH1F(Form("invMLambdadEtaCut_Cent%dRapidity%dPt%d",cent,rapidity,pt),Form("invMLambdadEtaCut_Cent%dRapidity%dPt%d",cent,rapidity,pt),340,1.05,1.18);
								pInvMassvsSinpoldEtaCut[cent][rapidity][pt] = new TProfile(Form("invMLambdavsSinpoldEtaCut_Cent%dRapidity%dPt%d",cent,rapidity,pt),Form("invMLambdavsSinpoldEtaCut_Cent%dRapidity%dPt%d",cent,rapidity,pt),30,1.10,1.13);
						}
				}
		}

		//No divided
		for(int dphi=0;dphi<10;dphi++){
				for(int cent=0;cent<8;cent++){
						hinvMCentLambda[cent][dphi] = new TH11F(Form("invMLambdaCent%ddPhi%d",cent,dphi),Form("invMLambdaCent%ddPhi%d",cent,dphi),340,1.05,1.18);
						pINvMvsSinpolCentLambda[cent][dphi] = new TProfile(Form());
				}
		}

		//Lambda bar PH
		for(int cent=0;cent<8;cent++){
				hinvMLambdaBarCent[cent] = new TH1F(Form("invMLambdaBar_Cent%d",cent),Form("invMLambdaBar_Cent%d",cent),340,1.05,1.18);
				pInvMassvsSinpolLambdaBarCent[cent] = new TProfile(Form("InvMassLambdaBarvsSinpol_Cent%d",cent),Form("InvMassLambdaBarvsSinpol_Cent%d",cent),340,1.05,1.18);
		}
		for(int y=0;y<6;y++){
				hinvMLambdaBarRap[y] = new TH1F(Form("invMLambdaBar_Rap%d",y),Form("invMLambdaBar_Rap%d",y),340,1.05,1.18);
				pInvMassvsSinpolLambdaBarRap[y] = new TProfile(Form("InvMassLambdaBarvsSinpol_Rap%d",y),Form("InvMassLambdaBarvsSinpol_Rap%d",y),340,1.05,1.18);
		}
		for(int pt=0;pt<10;pt++){
				hinvMLambdaBarPt[pt] = new TH1F(Form("invMLambdaBar_Pt%d",pt),Form("invMLambdaBar_Pt%d",pt),340,1.05,1.18);
				pInvMassvsSinpolLambdaBarPt[pt] = new TProfile(Form("InvMassLambdaBarvsSinpol_Pt%d",pt),Form("InvMassLambdaBarvsSinpol_Pt%d",pt),340,1.05,1.18);
		}

		hInvMassvsPhistar = new TH2F("InvMassvsPhistar","InvMassvsPhistar",50,1.09,1.14,100,-TMath::Pi(),TMath::Pi());
		hInvMassvsPhistarAfterCorrection = new TH2F("InvMassvsPhistarAfterCorrection","InvMassvsPhistarAfterCorrection",50,1.09,1.14,100,-TMath::Pi(),TMath::Pi());
		//some test for Polarization
		for(int y=0;y<6;y++){
				for(int particle=0;particle<2;particle++){
						for(int cut=0;cut<8;cut++){
								pInvMvsSinpolDaughterCut[y][particle][cut] = new TProfile(Form("invMvsSinpolDaugghterCut_Rapidity%dCutParticle%dCut%d",y,particle,cut),Form("invMvsSinpolDaugghterCut_Rapidity%dCutParticle%dCut%d",y,particle,cut),30,1.10,1.13);
						}
				}
				for(int particle=0;particle<4;particle++){
						pInvMvsSinpolPhiFlat[y][particle] = new TProfile(Form("invMvsSinpolAferPhiFlat_Rapidity%dParticle%d",y,particle),Form("invMvsSinpolAfterPhiFlat_Rapidity%dParticle%d",y,particle),30,1.10,1.13);
				}
		}
		//Polarization w.r.t production plane
		for(int y=0;y<15;y++){
				for(int pt=0;pt<10;pt++){
						hinvMassforProductionPlane[y][pt]  = new TH1F(Form("invMassForProductionPlane_Rapidity%dPt%d",y,pt),Form("invMassForProductionPlane_Rapidity%dPt%d",y,pt),340,1.10,1.13);
						pInvMassvsProductedPlanePol[y][pt] = new TProfile(Form("invMvsProducedPlanePol_Rapidity%dPt%d",y,pt),Form("invMvsProducedPlanePol_Rapidity%dPt%d",y,pt),30,1.10,1.13);
				}
		}

		hInvMvsDeltaPhi          = new TH2F("InvMvsDeltaPhi","InvMvsDeltaPhi",50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		hInvMvsDeltaPhiLambdabar = new TH2F("InvMvsDeltaPhiLambdabar","InvMvsDeltaPhiLambdabar",50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		pInvMvsSinPol            = new TProfile2D("InvMvsSinpol","InvMvsSinpol",50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		pInvMvsSinPolLambdabar   = new TProfile2D("InvMvsSinpolLambdabar","InvMvsSinpolLambdabar",50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		hInvMvsDeltaPhiFull = new TH2F("InvMvsDeltaPhiFull","InvMvsDeltaPhiFull",50,1.09,1.14,100,2*-TMath::Pi(),2*TMath::Pi());
		hInvMvsDeltaPhiTurn = new TH2F("InvMvsDeltaPhiTurn","InvMvsDeltaPhiTurn",50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		pInvMvsSinPolFull   = new TProfile2D("InvMvsSinpolFull","InvMvsSinpolFull",50,1.09,1.14,100,-2*TMath::Pi(),2*TMath::Pi());
		for(int pt=0;pt<5;pt++){
				for(int deta=0;deta<12;deta++){
						hInvMvsDeltaPhi_dEta[deta][pt] = new TH2F(Form("InvMvsDeltaPhi_dEta%dPt%d",deta,pt),Form("InvMvsDeltaPhi_dEta%dPt%d",deta,pt),50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
						pInvMvsSinPol_dEta[deta][pt]   = new TProfile2D(Form("InvMvsSinpol_dEta%dPt%d",deta,pt),Form("InvMvsSinpol_dEta%dPt%d",deta,pt),50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
				}
				for(int mass=0;mass<5;mass++){
						hdEtavsdPhi_invM[pt][mass] = new TH2F(Form("dEtavsdPhi_Pt%dMass%d",pt,mass),Form("dEtavsdPhi_Pt%dMass%d",pt,mass),100,-1.5,1.5,50,-TMath::Pi(),TMath::Pi());
				}
		}
		for(int deta=0;deta<8;deta++){
				for(int pt=0;pt<5;pt++){
						hInvMvsDeltaPhi_DivdEta[deta][pt] = new TH2F(Form("InvMvsDeltaPhi_DivdEta%dPt%d",deta,pt),Form("InvMvsDeltaPhi_DivdEta%dPt%d",deta,pt),50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
				}
		}
		hinvMvsDeltaPhi_dEta3D = new TH3F("InvMvsDeltaPhi_dEta3D","InvMvsDeltaPhi_dEta3D",50,-TMath::Pi(),TMath::Pi(),50,-1.0,1.0,50,1.09,1.14);

		for(int Lambda=0;Lambda<2;Lambda++){
				hInvMvsPhiProton[Lambda]      = new TH2F(Form("InvMvsPhiProton%d",Lambda),Form("InvMvsPhiProton%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsPhiPion[Lambda]        = new TH2F(Form("InvMvsPhiPion%d",Lambda),Form("InvMvsPhiPion%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsPhiLambda[Lambda]      = new TH2F(Form("InvMvsPhiLambda%d",Lambda),Form("InvMvsPhiLambda%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsPhistar[Lambda]        = new TH2F(Form("InvMvsPhistar%d",Lambda),Form("InvMvsPhistar%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsdPsiPhiProton[Lambda]  = new TH2F(Form("InvMvsPhidPsiPhiProton%d",Lambda),Form("InvMvsdPsiPhiProton%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsdPsiPhiPion[Lambda]    = new TH2F(Form("InvMvsPhidPsiPhiPion%d",Lambda),Form("InvMvsdPsiPhiPion%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsdPsiPhiLambda[Lambda]  = new TH2F(Form("InvMvsPhidPsiPhiLambda%d",Lambda),Form("InvMvsdPsiPhiLambda%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsdPsiPhistar[Lambda]    = new TH2F(Form("InvMvsPhidPsiPhistar%d",Lambda),Form("InvMvsdPsiPhistar%d",Lambda),200,1.08,1.18,100,-TMath::Pi(),TMath::Pi());
				hInvMvsDaughterdPhi[Lambda]   = new TH2F(Form("InvMvsDaughterdPhi%d",Lambda),Form("InvMvsDaughterdPhi%d",Lambda),200,1.08,1.18,200,-TMath::Pi(),TMath::Pi());
				hInvMvsdEta[Lambda]           = new TH2F(Form("InvMvsdEta%d",Lambda),Form("InvMvsdEta%d",Lambda),200,1.08,1.18,200,-1,1);
				hDaughterdPhiLaboRest[Lambda] = new TH2F(Form("DaughterdPhiLaboRest%d",Lambda),Form("DaughterdPhiLaboRest%d",Lambda),200,-TMath::Pi(),TMath::Pi(),200,-TMath::Pi(),TMath::Pi());
		}

		for(int cent=0;cent<8;cent++){
				for(int y=0;y<6;y++){
						for(int pt=0;pt<10;pt++){
								for(int dphi=0;dphi<50;dphi++){
										hInvMvsDeltaPhiWidthCorrect[cent][y][pt][dphi] = new TH1F(Form("InvMvsDeltaPhiWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),Form("InvMvsDeltaPhiWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),50,1.09,1.14);
										pInvMvsSinPolWidthCorrect[cent][y][pt][dphi]   = new TProfile(Form("InvMvsSinpolWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),Form("InvMvsSinpolWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),50,1.09,1.14);
										hInvMvsDeltaPhiNotWidthCorrect[cent][y][pt][dphi] = new TH1F(Form("InvMvsDeltaPhiNotWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),Form("InvMvsDeltaPhiNotWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),340,1.09,1.14);
										pInvMvsSinPolNotWidthCorrect[cent][y][pt][dphi]   = new TProfile(Form("InvMvsSinpolNotWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),Form("InvMvsSinpolNotWidthCorrect_Cent%dRapidity%dPt%ddPhi%d",cent,y,pt,dphi),50,1.09,1.14);
								}
						}
				}
		}
		for(int y=0;y<6;y++){
				hInvMvsDeltaPhiWidthCorrectYdivi[y] = new TH2F(Form("InvMvsDeltaPhiWidthCorrectRapidity%d",y),Form("InvMvsDeltaPhiWidthCorrectRapidity%d",y),50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
				pInvMvsSinPolWidthCorrectYdivi[y]   = new TProfile2D(Form("InvMvsSinpolWidthCorrectRapidity%d",y),Form("InvMvsSinpolWidthCorrectRapidity%d",y),50,1.09,1.14,50,-TMath::Pi(),TMath::Pi());
		}
		for(int y=0;y<6;y++){
				pDeltaPhivsSinPol[y] = new TProfile(Form("DeltaPhivsSinpol_Rapidity%d",y),Form("DeltaPhivsSinpol_Rapidity%d",y),30,-TMath::Pi(),TMath::Pi());
				pInvMvsSinPolMassWcorrect[y] = new TProfile(Form("InvMvsSinpolMassWcorrectRapidity%d",y),Form("InvMvsSinpolMassWcorrectRapidity%d",y),30,1.1,1.13);
		}

		//Change Topological cut for systematic
		for(int Lambda=0;Lambda<2;Lambda++){
				for(int topolo=0;topolo<10;topolo++){
						for(int dphi=0;dphi<10;dphi++){
								for(int cent=0;cent<8;cent++){
										hInvMCentTopoloSys[Lambda][topolo][cent][dphi] = new TH1F(Form("InvMLambda%dCent%ddPhi%dTopolo%d",Lambda,cent,dphi,topolo),Form("InvMLambda%dCent%ddPhi%dTopolo%d",Lambda,cent,dphi,topolo),340,1.05,1.18);
										pInvMvsSinpolCentTopoloSys[Lambda][topolo][cent][dphi] = new TProfile(Form("InvMvsSinpolLambda%dCent%ddPhi%dTopolo%d",Lambda,cent,dphi,topolo),Form("InvMvsSinpolLambda%dCent%ddPhi%dTopolo%d",Lambda,cent,dphi,topolo),30,1.1,1.13);
								}
								for(int y=0;y<6;y++){
										hInvMRapTopoloSys[Lambda][topolo][y][dphi] = new TH1F(Form("InvMLambda%dRapidity%ddPhi%dTopolo%d",Lambda,y,dphi,topolo),Form("InvMLambda%dRapidity%ddPhi%dTopolo%d",Lambda,y,dphi,topolo),340,1.05,1.18);
										pInvMvsSinpolRapTopoloSys[Lambda][topolo][y][dphi] = new TProfile(Form("InvMvsSinpolLambda%dRapidity%ddPhi%dTopolo%d",Lambda,y,dphi,topolo),Form("InvMvsSinpolLambda%dRapidity%ddPhi%dTopolo%d",Lambda,y,dphi,topolo),30,1.1,1.13);
								}
								for(int pt=0;pt<10;pt++){
										hInvMPtTopoloSys[Lambda][topolo][pt][dphi] = new TH1F(Form("InvMLambda%dPt%ddPhi%dTopolo%d",Lambda,pt,dphi,topolo),Form("InvMLambda%dPt%ddPhi%dTopolo%d",Lambda,pt,dphi,topolo),340,1.05,1.18);
										pInvMvsSinpolPtTopoloSys[Lambda][topolo][pt][dphi] = new TProfile(Form("InvMvsSinpolLambda%dPt%ddPhi%dTopolo%d",Lambda,pt,dphi,topolo),Form("InvMvsSinpolLambda%dPt%ddPhi%dTopolo%d",Lambda,pt,dphi,topolo),30,1.1,1.13);
								}
						}
				}
		}
		//DecayLvsPt
		for(int mass=0;mass<3;mass++){
				hDecayLvspt[mass] = new TH2F(Form("DecayLvsPt_Mass%d",mass),Form("DecayLvsPt_Mass%d",mass),100,0,10,100,0,5);
		}

		//Event plane method
		for(int dpsiphi=0;dpsiphi<6;dpsiphi++){
				for(int dphi=0;dphi<10;dphi++){
						for(int cent=0;cent<8;cent++){
								hInvMCentEPmethod[cent][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodCent%ddPsiphi%ddPhi%d",cent,dpsiphi,dphi),Form("InvMEPmethodCent%ddPsiphi%ddPhi%d",cent,dpsiphi,dphi),340,1.05,1.18);
						}
						for(int y=0;y<6;y++){
								hInvMRapEPmethod[y][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodRapidity%ddPsiphi%ddPhi%d",y,dpsiphi,dphi),Form("InvMEPmethodRapidity%ddPsiphi%ddPhi%d",y,dpsiphi,dphi),340,1.05,1.18);
						}
						for(int pt=0;pt<10;pt++){
								hInvMPtEPmethod[pt][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodPt%ddPsiphi%ddPhi%d",pt,dpsiphi,dphi),Form("InvMEPmethodPt%ddPsiphi%ddPhi%d",pt,dpsiphi,dphi),340,1.05,1.18);
						}
				}
				for(int dphi=0;dphi<2;dphi++){
						for(int cent=0;cent<8;cent++){
								hInvMLambdaBarCentEPmethod[cent][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodLambdaBarCent%ddPsiphi%ddPhi%d",cent,dpsiphi,dphi),Form("InvMEPmethodLambdaBarCent%ddPsiphi%ddPhi%d",cent,dpsiphi,dphi),340,1.05,1.18);
						}
						for(int y=0;y<6;y++){
								hInvMLambdaBarRapEPmethod[y][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodLambdaBarRapidity%ddPsiphi%ddPhi%d",y,dpsiphi,dphi),Form("InvMEPmethodLambdaBarRapidity%ddPsiphi%ddPhi%d",y,dpsiphi,dphi),340,1.05,1.18);
						}
						for(int pt=0;pt<10;pt++){
								hInvMLambdaBarPtEPmethod[pt][dpsiphi][dphi] = new TH1F(Form("InvMEPmethodLambdaBarPt%ddPsiphi%ddPhi%d",pt,dpsiphi,dphi),Form("InvMEPmethodLambdaBarPt%ddPsiphi%ddPhi%d",pt,dpsiphi,dphi),340,1.05,1.18);
						}
				}
		}

		//----------------  Get  -------------------
		TFile *CalibFile;
		if(QvWeight == 0)CalibFile = TFile::Open(Form("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/QvCalib/Calibration_%dmarge.root",runnumber));
		//if(QvWeight == 1)CalibFile = TFile::Open(Form("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/QvCalib/weight/Calibration_%dmarge.root",runnumber));
		if(QvWeight == 1)CalibFile = TFile::Open(Form("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/QvCalib/weight/removeEPD0/Calibration_%dmarge.root",runnumber));
		//for re-centering
		for(int rapidity=0;rapidity<3;rapidity++){
				getTpc[rapidity][0] = (TProfile*)CalibFile -> Get(Form("TpcQx_rapidity%d",rapidity));
				getTpc[rapidity][1] = (TProfile*)CalibFile -> Get(Form("TpcQy_rapidity%d",rapidity));
		}
		for(int nmipmax=0;nmipmax<4;nmipmax++){
				for(int iring=0;iring<4;iring++){
						for(int iring=0;iring<4;iring++){
								getEpd[nmipmax][iring][0] = (TProfile*)CalibFile -> Get(Form("EpdQx_nmipmax%d_ring%d",nmipmax+2,iring));
								getEpd[nmipmax][iring][1] = (TProfile*)CalibFile -> Get(Form("EpdQy_nmipmax%d_ring%d",nmipmax+2,iring));
						}
				}
		}
		//for flattening
		for(int rapidity=0;rapidity<3;rapidity++){
				mTPC_shift_s[rapidity] = (TProfile2D*)CalibFile->Get(Form("TPC_shiftsin_rapidity%d",rapidity));
				mTPC_shift_c[rapidity] = (TProfile2D*)CalibFile->Get(Form("TPC_shiftcos_rapidity%d",rapidity));
		}
		for(int nmip=0;nmip<4;nmip++){
				for(int iring=0;iring<4;iring++){
						mEPD_shift_s[nmip][iring] = (TProfile2D*)CalibFile->Get(Form("EPD_shiftsin_nmipmax%d_iring%d",nmip+2,iring));
						mEPD_shift_c[nmip][iring] = (TProfile2D*)CalibFile->Get(Form("EPD_shiftcos_nmipmax%d_iring%d",nmip+2,iring));
				}
		}
		TFile *V1weightFile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/QvCalib/weight/removeEPD0/V1EpdTpc.root");
		for(int cent=0;cent<8;cent++){
				for(int pt=0;pt<7;pt++){
						for(int particle=0;particle<2;particle++){
								wV1TpcPion[particle][cent][pt]   = (TProfile*)V1weightFile->Get(Form("CalcV1TPC_Psi3Cent%dPt%dPion%d",cent,pt,particle));
								wV1TpcProton[particle][cent][pt] = (TProfile*)V1weightFile->Get(Form("CalcV1TPC_Psi3Cent%dPt%dProton%d",cent,pt,particle));
						}
				}
				wV1Epd[cent] = (TProfile*)V1weightFile->Get(Form("CalcV1EPD_Psi0Cent%d",cent));
		}
		TFile *V1correctFile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/Correction/LambdaAnalisys_V1correction.root");
		for(int mass=0;mass<7;mass++){
				RapidityvsDeltaPhiBeforeCorrect[mass] = (TH2F*)V1correctFile->Get(Form("RapidityvsDeltaphi_Lambda0Mass%d",mass));
		}
		TFile *Phiflatfile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/Correction/PhiFlat.root");
		pinvMvsPhiforcorrection[0] = (TProfile*)Phiflatfile->Get("invMvsLambdaPhi0");
		pinvMvsPhiforcorrection[1] = (TProfile*)Phiflatfile->Get("invMvsDaughterPionPhi0");
		pinvMvsPhiforcorrection[2] = (TProfile*)Phiflatfile->Get("invMvsDaughterProtonPhi0");

		TFile *PhiRFflatfile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/Correction/PhiFlatRF.root");
		pinvMvsPhiRFforcorrection = (TProfile*)PhiRFflatfile->Get("invMvsDaughterPionPhiRF0");

		//TFile *PhiStarFile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/Correction/PhiStarCorrection.root");
		TFile *PhiStarFile = TFile::Open("/star/u/kokubo/ana/run18/Fxt7.2GeV/AnaPol/Correction/PhiStarCorrection2.root");//loosed bin cut
		pinvMvsPhiStarCorrection = (TH2F*)PhiStarFile->Get("InvMassvsPhistar");



		cout << "End of Histograms" << endl;
		return kStOK;
}
//__________________________________________________________________________________
void PicoAnalyzer::Clear(Option_t *opt)
{
		StMaker::Clear();
}

//__________________________________________________________________________________
Int_t PicoAnalyzer::Finish() {
		cout << "PicoAnalyzer::Finish()\n";
		cout << "\tProcessed " << mEventCounter << " events." << endl;
		//===============================
		//  Write Histograms
		//===============================


		LambdaAnalisysFile = new TFile(mLambdaAnalisysFileName.c_str(),"RECREATE","output");
		hVz                         -> Write();
		hVzAfter                    -> Write();
		hTPCmultiplicity                              -> Write(); 
		hTPCmultiplicityAfterPileup -> Write();
		hMultiplicityvsTofmult                        -> Write();
		hMultiplicityvsTofmultAfterPileup             -> Write();
		hMultiplicityvsTofMatch                       -> Write();
		hMultiplicityvsTofMatchAfterPileup            -> Write();
		hNumOfPrimaryvsNumOfGlobal                    -> Write();
		hNumOfPrimaryvsNumOfGlobalAfterPileup         -> Write();
		hNumOfPionProtonvsMultiplicity                -> Write();
		hNumOfPionProtonvsMultiplicityAfterPileup     -> Write();
		hCountCentrality                              -> Write();
		//pid information
		//hDedx                     -> Write();
		//hDedx_proton              -> Write();
		//hDedx_pion                -> Write();
		//hMM                       -> Write();
		//hMM_proton                -> Write();
		//hMM_pion                  -> Write();
		//for(int Lambda=0;Lambda<1;Lambda++){
		//hDaughterPionRapidityvsPt[Lambda]         -> Write();
		//hDaughterPionDeltaPhi[Lambda]             -> Write();
		//hDaughterPionPhi[Lambda]                  -> Write();
		//hDaughterProtonRapidityvsPt[Lambda]       -> Write();
		//hDaughterProtonDeltaPhi[Lambda]           -> Write();
		//hDaughterProtonPhi[Lambda]                -> Write();
		//hDaughterPionRapidityvsPtRF[Lambda]       -> Write();
		//hDaughterPionPhiRF[Lambda]                -> Write();
		//hDaughterProtonRapidityvsPtRF[Lambda]     -> Write();
		//hDaughterProtonPhiRF[Lambda]              -> Write();
		//hLambdaEta[Lambda]                        -> Write();
		//hLambdaPt[Lambda]                         -> Write();
		//hLambdaDeltaPhi[Lambda]                   -> Write();
		//hLambdaPhi[Lambda]                        -> Write();
		//hLambdaRapidity[Lambda]                   -> Write();
		//hLambdaRapidityvsPt[Lambda]               -> Write();
		//pInvMvsDaughterPionDeltaPhi[Lambda]       -> Write();
		//pInvMvsDaughterPionPhi[Lambda]            -> Write();
		//pInvMvsDaughterPionPhiFlat[Lambda]        -> Write();
		//pInvMvsDaughterPionMom[Lambda]            -> Write();
		//pInvMvsDaughterProtonDeltaPhi[Lambda]     -> Write();
		//pInvMvsDaughterProtonPhi[Lambda]          -> Write();
		//pInvMvsDaughterProtonPhiFlat[Lambda]      -> Write();
		//pInvMvsDaughterProtonMom[Lambda]          -> Write();
		//pInvMvsLambdaDeltaPhi[Lambda]             -> Write();
		//pInvMvsLambdaPhi[Lambda]                  -> Write();
		//pInvMvsLambdaPhiFlat[Lambda]              -> Write();
		//pInvMvsLambdaMom[Lambda]                  -> Write();
		//pInvMvsDaughterPionPhiRF[Lambda]          -> Write();
		//pInvMvsDaughterPionMomRF[Lambda]          -> Write();
		//pInvMvsDaughterProtonPhiRF[Lambda]          -> Write();
		//pInvMvsDaughterProtonDeltaPhiRF[Lambda]     -> Write();
		//pInvMvsDaughterProtonMomRF[Lambda]        -> Write();
		//pInvMvsDaughterProtonPhiRFFlat[Lambda]      -> Write();
		//pInvMvsDaughterProtonDeltaPhiRFFlat[Lambda] -> Write();
		//for(int y=0;y<6;y++){
		//		pInvMvsDaughterPionPt[Lambda][y]         -> Write();
		//		pInvMvsDaughterProtonPt[Lambda][y]       -> Write();
		//		pInvMvsLambdaPt[Lambda][y]               -> Write();
		//		pInvMvsDaughterPionPtRF[Lambda][y]         -> Write();
		//		pInvMvsDaughterProtonPtRF[Lambda][y]       -> Write();
		//		for(int pt=0;pt<5;pt++){
		//				pInvMvsDaughterPionPtDivide[Lambda][y][pt]   -> Write();
		//				pInvMvsDaughterProtonPtDivide[Lambda][y][pt] -> Write();
		//				pInvMvsLambdaPtDivide[Lambda][y][pt]         -> Write();
		//		}
		//}
		//for(int pt=0;pt<10;pt++){
		//		pInvMvsDaughterPionRapidity[Lambda][pt]   -> Write();
		//		pInvMvsDaughterProtonRapidity[Lambda][pt] -> Write();
		//		pInvMvsLambdaRapidity[Lambda][pt]         -> Write();
		//		pInvMvsDaughterPionRapidityRF[Lambda][pt]   -> Write();
		//		pInvMvsDaughterProtonRapidityRF[Lambda][pt] -> Write();
		//		for(int y=0;y<6;y++){
		//				pInvMvsDaughterPionRapidityDivide[Lambda][pt][y]   -> Write();
		//				pInvMvsDaughterProtonRapidityDivide[Lambda][pt][y] -> Write();
		//				pInvMvsLambdaRapidityDivide[Lambda][pt][y]         -> Write();
		//				pInvMvsDaughterPionDeltaPhiDivide[Lambda][pt][y]   -> Write();
		//				pInvMvsDaughterProtonDeltaPhiDivide[Lambda][pt][y] -> Write();
		//				pInvMvsLambdaDeltaPhiDivide[Lambda][pt][y]         -> Write();
		//		}
		//}
		//}
		//for(int Lambda=0;Lambda<2;Lambda++){
		//		for(int topolo=0;topolo<5;topolo++){
		//				hinvMLambda_topolo[Lambda][topolo] -> Write();
		//		}
		//		for(int cent=0;cent<8;cent++){
		//				for(int rap=0;rap<6;rap++){
		//						for(int pt=0;pt<12;pt++){
		//								//for(int psi=0;psi<7;psi++){
		//								//hinvMLambda[Lambda][cent][rap][pt][3] -> Write();
		//								//pInvMassvsSinpol[Lambda][cent][rap][pt][3] -> Write();
		//								pInvMassvsCospol[Lambda][cent][rap][pt][3] -> Write();
		//								pInvMassvsCos2pol[Lambda][cent][rap][pt][3] -> Write();
		//								//}
		//						}
		//				}
		//		}
		//}
		//for(int mass=0;mass<7;mass++){
		//		RapidityvsDeltaphiAfterCorrect[mass] -> Write();
		//}
		//for(int Lambda=0;Lambda<1;Lambda++){
		//		for(int cent=0;cent<8;cent++){
		//				for(int rapidity=0;rapidity<6;rapidity++){
		//						for(int pt=0;pt<10;pt++){
		//								pInvMassvsSinpolPhiCorrect[Lambda][cent][rapidity][pt] -> Write();
		//						}
		//				}
		//		}
		//}
		//for(int y=0;y<15;y++){
		//		for(int pt=0;pt<10;pt++){
		//				hinvMassforProductionPlane[y][pt]     -> Write();
		//				pInvMassvsProductedPlanePol[y][pt] -> Write();
		//		}
		//}
		//for(int y=0;y<6;y++){
		//		for(int particle=0;particle<2;particle++){//0 -> Daughter pion 1 -> Daughter proton
		//				for(int cut=0;cut<8;cut++){
		//						pInvMvsSinpolDaughterCut[y][particle][cut] -> Write();
		//				}
		//		}
		//		//for(int particle=0;particle<4;particle++){//0 -> all particle 1->Lambda 2-> Daughter pion 3-> Daughter proton
		//		//		pInvMvsSinpolPhiFlat[y][particle] -> Write();
		//		//}
		//}

		//----------   koko  -----------------
		//for(int Lambda=1;Lambda<2;Lambda++){
		//		for(int cent=0;cent<8;cent++){
		//				for(int rap=0;rap<6;rap++){
		//						for(int pt=0;pt<10;pt++){
		//								//hinvMLambdaMassWidthCorrect[Lambda][cent][rap][pt] -> Write();
		//								//pInvMassvsSinpolMassWidthCorrect[Lambda][cent][rap][pt] -> Write();
		//								for(int dphi=0;dphi<20;dphi++){
		//										hinvMLambda[Lambda][cent][rap][pt][dphi] -> Write();
		//										pInvMassvsSinpol[Lambda][cent][rap][pt][dphi] -> Write();
		//								}
		//								//for(int dPsiPhi=0;dPsiPhi<8;dPsiPhi++){
		//								//		hinvMLambdaDpsiphi[Lambda][cent][rap][pt][dPsiPhi] -> Write();
		//								//		pInvMassvsSinpolDpsiphi[Lambda][cent][rap][pt][dPsiPhi] -> Write();
		//								//}
		//						}
		//				}
		//		}
		//}

		//No divided
		for(int dphi=0;dphi<10;dphi++){
				for(int cent=0;cent<8;cent++){
						hinvMCentLambda[cent][dphi] -> Write();
						pInvMvsSinpolCentLambda[cent][dphi] -> Write();
				}
				for(int y=0：wq;y<6;y++){
						hinvMRapLambda[y][dphi] -> Write();
						pInvMvsSinpolRapLambda[y][dphi] -> Write();
				}
				for(int pt=0;pt<10;pt++){
						hinvMPtLambda[pt][dphi] -> Write();
						pInvMvsSinpolPtLambda[pt][dphi] -> Write();
				}
		}

		//Lambda bar PH
		//for(int cent=0;cent<8;cent++){
		//		hinvMLambdaBarCent[cent]->Write();
		//		pInvMassvsSinpolLambdaBarCent[cent]->Write();
		//}
		//for(int y=0;y<6;y++){
		//		hinvMLambdaBarRap[y]->Write();
		//		pInvMassvsSinpolLambdaBarRap[y]->Write();
		//}
		//for(int pt=0;pt<10;pt++){
		//		hinvMLambdaBarPt[pt]->Write();
		//		pInvMassvsSinpolLambdaBarPt[pt]->Write();
		//}

		//for(int cent=0;cent<8;cent++){
		//		for(int rap=0;rap<6;rap++){
		//				for(int pt=0;pt<12;pt++){
		//						hinvMLambdadEtaCut[cent][rap][pt] -> Write();
		//						pInvMassvsSinpoldEtaCut[cent][rap][pt] -> Write();
		//				}
		//		}
		//}

		//hInvMassvsPhistar -> Write();
		//hInvMassvsPhistarAfterCorrection -> Write();
		//hTPCphivsTOFphi        -> Write();
		//hTPCphivsTOFphiDproton -> Write();
		//hDeltaphiTpcTof        -> Write();
		//hTofXvsY               -> Write();
		//hTofXvsZ               -> Write();
		//hTofYvsZ               -> Write();
		//hDeltaPhi2D            -> Write();
		//for(int dphi=0;dphi<2;dphi++){//0:phiLambda-phi*p<0  1:phiLambda0phi*p>0
		//		hDeltaPhiPionvsProton[dphi] -> Write();
		//
		//for(int charge=0;charge<2;charge++){
		//		hDeltaPhivsPt[charge]->Write();
		//		hDeltaEtavsPt[charge]->Write();
		//}
		//hDeltaLambdaProtonPhi_LaboV0RF->Write();
		//hDeltaLambdaProtonPhivsPionPhi->Write();
		//hInvMvsDeltaPhi        -> Write();
		//hInvMvsDeltaPhiLambdabar -> Write();
		//pInvMvsSinPol     -> Write();
		//pInvMvsSinPolLambdabar     -> Write();
		//hInvMvsDeltaPhiFull -> Write();
		//pInvMvsSinPolFull     -> Write();
		//hInvMvsDeltaPhiTurn -> Write();
		//for(int cent=0;cent<8;cent++){
		//		for(int y=0;y<6;y++){
		//				for(int pt=0;pt<10;pt++){
		//						for(int dphi=0;dphi<50;dphi++){
		//								//hInvMvsDeltaPhiWidthCorrect[cent][y][pt][dphi] -> Write();
		//								//pInvMvsSinPolWidthCorrect[cent][y][pt][dphi]   -> Write();
		//								hInvMvsDeltaPhiNotWidthCorrect[cent][y][pt][dphi] -> Write();
		//								pInvMvsSinPolNotWidthCorrect[cent][y][pt][dphi]   -> Write();
		//						}
		//				}
		//		}
		//}
		//for(int y=0;y<6;y++){
		//		hInvMvsDeltaPhiWidthCorrectYdivi[y] -> Write();
		//		pInvMvsSinPolWidthCorrectYdivi[y]   -> Write();
		//		pDeltaPhivsSinPol[y]                -> Write();
		//		pInvMvsSinPolMassWcorrect[y]        -> Write();
		//}
		//for(int deta=0;deta<12;deta++){
		//		for(int pt=0;pt<5;pt++){
		//				hInvMvsDeltaPhi_dEta[deta][pt] -> Write();
		//				pInvMvsSinPol_dEta[deta][pt]   -> Write();
		//		}
		//}
		//for(int mass=0;mass<5;mass++){
		//		for(int pt=0;pt<5;pt++){
		//				hdEtavsdPhi_invM[pt][mass] -> Write();
		//		}
		//}
		//for(int deta=0;deta<8;deta++){
		//		for(int pt=0;pt<5;pt++){
		//				hInvMvsDeltaPhi_DivdEta[deta][pt] -> Write();
		//		}
		//}
		//for(int Lambda=0;Lambda<2;Lambda++){
		//		hInvMvsPhiProton[Lambda]      -> Write();
		//		hInvMvsPhiPion[Lambda]        -> Write();
		//		hInvMvsPhiLambda[Lambda]      -> Write();
		//		hInvMvsPhistar[Lambda]        -> Write();
		//		hInvMvsdPsiPhiProton[Lambda]  -> Write();
		//		hInvMvsdPsiPhiPion[Lambda]    -> Write();
		//		hInvMvsdPsiPhiLambda[Lambda]  -> Write();
		//		hInvMvsdPsiPhistar[Lambda]    -> Write();
		//		hInvMvsDaughterdPhi[Lambda]   -> Write();
		//		hInvMvsdEta[Lambda]           -> Write();
		//		hDaughterdPhiLaboRest[Lambda] -> Write();
		//}

		//Change Topological cut for systematic
		//for(int Lambda=0;Lambda<1;Lambda++){
		//		for(int topolo=0;topolo<10;topolo++){
		//				for(int dphi=0;dphi<10;dphi++){
		//						for(int cent=0;cent<8;cent++){
		//								hInvMCentTopoloSys[Lambda][topolo][cent][dphi]->Write();
		//								pInvMvsSinpolCentTopoloSys[Lambda][topolo][cent][dphi]->Write();
		//						}
		//						for(int y=0;y<6;y++){
		//								hInvMRapTopoloSys[Lambda][topolo][y][dphi]->Write();
		//								pInvMvsSinpolRapTopoloSys[Lambda][topolo][y][dphi]->Write();
		//						}
		//						for(int pt=0;pt<10;pt++){
		//								hInvMPtTopoloSys[Lambda][topolo][pt][dphi]->Write();
		//								pInvMvsSinpolPtTopoloSys[Lambda][topolo][pt][dphi]->Write();
		//						}
		//				}
		//		}
		//}

		//DecayLvsPt
		//for(int mass=0;mass<3;mass++){
		//		hDecayLvspt[mass] -> Write();
		//}

		//Event plane method
		//Lambda
		for(int dphi=0;dphi<10;dphi++){
				for(int dPsiphi=0;dPsiphi<6;dPsiphi++){
						for(int cent=0;cent<8;cent++){
								hInvMCentEPmethod[cent][dPsiphi][dphi]->Write();
						}
						for(int y=0;y<6;y++){
								hInvMRapEPmethod[y][dPsiphi][dphi]->Write();
						}
						for(int pt=0;pt<10;pt++){
								hInvMPtEPmethod[pt][dPsiphi][dphi]->Write();
						}
				}
		}
		//Anti-Lambda
		for(int dphi=0;dphi<2;dphi++){
				for(int dPsiphi=0;dPsiphi<6;dPsiphi++){
						for(int cent=0;cent<8;cent++){
								hInvMLambdaBarCentEPmethod[cent][dPsiphi][dphi]->Write();
						}
						for(int y=0;y<6;y++){
								hInvMLambdaBarRapEPmethod[y][dPsiphi][dphi]->Write();
						}
						for(int pt=0;pt<10;pt++){
								hInvMLambdaBarPtEPmethod[pt][dPsiphi][dphi]->Write();
						}
				}
		}

		LambdaAnalisysFile -> Write();
		LambdaAnalisysFile -> Close();




		return kStOK;
}
//__________________________________________________________________________________
Int_t PicoAnalyzer::GetRunIndex( const Int_t run ) {
		Int_t runindex = 9999;
		for(Int_t i=0; i<nRun; i++){
				if(run==RunNumber[i]) runindex = i; 
		}
		return runindex;
}
//__________________________________________________________________________________
Int_t PicoAnalyzer::Make() {
		mEventCounter++;
		//Begining of Event loop  

		//------------------------------------------------------------------
		if(!mPicoDstMaker) {
				LOG_WARN << " No PicoDstMaker! Skip! " << endm;
				return kStWarn;
		}
		mPicoDst = mPicoDstMaker->picoDst();
		if(!mPicoDst) {
				LOG_WARN << " No PicoDst! Skip! " << endm;
				return kStWarn;
		}
		picoEvent = (StPicoEvent*)mPicoDst->event();
		if( !picoEvent ){
				LOG_WARN << " No PicoEvent! Skip! " << endm;
				return kStWarn;
		}
		//------------------------------------------------------------------



		// check of triggerID ----------------------------------------------
		bool triggerId = 0;
		vector<unsigned int> trigIds = picoEvent->triggerIds();
		//  for(Int_t i=0; i<nTrigID; i++){
		//    if(picoEvent->isTrigger(TrigID[i])){
		//      triggerId = true;
		//    }
		//  }
		for(int i=0;i<(int)trigIds.size();i++) {
				//cout << trigIds[i] << endl;
				for(int j=0;j<nTrigID;j++) {
						if( TrigID[j]==trigIds[i] )  triggerId = true;
				}
				//cout << i << " " << trigIds[i] << " " << triggerId << endl;  
		}
		if(triggerId==false)  return kStOK;
		// get event info. -------------------------------------------------
		const int runnumber = picoEvent->runId();
		const int eventid = picoEvent->eventId();
		const int refmult = picoEvent->refMult();
		const int grefmult = picoEvent->grefMult();
		// Run index for QA
		const int runindex = GetRunIndex(runnumber); 
		// Vertex
		TVector3 pVertex = picoEvent->primaryVertex();
		StThreeVectorF PrimaryVertex(pVertex.x(),pVertex.y(),pVertex.z());
		const float Vzvpd = picoEvent->vzVpd();
		const float BField = picoEvent->bField();
		const Int_t nTrack = mPicoDst->numberOfTracks();
		const UShort_t tofmult = picoEvent->btofTrayMultiplicity();
		const float nBTofMatch = picoEvent->nBTOFMatch();
		const float pver_x = pVertex.x();
		const float pver_y = pVertex.y();
		const float pver_z = pVertex.z();
		const float pver_r = sqrt(pver_x*pver_x + (pver_y+2)*(pver_y+2));
		float ZDCx = picoEvent->ZDCx();
		float BBCx = picoEvent->BBCx();
		int tofmatch = 0;
		float VzRank = picoEvent->ranking();

		hVz -> Fill(pver_z);

		double QvTpcRaw[3][2]   = {0};
		double QvTpcRece[3][2]  = {0};
		double QvTpcmean[3][2]  = {0};
		double QvTpcsigma[3][2] = {0};
		double PsiTpcRaw[3]     = {0};
		double PsiTpcRece[3]    = {0};
		double PsiTpcFlat[3]    = {0};
		double WeightSumTpc[3]  = {0};
		int CountTpc[3] = {0};
		int NumGlobal = 0;
		int NumPrimary = 0;
		int NumPionProton = 0;

		for(int i=0;i<3;i++){
				WeightSumTpc[i]  = {0};
				CountTpc[i] = {0};
		}



		int Multiplicity =0;

		// primary track loop for determine refmult ----------------------------------------------
		for (Int_t itr=0;itr<nTrack;itr++) { 
				const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);

				//if(!isGoodTrack(ptrk))  continue;
				NumGlobal++;
				if(!ptrk)  continue;
				if(!ptrk->isPrimary())  continue;  // now selecting primary tracks
				NumPrimary++;

				const Float_t dca = ptrk->gDCA( pVertex ).Mag();
				const Int_t nHitsFit = ptrk->nHitsFit();
				const Int_t nHitsPoss = ptrk->nHitsMax();
				const Float_t nHitsR = (Float_t)nHitsFit/(Float_t)nHitsPoss;

				if(dca > 3) continue;
				if(nHitsFit < 10) continue;
				if(nHitsR < 0.52) continue;

				Multiplicity++;
		}

		int centnumber = mCentMaker -> CentDifine(Multiplicity);

		// primary track loop for calculate TPC Qvector ----------------------------------------------
		for (Int_t itr=0;itr<nTrack;itr++) { 
				const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);

				if(!ptrk)  continue;
				if(!ptrk->isPrimary())  continue;  // now selecting primary tracks

				const Float_t p_p = ptrk->pMom().Mag();
				const Float_t p_pt  = ptrk->pMom().Perp(); // zero for global tracks
				const Float_t p_phi = ptrk->pMom().Phi();  // zero for global track
				const Float_t p_eta = ptrk->pMom().PseudoRapidity();
				const Float_t dca = ptrk->gDCA( pVertex ).Mag();
				const Int_t nHitsFit = ptrk->nHitsFit();
				const Int_t nHitsPoss = ptrk->nHitsMax();
				const Float_t nHitsR = (Float_t)nHitsFit/(Float_t)nHitsPoss;
				const Float_t Dedx = ptrk->dEdx();  // [keV/cm] 
				const Short_t Charge = ptrk->charge();
				const Float_t nSigmaPion = ptrk->nSigmaPion();
				const Float_t nSigmaProton = ptrk->nSigmaProton();
				const Int_t bTofPidTraitsIndex = ptrk->bTofPidTraitsIndex();
				const StThreeVectorF origin(ptrk->origin().X(),ptrk->origin().Y(),ptrk->origin().Z());

				float beta = -1.;
				float mm = 0;

				if(dca > 3) continue;
				if(nHitsFit < 10) continue;
				if(nHitsR < 0.52) continue;

				//calculate Start time to remove pile up event
				if( bTofPidTraitsIndex>=0) {
						StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
						TVector3 tofPos = btofPidTraits->btofHitPos();
						const StThreeVectorF btofHitPos(tofPos.X(),tofPos.Y(),tofPos.Z());
						float L = tofPathLength( &origin, &btofHitPos, ptrk->helix(picoEvent->bField()).curvature())*0.01;  // cm -> m
						float tof = btofPidTraits->btof();
						NumPionProton = mCentMaker -> CountForT0cut(p_p, nSigmaPion, nSigmaProton, bTofPidTraitsIndex, tof, L, NumPionProton);
				}


				int pPionID = mMyPidMaker   -> PrimaryPion( bTofPidTraitsIndex, nSigmaPion, mm, Charge);
				int pProtonID = mMyPidMaker -> PrimaryProton( p_pt, bTofPidTraitsIndex, nSigmaProton, mm, Charge);
				int ptID  = mQvMaker -> PtId(p_pt);
				int etaID = mQvMaker -> TpcEtaGroup(p_eta);
				int etaIDforWeight = mQvMaker -> TpcEtaForWeight(p_eta);
				double WeightQvTpc;
				if(ptID != -1 && QvWeight == 1){
						if(pProtonID==0 || pProtonID==1) WeightQvTpc = wV1TpcProton[pProtonID][centnumber][ptID] -> GetBinContent(etaIDforWeight+1);
						if(pPionID==0   || pPionID==1  ) WeightQvTpc = wV1TpcPion[pPionID][centnumber][ptID] -> GetBinContent(etaIDforWeight+1);//if PionID != -1 && Proton != -1 -> it is regarded as Pion
				}
				bool QvTrackFlag = mQvMaker -> ContFlagQvTpc(QvWeight, etaID, ptID, pPionID, pProtonID);
				if(!QvTrackFlag) continue;
				double weight;
				if(QvWeight==0) weight = p_eta - midrapidity;
				if(QvWeight==1) weight = WeightQvTpc;
				QvTpcRaw[etaID][0] += weight*cos(p_phi);//Qx
				QvTpcRaw[etaID][1] += weight*sin(p_phi);//Qy
				WeightSumTpc[etaID] += weight;
				CountTpc[etaID]++;
		}
		// End of track loop ---------------------------------------------



		bool PileUpFlag = mCentMaker -> PileUpCut(Multiplicity, nBTofMatch, tofmult, NumGlobal, NumPrimary, 2, NumPionProton);
		//Event Cut 
		if( pver_z < 199. || 202. < pver_z  || pver_r > 2.) return kStOK;
		if(nBTofMatch <= 0) return kStOK;
		hTPCmultiplicity -> Fill(Multiplicity);
		hMultiplicityvsTofmult     -> Fill(Multiplicity,tofmult);
		hMultiplicityvsTofMatch    -> Fill(Multiplicity, nBTofMatch);
		hNumOfPrimaryvsNumOfGlobal -> Fill(NumPrimary, NumGlobal);
		hNumOfPionProtonvsMultiplicity   -> Fill(Multiplicity,NumPionProton);
		if(!PileUpFlag) return kStOK;
		hTPCmultiplicityAfterPileup -> Fill(Multiplicity);
		hMultiplicityvsTofmultAfterPileup     -> Fill(Multiplicity,tofmult);
		hMultiplicityvsTofMatchAfterPileup    -> Fill(Multiplicity, nBTofMatch);
		hNumOfPrimaryvsNumOfGlobalAfterPileup -> Fill(NumPrimary, NumGlobal);
		hNumOfPionProtonvsMultiplicityAfterPileup   -> Fill(Multiplicity,NumPionProton);
		hVzAfter -> Fill(pver_z);
		hCountCentrality                      -> Fill(centnumber);



		if(centnumber == -1) return kStOK;



		for(int etaId=0;etaId<3;etaId++){
				for(int ixy=0;ixy<2;ixy++){
						//if(CountTpc[etaId] > 1) QvTpcRaw[etaId][ixy] /=sqrt( fabs(WeightSumTpc[etaId]));
						//else QvTpcRaw[etaId][ixy] = -9999;
						if(CountTpc[etaId] <= 1)QvTpcRaw[etaId][ixy] = -9999;
				}
				PsiTpcRaw[etaId] = atan2(QvTpcRaw[etaId][1],QvTpcRaw[etaId][0]);
				if(QvTpcRaw[etaId][0] != -9999 && QvTpcRaw[etaId][1] != -9999){
						hQvTpcRaw[centnumber][etaId] -> Fill(QvTpcRaw[etaId][0],QvTpcRaw[etaId][1]);
						hPsiTpcRaw[centnumber][etaId] -> Fill(PsiTpcRaw[etaId]);
				}
				//for re-centring
				for(int ixy=0;ixy<2;ixy++){
						if(QvTpcRaw[etaId][ixy] != -9999)pTpcQv[etaId][ixy] -> Fill(centnumber,QvTpcRaw[etaId][ixy]);
				}
		}

		//re-centraing
		for(int rapidity=0;rapidity<3;rapidity++){
				for(int ixy=0;ixy<2;ixy++){
						QvTpcmean[rapidity][ixy]  = getTpc[rapidity][ixy] -> GetBinContent(centnumber+1);
						QvTpcsigma[rapidity][ixy] = getTpc[rapidity][ixy] -> GetBinError(centnumber+1);
						QvTpcRece[rapidity][ixy]  = (QvTpcRaw[rapidity][ixy] - QvTpcmean[rapidity][ixy])/QvTpcsigma[rapidity][ixy];
				}
				PsiTpcRece[rapidity] = atan2(QvTpcRece[rapidity][1],QvTpcRece[rapidity][0]);
				if(QvTpcRaw[rapidity][0] != -9999 && QvTpcRaw[rapidity][1] != -9999)hQvTpcRece[centnumber][rapidity] -> Fill(QvTpcRece[rapidity][0],QvTpcRece[rapidity][1]);
				if(QvTpcRaw[rapidity][0] != -9999 && QvTpcRaw[rapidity][1] != -9999)hPsiTpcRece[centnumber][rapidity] -> Fill(PsiTpcRece[rapidity]);
		}
		//flattening
		for(int nfill=1;nfill<9;nfill++){
				for(int rapidity=0;rapidity<3;rapidity++){
						if(QvTpcRaw[rapidity][0] != -9999 && QvTpcRaw[rapidity][1] != -9999)pTPC_shift_sin[rapidity]->Fill(nfill-1,centnumber,sin(nfill*PsiTpcRece[rapidity]));
						if(QvTpcRaw[rapidity][0] != -9999 && QvTpcRaw[rapidity][1] != -9999)pTPC_shift_cos[rapidity]->Fill(nfill-1,centnumber,cos(nfill*PsiTpcRece[rapidity]));
				}
		}
		float mTPC_shift_sin[8][3];
		float mTPC_shift_cos[8][3];
		for(int n=0;n<8;n++){
				for(int rapidity=0;rapidity<3;rapidity++){
						mTPC_shift_sin[n][rapidity] = mTPC_shift_s[rapidity]->GetBinContent(n+1, centnumber+1);
						mTPC_shift_cos[n][rapidity] = mTPC_shift_c[rapidity]->GetBinContent(n+1, centnumber+1);
				}
		}
		for(int rapidity=0;rapidity<3;rapidity++){
				PsiTpcFlat[rapidity] = PsiTpcRece[rapidity];
				for(int i=1;i<9;i++){
						PsiTpcFlat[rapidity] += 2*(-mTPC_shift_sin[i-1][rapidity]*cos(i*PsiTpcRece[rapidity])+mTPC_shift_cos[i-1][rapidity]*sin(i*PsiTpcRece[rapidity]))/(Float_t)(i);
				}
				if(QvTpcRaw[rapidity][0] != -9999 && QvTpcRaw[rapidity][1] != -9999)hPsiTpcFlat[centnumber][rapidity] -> Fill(PsiTpcFlat[rapidity]);
		}





		//-------------------Get EPD information------------------------
		Int_t nepdHits = mPicoDst->numberOfEpdHits();



		StPicoEpdHit *epdHit;
		TVector3 StraightLine_random[nepdHits];
		float phi_epd[nepdHits];
		float eta_epd[nepdHits];

		float mip       =0;
		float nepdMIPsE =0;
		float nepdMIPsW =0;

		float EpdnMipSumE = 0;
		float EpdnMipSumW = 0;

		double QvEpdRaw[4][4][2]       = {0};//[nMipmax=2,3,4,5][ringgroup][x or y]
		double QvEpdRece[4][4][2]      = {0};//[nMipmax=2,3,4,5][ringgroup][x or y] after re-centering
		double QvEpdmean[4][4][2]      = {0};
		double QvEpdsigma[4][4][2]     = {0};
		float PsiEpdRaw [4][4]        = {0};
		float PsiEpdRece[4][4]        = {0};
		float PsiEpdFlat[4][4]        = {0};
		float nepdMIPs[16]            = {0};
		float nepdMipsSumEast[5][16]  = {0};//[nMipmax(no,2,3,4,5)][EPD ring]
		float nepdMipsSumWest[5][16]  = {0};//[nMipmax(no,2,3,4,5)][EPD ring]
		float TileWeight[4][4] = {0};//[nMipMax = 2,3,4,5][EPDgroup]
		float TileWeight_sum[4][4] = {0};//[nMipMax = 2,3,4,5][EPDgroup]

		for(Int_t iHit=0; iHit<nepdHits; iHit++){
				epdHit = mPicoDst->epdHit(iHit);
				mip = epdHit->nMIP();
				int iring = epdHit->row() -1;//(1~16)-1 -> 0-15
				int ringgroup = mQvMaker -> RingGroup(iring);//0~3->0, 4~7->1, 8~11->2 12~15->3 0->most inner ring

				if( !epdHit) continue;
				Short_t side_EW = epdHit->side();// +1 for West  -1 for East


				//	StraightLine_center[iHit] = mEpdGeom->TileCenter(epdHit->id())        - picoEvent->primaryVertex();
				StraightLine_random[iHit] = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();

				phi_epd[iHit] = StraightLine_random[iHit].Phi();
				eta_epd[iHit] = StraightLine_random[iHit].Eta();




				//Fill TPC multiplicity vs EPD nMip ring by ring
				if(side_EW == -1)nepdMipsSumEast[0][iring] += mip;
				if(side_EW ==  1)nepdMipsSumWest[0][iring] += mip;
				for(int nmipmax=2;nmipmax<6;nmipmax++){//nMipmax = 2,3,4,5
						if(side_EW == -1) nepdMipsSumEast[nmipmax-1][iring] += (mip > nmipmax) ? nmipmax : mip;
						if(side_EW ==  1) nepdMipsSumWest[nmipmax-1][iring] += (mip > nmipmax) ? nmipmax : mip;
				}

				//sum epd nmip eash side
				if(side_EW == -1)nepdMIPsE += mip;
				if(side_EW ==  1)nepdMIPsW += mip;

				float WeightQvEpd = wV1Epd[centnumber] -> GetBinContent(iring+1);
				for(int nmipmax=0;nmipmax<4;nmipmax++){
						TileWeight[nmipmax][ringgroup] = (mip > nmipmax+2) ? nmipmax+2 : mip;
						if(QvWeight==1 && ringgroup != 0)TileWeight[nmipmax][ringgroup] *= WeightQvEpd;
						TileWeight_sum[nmipmax][ringgroup] += TileWeight[nmipmax][ringgroup];
				}

				//calculate Qvector
				for(int nmipmax=0;nmipmax<4;nmipmax++){
						QvEpdRaw[nmipmax][ringgroup][0] += TileWeight[nmipmax][ringgroup] * cos(phi_epd[iHit]);
						QvEpdRaw[nmipmax][ringgroup][1] += TileWeight[nmipmax][ringgroup] * sin(phi_epd[iHit]);
				}

		}//end of EPD roop


		//calculate EPD event plane
		for(int nMipMax=0;nMipMax<4;nMipMax++){
				for(int RingGroup=0;RingGroup<4;RingGroup++){
						for(int ixy=0;ixy<2;ixy++){
								if(TileWeight_sum[nMipMax][RingGroup] != 0)QvEpdRaw[nMipMax][RingGroup][ixy] /= sqrt(fabs(TileWeight_sum[nMipMax][RingGroup]));
								else QvEpdRaw[nMipMax][RingGroup][ixy] = -9999;
						}
						PsiEpdRaw[nMipMax][RingGroup] = atan2(QvEpdRaw[nMipMax][RingGroup][1],QvEpdRaw[nMipMax][RingGroup][0]);
						if(QvEpdRaw[nMipMax][RingGroup][0] != -9999 && QvEpdRaw[nMipMax][RingGroup][1] != -9999)hQvEpdRaw[centnumber][nMipMax][RingGroup] -> Fill(QvEpdRaw[nMipMax][RingGroup][0],QvEpdRaw[nMipMax][RingGroup][1]);
						if(QvEpdRaw[nMipMax][RingGroup][0] != -9999 && QvEpdRaw[nMipMax][RingGroup][1] != -9999)hPsiEpdRaw[centnumber][nMipMax][RingGroup] -> Fill(PsiEpdRaw[nMipMax][RingGroup]);
				}
		}

		for(int nmipmax=0;nmipmax<4;nmipmax++){
				for(int iring=0;iring<4;iring++){
						for(int ixy=0;ixy<2;ixy++){
								if(QvEpdRaw[nmipmax][iring][ixy] != -9999)pEpdQv[nmipmax][iring][ixy] -> Fill(centnumber,QvEpdRaw[nmipmax][iring][ixy]);
						}
				}
		}

		//for re-centering
		for(int nmip=0;nmip<4;nmip++){
				for(int iring=0;iring<4;iring++){
						for(int ixy=0;ixy<2;ixy++){
								QvEpdmean[nmip][iring][ixy]  = getEpd[nmip][iring][ixy] -> GetBinContent(centnumber+1);
								QvEpdsigma[nmip][iring][ixy] = getEpd[nmip][iring][ixy] -> GetBinError(centnumber+1);
								QvEpdRece[nmip][iring][ixy]  = (QvEpdRaw[nmip][iring][ixy] - QvEpdmean[nmip][iring][ixy])/QvEpdsigma[nmip][iring][ixy];
						}
						PsiEpdRece[nmip][iring] = atan2(QvEpdRece[nmip][iring][1],QvEpdRece[nmip][iring][0]);
						if(QvEpdRaw[nmip][iring][0] != -9999 && QvEpdRaw[nmip][iring][1] != -9999)hQvEpdRece[centnumber][nmip][iring]  -> Fill(QvEpdRece[nmip][iring][0],QvEpdRece[nmip][iring][1]);
						if(QvEpdRaw[nmip][iring][0] != -9999 && QvEpdRaw[nmip][iring][1] != -9999)hPsiEpdRece[centnumber][nmip][iring] -> Fill(PsiEpdRece[nmip][iring]);
				}
		}

		//for flattening
		for(int nfill=1;nfill<9;nfill++){
				for(int nmip=0;nmip<4;nmip++){
						for(int iring=0;iring<4;iring++){
								if(QvEpdRaw[nmip][iring][0] != -9999 && QvEpdRaw[nmip][iring][1] != -9999)pEPD_shift_sin[nmip][iring]->Fill(nfill-1,centnumber,sin(nfill*PsiEpdRece[nmip][iring]));
								if(QvEpdRaw[nmip][iring][0] != -9999 && QvEpdRaw[nmip][iring][1] != -9999)pEPD_shift_cos[nmip][iring]->Fill(nfill-1,centnumber,cos(nfill*PsiEpdRece[nmip][iring]));
						}
				}
		}
		float mEPD_shift_sin[8][4][4];
		float mEPD_shift_cos[8][4][4];
		for(int n=0;n<8;n++){
				for(int nmip=0;nmip<4;nmip++){
						for(int iring=0;iring<4;iring++){
								mEPD_shift_sin[n][nmip][iring] = mEPD_shift_s[nmip][iring]->GetBinContent(n+1, centnumber+1);
								mEPD_shift_cos[n][nmip][iring] = mEPD_shift_c[nmip][iring]->GetBinContent(n+1, centnumber+1);
						}
				}
		}
		for(int nmip=0;nmip<4;nmip++){
				for(int iring=0;iring<4;iring++){
						PsiEpdFlat[nmip][iring] = PsiEpdRece[nmip][iring];
						for(int i=1;i<9;i++){
								PsiEpdFlat[nmip][iring] += 2*(-mEPD_shift_sin[i-1][nmip][iring]*cos(i*PsiEpdRece[nmip][iring])+mEPD_shift_cos[i-1][nmip][iring]*sin(i*PsiEpdRece[nmip][iring]))/(Float_t)(i);
						}
						if(QvEpdRaw[nmip][iring][0] != -9999 && QvEpdRaw[nmip][iring][1] != -9999)hPsiEpdFlat[centnumber][nmip][iring] -> Fill(PsiEpdFlat[nmip][iring]);
				}
		}



		const int EPcorr_epdA[6]   = {0,0,0,1,1,2};
		const int EPcorr_epdB[6]   = {1,2,3,2,3,3};
		const int EPcorr_tpcA[3]   = {0,0,1};
		const int EPcorr_tpcB[3]   = {1,2,2};
		const int EPcorr_evstA[12] = {0,0,0,1,1,1,2,2,2,3,3,3};
		const int EPcorr_evstB[12] = {0,1,2,0,1,2,0,1,2,0,1,2};

		//Event plane correlation
		for(int tpc=0;tpc<3;tpc++){
				if(QvTpcRaw[EPcorr_tpcA[tpc]][0] != -9999 && QvTpcRaw[EPcorr_tpcA[tpc]][1] != -9999  && QvTpcRaw[EPcorr_tpcB[tpc]][0] != -9999 && QvTpcRaw[EPcorr_tpcB[tpc]][1] != -9999) {
						pEPcorr_TPC[tpc]    -> Fill(centnumber,cos(PsiTpcFlat[EPcorr_tpcA[tpc]]-PsiTpcFlat[EPcorr_tpcB[tpc]]));
						pEPcorrsin_TPC[tpc] -> Fill(centnumber,sin(PsiTpcFlat[EPcorr_tpcA[tpc]]-PsiTpcFlat[EPcorr_tpcB[tpc]]));
				}
		}
		for(int nmip=0;nmip<4;nmip++){
				for(int epd=0;epd<6;epd++){
						if(QvEpdRaw[nmip][EPcorr_epdA[epd]][0] != -9999 && QvEpdRaw[nmip][EPcorr_epdA[epd]][1] != -9999 && QvEpdRaw[nmip][EPcorr_epdB[epd]][0] != -9999 && QvEpdRaw[nmip][EPcorr_epdB[epd]][1] != -9999){
								pEPcorr_EPD[epd][nmip]    -> Fill(centnumber,cos(PsiEpdFlat[nmip][EPcorr_epdA[epd]]-PsiEpdFlat[nmip][EPcorr_epdB[epd]]));
								pEPcorrsin_EPD[epd][nmip] -> Fill(centnumber,sin(PsiEpdFlat[nmip][EPcorr_epdA[epd]]-PsiEpdFlat[nmip][EPcorr_epdB[epd]]));
						}
				}
				for(int evst=0;evst<12;evst++){
						if(QvEpdRaw[nmip][EPcorr_evstA[evst]][0] != -9999 && QvEpdRaw[nmip][EPcorr_evstA[evst]][1] != -9999 && QvTpcRaw[EPcorr_evstB[evst]][0] != -9999 && QvTpcRaw[EPcorr_evstB[evst]][1] != -9999){
								pEPcorr_EPDvsTPC[evst][nmip]    -> Fill(centnumber,cos(PsiEpdFlat[nmip][EPcorr_evstA[evst]]-PsiTpcFlat[EPcorr_evstB[evst]]));
								pEPcorrsin_EPDvsTPC[evst][nmip] -> Fill(centnumber,sin(PsiEpdFlat[nmip][EPcorr_evstA[evst]]-PsiTpcFlat[EPcorr_evstB[evst]]));
						}
				}
		}

		//Calculate v1 with EPD
		for(Int_t iHit=0; iHit<nepdHits; iHit++){
				epdHit = mPicoDst->epdHit(iHit);
				mip = epdHit->nMIP();
				int iring = epdHit->row() -1;//(1~16)-1 -> 0-15

				if( !epdHit) continue;
				Short_t side_EW = epdHit->side();// +1 for West  -1 for East


				//	StraightLine_center[iHit] = mEpdGeom->TileCenter(epdHit->id())        - picoEvent->primaryVertex();
				StraightLine_random[iHit] = mEpdGeom->RandomPointOnTile(epdHit->id()) - picoEvent->primaryVertex();

				phi_epd[iHit] = StraightLine_random[iHit].Phi();
				eta_epd[iHit] = StraightLine_random[iHit].Eta();

				if(mip > 5) mip = 5; //nMipMax = 5
				for(int psitpc=0;psitpc<3;psitpc++){
						if(QvTpcRaw[psitpc][0] != -9999 && QvTpcRaw[psitpc][1] != -9999 ) CalcV1EPD[psitpc][centnumber] -> Fill(eta_epd[iHit],cos(phi_epd[iHit]- PsiTpcFlat[psitpc]));
				}
				for(int psiepd=0;psiepd<4;psiepd++){
						if(QvEpdRaw[3][psiepd][0] != -9999 && QvEpdRaw[3][psiepd][1] != -9999 )CalcV1EPD[psiepd+3][centnumber] -> Fill(eta_epd[iHit],cos(phi_epd[iHit]- PsiEpdFlat[3][psiepd]));
				}

				if(QvTpcRaw[0][0] != -9999 && QvTpcRaw[0][1] != -9999) pV1RawFullRegion -> Fill(eta_epd[iHit], centnumber*10+5, cos(phi_epd[iHit] - PsiTpcFlat[0]));




		}//end of EPD roop


		//Calculate v1 with TPV
		// primary track loop for determine refmult ----------------------------------------------
		for (Int_t itr=0;itr<nTrack;itr++) { 
				const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);
				if(!ptrk)  continue;
				if(!ptrk->isPrimary())  continue;  // now selecting primary tracks

				const Float_t p_p = ptrk->pMom().Mag();
				const Float_t p_pt  = ptrk->pMom().Perp(); // zero for global tracks
				const Float_t p_phi = ptrk->pMom().Phi();  // zero for global track
				const Float_t p_eta = ptrk->pMom().PseudoRapidity();
				const Float_t p_pz = sqrt(p_p*p_p - p_pt*p_pt);
				const Float_t dca = ptrk->gDCA( pVertex ).Mag();
				const Int_t nHitsFit = ptrk->nHitsFit();
				const Int_t nHitsPoss = ptrk->nHitsMax();
				const Float_t nHitsR = (Float_t)nHitsFit/(Float_t)nHitsPoss;
				const Float_t Dedx = ptrk->dEdx();  // [keV/cm] 
				const Short_t Charge = ptrk->charge();
				const Float_t nSigmaPion = ptrk->nSigmaPion();
				const Float_t nSigmaProton = ptrk->nSigmaProton();
				const Int_t bTofPidTraitsIndex = ptrk->bTofPidTraitsIndex();
				const StThreeVectorF origin(ptrk->origin().X(),ptrk->origin().Y(),ptrk->origin().Z());
				const Float_t originx = origin.x();
				const Float_t originy = origin.y();
				const Float_t originz = origin.z();

				float beta = -1;
				float mm = 0;

				if( bTofPidTraitsIndex>=0 ) {
						StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
						beta = btofPidTraits->btofBeta();
						TVector3 tofPos = btofPidTraits->btofHitPos();
						const StThreeVectorF btofHitPos(tofPos.X(),tofPos.Y(),tofPos.Z());
						float L = tofPathLength( &origin, &btofHitPos, ptrk->helix(picoEvent->bField()).curvature())*0.01;  // cm -> m
						float tof = btofPidTraits->btof();
						if( tof>0.0 )  beta = L/(tof*(C_C_LIGHT/1.e9));
						else  beta = -1.0;
						mm = pow(p_p,2)*((1/pow(beta,2))-1);
				}

				if(dca > 3) continue;
				if(nHitsFit < 10) continue;
				if(nHitsR < 0.52) continue;

				int ptID  = mQvMaker -> PtId(p_pt);
				int pPionID = mMyPidMaker   -> PrimaryPion( bTofPidTraitsIndex, nSigmaPion, mm, Charge);
				int pProtonID = mMyPidMaker -> PrimaryProton( p_pt, bTofPidTraitsIndex, nSigmaProton, mm, Charge);
				hDedx -> Fill(p_p/Charge,Dedx);
				if( bTofPidTraitsIndex>=0 ) hMM   -> Fill(p_p/Charge,mm);
				if(ptID != -1){
						if(pPionID==0   || pPionID==1){
								hprimaryDedx_pion -> Fill(p_p/Charge,Dedx);
								if( bTofPidTraitsIndex>=0 )hprimaryMM_pion   -> Fill(p_p/Charge,mm);
								for(int psitpc=0;psitpc<3;psitpc++){
										if(QvTpcRaw[psitpc][0]!=-9999 && QvTpcRaw[psitpc][1]!=-9999) CalcV1PionTPC[psitpc][centnumber][ptID][pPionID]      -> Fill(p_eta,cos(p_phi - PsiTpcFlat[psitpc]));
								}
								for(int psiepd=0;psiepd<4;psiepd++){
										if(QvEpdRaw[3][psiepd][0]!=-9999 && QvEpdRaw[3][psiepd][1]!=-9999)CalcV1PionTPC[psiepd+3][centnumber][ptID][pPionID]     -> Fill(p_eta,cos(p_phi - PsiEpdFlat[3][psiepd]));
								}
						}
						if(pProtonID==0 || pProtonID==1){
								hprimaryDedx_proton -> Fill(p_p/Charge,Dedx);
								if( bTofPidTraitsIndex>=0 )hprimaryMM_proton   -> Fill(p_p/Charge,mm);
								for(int psitpc=0;psitpc<3;psitpc++){
										if(QvTpcRaw[psitpc][0]!=-9999 && QvTpcRaw[psitpc][1]!=-9999) CalcV1ProtonTPC[psitpc][centnumber][ptID][pProtonID] -> Fill(p_eta,cos(p_phi - PsiTpcFlat[psitpc]));
								}
								for(int psiepd=0;psiepd<4;psiepd++){
										if(QvEpdRaw[3][psiepd][0]!=-9999 && QvEpdRaw[3][psiepd][1]!=-9999)CalcV1ProtonTPC[psiepd+3][centnumber][ptID][pProtonID] -> Fill(p_eta,cos(p_phi - PsiEpdFlat[3][psiepd]));
								}
						}
				}
				if(QvEpdRaw[3][0][0] != -9999 && QvEpdRaw[3][0][1] != -9999) pV1RawFullRegion -> Fill(p_eta, centnumber*10+5, cos(p_phi - PsiEpdFlat[3][0]));


		}

		int CountPion[2]   = {0};
		int CountProton[2] = {0};
		//global track for PID of proton & pion
		for (Int_t itr=0;itr<nTrack;itr++) { 
				const StPicoTrack *ptrk = (StPicoTrack*)mPicoDst->track(itr);

				if(!isGoodTrack(ptrk)) continue;

				const Float_t p = ptrk->gMom().Mag();
				const Float_t pt = ptrk->gMom().Perp();
				const Float_t phi = ptrk->gMom().Phi();
				const Float_t eta = ptrk->gMom().PseudoRapidity();
				const Float_t px = ptrk->gMom().x();
				const Float_t py = ptrk->gMom().y();
				const Float_t pz = ptrk->gMom().z();
				const Float_t dca = ptrk->gDCA( pVertex ).Mag();
				const Int_t nHits = ptrk->nHits();
				const Int_t nHitsFit = ptrk->nHitsFit();
				const Int_t nHitsPoss = ptrk->nHitsMax();
				const Float_t nHitsDedx = ptrk->nHitsDedx();
				const Float_t nHitsR = (Float_t)nHitsFit/(Float_t)nHitsPoss;
				const Float_t Dedx = ptrk->dEdx();  // [keV/cm] 
				const Short_t Charge = ptrk->charge();
				const Float_t nSigmaPion = ptrk->nSigmaPion();
				const Float_t nSigmaProton = ptrk->nSigmaProton();
				const Int_t bTofPidTraitsIndex = ptrk->bTofPidTraitsIndex();
				const StThreeVectorF origin(ptrk->origin().X(),ptrk->origin().Y(),ptrk->origin().Z());
				const Float_t originx = origin.x();
				const Float_t originy = origin.y();
				const Float_t originz = origin.z();

				float beta = -1;
				float mm = 0;
				float tofphi = 0;

				if( bTofPidTraitsIndex>=0 ) {
						StPicoBTofPidTraits *btofPidTraits = mPicoDst->btofPidTraits(bTofPidTraitsIndex);
						beta = btofPidTraits->btofBeta();
						TVector3 tofPos = btofPidTraits->btofHitPos();
						const StThreeVectorF btofHitPos(tofPos.X(),tofPos.Y(),tofPos.Z());
						float L = tofPathLength( &origin, &btofHitPos, ptrk->helix(picoEvent->bField()).curvature())*0.01;  // cm -> m
						float tof = btofPidTraits->btof();
						float TOFx = tofPos.X();
						float TOFy = tofPos.Y();
						float TOFz = tofPos.Z();
						tofphi = atan2(TOFy,TOFx);
						hTPCphivsTOFphi -> Fill(phi,tofphi);
						hDeltaphiTpcTof -> Fill(tofphi - phi);
						hTofXvsY -> Fill(TOFx,TOFy);
						hTofXvsZ -> Fill(TOFx,TOFz);
						hTofYvsZ -> Fill(TOFy,TOFz);
						float dPhi = atan2(sin(tofphi-phi), cos(tofphi-phi));
						TVector3 TofTpcVector(tofPos.X()-pVertex.x(), tofPos.Y()-pVertex.y(), tofPos.Z()-pVertex.z());
						float CosTheta = TofTpcVector.Z()/(sqrt(TofTpcVector.X()*TofTpcVector.X() + TofTpcVector.Y()*TofTpcVector.Y() + TofTpcVector.Z()*TofTpcVector.Z()));
						float tofeta = -0.5*TMath::Log( (1.0-CosTheta)/(1.0+CosTheta) );
						if(Charge<0)hDeltaPhivsPt[0]->Fill(dPhi, pt);
						if(Charge>0)hDeltaPhivsPt[1]->Fill(dPhi, pt);
						if(Charge<0)hDeltaEtavsPt[0]->Fill(tofeta-eta, pt);
						if(Charge>0)hDeltaEtavsPt[1]->Fill(tofeta-eta, pt);
						if( tof>0.0 )  beta = L/(tof*(C_C_LIGHT/1.e9));

						else  beta = -1.0;
						mm = pow(p,2)*((1/pow(beta,2))-1);
						tofmatch++;
				}
				int pionID   = mMyPidMaker -> DaughterPion(bTofPidTraitsIndex, nSigmaPion, mm, Charge);
				int protonID = mMyPidMaker -> DaughterProton(bTofPidTraitsIndex, nSigmaProton, mm, Charge);
				if(pionID != -1){
						float energypi = sqrt(p*p*+mpi*mpi);
						float rapiditypi = log(fabs((energypi+pz)/(energypi-pz)))/2;
						//cout << "rapiditypi = " << rapiditypi << "  Energy = " << energypi << "  p = " << p  << "  px = " << px << "   py = " << py << "  pz = " << pz <<  "  log() = " << (energypi+pz)/(energypi-pz) << endl;
						pxpionpool[pionID][CountPion[pionID]]       = px;
						pypionpool[pionID][CountPion[pionID]]       = py;
						pzpionpool[pionID][CountPion[pionID]]       = pz;
						originxpionpool[pionID][CountPion[pionID]]  = originx; 
						originypionpool[pionID][CountPion[pionID]]  = originy; 
						originzpionpool[pionID][CountPion[pionID]]  = originz; 
						dcapionpool[pionID][CountPion[pionID]]      = dca;
						chargepionpool[pionID][CountPion[pionID]]   = Charge;
						etapionpool[pionID][CountPion[pionID]]      = eta;
						phipionpool[pionID][CountPion[pionID]]      = phi;
						rapiditypionpool[pionID][CountPion[pionID]] = rapiditypi;
						deltaphipionpool[pionID][CountPion[pionID]] = tofphi - phi; 
						if(bTofPidTraitsIndex>=0 ) tofflagpionpool[pionID][CountPion[pionID]] = true;
						else                       tofflagpionpool[pionID][CountPion[pionID]] = false;
						CountPion[pionID]++;
				}
				if(protonID != -1){
						float energyp = sqrt(p*p*+mp*mp);
						float rapidityp = log(fabs((energyp+pz)/(energyp-pz)))/2;
						//cout << "rapiditypro = " << rapidityp << "  Energy = " << energyp << "  pz = " << pz  <<  "  log() = " << (energyp+pz)/(energyp-pz) << endl;
						pxprotonpool[protonID][CountProton[protonID]]       = px;
						pyprotonpool[protonID][CountProton[protonID]]       = py;
						pzprotonpool[protonID][CountProton[protonID]]       = pz;
						originxprotonpool[protonID][CountProton[protonID]]  = originx; 
						originyprotonpool[protonID][CountProton[protonID]]  = originy; 
						originzprotonpool[protonID][CountProton[protonID]]  = originz; 
						dcaprotonpool[protonID][CountProton[protonID]]      = dca;
						chargeprotonpool[protonID][CountProton[protonID]]   = Charge;
						etaprotonpool[protonID][CountProton[protonID]]      = eta;
						phiprotonpool[protonID][CountProton[protonID]]      = phi;
						tofphiprotonpool[protonID][CountProton[protonID]]   = tofphi;
						rapidityprotonpool[protonID][CountProton[protonID]] = rapidityp;
						deltaphiprotonpool[protonID][CountProton[protonID]] = tofphi - phi; 
						if(bTofPidTraitsIndex>=0 ) tofflagprotonpool[protonID][CountProton[protonID]] = true;
						else                       tofflagprotonpool[protonID][CountProton[protonID]] = false;
						CountProton[protonID]++;
				}
				//hDedx -> Fill(p/Charge,Dedx);
				if(pionID == 0   || pionID ==1)    hDedx_pion   -> Fill(p/Charge,Dedx);
				if(protonID == 0 || protonID == 1) hDedx_proton -> Fill(p/Charge,Dedx);
				if( bTofPidTraitsIndex>=0 ) {
						//hMM   -> Fill(p/Charge,mm);
						if(pionID == 0   || pionID ==1)    hMM_pion  -> Fill(p/Charge,mm);
						if(protonID == 0 || protonID ==1)  hMM_proton  -> Fill(p/Charge,mm);
				}

		}//end of global track loop

		//for identify Lambda
		for(int proton=0;proton<CountProton[0];proton++){

				StThreeVectorF pmom(pxprotonpool[0][proton],pyprotonpool[0][proton],pzprotonpool[0][proton]); //momentum of proton
				StThreeVectorF porigin(originxprotonpool[0][proton],originyprotonpool[0][proton],originzprotonpool[0][proton]);//origin of proton

				for(int pion=0;pion<CountPion[0];pion++){

						StThreeVectorF pimom(pxpionpool[0][pion],pypionpool[0][pion],pzpionpool[0][pion]); //momentum of pion
						StThreeVectorF piorigin(originxpionpool[0][pion],originypionpool[0][pion],originzpionpool[0][pion]);//origin of pion

						StPhysicalHelixD helix_pos(pmom,porigin, BField*kilogauss, chargeprotonpool[0][proton]);
						StPhysicalHelixD helix_neg(pimom,piorigin, BField*kilogauss, chargepionpool[0][pion]);

						const double res = 3.0;
						float xc_pos = helix_pos.xcenter();
						float yc_pos = helix_pos.ycenter();
						float xc_neg = helix_neg.xcenter();
						float yc_neg = helix_neg.ycenter();
						float dd = sqrt((xc_pos-xc_neg)*(xc_pos-xc_neg) + (yc_pos-yc_neg)*(yc_pos-yc_neg));
						float r_pos = 1./helix_pos.curvature();
						float r_neg = 1./helix_neg.curvature();

						if(dd < fabs(r_pos-r_neg)-res || r_pos+r_neg+res < dd) continue;

						//daughter particle information
						pair<double,double> s = helix_pos.pathLengths(helix_neg);
						StThreeVectorF mom_pos = helix_pos.momentumAt(s.first,BField*kilogauss);
						StThreeVectorF mom_neg = helix_neg.momentumAt(s.second,BField*kilogauss);
						float ep = sqrt(mom_pos*mom_pos + mp*mp);//energy of proton
						float epi = sqrt(mom_neg*mom_neg + mpi*mpi);//energy of pion-

						//parent particle information
						StThreeVectorF mom_v0 = mom_pos + mom_neg; //momentum of Lambda
						float eV0 = ep + epi;//energy of Lambda
						float invM = sqrt(eV0*eV0 - mom_v0.mag()*mom_v0.mag());

						//  skip particle too large invM
						if(invM > 1.16) continue;
						hinvMLambda_topolo[0][0] -> Fill(invM);

						bool isGoodppiDCA = mMyPidMaker -> GoodppiDCA(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion]);
						//if(!isGoodppiDCA) continue;
						//hinvMLambda_topolo[0][1] -> Fill(invM);

						//topological cut information
						StThreeVectorF dca_pos = helix_pos.at(s.first);
						StThreeVectorF dca_neg = helix_neg.at(s.second);
						float dcaDaughters = (dca_pos - dca_neg).mag(); //p-pi DCA
						StThreeVectorF V0 = (dca_pos+dca_neg)*0.5;

						bool isGoodDaughterDCA = mMyPidMaker -> GoodDaughterDCA(0,centnumber, dcaDaughters);
						//if(!isGoodDaughterDCA) continue;
						//hinvMLambda_topolo[0][2] -> Fill(invM);

						StThreeVectorF v0toPV = V0 - PrimaryVertex;
						float angle   = (v0toPV).angle(mom_v0);
						float decayl  = (v0toPV).mag();//decay length
						float dca2vtx = (v0toPV).mag()*TMath::Sin(angle); //Lambda DCA

						bool isGoodLambdaDCA = mMyPidMaker -> GoodLambdaDCA(0,centnumber, dca2vtx);
						bool isGoodDecayLength = mMyPidMaker -> GoodDecayLength(0,centnumber, decayl);
						//if(!isGoodLambdaDCA) continue;
						//hinvMLambda_topolo[0][3] -> Fill(invM);
						//if(!isGoodDecayLength) continue;
						//hinvMLambda_topolo[0][4] -> Fill(invM);

						//Change Topological cut for finding optimal cut 
						//for(int i=0;i<10;i++){
						//		bool isGoodLambdaProtonDCA = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton]+0.2-i*0.1, dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl);//0.3~1.2
						//		if(isGoodLambdaProtonDCA)hinvMLambdaTopolo[0][centnumber][i] -> Fill(invM);
						//		bool isGoodLambdaPionDCA = mMyPidMaker -> ChangeTopological(0, centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion]+0.4-i*0.1, dcaDaughters, dca2vtx, decayl);//0.9~1.8
						//		if(isGoodLambdaPionDCA)hinvMLambdaTopolo[0][centnumber][i+10] -> Fill(invM);
						//		bool isGoodLambdaDaughterDCA = mMyPidMaker -> ChangeTopological(0, centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters-0.5+i*0.1, dca2vtx, decayl);//0.3~1.2
						//		if(isGoodLambdaDaughterDCA)hinvMLambdaTopolo[0][centnumber][i+20] -> Fill(invM);
						//		bool isGoodLambdaLambdaDCA = mMyPidMaker -> ChangeTopological(0, centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx-0.5+i*0.1, decayl);//0.3~1.2
						//		if(isGoodLambdaLambdaDCA)hinvMLambdaTopolo[0][centnumber][i+30] -> Fill(invM);
						//		bool isGoodLambdaDecayL= mMyPidMaker -> ChangeTopological(0, centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl+3.0-i*0.5);//5.0~9.5
						//		if(isGoodLambdaDecayL)hinvMLambdaTopolo[0][centnumber][i+40] -> Fill(invM);
						//}


						if(!isGoodppiDCA) continue;
						hinvMLambda_topolo[0][1] -> Fill(invM);
						if(!isGoodDaughterDCA) continue;
						hinvMLambda_topolo[0][2] -> Fill(invM);
						if(!isGoodLambdaDCA) continue;
						hinvMLambda_topolo[0][3] -> Fill(invM);

						TLorentzVector mom4d_v0( mom_v0.x(), mom_v0.y(), mom_v0.z(), eV0);
						float etav0 = mom4d_v0.Eta();
						float ptv0  = mom4d_v0.Pt();
						float phiv0 = mom4d_v0.Phi();
						float pzv0  = mom4d_v0.Pz();
						float momv0 = mom_v0.mag();
						float energyv0   = sqrt(momv0*momv0+Lambdamass*Lambdamass);
						float rapidityv0 = log((energyv0+pzv0)/(energyv0-pzv0))/2;
						int MassDecayBin = mMyPidMaker -> GetMassGroup(invM);
						hDecayLvspt[MassDecayBin] -> Fill(decayl,ptv0);

						if(!isGoodDecayLength) continue;
						hinvMLambda_topolo[0][4] -> Fill(invM);

						if(1.11 < invM && invM < 1.121){
								hDaughterPionRapidityvsPt[0]   -> Fill(rapiditypionpool[0][pion],sqrt(pxpionpool[0][pion]*pxpionpool[0][pion]+pypionpool[0][pion]*pypionpool[0][pion]));
								hDaughterProtonRapidityvsPt[0] -> Fill(rapidityprotonpool[0][proton],sqrt(pxprotonpool[0][proton]*pxprotonpool[0][proton]+pyprotonpool[0][proton]*pyprotonpool[0][proton]));
								hDaughterPionPhi[0]        -> Fill(phipionpool[0][pion]);
								hDaughterProtonPhi[0]      -> Fill(phiprotonpool[0][proton]);
								if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999) hDaughterPionDeltaPhi[0]   -> Fill(phipionpool[0][pion] - PsiEpdFlat[nMipMaxCut-2][0]);
								if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999) hDaughterProtonDeltaPhi[0] -> Fill(phiprotonpool[0][pion] - PsiEpdFlat[nMipMaxCut-2][0]);
								hLambdaEta[0]              -> Fill(etav0);
								hLambdaPt[0]               -> Fill(ptv0);
								hLambdaDeltaPhi[0]         -> Fill(phiv0 - PsiEpdFlat[nMipMaxCut-2][0]);
								hLambdaPhi[0]              -> Fill(phiv0);
								hLambdaRapidity[0]         -> Fill(rapidityv0);
								hLambdaRapidityvsPt[0]     -> Fill(rapidityv0,ptv0);
						}
						//phi flatting
						float phiave[3] = {0};
						float phiget[3][30] = {0};//0-> Lambda 1-> DaughterPion 2 -> DaughterProton
						float phiflat[3] = {0};
						int MassBinForPhiFlat = mMyPidMaker -> GetMassBin(invM);
						for(int particle=0;particle<3;particle++){
								for(int massbin=0;massbin<30;massbin++){
										phiave[particle] += pinvMvsPhiforcorrection[particle] -> GetBinContent(massbin+1)/30;
										phiget[particle][massbin]= pinvMvsPhiforcorrection[particle] -> GetBinContent(massbin+1);
								}
						}
						phiflat[0] = phiv0*phiave[0]/phiget[0][MassBinForPhiFlat];
						phiflat[1] = phipionpool[0][pion]*phiave[1]/phiget[1][MassBinForPhiFlat];
						phiflat[2] = phiprotonpool[0][proton]*phiave[2]/phiget[2][MassBinForPhiFlat];
						//phi flatting for Daughter proton in RF
						const float Ave = -0.00807147;
						float WforPhiRF = Ave/pinvMvsPhiRFforcorrection -> GetBinContent(MassBinForPhiFlat+1);
						int RapidityBin = mMyPidMaker -> GetYforLambda(rapidityv0);
						int PtBin       = mMyPidMaker -> GetPtforLambda(ptv0);
						int MassBin     = mMyPidMaker -> MassClass(invM);
						int DividePtBin = mMyPidMaker -> DividePtBin(ptv0);

						//DaughterPion
						if(PtBin != -1)pInvMvsDaughterPionRapidity[0][PtBin]   -> Fill(invM,rapiditypionpool[0][pion]);
						if(RapidityBin != -1)pInvMvsDaughterPionPt[0][RapidityBin] -> Fill(invM,sqrt(pxpionpool[0][pion]*pxpionpool[0][pion] + pypionpool[0][pion]*pypionpool[0][pion]));
						if(PtBin != -1 && RapidityBin != -1)pInvMvsDaughterPionRapidityDivide[0][PtBin][RapidityBin]   -> Fill(invM,rapiditypionpool[0][pion]);
						if(RapidityBin != -1 && DividePtBin != -1)pInvMvsDaughterPionPtDivide[0][RapidityBin][DividePtBin] -> Fill(invM,sqrt(pxpionpool[0][pion]*pxpionpool[0][pion] + pypionpool[0][pion]*pypionpool[0][pion]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999)pInvMvsDaughterPionDeltaPhi[0]   -> Fill(invM,phipionpool[0][pion] - PsiEpdFlat[nMipMaxCut-2][0]);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && PtBin != -1 && RapidityBin != -1)pInvMvsDaughterPionDeltaPhiDivide[0][PtBin][RapidityBin]   -> Fill(invM,phipionpool[0][pion] - PsiEpdFlat[nMipMaxCut-2][0]);
						pInvMvsDaughterPionPhi[0]        -> Fill(invM,phipionpool[0][pion]);
						pInvMvsDaughterPionPhiFlat[0]    -> Fill(invM,phiflat[1]);
						pInvMvsDaughterPionMom[0]        -> Fill(invM,pimom.mag());
						//DaughterProton
						if(PtBin != -1)pInvMvsDaughterProtonRapidity[0][PtBin] -> Fill(invM,rapidityprotonpool[0][proton]);
						if(RapidityBin != -1)pInvMvsDaughterProtonPt[0][RapidityBin] -> Fill(invM,sqrt(pxprotonpool[0][proton]*pxprotonpool[0][proton] + pyprotonpool[0][proton]*pyprotonpool[0][proton]));
						if(PtBin != -1 && RapidityBin != -1)pInvMvsDaughterProtonRapidityDivide[0][PtBin][RapidityBin] -> Fill(invM,rapidityprotonpool[0][proton]);
						if(RapidityBin != -1 && DividePtBin != -1)pInvMvsDaughterProtonPtDivide[0][RapidityBin][DividePtBin] -> Fill(invM,sqrt(pxprotonpool[0][proton]*pxprotonpool[0][proton] + pyprotonpool[0][proton]*pyprotonpool[0][proton]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999)pInvMvsDaughterProtonDeltaPhi[0] -> Fill(invM,phiprotonpool[0][proton] - PsiEpdFlat[nMipMaxCut-2][0]);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && PtBin != -1 && RapidityBin != -1)pInvMvsDaughterProtonDeltaPhiDivide[0][PtBin][RapidityBin]   -> Fill(invM,phiprotonpool[0][proton] - PsiEpdFlat[nMipMaxCut-2][0]);
						pInvMvsDaughterProtonPhi[0]      -> Fill(invM,phiprotonpool[0][proton]);
						pInvMvsDaughterProtonPhiFlat[0]  -> Fill(invM,phiflat[2]);
						pInvMvsDaughterProtonMom[0]      -> Fill(invM,pmom.mag());
						//Lambda
						if(PtBin != -1)pInvMvsLambdaRapidity[0][PtBin] -> Fill(invM,rapidityv0);
						if(RapidityBin != -1)pInvMvsLambdaPt[0][RapidityBin] -> Fill(invM,ptv0);
						if(PtBin != -1 && RapidityBin != -1)pInvMvsLambdaRapidityDivide[0][PtBin][RapidityBin] -> Fill(invM,rapidityv0);
						if(RapidityBin != -1 && DividePtBin != -1)pInvMvsLambdaPtDivide[0][RapidityBin][DividePtBin] -> Fill(invM,ptv0);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999)pInvMvsLambdaDeltaPhi[0]         -> Fill(invM,phiv0 - PsiEpdFlat[nMipMaxCut-2][0]);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && PtBin != -1 && RapidityBin != -1)pInvMvsLambdaDeltaPhiDivide[0][PtBin][RapidityBin]   -> Fill(invM,phiv0 - PsiEpdFlat[nMipMaxCut-2][0]);
						pInvMvsLambdaPhi[0]              -> Fill(invM,phiv0);
						pInvMvsLambdaPhiFlat[0]          -> Fill(invM,phiflat[0]);
						pInvMvsLambdaMom[0]              -> Fill(invM,momv0);


						float Deltaphi        = mMyPidMaker -> GetDeltaPhi(phiv0, PsiEpdFlat[nMipMaxCut-2][0]);
						if(MassBin != -1)RapidityvsDeltaphi[0][MassBin] -> Fill(rapidityv0,Deltaphi);
						int RapidityCorrectv1 = mMyPidMaker -> GetYforCorrectv1(rapidityv0);
						int DeltaphiBin       = mMyPidMaker -> GetDeltaphiBin(Deltaphi);
						float Average = RapidityvsDeltaPhiBeforeCorrect[MassBin] -> Integral(RapidityCorrectv1+1,RapidityCorrectv1+1,1,16)/16;
						float Entry = RapidityvsDeltaPhiBeforeCorrect[MassBin] -> GetBinContent(RapidityvsDeltaPhiBeforeCorrect[MassBin]->GetBin(RapidityCorrectv1+1,DeltaphiBin+1));
						float Weightforv1correct = (Float_t)(Average/Entry);

						if(MassBin != -1)RapidityvsDeltaphiAfterCorrect[MassBin] -> Fill(rapidityv0,Deltaphi,Weightforv1correct);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsCospol[0][centnumber][RapidityBin][PtBin][3] -> Fill(invM,cos(Deltaphi));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsCos2pol[0][centnumber][RapidityBin][PtBin][3] -> Fill(invM,cos(2*(Deltaphi)));

						//Lambda rest frame
						TVector3 v0beta = mom4d_v0.BoostVector();
						TLorentzVector mom4d_p( TVector3( mom_pos.x(), mom_pos.y(), mom_pos.z() ), ep );
						TLorentzVector mom4d_pi( TVector3( mom_neg.x(), mom_neg.y(), mom_neg.z() ), epi );
						mom4d_p.Boost( -v0beta);
						mom4d_pi.Boost(-v0beta);
						float phip_V0RF = mom4d_p.Phi();
						float Rapidityp_V0RF = mom4d_p.Rapidity();
						float Ptp_V0RF = mom4d_p.Pt();
						float Momp_V0RF = mom4d_p.Rho();
						float phipi_V0RF = mom4d_pi.Phi();
						float Rapiditypi_V0RF = mom4d_pi.Rapidity();
						float Ptpi_V0RF = mom4d_pi.Pt();
						float Mompi_V0RF = mom4d_pi.Rho();
						float dphiEpd[4]={0};
						float dphiTpc[3]={0};
						if(1.11 < invM && invM < 1.121){
								hDaughterPionRapidityvsPtRF[0]   -> Fill(Rapiditypi_V0RF,Ptpi_V0RF);
								hDaughterProtonRapidityvsPtRF[0] -> Fill(Rapidityp_V0RF,Ptp_V0RF);
								hDaughterPionPhiRF[0]    -> Fill(phipi_V0RF);
								hDaughterProtonPhiRF[0]  -> Fill(phip_V0RF);
						}
						hInvMassvsPhistar -> Fill(invM,phip_V0RF);
						int MassBinforPhiStarCorrection = mMyPidMaker -> GetMassBinforPhiStarCorrection(invM);
						int PhistarBin = mMyPidMaker -> GetPhiStarBin(phip_V0RF);
						float AveragePhiStar = pinvMvsPhiStarCorrection -> Integral(MassBinforPhiStarCorrection+1,MassBinforPhiStarCorrection+1,1,100)/100;
						float EntryPhiStar   = pinvMvsPhiStarCorrection -> GetBinContent(pinvMvsPhiStarCorrection -> GetBin(MassBinforPhiStarCorrection+1,PhistarBin+1));
						float WeightforPhiStar = AveragePhiStar/EntryPhiStar;
						//cout << "invM = " << invM << "    MassBin = " << MassBinforPhiStarCorrection << "   Ave = " << AveragePhiStar << "   Entry = " << EntryPhiStar << endl;
						if(MassBinforPhiStarCorrection != -1 && PhistarBin != -1)hInvMassvsPhistarAfterCorrection -> Fill(invM,phip_V0RF,WeightforPhiStar);
						//cout << "Boost  proton  phi = " << phip_V0RF << "   rapidity = " << Rapidityp_V0RF << "   pt = " << Ptp_V0RF << "   mom = " << Momp_V0RF << endl;
						//cout << "Boost  pion  phi = " << phipi_V0RF << "   rapidity = " << Rapiditypi_V0RF << "   pt = " << Ptpi_V0RF << "   mom = " << Mompi_V0RF << endl;
						int RapidityBinProtonRF = mMyPidMaker -> GetYforLambda(Rapidityp_V0RF);
						int PtBinProtonRF       = mMyPidMaker -> GetPtforLambda(Ptp_V0RF); 
						int RapidityBinPionRF = mMyPidMaker -> GetYforLambda(Rapiditypi_V0RF);
						int PtBinPionRF       = mMyPidMaker -> GetPtforLambda(Ptpi_V0RF); 
						if(PtBinProtonRF != -1)pInvMvsDaughterProtonRapidityRF[0][PtBinProtonRF] -> Fill(invM,Rapidityp_V0RF);
						if(RapidityBinProtonRF != -1)pInvMvsDaughterProtonPtRF[0][RapidityBinProtonRF] -> Fill(invM,Ptp_V0RF);
						pInvMvsDaughterPionPhiRF[0]        -> Fill(invM,phipi_V0RF);
						pInvMvsDaughterPionMomRF[0]        -> Fill(invM,Mompi_V0RF);
						if(PtBinPionRF != -1)pInvMvsDaughterPionRapidityRF[0][PtBinPionRF] -> Fill(invM,Rapiditypi_V0RF);
						if(RapidityBinPionRF != -1)pInvMvsDaughterPionPtRF[0][RapidityBinPionRF] -> Fill(invM,Ptpi_V0RF);
						pInvMvsDaughterProtonPhiRF[0]      -> Fill(invM,phip_V0RF);
						pInvMvsDaughterProtonMomRF[0]      -> Fill(invM,Momp_V0RF);
						pInvMvsDaughterProtonPhiRFFlat[0]  -> Fill(invM,phip_V0RF*WforPhiRF);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999)pInvMvsDaughterProtonDeltaPhiRF[0]     -> Fill(invM,PsiEpdFlat[nMipMaxCut-2][0] - phip_V0RF);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999)pInvMvsDaughterProtonDeltaPhiRFFlat[0] -> Fill(invM,PsiEpdFlat[nMipMaxCut-2][0] - phip_V0RF*WforPhiRF);

						//for(int etaId=0;etaId<3;etaId++){
						//		bool AutoCorrFlag = mMyPidMaker -> RemoveAutoCorr(etaId, RapidityBin, PtBin, etav0);
						//		dphiTpc[etaId] = PsiTpcFlat[etaId] - phip_V0RF;
						//		if(QvTpcRaw[etaId][0] != -9999 && QvTpcRaw[etaId][1] != -9999 && AutoCorrFlag)hinvMLambda[0][centnumber][RapidityBin][PtBin][etaId] -> Fill(invM);
						//		if(QvTpcRaw[etaId][0] != -9999 && QvTpcRaw[etaId][1] != -9999 && AutoCorrFlag)pInvMassvsSinpol[0][centnumber][RapidityBin][PtBin][etaId] -> Fill(invM,sin(dphiTpc[etaId]));
						//}

						//for(int iring=0;iring<4;iring++){
						//		dphiEpd[iring] = PsiEpdFlat[nMipMaxCut-2][iring] - phip_V0RF;
						//		if(QvEpdRaw[nMipMaxCut-2][iring][0] != -9999 && QvEpdRaw[nMipMaxCut-2][iring][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambda[0][centnumber][RapidityBin][PtBin][iring+3] -> Fill(invM);
						//		if(QvEpdRaw[nMipMaxCut-2][iring][0] != -9999 && QvEpdRaw[nMipMaxCut-2][iring][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpol[0][centnumber][RapidityBin][PtBin][iring+3] -> Fill(invM,sin(dphiEpd[iring]));
						//}
						dphiEpd[0] = PsiEpdFlat[nMipMaxCut-2][0] - phip_V0RF;
						int DeltaPhiBinforPol = mMyPidMaker -> GetDeltaPhiBinForPol(phiv0 - phip_V0RF);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambda[0][centnumber][RapidityBin][PtBin][DeltaPhiBinforPol] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && DeltaPhiBinforPol != -1)pInvMassvsSinpol[0][centnumber][RapidityBin][PtBin][DeltaPhiBinforPol] -> Fill(invM,sin(dphiEpd[0]));

						float dphiEpdPhiCorrect = PsiEpdFlat[nMipMaxCut-2][0] - phip_V0RF*WeightforPhiStar;
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && MassBinforPhiStarCorrection != -1 && PhistarBin != -1)pInvMassvsSinpolPhiCorrect[0][centnumber][RapidityBin][PtBin] -> Fill(invM,sin(dphiEpdPhiCorrect));
						if(RapidityBin != -1 && PtBin != -1 && PtBin != 10 && PtBin != 11)hinvMassforProductionPlane[RapidityBin][PtBin] -> Fill(invM);
						if(RapidityBin != -1 && PtBin != -1 && PtBin != 10 && PtBin != 11)pInvMassvsProductedPlanePol[RapidityBin][PtBin] -> Fill(invM,sin(phiv0-phip_V0RF));
						int DaughterPionYbin = mMyPidMaker -> DaughterRapidityBin(rapiditypionpool[0][pion]);
						int DaughterProtonYbin = mMyPidMaker -> DaughterRapidityBin(rapidityprotonpool[0][proton]);
						//rapidity acceptance cut
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1 && DaughterPionYbin != -1)pInvMvsSinpolDaughterCut[RapidityBin][0][DaughterPionYbin] -> Fill(invM,sin(dphiEpd[0]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1 && DaughterProtonYbin != -1)pInvMvsSinpolDaughterCut[RapidityBin][1][DaughterProtonYbin] -> Fill(invM,sin(dphiEpd[0]));

						float invMWidthCorrect = mMyPidMaker -> MassWidthCorrection(phiv0-phip_V0RF,invM);
						//cout << "DeltaPhi = " << phiv0-phip_V0RF << "   invM = " << invM << "   invM After Correct = " << invMWidthCorrect << endl;
						hInvMvsDeltaPhi -> Fill(invM,phiv0-phip_V0RF);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999) pInvMvsSinPol     -> Fill(invM,phiv0-phip_V0RF,-sin(dphiEpd[0]));
						hInvMvsDeltaPhiFull -> Fill(invM,phiv0-phip_V0RF);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999) pInvMvsSinPolFull     -> Fill(invM,phiv0-phip_V0RF,-sin(dphiEpd[0]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1)pDeltaPhivsSinPol[RapidityBin] -> Fill(phiv0-phip_V0RF,-sin(dphiEpd[0]));

						//Mass Width divided by delta Eta
						float dPhiCorrect = phiv0-phip_V0RF;
						if(-2*TMath::Pi() < dPhiCorrect && dPhiCorrect < -TMath::Pi()) dPhiCorrect += 2*TMath::Pi();
						if(TMath::Pi() < dPhiCorrect && dPhiCorrect < 2*TMath::Pi())    dPhiCorrect -= 2*TMath::Pi();
						int FilldPhiBin = mMyPidMaker -> FilldPhiBin(dPhiCorrect);
						int dEtaBin = mMyPidMaker -> DeltaEtaBin(etapionpool[0][pion],etaprotonpool[0][proton]);
						int PtBinForinvMdPhi = mMyPidMaker -> GetPtBinForInvMvsdPhi(ptv0);
						float dEta = etapionpool[0][pion]-etaprotonpool[0][proton];
						hInvMvsDeltaPhi_dEta[dEtaBin][PtBinForinvMdPhi] -> Fill(invM,dPhiCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999) pInvMvsSinPol_dEta[dEtaBin][PtBinForinvMdPhi]   -> Fill(invM,dPhiCorrect,-sin(dphiEpd[0]));
						if(1.08<invM && invM <1.16)  hdEtavsdPhi_invM[PtBinForinvMdPhi][0] -> Fill(dEta,dPhiCorrect);
						if(1.1<invM && invM <1.14)   hdEtavsdPhi_invM[PtBinForinvMdPhi][1] -> Fill(dEta,dPhiCorrect);
						if(1.11<invM && invM <1.122) hdEtavsdPhi_invM[PtBinForinvMdPhi][2] -> Fill(dEta,dPhiCorrect);
						if(1.112<invM && invM <1.12) hdEtavsdPhi_invM[PtBinForinvMdPhi][3] -> Fill(dEta,dPhiCorrect);
						if(1.114<invM && invM <1.118)hdEtavsdPhi_invM[PtBinForinvMdPhi][4] -> Fill(dEta,dPhiCorrect);
						int DeltaEtaBin = mMyPidMaker -> GetDeltaEtaRoughBin(etapionpool[0][pion],etaprotonpool[0][proton]);
						hInvMvsDeltaPhi_DivdEta[DeltaEtaBin][PtBinForinvMdPhi] -> Fill(invM,dPhiCorrect);
						hinvMvsDeltaPhi_dEta3D -> Fill(dPhiCorrect,dEta,invM);

						hInvMvsPhiProton[0]      -> Fill(invM, phiprotonpool[0][proton]);
						hInvMvsPhiPion[0]        -> Fill(invM, phipionpool[0][pion]);
						hInvMvsPhiLambda[0]      -> Fill(invM, phiv0);
						hInvMvsPhistar[0]        -> Fill(invM, phip_V0RF);
						hInvMvsdPsiPhiProton[0]  -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phiprotonpool[0][proton]),cos(PsiEpdFlat[nMipMaxCut-2][0]-phiprotonpool[0][proton])));
						hInvMvsdPsiPhiPion[0]    -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phipionpool[0][pion]),cos(PsiEpdFlat[nMipMaxCut-2][0]-phipionpool[0][pion])));
						hInvMvsdPsiPhiLambda[0]  -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phiv0),cos(PsiEpdFlat[nMipMaxCut-2][0]-phiv0)));
						hInvMvsdPsiPhistar[0]    -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phip_V0RF),cos(PsiEpdFlat[nMipMaxCut-2][0]-phip_V0RF)));
						hInvMvsDaughterdPhi[0]   -> Fill(invM, atan2(sin(phiprotonpool[0][proton]-phipionpool[0][pion]),cos(phiprotonpool[0][proton]-phipionpool[0][pion])));
						hInvMvsdEta[0]           -> Fill(invM, atan2(sin(etaprotonpool[0][proton]-etapionpool[0][pion]),cos(etaprotonpool[0][proton]-etapionpool[0][pion])));
						hDaughterdPhiLaboRest[0] -> Fill(atan2(sin(phiprotonpool[0][proton]-phipionpool[0][pion]),cos(phiprotonpool[0][proton]-phipionpool[0][pion])),atan2(sin(phip_V0RF-phiv0),cos(phip_V0RF-phiv0)));



						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && FilldPhiBin != -1)hInvMvsDeltaPhiNotWidthCorrect[centnumber][RapidityBin][PtBin][FilldPhiBin] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && FilldPhiBin != -1) pInvMvsSinPolNotWidthCorrect[centnumber][RapidityBin][PtBin][FilldPhiBin]  -> Fill(invM,-sin(dphiEpd[0]));
						//After Mass Width Correction
						hInvMvsDeltaPhiTurn -> Fill(invM,dPhiCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && FilldPhiBin != -1)hInvMvsDeltaPhiWidthCorrect[centnumber][RapidityBin][PtBin][FilldPhiBin] -> Fill(invMWidthCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && FilldPhiBin != -1) pInvMvsSinPolWidthCorrect[centnumber][RapidityBin][PtBin][FilldPhiBin]  -> Fill(invMWidthCorrect,-sin(dphiEpd[0]));
						if(RapidityBin != -1)hInvMvsDeltaPhiWidthCorrectYdivi[RapidityBin] -> Fill(invMWidthCorrect,dPhiCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1) pInvMvsSinPolWidthCorrectYdivi[RapidityBin]     -> Fill(invMWidthCorrect,dPhiCorrect,-sin(dphiEpd[0]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut][0][1] != -9999 && RapidityBin != -1) pInvMvsSinPolMassWcorrect[RapidityBin]     -> Fill(invMWidthCorrect,-sin(dphiEpd[0]));

						//Mass Width Correction
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambdaMassWidthCorrect[0][centnumber][RapidityBin][PtBin] -> Fill(invMWidthCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpolMassWidthCorrect[0][centnumber][RapidityBin][PtBin] -> Fill(invMWidthCorrect,sin(dphiEpd[0]));

						//Delta Eta cut
						if(fabs(deltaphipionpool[0][pion] - deltaphiprotonpool[0][proton])>0.2)
								if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambdadEtaCut[centnumber][RapidityBin][PtBin] -> Fill(invMWidthCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpoldEtaCut[centnumber][RapidityBin][PtBin] -> Fill(invMWidthCorrect,sin(dphiEpd[0]));


						int dPsiPhiBin = mMyPidMaker -> DeltaPsiLambdaPhi(PsiEpdFlat[nMipMaxCut-2][0], phiv0);//Psi phi
						//DeltaPsiPhi dependence
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && dPsiPhiBin != -1)hinvMLambdaDpsiphi[0][centnumber][RapidityBin][PtBin][dPsiPhiBin] -> Fill(invMWidthCorrect);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && dPsiPhiBin != -1)pInvMassvsSinpolDpsiphi[0][centnumber][RapidityBin][PtBin][dPsiPhiBin] -> Fill(invMWidthCorrect,sin(dphiEpd[0]));

						if(tofflagprotonpool[0][proton]){
								hDeltaPhi2D -> Fill(phiv0-phip_V0RF,deltaphiprotonpool[0][proton]);
								hTPCphivsTOFphiDproton -> Fill(phiprotonpool[0][proton],tofphiprotonpool[0][proton]);
								//cout << "phiv0 - phi* = " << phiv0-phip_V0RF <<"   phiTOF = " << tofphiprotonpool[0][proton] << "   phiTPC = " << phiprotonpool[0][proton] << "   phiTOF - phiTPC = " << deltaphiprotonpool[0][proton] << endl;
						}
						if(tofflagprotonpool[0][proton]&&tofflagpionpool[0][pion]){
								if(phiv0-phip_V0RF<0) hDeltaPhiPionvsProton[0] -> Fill(deltaphipionpool[0][pion],deltaphiprotonpool[0][proton]);
								if(phiv0-phip_V0RF>0) hDeltaPhiPionvsProton[1] -> Fill(deltaphipionpool[0][pion],deltaphiprotonpool[0][proton]);
						}
						hDeltaLambdaProtonPhi_LaboV0RF -> Fill(atan2(sin(phiv0-phiprotonpool[0][proton]),cos(phiv0-phiprotonpool[0][proton])), atan2(sin(phiv0-phip_V0RF), cos(phiv0-phip_V0RF)));
						hDeltaLambdaProtonPhivsPionPhi -> Fill(atan2(sin(phiv0-phiprotonpool[0][proton]),cos(phiv0-phiprotonpool[0][proton])), atan2(sin(phiv0-phipionpool[0][pion]),cos(phiv0-phipionpool[0][pion])));

						int EpMethodBin = mMyPidMaker -> GetdPsiphiBinforEPmethod(PsiEpdFlat[nMipMaxCut-2][0], phiv0);
						int DeltaPhiEPmethod = mMyPidMaker -> DeltaPhiforEPmethod(phiv0, phip_V0RF);
						//cout << "DeltaPsiPhi = " << atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phiv0),cos(PsiEpdFlat[nMipMaxCut-2][0]-phiv0)) << "  " << EpMethodBin << "  ||  " << "   DeltaPhi = " << atan2(sin(phiv0-phip_V0RF),cos(phiv0-phip_V0RF)) << "   " << DeltaPhiEPmethod << endl;

						//Event Plane method
						if(RapidityBin!=-1&&PtBin!=-1&&EpMethodBin!=-1&&DeltaPhiEPmethod!=-1){
								hInvMCentEPmethod[centnumber][EpMethodBin][DeltaPhiEPmethod] -> Fill(invM);
								if(centnumber==2||centnumber==3||centnumber==4||centnumber==5)hInvMRapEPmethod[RapidityBin][EpMethodBin][DeltaPhiEPmethod] -> Fill(invM);
								if(centnumber==2||centnumber==3||centnumber==4||centnumber==5)hInvMPtEPmethod[PtBin][EpMethodBin][DeltaPhiEPmethod]        -> Fill(invM);
						}

						//Change Topological cut for systematic
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1){
								bool isGoodLambdaProtonDCAminus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton]+0.1, dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl);//DCAproton >0.5 -> >0.4
								if(isGoodLambdaProtonDCAminus){
										hInvMCentTopoloSys[0][0][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][0][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][0][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][0][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][0][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][0][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaProtonDCAplus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton]-0.1, dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl);//DCAproton >0.5 -> >0.6
								if(isGoodLambdaProtonDCAplus){
										hInvMCentTopoloSys[0][1][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][1][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][1][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][1][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][1][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][1][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaPionDCAminus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion]+0.1, dcaDaughters, dca2vtx, decayl);//DCApion >1.7 -> >1.6
								if(isGoodLambdaPionDCAminus){
										hInvMCentTopoloSys[0][2][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][2][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][2][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][2][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][2][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][2][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaPionDCAplus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion]-0.1, dcaDaughters, dca2vtx, decayl);//DCApion >1.7 -> >1.8
								if(isGoodLambdaPionDCAplus){
										hInvMCentTopoloSys[0][3][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][3][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][3][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][3][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][3][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][3][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdappiDCAminus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters+0.1, dca2vtx, decayl);//Daughter DCA <1.1 -> <1.0
								if(isGoodLambdappiDCAminus){
										hInvMCentTopoloSys[0][4][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][4][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][4][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][4][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][4][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][4][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdappiDCAplus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters-0.1, dca2vtx, decayl);//Daughter DCA <1.1 -> <1.2
								if(isGoodLambdappiDCAplus){
										hInvMCentTopoloSys[0][5][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][5][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][5][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][5][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][5][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][5][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDCAminus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx+0.1, decayl);//Lambda DCA <0.8 -> <0.7
								if(isGoodLambdaDCAminus){
										hInvMCentTopoloSys[0][6][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][6][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][6][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][6][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][6][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][6][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDCAplus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx-0.1, decayl);//Lambda DCA <0.8 -> <0.9
								if(isGoodLambdaDCAplus){
										hInvMCentTopoloSys[0][7][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][7][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][7][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][7][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][7][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][7][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDecayLminus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl-0.5);//Decay L > 6.0 -> >5.5
								if(isGoodLambdaDecayLminus){
										hInvMCentTopoloSys[0][8][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][8][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][8][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][8][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][8][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][8][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDecayLplus = mMyPidMaker -> ChangeTopological(0,centnumber, dcaprotonpool[0][proton], dcapionpool[0][pion], dcaDaughters, dca2vtx, decayl+0.5);//Decay L > 6.0 -> >6.5
								if(isGoodLambdaDecayLplus){
										hInvMCentTopoloSys[0][9][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[0][9][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[0][9][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[0][9][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[0][9][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[0][9][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
						}// skip wrong event plane

				}

		}// end of Lambda loop



		//	for identify antiLambda

		for(int proton=0;proton<CountProton[1];proton++){

				StThreeVectorF pmom(pxprotonpool[1][proton],pyprotonpool[1][proton],pzprotonpool[1][proton]); //momentum of antiproton
				StThreeVectorF porigin(originxprotonpool[1][proton],originyprotonpool[1][proton],originzprotonpool[1][proton]);//origin of antiproton

				for(int pion=0;pion<CountPion[1];pion++){

						StThreeVectorF pimom(pxpionpool[1][pion],pypionpool[1][pion],pzpionpool[1][pion]); //momentum of pion
						StThreeVectorF piorigin(originxpionpool[1][pion],originypionpool[1][pion],originzpionpool[1][pion]);//origin of pion

						StPhysicalHelixD helix_pos(pimom,piorigin, BField*kilogauss, chargepionpool[1][proton]);
						StPhysicalHelixD helix_neg(pmom,porigin, BField*kilogauss, chargeprotonpool[1][pion]);

						const double res = 3.0;
						float xc_pos = helix_pos.xcenter();
						float yc_pos = helix_pos.ycenter();
						float xc_neg = helix_neg.xcenter();
						float yc_neg = helix_neg.ycenter();
						float dd = sqrt((xc_pos-xc_neg)*(xc_pos-xc_neg) + (yc_pos-yc_neg)*(yc_pos-yc_neg));
						float r_pos = 1./helix_pos.curvature();
						float r_neg = 1./helix_neg.curvature();

						if(dd < fabs(r_pos-r_neg)-res || r_pos+r_neg+res < dd) continue;

						//daughter particle information
						pair<double,double> s = helix_pos.pathLengths(helix_neg);
						StThreeVectorF mom_pos = helix_pos.momentumAt(s.first,BField*kilogauss);
						StThreeVectorF mom_neg = helix_neg.momentumAt(s.second,BField*kilogauss);
						float epi = sqrt(mom_pos*mom_pos + mpi*mpi);//energy of pion+
						float ep = sqrt(mom_neg*mom_neg + mp*mp);//energy of anti-proton

						//parent particle information
						StThreeVectorF mom_v0 = mom_pos + mom_neg; //momentum of Lambda
						float eV0 = ep + epi;//energy of Lambda
						float invM = sqrt(eV0*eV0 - mom_v0.mag()*mom_v0.mag());

						//  skip particle too large invM
						if(invM > 1.17) continue;
						hinvMLambda_topolo[1][0] -> Fill(invM);

						bool isGoodppiDCA = mMyPidMaker -> GoodppiDCA(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion]);
						//if(!isGoodppiDCA) continue;
						//hinvMLambda_topolo[1][1] -> Fill(invM);

						//topological cut information
						StThreeVectorF dca_pos = helix_pos.at(s.first);
						StThreeVectorF dca_neg = helix_neg.at(s.second);
						float dcaDaughters = (dca_pos - dca_neg).mag(); //p-pi DCA
						StThreeVectorF V0 = (dca_pos+dca_neg)*0.5;

						bool isGoodDaughterDCA = mMyPidMaker -> GoodDaughterDCA(1,centnumber, dcaDaughters);
						//if(!isGoodDaughterDCA) continue;
						//hinvMLambda_topolo[1][2] -> Fill(invM);

						StThreeVectorF v0toPV = V0 - PrimaryVertex;
						float angle   = (v0toPV).angle(mom_v0);
						float decayl  = (v0toPV).mag();//decay length
						float dca2vtx = (v0toPV).mag()*TMath::Sin(angle); //Lambda DCA

						bool isGoodLambdaDCA = mMyPidMaker -> GoodLambdaDCA(1,centnumber, dca2vtx);
						bool isGoodDecayLength = mMyPidMaker -> GoodDecayLength(1,centnumber, decayl);
						//if(!isGoodLambdaDCA) continue;
						//hinvMLambda_topolo[1][3] -> Fill(invM);
						//if(!isGoodDecayLength) continue;
						//hinvMLambda_topolo[1][4] -> Fill(invM);

						//Change Topological cut
						//for(int i=0;i<10;i++){
						//		bool isGoodLambdaProtonDCA = mMyPidMaker -> ChangeTopological(1, centnumber, dcaprotonpool[1][proton]+0.2-i*0.1, dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl);//0.3~1.2
						//		if(isGoodLambdaProtonDCA)hinvMLambdaTopolo[1][centnumber][i] -> Fill(invM);
						//		bool isGoodLambdaPionDCA = mMyPidMaker -> ChangeTopological(1, centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion]+0.4-i*0.1, dcaDaughters, dca2vtx, decayl);//1.6~2.5
						//		if(isGoodLambdaPionDCA)hinvMLambdaTopolo[1][centnumber][i+10] -> Fill(invM);
						//		bool isGoodLambdaDaughterDCA = mMyPidMaker -> ChangeTopological(1, centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters-0.5+i*0.1, dca2vtx, decayl);//0.3~1.2
						//		if(isGoodLambdaDaughterDCA)hinvMLambdaTopolo[1][centnumber][i+20] -> Fill(invM);
						//		bool isGoodLambdaLambdaDCA = mMyPidMaker -> ChangeTopological(1, centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx-0.5+i*0.1, decayl);//0.8~1.7
						//		if(isGoodLambdaLambdaDCA)hinvMLambdaTopolo[1][centnumber][i+30] -> Fill(invM);
						//		bool isGoodLambdaDecayL= mMyPidMaker -> ChangeTopological(1, centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl+3.0-i*0.5);//8.0~12.5
						//		if(isGoodLambdaDecayL)hinvMLambdaTopolo[1][centnumber][i+40] -> Fill(invM);
						//}
						if(!isGoodppiDCA) continue;
						hinvMLambda_topolo[1][1] -> Fill(invM);
						if(!isGoodDaughterDCA) continue;
						hinvMLambda_topolo[1][2] -> Fill(invM);
						if(!isGoodLambdaDCA) continue;
						hinvMLambda_topolo[1][3] -> Fill(invM);
						if(!isGoodDecayLength) continue;
						hinvMLambda_topolo[1][4] -> Fill(invM);


						TLorentzVector mom4d_v0( mom_v0.x(), mom_v0.y(), mom_v0.z(), eV0);
						float etav0 = mom4d_v0.Eta();
						float ptv0  = mom4d_v0.Pt();
						float phiv0 = mom4d_v0.Phi();
						float pzv0  = mom4d_v0.Pz();
						float momv0 = mom_v0.mag();
						float energyv0   = sqrt(momv0*momv0+Lambdamass*Lambdamass);
						float rapidityv0 = log((energyv0+pzv0)/(energyv0-pzv0))/2;

						if(1.11 < invM && invM < 1.121){
								hDaughterPionRapidityvsPt[1]   -> Fill(rapiditypionpool[1][pion],sqrt(pxpionpool[1][pion]*pxpionpool[1][pion]+pypionpool[1][pion]*pypionpool[1][pion]));
								hDaughterProtonRapidityvsPt[1] -> Fill(rapidityprotonpool[1][proton],sqrt(pxprotonpool[1][proton]*pxprotonpool[1][proton]+pyprotonpool[1][proton]*pyprotonpool[1][proton]));
								hDaughterPionPhi[1]    -> Fill(phipionpool[1][pion]);
								hDaughterProtonPhi[1]  -> Fill(phiprotonpool[1][proton]);
								hLambdaEta[1]          -> Fill(etav0);
								hLambdaPt[1]           -> Fill(ptv0);
								hLambdaPhi[1]          -> Fill(phiv0);
								hLambdaRapidity[1]     -> Fill(rapidityv0);
								hLambdaRapidityvsPt[1] -> Fill(rapidityv0,ptv0);
						}
						int RapidityBin = mMyPidMaker -> GetYforLambda(rapidityv0);
						int PtBin       = mMyPidMaker -> GetPtforLambda(ptv0);
						int MassBin     = mMyPidMaker -> MassClass(invM);
						if(PtBin != -1)pInvMvsDaughterPionRapidity[1][PtBin]   -> Fill(invM,rapiditypionpool[1][pion]);
						if(RapidityBin != -1)pInvMvsDaughterPionPt[1][RapidityBin] -> Fill(invM,sqrt(pxpionpool[1][pion]*pxpionpool[1][pion] + pypionpool[1][pion]*pypionpool[1][pion]));
						pInvMvsDaughterPionPhi[1]        -> Fill(invM,phipionpool[1][pion]);
						pInvMvsDaughterPionMom[1]        -> Fill(invM,pimom.mag());
						if(PtBin != -1)pInvMvsDaughterProtonRapidity[1][PtBin] -> Fill(invM,rapidityprotonpool[1][proton]);
						if(RapidityBin != -1)pInvMvsDaughterProtonPt[1][RapidityBin] -> Fill(invM,sqrt(pxprotonpool[1][proton]*pxprotonpool[1][proton] + pyprotonpool[1][proton]*pyprotonpool[1][proton]));
						pInvMvsDaughterProtonPhi[1]      -> Fill(invM,phiprotonpool[1][proton]);
						pInvMvsDaughterProtonMom[1]      -> Fill(invM,pmom.mag());
						if(PtBin != -1)pInvMvsLambdaRapidity[1][PtBin]  -> Fill(invM,rapidityv0);
						if(RapidityBin != -1)pInvMvsLambdaPt[1][RapidityBin]  -> Fill(invM,ptv0);
						pInvMvsLambdaPhi[1]              -> Fill(invM,phiv0);
						pInvMvsLambdaMom[1]              -> Fill(invM,momv0);

						float Deltaphi        = mMyPidMaker -> GetDeltaPhi(phiv0, PsiEpdFlat[nMipMaxCut-2][0]);

						RapidityvsDeltaphi[1][MassBin] -> Fill(rapidityv0,Deltaphi);

						//Lambda rest frame
						TVector3 v0beta = mom4d_v0.BoostVector();
						TLorentzVector mom4d_p( TVector3( mom_neg.x(), mom_neg.y(), mom_neg.z() ), ep );
						TLorentzVector mom4d_pi( TVector3( mom_pos.x(), mom_pos.y(), mom_pos.z() ), ep );
						mom4d_p.Boost( -v0beta);
						mom4d_pi.Boost( -v0beta);
						float phip_V0RF = mom4d_p.Phi();
						float Rapidityp_V0RF = mom4d_p.Rapidity();
						float Ptp_V0RF = mom4d_p.Pt();
						float Momp_V0RF = mom4d_p.Rho();
						float phipi_V0RF = mom4d_pi.Phi();
						float Rapiditypi_V0RF = mom4d_pi.Rapidity();
						float Ptpi_V0RF = mom4d_pi.Pt();
						float Mompi_V0RF = mom4d_pi.Rho();
						float dphiEpd[4]={0};
						float dphiTpc[3]={0};
						if(1.11 < invM && invM < 1.121){
								hDaughterPionRapidityvsPtRF[1]   -> Fill(Rapiditypi_V0RF,Ptpi_V0RF);
								hDaughterProtonRapidityvsPtRF[1] -> Fill(Rapidityp_V0RF,Ptp_V0RF);
								hDaughterPionPhiRF[1]    -> Fill(phipi_V0RF);
								hDaughterProtonPhiRF[1]  -> Fill(phip_V0RF);
						}
						int RapidityBinProtonRF = mMyPidMaker -> GetYforLambda(Rapidityp_V0RF);
						int PtBinProtonRF       = mMyPidMaker -> GetPtforLambda(Ptp_V0RF); 
						int RapidityBinPionRF = mMyPidMaker -> GetYforLambda(Rapiditypi_V0RF);
						int PtBinPionRF       = mMyPidMaker -> GetPtforLambda(Ptpi_V0RF); 
						if(PtBinProtonRF != -1)pInvMvsDaughterProtonRapidityRF[1][PtBinProtonRF] -> Fill(invM,Rapidityp_V0RF);
						if(RapidityBinProtonRF != -1)pInvMvsDaughterProtonPtRF[1][RapidityBinProtonRF] -> Fill(invM,Ptp_V0RF);
						pInvMvsDaughterPionPhiRF[1]        -> Fill(invM,phip_V0RF);
						pInvMvsDaughterPionMomRF[1]        -> Fill(invM,Momp_V0RF);
						if(PtBinPionRF != -1)pInvMvsDaughterPionRapidityRF[1][PtBinPionRF] -> Fill(invM,Rapiditypi_V0RF);
						if(RapidityBinPionRF != -1)pInvMvsDaughterPionPtRF[1][RapidityBinPionRF] -> Fill(invM,Ptpi_V0RF);
						pInvMvsDaughterProtonPhiRF[1]      -> Fill(invM,phipi_V0RF);
						pInvMvsDaughterProtonMomRF[1]      -> Fill(invM,Mompi_V0RF);

						hInvMvsDeltaPhiLambdabar -> Fill(invM,atan2(sin(phiv0-phip_V0RF),cos(phiv0-phip_V0RF)));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999) pInvMvsSinPolLambdabar     -> Fill(invM,atan2(sin(phiv0-phip_V0RF),cos(phiv0-phip_V0RF)),-sin(dphiEpd[0]));

						//for(int etaId=0;etaId<3;etaId++){
						//		bool AutoCorrFlag = mMyPidMaker -> RemoveAutoCorr(etaId, RapidityBin, PtBin, etav0);
						//		dphiTpc[etaId] = PsiTpcFlat[etaId] - phip_V0RF;
						//		if(QvTpcRaw[etaId][0] != -9999 && QvTpcRaw[etaId][1] != -9999 && AutoCorrFlag)hinvMLambda[1][centnumber][RapidityBin][PtBin][etaId] -> Fill(invM);
						//		if(QvTpcRaw[etaId][0] != -9999 && QvTpcRaw[etaId][1] != -9999 && AutoCorrFlag)pInvMassvsSinpol[1][centnumber][RapidityBin][PtBin][etaId] -> Fill(invM,sin(dphiTpc[etaId]));
						//}

						//for(int iring=0;iring<4;iring++){
						//		dphiEpd[iring] = PsiEpdFlat[nMipMaxCut-2][iring] - phip_V0RF;
						//		if(QvEpdRaw[nMipMaxCut-2][iring][0] != -9999 && QvEpdRaw[nMipMaxCut-2][iring][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambda[1][centnumber][RapidityBin][PtBin][iring+3] -> Fill(invM);
						//		if(QvEpdRaw[nMipMaxCut-2][iring][0] != -9999 && QvEpdRaw[nMipMaxCut-2][iring][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpol[1][centnumber][RapidityBin][PtBin][iring+3] -> Fill(invM,sin(dphiEpd[iring]));
						//}
						dphiEpd[0] = PsiEpdFlat[nMipMaxCut-2][0] - phip_V0RF;
						int DeltaPhiBinforPol = mMyPidMaker -> GetDeltaPhiBinForPol(phiv0 - phip_V0RF);
						//cout << "Qx=" << QvEpdRaw[nMipMaxCut-2][0][0] << "   Qy=" << QvEpdRaw[nMipMaxCut-2][0][1] << "   rapidity=" << rapidityv0 << " bin = "<<RapidityBin  << "   pt=" << ptv0 << "   bin = " << PtBin << "   phi=" << phiv0 << "   phip_V0RF=" << phip_V0RF << "   DeltaPhiBin="<< phiv0-phip_V0RF << "   " << DeltaPhiBinforPol << endl;
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && DeltaPhiBinforPol != -1)hinvMLambda[1][centnumber][RapidityBin][PtBin][DeltaPhiBinforPol] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1 && DeltaPhiBinforPol != -1)pInvMassvsSinpol[1][centnumber][RapidityBin][PtBin][DeltaPhiBinforPol] -> Fill(invM,sin(dphiEpd[0]));

						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambdaBarCent[centnumber] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpolLambdaBarCent[centnumber] -> Fill(invM,sin(dphiEpd[0]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambdaBarRap[RapidityBin] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpolLambdaBarRap[RapidityBin] -> Fill(invM,sin(dphiEpd[0]));
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)hinvMLambdaBarPt[PtBin] -> Fill(invM);
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1)pInvMassvsSinpolLambdaBarPt[PtBin] -> Fill(invM,sin(dphiEpd[0]));

						hInvMvsPhiProton[1]      -> Fill(invM, phiprotonpool[1][proton]);
						hInvMvsPhiPion[1]        -> Fill(invM, phipionpool[1][pion]);
						hInvMvsPhiLambda[1]      -> Fill(invM, phiv0);
						hInvMvsPhistar[1]        -> Fill(invM, phip_V0RF);
						hInvMvsdPsiPhiProton[1]  -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phiprotonpool[1][proton]),cos(PsiEpdFlat[nMipMaxCut-2][0]-phiprotonpool[1][proton])));
						hInvMvsdPsiPhiPion[1]    -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phipionpool[1][pion]),cos(PsiEpdFlat[nMipMaxCut-2][0]-phipionpool[1][pion])));
						hInvMvsdPsiPhiLambda[1]  -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phiv0),cos(PsiEpdFlat[nMipMaxCut-2][0]-phiv0)));
						hInvMvsdPsiPhistar[1]    -> Fill(invM, atan2(sin(PsiEpdFlat[nMipMaxCut-2][0]-phip_V0RF),cos(PsiEpdFlat[nMipMaxCut-2][0]-phip_V0RF)));
						hInvMvsDaughterdPhi[1]   -> Fill(invM, atan2(sin(phiprotonpool[1][proton]-phipionpool[1][pion]),cos(phiprotonpool[1][proton]-phipionpool[1][pion])));
						hInvMvsdEta[1]           -> Fill(invM, atan2(sin(etaprotonpool[1][proton]-etapionpool[1][pion]),cos(etaprotonpool[1][proton]-etapionpool[1][pion])));
						hDaughterdPhiLaboRest[1] -> Fill(atan2(sin(phiprotonpool[1][proton]-phipionpool[1][pion]),cos(phiprotonpool[1][proton]-phipionpool[1][pion])),atan2(sin(phip_V0RF-phiv0),cos(phip_V0RF-phiv0)));

						int EpMethodBin = mMyPidMaker -> GetdPsiphiBinforEPmethod(PsiEpdFlat[nMipMaxCut-2][0], phiv0);
						int DeltaPhiEPmethod = mMyPidMaker -> DeltaPhiforEPmethod(phiv0, phip_V0RF);

						//Event Plane method
						if(RapidityBin!=-1&&PtBin!=-1&&EpMethodBin!=-1&&DeltaPhiEPmethod!=-1){
								hInvMLambdaBarCentEPmethod[centnumber][EpMethodBin][DeltaPhiEPmethod] -> Fill(invM);
								if(centnumber==2||centnumber==3||centnumber==4||centnumber==5)hInvMLambdaBarRapEPmethod[RapidityBin][EpMethodBin][DeltaPhiEPmethod] -> Fill(invM);
								if(centnumber==2||centnumber==3||centnumber==4||centnumber==5)hInvMLambdaBarPtEPmethod[PtBin][EpMethodBin][DeltaPhiEPmethod]        -> Fill(invM);
						}

						//Change Topological cut for systematic
						if(QvEpdRaw[nMipMaxCut-2][0][0] != -9999 && QvEpdRaw[nMipMaxCut-2][0][1] != -9999 && RapidityBin != -1 && PtBin != -1){
								bool isGoodLambdaProtonDCAminus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton]+0.1, dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl);//DCAproton >0.6 -> >0.5
								if(isGoodLambdaProtonDCAminus){
										hInvMCentTopoloSys[1][0][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][0][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][0][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][0][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][0][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][0][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaProtonDCAplus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton]-0.1, dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl);//DCAproton >0.6 -> >0.7
								if(isGoodLambdaProtonDCAplus){
										hInvMCentTopoloSys[1][1][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][1][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][1][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][1][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][1][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][1][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaPionDCAminus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion]+0.1, dcaDaughters, dca2vtx, decayl);//DCApion >1.7 -> >1.6
								if(isGoodLambdaPionDCAminus){
										hInvMCentTopoloSys[1][2][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][2][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][2][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][2][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][2][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][2][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaPionDCAplus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion]-0.1, dcaDaughters, dca2vtx, decayl);//DCApion >1.7 -> >1.8
								if(isGoodLambdaPionDCAplus){
										hInvMCentTopoloSys[1][3][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][3][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][3][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][3][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][3][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][3][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdappiDCAminus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters+0.1, dca2vtx, decayl);//Daughter DCA <0.8 -> <0.7
								if(isGoodLambdappiDCAminus){
										hInvMCentTopoloSys[1][4][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][4][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][4][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][4][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][4][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][4][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdappiDCAplus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters-0.1, dca2vtx, decayl);//Daughter DCA <0.8 -> <0.9
								if(isGoodLambdappiDCAplus){
										hInvMCentTopoloSys[1][5][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][5][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][5][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][5][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][5][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][5][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDCAminus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx+0.1, decayl);//Lambda DCA <0.6 -> <0.5
								if(isGoodLambdaDCAminus){
										hInvMCentTopoloSys[1][6][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][6][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][6][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][6][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][6][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][6][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDCAplus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx-0.1, decayl);//Lambda DCA <0.6 -> <0.7
								if(isGoodLambdaDCAplus){
										hInvMCentTopoloSys[1][7][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][7][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][7][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][7][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][7][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][7][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDecayLminus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl-0.5);//Decay L > 7.0 -> >6.5
								if(isGoodLambdaDecayLminus){
										hInvMCentTopoloSys[1][8][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][8][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][8][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][8][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][8][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][8][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
								bool isGoodLambdaDecayLplus = mMyPidMaker -> ChangeTopological(1,centnumber, dcaprotonpool[1][proton], dcapionpool[1][pion], dcaDaughters, dca2vtx, decayl+0.5);//Decay L > 6.0 -> >7.5
								if(isGoodLambdaDecayLplus){
										hInvMCentTopoloSys[1][9][centnumber][DeltaPhiEPmethod] -> Fill(invM);
										pInvMvsSinpolCentTopoloSys[1][9][centnumber][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										if(centnumber==2 || centnumber==3 || centnumber==4 || centnumber==5){
												hInvMRapTopoloSys[1][9][RapidityBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolRapTopoloSys[1][9][RapidityBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
												hInvMPtTopoloSys[1][9][PtBin][DeltaPhiEPmethod] -> Fill(invM);
												pInvMvsSinpolPtTopoloSys[1][9][PtBin][DeltaPhiEPmethod] -> Fill(invM,sin(dphiEpd[0]));
										}
								}
						}// skip wrong event plane

				}

		}
		//end of antiLambda loop
		return kStOK;
}

//__________________________________________________________________________________
bool PicoAnalyzer::isGoodTrack(const StPicoTrack *ptrk) {
		const Float_t pt  = ptrk->gMom().Perp(); // zero for global tracks
		const Float_t mom = ptrk->gMom().Mag();
		const Float_t eta = ptrk->gMom().PseudoRapidity();
		const Int_t nHits = ptrk->nHits(); //TPCHits? 
		const Float_t dca = ptrk->gDCA( picoEvent->primaryVertex() ).Mag();
		const Int_t nHitsFit = ptrk->nHitsFit();
		const Int_t nHitsPoss = ptrk->nHitsMax();
		const Float_t nHitsDedx = ptrk->nHitsDedx();
		const Float_t quality = (Float_t)nHitsFit/(Float_t)nHitsPoss;

		if( pt < 0.15 )  return false;
		if( mom < 0 )   return false;
		if( eta < -1.5 || 0 < eta) return false;
		//if( fabs(dca)>3.0 ) return false;
		if( nHitsFit < 10 )  return false;
		if( quality < 0.52 )  return false;
		if(nHitsDedx < 0 ) return false;

		return true;
}
