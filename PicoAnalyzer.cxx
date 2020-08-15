/**
 * \brief Example of how to read a file (list of files) using StPicoEvent classes
 *
 * RunPicoDstAnalyzer.C is an example of reading STAR picoDst format.
 * One can use either picoDst file or a list of picoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov
 * \date May 29, 2018
 *
 * BBC Event Plane Builder
 * \author Ding Chen
 *
 * Updated to Calculate v2
 * \date Dec 23, 2019
 *
 * Tweak to analyze BES-II FXT 3GeV, 7.2GeV and more
 * \author Ding Chen
 * \date Feb 19, 2020
 */

// This is needed for calling standalone classes (not needed on RACF)
//#define _VANILLA_ROOT_

// C++ headers
#include <iostream>
#include <map>
#include <iterator>
#include <stdio.h>

// ROOT headers
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TMath.h"
#include "TString.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

// PicoDst headers
#include "StRoot/StPicoEvent/StPicoDstReader.h"
#include "StRoot/StPicoEvent/StPicoHelix.h"
#include "StRoot/StPicoEvent/StPicoDst.h"
#include "StRoot/StPicoEvent/StPicoEvent.h"
#include "StRoot/StPicoEvent/StPicoTrack.h"
#include "StRoot/StPicoEvent/StPicoBTofHit.h"
#include "StRoot/StPicoEvent/StPicoBTowHit.h"
#include "StRoot/StPicoEvent/StPicoEmcTrigger.h"
#include "StRoot/StPicoEvent/StPicoBTofPidTraits.h"
#include "StRoot/StPicoEvent/StPicoTrackCovMatrix.h"
//EPD header files
#include "StRoot/StPicoEvent/StPicoEpdHit.h"
#include "StRoot/StEpdUtil/StEpdGeom.h"
#include "StRoot/StEpdUtil/StEpdEpFinder.h"
#include "StRoot/StEpdUtil/StEpdEpInfo.h"

// Define global constants
// const Int_t daynumber     = 6;
const Int_t _Ncentralities = 9; // 9 centrality bins
const Int_t _EpTermsMaxIni = 20; // Shift Order
const Int_t _nEventTypeBins = 5; // 5 etaRange
const Int_t _nEventTypeBins_tpc = 2; // 2 etaRange for TPC
const Double_t _massPion     = 0.13957039;
const Double_t _massKaon     = 0.493677;
const Double_t _massProton   = 0.938272081;
const Double_t _massPhi = 1.019461;
const Double_t _y_mid = -2.03; // mid rapidity

// const Int_t order         = 20;
// const Int_t twoorder      = 2 * order;
Double_t GetPsi(Double_t Qx, Double_t Qy, Int_t order);

//////////////////////////////// Main Function /////////////////////////////////
void PicoAnalyzer(const Char_t *inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_EpdEP",
                      Int_t   inputp1 = 1, // event plane orders: 1st, 2nd order \psi
                      Int_t   inputp2 = 0, // sysErr cut Indexes 0-15
                      Int_t   inputp3 = 0, // sysErr cut variations, each systematic check has 2 or 3 vertions
                      Int_t   inputp4 = 1 // Iteration of the analysis is. In this analysis, 2 iterations is enough
                    )
{

  Int_t EpOrder = inputp1; // Event plane Fourier expansion order = 1, 2, 3
  Int_t sys_cutN = inputp2; // sysErr cut Indexes 0-15
  Int_t sys_varN = inputp3; // sysErr cut variations, each systematic check has 2 or 3 vertions
  Int_t sys_iterN = inputp4; // Iteration of the analysis is. In this analysis, 2 iterations is enough
  TString sys_object[16]  = {"primary", "etaGap", "etaRange", "binning",
                                "vtxDiff", "mthdDiff", "vz", "vr",
                                "dedx", "dca", "nHitsFit", "ratio",
                                "nSigK", "mass2", "pT", "dipAngle"};
  if(sys_cutN == 0) std::cout<< "sys_cutN == 0: " << sys_object[sys_cutN] << std::endl; //variationID: {0}
  if(sys_cutN == 1) std::cout<< "sys_cutN == 1: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 2) std::cout<< "sys_cutN == 2: " << sys_object[sys_cutN] << std::endl; //variationID: {0, 1, 2, 3, 4}
  if(sys_cutN == 3) std::cout<< "sys_cutN == 3: " << sys_object[sys_cutN] << std::endl; //variationID: {0, 1, 2, 3, 4}
  if(sys_cutN == 4) std::cout<< "sys_cutN == 4: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 5) std::cout<< "sys_cutN == 5: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 6) std::cout<< "sys_cutN == 6: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 7) std::cout<< "sys_cutN == 7: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 8) std::cout<< "sys_cutN == 8: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 9) std::cout<< "sys_cutN == 9: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 10) std::cout<< "sys_cutN == 10: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 11) std::cout<< "sys_cutN == 11: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 12) std::cout<< "sys_cutN == 12: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 13) std::cout<< "sys_cutN == 13: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 14) std::cout<< "sys_cutN == 14: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  if(sys_cutN == 15) std::cout<< "sys_cutN == 15: " << sys_object[sys_cutN] << std::endl; //variationID: {0-10}
  outFile.Prepend(Form("sys_%s_var%d_iter%d", sys_object[sys_cutN], sys_varN, sys_iterN);

  int mEvtcut[5] = {0};
  int mTrkcut[6] = {0};
  // (0) ================== Read input files and set status =====================
  StPicoDstReader* picoReader = new StPicoDstReader(inFile);
  picoReader->Init();
  std::cout << "Explicit read status for some branches" << std::endl;
  picoReader->SetStatus("*",0);
  picoReader->SetStatus("Event",1);
  picoReader->SetStatus("Track",1);
  picoReader->SetStatus("BTofHit",1);
  picoReader->SetStatus("BTofPidTraits",1);
  picoReader->SetStatus("EpdHit",1);

  if( !picoReader->chain() ) {
      std::cout << "No chain has been found." << std::endl;
  }

  Long64_t eventsInTree = picoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  // ------------------ Get event numbers from TChain --------------------------------
  Long64_t events2read = picoReader->chain()->GetEntries();
  std::cout << "Number of events to read: " << events2read << std::endl;

  // (1) ================= Set up EPD EP info to get EPD event plane ============
  TString EpdEpOutputName = "EpdEpCorrectionHistograms_OUTPUT_";
  EpdEpOutputName += outFile;
  EpdEpOutputName += ".root";
  StEpdGeom *mEpdGeom = new StEpdGeom();
  Double_t mThresh = 0.3; // EPD EP by hand
  Double_t mMax = 2.0; // EPD EP by hand
  Double_t etaRange[_nEventTypeBins] = {-5.0,-4.4,-4.35,-3.95,-2.60}; // EPD eta range to set 4 sub EPD EP -5.0,-4.4,-4.35,-3.95,-2.60
  // sys_cutN == 1;
  if(sys_cutN == 1 && sys_varN == 1)  etaRange[2] = -4.0; // EPD-2 as reference
  if(sys_cutN == 1 && sys_varN == 2){ // EPD-3 as reference
    etaRange[2] = -4.3;
    etaRange[3] = -3.9;
  }
  TH2D wt("Order1etaWeight","Order1etaWeight",500,1.5,6.5,5,0,5);
  for (int ix=1; ix<501; ix++){
    for (int iy=1; iy<6; iy++){
      double eta = wt.GetXaxis()->GetBinCenter(ix);
      if(iy==1) wt.SetBinContent(ix,iy,1);
      else {
        if(eta<=abs(etaRange[iy-2]) && eta>abs(etaRange[iy-1])) wt.SetBinContent(ix,iy,1.0);
        else wt.SetBinContent(ix,iy,0.0);
      }
    }
  }
  // resolution //{0.305527,0.346768,0.407968,0.452254,0.47444,0.486652,0.437499,0.276291,0.263857}
  double d_resolution[2][_Ncentralities] = { // EPD-1
    {0.305527,0.346768,0.407968,0.452254,0.47444,0.486652,0.437499,0.276291,0.263857},
    {0.0553539,0.058153,0.265089,0.30708,0.165371,0.162038,0.0392603,0.0485935,0.0441441}
    //{0.158,   0.1737,  0.2045,  0.2256,  0.2218,  0.2430,  0.20613,0.1312,0.1276}
  };
  double d_resolution_EPD_3[_Ncentralities] = {0.189401,0.196268,0.195405,0.189716,0.177785,0.163757,0.170117,0.272917,0.296757};
  // v1 eta weighting
  // double pr0[9] =  {-0.00146843,-0.00120124,-0.0016722,-0.00170085,-0.00198228,-0.00281638,-0.00343895,-0.00415811,-0.00537868};
  // double pr1[9] =  {-0.000426648,-0.000325884,-0.000581794,-0.000747585,-0.00103814,-0.00151215,-0.00188201,-0.00212849,-0.00257166};
  // double pr2[9] =  {0.00484939,0.00620883,0.00834932,0.0117018,0.0154964,0.0198621,0.0245629,0.0297327,0.0355797};
  // double pr3[9] =  {-0.00138009,-0.00193065,-0.00283936,-0.0042699,-0.00576823,-0.00747499,-0.00921599,-0.0111905,-0.0133779};

  double lin[9] =        {-0.000479872,-0.000468419,-0.000698331,-0.00136243,-0.00227147,-0.00314487,-0.00381054,-0.00416527,-0.00382669};
  double cub[9] =        {0.000453689,0.000550043,0.00072002,0.00100187,0.00129868,0.00160751,0.0018985,0.00218509,0.00234319};
  TH2D *v1WtaWt = new TH2D("v1WtaWt","v1WtaWt",200,-6.5,-1.5,_Ncentralities,0.5,0.5+_Ncentralities);
  for (int ix=1; ix<201; ix++){
    for (int iy=1; iy<10; iy++){
      double eta = v1WtaWt->GetXaxis()->GetBinCenter(ix);
      v1WtaWt->SetBinContent(ix,iy,lin[iy-1]*(eta-_y_mid)+cub[iy-1]*pow((eta-_y_mid),3));
      // if(eta>=-3.66){
      //   v1WtaWt->SetBinContent(ix,iy,pr3[iy-1]*pow(eta-_y_mid,3) + pr2[iy-1]*pow(eta-_y_mid,1));
      // } else {
      //   v1WtaWt->SetBinContent(ix,iy,pr1[iy-1]*pow(eta-_y_mid,1) + pr0[iy-1]);
      // }

    }
  }

  TClonesArray * mEpdHits = new TClonesArray("StPicoEpdHit");
  unsigned int found;
  // --------------------- Retrieve EpdHits TClonesArray ----------------------------
  TChain *mPicoDst = picoReader->chain();
  mPicoDst->SetBranchStatus("EpdHit*",1,&found);   // note you need the asterisk
  std::cout << "EpdHit Branch returned found= " << found << std::endl; // ? What is the EpdHit branch ? Check it on StRoot.
  mPicoDst->SetBranchAddress("EpdHit",&mEpdHits);
  // (2) ================ Output files and histograms ==========================
  outFile.Append(".picoDst.result.root");
  TFile *outputFile = new TFile(outFile,"recreate");
  // ------------------- Event cuts QA histograms ------------------------------
  TH1D *hist_runId = new TH1D("hist_runId","Event runId",20001,-0.5,20000.5);
  TH1D *hist_eventCuts = new TH1D("hist_eventCuts","# of Events after cuts",10,-0.5,9.5);
  TH1D *hist_trackCuts = new TH1D("hist_trackCuts","# of tracks after cuts",10,-0.5,9.5);
  TH1D *hist_Vz_pri = new TH1D("hist_Vz_pri","V_{Z} [cm]",6000,-300.0,300.0);
  TH2D *hist_VyVx_pri = new TH2D("hist_VyVx_pri","V_{Y} [cm] vs. V_{X} [cm]",500,-5.0,5.0,500,-5.0,5.0);
  TH1D *hist_Vr_pri = new TH1D("hist_Vr_pri","V_{R} [cm]",500,0.0,20.0);
  TH1D *hist_triggerID = new TH1D("hist_triggerID","Event TriggerId",20001,-0.5,20000.5);
  TH1D *hist_Vz_cut = new TH1D("hist_Vz_cut","V_{Z} after cut [cm]",6000,-300.0,300.0);
  TH1D *hist_Vr_cut = new TH1D("hist_Vr_cut","V_{R} after cut [cm]",500,0.0,20.0);
  TH2D *hist_VyVx_cut = new TH2D("hist_VyVx_cut","V_{Y} [cm] vs. V_{X} after cut [cm]",500,-5.0,5.0,500,-5.0,5.0);
  // -------------------- Track loop QA histograms --------------------------------
  TH2D *hist_px_py=new TH2D("hist_px_py","hist_px_py",4000,-10.0,10.0,4000,-10.0,10.0);
  TH1D *hist_pz = new TH1D("hist_pz","p_{z} [GeV/c]",4000,-10.0,10.0);
  TH1D *hist_pt = new TH1D("hist_pt","p_{T} [GeV/c]",2000,0.0,10.0);
  TH1D *hist_mom = new TH1D("hist_mom","p_{mom} [GeV/c]",2000,0.0,10.0);
  TH1D *hist_mass2 = new TH1D("hist_mass2","hist_mass2",4000,-10.0,10.0);
  TH1D *hist_eta = new TH1D("hist_eta","#eta",200,-3.0,0.5);
  TH1D *hist_phi = new TH1D("hist_phi","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH1D *hist_ratio = new TH1D("hist_ratio","hist_ratio",100,0,2);
  TH1D *hist_nHits = new TH1D("hist_nHits","hist_nHits",100,-0.5,99.5);
  TH1D *hist_ndEdx = new TH1D("hist_ndEdx","hist_ndEdx",100,-0.5,99.5);
  TH1D *hist_DCA = new TH1D("hist_DCA","hist_DCA",100,0,10.0);
  TH2D *hist_px_py_cut=new TH2D("hist_px_py_cut","hist_px_py_cut",4000,-10.0,10.0,4000,-10.0,10.0);
  TH1D *hist_pz_cut = new TH1D("hist_pz_cut","p_{z} [GeV/c]",4000,-10.0,10.0);
  TH1D *hist_pt_cut = new TH1D("hist_pt_cut","p_{T} [GeV/c]",2000,0.0,10.0);
  TH1D *hist_mom_cut = new TH1D("hist_mom_cut","p_{mom} [GeV/c]",2000,0.0,10.0);
  TH1D *hist_mass2_cut = new TH1D("hist_mass2_cut","hist_mass2_cut",4000,-10.0,10.0);
  TH1D *hist_eta_cut = new TH1D("hist_eta_cut","#eta",200,-3.0,0.5);
  TH1D *hist_phi_cut = new TH1D("hist_phi_cut","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH1D *hist_ratio_cut = new TH1D("hist_ratio_cut","hist_ratio_cut",100,0,2);
  TH1D *hist_nHits_cut = new TH1D("hist_nHits_cut","hist_nHits_cut",100,-0.5,99.5);
  TH1D *hist_ndEdx_cut = new TH1D("hist_ndEdx_cut","hist_ndEdx_cut",100,-0.5,99.5);
  TH1D *hist_DCA_cut = new TH1D("hist_DCA_cut","hist_DCA_cut",100,0,10.0);
  // ------------------ Centrality QA histograms ----------------------------------
  TH1D *hist_cent = new TH1D("hist_cent","Centrality",_Ncentralities+1,-0.5,_Ncentralities+0.5);
  TH1D *hist_realTrackMult = new TH1D("hist_realTrackMult","Actual track multiplicity",1001,-0.5,1000.5);
  TH2D *hist_realTrackMult_refmult = new TH2D("hist_realTrackMult_refmult","Actual track multiplicity vs. RefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  TH2D *hist_realTrackMult_grefmult = new TH2D("hist_realTrackMult_grefmult","Actual track multiplicity vs. gRefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  TH2D *hist_realTrackMult_tofmult = new TH2D("hist_realTrackMult_tofmult","Actual track multiplicity vs. TofMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  // ------------------ EPD event plane histograms ----------------------------------
  TH2D *hist2_Epd_east_Qy_Qx_raw_ini[_nEventTypeBins];
  TH1D *hist_Epd_Sub_psi_raw_ini = new TH1D("hist_Epd_Sub_psi_raw_ini","raw EPD-Sub EP for each & every EPD hit in EPD-1",1024,-1.0,7.0);
  TH1D *hist_Epd_Sub_psi_Shifted_ini = new TH1D("hist_Epd_Sub_psi_Shifted_ini","shifted EPD-Sub EP for each & every EPD hit in EPD-1",1024,-1.0,7.0);
  TH1D *hist_Epd_east_psi_raw_ini[_nEventTypeBins],/**hist_Epd_east_psi_Weighted_ini[_nEventTypeBins],*/*hist_Epd_east_psi_Shifted_ini[_nEventTypeBins];
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    hist2_Epd_east_Qy_Qx_raw_ini[EventTypeId]= new TH2D(Form("hist2_Epd_east_Qy_Qx_raw_ini_%d",EventTypeId),Form("EPD east Qy vs Qx EventTypeId%d",EventTypeId),600,-3.0,3.0,600,-3.0,3.0);
    hist_Epd_east_psi_raw_ini[EventTypeId] = new TH1D(Form("hist_Epd_east_psi_raw_ini_%d",EventTypeId),Form("EPD east EP EventTypeId%d",EventTypeId),1024,-1.0,7.0);
    // hist_Epd_east_psi_Weighted_ini[EventTypeId] = new TH1D(Form("hist_Epd_east_psi_Weighted_ini_%d",EventTypeId),Form("EPD east EP (Weighted) EventTypeId%d",EventTypeId),1024,-1.0,7.0);
    hist_Epd_east_psi_Shifted_ini[EventTypeId] = new TH1D(Form("hist_Epd_east_psi_Shifted_ini_%d",EventTypeId),Form("EPD east EP (Weighted & Shifted) EventTypeId%d",EventTypeId),1024,-1.0,7.0);
  }
  // ------------------ EPD event plane ab intio QA histograms ----------------------------------
  TH1D *hist_Epdeta = new TH1D("hist_Epdeta","epd eta",700,-6.5,0.5);
  TH1D *hist_Epdphi = new TH1D("hist_Epdphi","epd phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());

  TProfile2D *profile2D_PpVsEta = new TProfile2D("profile2D_PpVsEta","<TnMIP> vs. #eta vs. supersector",700,-6.5,0.5,12,0.5,12.5,0.3,3.0,"");
  profile2D_PpVsEta->Sumw2();
  TH2D *h2_hits_PpVsEta = new TH2D("h2_hits_PpVsEta","# of hits vs. #eta vs. supersector ",700,-6.5,0.5,12,0.5,12.5);
  TH1D *hist_nMip = new TH1D("hist_nMip","nMIP of tile: 0:1:1 ",64,-0.5,9.5);
  TH2D *h2_nMip_eta_cent = new TH2D("h2_nMip_eta_cent","Sum of nMIP VS. #eta VS. centrality ",20,-6.5,-1.5,_Ncentralities,0.5,_Ncentralities+0.5);
  TH2D *h2_TtVsPp[_nEventTypeBins], *h2_TtVsPpNmip[_nEventTypeBins], *h2_TtVsPpHit[_nEventTypeBins];
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    h2_TtVsPp[EventTypeId]= new TH2D(Form("h2_TtVsPp_%d",EventTypeId),Form("Tile vs Supersector of #eta range %d",EventTypeId),12,0.5,12.5,31,0.5,31.5);
    h2_TtVsPpNmip[EventTypeId] = new TH2D(Form("h2_TtVsPpNmip_%d",EventTypeId),Form("nMIP in Tile vs Supersector of #eta range %d",EventTypeId),12,0.5,12.5,31,0.5,31.5);
    h2_TtVsPpHit[EventTypeId] = new TH2D(Form("h2_TtVsPpHit_%d",EventTypeId),Form("Hits in Tile vs Supersector of #eta range %d",EventTypeId),12,0.5,12.5,31,0.5,31.5);
  }
  // --------------------- TPC event plane QA histograms ----------------------------------
  TH2D *h2_dEdxVsPq = new TH2D("h2_dEdxVsPq","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *h2_dEdxVspTq = new TH2D("h2_dEdxVspTq","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *h2_beta = new TH2D("h2_beta","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *h2_mass = new TH2D("hist_mass","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  // --------------------- TPC EP PID QA histograms ----------------------------------
  TH1D *  hist_trackmult_proton = new TH1D("hist_trackmult_proton","hist_trackmult_proton",100,-0.5,99.5);
  TH1D *  hist_trackmult_kaonPlus = new TH1D("hist_trackmult_kaonPlus","hist_trackmult_kaonPlus",100,-0.5,99.5);
  TH1D *  hist_trackmult_kaonMinus = new TH1D("hist_trackmult_kaonMinus","hist_trackmult_kaonMinus",100,-0.5,99.5);
  TH1D *  hist_trackmult_pionPlus = new TH1D("hist_trackmult_pionPlus","hist_trackmult_pionPlus",100,-0.5,99.5);
  TH1D *  hist_trackmult_pionMinus = new TH1D("hist_trackmult_pionMinus","hist_trackmult_pionMinus",100,-0.5,99.5);
  TH1D *hist_pt_proton = new TH1D("hist_pt_proton","p_{T} [GeV/c]",1000,0.0,5.0);
  TH1D *hist_eta_proton = new TH1D("hist_eta_proton","#eta",200,-3.0,0.5);
  TH1D *hist_y_proton = new TH1D("hist_y_proton","Rapidity y",200,-3.0,0.5);
  TH1D *hist_phi_proton = new TH1D("hist_phi_proton","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH2D *hist_rap_eta_proton = new TH2D("hist_rap_eta_proton","proton y versus #eta",250,-2.5,0,250,-2.5,0);
  TH2D *hist_pt_y_proton = new TH2D("hist_pt_y_proton","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_proton = new TH2D("hist_pt_eta_proton","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_dEdx_proton = new TH2D("hist_dEdx_proton","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *hist_beta_proton = new TH2D("hist_beta_proton","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *hist_mass_proton = new TH2D("hist_mass_proton","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  TH1D *hist_pt_kaonPlus = new TH1D("hist_pt_kaonPlus","p_{T} [GeV/c]",1000,0.0,5.0);
  TH1D *hist_eta_kaonPlus = new TH1D("hist_eta_kaonPlus","#eta",200,-3.0,0.5);
  TH1D *hist_y_kaonPlus = new TH1D("hist_y_kaonPlus","y",200,-3.0,0.5);
  TH1D *hist_phi_kaonPlus = new TH1D("hist_phi_kaonPlus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH2D *hist_rap_eta_kaonPlus = new TH2D("hist_rap_eta_kaonPlus","kaonPlus y versus #eta",250,-2.5,0,250,-2.5,0);
  TH2D *hist_pt_y_kaonPlus = new TH2D("hist_pt_y_kaonPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_kaonPlus = new TH2D("hist_pt_eta_kaonPlus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_dEdx_kaonPlus = new TH2D("hist_dEdx_kaonPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *hist_beta_kaonPlus = new TH2D("hist_beta_kaonPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *hist_mass_kaonPlus = new TH2D("hist_mass_kaonPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  TH1D *hist_pt_kaonMinus = new TH1D("hist_pt_kaonMinus","p_{T} [GeV/c]",1000,0.0,5.0);
  TH1D *hist_eta_kaonMinus = new TH1D("hist_eta_kaonMinus","#eta",200,-3.0,0.5);
  TH1D *hist_y_kaonMinus = new TH1D("hist_y_kaonMinus","y",200,-3.0,0.5);
  TH1D *hist_phi_kaonMinus = new TH1D("hist_phi_kaonMinus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH2D *hist_rap_eta_kaonMinus = new TH2D("hist_rap_eta_kaonMinus","kaonMinus y versus #eta",250,-2.5,0,250,-2.5,0);
  TH2D *hist_pt_y_kaonMinus = new TH2D("hist_pt_y_kaonMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_kaonMinus = new TH2D("hist_pt_eta_kaonMinus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_dEdx_kaonMinus = new TH2D("hist_dEdx_kaonMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *hist_beta_kaonMinus = new TH2D("hist_beta_kaonMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *hist_mass_kaonMinus = new TH2D("hist_mass_kaonMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  TH1D *hist_pt_pionPlus = new TH1D("hist_pt_pionPlus","p_{T} [GeV/c]",1000,0.0,5.0);
  TH1D *hist_eta_pionPlus = new TH1D("hist_eta_pionPlus","#eta",200,-3.0,0.5);
  TH1D *hist_y_pionPlus = new TH1D("hist_y_pionPlus","y",200,-3.0,0.5);
  TH1D *hist_phi_pionPlus = new TH1D("hist_phi_pionPlus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH2D *hist_rap_eta_pionPlus = new TH2D("hist_rap_eta_pionPlus","pionPlus y versus #eta",250,-2.5,0,250,-2.5,0);
  TH2D *hist_pt_y_pionPlus = new TH2D("hist_pt_y_pionPlus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_pionPlus = new TH2D("hist_pt_eta_pionPlus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_dEdx_pionPlus = new TH2D("hist_dEdx_pionPlus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *hist_beta_pionPlus = new TH2D("hist_beta_pionPlus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *hist_mass_pionPlus = new TH2D("hist_mass_pionPlus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  TH1D *hist_pt_pionMinus = new TH1D("hist_pt_pionMinus","p_{T} [GeV/c]",1000,0.0,5.0);
  TH1D *hist_eta_pionMinus = new TH1D("hist_eta_pionMinus","#eta",200,-3.0,0.5);
  TH1D *hist_y_pionMinus = new TH1D("hist_y_pionMinus","y",200,-3.0,0.5);
  TH1D *hist_phi_pionMinus = new TH1D("hist_phi_pionMinus","#phi [Radian]",1000,-0.5*TMath::Pi(),2.5*TMath::Pi());
  TH2D *hist_rap_eta_pionMinus = new TH2D("hist_rap_eta_pionMinus","pionMinus y versus #eta",250,-2.5,0,250,-2.5,0);
  TH2D *hist_pt_y_pionMinus = new TH2D("hist_pt_y_pionMinus","p_{T} [GeV/c] vs. y",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_pt_eta_pionMinus = new TH2D("hist_pt_eta_pionMinus","p_{T} [GeV/c] vs. #eta",500,-3.0,0.5,500,0.0,3.5);
  TH2D *hist_dEdx_pionMinus = new TH2D("hist_dEdx_pionMinus","dE/dx vs q*|p|",500,-3.0,3.0,500,0.0,10.0);
  TH2D *hist_beta_pionMinus = new TH2D("hist_beta_pionMinus","1/#beta vs q*|p|",1000,-5.0,5.0,500,0.0,5.0);
  TH2D *hist_mass_pionMinus = new TH2D("hist_mass_pionMinus","m^{2} vs q*|p|",1000,-5.0,5.0,1000,-0.6,4.0);
  // -------------------------- TPC event planes ----------------------------------
  Double_t etaRange_tpc[2] = {-0.6,0.}; // TPC eta range {-0.4, 0.0}
  TH2D wt_tpc("Order1etaWeight_tpc","Order1etaWeight_tpc",300,0,3.0,2,0,2);
  for (int ix=1; ix<301; ix++){
    for (int iy=1; iy<3; iy++){
      double eta = wt_tpc.GetXaxis()->GetBinCenter(ix);
      if(iy==1) wt_tpc.SetBinContent(ix,iy,1);
      else {
        if(eta<=abs(etaRange_tpc[iy-2]) && eta>abs(etaRange_tpc[iy-1])) wt_tpc.SetBinContent(ix,iy,1.0);
        else wt_tpc.SetBinContent(ix,iy,0.0);
      }
    }
  }
  TProfile2D *profile2D_v1VsEtaTpcOnly = new TProfile2D("profile2D_v1VsEtaTpcOnly","<( y - y_{CM} ) * cos ( #phi_{Track} - #psi_{EPD-full} ) > vs #eta vs centrality"
  ,64,-3.0,3.0,_Ncentralities,0.5,0.5+_Ncentralities,"");
  profile2D_v1VsEtaTpcOnly->Sumw2();
  TProfile2D *profile2D_v1VsEtaTpcOnly_1 = new TProfile2D("profile2D_v1VsEtaTpcOnly_1","< cos ( #phi_{Track} - #psi_{EPD-full} ) > vs #eta vs centrality"
  ,64,-3.0,3.0,_Ncentralities,0.5,0.5+_Ncentralities,"");
  profile2D_v1VsEtaTpcOnly_1->Sumw2();
  TH2D *hist_nTracksVsEta= new TH2D("hist_nTracksVsEta","# of good tracks VS #eta",64,-3.0,3.0,_Ncentralities,0.5,0.5+_Ncentralities);
  TH1D *hist_tpc_all_psi_raw[_nEventTypeBins_tpc] ;
  TH1D *hist_tpc_all_psi_shifted[_nEventTypeBins_tpc] ;
  for(int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
    hist_tpc_all_psi_raw[EventTypeId_tpc]= new TH1D(Form("hist_tpc_all_psi_raw_%d",EventTypeId_tpc),Form("TPC-sub%d event plane",EventTypeId_tpc),500,-0.5*TMath::Pi(),2.5*TMath::Pi());
    hist_tpc_all_psi_shifted[EventTypeId_tpc] = new TH1D(Form("hist_tpc_all_psi_shifted_%d",EventTypeId_tpc),Form("TPC-sub%d EP (shifted)",EventTypeId_tpc),500,-0.5*TMath::Pi(),2.5*TMath::Pi());
  }
  // "Shift correction" histograms that we INPUT and apply here
  TProfile2D *mEpdShiftInput_sin[_nEventTypeBins], *mEpdShiftInput_cos[_nEventTypeBins];
  TProfile2D *mTpcShiftInput_sin[_nEventTypeBins_tpc], *mTpcShiftInput_cos[_nEventTypeBins_tpc]; // TPC EP input
  // TH1D* mPhiWeightInput[_nEventTypeBins];
  TFile* mCorrectionInputFile = new TFile("EpCorrection_INPUT.root","READ");
  if (mCorrectionInputFile->IsZombie()) {
    std::cout << "Error opening file with Ab initio Correction Histograms" << std::endl;
    std::cout << "I will use no correction at all for my own EPD Ep." << std::endl;
    for (int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
      mEpdShiftInput_sin[EventTypeId] = 0;
    	mEpdShiftInput_cos[EventTypeId] = 0;
      // mPhiWeightInput[EventTypeId] = 0;
    }
    for (int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
      mTpcShiftInput_sin[EventTypeId_tpc] = 0;
    	mTpcShiftInput_cos[EventTypeId_tpc] = 0;
    }
  }
  else{
    for (int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
      mEpdShiftInput_sin[EventTypeId] = (TProfile2D*)mCorrectionInputFile->Get(Form("EpdShiftEW0Psi%d_sin",EventTypeId));
      mEpdShiftInput_cos[EventTypeId] = (TProfile2D*)mCorrectionInputFile->Get(Form("EpdShiftEW0Psi%d_cos",EventTypeId));
      // mPhiWeightInput[EventTypeId] = (TH1D*)mCorrectionInputFile->Get(Form("PhiWeight%d",EventTypeId));
      // mPhiWeightInput[EventTypeId]->Scale((double)12.0/((double)(mPhiWeightInput[EventTypeId]->GetEntries())));
    }
    for (int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
      mTpcShiftInput_sin[EventTypeId_tpc] = (TProfile2D*)mCorrectionInputFile->Get(Form("mTpcShiftOutput_%d_sin",EventTypeId_tpc));
      mTpcShiftInput_cos[EventTypeId_tpc] = (TProfile2D*)mCorrectionInputFile->Get(Form("mTpcShiftOutput_%d_cos",EventTypeId_tpc));
    }
  }

  // "Shift correction" histograms that we produce and OUTPUT
  TString EpOutputNameIni = "EpCorrection_OUTPUT_";
  EpOutputNameIni += outFile;
  EpOutputNameIni += ".root";
  TFile* mCorrectionOutputFile = new TFile(EpOutputNameIni,"RECREATE");
  TProfile2D *mEpdShiftOutput_sin[_nEventTypeBins], *mEpdShiftOutput_cos[_nEventTypeBins]; // EPD EP output
  TProfile2D *mTpcShiftOutput_sin[_nEventTypeBins_tpc], *mTpcShiftOutput_cos[_nEventTypeBins_tpc]; // TPC EP output
  // TH1D* mPhiWeightOutput[_nEventTypeBins];     // the array index is for EPD sub 0,1,2,3,4
  // TH1D* mPhiAveraged[_nEventTypeBins];         // the bins are (Phi bin) Sum of TnMIP vs Phi bin
  TProfile2D *profile2D_v1VsCentVsEta = new TProfile2D("profile2D_v1VsCentVsEta","v_{1} vs. #eta vs. centrality",
          40,-7.0,3.0, // total eta range
          _Ncentralities,0.5,_Ncentralities+0.5, // Centrality
          -1.0,1.0,"");//Use EPD-3 as primary event plane
  profile2D_v1VsCentVsEta->Sumw2();
  TProfile *profile_v1VsEta[_Ncentralities]; // [] is from 0 to 8, centrality is from 1 to 9.
  for(int cent=0; cent<_Ncentralities; cent++){
    profile_v1VsEta[cent]   = new TProfile(Form("profile_v1VsEta_cent%d",cent),Form("Directed flow VS. #eta in cent bin %d",cent),40,-7.0,3.0,-1.0,1.0,"");
    profile_v1VsEta[cent]->Sumw2();
  }
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    // mPhiWeightOutput[EventTypeId]   = new TH1D(Form("PhiWeight%d",EventTypeId),Form("Phi Weight divided by Averaged EPD-%d",EventTypeId),12,0.,2.0*TMath::Pi()); // bins are Phi bin
    // mPhiAveraged[EventTypeId]       = new TH1D(Form("PhiAveraged%d",EventTypeId),Form("Average for this phi EPD-%d",EventTypeId),12,0.,2.0*TMath::Pi()); // just for normalization. discard after use
    mEpdShiftOutput_sin[EventTypeId] = new TProfile2D(Form("EpdShiftEW0Psi%d_sin",EventTypeId),Form("EpdShiftEW0Psi%d_sin",EventTypeId),
            _EpTermsMaxIni,0.5,1.0*_EpTermsMaxIni+.5, // Shift order
            _Ncentralities,0.5,_Ncentralities+0.5, // Centrality
            -1.0,1.0);
    mEpdShiftOutput_cos[EventTypeId] = new TProfile2D(Form("EpdShiftEW0Psi%d_cos",EventTypeId),Form("EpdShiftEW0Psi%d_cos",EventTypeId),
            _EpTermsMaxIni,0.5,1.0*_EpTermsMaxIni+.5, // Shift order
            _Ncentralities,0.5,_Ncentralities+0.5, // Centrality
            -1.0,1.0);
  }
  for(int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
    mTpcShiftOutput_sin[EventTypeId_tpc] = new TProfile2D(Form("mTpcShiftOutput_%d_sin",EventTypeId_tpc),Form("mTpcShiftOutput_%d_sin",EventTypeId_tpc),
            _EpTermsMaxIni,0.5,1.0*_EpTermsMaxIni+.5, // Shift order
            _Ncentralities,0.5,_Ncentralities+0.5, // Centrality
            -1.0,1.0);
    mTpcShiftOutput_cos[EventTypeId_tpc] = new TProfile2D(Form("mTpcShiftOutput_%d_cos",EventTypeId_tpc),Form("mTpcShiftOutput_%d_cos",EventTypeId_tpc),
            _EpTermsMaxIni,0.5,1.0*_EpTermsMaxIni+.5, // Shift order
            _Ncentralities,0.5,_Ncentralities+0.5, // Centrality
            -1.0,1.0);
  }
  // ------------------ TPC event plane ab intio Correlations histograms ----------------------------------
  TProfile *profile_correlation_epd_east[2][6], *profile_correlation_epd_tpc[2][4], *profile_correlation_epd_tpc_all[2];
  TH2D *correlation2D_epd_east[6],*correlation2D_epd_tpc[4], *correlation2D_epd_tpc_all;
  int pairs =0;
  for(int i = 0; i<3;i++){ // Correlations between EPD EP 1, 2, 3, 4. 6 pairs of correlations
    for(int j=i+1;j<4;j++){
      for(int n=0; n<2; n++){
        profile_correlation_epd_east[n][pairs]  =
        new TProfile(Form("profile_correlation_n%d_epd_east%d",n+1,pairs),
        Form("<cos(%d * (#psi^{EPD east}[%d] #minus #psi^{EPD east}[%d]))>",n+1,i+1,j+1),
        _Ncentralities,0.5,_Ncentralities+0.5,-1.0,1.0,"");
      }
      correlation2D_epd_east[pairs]   =
      new TH2D(Form("correlation2D_epd_east%d",pairs),
      Form("#psi^{EPD east}[%d] vs. #psi^{EPD east}[%d]",i+1,j+1),
      50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
      pairs++;
    }
  }
  for(int i=0;i<4;i++){// Correlaitons between TPC and 4 EPD event planes 1,2,3,4
    for(int n=0; n<2; n++){
      profile_correlation_epd_tpc[n][i]  =
      new TProfile(Form("profile_correlation_n%d_epd%d_tpc",n+1,i+1),
      Form("<cos(%d * (#psi^{EPD east}[%d] #minus #psi^{TPC}))>",n+1 ,i+1),
      _Ncentralities,0.5,_Ncentralities+0.5,-1.0,1.0,"");
    }
    correlation2D_epd_tpc[i]   =
    new TH2D(Form("correlation2D_epd%d_tpc",i+1),
    Form("#psi^{EPD east}[%d] vs. #psi^{TPC}",i+1),
    50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  }
  for(int n=0; n<2; n++){
    profile_correlation_epd_tpc_all[n]  =
    new TProfile(Form("profile_correlation_n%d_epd_tpc_all",n+1),
    Form("<cos(%d * (#psi^{EPD east}[full] #minus #psi^{TPC}))>", n+1),
    _Ncentralities,0.5,_Ncentralities+0.5,-1.0,1.0,"");
  }
  correlation2D_epd_tpc_all   =
  new TH2D("correlation2D_epd_tpc_all",
  "#psi^{EPD east}[full] vs. #psi^{TPC}",
  50,-0.5*TMath::Pi(),2.5*TMath::Pi(),50,-0.5*TMath::Pi(),2.5*TMath::Pi());
  // ------------- phi-meson output file and plots -----------------------------
  double ptSetA[3]  = {0.6, 1.2, 2.4};
  double ptSetB[5]  = {0.4, 0.7, 1.0, 1.4, 2.0};
  double ptSetC[11] = {0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0, 4.0};

  double rapSetA[5]  = {-2.0, -1.5, -1.0, -0.5, 0};

  double centSetA[5]  = {0, 10, 40, 60, 80}; // %
  double centSetB[10]  = {0, 5, 10, 20, 30, 40, 50, 60, 70, 80}; // %
  TString PhiOutputName = "PhiMesonAna_OUTPUT_";
  PhiOutputName += outFile;
  PhiOutputName += ".root";
  TFile* PhiMesonAnaOutputFile = new TFile(PhiOutputName,"RECREATE");
  TH1D * hist_dip_angle = new TH1D("hist_dip_angle","hist_dip_angle",1000,-1,1.0);
  TH1D * hist_mother_decay_length = new TH1D("hist_mother_decay_length","hist_mother_decay_length",1000,-1.0,4.0);
  TH1D * hist_SE_mass_Phi     = new TH1D("hist_SE_mass_Phi","Same event invariant mass",200,0.9,1.1);
  TH1D *hist_SE_PhiMeson_pT  = new TH1D("hist_SE_PhiMeson_pT","pT distribution of #phi",200,0.0,10);
  TH1D *hist_SE_PhiMeson_mT  = new TH1D("hist_SE_PhiMeson_mT","mT distribution of #phi",200,0.0,10);
  TH1D *hist_SE_PhiMeson_rap  = new TH1D("hist_SE_PhiMeson_rap","y distribution of #phi",200,-10.,10);
  TH2D *hist_SE_pt_y_PhiMeson = new TH2D("hist_SE_pt_y_PhiMeson","p_{T} [GeV/c] vs. y of #phi",500,-3.0,0.5,500,0.0,3.5);
  TH2D * h2_TOF_beta_pq       = new TH2D("h2_TOF_beta_pq","1/#beta vs. pq",500,-3,3,500,0,3);
// pt SetA, cent SetA
  TH1D *mHist_SE_InvM_ptSetA_centSetA[2][6];
  TH2D *mHist_v1_raw_ptSetA_centSetA[2][6];
  TH2D *mHist_v1_reso_ptSetA_centSetA[2][6];
  TH2D *mHist_v2_raw_ptSetA_centSetA[2][6];
  TH2D *mHist_v2_reso_ptSetA_centSetA[2][6];
  TProfile *mProfile_v1_raw_ptSetA_centSetA[2][6];
  TProfile *mProfile_v1_reso_ptSetA_centSetA[2][6];
  TProfile *mProfile_v2_raw_ptSetA_centSetA[2][6];
  TProfile *mProfile_v2_reso_ptSetA_centSetA[2][6];
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_ptSetA_centSetA[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetA_centSetA[pt][cent] = new TH2D(Form("Hist_v1_raw_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_v1_raw_ptSetA%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetA_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");
      mHist_v1_reso_ptSetA_centSetA[pt][cent] = new TH2D(Form("Hist_v1_reso_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_v1_reso_ptSetA%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_ptSetA_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>/R_{1}^{EPD}");
      mHist_v2_raw_ptSetA_centSetA[pt][cent] = new TH2D(Form("Hist_v2_raw_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_v2_raw_ptSetA%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_ptSetA_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");
      mHist_v2_reso_ptSetA_centSetA[pt][cent] = new TH2D(Form("Hist_v2_reso_ptSetA%d_centSetA%d",pt,cent),
      Form("Hist_v2_reso_ptSetA%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_ptSetA_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_ptSetA_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
  }
  // pt SetA, cent SetB
  TH1D *mHist_SE_InvM_ptSetA_centSetB[2][9];
  TH2D *mHist_v1_raw_ptSetA_centSetB[2][9];
  TH2D *mHist_v1_reso_ptSetA_centSetB[2][9];
  TH2D *mHist_v2_raw_ptSetA_centSetB[2][9];
  TH2D *mHist_v2_reso_ptSetA_centSetB[2][9];
  TProfile *mProfile_v1_raw_ptSetA_centSetB[2][9];
  TProfile *mProfile_v1_reso_ptSetA_centSetB[2][9];
  TProfile *mProfile_v2_raw_ptSetA_centSetB[2][9];
  TProfile *mProfile_v2_reso_ptSetA_centSetB[2][9];
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_ptSetA_centSetB[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetA%d_centSetB%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_ptSetA_centSetB[pt][cent] = new TH2D(Form("Hist_v1_raw_ptSetA%d_centSetB%d",pt,cent),
      Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetA_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_ptSetA_centSetB[pt][cent] = new TH2D(Form("Hist_v1_reso_ptSetA%d_centSetB%d",pt,cent),
      Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_ptSetA_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>/R_{1}^{EPD}");

      mHist_v2_raw_ptSetA_centSetB[pt][cent] = new TH2D(Form("Hist_v2_raw_ptSetA%d_centSetB%d",pt,cent),
      Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_ptSetA_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");

      mHist_v2_reso_ptSetA_centSetB[pt][cent] = new TH2D(Form("Hist_v2_reso_ptSetA%d_centSetB%d",pt,cent),
      Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_ptSetA_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_ptSetA_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
  }
  // pt SetB, cent SetA
  TH1D *mHist_SE_InvM_ptSetB_centSetA[4][6];
  TH2D *mHist_v1_raw_ptSetB_centSetA[4][6];
  TH2D *mHist_v1_reso_ptSetB_centSetA[4][6];
  TH2D *mHist_v2_raw_ptSetB_centSetA[4][6];
  TH2D *mHist_v2_reso_ptSetB_centSetA[4][6];
  TProfile *mProfile_v1_raw_ptSetB_centSetA[4][6];
  TProfile *mProfile_v1_reso_ptSetB_centSetA[4][6];
  TProfile *mProfile_v2_raw_ptSetB_centSetA[4][6];
  TProfile *mProfile_v2_reso_ptSetB_centSetA[4][6];
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_ptSetB_centSetA[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_SE_InvM_ptSetA%d_centSetA%d",pt,cent),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_ptSetB_centSetA[pt][cent] = new TH2D(Form("Hist_v1_raw_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_v1_raw_ptSetB%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetB_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_ptSetB_centSetA[pt][cent] = new TH2D(Form("Hist_v1_reso_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_v1_reso_ptSetB%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_ptSetB_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");

      mHist_v2_raw_ptSetB_centSetA[pt][cent] = new TH2D(Form("Hist_v2_raw_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_v2_raw_ptSetB%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_ptSetB_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v2_reso_ptSetB_centSetA[pt][cent] = new TH2D(Form("Hist_v2_reso_ptSetB%d_centSetA%d",pt,cent),
      Form("Hist_v2_reso_ptSetB%d_centSetA%d",pt,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_ptSetB_centSetA[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_ptSetB_centSetA[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
  }
  // pt SetB, cent SetB
  TH1D *mHist_SE_InvM_ptSetB_centSetB[4][9];
  TH2D *mHist_v1_raw_ptSetB_centSetB[4][9];
  TH2D *mHist_v1_reso_ptSetB_centSetB[4][9];
  TH2D *mHist_v2_raw_ptSetB_centSetB[4][9];
  TH2D *mHist_v2_reso_ptSetB_centSetB[4][9];
  TProfile *mProfile_v1_raw_ptSetB_centSetB[4][9];
  TProfile *mProfile_v1_reso_ptSetB_centSetB[4][9];
  TProfile *mProfile_v2_raw_ptSetB_centSetB[4][9];
  TProfile *mProfile_v2_reso_ptSetB_centSetB[4][9];
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_ptSetB_centSetB[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetB%d_centSetB%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_ptSetB_centSetB[pt][cent] = new TH2D(Form("Hist_v1_raw_ptSetB%d_centSetB%d",pt,cent),
      Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetB_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_ptSetB_centSetB[pt][cent] = new TH2D(Form("Hist_v1_reso_ptSetB%d_centSetB%d",pt,cent),
      Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_ptSetB_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>/R_{1}^{EPD}");

      mHist_v2_raw_ptSetB_centSetB[pt][cent] = new TH2D(Form("Hist_v2_raw_ptSetB%d_centSetB%d",pt,cent),
      Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_ptSetB_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");

      mHist_v2_reso_ptSetB_centSetB[pt][cent] = new TH2D(Form("Hist_v2_reso_ptSetB%d_centSetB%d",pt,cent),
      Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[pt],ptSetB[pt+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_ptSetB_centSetB[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_ptSetB_centSetB[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
  }
  // pt SetC, cent 0-60%, 0-80%
  TH1D *mHist_SE_InvM_ptSetC_centAll[10][2];
  TH2D *mHist_v1_raw_ptSetC_centAll[10][2];
  TH2D *mHist_v1_reso_ptSetC_centAll[10][2];
  TH2D *mHist_v2_raw_ptSetC_centAll[10][2];
  TH2D *mHist_v2_reso_ptSetC_centAll[10][2];
  TProfile *mProfile_v1_raw_ptSetC_centAll[10][2];
  TProfile *mProfile_v1_reso_ptSetC_centAll[10][2];
  TProfile *mProfile_v2_raw_ptSetC_centAll[10][2];
  TProfile *mProfile_v2_reso_ptSetC_centAll[10][2];
  for(int pt=0; pt<10; pt++)
  {
    for(int cent=0; cent<2;cent++){
      mHist_SE_InvM_ptSetC_centAll[pt][cent] = new TH1D(Form("Hist_SE_InvM_ptSetC%d_centAll%d",pt,cent),
      Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      200,0.9,1.1);
      mHist_SE_InvM_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_ptSetC_centAll[pt][cent] = new TH2D(Form("Hist_v1_raw_ptSetC%d_centAll%d",pt,cent),
      Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_ptSetC_centAll[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_ptSetC_centAll[pt][cent] = new TH2D(Form("Hist_v1_reso_ptSetC%d_centAll%d",pt,cent),
      Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_ptSetC_centAll[pt][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>/R_{1}^{EPD}");

      mHist_v2_raw_ptSetC_centAll[pt][cent] = new TH2D(Form("Hist_v2_raw_ptSetC%d_centAll%d",pt,cent),
      Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_ptSetC_centAll[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");

      mHist_v2_reso_ptSetC_centAll[pt][cent] = new TH2D(Form("Hist_v2_reso_ptSetC%d_centAll%d",pt,cent),
      Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetC[pt],ptSetC[pt+1],centSetA[0],centSetA[cent+3]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_ptSetC_centAll[pt][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_ptSetC_centAll[pt][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>/R_{1}^{EPD}");
    }
  }
  // rap SetA, cent SetA
  TH1D *mHist_SE_InvM_rapSetA_centSetA[4][6];
  TH2D *mHist_v1_raw_rapSetA_centSetA[4][6];
  TH2D *mHist_v1_reso_rapSetA_centSetA[4][6];
  TH2D *mHist_v2_raw_rapSetA_centSetA[4][6];
  TH2D *mHist_v2_reso_rapSetA_centSetA[4][6];
  TProfile *mProfile_v1_raw_rapSetA_centSetA[4][6];
  TProfile *mProfile_v1_reso_rapSetA_centSetA[4][6];
  TProfile *mProfile_v2_raw_rapSetA_centSetA[4][6];
  TProfile *mProfile_v2_reso_rapSetA_centSetA[4][6];
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<6;cent++){
      mHist_SE_InvM_rapSetA_centSetA[rap][cent] = new TH1D(Form("Hist_SE_InvM_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_SE_InvM_rapSetA%d_centSetA%d",rap,cent),
      200,0.9,1.1);
      mHist_SE_InvM_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_rapSetA_centSetA[rap][cent] = new TH2D(Form("Hist_v1_raw_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_v1_raw_rapSetA%d_centSetA%d",rap,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_rapSetA_centSetA[rap][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_rapSetA_centSetA[rap][cent] = new TH2D(Form("Hist_v1_reso_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_v1_reso_rapSetA%d_centSetA%d",rap,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_rapSetA_centSetA[rap][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v2_raw_rapSetA_centSetA[rap][cent] = new TH2D(Form("Hist_v2_raw_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_v2_raw_rapSetA%d_centSetA%d",rap,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_rapSetA_centSetA[rap][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");

      mHist_v2_reso_rapSetA_centSetA[rap][cent] = new TH2D(Form("Hist_v2_reso_rapSetA%d_centSetA%d",rap,cent),
      Form("Hist_v2_reso_rapSetA%d_centSetA%d",rap,cent),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_rapSetA_centSetA[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_rapSetA_centSetA[rap][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");
    }
  }
  // rap SetA, cent SetB
  TH1D *mHist_SE_InvM_rapSetA_centSetB[4][9];
  TH2D *mHist_v1_raw_rapSetA_centSetB[4][9];
  TH2D *mHist_v1_reso_rapSetA_centSetB[4][9];
  TH2D *mHist_v2_raw_rapSetA_centSetB[4][9];
  TH2D *mHist_v2_reso_rapSetA_centSetB[4][9];
  TProfile *mProfile_v1_raw_rapSetA_centSetB[4][9];
  TProfile *mProfile_v1_reso_rapSetA_centSetB[4][9];
  TProfile *mProfile_v2_raw_rapSetA_centSetB[4][9];
  TProfile *mProfile_v2_reso_rapSetA_centSetB[4][9];
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<9;cent++){
      mHist_SE_InvM_rapSetA_centSetB[rap][cent] = new TH1D(Form("Hist_SE_InvM_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      200,0.9,1.1);
      mHist_SE_InvM_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");

      mHist_v1_raw_rapSetA_centSetB[rap][cent] = new TH2D(Form("Hist_v1_raw_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_raw_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_raw_rapSetA_centSetB[rap][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v1_reso_rapSetA_centSetB[rap][cent] = new TH2D(Form("Hist_v1_reso_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v1_reso_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v1_reso_rapSetA_centSetB[rap][cent]->GetYaxis()->SetTitle("<cos(#phi - #psi_{1})>");

      mHist_v2_raw_rapSetA_centSetB[rap][cent] = new TH2D(Form("Hist_v2_raw_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_raw_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_raw_rapSetA_centSetB[rap][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");

      mHist_v2_reso_rapSetA_centSetB[rap][cent] = new TH2D(Form("Hist_v2_reso_rapSetA%d_centSetB%d",rap,cent),
      Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[rap],rapSetA[rap+1],centSetB[cent],centSetB[cent+1]),
      100,0.9,1.1,
      1000,-1.0,1.0);
      mHist_v2_reso_rapSetA_centSetB[rap][cent]->GetXaxis()->SetTitle("mass [GeV/c^{2}]");
      mHist_v2_reso_rapSetA_centSetB[rap][cent]->GetYaxis()->SetTitle("<cos(2(#phi - #psi_{1}))>");
    }
  }
  // ------------------ EPD & TPC event plane ab intio Correlations histograms ----------------------------------
  // (3) =========================== Event loop ====================================
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
  {
    // ---------------------- Event reading quality assurance ----------------------
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1)
    << "/" << events2read << "]" << std::endl;
    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
        std::cout << "Something went wrong, my Lord! Nothing to analyze..."
        << std::endl;
        break;
    }
    StPicoDst     *dst = picoReader->picoDst();
    StPicoEvent *event = dst->event();
    if( !event ) {
        std::cout << "Something went wrong, my Lord! Event is hiding from me..."
        << std::endl;
        break;
    }
    mEvtcut[0]++;// No event cut yet
    // (4) =================== Get event parameters ================================
    Int_t runId       = event->runId();
    Int_t nTracks     = dst->numberOfTracks();

    const Float_t   f_MagField = event->bField(); // Magnetic field
    Double_t Day      = (Double_t)runId - 19151028.0; // a day bin
    hist_runId->Fill(Day);

    Double_t primaryVertex_X    = (Double_t)event->primaryVertex().X();
    Double_t primaryVertex_Y    = (Double_t)event->primaryVertex().Y();
    Double_t primaryVertex_Z    = (Double_t)event->primaryVertex().Z();
    Double_t primaryVertex_perp = (Double_t)event->primaryVertex().Perp();
    hist_Vz_pri  ->Fill(primaryVertex_Z);
    hist_VyVx_pri->Fill(primaryVertex_X,primaryVertex_Y);
    hist_Vr_pri  ->Fill(primaryVertex_perp);

    // ---------------------- trigger selection ---------------------------------
    std::vector <unsigned int> triggerIDs;
    triggerIDs.clear();
    triggerIDs      = event->triggerIds();
    bool b_bad_trig = true;
    for(unsigned int i=0; i < triggerIDs.size(); i++)
      {
        Double_t d_trigger = (Double_t)triggerIDs[i] - 620050.0;
        hist_triggerID->Fill(d_trigger);
        if(triggerIDs[i] == 630052) b_bad_trig = false; // bbce_tofmult1 7.2GeV
      }

    // --------------------------- Vertex cut -----------------------------------
    double      d_zvtx  = -9999.0;
    double      d_xvtx  = -9999.0;
    double      d_yvtx  = -9999.0;
    double  d_vtx_perp  = -9999.0;
    TVector3       pVtx = event->primaryVertex();
    d_zvtx     = pVtx.z();
    d_xvtx     = pVtx.x();
    d_yvtx     = pVtx.y();
    d_vtx_perp = pVtx.Perp();
    bool b_bad_zvtx   =  ((d_zvtx < 198.0) || (d_zvtx > 202.0)); //FXT_26p5_2018
    bool b_bad_xvtx   =  ((d_xvtx < -1.0) || (d_xvtx > 1.0)); //FXT_26p5_2018
    bool b_bad_yvtx   =  ((d_yvtx < -3.0) || (d_yvtx > -0.5)); //FXT_26p5_2018
    bool b_bad_rvtx   =   sqrt(pow(d_xvtx,2)+pow(d_yvtx+2,2))> 2.0;
    bool b_bad_evt  = b_bad_zvtx || b_bad_trig /*|| b_bad_xvtx || b_bad_yvtx */|| b_bad_rvtx;
    if(b_bad_evt) continue;
    hist_Vz_cut->Fill(primaryVertex_Z);
    hist_Vr_cut->Fill(primaryVertex_perp);
    hist_VyVx_cut->Fill(primaryVertex_X,primaryVertex_Y);
    mEvtcut[1]++; // 1. vertex event cut

    Int_t  refMult = event->refMult(); // refMult
    Int_t grefMult = event->grefMult();
    Int_t  tofMult =(Int_t)event->nBTOFMatch();
    // (5) =============== Track loop to determine good tracks =================
    int nGoodTracks = 0;
    std::vector<StPicoTrack *> vGoodTracks; // vector of good tracks for TPC event plane Q-vector loop
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
      StPicoTrack *picoTrack = dst->track(iTrk);
      mTrkcut[0]++; // 0. No track cut
      if(!picoTrack) continue;
      mTrkcut[1]++; // 1. pico track cut
      StPicoBTofPidTraits *trait = NULL;
      // ----------------------- Physics values of tracks --------------------------
      double        tofBeta    = -999.;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      if(trait)        tofBeta = trait->btofBeta();
      double d_px  = picoTrack->pMom().x();
      double d_py  = picoTrack->pMom().y();
      double d_pz  = picoTrack->pMom().z();
      double d_pT  = picoTrack->pPt();
      double d_mom = sqrt(d_pT*d_pT + d_pz*d_pz);
      double mass2 = d_mom*d_mom*((1.0/(tofBeta*tofBeta))-1.0);
      Double_t eta = picoTrack->pMom().Eta();
      Double_t phi    = picoTrack->pMom().Phi();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();
      // --------------- QA plots before major track cuts ----------------------
      hist_px_py->Fill(d_px,d_py);
      hist_pz   ->Fill(d_pz);
      hist_pt   ->Fill(d_pT);
      hist_mom  ->Fill(d_mom);
      hist_mass2->Fill(mass2);
      hist_eta  ->Fill(eta);
      hist_phi  ->Fill(phi);
      hist_ratio->Fill(((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()));
      hist_nHits->Fill((double)picoTrack->nHitsFit());
      hist_ndEdx->Fill(picoTrack->nHitsDedx());
      hist_DCA  ->Fill(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z));

      if(!picoTrack->isPrimary()) continue;
      mTrkcut[2]++; // 2. Primary track cut
      bool    b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool    b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.51);
      bool b_not_enough_hits = ((double)picoTrack->nHitsFit()) < 15;
      bool    b_bad_DCA      = (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);
      bool    b_bad_track    = b_bad_dEdx || b_bad_tracking || b_not_enough_hits || b_bad_DCA;
      if(b_bad_track) continue;
      mTrkcut[3]++; // 3. Bad track cuts
      nGoodTracks++; // nGoodTracks is used to determine centrality later in the event loop
      vGoodTracks.push_back(picoTrack);
      // --------------- QA plots after major track cuts ----------------------
      hist_px_py_cut->Fill(d_px,d_py);
      hist_pz_cut   ->Fill(d_pz);
      hist_pt_cut   ->Fill(d_pT);
      hist_mom_cut  ->Fill(d_mom);
      hist_mass2_cut->Fill(mass2);
      hist_eta_cut  ->Fill(eta);
      hist_phi_cut  ->Fill(phi);
      hist_ratio_cut->Fill(((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()));
      hist_nHits_cut->Fill((double)picoTrack->nHitsFit());
      hist_ndEdx_cut->Fill(picoTrack->nHitsDedx());
      hist_DCA_cut  ->Fill(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z));
      if(tofBeta == -999) continue;
      mTrkcut[4]++; // 4. Bad tof track cut, to see how many tracks with tof information
    } // Track loop to determine good tracks
    for(int i=0;i<5;i++){ // fill the tracks after cut
      hist_trackCuts->SetBinContent(i+1,mTrkcut[i]);
    }
    // (6) ================ Centrality definition ===============================
    Int_t centrality = 0;
    bool a_b_cent[9]={false};
    Int_t cenSection[9]={11,22,37,57,82,113,151,174,245};//10,17,28,41,57,77,100,127,160,245 version 0 cent
    bool b_pileup   = (nGoodTracks > 245);
    bool b_low_mult = (nGoodTracks < 5);
    a_b_cent[0]     = (nGoodTracks >= cenSection[7] && nGoodTracks < cenSection[8]); // 0 - 5%
    a_b_cent[1]     = (nGoodTracks >= cenSection[6] && nGoodTracks < cenSection[7]); // 5 - 10%
    a_b_cent[2]     = (nGoodTracks >= cenSection[5] && nGoodTracks < cenSection[6]); // 10 - 20%
    a_b_cent[3]     = (nGoodTracks >= cenSection[4]  && nGoodTracks < cenSection[5]); // 20 - 30%
    a_b_cent[4]     = (nGoodTracks >= cenSection[3]  && nGoodTracks < cenSection[4]); // 30 - 40%
    a_b_cent[5]     = (nGoodTracks >= cenSection[2]  && nGoodTracks < cenSection[3]); // 40 - 50%
    a_b_cent[6]     = (nGoodTracks >= cenSection[1]  && nGoodTracks < cenSection[2]); // 50 - 60%
    a_b_cent[7]     = (nGoodTracks >= cenSection[0]  && nGoodTracks < cenSection[1]); // 60 - 70%
    a_b_cent[8]     = (nGoodTracks >= 5  && nGoodTracks < cenSection[0]); // 70 - 80%
    for(int i=0;i<_Ncentralities;i++){
      if(a_b_cent[i]) centrality = i+1;
    }
    hist_cent->Fill(centrality);
    hist_realTrackMult->Fill(nGoodTracks);
    hist_realTrackMult_refmult->Fill(nGoodTracks,refMult);
    hist_realTrackMult_grefmult->Fill(nGoodTracks,grefMult);
    hist_realTrackMult_tofmult->Fill(nGoodTracks,tofMult);
    if(b_pileup||b_low_mult) continue; //Pile/lowMult cut
    mEvtcut[2]++; // 2. Pile Up event cut

    // (7) ================ EPD event plane ====================================
    // (7.1) ------------- EPD ep from Mike Lisa's class StEpdEpFinder // removed due to redundancy
    // (7.2) ------------------- EPD EP by hand ---------------------------------
    // refer to Mike's StEpdEpFinder and Yang's BBC Ep
    Int_t N_Epd_east[5]={0}; //Count # of hits in each eta region /// indices: [etaBin]
    Double_t QrawEastSide[5][2]={0};       /// indices: [etaBin][x,y]
    // Double_t QphiWeightedEastSide[5][2]={0};       /// indices: [etaBin][x,y]
    Double_t PsiEastRaw[5]={-999.0,-999.0,-999.0,-999.0,-999.0};           /// indices: [etaBin]
    // Double_t PsiEastPhiWeighted[5]={-999.0,-999.0,-999.0,-999.0,-999.0};       /// indices: [etaBin]
    Double_t PsiEastShifted[5]={-999.0,-999.0,-999.0,-999.0,-999.0};       /// indices: [etaBin]
    for (int iEpdHit = 0; iEpdHit < mEpdHits->GetEntries(); iEpdHit++){
      StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[iEpdHit]);
      int tileId,ring,TT,PP,EW,ADC;
      float nMip;
    	tileId = epdHit->id();
    	EW = (tileId<0)?0:1;
      if(EW!=0) continue; // EPD east event plane needed
    	ring = epdHit->row();
    	TT = epdHit->tile();
    	PP = epdHit->position();
    	ADC = epdHit->adc();
    	nMip = epdHit->nMIP();   // malisa 20feb2019 - I have finally made the transition from ADC (next line) to truly nMip, now that calibrations are done.
      //      nMip = (TT<10)?(double)ADC/160.0:(double)ADC/115.0;
      if(PP==1 && TT==1) hist_nMip->Fill(nMip);
      if (nMip<mThresh) continue;
      double TileWeight = (nMip<mMax)?nMip:mMax;
      TVector3 StraightLine = mEpdGeom->TileCenter(tileId) - event->primaryVertex();
      double phi = StraightLine.Phi();
      double eta = StraightLine.Eta();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();
      hist_Epdeta->Fill(eta);
      hist_Epdphi->Fill(phi);
      profile2D_PpVsEta->Fill(eta,PP,TileWeight);
      h2_hits_PpVsEta->Fill(eta,PP);
      h2_nMip_eta_cent->Fill(eta,centrality,TileWeight);

      //---------------------------------
      // fill Phi Weight histograms to be used in next iteration (if desired)
      // Obviously, do this BEFORE phi weighting!
      //---------------------------------
      // for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){
      //   int etaBin = (int)wt.GetXaxis()->FindBin(fabs(eta));
      //   double etaWeight = (double)wt.GetBinContent(etaBin,EventTypeId+1);
        // if(etaWeight==1){
          // mPhiWeightOutput[EventTypeId]->Fill(phi,TileWeight);
          // for(int bin=1;bin<13;bin++) mPhiAveraged[EventTypeId]->Fill((double)bin*TMath::Pi()/6.0-0.1,TileWeight/12.0);
        // }
      // }
      //--------------------------------
      // now calculate Q-vectors
      //--------------------------------
      // double PhiWeightedTileWeight = TileWeight;
      for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){
        // if (mPhiWeightInput[EventTypeId]){
          // int phiBin = (int)mPhiWeightInput[EventTypeId]->GetXaxis()->FindBin(phi);
          // PhiWeightedTileWeight /= mPhiWeightInput[EventTypeId]->GetBinContent(phiBin); // Phi weighting :https://drupal.star.bnl.gov/STAR/blog/lisa/phi-weighting-and-optimizing-ring-weights-auau-27-gev
          // std::cout<<"Tile weight: "<< TileWeight ;
          // std::cout<<" Phi weighted tile weight: "<< PhiWeightedTileWeight<<std::endl;
        // }
        int etaBin = (int)wt.GetXaxis()->FindBin(fabs(eta));
        double etaWeight = (double)wt.GetBinContent(etaBin,EventTypeId+1);
        int v1etaBin = (int)v1WtaWt->GetXaxis()->FindBin(eta);
        double v1EtaWeight = (double)v1WtaWt->GetBinContent(v1etaBin,centrality);
        v1EtaWeight = 1.0; // disable v1 eta weighting
        if(v1EtaWeight == 0){
          std::cout<<"Centality is "<<centrality<<"\t"<< "eta : " << eta<<"\t"<<"eta weighting: " << v1EtaWeight << std::endl;
        }
        if(etaWeight>0.0) N_Epd_east[EventTypeId]++;
        double Cosine = cos(phi*(double)EpOrder);
        double Sine   = sin(phi*(double)EpOrder);
        QrawEastSide[EventTypeId][0] += etaWeight * v1EtaWeight * TileWeight * Cosine;
        QrawEastSide[EventTypeId][1] += etaWeight * v1EtaWeight * TileWeight * Sine;

        // QphiWeightedEastSide[EventTypeId][0]      += etaWeight * PhiWeightedTileWeight * Cosine;
        // QphiWeightedEastSide[EventTypeId][1]      += etaWeight * PhiWeightedTileWeight * Sine;
        if(etaWeight==1){
          h2_TtVsPp[EventTypeId]->Fill(PP,TT);
          h2_TtVsPpNmip[EventTypeId]->Fill(PP,TT,TileWeight);
          h2_TtVsPpHit[EventTypeId]->Fill(PP,TT);

        }
      }
    } // loop over EPD hits
    // Before going any farther, flip the sign of the 1st-order Q-vector on the East side.
    //  I want the rapidity-odd first-order event plane.
    // Comment this out if v1 eta weighting used
    for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){// Comment this out if v1 eta weighting used
      for (int xy=0; xy<2; xy++){
        QrawEastSide[EventTypeId][xy]           *= -1.0;
        // QphiWeightedEastSide[EventTypeId][xy]           *= -1.0;
      }
    }

    // To remove autocorrelation in EPD-3, calculate Qvector for each epd hit in EPD-3: -5.16 <= eta < -3.82
    std::map<int,TVector2> mpQvctrEpdSub;
    for (int iEpdHit = 0; iEpdHit < mEpdHits->GetEntries(); iEpdHit++){
      StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[iEpdHit]);
      int tileId,ring,TT,PP,EW,ADC;
      float nMip;
      tileId = epdHit->id();
      EW = (tileId<0)?0:1;
      if(EW!=0) continue; // EPD east event plane needed
      ring = epdHit->row();
      TT = epdHit->tile();
      PP = epdHit->position();
      ADC = epdHit->adc();
      nMip = epdHit->nMIP();   // malisa 20feb2019 - I have finally made the transition from ADC (next line) to truly nMip, now that calibrations are done.
      //      nMip = (TT<10)?(double)ADC/160.0:(double)ADC/115.0;
      if (nMip<mThresh) continue;
      double TileWeight = (nMip<mMax)?nMip:mMax;
      TVector3 StraightLine = mEpdGeom->TileCenter(tileId) - event->primaryVertex();
      double phi = StraightLine.Phi();
      double eta = StraightLine.Eta();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();
      //--------------------------------
      // now calculate Q-vectors for each hit in EPD-3
      //--------------------------------
      if(eta>=etaRange[2] && eta < etaRange[3]){ // EPD-3
        double QxEpdSub, QyEpdSub;
        TVector2 Qvec;
        double Cosine = cos(phi*(double)EpOrder);
        double Sine   = sin(phi*(double)EpOrder);
        QxEpdSub = QrawEastSide[3][0] + TileWeight * Cosine; // Since QrawEastSide[EventTypeId][xy] already times -1.0, here shoud "+ Qx_i" to remove autocorrelation
        QyEpdSub = QrawEastSide[3][1] + TileWeight * Sine; // Since QrawEastSide[EventTypeId][xy] already times -1.0, here should "+ Qy_i" to remove autocorrelation
        Qvec=TVector2(QxEpdSub,QyEpdSub);

        mpQvctrEpdSub.insert(pair<int, TVector2>(iEpdHit, Qvec));
      }
    } // loop over EPD hits
    //Print out the map and fill the PsiRawEpdSub map
    std::map<int, TVector2>::iterator itr;
    // std::cout << "\nThe map mpPsiRawEpdSub is : \n";
    // cout << "\tKEY\tELEMENT\n";
    std::map<int,double> mpPsiRawEpdSub;
    for (itr = mpQvctrEpdSub.begin(); itr != mpQvctrEpdSub.end(); itr++) { // insert a map of key: iEpdHit, value: PsiRawEpdSub
        // std::cout << '\t' << itr->first
        //      << '\t' << (double)(itr->second).X()
        //      << '\t' << (double)(itr->second).Y() << '\n';
        if(N_Epd_east[3]<5) continue; // EPD-3
        Double_t PsiRawEpdSub;
        if(QrawEastSide[3][0] || QrawEastSide[3][1] ){ // EPD-3
          PsiRawEpdSub = GetPsi((double)(itr->second).X(),(double)(itr->second).Y(),EpOrder);
          mpPsiRawEpdSub.insert(pair<int, double>(itr->first, PsiRawEpdSub));
          hist_Epd_Sub_psi_raw_ini->Fill(PsiRawEpdSub);
          // std::cout << '\t' << itr->first
          //      << '\t' << PsiRawEpdSub << '\n';
        }
    }
    // std::cout << std::endl;
    //---------------------------------
    // Calculate unshifted EP angles
    //---------------------------------
    for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){
      if(N_Epd_east[EventTypeId]<5) continue;
      if(QrawEastSide[EventTypeId][0] || QrawEastSide[EventTypeId][1] ){
        PsiEastRaw[EventTypeId] = GetPsi(QrawEastSide[EventTypeId][0],QrawEastSide[EventTypeId][1],EpOrder);
        // PsiEastPhiWeighted[EventTypeId] = GetPsi(QphiWeightedEastSide[EventTypeId][0],QphiWeightedEastSide[EventTypeId][1],EpOrder);
      }
    }
    for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){
      if(QrawEastSide[EventTypeId][0] || QrawEastSide[EventTypeId][1] )
      {
        hist2_Epd_east_Qy_Qx_raw_ini[EventTypeId]->Fill(QrawEastSide[EventTypeId][0],QrawEastSide[EventTypeId][1]);
        if(PsiEastRaw[EventTypeId]!=-999.0){
          hist_Epd_east_psi_raw_ini[EventTypeId]->Fill(PsiEastRaw[EventTypeId]);
          // hist_Epd_east_psi_Weighted_ini[EventTypeId]->Fill(PsiEastPhiWeighted[EventTypeId]);
        }
      }
    }
    // --------------------------- " Do the SHIFT thing " ------------------------
    // Fill the PsiShiftedEpdSub map: Key: iEpdHit, value: PsiShiftedEpdSub
    std::map<int, double>::iterator itr1;
    // std::cout << "\nThe map mpPsiShiftedEpdSub is : \n";
    // cout << "\tKEY\tELEMENT\n";
    std::map<int,double> mpPsiShiftedEpdSub;
    for (itr1 = mpPsiRawEpdSub.begin(); itr1 != mpPsiRawEpdSub.end(); itr1++) { // insert a map of key: iEpdHit, value: PsiRawEpdSub
        Double_t PsiRawEpdSub = (double)itr1->second ;
        Double_t PsiShiftedEpdSub = PsiRawEpdSub ;
        if(PsiShiftedEpdSub==-999.0) continue;
        if (mEpdShiftInput_sin[3] != 0 && mEpdShiftInput_cos[3]!= 0){
          if(QrawEastSide[3][0] || QrawEastSide[3][1] ){
            for (int i=1; i<=_EpTermsMaxIni; i++){
          	  double tmp = (double)(EpOrder*i);
              double sinAve = mEpdShiftInput_sin[3]->GetBinContent(i,centrality);
          	  double cosAve = mEpdShiftInput_cos[3]->GetBinContent(i,centrality);
          	  PsiShiftedEpdSub +=
          	    2.0*(cosAve*sin(tmp*PsiRawEpdSub) - sinAve*cos(tmp*PsiRawEpdSub))/tmp; // use raw EP rather than Phi weighing EP
          	}
            double AngleWrapAround = 2.0*TMath::Pi()/(double)EpOrder;
             if (PsiShiftedEpdSub<0) PsiShiftedEpdSub += AngleWrapAround;
              else if (PsiShiftedEpdSub>AngleWrapAround) PsiShiftedEpdSub -= AngleWrapAround;

            mpPsiShiftedEpdSub.insert(pair<int, double>(itr1->first, PsiShiftedEpdSub));
            hist_Epd_Sub_psi_Shifted_ini->Fill(PsiShiftedEpdSub);
            // std::cout << '\t' << itr1->first
            //      << '\t' << PsiShiftedEpdSub << '\n';
          }
        }
    }
    // std::cout << std::endl;
    for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){ //etaRange {-5.1,-4.2,-3.28,-2.87,-2.60}
        PsiEastShifted[EventTypeId] = PsiEastRaw[EventTypeId]; // use raw EP rather than Phi weighing EP
        if(PsiEastShifted[EventTypeId]==-999.0) continue;
        if (mEpdShiftInput_sin[EventTypeId] != 0 && mEpdShiftInput_cos[EventTypeId]!= 0){
          for (int i=1; i<=_EpTermsMaxIni; i++){
        	  double tmp = (double)(EpOrder*i);
            double sinAve = mEpdShiftInput_sin[EventTypeId]->GetBinContent(i,centrality);
        	  double cosAve = mEpdShiftInput_cos[EventTypeId]->GetBinContent(i,centrality);
        	  PsiEastShifted[EventTypeId] +=
        	    2.0*(cosAve*sin(tmp*PsiEastRaw[EventTypeId]) - sinAve*cos(tmp*PsiEastRaw[EventTypeId]))/tmp; // use raw EP rather than Phi weighing EP
        	}
	         double AngleWrapAround = 2.0*TMath::Pi()/(double)EpOrder;
  	        if (PsiEastShifted[EventTypeId]<0) PsiEastShifted[EventTypeId] += AngleWrapAround;
  	         else if (PsiEastShifted[EventTypeId]>AngleWrapAround) PsiEastShifted[EventTypeId] -= AngleWrapAround;
        }
        hist_Epd_east_psi_Shifted_ini[EventTypeId]->Fill(PsiEastShifted[EventTypeId]);
      }
      // --------------------------- Fill the Correlations among EPD sub EPs ------------------------
      pairs = -1;
      for(int i = 0; i<3;i++){ // Correlations between EPD EP 1, 2, 3, 4. 6 pairs of correlations
        for(int j=i+1;j<4;j++){
          pairs++;
          if(PsiEastShifted[i+1]!=-999.0&&PsiEastShifted[j+1]!=-999.0){
            for(int n=0; n<2; n++){
              profile_correlation_epd_east[n][pairs]->Fill(centrality,TMath::Cos((double)(n+1) * (PsiEastShifted[i+1] - PsiEastShifted[j+1] )));
            }
            correlation2D_epd_east[pairs]->Fill(PsiEastShifted[i+1],PsiEastShifted[j+1]);
          }
        }
      }

    // -------------------- "Shift correction histograms Output" ----------------
    // -------------------- "calculate shift histograms for a future run" ----------------
    for (int i=1; i<=_EpTermsMaxIni; i++){
      for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){//etaRange {-5.1,-4.2,-3.28,-2.87,-2.60}
        double tmp = (double)(EpOrder*i);
        if(PsiEastRaw[EventTypeId]==-999.0) continue;
        mEpdShiftOutput_sin[EventTypeId]->Fill(i,centrality,sin(tmp*PsiEastRaw[EventTypeId]));// use raw EP rather than Phi weighing EP
        mEpdShiftOutput_cos[EventTypeId]->Fill(i,centrality,cos(tmp*PsiEastRaw[EventTypeId]));// use raw EP rather than Phi weighing EP
      }
    }
    // (8) ================ TPC event plane : use identedfied particles ====================================
    // Fill kaon tracks for phi meson analysis
    std::vector<StPicoTrack *> v_KaonPlus_tracks;
    std::vector<StPicoTrack *> v_KaonMinus_tracks;
    // Define TPC EP parameters
    Int_t NTpcAll[2] = {0};
    Double_t QrawTpcAll[2][2]={0.0};       /// indices:[TPCetaRange] [x,y]
    Double_t PsiTpcAllRaw[2]={-999.0,-999.0};
    Double_t PsiTpcAllShifted[2]={-999.0,-999.0};
    std::vector<std::vector<Double_t> > vQrawTpcAll; // {{Xn,Yn},...} For TPC EP flow, one want to remove current track
    Int_t nProtons=0,nKaonPlus=0,nKaonMinus=0,nPionPlus=0,nPionMinus=0; // PID parameters
    // TPC Q-vector loop
    for(unsigned int i=0; i<vGoodTracks.size();i++){
      StPicoTrack* picoTrack = vGoodTracks[i];
      StPicoBTofPidTraits *trait        = NULL;
      Int_t charge;
      Double_t pt,pz,eta,ptot,phi;
      Double_t mass2 =-999.0,tofBeta =-999.0;
      Double_t rapWeight = 0.0; // Weight based on rapidity
      charge = picoTrack->charge();
      pt     = picoTrack->pPt();
      pz     = picoTrack->pMom().Z();
      eta    = picoTrack->pMom().Eta();
      ptot = picoTrack->pPtot();
      phi    = picoTrack->pMom().Phi();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();
      // ---------------- Check if TOF info available --------------------------
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      if(trait) tofBeta               = trait->btofBeta();
      if(tofBeta != -999.0) mass2 = ptot * ptot *( ( 1.0 / ( tofBeta*tofBeta ) ) - 1.0 );
      // ---------------- Particle Physics Variables presumed --------------------------
      Double_t energyProton,energyKaon,energyPion,rapProton,rapKaon,rapPion,mtProton,mtKaon,mtPion;
      Int_t particleType = -999; //default -999. 0,1,2,3,4 indicate p, K+, K-, \Pi+, \Pi-
      energyProton = TMath::Sqrt(ptot*ptot + _massProton*_massProton);
      energyKaon = TMath::Sqrt(ptot*ptot + _massKaon*_massKaon);
      energyPion = TMath::Sqrt(ptot*ptot + _massPion*_massPion);
      rapProton    = 0.5*TMath::Log( (energyProton + pz) / (energyProton - pz) );
      rapKaon    = 0.5*TMath::Log( (energyKaon + pz) / (energyKaon - pz) );
      rapPion    = 0.5*TMath::Log( (energyPion + pz) / (energyPion - pz) );
      mtProton   = TMath::Sqrt(pt*pt + _massProton*_massProton);
      mtKaon   = TMath::Sqrt(pt*pt + _massProton*_massProton);
      mtPion   = TMath::Sqrt(pt*pt + _massProton*_massProton);
      // ------------------------ TPC EP QA plots ------------------------------
      h2_dEdxVsPq->Fill(charge*ptot,picoTrack->dEdx());
      h2_dEdxVspTq->Fill(charge*pt,picoTrack->dEdx());
      if(tofBeta!=-999.0){
        h2_beta->Fill(charge*ptot,1.0/tofBeta);
        h2_mass->Fill(charge*ptot,mass2);
      }
      // ------------------------ Particle identifications ------------------------------
      if( // Proton PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaProton()) < 2.0 &&
        tofBeta != -999.0 && mass2 > 0.7 && mass2 < 1.1 &&
        pt > 0.2 &&
        charge > 0
      ){
        particleType=0;// Proton
        nProtons++;
        // Fill histograms
        hist_pt_proton->Fill(pt);
        hist_eta_proton->Fill(eta);
        hist_y_proton->Fill(rapProton);
        hist_phi_proton->Fill(phi);
        hist_rap_eta_proton->Fill(eta,rapProton);
        hist_pt_y_proton->Fill(rapProton,pt,1);
        hist_pt_eta_proton->Fill(eta,pt,1);
        hist_dEdx_proton->Fill(charge*ptot,picoTrack->dEdx());
        hist_beta_proton->Fill(charge*ptot,1.0/tofBeta);
        hist_mass_proton->Fill(charge*ptot,mass2);
      } else if( // Kaons PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaKaon()) < 2.0 &&
        tofBeta != -999.0 && mass2 > 0.16 && mass2 < 0.32
        && pt > 0.2
      ){
        if(charge > 0){
          particleType=1;// K+
          nKaonPlus++;
          v_KaonPlus_tracks.push_back(picoTrack); // push back K+ tracks
          // Fill histograms
          hist_pt_kaonPlus->Fill(pt);
          hist_eta_kaonPlus->Fill(eta);
          hist_y_kaonPlus->Fill(rapKaon);
          hist_phi_kaonPlus->Fill(phi);
          hist_rap_eta_kaonPlus->Fill(eta,rapKaon);
          hist_pt_y_kaonPlus->Fill(rapKaon,pt,1);
          hist_pt_eta_kaonPlus->Fill(eta,pt,1);
          hist_dEdx_kaonPlus->Fill(charge*ptot,picoTrack->dEdx());
          hist_beta_kaonPlus->Fill(charge*ptot,1.0/tofBeta);
          hist_mass_kaonPlus->Fill(charge*ptot,mass2);
        } else { // charge < 0
          particleType=2;// K-
          nKaonMinus++;
          v_KaonMinus_tracks.push_back(picoTrack); // push back K+ tracks
          // Fill histograms
          hist_pt_kaonMinus->Fill(pt);
          hist_eta_kaonMinus->Fill(eta);
          hist_y_kaonMinus->Fill(rapKaon);
          hist_phi_kaonMinus->Fill(phi);
          hist_rap_eta_kaonMinus->Fill(eta,rapKaon);
          hist_pt_y_kaonMinus->Fill(rapKaon,pt,1);
          hist_pt_eta_kaonMinus->Fill(eta,pt,1);
          hist_dEdx_kaonMinus->Fill(charge*ptot,picoTrack->dEdx());
          hist_beta_kaonMinus->Fill(charge*ptot,1.0/tofBeta);
          hist_mass_kaonMinus->Fill(charge*ptot,mass2);
        }
      } else if( // Pions PID: require both TPC and TOF
        TMath::Abs(picoTrack->nSigmaPion()) <  2.0 &&
        tofBeta != -999.0 && mass2 > -0.01 && mass2 < 0.05 &&
        pt > 0.2  &&
        !(TMath::Abs(mass2)<0.005 && ptot<0.25) // Remove electron influence
      ){
        if(charge > 0){
          particleType=3;// \Pi+
          nPionPlus++;
          // Fill histograms
          hist_pt_pionPlus->Fill(pt);
          hist_eta_pionPlus->Fill(eta);
          hist_y_pionPlus->Fill(rapPion);
          hist_phi_pionPlus->Fill(phi);
          hist_rap_eta_pionPlus->Fill(eta,rapPion);
          hist_pt_y_pionPlus->Fill(rapPion,pt,1);
          hist_pt_eta_pionPlus->Fill(eta,pt,1);
          hist_dEdx_pionPlus->Fill(charge*ptot,picoTrack->dEdx());
          hist_beta_pionPlus->Fill(charge*ptot,1.0/tofBeta);
          hist_mass_pionPlus->Fill(charge*ptot,mass2);
        } else { // charge < 0
          particleType=4;// \Pi-
          nPionMinus++;
          // Fill histograms
          hist_pt_pionMinus->Fill(pt);
          hist_eta_pionMinus->Fill(eta);
          hist_y_pionMinus->Fill(rapPion);
          hist_phi_pionMinus->Fill(phi);
          hist_rap_eta_pionMinus->Fill(eta,rapPion);
          hist_pt_y_pionMinus->Fill(rapPion,pt,1);
          hist_pt_eta_pionMinus->Fill(eta,pt,1);
          hist_dEdx_pionMinus->Fill(charge*ptot,picoTrack->dEdx());
          hist_beta_pionMinus->Fill(charge*ptot,1.0/tofBeta);
          hist_mass_pionMinus->Fill(charge*ptot,mass2);
        }
      }
      // if(particleType==-999) continue; // No particle identified
      // if(particleType==0) rapWeight= rapProton + 2.02; // y_CM = -2.02, COM rapidity
      // if(particleType==1||particleType==2) rapWeight= rapKaon + 2.02; // y_CM = -2.02, COM rapidity
      // if(particleType==3||particleType==4) rapWeight= rapPion + 2.02;// y_CM = -2.02, COM rapidity
      // Use all the good tracks to determine TPC EP
      double etaTrkWeight = 0.;
      if(eta>=_y_mid) {etaTrkWeight = 1.;} else{
        etaTrkWeight = -1;
      }
      for(int EventTypeId_tpc=0;EventTypeId_tpc<_nEventTypeBins_tpc;EventTypeId_tpc++){
        int etaBin = (int)wt_tpc.GetXaxis()->FindBin(fabs(eta));
        double etaWeight = (double)wt_tpc.GetBinContent(etaBin,EventTypeId_tpc+1);
        if(EpOrder == 1){ // \psi_1^{TPC}
          if(etaWeight>0.0 && etaTrkWeight /*rapWeight*/!=0) NTpcAll[EventTypeId_tpc]++;
          double Cosine = cos(phi*(double)EpOrder);
          double Sine   = sin(phi*(double)EpOrder);
          QrawTpcAll[EventTypeId_tpc][0] += etaWeight * etaTrkWeight /*rapWeight*/ * Cosine;
          QrawTpcAll[EventTypeId_tpc][1] += etaWeight * etaTrkWeight /*rapWeight*/ * Sine;
        } else { // \psi_2^{TPC}
          if(etaWeight>0.0) NTpcAll[EventTypeId_tpc]++;
          double Cosine = cos(phi*(double)EpOrder);
          double Sine   = sin(phi*(double)EpOrder);
          QrawTpcAll[EventTypeId_tpc][0] += etaWeight * pt * Cosine;
          QrawTpcAll[EventTypeId_tpc][1] += etaWeight * pt * Sine;
        }
      }
      // calculate the v1 in TPC region using EPD EP
      if(PsiEastShifted[1]!=-999.0){// Using EPD-1
        // ------------- Fill histograms for the determination of TPC eta range -----
        profile2D_v1VsEtaTpcOnly->Fill(eta,centrality,etaTrkWeight /*rapWeight*/ * TMath::Cos((phi-PsiEastShifted[1])*(Double_t)EpOrder));
        profile2D_v1VsEtaTpcOnly_1->Fill(eta,centrality,TMath::Cos((phi-PsiEastShifted[1])*(Double_t)EpOrder));
      // ------------------- Fill the eta weighting histograms --------------------------
        profile2D_v1VsCentVsEta->Fill(eta,centrality,TMath::Cos(phi-PsiEastShifted[1])/d_resolution[0][centrality-1]);//Use EPD-3 as primary event plane
        profile_v1VsEta[centrality-1]->Fill(eta,TMath::Cos(phi-PsiEastShifted[1])/d_resolution[0][centrality-1]); // [] is from 0 to 8, centrality is from 1 to 9.
      }
      hist_nTracksVsEta->Fill(eta,centrality);//histograms for the determination of TPC eta range
    } // TPC Q-vector loop


    // Track multiplicity for each particle
    hist_trackmult_proton->Fill(nProtons);
    hist_trackmult_pionPlus->Fill(nPionPlus);
    hist_trackmult_pionMinus->Fill(nPionMinus);
    hist_trackmult_kaonPlus->Fill(nKaonPlus);
    hist_trackmult_kaonMinus->Fill(nKaonMinus);
    //---------------------------------
    // Calculate unshifted EP angles
    //---------------------------------
    for(int EventTypeId_tpc=0;EventTypeId_tpc<_nEventTypeBins_tpc;EventTypeId_tpc++){
      if(NTpcAll[EventTypeId_tpc]<5) continue; // at least 5 tracks to get TPC event plane
      if(QrawTpcAll[EventTypeId_tpc][0] || QrawTpcAll[EventTypeId_tpc][1] ){ // Qx, Qy cannot be 0 at the same time
        PsiTpcAllRaw[EventTypeId_tpc] = (1./(Double_t)EpOrder)*TMath::ATan2(QrawTpcAll[EventTypeId_tpc][1],QrawTpcAll[EventTypeId_tpc][0]);
        if(PsiTpcAllRaw[EventTypeId_tpc] < 0.0                             )         PsiTpcAllRaw[EventTypeId_tpc] += (1. / (double)EpOrder) * 2.0*TMath::Pi();
        if(PsiTpcAllRaw[EventTypeId_tpc] > (1. / (double)EpOrder) * 2.0*TMath::Pi()) PsiTpcAllRaw[EventTypeId_tpc] -= (1. / (double)EpOrder) * 2.0*TMath::Pi();
        if(PsiTpcAllRaw[EventTypeId_tpc]!=-999.0) hist_tpc_all_psi_raw[EventTypeId_tpc]->Fill(PsiTpcAllRaw[EventTypeId_tpc]);
      }
    }
    // --------------------------- " Do the SHIFT thing (TPC) " ------------------------
    for(int EventTypeId_tpc=0;EventTypeId_tpc<_nEventTypeBins_tpc;EventTypeId_tpc++){
      PsiTpcAllShifted[EventTypeId_tpc] = PsiTpcAllRaw[EventTypeId_tpc];
      if(PsiTpcAllShifted[EventTypeId_tpc]==-999.0) continue; // Bad PsiTpcAllRaw
      if (mTpcShiftInput_sin[EventTypeId_tpc] != 0 && mTpcShiftInput_cos[EventTypeId_tpc]!= 0){
        for (int i=1; i<=_EpTermsMaxIni; i++){
            double tmp = (double)(EpOrder*i);
            double sinAve = mTpcShiftInput_sin[EventTypeId_tpc]->GetBinContent(i,centrality);
            double cosAve = mTpcShiftInput_cos[EventTypeId_tpc]->GetBinContent(i,centrality);
            PsiTpcAllShifted[EventTypeId_tpc] +=
              2.0*(cosAve*sin(tmp*PsiTpcAllRaw[EventTypeId_tpc]) - sinAve*cos(tmp*PsiTpcAllRaw[EventTypeId_tpc]))/tmp;
        }
           double AngleWrapAround = 2.0*TMath::Pi()/(double)EpOrder;
            if (PsiTpcAllShifted[EventTypeId_tpc]<0) PsiTpcAllShifted[EventTypeId_tpc] += AngleWrapAround;
             else if (PsiTpcAllShifted[EventTypeId_tpc]>AngleWrapAround) PsiTpcAllShifted[EventTypeId_tpc] -= AngleWrapAround;
      }
      hist_tpc_all_psi_shifted[EventTypeId_tpc]->Fill(PsiTpcAllShifted[EventTypeId_tpc]);
    }
    //---------------------------- Fill the directed flow from EPD (forward) region -----
    for (int iEpdHit = 0; iEpdHit < mEpdHits->GetEntries(); iEpdHit++){
      StPicoEpdHit* epdHit = (StPicoEpdHit*)((*mEpdHits)[iEpdHit]);
      int tileId,ring,TT,PP,EW,ADC;
      float nMip;
      tileId = epdHit->id();
      EW = (tileId<0)?0:1;
      if(EW!=0) continue; // EPD east event plane needed
      ring = epdHit->row();
      TT = epdHit->tile();
      PP = epdHit->position();
      ADC = epdHit->adc();
      nMip = epdHit->nMIP();   // malisa 20feb2019 - I have finally made the transition from ADC (next line) to truly nMip, now that calibrations are done.
      //      nMip = (TT<10)?(double)ADC/160.0:(double)ADC/115.0;
      if (nMip<mThresh) continue;
      double TileWeight = (nMip<mMax)?nMip:mMax;
      TVector3 StraightLine = mEpdGeom->TileCenter(tileId) - event->primaryVertex();
      double phi = StraightLine.Phi();
      double eta = StraightLine.Eta();
      if(phi < 0.0            ) phi += 2.0*TMath::Pi();
      if(phi > 2.0*TMath::Pi()) phi -= 2.0*TMath::Pi();

      //--------------------------------
      // Fill the directed flow into the TProfile2D and TProfile
      //--------------------------------
      // if(PsiTpcAllShifted[1]!=-999.0){//Use TPC EP for EPD v1 <(y-y_CM)*Cos(\phi - \psi_1)>
      //     profile2D_v1VsCentVsEta->Fill(eta,centrality,/*-(eta-_y_mid)**/TMath::Cos(phi-PsiTpcAllShifted[1]));//Use TPC
      //     profile_v1VsEta[centrality-1]->Fill(eta,/*-(eta-_y_mid)**/TMath::Cos(phi-PsiTpcAllShifted[1])); // [] is from 0 to 8, centrality is from 1 to 9.
      // }
      if( eta > etaRange[0] && eta < etaRange[1]){// Using EPD-1
        if(PsiEastShifted[3]!=-999.0){
          // ------------------- Fill the eta weighting histograms --------------------------
            profile2D_v1VsCentVsEta->Fill(eta,centrality,TMath::Cos(phi-PsiEastShifted[3])/d_resolution_EPD_3[centrality-1]);//Use EPD-3 as primary event plane
            profile_v1VsEta[centrality-1]->Fill(eta,TMath::Cos(phi-PsiEastShifted[3])/d_resolution_EPD_3[centrality-1]); // [] is from 0 to 8, centrality is from 1 to 9.
        }
      } else if(PsiEastShifted[1]!=-999.0){
        profile2D_v1VsCentVsEta->Fill(eta,centrality,TMath::Cos(phi-PsiEastShifted[1])/d_resolution[0][centrality-1]);//Use EPD-3 as primary event plane
        profile_v1VsEta[centrality-1]->Fill(eta,TMath::Cos(phi-PsiEastShifted[1])/d_resolution[0][centrality-1]); // [] is from 0 to 8, centrality is from 1 to 9.
      }
    } // loop over EPD hits
    // std::cout << std::endl;

    // ------------------- Fill the Correlations among TPC EP and EPD sub EPs ------------------------
    for(int n=0; n<2; n++){
      profile_correlation_epd_tpc_all[n]->Fill(centrality,TMath::Cos((double)(n+1) * (PsiEastShifted[0] - PsiTpcAllShifted[1] /*- TMath::Pi()/(double)EpOrder*/ )));
    }
    correlation2D_epd_tpc_all->Fill(PsiTpcAllShifted[1],PsiEastShifted[0]);
    for(int i=0;i<4;i++){// Correlaitons between TPC and EPD sub event planes 1,2,3,4
      if(PsiEastShifted[i+1]!=-999.0&&PsiTpcAllShifted[1]!=-999.0){
        for(int n=0; n<2; n++){
          profile_correlation_epd_tpc[n][i]->Fill(centrality,TMath::Cos((double)(n+1) * (PsiEastShifted[i+1] - PsiTpcAllShifted[1] /*- TMath::Pi()/(double)EpOrder*/ )));
        }
        correlation2D_epd_tpc[i]->Fill(PsiTpcAllShifted[1],PsiEastShifted[i+1]);
      }
    }
    // -------------------- "Shift correction histograms (TPC) Output" ----------------
    // -------------------- "calculate shift histograms for a future run" ----------------
    for(int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
      for (int i=1; i<=_EpTermsMaxIni; i++){ // TPC shifted Output
        double tmp = (double)(EpOrder*i);
        if(PsiTpcAllRaw[EventTypeId_tpc]==-999.0) break;
        mTpcShiftOutput_sin[EventTypeId_tpc]->Fill(i,centrality,sin(tmp*PsiTpcAllRaw[EventTypeId_tpc]));
        mTpcShiftOutput_cos[EventTypeId_tpc]->Fill(i,centrality,cos(tmp*PsiTpcAllRaw[EventTypeId_tpc]));
      }
    }
    // (9) ======================= Phi meson analysis  =========================
    double d_cut_mother_decay_length_PHI = 0.5; // must be LESS than this
    for(unsigned int i = 0; i < v_KaonPlus_tracks.size(); i++){
      StPicoTrack * picoTrack0 = v_KaonPlus_tracks.at(i); // i-th K+ track
      if(!picoTrack0) continue;
      // K+ Variables
      double d_charge0  = picoTrack0->charge();
      double d_px0      = picoTrack0->pMom().x();
      double d_py0      = picoTrack0->pMom().y();
      double d_pz0      = picoTrack0->pMom().z();
      double d_pT0      = picoTrack0->pPt();
      double d_mom0     = sqrt(d_pT0*d_pT0 + d_pz0*d_pz0);
      StPicoBTofPidTraits *trait0 = NULL;
      double d_tofBeta0    = -999.;
      double d_inv_tofBeta0    = -999.;
      if(picoTrack0->isTofTrack()) trait0 = dst->btofPidTraits( picoTrack0->bTofPidTraitsIndex() );
      if(trait0)        d_tofBeta0 = trait0->btofBeta();
      double d_M0   = _massKaon;
      double d_E0   = sqrt((d_px0*d_px0+d_py0*d_py0+d_pz0*d_pz0)+_massKaon*_massKaon);
      double d_y0   = ((d_E0-d_pz0) != 0.0) ? 0.5*TMath::Log( (d_E0 + d_pz0) / (d_E0 - d_pz0) ) : -999.0;
      double eta0   = ((d_mom0 - d_pz0) != 0.0) ? 0.5*TMath::Log( (d_mom0 + d_pz0) / (d_mom0 - d_pz0) ) : -999.0;
      double d_mT0  = sqrt(d_pT0*d_pT0 + d_M0*d_M0);
      double mass2_0 = d_mom0*d_mom0*((1.0/(d_tofBeta0*d_tofBeta0))-1.0);
      double d_pq0   = fabs(d_mom0) * d_charge0;
      if(d_tofBeta0 != -999. && d_tofBeta0 != 0.){
        d_inv_tofBeta0 = 1.0 / d_tofBeta0;
        h2_TOF_beta_pq  -> Fill(d_pq0,d_inv_tofBeta0);
      }
      for(unsigned int j = 0; j < v_KaonMinus_tracks.size(); j++){
        StPicoTrack * picoTrack1 = v_KaonMinus_tracks.at(j); // j-th K- track
        if(!picoTrack1) continue;
        // K- Variables
        double d_charge1  = picoTrack1->charge();
        if(d_charge0 == d_charge1) continue; // same charge cut
        double d_px1      = picoTrack1->pMom().x();
        double d_py1      = picoTrack1->pMom().y();
        double d_pz1      = picoTrack1->pMom().z();
        double d_pT1      = picoTrack1->pPt();
        double d_mom1     = sqrt(d_pT1*d_pT1 + d_pz1*d_pz1);
        StPicoBTofPidTraits *trait1 = NULL;
        double d_tofBeta1    = -999.;
        double d_inv_tofBeta1    = -999.;
        if(picoTrack1->isTofTrack()) trait1 = dst->btofPidTraits( picoTrack1->bTofPidTraitsIndex() );
        if(trait1)        d_tofBeta1 = trait1->btofBeta();
        double d_M1   = _massKaon;
        double d_E1   = sqrt((d_px1*d_px1+d_py1*d_py1+d_pz1*d_pz1)+_massKaon*_massKaon);
        double d_y1   = ((d_E1-d_pz1) != 0.0) ? 0.5*TMath::Log( (d_E1 + d_pz1) / (d_E1 - d_pz1) ) : -999.0;
        double eta1   = ((d_mom1 - d_pz1) != 0.0) ? 0.5*TMath::Log( (d_mom1 + d_pz1) / (d_mom1 - d_pz1) ) : -999.0;
        double d_mT1  = sqrt(d_pT1*d_pT1 + d_M1*d_M1);
        double mass2_1 = d_mom1*d_mom1*((1.0/(d_tofBeta1*d_tofBeta1))-1.0);
        double d_pq1   = fabs(d_mom0) * d_charge0;
        if(d_tofBeta1 != -999. && d_tofBeta1 != 0.){
          d_inv_tofBeta1 = 1.0 / d_tofBeta1;
          h2_TOF_beta_pq  -> Fill(d_pq1,d_inv_tofBeta1);
        }
        // phi Variables
        double d_dip_angle = TMath::ACos((d_pT0*d_pT1+d_pz0*d_pz1) / (d_mom0*d_mom1) );
        double d_Phi_pT = sqrt(d_px0*d_px0 + d_py0*d_py0 +d_px1*d_px1 +d_py1+d_py1 + 2.*d_px0*d_px1 + 2.*d_py0*d_py1);
        double d_mT_phi = sqrt(d_Phi_pT*d_Phi_pT + _massPhi*_massPhi );
        double d_phi_pz = d_pz0+d_pz1;
        double d_phi_E  = d_E0+d_E1;
        double d_phi_y  = ((d_phi_E - d_phi_pz) != 0.0) ?  0.5*TMath::Log( (d_phi_E + d_phi_pz) / (d_phi_E - d_phi_pz) ) : -9999;
        double d_inv_m  = sqrt(  d_M0*d_M0
                              + d_M1*d_M1
                              + 2.0 *d_E0*d_E1
                              - 2.0 *(d_px0*d_px1+d_py0*d_py1+d_pz0*d_pz1) );

        hist_dip_angle         ->Fill(d_dip_angle);
        hist_SE_mass_Phi    ->Fill(d_inv_m);
        hist_SE_PhiMeson_pT ->Fill(d_Phi_pT);
        hist_SE_PhiMeson_mT ->Fill(d_mT_phi);
        hist_SE_PhiMeson_rap ->Fill(d_phi_y);
        hist_SE_pt_y_PhiMeson ->Fill(d_phi_y,d_Phi_pT);
        // ---------------- phi-meson cuts: decay length, dip angle ------------
        StPicoPhysicalHelix    trackhelix0 = picoTrack0->helix(f_MagField);
        StPicoPhysicalHelix    trackhelix1 = picoTrack1->helix(f_MagField);
        pair<double,double> pairLengths = trackhelix0.pathLengths(trackhelix1);
        TVector3 v3D_x_daughter0 = trackhelix0.at(pairLengths.first);
        TVector3 v3D_x_daughter1 = trackhelix1.at(pairLengths.second);
        TVector3 v3D_x_mother    = (v3D_x_daughter0+v3D_x_daughter1)*0.5;
        TVector3 v3D_xvec_decayl = v3D_x_mother - pVtx;
        double d_mother_decay_length =  v3D_xvec_decayl.Mag();
        hist_mother_decay_length->Fill(d_mother_decay_length);
        // if(d_mother_decay_length > d_cut_mother_decay_length_PHI) continue; //decay length cut
        if(d_dip_angle<0.04) continue; // dip-angle cut
        // --------------------- phi-meson flows -------------------------------
        TVector3 v3D_p_daughter0 = trackhelix0.momentumAt(pairLengths.first, f_MagField*kilogauss);
        TVector3 v3D_p_daughter1 = trackhelix1.momentumAt(pairLengths.second, f_MagField*kilogauss);
        TVector3 v3D_p_mother    = v3D_p_daughter0+v3D_p_daughter1;
        double d_pmom = v3D_xvec_decayl.Dot(v3D_p_mother);
        double d_dca_mother = sqrt(v3D_xvec_decayl.Mag2() - (d_pmom*d_pmom/v3D_p_mother.Mag2()) );
        double d_phi_azimuth = v3D_p_mother.Phi();
        if(d_phi_azimuth < 0.0            ) d_phi_azimuth += 2.0*TMath::Pi();
        if(d_phi_azimuth > 2.0*TMath::Pi()) d_phi_azimuth -= 2.0*TMath::Pi();
        double d_flow_PHI_raw[2] = {-999.0,-999.0}; // v1, v2 raw flow
        double d_flow_PHI_resolution[2] = {-999.0,-999.0}; // v1, v2 flow corrected by resolution
        if(PsiEastShifted[1]!=-999.0){// Using EPD-1
          for(int km=0;km<2;km++){ // km - flow order
            d_flow_PHI_raw[km]        = TMath::Cos((double)(km+1.) * (d_phi_azimuth - PsiEastShifted[1]));
            d_flow_PHI_resolution[km] = TMath::Cos((double)(km+1.) * (d_phi_azimuth - PsiEastShifted[1]))/(d_resolution[km][centrality-1]); // km {0,1}, centrality [1,9]
          }
        }
        // -------------------- (9.1) Fill SE InvM plots -------------------------
        for(int pt=0; pt<2; pt++)
        {// pt SetA, cent SetA
          if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
            if(centrality >= 1 && centrality <= 2){
              mHist_SE_InvM_ptSetA_centSetA[pt][0]->Fill(d_inv_m); // 0-10%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][0]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-10%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][0]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-10%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][0]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-10%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][0]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-10%
            }
            if(centrality >= 3 && centrality <= 5){
              mHist_SE_InvM_ptSetA_centSetA[pt][1]->Fill(d_inv_m); // 10-40%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][1]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 10-40%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][1]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 10-40%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][1]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-10%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][1]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 10-40%
            }
            if(centrality >= 6 && centrality <= 7){
              mHist_SE_InvM_ptSetA_centSetA[pt][2]->Fill(d_inv_m); // 40-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][2]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][2]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][2]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][2]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-60%
            }
            if(centrality >= 6 && centrality <= 9){
              mHist_SE_InvM_ptSetA_centSetA[pt][3]->Fill(d_inv_m); // 40-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][3]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][3]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][3]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][3]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-80%
            }
            if(centrality >= 1 && centrality <= 7){
              mHist_SE_InvM_ptSetA_centSetA[pt][4]->Fill(d_inv_m); // 0-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][4]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][4]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][4]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][4]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-60%
            }
            if(centrality >= 1 && centrality <= 9){
              mHist_SE_InvM_ptSetA_centSetA[pt][5]->Fill(d_inv_m); // 0-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetA[pt][5]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetA[pt][5]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetA[pt][5]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetA[pt][5]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-80%
            }
          }
        }
        for(int pt=0; pt<2; pt++)
        {// pt SetA, cent SetB
          for(int cent=0; cent<9;cent++){
            if(d_Phi_pT >= ptSetA[pt] && d_Phi_pT <= ptSetA[pt+1]){
              if(centrality == cent+1 ){
                mHist_SE_InvM_ptSetA_centSetB[pt][cent]->Fill(d_inv_m);
                if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetA_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[0]);
                if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetA_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[0]);
                if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetA_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[1]);
                if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetA_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[1]);
              }
            }
          }
        }
        for(int i=0; i<4; i++)
        {// pt SetB, cent SetA
          if(d_Phi_pT >= ptSetB[i] && d_Phi_pT <= ptSetB[i+1]){
            if(centrality >= 1 && centrality <= 2){
              mHist_SE_InvM_ptSetB_centSetA[i][0]->Fill(d_inv_m); // 0-10%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-10%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-10%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-10%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-10%
            }
            if(centrality >= 3 && centrality <= 5){
              mHist_SE_InvM_ptSetB_centSetA[i][1]->Fill(d_inv_m); // 10-40%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 10-40%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 10-40%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 10-40%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 10-40%
            }
            if(centrality >= 6 && centrality <= 7){
              mHist_SE_InvM_ptSetB_centSetA[i][2]->Fill(d_inv_m); // 40-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-60%
            }
            if(centrality >= 6 && centrality <= 9){
              mHist_SE_InvM_ptSetB_centSetA[i][3]->Fill(d_inv_m); // 40-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-80%
            }
            if(centrality >= 1 && centrality <= 7){
              mHist_SE_InvM_ptSetB_centSetA[i][4]->Fill(d_inv_m); // 0-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-60%
            }
            if(centrality >= 1 && centrality <= 9){
              mHist_SE_InvM_ptSetB_centSetA[i][5]->Fill(d_inv_m); // 0-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-80%
            }
          }

          // rap SetA, cent SetA
          if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
            if(centrality >= 1 && centrality <= 2){
              mHist_SE_InvM_rapSetA_centSetA[i][0]->Fill(d_inv_m); // 0-10%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-10%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-10%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-10%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][0]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-10%
            }
            if(centrality >= 3 && centrality <= 5){
              mHist_SE_InvM_rapSetA_centSetA[i][1]->Fill(d_inv_m); // 10-40%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 10-40%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 10-40%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 10-40%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][1]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 10-40%
            }
            if(centrality >= 6 && centrality <= 7){
              mHist_SE_InvM_rapSetA_centSetA[i][2]->Fill(d_inv_m); // 40-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][2]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-60%
            }
            if(centrality >= 6 && centrality <= 9){
              mHist_SE_InvM_rapSetA_centSetA[i][3]->Fill(d_inv_m); // 40-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 40-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 40-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 40-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][3]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 40-80%
            }
            if(centrality >= 1 && centrality <= 7){
              mHist_SE_InvM_rapSetA_centSetA[i][4]->Fill(d_inv_m); // 0-60%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-60%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-60%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-60%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][4]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-60%
            }
            if(centrality >= 1 && centrality <= 9){
              mHist_SE_InvM_rapSetA_centSetA[i][5]->Fill(d_inv_m); // 0-80%
              if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_raw[0]); // 0-80%
              if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_resolution[0]); // 0-80%
              if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_raw[1]); // 0-80%
              if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetA[i][5]->Fill(d_inv_m,d_flow_PHI_resolution[1]); // 0-80%
            }
          }
          // rap SetA, cent SetB
          for(int cent=0; cent<9;cent++){
            if(d_phi_y >= rapSetA[i] && d_phi_y <= rapSetA[i+1]){
              if(centrality == cent+1 ){
                mHist_SE_InvM_rapSetA_centSetB[i][cent]->Fill(d_inv_m);
                if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_rapSetA_centSetB[i][cent]->Fill(d_inv_m,d_flow_PHI_raw[0]);
                if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_rapSetA_centSetB[i][cent]->Fill(d_inv_m,d_flow_PHI_resolution[0]);
                if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_rapSetA_centSetB[i][cent]->Fill(d_inv_m,d_flow_PHI_raw[1]);
                if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_rapSetA_centSetB[i][cent]->Fill(d_inv_m,d_flow_PHI_resolution[1]);
              }
            }
          }
        }
        for(int pt=0; pt<4; pt++)
        {// pt SetB, cent SetB
          for(int cent=0; cent<9;cent++){
            if(d_Phi_pT >= ptSetB[pt] && d_Phi_pT <= ptSetB[pt+1]){
              if(centrality == cent+1 ){
                mHist_SE_InvM_ptSetB_centSetB[pt][cent]->Fill(d_inv_m);
                if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetB_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[0]);
                if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetB_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[0]);
                if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetB_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[1]);
                if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetB_centSetB[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[1]);
              }
            }
          }
        }
        for(int pt=0; pt<10; pt++)
        {// pt SetC, cent 0-60%, 0-80%
          for(int cent=0; cent<2;cent++){
            if(d_Phi_pT >= ptSetC[pt] && d_Phi_pT <= ptSetC[pt+1]){
              if(centrality >= 1 && centrality <= cent*2 + 7 ){
                mHist_SE_InvM_ptSetC_centAll[pt][cent]->Fill(d_inv_m);
                if(d_flow_PHI_raw[0]!=-999.0)        mHist_v1_raw_ptSetC_centAll[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[0]);
                if(d_flow_PHI_resolution[0]!=-999.0) mHist_v1_reso_ptSetC_centAll[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[0]);
                if(d_flow_PHI_raw[1]!=-999.0)        mHist_v2_raw_ptSetC_centAll[pt][cent]->Fill(d_inv_m,d_flow_PHI_raw[1]);
                if(d_flow_PHI_resolution[1]!=-999.0) mHist_v2_reso_ptSetC_centAll[pt][cent]->Fill(d_inv_m,d_flow_PHI_resolution[1]);
              }
            }
          }
        }


      }
    }
    v_KaonPlus_tracks.clear();
    v_KaonMinus_tracks.clear();
  }  // Event Loop
  // --------------------- Set histograms axises titles --------------------------------
  hist_runId->GetXaxis()->SetTitle("RunId");
  hist_runId->GetYaxis()->SetTitle("# of events");
  hist_eventCuts->SetBinContent(1,mEvtcut[0]);
  hist_eventCuts->SetBinContent(2,mEvtcut[1]);
  hist_eventCuts->SetBinContent(3,mEvtcut[2]);
  hist_eventCuts->GetXaxis()->SetBinLabel(1,"no cuts");
  hist_eventCuts->GetXaxis()->SetBinLabel(2,"Vertex cuts");
  hist_eventCuts->GetXaxis()->SetBinLabel(2,"Pile up/lowMult cut");
  hist_trackCuts->GetXaxis()->SetBinLabel(1,"no cuts");
  hist_trackCuts->GetXaxis()->SetBinLabel(2,"picoTrack");
  hist_trackCuts->GetXaxis()->SetBinLabel(3,"primary track");
  hist_trackCuts->GetXaxis()->SetBinLabel(4,"Good track");
  hist_trackCuts->GetXaxis()->SetBinLabel(4,"With TOF");
  hist_Vz_pri->GetXaxis()->SetTitle("V_{Z} [cm]");
  hist_Vz_pri->GetYaxis()->SetTitle("# of events");
  hist_VyVx_pri->GetXaxis()->SetTitle("V_{X} [cm]");
  hist_VyVx_pri->GetYaxis()->SetTitle("V_{Y} [cm]");
  hist_Vr_pri->GetXaxis()->SetTitle("V_{R} [cm]");
  hist_Vr_pri->GetYaxis()->SetTitle("# of events");
  hist_triggerID->GetXaxis()->SetTitle("TriggerID");
  hist_triggerID->GetYaxis()->SetTitle("# of events");
  hist_Vz_cut->GetXaxis()->SetTitle("V_{Z} [cm]");
  hist_Vz_cut->GetYaxis()->SetTitle("# of events");
  hist_Vr_cut->GetXaxis()->SetTitle("V_{R} [cm]");
  hist_Vr_cut->GetYaxis()->SetTitle("# of events");
  hist_VyVx_cut->GetXaxis()->SetTitle("V_{X} [cm]");
  hist_VyVx_cut->GetYaxis()->SetTitle("V_{Y} [cm]");
  hist_px_py->GetXaxis()->SetTitle("p_{x} [GeV/c]");
  hist_px_py->GetYaxis()->SetTitle("p_{y} [GeV/c]");
  hist_pz->GetXaxis()->SetTitle("p_{z} [GeV/c]");
  hist_pz->GetYaxis()->SetTitle("# of tracks");
  hist_pt->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt->GetYaxis()->SetTitle("# of tracks");
  hist_mom->GetXaxis()->SetTitle("p_{mom} [GeV/c]");
  hist_mom->GetYaxis()->SetTitle("# of tracks");
  hist_mass2->GetXaxis()->SetTitle("M^{2} [(GeV/c^{2})^{2}]");
  hist_mass2->GetYaxis()->SetTitle("# of tracks");
  hist_eta->GetXaxis()->SetTitle("#eta");
  hist_eta->GetYaxis()->SetTitle("# of tracks");
  hist_phi->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi->GetYaxis()->SetTitle("# of tracks");
  hist_ratio->GetXaxis()->SetTitle("nHitsFit/nHitsPoss");
  hist_ratio->GetYaxis()->SetTitle("# of Tracks");
  hist_nHits->GetXaxis()->SetTitle("nHits");
  hist_nHits->GetYaxis()->SetTitle("# of Tracks");
  hist_ndEdx->GetXaxis()->SetTitle("nDedx");
  hist_ndEdx->GetYaxis()->SetTitle("# of Tracks");
  hist_DCA->GetXaxis()->SetTitle("DCA [cm]");
  hist_DCA->GetYaxis()->SetTitle("# of Tracks");
  hist_px_py_cut->GetXaxis()->SetTitle("p_{x} [GeV/c]");
  hist_px_py_cut->GetYaxis()->SetTitle("p_{y} [GeV/c]");
  hist_pz_cut->GetXaxis()->SetTitle("p_{z} [GeV/c]");
  hist_pz_cut->GetYaxis()->SetTitle("# of tracks");
  hist_pt_cut->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_cut->GetYaxis()->SetTitle("# of tracks");
  hist_mom_cut->GetXaxis()->SetTitle("p_{mom} [GeV/c]");
  hist_mom_cut->GetYaxis()->SetTitle("# of tracks");
  hist_mass2_cut->GetXaxis()->SetTitle("M^{2} [(GeV/c^{2})^{2}]");
  hist_mass2_cut->GetYaxis()->SetTitle("# of tracks");
  hist_eta_cut->GetXaxis()->SetTitle("#eta");
  hist_eta_cut->GetYaxis()->SetTitle("# of tracks");
  hist_phi_cut->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_cut->GetYaxis()->SetTitle("# of tracks");
  hist_ratio_cut->GetXaxis()->SetTitle("nHitsFit/nHitsPoss");
  hist_ratio_cut->GetYaxis()->SetTitle("# of Tracks");
  hist_nHits_cut->GetXaxis()->SetTitle("nHits");
  hist_nHits_cut->GetYaxis()->SetTitle("# of Tracks");
  hist_ndEdx_cut->GetXaxis()->SetTitle("nDedx");
  hist_ndEdx_cut->GetYaxis()->SetTitle("# of Tracks");
  hist_DCA_cut->GetXaxis()->SetTitle("DCA [cm]");
  hist_DCA_cut->GetYaxis()->SetTitle("# of Tracks");
  hist_cent->GetXaxis()->SetTitle("Centrality bin");
  hist_cent->GetYaxis()->SetTitle("# of events");
  hist_realTrackMult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult->GetYaxis()->SetTitle("# of events");
  hist_realTrackMult_refmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_refmult->GetYaxis()->SetTitle("RefMult");
  hist_realTrackMult_grefmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_grefmult->GetYaxis()->SetTitle("gRefMult");
  hist_realTrackMult_tofmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_tofmult->GetYaxis()->SetTitle("tofMult");
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    hist2_Epd_east_Qy_Qx_raw_ini[EventTypeId]->GetXaxis()->SetTitle("Q_x^{EPD east}_{1} ");
    hist2_Epd_east_Qy_Qx_raw_ini[EventTypeId]->GetYaxis()->SetTitle("Q_y^{EPD east}_{1} ");
    hist_Epd_east_psi_raw_ini[EventTypeId]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_raw_ini[EventTypeId]->GetYaxis()->SetTitle("# of events");
    // hist_Epd_east_psi_Weighted_ini[EventTypeId]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    // hist_Epd_east_psi_Weighted_ini[EventTypeId]->GetYaxis()->SetTitle("# of events");
    hist_Epd_east_psi_Shifted_ini[EventTypeId]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_Shifted_ini[EventTypeId]->GetYaxis()->SetTitle("# of events");
  }
  hist_Epdeta->GetXaxis()->SetTitle("#eta");
  hist_Epdeta->GetYaxis()->SetTitle("# of hits");
  hist_Epdphi->GetXaxis()->SetTitle("#phi [Radian]");
  hist_Epdphi->GetYaxis()->SetTitle("# of hits");
  profile2D_PpVsEta->GetXaxis()->SetTitle("#eta");
  profile2D_PpVsEta->GetYaxis()->SetTitle("Supersector");
  hist_nMip->GetXaxis()->SetTitle("nMIP");
  hist_nMip->GetYaxis()->SetTitle("# of hits");
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    h2_TtVsPp[EventTypeId]->GetXaxis()->SetTitle("Supersector");
    h2_TtVsPp[EventTypeId]->GetYaxis()->SetTitle("Tile");
    h2_TtVsPpNmip[EventTypeId]->GetXaxis()->SetTitle("Supersector");
    h2_TtVsPpNmip[EventTypeId]->GetYaxis()->SetTitle("Tile");
  }
  h2_dEdxVsPq->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  h2_dEdxVsPq->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  h2_dEdxVspTq->GetXaxis()->SetTitle("q*|pT| (GeV/c)");
  h2_dEdxVspTq->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  h2_beta->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  h2_beta->GetYaxis()->SetTitle("1/#beta");
  h2_mass->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  h2_mass->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_trackmult_proton->GetXaxis()->SetTitle("# of Protons");
  hist_trackmult_proton->GetYaxis()->SetTitle("# of Events");
  hist_trackmult_pionPlus->GetXaxis()->SetTitle("# of #pi^{#plus}");
  hist_trackmult_pionPlus->GetYaxis()->SetTitle("# of Events");
  hist_trackmult_pionMinus->GetXaxis()->SetTitle("# of #pi^{#minus}");
  hist_trackmult_pionMinus->GetYaxis()->SetTitle("# of Events");
  hist_trackmult_kaonPlus->GetXaxis()->SetTitle("# of K^{#plus}");
  hist_trackmult_kaonPlus->GetYaxis()->SetTitle("# of Events");
  hist_trackmult_kaonMinus->GetXaxis()->SetTitle("# of K^{#minus}");
  hist_trackmult_kaonMinus->GetYaxis()->SetTitle("# of Events");
  hist_pt_proton->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_proton->GetYaxis()->SetTitle("# of tracks");
  hist_eta_proton->GetXaxis()->SetTitle("#eta");
  hist_eta_proton->GetYaxis()->SetTitle("# of tracks");
  hist_y_proton->GetXaxis()->SetTitle("Rapidity y");
  hist_y_proton->GetYaxis()->SetTitle("# of tracks");
  hist_phi_proton->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_proton->GetYaxis()->SetTitle("# of tracks");
  hist_rap_eta_proton->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_proton->GetYaxis()->SetTitle("Rapidity y");
  hist_pt_y_proton->GetXaxis()->SetTitle("y");
  hist_pt_y_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_eta_proton->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_proton->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_dEdx_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_proton->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  hist_beta_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_proton->GetYaxis()->SetTitle("1/#beta");
  hist_mass_proton->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_proton->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_mass_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_kaonPlus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_beta_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_kaonPlus->GetYaxis()->SetTitle("1/#beta");
  hist_dEdx_kaonPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_kaonPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  hist_pt_eta_kaonPlus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_y_kaonPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_rap_eta_kaonPlus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_kaonPlus->GetYaxis()->SetTitle("Rapidity y");
  hist_phi_kaonPlus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_kaonPlus->GetYaxis()->SetTitle("# of tracks");
  hist_y_kaonPlus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_kaonPlus->GetYaxis()->SetTitle("# of tracks");
  hist_eta_kaonPlus->GetXaxis()->SetTitle("#eta");
  hist_eta_kaonPlus->GetYaxis()->SetTitle("# of tracks");
  hist_pt_kaonPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_kaonPlus->GetYaxis()->SetTitle("# of tracks");
  hist_mass_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_kaonMinus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_beta_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_kaonMinus->GetYaxis()->SetTitle("1/#beta");
  hist_dEdx_kaonMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_kaonMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  hist_pt_eta_kaonMinus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_y_kaonMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_kaonMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_rap_eta_kaonMinus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_kaonMinus->GetYaxis()->SetTitle("Rapidity y");
  hist_phi_kaonMinus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_kaonMinus->GetYaxis()->SetTitle("# of tracks");
  hist_y_kaonMinus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_kaonMinus->GetYaxis()->SetTitle("# of tracks");
  hist_eta_kaonMinus->GetXaxis()->SetTitle("#eta");
  hist_eta_kaonMinus->GetYaxis()->SetTitle("# of tracks");
  hist_pt_kaonMinus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_kaonMinus->GetYaxis()->SetTitle("# of tracks");
  hist_pt_pionPlus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_pionPlus->GetYaxis()->SetTitle("# of tracks");
  hist_eta_pionPlus->GetXaxis()->SetTitle("#eta");
  hist_eta_pionPlus->GetYaxis()->SetTitle("# of tracks");
  hist_y_pionPlus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_pionPlus->GetYaxis()->SetTitle("# of tracks");
  hist_phi_pionPlus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_pionPlus->GetYaxis()->SetTitle("# of tracks");
  hist_rap_eta_pionPlus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_pionPlus->GetYaxis()->SetTitle("Rapidity y");
  hist_pt_y_pionPlus->GetXaxis()->SetTitle("y");
  hist_pt_y_pionPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_eta_pionPlus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_pionPlus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_dEdx_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_pionPlus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  hist_beta_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_pionPlus->GetYaxis()->SetTitle("1/#beta");
  hist_mass_pionPlus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_pionPlus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_mass_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_mass_pionMinus->GetYaxis()->SetTitle("m^{2} (GeV/c^{2})^{2}");
  hist_beta_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_beta_pionMinus->GetYaxis()->SetTitle("1/#beta");
  hist_dEdx_pionMinus->GetXaxis()->SetTitle("q*|p| (GeV/c)");
  hist_dEdx_pionMinus->GetYaxis()->SetTitle("dE/dx (keV/cm)");
  hist_pt_eta_pionMinus->GetXaxis()->SetTitle("#eta");
  hist_pt_eta_pionMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_y_pionMinus->GetXaxis()->SetTitle("y");
  hist_pt_y_pionMinus->GetYaxis()->SetTitle("p_{T} [GeV/c]");
  hist_rap_eta_pionMinus->GetXaxis()->SetTitle("Pseudorapidity #eta");
  hist_rap_eta_pionMinus->GetYaxis()->SetTitle("Rapidity y");
  hist_phi_pionMinus->GetXaxis()->SetTitle("#phi [Radian]");
  hist_phi_pionMinus->GetYaxis()->SetTitle("# of tracks");
  hist_y_pionMinus->GetXaxis()->SetTitle("Rapidity y");
  hist_y_pionMinus->GetYaxis()->SetTitle("# of tracks");
  hist_eta_pionMinus->GetXaxis()->SetTitle("#eta");
  hist_eta_pionMinus->GetYaxis()->SetTitle("# of tracks");
  hist_pt_pionMinus->GetXaxis()->SetTitle("p_{T} [GeV/c]");
  hist_pt_pionMinus->GetYaxis()->SetTitle("# of tracks");
  hist_nTracksVsEta->GetXaxis()->SetTitle("#eta");
  hist_nTracksVsEta->GetYaxis()->SetTitle("centrality");
  profile2D_v1VsEtaTpcOnly->GetXaxis()->SetTitle("#eta");
  profile2D_v1VsEtaTpcOnly->GetYaxis()->SetTitle("centrality");
  profile2D_v1VsEtaTpcOnly_1->GetXaxis()->SetTitle("#eta");
  profile2D_v1VsEtaTpcOnly_1->GetYaxis()->SetTitle("centrality");
  pairs =0;
  for(int i = 0; i<3;i++){ // Correlations between EPD EP 1, 2, 3, 4. 6 pairs of correlations
    for(int j=i+1;j<4;j++){
      correlation2D_epd_east[pairs] ->GetXaxis()->SetTitle(Form("#phi of EPD%d",i+1));
      correlation2D_epd_east[pairs] ->GetYaxis()->SetTitle(Form("#phi of EPD%d",j+1));
      pairs++;
    }
  }
  for(int i=0;i<4;i++){// Correlaitons between TPC and EPD sub event planes 1,2,3,4
    correlation2D_epd_tpc[i] ->GetXaxis()->SetTitle("#phi of TPC");
    correlation2D_epd_tpc[i] ->GetYaxis()->SetTitle(Form("#phi of EPD%d",i+1));
  }
  for(int EventTypeId=0; EventTypeId<_nEventTypeBins; EventTypeId++){
    mEpdShiftOutput_sin[EventTypeId]->GetXaxis()->SetTitle("Shift order");
    mEpdShiftOutput_sin[EventTypeId]->GetYaxis()->SetTitle("Centrality");
    mEpdShiftOutput_cos[EventTypeId]->GetXaxis()->SetTitle("Shift order");
    mEpdShiftOutput_cos[EventTypeId]->GetYaxis()->SetTitle("Centrality");
  }
  for(int EventTypeId_tpc=0; EventTypeId_tpc<_nEventTypeBins_tpc; EventTypeId_tpc++){
    mTpcShiftOutput_sin[EventTypeId_tpc]->GetXaxis()->SetTitle("Shift order");
    mTpcShiftOutput_sin[EventTypeId_tpc]->GetYaxis()->SetTitle("Centrality");
    mTpcShiftOutput_cos[EventTypeId_tpc]->GetXaxis()->SetTitle("Shift order");
    mTpcShiftOutput_cos[EventTypeId_tpc]->GetYaxis()->SetTitle("Centrality");
  }
  profile2D_v1VsCentVsEta->GetXaxis()->SetTitle("#eta");
  profile2D_v1VsCentVsEta->GetYaxis()->SetTitle("centrality (%)");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(1,"0-5");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(2,"5-10");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(3,"10-20");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(4,"20-30");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(5,"30-40");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(6,"40-50");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(7,"50-60");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(8,"60-70");
  profile2D_v1VsCentVsEta->GetYaxis()->SetBinLabel(9,"70-80");

  // -------------------------------- Set titles -------------------------------
  for(int i=0; i<4; i++)
  {// pt SetB, cent SetA
    mHist_SE_InvM_ptSetB_centSetA[i][0]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_ptSetB_centSetA[i][1]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_ptSetB_centSetA[i][2]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_ptSetB_centSetA[i][3]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_ptSetB_centSetA[i][4]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_ptSetB_centSetA[i][5]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));

    mHist_v1_raw_ptSetB_centSetA[i][0]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_v1_raw_ptSetB_centSetA[i][1]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_v1_raw_ptSetB_centSetA[i][2]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_v1_raw_ptSetB_centSetA[i][3]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_v1_raw_ptSetB_centSetA[i][4]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_v1_raw_ptSetB_centSetA[i][5]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));

    mHist_v1_reso_ptSetB_centSetA[i][0]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_v1_reso_ptSetB_centSetA[i][1]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_v1_reso_ptSetB_centSetA[i][2]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_v1_reso_ptSetB_centSetA[i][3]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_v1_reso_ptSetB_centSetA[i][4]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_v1_reso_ptSetB_centSetA[i][5]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));

    mHist_v2_raw_ptSetB_centSetA[i][0]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_v2_raw_ptSetB_centSetA[i][1]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_v2_raw_ptSetB_centSetA[i][2]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_v2_raw_ptSetB_centSetA[i][3]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_v2_raw_ptSetB_centSetA[i][4]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_v2_raw_ptSetB_centSetA[i][5]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));

    mHist_v2_reso_ptSetB_centSetA[i][0]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[1]));
    mHist_v2_reso_ptSetB_centSetA[i][1]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[1],centSetA[2]));
    mHist_v2_reso_ptSetB_centSetA[i][2]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[3]));
    mHist_v2_reso_ptSetB_centSetA[i][3]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[2],centSetA[4]));
    mHist_v2_reso_ptSetB_centSetA[i][4]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[3]));
    mHist_v2_reso_ptSetB_centSetA[i][5]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetB[i],ptSetB[i+1],centSetA[0],centSetA[4]));
    // rap SetA, cent SetA
    mHist_SE_InvM_rapSetA_centSetA[i][0]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_rapSetA_centSetA[i][1]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_rapSetA_centSetA[i][2]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_rapSetA_centSetA[i][3]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_rapSetA_centSetA[i][4]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_rapSetA_centSetA[i][5]->SetTitle(Form("SE, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));

    mHist_v1_raw_rapSetA_centSetA[i][0]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_v1_raw_rapSetA_centSetA[i][1]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_v1_raw_rapSetA_centSetA[i][2]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_v1_raw_rapSetA_centSetA[i][3]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_v1_raw_rapSetA_centSetA[i][4]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_v1_raw_rapSetA_centSetA[i][5]->SetTitle(Form("v_{1}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));

    mHist_v1_reso_rapSetA_centSetA[i][0]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_v1_reso_rapSetA_centSetA[i][1]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_v1_reso_rapSetA_centSetA[i][2]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_v1_reso_rapSetA_centSetA[i][3]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_v1_reso_rapSetA_centSetA[i][4]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_v1_reso_rapSetA_centSetA[i][5]->SetTitle(Form("v_{1}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));

    mHist_v2_raw_rapSetA_centSetA[i][0]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_v2_raw_rapSetA_centSetA[i][1]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_v2_raw_rapSetA_centSetA[i][2]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_v2_raw_rapSetA_centSetA[i][3]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_v2_raw_rapSetA_centSetA[i][4]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_v2_raw_rapSetA_centSetA[i][5]->SetTitle(Form("v_{2}^{raw}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));

    mHist_v2_reso_rapSetA_centSetA[i][0]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[1]));
    mHist_v2_reso_rapSetA_centSetA[i][1]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[1],centSetA[2]));
    mHist_v2_reso_rapSetA_centSetA[i][2]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[3]));
    mHist_v2_reso_rapSetA_centSetA[i][3]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[2],centSetA[4]));
    mHist_v2_reso_rapSetA_centSetA[i][4]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[3]));
    mHist_v2_reso_rapSetA_centSetA[i][5]->SetTitle(Form("v_{2}^{resolution}, %3.1f<y<%3.1f, %3.f -%3.f%%",rapSetA[i],rapSetA[i+1],centSetA[0],centSetA[4]));
  }
  for(int pt=0; pt<2; pt++)
  {// pt SetA, cent SetA
    mHist_SE_InvM_ptSetA_centSetA[pt][0]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_SE_InvM_ptSetA_centSetA[pt][1]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_SE_InvM_ptSetA_centSetA[pt][2]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_SE_InvM_ptSetA_centSetA[pt][3]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_SE_InvM_ptSetA_centSetA[pt][4]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_SE_InvM_ptSetA_centSetA[pt][5]->SetTitle(Form("SE, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));

    mHist_v1_raw_ptSetA_centSetA[pt][0]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_v1_raw_ptSetA_centSetA[pt][1]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_v1_raw_ptSetA_centSetA[pt][2]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_v1_raw_ptSetA_centSetA[pt][3]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_v1_raw_ptSetA_centSetA[pt][4]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_v1_raw_ptSetA_centSetA[pt][5]->SetTitle(Form("v_{1}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));

    mHist_v1_reso_ptSetA_centSetA[pt][0]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_v1_reso_ptSetA_centSetA[pt][1]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_v1_reso_ptSetA_centSetA[pt][2]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_v1_reso_ptSetA_centSetA[pt][3]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_v1_reso_ptSetA_centSetA[pt][4]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_v1_reso_ptSetA_centSetA[pt][5]->SetTitle(Form("v_{1}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));

    mHist_v2_raw_ptSetA_centSetA[pt][0]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_v2_raw_ptSetA_centSetA[pt][1]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_v2_raw_ptSetA_centSetA[pt][2]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_v2_raw_ptSetA_centSetA[pt][3]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_v2_raw_ptSetA_centSetA[pt][4]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_v2_raw_ptSetA_centSetA[pt][5]->SetTitle(Form("v_{2}^{raw}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));

    mHist_v2_reso_ptSetA_centSetA[pt][0]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[1]));
    mHist_v2_reso_ptSetA_centSetA[pt][1]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[1],centSetA[2]));
    mHist_v2_reso_ptSetA_centSetA[pt][2]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[3]));
    mHist_v2_reso_ptSetA_centSetA[pt][3]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[2],centSetA[4]));
    mHist_v2_reso_ptSetA_centSetA[pt][4]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[3]));
    mHist_v2_reso_ptSetA_centSetA[pt][5]->SetTitle(Form("v_{2}^{resolution}, %3.1f<pt<%3.1f, %3.f -%3.f%%",ptSetA[pt],ptSetA[pt+1],centSetA[0],centSetA[4]));
  }
  // pt SetA, cent SetA
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mProfile_v1_raw_ptSetA_centSetA[pt][cent]  = mHist_v1_raw_ptSetA_centSetA[pt][cent]->ProfileX();
      mProfile_v1_reso_ptSetA_centSetA[pt][cent] = mHist_v1_reso_ptSetA_centSetA[pt][cent]->ProfileX();;
      mProfile_v2_raw_ptSetA_centSetA[pt][cent]  = mHist_v2_raw_ptSetA_centSetA[pt][cent]->ProfileX();;
      mProfile_v2_reso_ptSetA_centSetA[pt][cent] = mHist_v2_reso_ptSetA_centSetA[pt][cent]->ProfileX();;
    }
  }
  // pt SetA, cent SetB
  for(int pt=0; pt<2; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mProfile_v1_raw_ptSetA_centSetB[pt][cent]  = mHist_v1_raw_ptSetA_centSetB[pt][cent]->ProfileX();
      mProfile_v1_reso_ptSetA_centSetB[pt][cent] = mHist_v1_reso_ptSetA_centSetB[pt][cent]->ProfileX();;
      mProfile_v2_raw_ptSetA_centSetB[pt][cent]  = mHist_v2_raw_ptSetA_centSetB[pt][cent]->ProfileX();;
      mProfile_v2_reso_ptSetA_centSetB[pt][cent] = mHist_v2_reso_ptSetA_centSetB[pt][cent]->ProfileX();;
    }
  }
  // pt SetB, cent SetA
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<6;cent++){
      mProfile_v1_raw_ptSetB_centSetA[pt][cent]  = mHist_v1_raw_ptSetB_centSetA[pt][cent]->ProfileX();
      mProfile_v1_reso_ptSetB_centSetA[pt][cent] = mHist_v1_reso_ptSetB_centSetA[pt][cent]->ProfileX();;
      mProfile_v2_raw_ptSetB_centSetA[pt][cent]  = mHist_v2_raw_ptSetB_centSetA[pt][cent]->ProfileX();;
      mProfile_v2_reso_ptSetB_centSetA[pt][cent] = mHist_v2_reso_ptSetB_centSetA[pt][cent]->ProfileX();;
    }
  }
  // pt SetB, cent SetB
  for(int pt=0; pt<4; pt++)
  {
    for(int cent=0; cent<9;cent++){
      mProfile_v1_raw_ptSetB_centSetB[pt][cent]  = mHist_v1_raw_ptSetB_centSetB[pt][cent]->ProfileX();
      mProfile_v1_reso_ptSetB_centSetB[pt][cent] = mHist_v1_reso_ptSetB_centSetB[pt][cent]->ProfileX();;
      mProfile_v2_raw_ptSetB_centSetB[pt][cent]  = mHist_v2_raw_ptSetB_centSetB[pt][cent]->ProfileX();;
      mProfile_v2_reso_ptSetB_centSetB[pt][cent] = mHist_v2_reso_ptSetB_centSetB[pt][cent]->ProfileX();;
    }
  }
  // pt SetC, cent 0-60%, 0-80%
  for(int pt=0; pt<10; pt++)
  {
    for(int cent=0; cent<2;cent++){
      mProfile_v1_raw_ptSetC_centAll[pt][cent]  = mHist_v1_raw_ptSetC_centAll[pt][cent]->ProfileX();
      mProfile_v1_reso_ptSetC_centAll[pt][cent] = mHist_v1_reso_ptSetC_centAll[pt][cent]->ProfileX();;
      mProfile_v2_raw_ptSetC_centAll[pt][cent]  = mHist_v2_raw_ptSetC_centAll[pt][cent]->ProfileX();;
      mProfile_v2_reso_ptSetC_centAll[pt][cent] = mHist_v2_reso_ptSetC_centAll[pt][cent]->ProfileX();;
    }
  }
  // rap SetA, cent SetA
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<6;cent++){
      mProfile_v1_raw_rapSetA_centSetA[rap][cent]  = mHist_v1_raw_rapSetA_centSetA[rap][cent]->ProfileX();
      mProfile_v1_reso_rapSetA_centSetA[rap][cent] = mHist_v1_reso_rapSetA_centSetA[rap][cent]->ProfileX();;
      mProfile_v2_raw_rapSetA_centSetA[rap][cent]  = mHist_v2_raw_rapSetA_centSetA[rap][cent]->ProfileX();;
      mProfile_v2_reso_rapSetA_centSetA[rap][cent] = mHist_v2_reso_rapSetA_centSetA[rap][cent]->ProfileX();;
    }
  }
  // rap SetA, cent SetB
  for(int rap=0; rap<4; rap++)
  {
    for(int cent=0; cent<9;cent++){
      mProfile_v1_raw_rapSetA_centSetB[rap][cent]  = mHist_v1_raw_rapSetA_centSetB[rap][cent]->ProfileX();
      mProfile_v1_reso_rapSetA_centSetB[rap][cent] = mHist_v1_reso_rapSetA_centSetB[rap][cent]->ProfileX();;
      mProfile_v2_raw_rapSetA_centSetB[rap][cent]  = mHist_v2_raw_rapSetA_centSetB[rap][cent]->ProfileX();;
      mProfile_v2_reso_rapSetA_centSetB[rap][cent] = mHist_v2_reso_rapSetA_centSetB[rap][cent]->ProfileX();;
    }
  }
  outputFile->cd();
  wt.Write();
  // wt_tpc.Write();
  v1WtaWt->Write();
  outputFile->Write();
  // for(int EventTypeId=0;EventTypeId<_nEventTypeBins;EventTypeId++){
  //   mPhiWeightOutput[EventTypeId]->Divide(mPhiAveraged[EventTypeId]);
  //   delete mPhiAveraged[EventTypeId];
  // }
  mCorrectionOutputFile->Write();
  PhiMesonAnaOutputFile->Write();
}

// =========================== Get Psi from Q vector =============================================
Double_t GetPsi(Double_t Qx, Double_t Qy, Int_t order){
  Double_t temp;
  if ((Qx==0.0) && (Qy==0.0)) temp=-999.0;
  else{
    temp =  TMath::ATan2(Qy,Qx)/((Double_t)order);
    Double_t AngleWrapAround = 2.0*TMath::Pi()/(Double_t)order;
    if (temp<0.0) temp+= AngleWrapAround;
    else if (temp>AngleWrapAround) temp -= AngleWrapAround;
  }
  return temp;
}
