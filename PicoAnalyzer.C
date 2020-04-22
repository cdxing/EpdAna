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
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TRandom.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TMath.h"

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
const Int_t Ncentralities = 10;
// const Int_t order         = 20;
// const Int_t twoorder      = 2 * order;
Double_t GetPsi(Double_t Qx, Double_t Qy, Int_t order);

//////////////////////////////// Main Function /////////////////////////////////
void PicoAnalyzer(const Char_t *inFile = "/star/data01/pwg/dchen/Ana/fxtPicoAna/files/PicoDst/st_physics_16140033_raw_0000002.picoDst.root",
                      TString outFile = "test_EpdEP",
                      Int_t   inputp1 = 1
                    )
{

  Int_t EpOrder = inputp1; // Event plane Fourier expansion order = 1, 2, 3
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
  StEpdEpFinder *mEpFinder = new StEpdEpFinder(5,EpdEpOutputName,"/star/u/dchen/GitHub/EpdAna/EpdEpCorrectionHistograms_INPUT.root");
  StEpdGeom *mEpdGeom = new StEpdGeom();
  int format = 2;
  mEpFinder->SetEpdHitFormat(format);    // format=0/1/2 for StEpdHit/StMuEpdHit/StPicoEpdHit
  Double_t mThresh = 0.3; // EPD EP by hand
  Double_t mMax = 2.0; // EPD EP by hand
  Double_t etaRange[5] = {-5.16,-3.82,-3.28,-2.87,-2.60}; // EPD eta range to set 4 sub EPD EP
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
  mEpFinder->SetEtaWeights(1,wt);
  mEpFinder->SetnMipThreshold(mThresh);    // recommended by EPD group
  mEpFinder->SetMaxTileWeight(mMax);     // recommended by EPD group 3.0
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
  TH1D  *hist_runId = new TH1D("hist_runId","Event runId",20001,-0.5,20000.5);
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
  TH1D *hist_ratio_cut = new TH1D("hist_ratio_cut","hist_ratio_cut",100,0,2);
  TH1D *hist_nHits_cut = new TH1D("hist_nHits_cut","hist_nHits_cut",100,-0.5,99.5);
  TH1D *hist_ndEdx_cut = new TH1D("hist_ndEdx_cut","hist_ndEdx_cut",100,-0.5,99.5);
  TH1D *hist_DCA_cut = new TH1D("hist_DCA_cut","hist_DCA_cut",100,0,10.0);
  // ------------------ Centrality QA histograms ----------------------------------
  TH1D *hist_cent = new TH1D("hist_cent","Centrality",Ncentralities+1,-0.5,Ncentralities+0.5);
  TH1D *hist_realTrackMult = new TH1D("hist_realTrackMult","Actual track multiplicity",1001,-0.5,1000.5);
  TH2D *hist_realTrackMult_refmult = new TH2D("hist_realTrackMult_refmult","Actual track multiplicity vs. RefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  TH2D *hist_realTrackMult_grefmult = new TH2D("hist_realTrackMult_grefmult","Actual track multiplicity vs. gRefMult",1001,-0.5,1000.5,1001,-0.5,1000.5);
  // ------------------ EPD event plane histograms ----------------------------------
  TH2D *hist2_Epd_east_Qy_Qx_raw[5],*hist2_Epd_east_Qy_Qx_Weighted[5];
  TH1D *hist_Epd_east_psi_raw[5],*hist_Epd_east_psi_Weighted[5],*hist_Epd_east_psi_Shifted[5];
  TH2D *hist2_Epd_east_Qy_Qx_raw_ini[5];
  TH1D *hist_Epd_east_psi_raw_ini[5],*hist_Epd_east_psi_Shifted_ini[5];
  for(int sub=0; sub<5; sub++){
    hist2_Epd_east_Qy_Qx_raw[sub]= new TH2D(Form("hist2_Epd_east_Qy_Qx_raw_%d",sub),Form("EPD east Qy vs Qx sub%d",sub),600,-3.0,3.0,600,-3.0,3.0);
    hist2_Epd_east_Qy_Qx_Weighted[sub]= new TH2D(Form("hist2_Epd_east_Qy_Qx_Weighted_%d",sub),Form("EPD east Qy vs Qx (Weighted) sub%d",sub),600,-3.0,3.0,600,-3.0,3.0);
    hist_Epd_east_psi_raw[sub] = new TH1D(Form("hist_Epd_east_psi_raw_%d",sub),Form("EPD east EP sub%d",sub),1024,-1.0,7.0);
    hist_Epd_east_psi_Weighted[sub] = new TH1D(Form("hist_Epd_east_psi_Weighted_%d",sub),Form("EPD east EP (Weighted) sub%d",sub),1024,-1.0,7.0);
    hist_Epd_east_psi_Shifted[sub] = new TH1D(Form("hist_Epd_east_psi_Shifted_%d",sub),Form("EPD east EP (Weighted & Shifted) sub%d",sub),1024,-1.0,7.0);

    hist2_Epd_east_Qy_Qx_raw_ini[sub]= new TH2D(Form("hist2_Epd_east_Qy_Qx_raw_ini_%d",sub),Form("EPD east Qy vs Qx sub%d",sub),600,-3.0,3.0,600,-3.0,3.0);
    hist_Epd_east_psi_raw_ini[sub] = new TH1D(Form("hist_Epd_east_psi_raw_ini_%d",sub),Form("EPD east EP sub%d",sub),1024,-1.0,7.0);
    hist_Epd_east_psi_Shifted_ini[sub] = new TH1D(Form("hist_Epd_east_psi_Shifted_ini_%d",sub),Form("EPD east EP (Weighted & Shifted) sub%d",sub),1024,-1.0,7.0);
  }
  // ------------------ EPD event plane ab intio QA histograms ----------------------------------
  TH1D *hist_Epdeta = new TH1D("hist_Epdeta","epd eta",700,-6.5,0.5);
  TH1D *hist_nMip = new TH1D("hist_nMip","nMIP of tile: 0:1:1 ",64,-0.5,9.5);
  // (3) =========================== Event loop ====================================
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++)
  {
    // ---------------------- Event reading quality assurance ----------------------
    if((iEvent+1)%100 == 0) std::cout << "Working on event #[" << (iEvent+1)
    << "/" << events2read << "]" << std::endl;
    Bool_t readEvent = picoReader->readPicoEvent(iEvent);
    if( !readEvent ) {
        std::cout << "Something went wrong, Master! Nothing to analyze..."
        << std::endl;
        break;
    }
    StPicoDst     *dst = picoReader->picoDst();
    StPicoEvent *event = dst->event();
    if( !event ) {
        std::cout << "Something went wrong, Master! Event is hiding from me..."
        << std::endl;
        break;
    }
    mEvtcut[0]++;// No event cut yet

    // (4) =================== Get event parameters ================================
    Int_t runId       = event->runId();
    Int_t nTracks     = dst->numberOfTracks();

    const Float_t   B = event->bField(); // Magnetic field
    double d_MagField = event->bField();
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
        if(triggerIDs[i] == 630802) b_bad_trig = false; // hlt_fixedTargetGood 7.2GeV
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
    bool b_bad_zvtx   =  ((d_zvtx < 199.0) || (d_zvtx > 202.0)); //FXT_26p5_2018
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

    // (5) =============== Track loop to determine good tracks =================
    int nGoodTracks = 0;
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++){
      StPicoTrack *picoTrack = dst->track(iTrk);
      mTrkcut[0]++; // 0. No track cut
      if(!picoTrack) continue;
      mTrkcut[1]++; // 1. pico track cut
      StPicoBTofPidTraits *trait = NULL;
      // ----------------------- Physics values of tracks --------------------------
      double        d_tofBeta    = -999.;
      if(picoTrack->isTofTrack()) trait = dst->btofPidTraits( picoTrack->bTofPidTraitsIndex() );
      if(trait)        d_tofBeta = trait->btofBeta();
      double d_px  = picoTrack->gMom().x();
      double d_py  = picoTrack->gMom().y();
      double d_pz  = picoTrack->gMom().z();
      double d_pT  = picoTrack->gPt();
      double d_mom = sqrt(d_pT*d_pT + d_pz*d_pz);
      double mass2 = d_mom*d_mom*((1.0/(d_tofBeta*d_tofBeta))-1.0);
      Double_t eta = picoTrack->pMom().Eta();
      // --------------- QA plots before major track cuts ----------------------
      hist_px_py->Fill(d_px,d_py);
      hist_pz   ->Fill(d_pz);
      hist_pt   ->Fill(d_pT);
      hist_mom  ->Fill(d_mom);
      hist_mass2->Fill(mass2);
      hist_eta  ->Fill(eta);
      hist_ratio->Fill(((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()));
      hist_nHits->Fill((double)picoTrack->nHitsFit());
      hist_ndEdx->Fill(picoTrack->nHitsDedx());
      hist_DCA  ->Fill(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z));

      if(!picoTrack->isPrimary()) continue;
      mTrkcut[2]++; // 2. Primary track cut
      bool    b_bad_dEdx     = (picoTrack->nHitsDedx() <= 0);
      bool    b_bad_tracking = (((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()) < 0.52);
      bool b_not_enough_hits = ((double)picoTrack->nHitsFit()) < 15;
      bool    b_bad_DCA      = (picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z) >= 3.0);
      bool    b_bad_track    = b_bad_dEdx || b_bad_tracking || b_not_enough_hits || b_bad_DCA;
      if(b_bad_track) continue;
      mTrkcut[3]++; // 3. Bad track cuts
      nGoodTracks++; // nGoodTracks is used to determine centrality later in the event loop
      // --------------- QA plots after major track cuts ----------------------
      hist_px_py_cut->Fill(d_px,d_py);
      hist_pz_cut   ->Fill(d_pz);
      hist_pt_cut   ->Fill(d_pT);
      hist_mom_cut  ->Fill(d_mom);
      hist_mass2_cut->Fill(mass2);
      hist_eta_cut  ->Fill(eta);
      hist_ratio_cut->Fill(((double)picoTrack->nHitsFit() / (double)picoTrack->nHitsPoss()));
      hist_nHits_cut->Fill((double)picoTrack->nHitsFit());
      hist_ndEdx_cut->Fill(picoTrack->nHitsDedx());
      hist_DCA_cut  ->Fill(picoTrack->gDCA(primaryVertex_X,primaryVertex_Y,primaryVertex_Z));
      if(d_tofBeta == -999) continue;
      mTrkcut[4]++; // 4. Bad tof track cut, to see how many tracks with tof information
    } // Track loop to determine good tracks

    // (6) ================ Centrality definition ===============================
    Int_t centrality = 0;
    bool a_b_cent[10]={false};
    bool b_pileup   = (nGoodTracks > 270);
    bool b_low_mult = (nGoodTracks < 10);
    a_b_cent[0]     = (nGoodTracks >= 200); // 0 - 10%
    a_b_cent[1]     = (nGoodTracks >= 150 && nGoodTracks < 200); // 10 - 20%
    a_b_cent[2]     = (nGoodTracks >= 124 && nGoodTracks < 150); // 20 - 30%
    a_b_cent[3]     = (nGoodTracks >= 100 && nGoodTracks < 124); // 30 - 40%
    a_b_cent[4]     = (nGoodTracks >= 72  && nGoodTracks < 100); // 40 - 50%
    a_b_cent[5]     = (nGoodTracks >= 50  && nGoodTracks < 72); // 50 - 60%
    a_b_cent[6]     = (nGoodTracks >= 40  && nGoodTracks < 50); // 60 - 70%
    a_b_cent[7]     = (nGoodTracks >= 30  && nGoodTracks < 40); // 70 - 80%
    a_b_cent[8]     = (nGoodTracks >= 20  && nGoodTracks < 30); // 80 - 90%
    a_b_cent[9]     = (nGoodTracks >= 10  && nGoodTracks < 20); // >90%
    for(int i=0;i<10;i++){
      if(a_b_cent[i]) centrality = i+1;
    }
    hist_cent->Fill(centrality);
    hist_realTrackMult->Fill(nGoodTracks);
    hist_realTrackMult_refmult->Fill(nGoodTracks,refMult);
    hist_realTrackMult_grefmult->Fill(nGoodTracks,grefMult);

    // (7) ================ EPD event plane ====================================
    // (7.1) ------------- EPD ep from Mike Lisa's class StEpdEpFinder -----------------
    StEpdEpInfo mResult[5];
    Double_t EastRawQx[5],EastRawQy[5],EastWeightedQx[5],EastWeightedQy[5];
    for(int i=0;i<5;i++){
      mResult[i] = (StEpdEpInfo)mEpFinder->Results(mEpdHits,pVtx,i);  // and now you have all the EP info you could ever want :-)
      EastRawQx[i] = (Double_t) mResult[i].EastRawQ(EpOrder).X();
      // std::cout<<Form("east raw qx sub %d =",i)<<EastRawQx[i]<<std::endl;
      EastRawQy[i] = (Double_t) mResult[i].EastRawQ(EpOrder).Y();
      EastWeightedQx[i] = (Double_t) mResult[i].EastPhiWeightedQ(EpOrder).X();
      EastWeightedQy[i] = (Double_t) mResult[i].EastPhiWeightedQ(EpOrder).Y();
      if(EastRawQx[i]!=0 || EastRawQy[i]!=0 )
      {
        hist2_Epd_east_Qy_Qx_raw[i]->Fill(EastRawQx[i],EastRawQy[i]);
        hist2_Epd_east_Qy_Qx_Weighted[i]->Fill(EastWeightedQx[i],EastWeightedQy[i]);
        hist_Epd_east_psi_raw[i]->Fill(mResult[i].EastRawPsi(EpOrder));
        hist_Epd_east_psi_Weighted[i]->Fill(mResult[i].EastPhiWeightedPsi(EpOrder));
        hist_Epd_east_psi_Shifted[i]->Fill(mResult[i].EastPhiWeightedAndShiftedPsi(EpOrder));
      }
    }
    // (7.2) ------------------- EPD EP by hand ---------------------------------
    // refer to Mike's StEpdEpFinder and Yang's BBC Ep
    Int_t N_Epd_east[5]={0}; //Count # of hits in each eta region /// indices: [etaBin]
    Double_t QrawEastSide[5][2]={0.0};       /// indices: [etaBin][x,y]
    Double_t PsiEastRaw[5]={-999.0};           /// indices: [etaBin]
    Double_t PsiEastShifted[5]={-999.0};       /// indices: [etaBin]
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
      hist_Epdeta->Fill(eta);
      //--------------------------------
      // now calculate Q-vectors
      //--------------------------------
      for(int EventTypeId=0;EventTypeId<5;EventTypeId++){
        int etaBin = (int)wt.GetXaxis()->FindBin(fabs(eta));
        double etaWeight = (double)wt.GetBinContent(etaBin,EventTypeId+1);
        if(etaWeight>0.0) N_Epd_east[EventTypeId]++;
        double Cosine = cos(phi*(double)EpOrder);
        double Sine   = sin(phi*(double)EpOrder);
        QrawEastSide[EventTypeId][0] += etaWeight * TileWeight * Cosine;
        QrawEastSide[EventTypeId][1] += etaWeight * TileWeight * Sine;
      }
    } // loop over EPD hits
    // Before going any farther, flip the sign of the 1st-order Q-vector on the East side.
    //  I want the rapidity-odd first-order event plane.
    for(int EventTypeId=0;EventTypeId<5;EventTypeId++){
      for (int xy=0; xy<2; xy++){
        QrawEastSide[EventTypeId][xy]           *= -1.0;
      }
    }
    //---------------------------------
    // Calculate unshifted EP angles
    //---------------------------------
    for(int EventTypeId=0;EventTypeId<5;EventTypeId++){
      if(N_Epd_east[EventTypeId]<5) continue;
      if(QrawEastSide[EventTypeId][0] || QrawEastSide[EventTypeId][1] )PsiEastRaw[EventTypeId] = GetPsi(QrawEastSide[EventTypeId][0],QrawEastSide[EventTypeId][1],EpOrder);
    }

    for(int EventTypeId=0;EventTypeId<5;EventTypeId++){
      if(QrawEastSide[EventTypeId][0] || QrawEastSide[EventTypeId][1] )
      {
        hist2_Epd_east_Qy_Qx_raw_ini[EventTypeId]->Fill(EastRawQx[EventTypeId],EastRawQy[EventTypeId]);
        if(PsiEastRaw[EventTypeId]!=-999.0) hist_Epd_east_psi_raw_ini[EventTypeId]->Fill(PsiEastRaw[EventTypeId]);
        if(PsiEastRaw[EventTypeId]==0.0){
          std::cout<<"Psi = " <<PsiEastRaw[EventTypeId]<<endl;
          std::cout<<"Qx = " <<QrawEastSide[EventTypeId][0]<<endl;
          std::cout<<"Qy = " <<QrawEastSide[EventTypeId][1]<<endl;
        }
        // hist_Epd_east_psi_Shifted_ini[EventTypeId]->Fill(PsiEastRaw[EventTypeId]);
      }

    }

  }  // Event Loop
  // --------------------- Set histograms axises titles --------------------------------
  hist_runId->GetXaxis()->SetTitle("RunId");
  hist_runId->GetYaxis()->SetTitle("# of events");
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
  hist_realTrackMult->GetXaxis()->SetTitle("# of events");
  hist_realTrackMult_refmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_refmult->GetXaxis()->SetTitle("RefMult");
  hist_realTrackMult_grefmult->GetXaxis()->SetTitle("TrackMult");
  hist_realTrackMult_grefmult->GetXaxis()->SetTitle("gRefMult");
  for(int sub=0; sub<5; sub++){
    hist2_Epd_east_Qy_Qx_raw[sub]->GetXaxis()->SetTitle("Q_x^{EPD east}_{1} ");
    hist2_Epd_east_Qy_Qx_raw[sub]->GetYaxis()->SetTitle("Q_y^{EPD east}_{1} ");
    hist2_Epd_east_Qy_Qx_Weighted[sub]->GetXaxis()->SetTitle("Q_x^{EPD east}_{1} ");
    hist2_Epd_east_Qy_Qx_Weighted[sub]->GetYaxis()->SetTitle("Q_y^{EPD east}_{1} ");
    hist_Epd_east_psi_raw[sub]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_raw[sub]->GetYaxis()->SetTitle("# of events");
    hist_Epd_east_psi_Weighted[sub]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_Weighted[sub]->GetYaxis()->SetTitle("# of events");
    hist_Epd_east_psi_Shifted[sub]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_Shifted[sub]->GetYaxis()->SetTitle("# of events");
    hist2_Epd_east_Qy_Qx_raw_ini[sub]->GetXaxis()->SetTitle("Q_x^{EPD east}_{1} ");
    hist2_Epd_east_Qy_Qx_raw_ini[sub]->GetYaxis()->SetTitle("Q_y^{EPD east}_{1} ");
    hist_Epd_east_psi_raw_ini[sub]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_raw_ini[sub]->GetYaxis()->SetTitle("# of events");
    hist_Epd_east_psi_Shifted_ini[sub]->GetXaxis()->SetTitle("#psi^{EPD east}_{1} [Radian]");
    hist_Epd_east_psi_Shifted_ini[sub]->GetYaxis()->SetTitle("# of events");
  }
  hist_Epdeta->GetXaxis()->SetTitle("#eta");
  hist_Epdeta->GetYaxis()->SetTitle("# of hits");
  hist_nMip->GetXaxis()->SetTitle("nMIP");
  hist_nMip->GetYaxis()->SetTitle("# of hits");
  wt.Write();
  outputFile->Write();
  mEpFinder->Finish();
}

// =========================== Get Psi from Q vector =============================================
Double_t GetPsi(Double_t Qx, Double_t Qy, Int_t order){
  Double_t temp;
  if ((Qx==0.0)||(Qy==0.0)) temp=-999.0;
  else{
    temp = atan2(Qy,Qx)/((double)order);
    Double_t AngleWrapAround = 2.0*TMath::Pi()/(Double_t)order;
    if (temp<0.0) temp+= AngleWrapAround;
    else if (temp>AngleWrapAround) temp -= AngleWrapAround;
  }
  return temp;
}
