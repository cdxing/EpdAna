#include <iostream>
#include <fstream>
#include <cmath>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"

#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"

Double_t resoErr(Double_t corrAB, Double_t corrAC, Double_t corrBC,
                 Double_t errAB, Double_t errAC, Double_t errBC);
Double_t resoVal(Double_t corrAB, Double_t corrAC, Double_t corrBC);
const int _Ncentralities = 10;
int corrFinder(){
  TFile* inFile = new TFile("corrINPUT_26p5_2ndEp.root","READ");

  TProfile *profile_correlation_epd_east[6], *profile_correlation_epd_tpc[4];
  Double_t binContentEpd[6][_Ncentralities];
  Double_t binErrorEpd[6][_Ncentralities];
  Double_t binContentTpc[4][_Ncentralities];
  Double_t binErrorTpc[4][_Ncentralities];
  for(int i=0;i<6;i++){
    profile_correlation_epd_east[i] = (TProfile*)inFile->Get(Form("profile_correlation_epd_east%d",i));
  }
  for(int i=0;i<4;i++){
    profile_correlation_epd_tpc[i] = (TProfile*)inFile->Get(Form("profile_correlation_epd%d_tpc",i+1));
  }
  for(int i=0;i<6;i++){
    for(int j=0;j<_Ncentralities;j++){
      binContentEpd[i][j]=profile_correlation_epd_east[i]->GetBinContent(j+1);
      std::cout<<Form("profile_correlation_epd_east%d binNumber%d binContent: ",i,j+1) << binContentEpd[i][j];
      binErrorEpd[i][j]=profile_correlation_epd_east[i]->GetBinError(j+1);
      std::cout<<Form("; binError: ") << binErrorEpd[i][j]<<std::endl;
    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<_Ncentralities;j++){
      binContentTpc[i][j]=profile_correlation_epd_tpc[i]->GetBinContent(j+1);
      std::cout<<Form("profile_correlation_epd%d_tpc binNumber%d binContent: ",i+1,j+1) << binContentTpc[i][j];
      binErrorTpc[i][j]=profile_correlation_epd_tpc[i]->GetBinError(j+1);
      std::cout<<Form("; binError: ") << binErrorTpc[i][j]<<std::endl;
    }
  }

  // --------------------- sqrt bin Contents ---------------------------------
  for(int i=0;i<6;i++){
    for(int j=0;j<_Ncentralities;j++){
      if(binContentEpd[i][j]>=0){ // make sure bin content is not negative to sqrt
        binErrorEpd[i][j] *=1./(2.*sqrt(binContentEpd[i][j])); // Error proagation;
        std::cout<<Form("profile_correlation_epd_east%d binNumber%d sqrt binError: ",i,j+1) << binErrorEpd[i][j];
        binContentEpd[i][j]= sqrt(binContentEpd[i][j]);
        std::cout<<Form("; sqrt binContent: ") << binContentEpd[i][j]<<std::endl;
      } else{
        binErrorEpd[i][j] =0; // Error proagation;
        std::cout<<Form("profile_correlation_epd_east%d binNumber%d sqrt binError: ",i,j+1) << binErrorEpd[i][j];
        binContentEpd[i][j]=0;
        std::cout<<Form("; sqrt binContent: ") << binContentEpd[i][j]<<std::endl;
      }

    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<_Ncentralities;j++){
      if(binContentTpc[i][j]>0){ // make sure bin content is not negative to sqrt
        binErrorTpc[i][j] *=1./(2.*sqrt(binContentTpc[i][j])); // Error proagation;
        std::cout<<Form("profile_correlation_epd%d_tpc binNumber%d sqrt binError: ",i+1,j+1) << binErrorTpc[i][j];
        binContentTpc[i][j]=sqrt(binContentTpc[i][j]);
        std::cout<<Form("; sqrt binContent: ") << binContentTpc[i][j]<<std::endl;
      } else{
        binErrorTpc[i][j] =0; // Error proagation;
        std::cout<<Form("profile_correlation_epd%d_tpc binNumber%d sqrt binError: ",i+1,j+1) << binErrorTpc[i][j];
        binContentTpc[i][j]=0;
        std::cout<<Form("; sqrt binContent: ") << binContentTpc[i][j]<<std::endl;
      }

    }
  }

  // ------------------ Fill into historgrams ----------------------------------
  TFile* outFile = new TFile("corrResult_26p5_2ndEp.root","RECREATE");
  TH1D *hist_correlation_epd_east[6], *hist_correlation_epd_tpc[4];
  TCanvas* c1[5];
  TLegend *legend[5];
  for(int i=0;i<5;i++){
    c1[i] = new TCanvas(Form("c1_%d",i),Form("Event Plane Resolutions %d",i+1),200,10,700,500);
    c1[i]->GetFrame()->SetFillColor(21);
    c1[i]->GetFrame()->SetBorderSize(12);
    c1[i]->DrawFrame(0., 0.01, 100., 0.6);
    legend[i] = new TLegend(0.1,0.7,0.48,0.9);
  }
  const Int_t n = _Ncentralities;
  Double_t x[n]  = {5,15,25,35,45,55,65,75,85,95};
  Double_t ex[n] = {5,5,5,5,5,5,5,5,5,5};
  Double_t y1_1[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y1_2[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y1_3[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey1_1[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey1_2[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey1_3[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t y2_1[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y2_2[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y2_3[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey2_1[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey2_2[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey2_3[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t y3_1[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y3_2[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y3_3[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey3_1[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey3_2[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey3_3[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t y4_1[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y4_2[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y4_3[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey4_1[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey4_2[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey4_3[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t y5_1[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y5_2[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t y5_3[n]  = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey5_1[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey5_2[n] = {0,0,0,0,0,0,0,0,0,0};
  Double_t ey5_3[n] = {0,0,0,0,0,0,0,0,0,0};

  TGraphErrors *gr1[3],*gr2[3],*gr3[3],*gr4[3],*gr5[4];
  for(int i=0;i<6;i++){
    hist_correlation_epd_east[i] = new TH1D(Form("hist_correlation_epd_east%d",i),Form("hist_correlation_epd_east%d",i),_Ncentralities,0.5,_Ncentralities+.5);
  }
  for(int i=0;i<4;i++){
    hist_correlation_epd_tpc[i] = new TH1D(Form("hist_correlation_epd_tpc%d",i),Form("hist_correlation_epd_tpc%d",i),_Ncentralities,0.5,_Ncentralities+.5);
  }
  hist_correlation_epd_east[0]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-1} #minus #psi_{1}^{EPD-2})>}");
  hist_correlation_epd_east[1]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-1} #minus #psi_{1}^{EPD-3})>}");
  hist_correlation_epd_east[2]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-1} #minus #psi_{1}^{EPD-4})>}");
  hist_correlation_epd_east[3]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-2} #minus #psi_{1}^{EPD-3})>}");
  hist_correlation_epd_east[4]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-2} #minus #psi_{1}^{EPD-4})>}");
  hist_correlation_epd_east[5]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-3} #minus #psi_{1}^{EPD-4})>}");
  hist_correlation_epd_tpc[0]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-1} #minus #psi_{1}^{TPC})>}");
  hist_correlation_epd_tpc[1]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-2} #minus #psi_{1}^{TPC})>}");
  hist_correlation_epd_tpc[2]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-3} #minus #psi_{1}^{TPC})>}");
  hist_correlation_epd_tpc[3]->SetTitle("#sqrt{<cos(#psi_{1}^{EPD-4} #minus #psi_{1}^{TPC})>}");
  for(int i=0;i<6;i++){
    for(int j=0;j<_Ncentralities;j++){
      if(binErrorEpd[i][j]!=0){
        hist_correlation_epd_east[i]->SetBinContent(j+1,binContentEpd[i][j]);
        hist_correlation_epd_east[i]->SetBinError(j+1,binErrorEpd[i][j]);
      }
    }
  }
  for(int i=0;i<4;i++){
    for(int j=0;j<_Ncentralities;j++){
      if(binErrorTpc[i][j]!=0){
        hist_correlation_epd_tpc[i]->SetBinContent(j+1,binContentTpc[i][j]);
        hist_correlation_epd_tpc[i]->SetBinError(j+1,binErrorTpc[i][j]);
      }
    }
  }
  // EPD ep  resolutions
  //binContentEpd[0,1,2,3,4,5] corresponds to correlations of epd subs {1 2,1 3,1 4,2 3,2 4,3 4}
  //binContentTpc[0,1,2,3] corresponds to correlations of epd subs & tpc {1 t, 2 t, 3 t, 4 t}

  Double_t ep1Reso[3][n]; // binContentEpd[0,1,2][j] * binContentTpc[0][j] / binContentTpc[1,2,3][j]
  Double_t ep1Err[3][n]; // binErrorEpd[0,1,2][j] * binErrorTpc[0][j] / binErrorTpc[1,2,3][j]
  Double_t ep2Reso[3][n]; // binContentEpd[0,3,4][j] * binContentTpc[1][j] / binContentTpc[0,2,3][j]
  Double_t ep2Err[3][n]; // binErrorEpd[0,3,4][j] * binErrorTpc[1][j] / binErrorTpc[0,2,3][j]
  Double_t ep3Reso[3][n]; // binContentEpd[1,3,5][j] * binContentTpc[2][j] / binContentTpc[0,1,3][j]
  Double_t ep3Err[3][n]; // binErrorEpd[1,3,5][j] * binErrorTpc[2][j] / binErrorTpc[0,1,3][j]
  Double_t ep4Reso[3][n]; // binContentEpd[2,4,5][j] * binContentTpc[3][j] / binContentTpc[0,1,2][j]
  Double_t ep4Err[3][n]; // binErrorEpd[2,4,5][j] * binErrorTpc[3][j] / binErrorTpc[0,1,2][j]
  Double_t ep5Reso[3][n]; // binContentTpc[0,1,2][j] * binContentTpc[1,2,3][j] / binContentEpd[0,3,5][j]
  Double_t ep5Err[3][n]; // binErrorTpc[0,1,2][j] * binErrorTpc[1,2,3][j] / binErrorEpd[0,3,5][j]
  for(int j=0;j<_Ncentralities;j++){
    for(int i=0;i<3;i++){
      int order2_1 = 0; // EPD 2 reso
      int order2_2 = 0; // EPD 2 reso
      int order3_1 = 2 * i + 1; // 1, 3, 5 // EPD 3 reso
      int order3_2 = 0; // EPD 3 reso
      int order4_1 = 2; // EPD 4 reso
      int order4_2 = i; // EPD 4 reso // 0, 1, 2
      int order5_1 = i; // TPC reso // 0, 1, 2
      int order5_2 = i+1; // TPC reso // 1, 2, 3
      int order5_3 = 0; // TPC reso
      if(i>0){
        order2_1 = i+2; // 0, 3, 4
        order2_2 = i+1; // 0, 2, 3
        order3_2 = 2*i-1; // 0, 1, 3
        order4_1 += i+1; // 2, 4, 5
        order5_3 = 2*i+1; // 0, 3, 5
      }

      if(binErrorEpd[i][j]!=0&&binErrorTpc[0][j]!=0&&binErrorTpc[i+1][j]!=0){ // EPD 1 VS EPD {2,3,4} & TPC
        ep1Reso[i][j] = resoVal(binContentEpd[i][j],binContentTpc[0][j],binContentTpc[i+1][j]);
        ep1Err[i][j]= resoErr(binContentEpd[i][j],binContentTpc[0][j],binContentTpc[i+1][j],
                              binErrorEpd[i][j],binErrorTpc[0][j],binErrorTpc[i+1][j]);
      }
      if(binErrorEpd[order2_1][j]!=0&&binErrorTpc[1][j]!=0&&binErrorTpc[order2_2][j]!=0){ // EPD 2 VS EPD {1,3,4} & TPC
        ep2Reso[i][j] = resoVal(binContentEpd[order2_1][j],binContentTpc[1][j],binContentTpc[order2_2][j]);
        ep2Err[i][j]= resoErr(binContentEpd[order2_1][j],binContentTpc[1][j],binContentTpc[order2_2][j],
                              binErrorEpd[order2_1][j],binErrorTpc[1][j],binErrorTpc[order2_2][j]);
      }
      if(binErrorEpd[order3_1][j]!=0&&binErrorTpc[2][j]!=0&&binErrorTpc[order2_2][j]!=0){ // EPD 3 VS EPD {1,2,4} & TPC
        ep3Reso[i][j] = resoVal(binContentEpd[order3_1][j],binContentTpc[2][j],binContentTpc[order2_2][j]);
        ep3Err[i][j]= resoErr(binContentEpd[order3_1][j],binContentTpc[2][j],binContentTpc[order2_2][j],
                              binErrorEpd[order3_1][j],binErrorTpc[2][j],binErrorTpc[order2_2][j]);
      }
      if(binErrorEpd[order4_1][j]!=0&&binErrorTpc[3][j]!=0&&binErrorTpc[order4_2][j]!=0){ // EPD 4 VS EPD {1, 2, 3} & TPC
        ep4Reso[i][j] = resoVal(binContentEpd[order4_1][j],binContentTpc[3][j],binContentTpc[order4_2][j]);
        ep4Err[i][j]= resoErr(binContentEpd[order4_1][j],binContentTpc[3][j],binContentTpc[order4_2][j],
                              binErrorEpd[order4_1][j],binErrorTpc[3][j],binErrorTpc[order4_2][j]);
      }
      if(binErrorTpc[order5_1][j]!=0&&binErrorTpc[order5_2][j]!=0&&binErrorEpd[order5_3][j]!=0){ // TPC VS EPD {1, 2, 3} & EPD {2, 3, 4}
        ep5Reso[i][j] = resoVal(binContentTpc[order5_1][j],binContentTpc[order5_2][j],binContentEpd[order5_3][j]);
        ep5Err[i][j]= resoErr(binContentTpc[order5_1][j],binContentTpc[order5_2][j],binContentEpd[order5_3][j],
                              binErrorTpc[order5_1][j],binErrorTpc[order5_2][j],binErrorEpd[order5_3][j]);
      }
    }
  }
  for(int j=0;j<_Ncentralities;j++){
    y1_1[j]=ep1Reso[0][j];
    ey1_1[j]=ep1Err[0][j];
    y2_1[j]=ep2Reso[0][j];
    ey2_1[j]=ep2Err[0][j];
    y3_1[j]=ep3Reso[0][j];
    ey3_1[j]=ep3Err[0][j];
    y4_1[j]=ep4Reso[0][j];
    ey4_1[j]=ep4Err[0][j];
    y5_1[j]=ep5Reso[0][j];
    ey5_1[j]=ep5Err[0][j];

    y1_2[j]=ep1Reso[1][j];
    ey1_2[j]=ep1Err[1][j];
    y2_2[j]=ep2Reso[1][j];
    ey2_2[j]=ep2Err[1][j];
    y3_2[j]=ep3Reso[1][j];
    ey3_2[j]=ep3Err[1][j];
    y4_2[j]=ep4Reso[1][j];
    ey4_2[j]=ep4Err[1][j];
    y5_2[j]=ep5Reso[1][j];
    ey5_2[j]=ep5Err[1][j];

    y1_3[j]=ep1Reso[2][j];
    ey1_3[j]=ep1Err[2][j];
    y2_3[j]=ep2Reso[2][j];
    ey2_3[j]=ep2Err[2][j];
    y3_3[j]=ep3Reso[2][j];
    ey3_3[j]=ep3Err[2][j];
    y4_3[j]=ep4Reso[2][j];
    ey4_3[j]=ep4Err[2][j];
    y5_3[j]=ep5Reso[2][j];
    ey5_3[j]=ep5Err[2][j];
  }
  // ey3_1[5]=0;
  gr1[0] = new TGraphErrors(n,x,y1_1,ex,ey1_1);
  gr1[0]->SetTitle("EPD sub1 event plane ");
  gr1[0]->SetMarkerColor(1);
  gr1[0]->SetMarkerStyle(21);
  gr1[1] = new TGraphErrors(n,x,y1_2,ex,ey1_2);
  gr1[1]->SetTitle("EPD sub1 event plane ");
  gr1[1]->SetMarkerColor(2);
  gr1[1]->SetMarkerStyle(22);
  gr1[2] = new TGraphErrors(n,x,y1_3,ex,ey1_3);
  gr1[2]->SetTitle("EPD sub1 event plane ");
  gr1[2]->SetMarkerColor(3);
  gr1[2]->SetMarkerStyle(23);
  legend[0]->AddEntry(gr1[0],"EPD-1 VS EPD-2 & TPC","p");
  legend[0]->AddEntry(gr1[1],"EPD-1 VS EPD-3 & TPC","p");
  legend[0]->AddEntry(gr1[2],"EPD-1 VS EPD-4 & TPC","p");
  gr2[0] = new TGraphErrors(n,x,y2_1,ex,ey2_1);
  gr2[0]->SetTitle("EPD sub2 event plane ");
  gr2[0]->SetMarkerColor(1);
  gr2[0]->SetMarkerStyle(21);
  gr2[1] = new TGraphErrors(n,x,y2_2,ex,ey2_2);
  gr2[1]->SetTitle("EPD sub2 event plane ");
  gr2[1]->SetMarkerColor(2);
  gr2[1]->SetMarkerStyle(22);
  gr2[2] = new TGraphErrors(n,x,y2_3,ex,ey2_3);
  gr2[2]->SetTitle("EPD sub2 event plane ");
  gr2[2]->SetMarkerColor(3);
  gr2[2]->SetMarkerStyle(23);
  legend[1]->AddEntry(gr2[0],"EPD-2 VS EPD-1 & TPC","p");
  legend[1]->AddEntry(gr2[1],"EPD-2 VS EPD-3 & TPC","p");
  legend[1]->AddEntry(gr2[2],"EPD-2 VS EPD-4 & TPC","p");
  gr3[0] = new TGraphErrors(n,x,y3_1,ex,ey3_1);
  gr3[0]->SetTitle("EPD sub3 event plane ");
  gr3[0]->SetMarkerColor(1);
  gr3[0]->SetMarkerStyle(21);
  gr3[1] = new TGraphErrors(n,x,y3_2,ex,ey3_2);
  gr3[1]->SetTitle("EPD sub3 event plane ");
  gr3[1]->SetMarkerColor(2);
  gr3[1]->SetMarkerStyle(22);
  gr3[2] = new TGraphErrors(n,x,y3_3,ex,ey3_3);
  gr3[2]->SetTitle("EPD sub3 event plane ");
  gr3[2]->SetMarkerColor(3);
  gr3[2]->SetMarkerStyle(23);
  legend[2]->AddEntry(gr3[0],"EPD-3 VS EPD-1 & TPC","p");
  legend[2]->AddEntry(gr3[1],"EPD-3 VS EPD-2 & TPC","p");
  legend[2]->AddEntry(gr3[2],"EPD-3 VS EPD-4 & TPC","p");
  gr4[0] = new TGraphErrors(n,x,y4_1,ex,ey4_1);
  gr4[0]->SetTitle("EPD sub4 event plane ");
  gr4[0]->SetMarkerColor(1);
  gr4[0]->SetMarkerStyle(21);
  gr4[1] = new TGraphErrors(n,x,y4_2,ex,ey4_2);
  gr4[1]->SetTitle("EPD sub4 event plane ");
  gr4[1]->SetMarkerColor(2);
  gr4[1]->SetMarkerStyle(22);
  gr4[2] = new TGraphErrors(n,x,y4_3,ex,ey4_3);
  gr4[2]->SetTitle("EPD sub4 event plane ");
  gr4[2]->SetMarkerColor(3);
  gr4[2]->SetMarkerStyle(23);
  legend[3]->AddEntry(gr4[0],"EPD-4 VS EPD-1 & TPC","p");
  legend[3]->AddEntry(gr4[1],"EPD-4 VS EPD-2 & TPC","p");
  legend[3]->AddEntry(gr4[2],"EPD-4 VS EPD-3 & TPC","p");
  gr5[0] = new TGraphErrors(n,x,y5_1,ex,ey5_1);
  gr5[0]->SetTitle("TPC event plane ");
  gr5[0]->SetMarkerColor(1);
  gr5[0]->SetMarkerStyle(21);
  gr5[1] = new TGraphErrors(n,x,y5_2,ex,ey5_2);
  gr5[1]->SetTitle("TPC event plane ");
  gr5[1]->SetMarkerColor(2);
  gr5[1]->SetMarkerStyle(22);
  gr5[2] = new TGraphErrors(n,x,y5_3,ex,ey5_3);
  gr5[2]->SetTitle("TPC event plane ");
  gr5[2]->SetMarkerColor(3);
  gr5[2]->SetMarkerStyle(23);
  legend[4]->AddEntry(gr5[0],"TPC VS EPD-1 & EPD-2","p");
  legend[4]->AddEntry(gr5[1],"TPC VS EPD-2 & EPD-3","p");
  legend[4]->AddEntry(gr5[2],"TPC VS EPD-3 & EPD-4","p");
  c1[0]->cd();
  gr1[0]->Draw("P");
  gr1[1]->Draw("P");
  gr1[2]->Draw("P");
  legend[0]->Draw();
  c1[1]->cd();
  gr2[0]->Draw("P");
  gr2[1]->Draw("P");
  gr2[2]->Draw("P");
  legend[1]->Draw();
  c1[2]->cd();
  gr3[0]->Draw("P");
  gr3[1]->Draw("P");
  gr3[2]->Draw("P");
  legend[2]->Draw();
  c1[3]->cd();
  gr4[0]->Draw("P");
  gr4[1]->Draw("P");
  gr4[2]->Draw("P");
  legend[3]->Draw();
  c1[4]->cd();
  gr5[0]->Draw("P");
  gr5[1]->Draw("P");
  gr5[2]->Draw("P");
  legend[4]->Draw();
  for(int i=0;i<5;i++){
    c1[i]->Write();
  }
  outFile->Write();
  return 0;
}

Double_t resoErr(Double_t corrAB, Double_t corrAC, Double_t corrBC, Double_t errAB, Double_t errAC, Double_t errBC){
  Double_t error = -999.;
  if(corrBC!=0){
    error = corrAC / corrBC * errAB
                 + corrAB / corrBC * errAC
                 - corrAB * corrAC / (corrBC * corrBC) * errBC;
  }
  return error;
}
Double_t resoVal(Double_t corrAB, Double_t corrAC, Double_t corrBC){
  Double_t resolution = -999.;
  if(corrBC!=0){
    resolution = corrAB * corrAC / corrBC;
  }
  return resolution;
}
