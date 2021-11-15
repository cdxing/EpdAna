#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include <iostream>
#include "Math/MinimizerOptions.h"
#include "TLatex.h"

Double_t Breit_Wigner_Poly_FitFunc(Double_t* x_val, Double_t* par)
{
  Double_t x, y, amp, mean, gamma;
  Double_t pol0, pol1, pol2, pol3, pol4, pol5;
  amp   = par[0];
  mean  = par[1];
  gamma = par[2];
  pol0  = par[3];
  pol1  = par[4];
  pol2  = par[5];
  pol3  = par[6];
  pol4  = par[7];
  pol5  = par[8];
  x = x_val[0];
  y = (amp/(2.0*TMath::Pi()))*gamma/((x-mean)*(x-mean)+gamma*gamma/4.0) + pol0 + pol1*x + pol2*x*x + pol3*x*x*x + pol4*x*x*x*x + pol5*x*x*x*x*x;
  return y;
}

TLatex* plotTopLegend(char* label,Float_t x=-1,Float_t y=-1,Float_t size=0.06,Int_t color=1,Float_t angle=0.0, Int_t font = 42, Int_t NDC = 1)
{
    // coordinates in NDC!
    // plots the string label in position x and y in NDC coordinates
    // size is the text size
    // color is the text color

//    if(x<0||y<0)
//    {   // defaults
//      x=gPad->GetLeftMargin()*1.15;
//      y=(1-gPad->GetTopMargin())*1.04;
//    }
    TLatex* text=new TLatex(x,y,label);
    text->SetTextFont(font);
    text->SetTextSize(size);
    if(NDC == 1) text->SetNDC();
    text->SetTextColor(color);
    text->SetTextAngle(angle);
    text->Draw();
    return text;
}

void FitPhi(TString inputfile = "/project/projectdirs/star/xusun/OutPut/AuAu39GeV/Phi/file_39GeV_Phi_0080_etagap_00_511_sig.root")
{
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(50000);

  TFile *File_InPut = TFile::Open(inputfile.Data());
  TH2F *h_Mass2[10];
  for(Int_t i = 0; i < 10; i ++)
  {
    TString HistName = Form("Mass2_pt_%d",i);
    h_Mass2[i] = (TH2F*)File_InPut->Get(HistName.Data());
  }

  TH1F *h_Mass2_proj[10];
  for(Int_t i = 0; i < 10; i++)
  {
    TString HistName = Form("Mass2_pt_proj_%d",i);
    h_Mass2_proj[i] = (TH1F*)h_Mass2[i]->ProjectionY(HistName.Data());
    h_Mass2_proj[i]->SetTitle("");
    h_Mass2_proj[i]->SetStats(0);
    h_Mass2_proj[i]->GetXaxis()->SetTitle("Iniv Mass");
    h_Mass2_proj[i]->GetYaxis()->SetTitle("Counts");
    h_Mass2_proj[i]->GetXaxis()->CenterTitle();
    h_Mass2_proj[i]->GetYaxis()->CenterTitle();
    h_Mass2_proj[i]->GetXaxis()->SetNdivisions(505,'N');
    h_Mass2_proj[i]->GetYaxis()->SetNdivisions(505,'N');
  }

  TF1 *f_BWP[10];
  TF1 *f_bg[10];
  for(Int_t i = 0; i < 10; i++)
  {
    TString FuncName = Form("f_Mass2_%d",i);
    f_BWP[i] = new TF1(FuncName.Data(),Breit_Wigner_Poly_FitFunc,0.99,1.08,9);
    f_BWP[i]->SetParameter(0,1000);
    f_BWP[i]->SetParameter(1,1.0193);
    //f_BWP[i]->SetParLimits(1,1.01,1.03);
    f_BWP[i]->SetParameter(2,0.002);
    f_BWP[i]->SetParameter(3,2500.0);
    f_BWP[i]->SetParameter(4,1000.0/0.06);
    //f_BWP[i]->SetParameter(5, 3.54347e+06);
    //f_BWP[i]->SetParameter(6, 3.63963e+06);
    //f_BWP[i]->SetParameter(7, 1.00524e+06);
    //f_BWP[i]->SetParameter(8,-4.44208e+06);
//    f_BWP[i]->Draw("l");

    f_BWP[i]->SetRange(1.0,1.05);
    h_Mass2_proj[i]->Fit(f_BWP[i],"NR");
    FuncName = Form("f_BG_%d",i);
    f_bg[i] = new TF1(FuncName.Data(),Breit_Wigner_Poly_FitFunc,0.99,1.08,9);
    f_bg[i]->FixParameter(0,0);
    for(Int_t j = 1; j < 9; j++)
    {
      f_bg[i]->FixParameter(j,f_BWP[i]->GetParameter(j));
    }
    f_bg[i]->SetLineColor(2);
    f_bg[i]->SetLineStyle(2);
    f_bg[i]->SetLineWidth(2);
  }

  TH1F *h_sig[10];
  Float_t significance[10];
  for(Int_t i = 0; i < 10; i++)
  {
    TH1F *h_Clone = h_Mass2_proj[i]->Clone("h_Clone");
    h_sig[i] = h_Clone;
    h_sig[i]->Add(f_bg[i],-1.0); 

    Float_t mean = f_BWP[i]->GetParameter(1);
    Float_t sigma = f_BWP[i]->GetParameter(2);
    Int_t bin_start = h_Mass2_proj[i]->FindBin(mean-2.5*sigma);
    Int_t bin_stop = h_Mass2_proj[i]->FindBin(mean+2.5*sigma);

    Float_t sig_back = h_Mass2_proj[i]->Integral(bin_start,bin_stop);
    Float_t sig = h_sig[i]->Integral(bin_start,bin_stop);
    significance[i] = sig/TMath::Sqrt(sig_back);
  }

  TCanvas *c_Mass2[10];
  for(Int_t i = 0; i < 10; i++)
  {
    TString CanName = Form("c_Mass2_%d",i);
    c_Mass2[i] = new TCanvas(CanName.Data(),CanName.Data(),10,10,800,800);
    c_Mass2[i]->cd();
    c_Mass2[i]->cd()->SetLeftMargin(0.15);
    c_Mass2[i]->cd()->SetBottomMargin(0.15);
    c_Mass2[i]->cd()->SetTicks(1,1);
    c_Mass2[i]->cd()->SetGrid(0,0);
    h_Mass2_proj[i]->DrawCopy("h");
    f_bg[i]->Draw("l same");
    f_BWP[i]->Draw("l same");
    h_sig[i]->DrawCopy("h same");

    TString legend = Form("significance = %2.2f",significance[i]);
    plotTopLegend((char*)legend.Data(),0.45,0.5,0.05,1,0.0,42,1);
  }
}
