#include "TH1F.h"
#include "TCanvas.h"
#include <iostream>

void FitKaonMass2()
{
  Float_t Mass2_kaon = 0.493677*0.493677;
  Float_t pt_bin[4] = {0.55,1.05,1.55,2.05};
  Float_t width[4] = {0.008365,0.01799,0.03305,0.05521};
  Float_t nSigma = 2.5;

  TH1F *h_low = new TH1F("h_low","h_low",300,0,3.0);
  TH1F *h_up = new TH1F("h_up","h_up",300,0,3.0);
  for(Int_t i = 0; i < 4; i++)
  {
    h_low->SetBinContent(h_low->FindBin(pt_bin[i]),(Mass2_kaon-nSigma*width[i]));
    h_up->SetBinContent(h_up->FindBin(pt_bin[i]),(Mass2_kaon+nSigma*width[i]));
  }
  h_low->SetMarkerStyle(20);
//  h_low->Fit("pol1");
  h_up->SetLineColor(2);
  h_up->SetMarkerStyle(21);
  h_up->Draw("p");
  h_low->Draw("p same");
}
