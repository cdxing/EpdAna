#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TH1.h"
#include "TMath.h"

int centFinder(Int_t Ncentralities = 16,
               Int_t pileUp = 245,
               Int_t lowMult = 5){
  TFile* inFile = new TFile("./plots/centAllEvents.root","READ");
  TH1D *h1;
  h1=(TH1D*)inFile->Get("hist_realTrackMult");
  Int_t binCenter[Ncentralities];
  for(int i =0;i<Ncentralities;i++){
    binCenter[i]=0;
  }
  // Int_t binNumber[Ncentralities];
  Double_t sections[Ncentralities];
  Int_t Nentries = 0;
  Int_t Ncounts = 0;
  TLine* lines[Ncentralities];
  for(Int_t i=(Int_t)h1->FindBin(lowMult);i<=(Int_t)h1->FindBin(pileUp);i++){
    Nentries+=(Int_t)h1->GetBinContent(i);
  }
  for(Int_t cent =0; cent<Ncentralities;cent++){
    sections[cent] = (Double_t)(cent+1)*Nentries/Ncentralities;
  }

  for(Int_t i=(Int_t)h1->FindBin(lowMult);i<=(Int_t)h1->FindBin(pileUp);i++){
    Ncounts+=(Int_t)h1->GetBinContent(i);
    for(Int_t cent =Ncentralities-1; cent>=0;cent--){
      if(Ncounts>=sections[cent]){
        if(binCenter[cent]!=0)break;
        binCenter[cent]=i-1;//(Int_t)h1->GetBinCenter(i);
        std::cout<<Form("Ncounts%d Section: ",cent)<< Ncounts <<std::endl; // Array to string?

        break;
      }
    }

  }
  std::string str;
  std::string str1;

  for(Int_t i : binCenter){
    str += std::to_string(i);
    str += ",";
  }
  for(Int_t i : sections){
    str1 += std::to_string(i);
    str1 += ",";
  }
  std::cout<<"Bin number of pileUp: "<< (Int_t)h1->FindBin(pileUp) <<std::endl; // Array to string?
  std::cout<<"Total entries: "<< Nentries <<std::endl; // Array to string?
  std::cout<<"Entries Section: "<< str1 <<std::endl; // Array to string?
  std::cout<<"Centrality Section: "<< str <<std::endl; // Array to string?

  h1->Draw();
  for(Int_t cent = 0; cent<Ncentralities;cent++){
    lines[cent] = new TLine(binCenter[cent],0,binCenter[cent],h1->GetMaximum());
    lines[cent]->SetLineColor(cent+1);
    lines[cent]->Draw("same");
  }


  return 0;
}
