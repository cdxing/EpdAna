void PutLabel(TGraph* g,int ring, double x);

void DrawTileEtaRange(){
  gStyle->SetOptStat(0);

  double vzmax=200.0;  // change this to control how far out in Vz the graph goes.  Other things adjust automatically

  double edge[17] = {4.6, 9.0, 13.4, 17.8, 23.33, 28.86, 34.39, 39.92, 45.45, 50.98, 56.51, 62.05, 67.58, 73.11, 78.64, 84.17, 89.70};  // radii of ring boundaries
  double zEPD=375.0;  // epd distance to center of STAR

  TGraph* Boundary[16];    // this will be filled in, to show the region in Vz-eta space for the 16 rings.
  for (int i=0; i<16; i++){Boundary[i] = new TGraph();}

  double vzStep=1.0;
  // this goes "from left to right on the inner radius of the ring"...
  for (double vz=-vzmax; vz<vzmax; vz+=vzStep){
    TVector3 Primary(-2.0,0.0,vz); // FXT set up vz (0,-2.,200)
    for (int ring=1; ring<17; ring++){
      TVector3 End(edge[ring-1],0.0,zEPD);
      TVector3 line = End-Primary;
      Boundary[ring-1]->SetPoint(Boundary[ring-1]->GetN(),vz,line.Eta());
      cout << ring << " " << vz << " " << line.Eta() << endl;
    }
  }
  // ...and now "from right to left on the outer radius of the ring"
  for (double vz=vzmax; vz>-vzmax; vz-=vzStep){
    TVector3 Primary(-2.0,0.0,vz);
    for (int ring=1; ring<17; ring++){
      TVector3 End(edge[ring],0.0,zEPD);
      TVector3 line = End-Primary;
      Boundary[ring-1]->SetPoint(Boundary[ring-1]->GetN(),vz,line.Eta());
    }
  }
  // and finally, close the boundary
  for (int ring=1; ring<17; ring++){
    double x,y;
    Boundary[ring-1]->GetPoint(0,x,y);
    Boundary[ring-1]->SetPoint(Boundary[ring-1]->GetN(),x,y);
  }

  TCanvas* tc = new TCanvas("RingBoundaries","Ring acceptances",800,800);
  TH2D* Frame = new TH2D("Bound","Ring acceptance versus vertex position",100,-vzmax*1.2,vzmax*1.2,100,1.,6.);
  Frame->GetXaxis()->SetTitle("V_{z} (cm)");
  Frame->GetYaxis()->SetTitle("#eta");
  Frame->DrawCopy();

  double fcolor[16] = {2,3,4,5,6,7,2,3,4,5,6,7,2,3,4,5};
  for (int ring=1; ring<17; ring++){
    Boundary[ring-1]->SetFillColor(fcolor[ring-1]);
    Boundary[ring-1]->Draw("same,f");
    PutLabel(Boundary[ring-1],ring,vzmax);
  }

  /*
  for (double eta=1.9; eta<5.5; eta+=0.2){
    TLine* tl = new TLine(-vzmax*1.2,eta,vzmax*1.2,eta);
    tl->Draw();
  }
  */

}
//=============================================================
void PutLabel(TGraph* g,int ring, double xedge){
  double *x,*y;
  x = g->GetX();
  y = g->GetY();
  double xlab = 1.05*xedge;
  double ylab;
  if (ring%2==1){    // odd number ring-- put label on left
    xlab *= -1;
    ylab = 0.5*(y[0]+y[g->GetN()-2]);
  }
  else{              // even number ring-- put label on right
    ylab = 0.5*(y[g->GetN()/2-1] + y[g->GetN()/2]);
  }
  TText* t = new TText(xlab,ylab,Form("%d",ring));
  t->SetTextSize(0.02);
  t->SetTextAlign(22);
  t->Draw();
}
