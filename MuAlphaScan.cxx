#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void MuAlphaScan(){

  double Lumi_scalefactor = 1.0;
  int Irebin = 1;

  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-- rebin histograms if necessary
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);


  TH2D *h_MuAlphaScan = new TH2D("h_MuAlphaScan", "", 1000, 0.3, 3.0, 1000, 0.5, 2.5);
  TH1D *h_MuRange = new TH1D("h_MuRange", "", 1000, 0.3, 3.0);
  TH1D *h_AlphaRange = new TH1D("h_AlphaRange", "", 1000, 0.5, 2.5);

  for (int i_bin = 1; i_bin<=h_MuAlphaScan->GetNbinsX(); i_bin++ ){
    double Mu = h_MuRange->GetBinCenter(i_bin);
    for (int j_bin = 1; j_bin<=h_MuAlphaScan->GetNbinsY(); j_bin++ ){
      double Alpha = h_AlphaRange->GetBinCenter(j_bin);
      double loglik = 1e-10;
      for (int l_bin = 1; l_bin <= h_bgr->GetNbinsX(); l_bin++){
        double si = h_sig->GetBinContent(l_bin);
        double bi = h_bgr->GetBinContent(l_bin);
        double ni = h_data->GetBinContent(l_bin);
        loglik += TMath::Log(TMath::Poisson(ni,(Mu*si + Alpha*bi))); 
      } // end loop over bins

      h_MuAlphaScan->SetBinContent(i_bin, j_bin, -2.*loglik);
    }
  }

  double min_val, max_val;
  int locminx, locminy, locminz;
  double mu_opti, alpha_opti;

  h_MuAlphaScan->GetMinimumBin(locminx, locminy, locminz);
  mu_opti = h_MuRange->GetBinCenter(locminx);
  alpha_opti = h_AlphaRange->GetBinCenter(locminy);
  min_val = h_MuAlphaScan->GetBinContent(locminx, locminy);

  cout << mu_opti << " " << alpha_opti << endl;
  cout << min_val << endl;
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetGrid();

  gStyle->SetOptStat(0);
  
  h_MuAlphaScan->Draw("COLZ"); 


  TH2D *h_MuAlphaScan_cont = (TH2D*) h_MuAlphaScan->Clone();

  double contours[2];

  contours[0] = min_val+1;
  contours[1] = min_val+2;

  h_MuAlphaScan_cont->SetContour(2,contours);

  h_MuAlphaScan_cont->SetLineColorAlpha(1, 0.00);;
  h_MuAlphaScan_cont->SetLineWidth(2);;
  h_MuAlphaScan_cont->Draw("same CONT3"); 

  double x[1], y[1];
  x[0] = mu_opti;
  y[0] = alpha_opti;
  TGraph *h_Min = new TGraph(1,x,y);
  
  h_Min->SetMarkerStyle(34); 
  h_Min->SetMarkerSize(1.2); 
  h_Min->Draw("same P"); 
  
  AddText( 0.900, 0.035, "mu_s",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, "alpha_bgr" ,0.060,90.,"right"); // Y-axis
  AddText( 0.920, 0.950, "-2ln(L)" ,0.060,0.,"right"); // Z-axis

  return;
}
