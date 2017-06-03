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

  //TAxis *xaxis = h->GetXaxis();
  //TAxis *yaxis = h->GetYaxis();
  //Int_t binx = xaxis->FindBin(x);
  //Int_t biny = yaxis->FindBin(y);

  cout << mu_opti << " " << alpha_opti << endl;
  cout << min_val << endl;

  //int i_mu_1sig, i_mu_2sig;
  //int dump1, dump2;

  //double mu_1sig = h_MuAlphaScan->GetBinWithContent2(min_val+1,i_mu_1sig,dump1,1, locminx, -1, 1, 0.2);

  //cout << mu_1sig << endl;
  //cout << i_mu_1sig << endl;

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



//  int Nbins = 200;
//    
//  TH1D *h_bgr_TS = new TH1D("h_bgr_TS","", Nbins, -40, 40);
//  TH1D *h_sigbgr_TS = new TH1D("h_sigbgr_TS","", Nbins, -40, 40);
//
//  for (int i_bin = 1; i_bin<=h_LumiScan_spb->GetNbinsX(); i_bin++ ){
//    double Lumi_scalefactor = h_LumiScan_spb->GetBinLowEdge(i_bin);
//  
//    //-- Get histograms from the files (higgs, zz and data)
//    TH1D *h_sig, *h_bgr, *h_data;
//    h_sig  = GetMassDistribution(125, Lumi_scalefactor);
//    h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
//    h_data = GetMassDistribution(2, Lumi_scalefactor);
//
////    /////// Scaled:
////    
////    double alpha_bgr_nom = 1.11;
////    double alpha_bgr_sig = 0.07;
////    
////    h_bgr->Scale(alpha_bgr_nom);
//
//    /////////////////
//    
//    TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
//    h_sigbgr->Add(h_bgr);
//
//    for (int toy = 1; toy <= NToys; toy++) {
//      TRandom3 *R = new TRandom3(0);
//
//      //-- Create new histogram for the data-set
//      TH1D *h_bgr_toy = (TH1D*) h_bgr->Clone("h_mass_toydataset"); h_bgr_toy->Reset(); 
//      TH1D *h_sigbgr_toy = (TH1D*) h_sigbgr->Clone("h_mass_toydataset"); h_sigbgr_toy->Reset(); 
//
//      //-- Loop over bins and draw Poisson number of event in each bin
//      for (int i_bin = 1; i_bin <= h_bgr->GetNbinsX(); i_bin++){
//        double toy_bgr = R->Poisson(h_bgr->GetBinContent(i_bin));
//        double toy_sig = R->Poisson(h_sig->GetBinContent(i_bin));
//        h_bgr_toy->SetBinContent(i_bin, toy_bgr);
//        h_sigbgr_toy->SetBinContent(i_bin, toy_bgr+toy_sig);
//      }
//
//
//      h_bgr_TS->Fill(Get_TestStatistic(h_bgr_toy,h_bgr,h_sig));
//      h_sigbgr_TS->Fill(Get_TestStatistic(h_sigbgr_toy,h_bgr,h_sig));
//    }
//
//    vector<double> sigbgr_t_quant = Get_Quantiles(h_sigbgr_TS);
//    vector<double> bgr_t_quant = Get_Quantiles(h_bgr_TS);
//
//    double t_sigbgr_median = sigbgr_t_quant[2];
//    double t_bgr_median = bgr_t_quant[2];
//    double t_data = Get_TestStatistic(h_data,h_bgr,h_sig);
//
//    TH1D *h_bgr_TS_norm = (TH1D*) h_bgr_TS->Clone();
//    TH1D *h_sigbgr_TS_norm = (TH1D*) h_sigbgr_TS->Clone();
//
//    h_bgr_TS_norm->Scale(1/h_bgr_TS_norm->Integral());
//    h_sigbgr_TS_norm->Scale(1/h_sigbgr_TS_norm->Integral());
//    
//    /// 1 - CL_b ///
//    
//    double p1mCLb_sb = 1 - IntegrateFromRight(h_bgr_TS_norm, t_sigbgr_median);
//    double p1mCLb_b = 1 - IntegrateFromRight(h_bgr_TS_norm, t_bgr_median);
//    double p1mCLb_n = 1 - IntegrateFromRight(h_bgr_TS_norm, t_data);
//
//    h_LumiScan_spb->SetBinContent(i_bin,p1mCLb_sb);
//    h_LumiScan_b->SetBinContent(i_bin,p1mCLb_b);
//    h_LumiScan_obs->SetBinContent(i_bin,p1mCLb_n);
//  }
//
//  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
//  canvas2->SetGrid();
//  
//  h_LumiScan_spb->SetLineWidth(2);
//  h_LumiScan_b->SetLineWidth(2);
//  h_LumiScan_obs->SetLineWidth(2);
//  h_LumiScan_spb->SetLineColor(2);
//  h_LumiScan_b->SetLineColor(3);
//  h_LumiScan_obs->SetLineColor(4);
//
//  THStack* stack = new THStack("mystack1", "");
//  stack->Add(h_LumiScan_spb);
//  stack->Add(h_LumiScan_b);
//  stack->Add(h_LumiScan_obs);
//  stack->Draw("l nostack");
//  stack->SetMinimum(ROOT::Math::normal_cdf_c(6));
//  stack->SetMaximum(1);
//  stack->Draw("l nostack");
//
//  canvas2->SetLogy(1);
//
//  TLine *line1 = new TLine(1,ROOT::Math::normal_cdf_c(3),10,ROOT::Math::normal_cdf_c(3));
//  line1->SetLineStyle(9);
//  line1->Draw();
//
//  TLine *line2 = new TLine(1,ROOT::Math::normal_cdf_c(4),10,ROOT::Math::normal_cdf_c(4));
//  line2->SetLineStyle(9);
//  line2->Draw();
//
//  TLine *line3 = new TLine(1,ROOT::Math::normal_cdf_c(5),10,ROOT::Math::normal_cdf_c(5));
//  line3->SetLineStyle(9);
//  line3->Draw();

  return;
}
