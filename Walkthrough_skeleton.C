/*
 * Project:        Exercises 11.1-11.3
 * File:           Walkthrough_skeleton.C
 * Author:         Ivo van Vulpen, Aart Heijboer
 * Version (date): 1.0 (23.06.2013)
 *
 * Copyright (C) 2013, Ivo van Vulpen, Aart Heijboer
 * All rights reserved.
 *
 * Description:
 * A code skeleton for the searches part.
 *
 * This code is distributed with the solution manual to the book
 *
 * Data Analysis in High Energy Physics: A Practical Guide to Statistical Methods,
 * Wiley-VCH (2013),
 * O. Behnke, K. Kroeninger, G. Schott, Th. Schoerner-Sadenius (editors)
 */


#include "TBox.h"
#include "TArrow.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLine.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "Math/QuantFuncMathCore.h"
#include "Math/ProbFuncMathCore.h" // for ROOT::Math::gaussian_cdf

#include <tuple>

#include <iostream>
using namespace std;

//-------------------------------------------------------------------------------------------------------------------------
//-- full functions
TH1D * GetMassDistribution(int Itype = 1, double scalefactor = 1.00);
void MassPlot(int Irebin = 20);

//-- skeleton functions
tuple<double,double,double> SideBandFit(TH1D* h_data, int Irebin = 10, double SB_lower = 150.0, double SB_upper = 400.0, bool plot=true);
tuple<double,double,double> SideBandFit_Toy(int Irebin = 10, double SB_lower = 150.0, double SB_upper = 400.0, bool plot=true, bool printout=true);
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig);
TH1D * GenerateToyDataSet(TH1D *h_mass_template);

//-- some usefull functions
double IntegratePoissonFromRight(double mu, int N_obs);
double IntegrateFromRight(TH1D * h_X_bgr, double X_value);
vector<double> Get_Quantiles( TH1D* hist );
void AddText( Double_t txt_x = 0.50, Double_t txt_y = 0.50, const char * txt = "dummy", Double_t txt_size = 0.045, 
	      Double_t txt_angle = 0., const char * Alignment = "left", Int_t UseNormalizedSize = 1, Int_t txt_color =1 );
//-------------------------------------------------------------------------------------------------------------------------



//========================================
// S O M E   F I N A L   F U N C T I O N S 
//========================================


//========================================================
TH1D * GetMassDistribution(int Itype, double scalefactor){
//========================================================
 //----------------------------------------------------------
 // Goal: return the histogram of the 4-lepton invariant mass
 //       for given type with an optional scale factor
 //
 //       Itype 1 = ZZ SM background
 //             2 = data
 //           125 = Higgs 125
 //           200 = Higgs 200
 //
 //      scalefactor: histograms will be scaled with this number
 //
 //  Note: Histograms have ~200 MeV bins, so need to rebin
 //---------------------------------------------------------

  //-- [1] Get histogram from the file
  TH1D *h_mass = 0;
  TDirectory* dir = gDirectory;   
  TFile *file = new TFile("Histograms_fake.root", "READ");
  dir->cd();

  //-- Higgs 125
  if(Itype == 125){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs125_fake")->Clone("h_mass");     
  }
  //-- Higgs 200
  if(Itype == 200){
    h_mass  = (TH1D*) file->Get("h_m4l_Higgs200_fake")->Clone("h_mass");     
  }
  //-- ZZ SM background
  if(Itype == 1){
    h_mass  = (TH1D*) file->Get("h_m4l_ZZ_fake")->Clone("h_mass");     
  }
  //-- data
  if(Itype == 2){
    h_mass  = (TH1D*) file->Get("h_m4l_data_fake")->Clone("h_mass");     
  }
  
  //-- [2] scale histograms
  int Nbins = h_mass->GetNbinsX();
  for (int i_bin = 1; i_bin < Nbins; i_bin++){
    double mu_bin = h_mass->GetBinContent(i_bin);
    h_mass -> SetBinContent( i_bin, scalefactor * mu_bin);
  }


  file->Close();
  //-- [3] return histogram   
  return h_mass;

  //===========================
} // end GetMassDistribution()
  //===========================




//========================
void MassPlot(int Irebin){
//========================
  // ------------------------------------------
  // Goal: produce SM+Higgs+data plot
  //       Note: rebinning is only for plotting
  // ------------------------------------------

  //------------------------------------
  //-- Standard stuff and prepare canvas
  //------------------------------------
  gROOT->Clear();
  gROOT->Delete();

  //-- Prepare canvas and plot histograms
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o Make cumulative histograms (for signal and background)
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125);
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);

  //-----------------------------------
  //-- [2] Plot histograms and make gif
  //--     o rebin histograms 
  //--     o prepare cumulative histogram
  //--     o make plot + opsmuk + gif
  //-----------------------------------

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);

  //-- Prepare cumulative histogram for signal + background 
  TH1D *h_sig_plus_bgr = (TH1D* ) h_bgr->Clone("h_sig_plus_bgr");
  h_sig_plus_bgr->Reset();
  for (int i_bin = 1; i_bin < h_bgr->GetNbinsX(); i_bin++){
       h_sig_plus_bgr->SetBinContent( i_bin, h_sig->GetBinContent(i_bin) + h_bgr->GetBinContent(i_bin));
       printf("  REBINNED HISTOGRAM:  bin %d, Ndata = %d\n",i_bin,(int)h_data->GetBinContent(i_bin));
  }

  //-- prepare histograms and plot them on canvas
  double Data_max = h_data->GetBinContent(h_data->GetMaximumBin());
  double Ymax_plot = 1.10* (Data_max + TMath::Sqrt(Data_max));
  h_sig_plus_bgr->SetFillColor(7); 
  h_sig_plus_bgr->SetAxisRange(0.,Ymax_plot,"Y");
  h_sig_plus_bgr->SetAxisRange(0.,400.,"X");
  h_bgr->SetFillColor(2); 
  h_sig_plus_bgr->Draw("hist");  
  h_bgr->Draw("same");  
  h_bgr->Draw("axis same");  
  h_data->Draw("e same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-lepton invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.85);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_bgr,          " SM(ZZ)", "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sig_plus_bgr, " Higgs" , "f");  leg1b->SetTextSize(0.04);
  leg1->Draw();

  //-- prepare gif
  canvas1->Print(Form("./MassPlot_rebin%d.gif",Irebin));

  return;

   //===============
 } // end MassPlot()
   //===============




//===============================================
// S O M E   S K E L E T O N    F U N C T I O N S 
//===============================================


//=============================================================
void Significance_Optimization(double Lumi_scalefactor = 1.00, double Higgs_mass = 125){
//=============================================================

  printf("\n Significance_Optimization()\n\n");

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o scale to correct luminosity
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  //h_sig  = GetMassDistribution(Higgs_mass, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-------------------------------------------------------
  //-- [2] Compute significance for various mass windows
  //--     o try various options for window (use histogram)
  //--     o compute expected and observed significance 
  //-------------------------------------------------------

  //-- Define histogram that defines mass windows to try
  TH1D *h_masswindow          = new TH1D("h_masswindow","",250,0.,25.);           // make a mass window - full width between 0 and 25 GeV
  TH1D *h_masswindow_expected = new TH1D("h_masswindow_expected","",250,0.,25.);  // histogram to hold results for expected
  TH1D *h_masswindow_observed = new TH1D("h_masswindow_observed","",250,0.,25.);  // histogram to hold results for observed

  //---------------------------------
  //-- Loop over various mass windows
  //---------------------------------
  for (int i_bin = 1; i_bin<h_masswindow->GetNbinsX(); i_bin++ ){
  
    //-- get full width of mass window (bin center of our histogram) and the number of events in mass window for each event type
    double masswindow_fullwidth = h_masswindow->GetBinCenter(i_bin);    
    printf("   Trying as mass window: %5.2f GeV\n", masswindow_fullwidth);

    double masswindow_lower = Higgs_mass-(masswindow_fullwidth/2);
    double masswindow_upper = Higgs_mass+(masswindow_fullwidth/2);

    //-- [a] determine the number of events in the mass window for each event type
    //       Ndata_win, Nbgr_win and Nsig_win
    int lower_bin_sig = h_sig->FindBin(masswindow_lower); 
    int upper_bin_sig = h_sig->FindBin(masswindow_upper); 
    int lower_bin_bgr = h_bgr->FindBin(masswindow_lower); 
    int upper_bin_bgr = h_bgr->FindBin(masswindow_upper); 
    int lower_bin_data = h_data->FindBin(masswindow_lower); 
    int upper_bin_data = h_data->FindBin(masswindow_upper); 

    double Nsig_win = h_sig->Integral(lower_bin_sig, upper_bin_sig);
    double Nbgr_win = h_bgr->Integral(lower_bin_bgr, upper_bin_bgr);
    double Ndata_win = h_data->Integral(lower_bin_data, upper_bin_data);

    //double p_exp = IntegratePoissonFromRight(Nbgr_win, Nsig_win+Nbgr_win);
    //double p_obs = IntegratePoissonFromRight(Nbgr_win, Ndata_win);

    //printf(" %f\n", Nsig_win);
    //printf(" %f\n", Nbgr_win);
    //printf(" %f\n", Ndata_win);

    //printf(" %f\n", p_exp);
    //printf(" %f\n", p_obs);

    //-- [b] compute EXPECTED significance and save in histogram
    //double pvalue_expected       = 0.75; // you need to do this yourself
    double pvalue_expected       = IntegratePoissonFromRight(Nbgr_win, Nsig_win+Nbgr_win);
    double significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    if ( !isfinite(significance_expected)) {
      significance_expected = 0;
    }
    h_masswindow_expected->SetBinContent(i_bin, significance_expected);

    //-- [c] compute OBSERVED significance and save in histogram
    //double pvalue_observed       = 0.25; // you need to do this yourself
    double pvalue_observed       = IntegratePoissonFromRight(Nbgr_win, Ndata_win);
    double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
    if ( !isfinite(significance_observed)) {
      significance_observed = 0;
    }
    h_masswindow_observed->SetBinContent(i_bin, significance_observed);

    printf("           expected significance = %.3f (p = %.6f) || observed significance = %.3f (p = %.6f) \n", significance_expected, pvalue_expected, significance_observed, pvalue_observed);

  } // end loop over width mass window

  int binmax_expected = h_masswindow_expected->GetMaximumBin();
  double optimum_exp_significance = h_masswindow_expected->GetBinContent(binmax_expected);
  double optimum_exp_masswindow = h_masswindow_expected->GetXaxis()->GetBinCenter(binmax_expected);

  printf(" Lumi factor: %2.2f \n", Lumi_scalefactor);
  //-- print optimum to the screen
  printf(" Optimal mass window (expected): \n"); 
  printf("    Expected significance: %.4f \n", optimum_exp_significance); 
  printf("    Mass window width: %.2f GeV \n", optimum_exp_masswindow); 
  printf("    Mass window range: [%.2f, %.2f] GeV \n", (Higgs_mass-optimum_exp_masswindow), (Higgs_mass+optimum_exp_masswindow)); 
  
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    int binmax_observed = h_masswindow_observed->GetMaximumBin();
    double optimum_obs_significance = h_masswindow_observed->GetBinContent(binmax_observed);
    double optimum_obs_masswindow = h_masswindow_observed->GetXaxis()->GetBinCenter(binmax_observed);

    //-- print optimum to the screen
    printf(" Optimal mass window (observed): \n"); 
    printf("    Observed significance: %.4f \n", optimum_obs_significance); 
    printf("    Mass window width: %.2f GeV \n", optimum_obs_masswindow); 
    printf("    Mass window range: [%.2f, %.2f] GeV \n", (Higgs_mass-optimum_obs_masswindow), (Higgs_mass+optimum_obs_masswindow)); 
  }

  //----------------------------------
  //-- [3] Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125); 
  canvas1->cd(); 
 
  h_masswindow_expected->SetLineColor(1);
  h_masswindow_expected->SetLineWidth(2);
  h_masswindow_observed->SetLineColor(4);
  h_masswindow_observed->SetLineWidth(2);

  h_masswindow_expected->SetAxisRange(-1.,6.,"Y");
  h_masswindow_expected->Draw("l");
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    h_masswindow_observed->Draw("l same");
  }
  //-- axes
  AddText( 0.900, 0.035, "Mass window GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    
  AddText( 0.225, 0.825, Form("Luminosity scalefactor = %5.1f",Lumi_scalefactor),0.050, 0.,"left");            

  AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);                        
  if(fabs(Lumi_scalefactor-1.00)<0.01){
    AddText( 0.700, 0.300, "Observed significance",0.050, 0.,"right",1,4);                        
  }
  //-- prepare gif
  canvas1->Print(Form("./Significance_Optimization_lumiscalefactor%d.gif",int(Lumi_scalefactor)));

  return;

  //================================
} // end Significance_Optimization()
  //================================


//=============================================================
double Significance(double Lumi_scalefactor = 1.00, double masswindow_fullwidth = 10.0, double rel_sigma_b = 0.0, double Higgs_mass = 125){
//=============================================================

  //printf("\n Significance()\n\n");

  //------------------------------------------------------------------
  //-- [1] Prepare histograms
  //--     o Get histograms from the files (signal, background and data)
  //--     o scale to correct luminosity
  //------------------------------------------------------------------

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  //printf ("\n  INFO: Mass distribution in the 4 lepton channel\n");
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  //h_sig  = GetMassDistribution(Higgs_mass, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-------------------------------------------------------
  //-- [2] Compute significance of mass windows
  //--     o compute expected and observed significance 
  //-------------------------------------------------------

  //printf("   Mass window: %5.2f GeV\n", masswindow_fullwidth);

  double masswindow_lower = Higgs_mass-(masswindow_fullwidth/2);
  double masswindow_upper = Higgs_mass+(masswindow_fullwidth/2);

  //-- [a] determine the number of events in the mass window for each event type
  //       Ndata_win, Nbgr_win and Nsig_win
  int lower_bin_sig = h_sig->FindBin(masswindow_lower); 
  int upper_bin_sig = h_sig->FindBin(masswindow_upper); 
  int lower_bin_bgr = h_bgr->FindBin(masswindow_lower); 
  int upper_bin_bgr = h_bgr->FindBin(masswindow_upper); 
  int lower_bin_data = h_data->FindBin(masswindow_lower); 
  int upper_bin_data = h_data->FindBin(masswindow_upper); 

  double Nsig_win = h_sig->Integral(lower_bin_sig, upper_bin_sig);
  double Nbgr_win = h_bgr->Integral(lower_bin_bgr, upper_bin_bgr);
  double Ndata_win = h_data->Integral(lower_bin_data, upper_bin_data);

  //cout << rel_sigma_b << endl;
  double significance_expected;
  if ( rel_sigma_b == 0 ) {
    //-- [b] compute EXPECTED significance and save in histogram
    //double pvalue_expected       = 0.75; // you need to do this yourself
    double pvalue_expected       = IntegratePoissonFromRight(Nbgr_win, Nsig_win+Nbgr_win);
    significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    if ( !isfinite(significance_expected)) {
      significance_expected = 0;
    }

    //-- [c] compute OBSERVED significance and save in histogram
    //double pvalue_observed       = 0.25; // you need to do this yourself
    double pvalue_observed       = IntegratePoissonFromRight(Nbgr_win, Ndata_win);
    double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
    if ( !isfinite(significance_observed)) {
      significance_observed = 0;
    }
  }
  else {
    TH1D *h_bkg_dist = new TH1D("h_bkg_dist","",201,0,201);           // make a mass window - full width between 0 and 25 GeV
    int n;
    double upper_bound = Nbgr_win + 10*(Nbgr_win*rel_sigma_b);
    double step = 0.01*(Nbgr_win*rel_sigma_b);
    for (int i_bin = 1; i_bin <= h_bkg_dist->GetNbinsX(); i_bin++){
         n = h_bkg_dist->GetBinLowEdge(i_bin);
         //cout << i_bin << " " << TMath::Poisson(n, Nbgr_win) << endl;
         //cout << i_bin << " " << TMath::Gaus(n, Nbgr_win, Nbgr_win*rel_sigma_b) << endl;
         double n_term = 0;
         for (double nu_b = 0; nu_b<= upper_bound; nu_b+=step){
            n_term += TMath::Poisson(n, nu_b)*TMath::Gaus(nu_b, Nbgr_win, Nbgr_win*rel_sigma_b);
         }
         h_bkg_dist->SetBinContent(i_bin, n_term);
         //h_bkg_dist->SetBinContent(i_bin, TMath::Poisson(n, Nbgr_win)*TMath::Gaus(n, Nbgr_win, Nbgr_win*rel_sigma_b));
    }
    h_bkg_dist->Scale(1/h_bkg_dist->Integral());
    
    double pvalue_expected = IntegrateFromRight(h_bkg_dist, Nsig_win+Nbgr_win); //IntegratePoissonFromRight(Nbgr_win, Nsig_win+Nbgr_win);
    significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
    double S = Nsig_win;
    double B = Nbgr_win;
    double D_B = Nbgr_win*rel_sigma_b;
    cout << 2*(TMath::Sqrt(S + B) - TMath::Sqrt(B))*B/(B + TMath::Power(D_B, 2)) << endl;
    //significance_expected = 2*(TMath::Sqrt(S + B) - TMath::Sqrt(B))*B/(B + TMath::Power(D_B, 2));
    if ( !isfinite(significance_expected)) {
      significance_expected = 0;
    }
  }

//  printf("           expected significance = %.3f (p = %.6f) || observed significance = %.3f (p = %.6f) \n", significance_expected, pvalue_expected, significance_observed, pvalue_observed);
//
//  printf(" Lumi factor: %2.2f \n", Lumi_scalefactor);
//  //-- print optimum to the screen
//  printf(" Mass window (expected): \n"); 
//  printf("    Expected significance: %.4f \n", significance_expected); 
//  printf("    Mass window width: %.2f GeV \n", masswindow_fullwidth); 
//  printf("    Mass window range: [%.2f, %.2f] GeV \n", (Higgs_mass-masswindow_fullwidth), (Higgs_mass+masswindow_fullwidth)); 
//  
//  if(fabs(Lumi_scalefactor-1.00)<0.01){
//    //-- print optimum to the screen
//    printf(" Mass window (observed): \n"); 
//    printf("    Observed significance: %.4f \n", significance_observed); 
//    printf("    Mass window width: %.2f GeV \n", masswindow_fullwidth); 
//    printf("    Mass window range: [%.2f, %.2f] GeV \n", (Higgs_mass-masswindow_fullwidth), (Higgs_mass+masswindow_fullwidth)); 
//  }


  return significance_expected;

  //================================
} // end Significance()
  //================================




//===========================
tuple<double,double,double> SideBandFit(int Irebin, double SB_lower, double SB_upper, bool plot){
//===========================

  printf("\n SideBandFit()\n\n");

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr, *h_data;
  h_bgr  = GetMassDistribution(1);
  h_data = GetMassDistribution(2);
 
  //-- rebin histograms if necessary
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);
  printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.2f GeV\n", Irebin, h_data->GetBinWidth(1));

  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",10000,0.1,3.1); // you'll want to put more here I guess
  double scalefactor_bgr  = 1.00; 

  //---------------------------------------------
  // [2a] Loop 1: loop over scalefactors in alpha
  //---------------------------------------------
  int lower_bin = h_data->FindBin(SB_lower); 
  int upper_bin = h_data->FindBin(SB_upper); 
  for (int i_bin_sf = 1; i_bin_sf <=h_scalefactor_bgr->GetNbinsX(); i_bin_sf++){

    //-- determine the scale factor for the background
    scalefactor_bgr = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
    printf(" Loop 1: I am now trying alpha = %5.3f\n",scalefactor_bgr);
  
    //-----------------------------------------------------------------------------------
    // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
    //-----------------------------------------------------------------------------------
    double loglik = 1e-10;
    for (int i_bin = lower_bin; i_bin <= upper_bin; i_bin++){
      //printf("        bin %d, m = %5.2f, Ndata = %5.2f, Nbgr = %5.2f\n",i_bin, h_bgr->GetBinCenter(i_bin),h_data->GetBinContent(i_bin),h_bgr->GetBinContent(i_bin));
      double ni = h_data->GetBinContent(i_bin);
      double bi = h_bgr->GetBinContent(i_bin);
      //loglik += ni*TMath::Log(scalefactor_bgr*bi) - scalefactor_bgr*bi; // - TMath::Log(TMath::Factorial(ni));
      loglik += TMath::Log(TMath::Poisson(ni,scalefactor_bgr*bi)); // - TMath::Log(TMath::Factorial(ni));
    } // end loop over bins

      h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik);   
  } // end loop over scale factors for the background


  //----------------------------------------------------
  //-- [3] Interpret the -2Log (Likelihood distribution)
  //----------------------------------------------------
  TH1D *h_scalefactor_bgr_rescaled  = (TH1D*) h_scalefactor_bgr->Clone("h_scalefactor_bgr_rescaled");  
  //h_scalefactor_bgr_rescaled->Reset();   
  
  //-- Find minimum
  int binmin = h_scalefactor_bgr_rescaled->GetMinimumBin();
  double min = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  double alpha_fit = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(binmin);

  //--Rescale and find +/- 1 sigma errors
  for (int i_bin = 1; i_bin <=h_scalefactor_bgr_rescaled->GetNbinsX(); i_bin++){
    h_scalefactor_bgr_rescaled->AddBinContent(i_bin, -min);
  }

  int bin_uncert_plus = binmin;
  double val_uncert_plus = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  while ( val_uncert_plus <= 1 ) {
    bin_uncert_plus++;
    val_uncert_plus = h_scalefactor_bgr_rescaled->GetBinContent(bin_uncert_plus);
  }
  double alpha_fit_sigma_max = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(bin_uncert_plus);
  
  int bin_uncert_minus = binmin;
  double val_uncert_minus = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  while ( val_uncert_minus <= 1 ) {
    bin_uncert_minus--;
    val_uncert_minus = h_scalefactor_bgr_rescaled->GetBinContent(bin_uncert_minus);
  }
  double alpha_fit_sigma_min = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(bin_uncert_minus);

  double alpha_fit_uncert_plus = alpha_fit_sigma_max-alpha_fit;
  double alpha_fit_uncert_minus = alpha_fit-alpha_fit_sigma_min;

  //-- print summary to screen
  printf("   Summary\n");
  printf("      alpha_fit = %5.2f  [%.2f - %.2f] \n", alpha_fit, alpha_fit_sigma_min, alpha_fit_sigma_max);
  printf("      Uncertainty: + %5.2f \n", alpha_fit_uncert_plus);
  printf("                   - %5.2f \n", alpha_fit_uncert_minus);
  

  if ( plot ) {
    TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
    canvas1->SetGrid();
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125); 
    canvas1->cd(); 

    gStyle->SetOptStat(0);
    
    h_scalefactor_bgr_rescaled->SetLineWidth(2);
    h_scalefactor_bgr_rescaled->SetAxisRange(0.0,5.0,"Y");
    h_scalefactor_bgr_rescaled->SetAxisRange(0.95,1.25,"X");
    h_scalefactor_bgr_rescaled->Draw("l");

    //-- axes
    AddText( 0.900, 0.035, "Background scale factor",0.060, 0.,"right"); // X-axis
    AddText( 0.040, 0.900, "(-2ln(L)) - min(-2ln(L))" ,0.060,90.,"right");   // Y-axis
    
    TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);  
    canvas2->SetGrid();
    canvas2->SetLeftMargin(0.125);
    canvas2->SetBottomMargin(0.125); 
    canvas2->cd(); 

    h_scalefactor_bgr->SetLineWidth(2);
    h_scalefactor_bgr->Draw("l");
    
    //-- axes
    AddText( 0.900, 0.035, "Background scale factor",0.060, 0.,"right"); // X-axis
    AddText( 0.040, 0.900, "-2ln(L)" ,0.060,90.,"right");   // Y-axis
  }

  tuple<double,double,double> results = make_tuple(alpha_fit, alpha_fit_sigma_min, alpha_fit_sigma_max);

  return results;

  //==================
} // end SideBandFit()
  //==================

  
//===========================
tuple<double,double,double> SideBandFit_Toy(TH1D* h_data, int Irebin, double SB_lower, double SB_upper, bool plot, bool printout){
//===========================

  if (printout){
    printf("\n SideBandFit()\n\n");
  }

  //-------------------------
  //-- [1] Prepare histograms
  //-------------------------
  TH1D *h_bgr;//, *h_data;
  h_bgr  = GetMassDistribution(1);
  //h_data = GetMassDistribution(2);
 
  //-- rebin histograms if necessary
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);
  if (printout) {
    printf(" INFO: Rebinning the histograms with a factor %d. Binwidth is now %5.2f GeV\n", Irebin, h_data->GetBinWidth(1));
  }

  //-----------------------------------------
  //-- [2] Loop over scale factor (alpha_bgr)
  //-----------------------------------------
  TH1D *h_scalefactor_bgr = new TH1D("h_scalefactor_bgr","",2000,0.1,3.1); // you'll want to put more here I guess
  double scalefactor_bgr  = 1.00; 

  //---------------------------------------------
  // [2a] Loop 1: loop over scalefactors in alpha
  //---------------------------------------------
  int lower_bin = h_data->FindBin(SB_lower); 
  int upper_bin = h_data->FindBin(SB_upper); 
  for (int i_bin_sf = 1; i_bin_sf <=h_scalefactor_bgr->GetNbinsX(); i_bin_sf++){

    //-- determine the scale factor for the background
    scalefactor_bgr = h_scalefactor_bgr->GetBinCenter(i_bin_sf);
    if (printout) {
      printf(" Loop 1: I am now trying alpha = %5.3f\n",scalefactor_bgr);
    }
  
    //-----------------------------------------------------------------------------------
    // [2b] Loop 2: loop over bins in the histogram, compute loglik and save in histogram
    //-----------------------------------------------------------------------------------
    double loglik = 1e-10;
    for (int i_bin = lower_bin; i_bin <= upper_bin; i_bin++){
      //printf("        bin %d, m = %5.2f, Ndata = %5.2f, Nbgr = %5.2f\n",i_bin, h_bgr->GetBinCenter(i_bin),h_data->GetBinContent(i_bin),h_bgr->GetBinContent(i_bin));
      double ni = h_data->GetBinContent(i_bin);
      double bi = h_bgr->GetBinContent(i_bin);
      //loglik += ni*TMath::Log(scalefactor_bgr*bi) - scalefactor_bgr*bi; // - TMath::Log(TMath::Factorial(ni));
      loglik += TMath::Log(TMath::Poisson(ni,scalefactor_bgr*bi)); // - TMath::Log(TMath::Factorial(ni));
    } // end loop over bins

      h_scalefactor_bgr->SetBinContent(i_bin_sf,-2.*loglik);   
  } // end loop over scale factors for the background


  //----------------------------------------------------
  //-- [3] Interpret the -2Log (Likelihood distribution)
  //----------------------------------------------------
  TH1D *h_scalefactor_bgr_rescaled  = (TH1D*) h_scalefactor_bgr->Clone("h_scalefactor_bgr_rescaled");  
  //h_scalefactor_bgr_rescaled->Reset();   
  
  //-- Find minimum
  int binmin = h_scalefactor_bgr_rescaled->GetMinimumBin();
  double min = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  double alpha_fit = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(binmin);

  //--Rescale and find +/- 1 sigma errors
  for (int i_bin = 1; i_bin <=h_scalefactor_bgr_rescaled->GetNbinsX(); i_bin++){
    h_scalefactor_bgr_rescaled->AddBinContent(i_bin, -min);
  }

  int bin_uncert_plus = binmin;
  double val_uncert_plus = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  while ( val_uncert_plus <= 1 ) {
    bin_uncert_plus++;
    val_uncert_plus = h_scalefactor_bgr_rescaled->GetBinContent(bin_uncert_plus);
  }
  double alpha_fit_sigma_max = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(bin_uncert_plus);
  
  int bin_uncert_minus = binmin;
  double val_uncert_minus = h_scalefactor_bgr_rescaled->GetBinContent(binmin);
  while ( val_uncert_minus <= 1 ) {
    bin_uncert_minus--;
    val_uncert_minus = h_scalefactor_bgr_rescaled->GetBinContent(bin_uncert_minus);
  }
  double alpha_fit_sigma_min = h_scalefactor_bgr_rescaled->GetXaxis()->GetBinCenter(bin_uncert_minus);

  double alpha_fit_uncert_plus = alpha_fit_sigma_max-alpha_fit;
  double alpha_fit_uncert_minus = alpha_fit-alpha_fit_sigma_min;

  if (printout) {
    //-- print summary to screen
    printf("   Summary\n");
    printf("      alpha_fit = %5.2f  [%.2f - %.2f] \n", alpha_fit, alpha_fit_sigma_min, alpha_fit_sigma_max);
    printf("      Uncertainty: + %5.2f \n", alpha_fit_uncert_plus);
    printf("                   - %5.2f \n", alpha_fit_uncert_minus);
  }
  

  if ( plot ) {
    TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);  
    canvas1->SetGrid();
    canvas1->SetLeftMargin(0.125);
    canvas1->SetBottomMargin(0.125); 
    canvas1->cd(); 

    h_scalefactor_bgr_rescaled->SetAxisRange(0.0,5.0,"Y");
    h_scalefactor_bgr_rescaled->SetAxisRange(0.95,1.25,"X");
    h_scalefactor_bgr_rescaled->Draw("l");
    
    TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);  
    canvas2->SetGrid();
    canvas2->SetLeftMargin(0.125);
    canvas2->SetBottomMargin(0.125); 
    canvas2->cd(); 

    h_scalefactor_bgr->Draw("l");
  }

  tuple<double,double,double> results = make_tuple(alpha_fit, alpha_fit_sigma_min, alpha_fit_sigma_max);

  delete h_scalefactor_bgr;
  h_scalefactor_bgr = NULL;

  return results;

  //==================
} // end SideBandFit()
  //==================



//=========================================================================================
double Get_TestStatistic(TH1D *h_mass_dataset, TH1D *h_template_bgr, TH1D *h_template_sig){
//=========================================================================================
   
  double loglik_bgr = 1e-10;
  double loglik_sig_plus_bgr = 1e-10;

  int Nbins = h_mass_dataset->GetNbinsX();
  double n, b, s;
  //-- do likelihood fit
  for (int i_bin = 1; i_bin <= Nbins; i_bin++){
    n = h_mass_dataset->GetBinContent(i_bin);
    b = h_template_bgr->GetBinContent(i_bin);
    s = h_template_sig->GetBinContent(i_bin);
    //loglik_bgr += n*TMath::Log(b) - b;  // likelihood for the mu=0 (b-only) scenario
    //loglik_sig_plus_bgr += b*TMath::Log(b+s) - (b+s);  // likelihood for the mu=1 (s+b) scenario
    loglik_bgr += TMath::Log(TMath::Poisson(n,b));  // likelihood for the mu=0 (b-only) scenario
    loglik_sig_plus_bgr += TMath::Log(TMath::Poisson(n,s+b));  // likelihood for the mu=1 (s+b) scenario
  } // end loop over bins
 
  //-- compute test statistic
  double test_statistic  = -2*(loglik_sig_plus_bgr - loglik_bgr); // find this out yourself

  //-- return test_statistic
  return test_statistic;

} // end Get_TestStatistic()





//===============================================
TH1D * GenerateToyDataSet(TH1D *h_mass_template){
//===============================================
  //-------------------------------------------------
  // Goal: Generate Toy data set from input histogram
  // How:  dumb way -> draw random Poisson in each bin       
  //--------------------------------------------------
  TRandom3 *R = new TRandom3(0);

  //-- Create new histogram for the data-set
  TH1D *h_mass_toydataset = (TH1D*) h_mass_template->Clone("h_mass_toydataset"); h_mass_toydataset->Reset();
 
  //-- Loop over bins and draw Poisson number of event in each bin
  for (int i_bin = 1; i_bin <= h_mass_toydataset->GetNbinsX(); i_bin++){
    h_mass_toydataset->SetBinContent(i_bin, R->Poisson(h_mass_template->GetBinContent(i_bin)));
  }

  //-- return histogram of toy data-set
  return h_mass_toydataset;
  
} // end GenerateToyDataSet()

//=============================================================
void ExpectedObservedSignificance_ToyMC(double b = 5.15, double s = 5.41, double sigma_b = 0.30, int Ntoys = 1e3, int Ndata = 13){
//=============================================================

  TRandom3 *R = new TRandom3(0);

  int Nbins = 60;

  TH1D *h_bgr_toys = new TH1D("h_bgr_toys","", Nbins, 0, Nbins); 
  TH1D *h_sigbgr_toys = new TH1D("h_sigbgr_toys","", Nbins, 0, Nbins); 

  for (int toy = 1; toy <= Ntoys; toy++){
    h_bgr_toys->Fill(R->Poisson(R->Gaus(b, sigma_b)));
    h_sigbgr_toys->Fill(R->Poisson(s+R->Gaus(b, sigma_b)));
  }

  vector<double> bgr_toys_quant = Get_Quantiles(h_bgr_toys);
  vector<double> sigbgr_toys_quant = Get_Quantiles(h_sigbgr_toys);

  cout << bgr_toys_quant[2] << endl;
  cout << h_bgr_toys->GetMean() << endl;
  cout << sigbgr_toys_quant[2] << endl;
  cout << h_sigbgr_toys->GetMean() << endl;



  h_bgr_toys->Draw();
  h_sigbgr_toys->Draw("same");

  //-- compute EXPECTED significance
  //double pvalue_expected       = IntegratePoissonFromRight(bgr_toys_quant[2], sigbgr_toys_quant[2]);
  double pvalue_expected       = IntegratePoissonFromRight(h_bgr_toys->GetMean(), h_sigbgr_toys->GetMean());
  double significance_expected = ROOT::Math::gaussian_quantile_c(pvalue_expected,1);
  if ( !isfinite(significance_expected)) {
    significance_expected = 0;
  }

  //-- compute OBSERVED significance
  //double pvalue_observed       = IntegratePoissonFromRight(bgr_toys_quant[2], Ndata);
  double pvalue_observed       = IntegratePoissonFromRight(h_bgr_toys->GetMean(), Ndata);
  double significance_observed = ROOT::Math::gaussian_quantile_c(pvalue_observed,1);
  if ( !isfinite(significance_observed)) {
    significance_observed = 0;
  }

  printf(" Summary: \n"); 
  printf("  Expected background: %.2f \n", b); 
  printf("  Expected background uncertainty (gauss): %.2f \n", sigma_b); 
  printf("  Expected signal: %.2f \n", s); 
  printf("  Observed data: %d \n", Ndata); 
  printf("  Number of toys: %.1e \n", double(Ntoys)); 
  printf("  \n"); 
  printf("    Expected significance: %.4f \n", significance_expected); 
  printf("    Observed significance: %.4f \n", significance_observed); 

  cout << "whoho!" << endl;
  
} // end ExpectedObservedSignificance_ToyMC()
















//================================================
// S O M E   U S E F U L L   F U N C T I O N S 
//================================================



//====================================================
double IntegratePoissonFromRight(double mu, int N_obs){
//====================================================
// --------------------------------------------------------------------
// Compute p-value for case zero background uncertainty, i.e. 
//         just integrate Poisson from the right from N_obs to infinity
// --------------------------------------------------------------------

  double integral = 1.; 
  for(int i_obs = 0; i_obs < N_obs; i_obs++){
    integral -= TMath::Poisson(i_obs,mu);
  }
  
  return integral;
  
} // end IntegratePoissonFromRight()


//========================================================
double IntegrateFromRight(TH1D * h_X_bgr, double X_value){
//========================================================
// --------------------------------------------------------------------
// Compute p-value: integrate number of events from X_value to infinity
// --------------------------------------------------------------------

  //-- Integrate distributions
  int Nbins = h_X_bgr->GetNbinsX();
  int X_bin = h_X_bgr->FindBin(X_value); 

  //-- Compute integral from X-value to infinity
  double pvalue = h_X_bgr->Integral(X_bin,Nbins) / h_X_bgr->Integral();
  
  return pvalue;

} // end IntegrateFrom Right()




//=========================================
vector<double> Get_Quantiles( TH1D* hist ){
//=========================================
// Quantiles returns a vector<double> with 5 entries.
// Entries 0 and 4 are the values on the histogram x-axis
// so that 95% of the content lies between these values.
// Entries 1 and 3 bracket 68% in the same way.
// Entry 2 is the median of the histogram.

  //-- define quantiles
  double fraction_1sigma = ROOT::Math::gaussian_cdf(-1.,1.,0.); // 15.8655 %
  double fraction_2sigma = ROOT::Math::gaussian_cdf(-2.,1.,0.); //  2.2750 %
  double probs[5] = {fraction_2sigma, fraction_1sigma, 0.50, 1.00-fraction_1sigma, 1.00-fraction_2sigma };

  //-- output of the quantiles
  double Xvalues[5];

  //-- extract quantiles
  hist->GetQuantiles( 5, Xvalues, probs );
  
  vector<double> Xvalues_output(5);
  for (int i=0; i<5; i++) 
    {
      Xvalues_output[i] = Xvalues[i];
    }

  return Xvalues_output;
} // end Get_Quantiles()





//=======================================================================================================================
void AddText( Double_t txt_x, Double_t txt_y, const char * txt, Double_t txt_size,
              Double_t txt_angle, const char * Alignment, Int_t UseNormalizedSize, Int_t txt_color)
//=======================================================================================================================
{
  Int_t txt_align = 12;
  if ( !strcmp(Alignment, "left"))   { txt_align = 12; } // left 
  if ( !strcmp(Alignment, "right"))  { txt_align = 32; } // right
  if ( !strcmp(Alignment, "center")) { txt_align = 22; } // center

  TLatex* t1 = new TLatex( txt_x, txt_y, txt);
  if(UseNormalizedSize) {t1->SetNDC(kTRUE);} // <- use NDC coordinate
  t1->SetTextSize(txt_size);
  t1->SetTextAlign(txt_align);
  t1->SetTextAngle(txt_angle);
  t1->SetTextColor(txt_color);
  t1->Draw();

} // end AddText()

