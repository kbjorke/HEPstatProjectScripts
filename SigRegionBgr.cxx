#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void SigRegionBgr(){

  tuple<double,double,double> alpha_fit = SideBandFit(1, 150.0, 400.0, false);
  
  double alpha_fit_center = get<0>(alpha_fit);
  double alpha_fit_lower = get<1>(alpha_fit);
  double alpha_fit_upper = get<2>(alpha_fit);

  double Higgs_mass = 125.0;
  double masswindow_width = 7.15;
  double Lumi_scalefactor = 1.0;

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  double masswindow_lower = Higgs_mass-(masswindow_width/2);
  double masswindow_upper = Higgs_mass+(masswindow_width/2);

  int lower_bin_sig = h_sig->FindBin(masswindow_lower);
  int upper_bin_sig = h_sig->FindBin(masswindow_upper);
  int lower_bin_bgr = h_bgr->FindBin(masswindow_lower);
  int upper_bin_bgr = h_bgr->FindBin(masswindow_upper);
  int lower_bin_data = h_data->FindBin(masswindow_lower);
  int upper_bin_data = h_data->FindBin(masswindow_upper);

  double Nsig_win = h_sig->Integral(lower_bin_sig, upper_bin_sig);
  double Nbgr_win = h_bgr->Integral(lower_bin_bgr, upper_bin_bgr);
  double Ndata_win = h_data->Integral(lower_bin_data, upper_bin_data);

  double Nbgr_est_win_center = alpha_fit_center*Nbgr_win;
  double Nbgr_est_win_lower =  alpha_fit_lower*Nbgr_win;
  double Nbgr_est_win_upper =  alpha_fit_upper*Nbgr_win;
  
  double Nbgr_est_win_uncert_plus = Nbgr_est_win_upper-Nbgr_est_win_center;
  double Nbgr_est_win_uncert_minus = Nbgr_est_win_center-Nbgr_est_win_lower;

  printf(" Background estimates: \n");
  printf("   Higgs mass:  %.1f GeV\n", Higgs_mass);
  printf("   Mass Window width:  %.2f GeV\n", masswindow_width);
  printf("   Luminosity scale factor:  %.1f \n", Lumi_scalefactor);
  printf("   \n");
  printf("   Expected background (unscaled):  %5.2f \n", Nbgr_win);
  printf("   Estimated expected background (scaled):  %5.2f   [%.2f,%.2f] \n", Nbgr_est_win_center, Nbgr_est_win_lower, Nbgr_est_win_upper);
  printf("      Uncertanity:  + %.2f \n", Nbgr_est_win_uncert_plus);
  printf("                    - %.2f \n", Nbgr_est_win_uncert_minus);
  printf("   Expected signal:  %5.2f \n", Nsig_win);
  printf("   Observed data:  %5d \n", int(Ndata_win));

  return;
}
