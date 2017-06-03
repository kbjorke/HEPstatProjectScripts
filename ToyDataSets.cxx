#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void ToyDataSets(){

  int NToys = 10000;
  double Lumi_scalefactor = 1.0;
  
  double Higgs_mass = 125.0;
  double masswindow_width = 7.15;
  
  double masswindow_lower = Higgs_mass-(masswindow_width/2);
  double masswindow_upper = Higgs_mass+(masswindow_width/2);
  
  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);
  
  int lower_bin_sig = h_sig->FindBin(masswindow_lower);
  int upper_bin_sig = h_sig->FindBin(masswindow_upper);

  double N_bgr_mw = h_bgr->Integral(lower_bin_sig,upper_bin_sig);

  int Nbins = 200;

  /////// Scaled:

  double alpha_bgr_nom = 1.11;
  double alpha_bgr_sig = 0.07;

  h_bgr->Scale(alpha_bgr_nom);

  ////////////////////
  
  TH1D *h_bgr_TS = new TH1D("h_bgr_TS","", Nbins, -35, 20);
  TH1D *h_sigbgr_TS = new TH1D("h_sigbgr_TS","", Nbins, -35, 20);

  TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
  h_sigbgr->Add(h_bgr);
  
  TH1D *h_Pull = new TH1D("Pull distribution","Pull distribution", 40, -5, 5);
  
  TH1D *h_unitGaus = new TH1D("h_unitGaus","", 40, -5, 5);

  double N_MC_truth = h_bgr->Integral(lower_bin_sig,upper_bin_sig);

  for (int toy = 1; toy <= NToys; toy++) {
    TRandom3 *R = new TRandom3(0);
    
    //-- Create new histogram for the data-set
    TH1D *h_bgr_toy = (TH1D*) h_bgr->Clone("h_bgr_toy"); h_bgr_toy->Reset(); 
    TH1D *h_sigbgr_toy = (TH1D*) h_sigbgr->Clone("h_sigbgr_toy"); h_sigbgr_toy->Reset(); 

    //-- Loop over bins and draw Poisson number of event in each bin
    for (int i_bin = 1; i_bin <= h_bgr->GetNbinsX(); i_bin++){
      double toy_bgr = R->Poisson(h_bgr->GetBinContent(i_bin));
      double toy_sig = R->Poisson(h_sig->GetBinContent(i_bin));
      h_bgr_toy->SetBinContent(i_bin, toy_bgr);
      h_sigbgr_toy->SetBinContent(i_bin, toy_bgr+toy_sig);
    }


    h_bgr_TS->Fill(Get_TestStatistic(h_bgr_toy,h_bgr,h_sig));
    h_sigbgr_TS->Fill(Get_TestStatistic(h_sigbgr_toy,h_bgr,h_sig));

    tuple<double,double,double> alpha_fit = SideBandFit_Toy(h_sigbgr_toy,1,150,400,false,false);

    double alpha_fit_center = get<0>(alpha_fit);
    double alpha_fit_lower = get<1>(alpha_fit);
    double alpha_fit_upper = get<2>(alpha_fit);

    double N_fit_predicted = alpha_fit_center*N_bgr_mw;
    double N_fit_predicted_lower =  alpha_fit_lower*N_bgr_mw;
    double N_fit_predicted_upper =  alpha_fit_upper*N_bgr_mw;


    double sigma_N_fit_predicted = 1;
    if ( (N_fit_predicted_upper - N_fit_predicted) >= (N_fit_predicted - N_fit_predicted_lower) ) {
      sigma_N_fit_predicted = (N_fit_predicted_upper - N_fit_predicted);
    }
    else {
      sigma_N_fit_predicted = (N_fit_predicted - N_fit_predicted_lower);
    }

    
    double Pull = (N_fit_predicted - N_MC_truth)/sigma_N_fit_predicted;

    h_Pull->Fill(Pull);
    
    delete h_bgr_toy;
    h_bgr_toy = NULL;
    delete h_sigbgr_toy;
    h_sigbgr_toy = NULL;
    
//    cout << Get_TestStatistic(h_data,h_bgr,h_sig) << endl;
//    cout << Get_TestStatistic(h_bgr_toy_unscaled,h_bgr,h_sig) << endl;
//    cout << Get_TestStatistic(h_sigbgr_toy_unscaled,h_bgr,h_sig) << endl;
//    cout << endl;
//
//    double alpha_bgr = R->Gaus(alpha_bgr_nom, alpha_bgr_sig);
//
//    //cout << alpha_bgr << endl;
//    
//    TH1D *h_bgr_scaled = (TH1D*) h_bgr->Clone("h_bgr");
//    TH1D *h_sigbgr_scaled = (TH1D*) h_sig->Clone("h_sig");
//    h_bgr_scaled->Scale(alpha_bgr);
//    h_sigbgr_scaled->Add(h_bgr_scaled);
//
//    TH1D *h_bgr_toy_scaled = GenerateToyDataSet(h_bgr_scaled);
//    TH1D *h_sigbgr_toy_scaled = GenerateToyDataSet(h_sigbgr_scaled);
//                        
//    h_data_scaled_TS->Fill(Get_TestStatistic(h_data,h_bgr_scaled,h_sig));
//    h_bgr_scaled_TS->Fill(Get_TestStatistic(h_bgr_toy_scaled,h_bgr_scaled,h_sig));
//    h_sigbgr_scaled_TS->Fill(Get_TestStatistic(h_sigbgr_toy_scaled,h_bgr_scaled,h_sig));
//
//    //cout << Get_TestStatistic(h_data,h_bgr_scaled,h_sig) << endl;
//    //cout << Get_TestStatistic(h_bgr_toy_scaled,h_bgr_scaled,h_sig) << endl;
//    //cout << Get_TestStatistic(h_sigbgr_toy_scaled,h_bgr_scaled,h_sig) << endl;
//    //cout << endl;

  }

  vector<double> sigbgr_t_quant = Get_Quantiles(h_sigbgr_TS);
  vector<double> bgr_t_quant = Get_Quantiles(h_bgr_TS);

  double t_sigbgr = sigbgr_t_quant[2];
  double t_bgr = bgr_t_quant[2];
  double t_data = Get_TestStatistic(h_data,h_bgr,h_sig);


  cout << t_sigbgr << endl;
  cout << t_bgr << endl;
  cout << t_data << endl;
  
  for (int i= 1; i <= h_unitGaus->GetNbinsX(); i++) {
    h_unitGaus->AddBinContent(i, TMath::Gaus(h_unitGaus->GetBinCenter(i),0,1, true));
  }

  h_unitGaus->Scale(NToys*h_unitGaus->GetBinWidth(1));

  
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetGrid();
  
  h_bgr_TS->SetLineWidth(2);
  h_sigbgr_TS->SetLineWidth(2);
  h_bgr_TS->SetLineColor(2);
  h_sigbgr_TS->SetLineColor(3);
  h_bgr_TS->SetStats(false);
  h_sigbgr_TS->SetStats(false);

  h_bgr_TS->Draw();
  h_sigbgr_TS->Draw("same");
  
  //-- axes
  AddText( 0.900, 0.035, "Test statistics t",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.2f ",h_bgr_TS->GetBinWidth(1)),0.060,90.,"right");   // Y-axis    
  
  auto legend1 = new TLegend(0.5,0.5,0.7,0.7);
  legend1->AddEntry(h_bgr_TS,"b-only");
  legend1->AddEntry(h_sigbgr_TS, "s + b");
  legend1->Draw();
  
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
  canvas2->SetGrid();
  
  h_Pull->SetLineWidth(2);
  h_unitGaus->SetLineWidth(2);
  
  h_Pull->SetLineColor(4);
  h_unitGaus->SetLineColor(3);

  h_Pull->Draw();
  h_unitGaus->Draw("l same");
  
  //-- axes
  AddText( 0.900, 0.035, "Pull",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.2f ",h_Pull->GetBinWidth(1)) ,0.060,90.,"right");   // Y-axis    

  auto legend2 = new TLegend(0.5,0.5,0.7,0.7);
  legend2->AddEntry(h_Pull,"Pull distribution");
  legend2->AddEntry(h_unitGaus, "Unit Gaussian");
  legend2->Draw();

//  h_data_scaled_TS->Draw();
//  h_bgr_scaled_TS->Draw("same");
//  h_sigbgr_scaled_TS->Draw("same");

  //double alpha_fit_center = get<0>(alpha_fit);
  //double alpha_fit_lower = get<1>(alpha_fit);
  //double alpha_fit_upper = get<2>(alpha_fit);

  //double Higgs_mass = 125.0;
  //double masswindow_width = 7.15;
  //double Lumi_scalefactor = 1.0;

  ////-- Get histograms from the files (higgs, zz and data)
  //TH1D *h_sig, *h_bgr, *h_data;
  //h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  //h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  //h_data = GetMassDistribution(2, Lumi_scalefactor);

  //double masswindow_lower = Higgs_mass-(masswindow_width/2);
  //double masswindow_upper = Higgs_mass+(masswindow_width/2);

  //int lower_bin_sig = h_sig->FindBin(masswindow_lower);
  //int upper_bin_sig = h_sig->FindBin(masswindow_upper);

  //double Nsig_win = h_sig->Integral(lower_bin_sig, upper_bin_sig);
  //double Nbgr_win = h_bgr->Integral(lower_bin_bgr, upper_bin_bgr);
  //double Ndata_win = h_data->Integral(lower_bin_data, upper_bin_data);

  //double Nbgr_est_win_center = alpha_fit_center*Nbgr_win;
  //double Nbgr_est_win_lower =  alpha_fit_lower*Nbgr_win;
  //double Nbgr_est_win_upper =  alpha_fit_upper*Nbgr_win;
  //
  //double Nbgr_est_win_uncert_plus = Nbgr_est_win_upper-Nbgr_est_win_center;
  //double Nbgr_est_win_uncert_minus = Nbgr_est_win_center-Nbgr_est_win_lower;

  //printf(" Background estimates: \n");
  //printf("   Higgs mass:  %.1f GeV\n", Higgs_mass);
  //printf("   Mass Window width:  %.2f GeV\n", masswindow_width);
  //printf("   Luminosity scale factor:  %.1f \n", Lumi_scalefactor);
  //printf("   \n");
  //printf("   Expected background (unscaled):  %5.2f \n", Nbgr_win);
  //printf("   Estimated expected background (scaled):  %5.2f   [%.2f,%.2f] \n", Nbgr_est_win_center, Nbgr_est_win_lower, Nbgr_est_win_upper);
  //printf("      Uncertanity:  + %.2f \n", Nbgr_est_win_uncert_plus);
  //printf("                    - %.2f \n", Nbgr_est_win_uncert_minus);
  //printf("   Expected signal:  %5.2f \n", Nsig_win);
  //printf("   Observed data:  %5d \n", int(Ndata_win));

  return;
}
