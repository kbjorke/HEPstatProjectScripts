#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

TH1D *h_MassScan_exp_sig = new TH1D("h_MassScan_exp_sig", "", 600, 40.0, 400.0);
TH1D *h_MassScan_exp_win = new TH1D("h_MassScan_exp_win", "", 600, 40.0, 400.0);

TH1D *h_MassScan_obs_sig = new TH1D("h_MassScan_obs_sig", "", 600, 40.0, 400.0);
TH1D *h_MassScan_obs_win = new TH1D("h_MassScan_obs_win", "", 600, 40.0, 400.0);

//double masswindow_width = 7.15;
//double rel_sigma_b = 0.10;

void mass_scan(){

  for (int i_bin = 1; i_bin<=h_MassScan_exp_sig->GetNbinsX(); i_bin++ ){
    double Mass = h_MassScan_exp_sig->GetBinLowEdge(i_bin);
    cout << Mass << endl;
    tuple<double,double,double,double> results = Significance_Optimization(1.0, Mass, false);
    double significance_exp = get<0>(results);
    double optimal_mass_window_width_exp = get<1>(results);
    double significance_obs = get<2>(results);
    double optimal_mass_window_width_obs = get<3>(results);
    h_MassScan_exp_sig->SetBinContent(i_bin, significance_exp);
    h_MassScan_exp_win->SetBinContent(i_bin, optimal_mass_window_width_exp);
    h_MassScan_obs_sig->SetBinContent(i_bin, significance_obs);
    h_MassScan_obs_win->SetBinContent(i_bin, optimal_mass_window_width_obs);
  }

  //----------------------------------
  //-- [3] Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
  canvas1->SetGrid();
  canvas1->cd();

  gStyle->SetOptStat(0);

  h_MassScan_exp_sig->SetLineColor(1);
  h_MassScan_exp_sig->SetLineWidth(2);
  //h_masswindow_observed->SetLineColor(4);
  //h_masswindow_observed->SetLineWidth(2);

  h_MassScan_exp_sig->SetAxisRange(-1.,7.,"Y");
  h_MassScan_exp_sig->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Mass GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);
  
//-- prepare gif
  canvas1->Print(Form("./Significance_massscalefactor-scan.gif"));
  
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
  canvas2->SetLeftMargin(0.125);
  canvas2->SetBottomMargin(0.125);
  canvas2->SetGrid();
  canvas2->cd();

  h_MassScan_exp_win->SetLineColor(1);
  h_MassScan_exp_win->SetLineWidth(2);
  //h_masswin_expdow_observed->SetLineColor(4);
  //h_masswin_expdow_observed->SetLineWidth(2);

  h_MassScan_exp_win->SetAxisRange(0.,25.,"Y");
  h_MassScan_exp_win->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Mass GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Mass window width GeV" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Optimal mass window width GeV",0.050, 0.,"right",1,1);
  
  TCanvas * canvas3 = new TCanvas( "canvas3","Standard Canvas",600,400);
  canvas3->SetLeftMargin(0.125);
  canvas3->SetBottomMargin(0.125);
  canvas3->SetGrid();
  canvas3->cd();

  h_MassScan_obs_sig->SetLineColor(1);
  h_MassScan_obs_sig->SetLineWidth(2);
  //h_massobs_sigdow_observed->SetLineColor(4);
  //h_massobs_sigdow_observed->SetLineWidth(2);

  h_MassScan_obs_sig->SetAxisRange(-1.,7.,"Y");
  h_MassScan_obs_sig->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Mass GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Observed significance",0.050, 0.,"right",1,1);
  
  TCanvas * canvas4 = new TCanvas( "canvas4","Standard Canvas",600,400);
  canvas4->SetLeftMargin(0.125);
  canvas4->SetBottomMargin(0.125);
  canvas4->SetGrid();
  canvas4->cd();

  h_MassScan_obs_win->SetLineColor(1);
  h_MassScan_obs_win->SetLineWidth(2);
  //h_massobs_window_observed->SetLineColor(4);
  //h_massobs_window_observed->SetLineWidth(2);

  h_MassScan_obs_win->SetAxisRange(0.,25.,"Y");
  h_MassScan_obs_win->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Mass GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Mass window width GeV" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Optimal mass winddow width GeV",0.050, 0.,"right",1,1);


  THStack* stack = new THStack("mystack", "");
  stack->Add(h_MassScan_exp_sig);
  stack->Add(h_MassScan_obs_sig);

  h_MassScan_exp_sig->SetLineColor(2);
  h_MassScan_exp_sig->SetLineWidth(2);

  h_MassScan_obs_sig->SetLineColor(3);
  h_MassScan_obs_sig->SetLineWidth(2); 

  stack->Draw("l nostack");
  stack->SetMinimum(-0.5);
  stack->SetMaximum(5);
  stack->Draw("l nostack");
  
//-- axes
  AddText( 0.900, 0.035, "Higgs Mass GeV",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    

  auto legend = new TLegend(0.5,0.5,0.7,0.7);
  legend->AddEntry(h_MassScan_exp_sig,"Expected Significance");
  legend->AddEntry(h_MassScan_obs_sig, "Observed Significance");
  legend->Draw();


  return;
}
