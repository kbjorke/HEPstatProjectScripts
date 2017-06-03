#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void SignalScaleScan(){

  int NToys = 1e5;
  double Lumi_scalefactor = 1.0;
  int Irebin = 10;

  TH1D *h_SignalScaleScan_spb = new TH1D("h_SignalScaleScan_spb", "", 50, 0.5, 6.25);
  TH1D *h_SignalScaleScan_b = new TH1D("h_SignalScaleScan_b", "", 50, 0.5, 6.25);
  TH1D *h_SignalScaleScan_obs = new TH1D("h_SignalScaleScan_obs", "", 50, 0.5, 6.25);
  
  int Nbins = 400;

  for (int i_bin = 1; i_bin<=h_SignalScaleScan_spb->GetNbinsX(); i_bin++ ){
    cout << i_bin << endl;
  
    TH1D *h_bgr_TS = new TH1D("h_bgr_TS","", Nbins, -150, 150);
    TH1D *h_sigbgr_TS = new TH1D("h_sigbgr_TS","", Nbins, -150, 150);
  
    //-- Get histograms from the files (higgs, zz and data)
    TH1D *h_sig, *h_bgr, *h_data;
    h_sig  = GetMassDistribution(125, Lumi_scalefactor);
    h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
    h_data = GetMassDistribution(2, Lumi_scalefactor);
    
    h_sig->Rebin(Irebin);  
    h_bgr->Rebin(Irebin);  
    h_data->Rebin(Irebin); 

    /////// Scaled:
    
    double alpha_bgr_nom = 1.11;
    double alpha_bgr_sig = 0.07;
    
    h_bgr->Scale(alpha_bgr_nom);

    /////////////////
    
    h_sig->Scale(h_SignalScaleScan_spb->GetBinCenter(i_bin));
    
    TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
    h_sigbgr->Add(h_bgr);

    for (int toy = 1; toy <= NToys; toy++) {
      TRandom3 *R = new TRandom3(0);

      //-- Create new histogram for the data-set
      TH1D *h_bgr_toy = (TH1D*) h_bgr->Clone("h_mass_toydataset"); h_bgr_toy->Reset(); 
      TH1D *h_sigbgr_toy = (TH1D*) h_sigbgr->Clone("h_mass_toydataset"); h_sigbgr_toy->Reset(); 

      //-- Loop over bins and draw Poisson number of event in each bin
      for (int i_bin = 1; i_bin <= h_bgr->GetNbinsX(); i_bin++){
        double toy_bgr = R->Poisson(h_bgr->GetBinContent(i_bin));
        double toy_sig = R->Poisson(h_sig->GetBinContent(i_bin));
        h_bgr_toy->SetBinContent(i_bin, toy_bgr);
        h_sigbgr_toy->SetBinContent(i_bin, toy_bgr+toy_sig);
      }


      h_bgr_TS->Fill(Get_TestStatistic(h_bgr_toy,h_bgr,h_sig));
      h_sigbgr_TS->Fill(Get_TestStatistic(h_sigbgr_toy,h_bgr,h_sig));

      delete h_bgr_toy;
      h_bgr_toy = NULL;
      delete h_sigbgr_toy;
      h_sigbgr_toy = NULL;

      delete R;
      R = NULL;
    }

    vector<double> sigbgr_t_quant = Get_Quantiles(h_sigbgr_TS);
    vector<double> bgr_t_quant = Get_Quantiles(h_bgr_TS);

    double t_sigbgr_median = sigbgr_t_quant[2];
    double t_bgr_median = bgr_t_quant[2];
    double t_data = Get_TestStatistic(h_data,h_bgr,h_sig);
    
    TH1D *h_bgr_TS_norm = (TH1D*) h_bgr_TS->Clone();
    TH1D *h_sigbgr_TS_norm = (TH1D*) h_sigbgr_TS->Clone();

    h_bgr_TS_norm->Scale(1/h_bgr_TS_norm->Integral());
    h_sigbgr_TS_norm->Scale(1/h_sigbgr_TS_norm->Integral());
    
    /// CL_s+b ///
    
    double CLspb_sb = IntegrateFromRight(h_sigbgr_TS_norm, t_sigbgr_median);
    double CLspb_b = IntegrateFromRight(h_sigbgr_TS_norm, t_bgr_median);
    double CLspb_n = IntegrateFromRight(h_sigbgr_TS_norm, t_data);

    h_SignalScaleScan_spb->SetBinContent(i_bin,CLspb_sb);
    h_SignalScaleScan_b->SetBinContent(i_bin,CLspb_b);
    h_SignalScaleScan_obs->SetBinContent(i_bin,CLspb_n);

    delete h_bgr_TS;
    h_bgr_TS = NULL;
    delete h_sigbgr_TS;
    h_sigbgr_TS = NULL;
    delete h_bgr_TS_norm;
    h_bgr_TS_norm = NULL;
    delete h_sigbgr_TS_norm;
    h_sigbgr_TS_norm = NULL;
  }
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
  canvas2->SetGrid();

  h_SignalScaleScan_spb->SetLineWidth(2);
  h_SignalScaleScan_b->SetLineWidth(2);
  h_SignalScaleScan_obs->SetLineWidth(2);
  h_SignalScaleScan_spb->SetLineColor(2);
  h_SignalScaleScan_b->SetLineColor(3);
  h_SignalScaleScan_obs->SetLineColor(4);

  THStack* stack = new THStack("mystack1", "");
  stack->Add(h_SignalScaleScan_spb);
  stack->Add(h_SignalScaleScan_b);
  stack->Add(h_SignalScaleScan_obs);
  stack->Draw("l nostack");
  stack->SetMinimum(0);
  stack->SetMaximum(1);
  stack->Draw("l nostack");

  //canvas2->SetLogy(1);

  TLine *line1 = new TLine(0.5,0.05,6.0,0.05);
  line1->SetLineStyle(9);
  line1->SetLinewidth(3);
  line1->Draw();

  //-- axes
  AddText( 0.900, 0.035, "Cross-section scale factor",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "CL_s+b",0.060,90.,"right");   // Y-axis   
  AddText( 0.890, 0.170, "exclusion (CL_s+b = 0.05)", 0.040, 0., "right");


  auto legend = new TLegend(0.5,0.5,0.7,0.7);
  legend->AddEntry(h_SignalScaleScan_spb,"median s+b");
  legend->AddEntry(h_SignalScaleScan_b,"median b-only");
  legend->AddEntry(h_SignalScaleScan_obs,"observed");
  legend->Draw();

  return;
}
