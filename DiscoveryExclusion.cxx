#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void DiscoveryExclusion(){

  int NToys = 5e6;
  double Lumi_scalefactor = 6.0;

  int Irebin = 40;
  
  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);
  
  int Nbins = 400;
  
  TH1D *h_bgr_TS = new TH1D("h_bgr_TS","", Nbins, -150, 150);
  TH1D *h_sigbgr_TS = new TH1D("h_sigbgr_TS","", Nbins, -150, 150);

  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin); 

  /////// Scaled:
  
  double alpha_bgr_nom = 1.11;
  double alpha_bgr_sig = 0.07;
  
  h_bgr->Scale(alpha_bgr_nom);

  /////////////////
  
  TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
  h_sigbgr->Add(h_bgr);

  for (int toy = 1; toy <= NToys; toy++) {
    cout << toy << endl;
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
  }

  vector<double> sigbgr_t_quant = Get_Quantiles(h_sigbgr_TS);
  vector<double> bgr_t_quant = Get_Quantiles(h_bgr_TS);

  double t_sigbgr_median = sigbgr_t_quant[2];
  double t_bgr_median = bgr_t_quant[2];
  double t_data = Get_TestStatistic(h_data,h_bgr,h_sig);


  cout << t_sigbgr_median << endl;
  cout << t_bgr_median << endl;
  cout << t_data << endl;
  

  
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
  
  
  TH1D *h_bgr_TS_norm = (TH1D*) h_bgr_TS->Clone();
  TH1D *h_sigbgr_TS_norm = (TH1D*) h_sigbgr_TS->Clone();


  h_bgr_TS_norm->Scale(1/h_bgr_TS_norm->Integral());
  h_sigbgr_TS_norm->Scale(1/h_sigbgr_TS_norm->Integral());
  
  printf("Number of toys:  %6.1e \n", double (NToys));
  printf("Luminosity scale factor:  %.1f \n", Lumi_scalefactor);
  printf("\n");

  /// 1 - CL_b ///
  
  double p1mCLb_sb = 1 - IntegrateFromRight(h_bgr_TS_norm, t_sigbgr_median);
  double p1mCLb_b = 1 - IntegrateFromRight(h_bgr_TS_norm, t_bgr_median);
  double p1mCLb_n = 1 - IntegrateFromRight(h_bgr_TS_norm, t_data);
  
  printf(" 1 - CL_b values: \n");
  printf("     median s+b experiment:     1 - CL_b = %8.6f (%5.3f sigma) \n", p1mCLb_sb, ROOT::Math::gaussian_quantile_c(p1mCLb_sb,1));
  printf("     median b-only experiment:  1 - CL_b = %8.6f (%5.3f sigma) \n", p1mCLb_b, ROOT::Math::gaussian_quantile_c(p1mCLb_b,1));
  printf("     data                       1 - CL_b = %8.6f (%5.3f sigma) \n", p1mCLb_n, ROOT::Math::gaussian_quantile_c(p1mCLb_n,1));

  /// CL_s+b ///
  
  double CLspb_sb = IntegrateFromRight(h_sigbgr_TS_norm, t_sigbgr_median);
  double CLspb_b = IntegrateFromRight(h_sigbgr_TS_norm, t_bgr_median);
  double CLspb_n = IntegrateFromRight(h_sigbgr_TS_norm, t_data);
  
  printf("   CL_s+b values: \n");
  printf("     median s+b experiment:       CL_s+b = %8.6f \n", CLspb_sb);
  printf("     median b-only experiment:    CL_s+b = %8.6f \n", CLspb_b);
  printf("     data                         CL_s+b = %8.6f \n", CLspb_n);
  
  
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
  canvas2->SetGrid();
  
  h_bgr_TS_norm->Draw();
  h_sigbgr_TS_norm->Draw("same");

  return;
}
