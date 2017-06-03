#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void LuminosityScan(){

  int NToys = 1e6;
  
  int Irebin = 20;
    
  TH1D *h_LumiScan_spb = new TH1D("h_LumiScan_spb", "", 7, 1.0, 8.0);
  TH1D *h_LumiScan_b = new TH1D("h_LumiScan_b", "", 7, 1.0, 8.0);
  
  int Nbins = 400;

  for (int i_bin = 1; i_bin<=h_LumiScan_spb->GetNbinsX(); i_bin++ ){
    double Lumi_scalefactor = h_LumiScan_spb->GetBinCenter(i_bin);
  
    //-- Get histograms from the files (higgs and zz)
    TH1D *h_sig, *h_bgr;
    h_sig  = GetMassDistribution(125, Lumi_scalefactor);
    h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
    
    h_sig->Rebin(Irebin);
    h_bgr->Rebin(Irebin);
    
    TH1D *h_bgr_TS = new TH1D("h_bgr_TS","", Nbins, -150, 150);
    TH1D *h_sigbgr_TS = new TH1D("h_sigbgr_TS","", Nbins, -150, 150);

    /////// Scaled:
    
    double alpha_bgr_nom = 1.11;
    double alpha_bgr_sig = 0.07;
    
    h_bgr->Scale(alpha_bgr_nom);

    /////////////////
    
    TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
    h_sigbgr->Add(h_bgr);

    cout << Lumi_scalefactor << endl;
    cout << int(Lumi_scalefactor*NToys) << endl;

    for (int toy = 1; toy <= int(Lumi_scalefactor*NToys); toy++) {
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

    TH1D *h_bgr_TS_norm = (TH1D*) h_bgr_TS->Clone();
    TH1D *h_sigbgr_TS_norm = (TH1D*) h_sigbgr_TS->Clone();

    h_bgr_TS_norm->Scale(1/h_bgr_TS_norm->Integral());
    h_sigbgr_TS_norm->Scale(1/h_sigbgr_TS_norm->Integral());
    
    /// 1 - CL_b ///
    
    double p1mCLb_sb = 1 - IntegrateFromRight(h_bgr_TS_norm, t_sigbgr_median);
    double p1mCLb_b = 1 - IntegrateFromRight(h_bgr_TS_norm, t_bgr_median);

    h_LumiScan_spb->SetBinContent(i_bin,p1mCLb_sb);
    h_LumiScan_b->SetBinContent(i_bin,p1mCLb_b);

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
  
  h_LumiScan_spb->SetLineWidth(2);
  h_LumiScan_b->SetLineWidth(2);
  h_LumiScan_spb->SetLineColor(2);
  h_LumiScan_b->SetLineColor(3);

  THStack* stack = new THStack("mystack1", "");
  stack->Add(h_LumiScan_spb);
  stack->Add(h_LumiScan_b);
  stack->Draw("l nostack");
  stack->SetMinimum(ROOT::Math::normal_cdf_c(6));
  stack->SetMaximum(1);
  stack->Draw("l nostack");

  canvas2->SetLogy(1);

  TLine *line1 = new TLine(1,ROOT::Math::normal_cdf_c(3),7.5,ROOT::Math::normal_cdf_c(3));
  line1->SetLineWidth(3);
  line1->SetLineStyle(9);
  line1->Draw();

  TLine *line2 = new TLine(1,ROOT::Math::normal_cdf_c(4),7.5,ROOT::Math::normal_cdf_c(4));
  line2->SetLineWidth(3);
  line2->SetLineStyle(9);
  line2->Draw();

  TLine *line3 = new TLine(1,ROOT::Math::normal_cdf_c(5),7.5,ROOT::Math::normal_cdf_c(5));
  line3->SetLineWidth(3);
  line3->SetLineStyle(9);
  line3->Draw();

  //-- axes
  AddText( 0.900, 0.035, "Luminosity scale factor",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "1 -CL_b",0.060,90.,"right");   // Y-axis   
  AddText( 0.800, 0.700, "3 sigma", 0.040, 0., "right");
  AddText( 0.800, 0.530, "4 sigma", 0.040, 0., "right");
  AddText( 0.800, 0.350, "5 sigma", 0.040, 0., "right");


  auto legend = new TLegend(0.5,0.5,0.7,0.7);
  legend->AddEntry(h_LumiScan_spb,"median s+b");
  legend->AddEntry(h_LumiScan_b,"median b-only");
  legend->Draw();


  return;
}
