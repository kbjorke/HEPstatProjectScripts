#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

void ToyTest(){

  double Lumi_scalefactor = 3.0;
  int Irebin = 20;

  //-- Get histograms from the files (higgs, zz and data)
  TH1D *h_sig, *h_bgr, *h_data;
  h_sig  = GetMassDistribution(125, Lumi_scalefactor);
  h_bgr  = GetMassDistribution(1, Lumi_scalefactor);
  h_data = GetMassDistribution(2, Lumi_scalefactor);

  //-- Rebin histograms (only for plotting)
  h_sig->Rebin(Irebin);
  h_bgr->Rebin(Irebin);
  h_data->Rebin(Irebin);

  TH1D *h_sigbgr = (TH1D*) h_sig->Clone("h_sig");
  h_sigbgr->Add(h_bgr);

  TH1D *h_bgr_toy = GenerateToyDataSet(h_bgr);
  TH1D *h_sigbgr_toy = GenerateToyDataSet(h_sigbgr);

  //-- Prepare canvas and plot histograms
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
  canvas1->cd();

  //-- prepare histograms and plot them on canvas
  h_sigbgr_toy->SetFillColor(7);
  h_sigbgr_toy->SetAxisRange(0.,60,"Y");
  h_sigbgr_toy->SetAxisRange(0.,400.,"X");
  h_bgr_toy->SetFillColor(2);
  h_sigbgr_toy->Draw("hist");
  //h_bgr_toy->Draw("same");
  //h_bgr_toy->Draw("axis same");
  h_sigbgr->SetLineColor(2);
  h_sigbgr->SetLineWidth(2);
  h_sigbgr->Draw("l same");

  //-- some nice axes and add legend
  AddText( 0.900, 0.035, "4-lepton invariant mass [GeV]",0.060, 0.,"right");                             // X-axis
  AddText( 0.040, 0.900, Form("Number of events / %3.1f GeV",h_bgr_toy->GetBinWidth(1)) ,0.060,90.,"right"); // Y-axis
  AddText( 0.225, 0.825, Form("Luminosity scalefactor = %5.1f",Lumi_scalefactor),0.050, 0.,"left");
  TLegend *leg1 = new TLegend(0.65,0.65,0.90,0.75);
  leg1->SetBorderSize(0); leg1->SetFillColor(0);
  TLegendEntry *leg1a = leg1->AddEntry(h_sigbgr_toy, "Pseudo data set" , "f");  leg1a->SetTextSize(0.04);
  TLegendEntry *leg1b = leg1->AddEntry(h_sigbgr, "Expected s+b" , "l");
  leg1->Draw();


  return;
}
