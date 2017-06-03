#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include <tuple>

#include "Walkthrough_skeleton.C"

using namespace std;

TH1D *h_LumiScan_sig = new TH1D("h_LumiScan_sig", "", 72, 1.0, 8.2);
TH1D *h_LumiScan_win = new TH1D("h_LumiScan_win", "", 72, 1.0, 8.2);

void lumi_scan_optimasswin(){

  for (int i_bin = 1; i_bin<=h_LumiScan_sig->GetNbinsX(); i_bin++ ){
    double Lumi = h_LumiScan_sig->GetBinLowEdge(i_bin);
    cout << Lumi << endl;
    tuple<double,double,double,double> results = Significance_Optimization(Lumi, 125, false);
    double significance = get<0>(results);
    double optimal_mass_window_width = get<1>(results);
    cout << significance << endl;
    h_LumiScan_sig->SetBinContent(i_bin, significance);
    h_LumiScan_win->SetBinContent(i_bin, optimal_mass_window_width);
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

  h_LumiScan_sig->SetLineColor(2);
  h_LumiScan_sig->SetLineWidth(2);

  h_LumiScan_sig->SetAxisRange(1.,7.,"Y");
  h_LumiScan_sig->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Luminosity scalefactor",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);
  
//-- prepare gif
  canvas1->Print(Form("./Significance_lumiscalefactor-scan.gif"));
  
  TCanvas * canvas2 = new TCanvas( "canvas2","Standard Canvas",600,400);
  canvas2->SetLeftMargin(0.125);
  canvas2->SetBottomMargin(0.125);
  canvas2->SetGrid();
  canvas2->cd();

  h_LumiScan_win->SetLineColor(2);
  h_LumiScan_win->SetLineWidth(2);

  h_LumiScan_win->SetAxisRange(0.,10.,"Y");
  h_LumiScan_win->Draw("l");

  //-- axes
  AddText( 0.900, 0.035, "Luminosity scalefactor",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Mass window width GeV" ,0.060,90.,"right");   // Y-axis    

  AddText( 0.700, 0.200, "Optimal mass window width GeV",0.050, 0.,"right",1,1);
  
//-- prepare gif
  canvas1->Print(Form("./Significance_lumiscalefactor-scan.gif"));

  return;
}
