#include <iostream>
#include <TH1D.h>
#include <TCanvas.h>

#include "Walkthrough_skeleton.C"

using namespace std;

TH1D *h_LumiScan = new TH1D("h_LumiScan", "", 72, 1.0, 8.2);
TH1D *h_LumiScan_sigma = new TH1D("h_LumiScan_sigma", "", 72, 1.0, 8.2);

double masswindow_width = 7.15;
double rel_sigma_b = 0.10;

void lumi_scan(){

  for (int i_bin = 1; i_bin<=h_LumiScan->GetNbinsX(); i_bin++ ){
    double Lumi = h_LumiScan->GetBinLowEdge(i_bin);
    cout << Lumi << endl;
    double significance = Significance(Lumi, masswindow_width, 0, 125);
    double significance_sigma = Significance(Lumi, masswindow_width, rel_sigma_b, 125);
    cout << significance << endl;
    h_LumiScan->SetBinContent(i_bin, significance);
    h_LumiScan_sigma->SetBinContent(i_bin, significance_sigma);
  }

  //----------------------------------
  //-- [3] Plot histogram and make gif
  //----------------------------------
  TCanvas * canvas1 = new TCanvas( "canvas1","Standard Canvas",600,400);
  canvas1->SetLeftMargin(0.125);
  canvas1->SetBottomMargin(0.125);
  canvas1->SetGrid();
  canvas1->cd();

  THStack* stack = new THStack("mystack", "");
  stack->Add(h_LumiScan);
  stack->Add(h_LumiScan_sigma);

  gStyle->SetOptStat(0);

  h_LumiScan->SetLineColor(2);
  h_LumiScan->SetLineWidth(2);
  
  h_LumiScan_sigma->SetLineColor(3);
  h_LumiScan_sigma->SetLineWidth(2);

  h_LumiScan->SetAxisRange(2.,6.5,"Y");
  h_LumiScan_sigma->SetAxisRange(2.,6.5,"Y");
  stack->Draw("l nostack");
  stack->SetMinimum(2);
  stack->SetMaximum(6.5);
  stack->Draw("l nostack");

  //-- axes
  AddText( 0.900, 0.035, "Luminosity scalefactor",0.060, 0.,"right"); // X-axis
  AddText( 0.040, 0.900, "Significance" ,0.060,90.,"right");   // Y-axis    
  AddText( 0.225, 0.825, Form("Mass window width = %2.2f", masswindow_width),0.050, 0.,"left");

  AddText( 0.700, 0.200, "Expected significance",0.050, 0.,"right",1,1);

  auto legend = new TLegend(0.5,0.5,0.7,0.7);
  legend->AddEntry(h_LumiScan,"No uncertainty");
  legend->AddEntry(h_LumiScan_sigma, "sigma b = 10%");
  legend->Draw();
  
//-- prepare gif
  canvas1->Print(Form("./Significance_lumiscalefactor-scan.gif"));

  return;
}
