#include "TString.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMarker.h"
#include "TColor.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TGraph.h"

#include <iostream>

using namespace std;

void compare2Histos(TString histo1, TString histo2, const char* Filename="compare2Histos.root", const char* OutputPath="/home/gu74yub/SimAna/Results") {
  
  TFile *File1 = TFile::Open(histo1.Data());
  TFile *File2  = TFile::Open(histo2.Data());

  TGraph *Fraction = new TGraph();

  if (!File1) {
    cout << "Histo 1 Not Open" << endl;
    return;
  }

  if (!File2) {
    cout << "Histo 2 Not Open" << endl;
    return;
  }

  TH1F *Histo1             = (TH1F*)(File1->FindObjectAny("hNElectronsPerCM"));
  TH1F *Histo2             = (TH1F*)(File2->FindObjectAny("NElePerCM"));


  for(int i=0 ; i<75; ++i) {
  	if (Histo2->GetBinContent(i) == 0 || Histo1->GetBinContent(i) == 0) {
      Fraction->SetPoint(i,i,0);
  	}
  	else {
      Fraction->SetPoint(i,i,float(Histo1->GetBinContent(i))/float(Histo2->GetBinContent(i)));
    }
  }

  short Color1  = kPink;
  short Color2  = kBlue-2;

  TCanvas *cHistos = new TCanvas();
  float maxHisto1  = Histo1->GetBinContent(Histo1->GetMaximumBin());
  float maxHisto2  = Histo2->GetBinContent(Histo2->GetMaximumBin());
  float maxHisto   = max(maxHisto1,maxHisto2);
  Histo1->GetYaxis()->SetRangeUser(0,1.2*maxHisto);
  Histo1->GetXaxis()->SetRangeUser(0,150);
  Histo1->SetLineColor(Color1);
  Histo1->SetMarkerStyle(20);
  Histo1->SetMarkerColor(Color1);
  Histo1->GetYaxis()->SetTitle("normalized counts");
  Histo1->Draw("hist");
  Histo2->SetLineColor(Color2);
  Histo2->SetMarkerStyle(20);
  Histo2->SetMarkerColor(Color2);
  Histo2->Draw("same, hist");
  TLegend *leg1 = new TLegend(0.75,0.75,0.9,0.9);
  leg1->AddEntry(Histo1,"O2","p");
  leg1->AddEntry(Histo2,"AliRoot","p");
  leg1->SetTextSize(0.05);
  leg1->SetBorderSize(0);
  leg1->Draw("same");

  TCanvas *cFrac = new TCanvas();
  Fraction->SetMarkerStyle(20);
  Fraction->GetYaxis()->SetTitle("O2/AliRoot");
  Fraction->GetXaxis()->SetTitle("bin number");
  Fraction->Draw("ap");

  TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");
  OutFile->WriteObject(cHistos, "cComparison");
  OutFile->WriteObject(cFrac, "cFraction");




}