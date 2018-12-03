#include <iostream>
#include <math.h>

#include "TPCSimulation/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/TrackTPC.h"
#include "TPCBase/CRU.h"
#include "TPCBase/ROC.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/Painter.h"

#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TF1.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TPaveText.h"
#include "TStyle.h"

void ConvertGainMapToNiceCPad(TString GainMapFile)
{
  using namespace o2::TPC;

  gStyle->SetOptStat(0);
  Mapper& mapper = Mapper::instance();

  TFile f(GainMapFile);
//  gROOT->cd();

  TH1D *hGainDist = nullptr;
  f.GetObject("GainDistribution", hGainDist);

  TCanvas *GainDist = new TCanvas("c2","gain distribution");
  float hMean = hGainDist->GetMean();
  float hRMS = hGainDist->GetRMS();
  TPaveText *pave = new TPaveText(0.6,.7,.9,.9,"NDC");
  pave->SetBorderSize(1);
  pave->AddText(Form("mean: %.2f", hMean));
  pave->AddText(Form("RMS: %.2f", hRMS));
  hGainDist->Draw();
//  GainDist->Update();
  pave->Draw("same");
//  GainDist->Update();

  TH2D *hGainMap = nullptr;
  f.GetObject("GainMap", hGainMap);

  TH2D *hGainMapCPad = new TH2D("GainMap", "; Row; Pad", 63,0,62,101,-50,50);

  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<90; ++ipad){
      double gain = hGainMap->GetBinContent(irow,ipad);
      float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;
      hGainMapCPad->Fill(irow, cpad, gain);
    }
  }

  TCanvas *GainMap = new TCanvas("c1","gain map");
  hGainMapCPad->GetZaxis()->SetRangeUser(0.0,1.3);
  hGainMapCPad->GetYaxis()->SetTitleSize(24);
  hGainMapCPad->GetYaxis()->SetLabelSize(21);
  hGainMapCPad->GetYaxis()->SetTitleFont(43);
  hGainMapCPad->GetYaxis()->SetLabelFont(43);
  hGainMapCPad->GetYaxis()->SetTitleOffset(1.1);

  hGainMapCPad->GetXaxis()->SetTitleSize(24);
  hGainMapCPad->GetXaxis()->SetLabelSize(21);
  hGainMapCPad->GetXaxis()->SetTitleFont(43);
  hGainMapCPad->GetXaxis()->SetLabelFont(43);
  hGainMapCPad->GetXaxis()->SetTitleOffset(1.1);
  hGainMapCPad->Draw("colz");

  return;
}
