#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <vector>

#include <fstream>

#include <iostream>



#include "TROOT.h"

#include "TLine.h"

#include "TMath.h"

#include "TFile.h"

#include "TTree.h"

#include "TClonesArray.h"

#include "TCanvas.h"

#include "TStyle.h"

#include "TGraph.h"

#include "TAxis.h"

#include "TH1F.h"

#include "TPCSimulation/Point.h"

#include "TPCBase/Digit.h"

#include "TPCSimulation/Digitizer.h"

#include "DataFormatsTPC/TrackTPC.h"

#include "ReconstructionDataFormats/Track.h"

#include "DataFormatsTPC/Cluster.h"

#include "TPCBase/Mapper.h"

#endif



using namespace o2::TPC;



void drawSectorBoundaries();



void drawHits(std::string simFile="", int iEv = -1, int sector =0)

{

  gStyle->SetMarkerStyle(20);

  gStyle->SetMarkerSize(0.5);

  gStyle->SetTitleSize(24);



  // ===| process the hits |====================================================

  TFile *hitFile   = TFile::Open(simFile.data());

  TTree *hitTree = (TTree *)gDirectory->Get("o2sim");



  std::vector<o2::TPC::HitGroup> *sectorHitsArray[Sector::MAXSECTOR];

  for (int s=0;s<Sector::MAXSECTOR;++s){

    if(sector != s) continue;

    sectorHitsArray[s] = nullptr;

    std::stringstream sectornamestr;

    sectornamestr << "TPCHitsSector" << s;

    hitTree->SetBranchAddress(sectornamestr.str().c_str(), &sectorHitsArray[s]);

  }



  TGraph *grHitsA = new TGraph();

  grHitsA->SetTitle(Form("Hits A-Side Event %d;x (cm);y (cm)", iEv));

  grHitsA->SetMarkerColor(kBlue+2);

  grHitsA->SetPoint(0,-1000,-1000);



  TGraph *grHitsC = new TGraph();

  grHitsC->SetTitle(Form("Hits C-Side Event %d;x (cm);y (cm)", iEv));

  grHitsC->SetMarkerColor(kBlue+2);

  grHitsC->SetPoint(0,-1000,-1000);



  TGraph *grHitsAzr = new TGraph();

  grHitsAzr->SetTitle(Form("Hits comparison A-Side Event %d;z (cm);r (cm)", iEv));

  grHitsAzr->SetMarkerColor(kBlue+2);

  grHitsAzr->SetPoint(0,-1000,-1000);



  TGraph *grHitsCzr = new TGraph();

  grHitsCzr->SetTitle(Form("Hits comparison C-Side Event %d;z (cm);r (cm)", iEv));

  grHitsCzr->SetMarkerColor(kBlue+2);

  grHitsCzr->SetPoint(0,-1000,-1000);





  int hitCounterA = 1;

  int hitCounterC = 1;

  for(int event=0; event<hitTree->GetEntriesFast(); ++event) {

    if(iEv != -1 && event != iEv) continue;

    hitTree->GetEntry(event);

    for(auto& inputgroup : *sectorHitsArray[sector]) {

      const int MCTrackID = inputgroup.GetTrackID();

      for(size_t hitindex = 0; hitindex < inputgroup.getSize(); ++hitindex){

	const auto& eh = inputgroup.getHit(hitindex);



	// A side

	if(eh.GetZ() > 0 ) {

	  grHitsA->SetPoint(hitCounterA, eh.GetX(), eh.GetY());

	  grHitsAzr->SetPoint(hitCounterA++, eh.GetZ(), TMath::Sqrt(eh.GetX()*eh.GetX() + eh.GetY()*eh.GetY()));

	}

	// C side

	if(eh.GetZ() < 0 ) {

	  grHitsC->SetPoint(hitCounterC, eh.GetX(), eh.GetY());

	  grHitsCzr->SetPoint(hitCounterC++, eh.GetZ(), TMath::Sqrt(eh.GetX()*eh.GetX() + eh.GetY()*eh.GetY()));

	}

      }

    }

  }

  



  // |=============================================================

  auto CDigits = new TCanvas("CDigits", "Compare Digits - Hits on A & C side", 1200, 600);

  CDigits->Divide(2,1);

  CDigits->cd(1);

  grHitsA->Draw("ap");

  grHitsA->GetXaxis()->SetLimits(-250, 250);

  grHitsA->SetMinimum(-250);

  grHitsA->SetMaximum(250);

  drawSectorBoundaries();

  

  CDigits->cd(2);

  grHitsC->Draw("ap");

  grHitsC->GetXaxis()->SetLimits(-250, 250);

  grHitsC->SetMinimum(-250);

  grHitsC->SetMaximum(250);

  drawSectorBoundaries();

  

  auto CDigitsXZ = new TCanvas("CDigitsXZ", "Compare Digits - Hits on A & C side", 600, 600);

  grHitsAzr->Draw("ap");

  grHitsAzr->GetXaxis()->SetLimits(-250, 250);

  grHitsAzr->SetMinimum(-250);

  grHitsAzr->SetMaximum(250);

  grHitsCzr->Draw("p");

  grHitsCzr->GetXaxis()->SetLimits(-250, 250);

  grHitsCzr->SetMinimum(-250);

  grHitsCzr->SetMaximum(250);



}



void drawSectorBoundaries()

{

  TLine *sectorBoundary[18];

  for(float i = 0; i<18; ++i) {

    const float angle = i*20.f*TMath::DegToRad();

    sectorBoundary[int(i)] = new TLine(80.f*std::cos(angle), 80.f*std::sin(angle), 250.f*std::cos(angle), 250.f*std::sin(angle));

    sectorBoundary[int(i)]->SetLineStyle(2);

    sectorBoundary[int(i)]->SetLineColor(kGray);

    sectorBoundary[int(i)]->Draw("same");

  }

}
