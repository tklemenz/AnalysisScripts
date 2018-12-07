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
#include "TGraph2D.h"
#include "TAxis.h"
#include "TH1F.h"
#include "TPCSimulation/Point.h"
#include "TPCBase/Digit.h"
#include "TPCSimulation/Digitizer.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/CDBInterface.h"
#endif

using namespace o2::TPC;

void drawSectorBoundaries();

void testHits(std::string simFile="o2sim.root")
{
  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetTitleSize(24);

  auto& cdb = CDBInterface::instance();
  cdb.setUseDefaults();
  auto& sampa = SAMPAProcessing::instance();

  // ===| process the hits |====================================================
  TFile *hitFile   = TFile::Open(simFile.data());
  TTree *hitTree = (TTree *)gDirectory->Get("o2sim");

  std::vector<o2::TPC::HitGroup> *sectorHitsArray[Sector::MAXSECTOR];
  for (int s=0;s<Sector::MAXSECTOR;++s){
    sectorHitsArray[s] = nullptr;
    std::stringstream sectornamestr;
    sectornamestr << "TPCHitsShiftedSector" << s;
    hitTree->SetBranchAddress(sectornamestr.str().c_str(), &sectorHitsArray[s]);
  }

  TGraph *grHitsA = new TGraph();
  grHitsA->SetTitle("Hits - Cluster comparison A-Side;x (cm);y (cm)");
  grHitsA->SetMarkerColor(kBlue+2);

  TGraph *grHitsC = new TGraph();
  grHitsC->SetTitle("Hits - Cluster comparison C-Side;x (cm);y (cm)");
  grHitsC->SetMarkerColor(kBlue+2);

  TGraph *grHitsAzr = new TGraph();
  grHitsAzr->SetTitle("Hits - Cluster comparison A-Side;z (cm);r (cm)");
  grHitsAzr->SetMarkerColor(kBlue+2);

  TGraph *grHitsCzr = new TGraph();
  grHitsCzr->SetTitle("Hits - Cluster comparison C-Side Event %d;z (cm);r (cm)");
  grHitsCzr->SetMarkerColor(kBlue+2);

  TGraph2D *grHits2D = new TGraph2D();
  grHits2D->SetTitle("Hits - 3D;x (cm);y (cm);z (cm)");
  grHits2D->SetMarkerColor(kBlue+2);


  int hitCounterA   = 0;
  int hitCounterC   = 0;
  int hitCounter2D  = 0;

  for(int event=0; event<hitTree->GetEntriesFast(); ++event) {
  hitTree->GetEntry(event);
  for (auto hits : sectorHitsArray) { // loop over sectors
    for(auto& inputgroup : *hits) {
      const int MCTrackID = inputgroup.GetTrackID();
      for(size_t hitindex = 0; hitindex < inputgroup.getSize(); ++hitindex){
	const auto& eh = inputgroup.getHit(hitindex);

	grHits2D->SetPoint(hitCounter2D++, eh.GetX(),eh.GetY(),eh.GetZ());

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
}
  
  // |=============================================================
  auto CHits = new TCanvas("CHits", "Hits on A & C side", 1200, 600);
  CHits->Divide(2,1);
  CHits->cd(1);
  grHitsA->Draw("ap");
  grHitsA->GetXaxis()->SetLimits(-250, 250);
  grHitsA->SetMinimum(-250);
  grHitsA->SetMaximum(250);
  drawSectorBoundaries();
  
  CHits->cd(2);
  grHitsC->Draw("ap");
  grHitsC->GetXaxis()->SetLimits(-250, 250);
  grHitsC->SetMinimum(-250);
  grHitsC->SetMaximum(250);
  drawSectorBoundaries();
  
  auto CDigitsXZ = new TCanvas("CDigitsXZ", "Compare Digits - Hits on A & C side", 600, 600);
  grHitsAzr->Draw("ap");
  grHitsAzr->GetXaxis()->SetLimits(-250, 250);
  grHitsAzr->SetMinimum(80);
  grHitsAzr->SetMaximum(250);

  auto CHits2D = new TCanvas("CHits2D", "Hits", 600, 600);
  grHits2D->Draw();
  grHits2D->GetXaxis()->SetLimits(-250, 250);
  grHits2D->GetXaxis()->SetLimits(-250, 250);
  grHits2D->GetYaxis()->SetLimits(-250, 250);
  grHits2D->GetYaxis()->SetLimits(-250, 250);
  grHits2D->SetMinimum(-250);
  grHits2D->SetMaximum(250);
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
