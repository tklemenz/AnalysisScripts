// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See http://alice-o2.web.cern.ch/license for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

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
#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/CDBInterface.h"
#endif

using namespace o2::TPC;

void drawSectorBoundaries();

void testHitsDigits(std::string simFile="o2sim.root",
		    std::string digiFile="tpcdigits.root")
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


  int hitCounterA = 0;
  int hitCounterC = 0;

  for(int event=0; event<hitTree->GetEntriesFast(); ++event) {
  hitTree->GetEntry(event);
  for (auto hits : sectorHitsArray) { // loop over sectors
    for(auto& inputgroup : *hits) {
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
}

  // ===| process the digits |==================================================
  TFile *digitFile = TFile::Open(digiFile.data());
  TTree *digitTree = (TTree *)gDirectory->Get("o2sim");
  
  std::vector<o2::TPC::Digit> *digitsArray[Sector::MAXSECTOR];

  for (int s=0;s<Sector::MAXSECTOR;++s){
    digitsArray[s] = nullptr;
    std::stringstream sectornamestr;
    sectornamestr << "TPCDigit_" << s;
    digitTree->SetBranchAddress(sectornamestr.str().c_str(), &digitsArray[s]);
  }
  
  const Mapper& mapper = Mapper::instance();

  TGraph *grDigitsA = new TGraph();
  grDigitsA->SetMarkerColor(kGreen+2);

  TGraph *grDigitsC = new TGraph();
  grDigitsC->SetMarkerColor(kGreen+2);

  TGraph *grDigitsAzr = new TGraph();
  grDigitsAzr->SetMarkerColor(kGreen+2);

  TGraph *grDigitsCzr = new TGraph();
  grDigitsCzr->SetMarkerColor(kGreen+2);

  int digiCounterA = 0;
  int digiCounterC = 0;
  for(int event=0; event<digitTree->GetEntriesFast(); ++event) {
  digitTree->GetEntry(event);
  for(auto digits : digitsArray) {
  if(!digits) continue;
    for(auto& digit : *digits) { 
      const CRU cru(digit.getCRU());
      const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
      const int row               = digit.getRow();
      const GlobalPadNumber pad   = mapper.globalPadNumber(PadPos(row, digit.getPad()));
      const PadCentre& padCentre  = mapper.padCentre(pad);
      const float localYfactor    = (cru.side()==Side::A)?-1.f:1.f;
      float zPosition       = sampa.getZfromTimeBin(digit.getTimeStamp(), cru.side());

      LocalPosition3D posLoc(padCentre.X(), localYfactor*padCentre.Y(), zPosition);
      GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

      const float digiX = posGlob.X();
      const float digiY = posGlob.Y();
      const float digiZ = zPosition;

      if(cru.side() == Side::A) {
        grDigitsA->SetPoint(digiCounterA, digiX, digiY);
        grDigitsAzr->SetPoint(digiCounterA++, digiZ, TMath::Sqrt(digiX*digiX + digiY*digiY));
      }
      if(cru.side() == Side::C) {
        grDigitsC->SetPoint(digiCounterC, digiX, digiY);
        grDigitsCzr->SetPoint(digiCounterC++, digiZ, TMath::Sqrt(digiX*digiX + digiY*digiY));
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

  auto CDigits = new TCanvas("CDigits", "Digits on A & C side", 1200, 600);
  CDigits->Divide(2,1);
  CDigits->cd(1);
  grDigitsA->Draw("ap");
  grDigitsA->GetXaxis()->SetLimits(-250, 250);
  grDigitsA->SetMinimum(-250);
  grDigitsA->SetMaximum(250);
  drawSectorBoundaries();
  
  CDigits->cd(2);
  grDigitsC->Draw("ap");
  grDigitsC->GetXaxis()->SetLimits(-250, 250);
  grDigitsC->SetMinimum(-250);
  grDigitsC->SetMaximum(250);
  drawSectorBoundaries();

  auto CCompare = new TCanvas("CCompare", "Compare Digits - Hits on A & C side", 1200, 600);
  CCompare->Divide(2,1);
  CCompare->cd(1);
  grHitsA->Draw("ap");
  grHitsA->GetXaxis()->SetLimits(-250, 250);
  grHitsA->SetMinimum(-250);
  grHitsA->SetMaximum(250);
  grDigitsA->Draw("p");
  drawSectorBoundaries();
  
  CCompare->cd(2);
  grHitsC->Draw("ap");
  grHitsC->GetXaxis()->SetLimits(-250, 250);
  grHitsC->SetMinimum(-250);
  grHitsC->SetMaximum(250);
  grDigitsC->Draw("p");
  drawSectorBoundaries();
  
  auto CDigitsXZ = new TCanvas("CDigitsXZ", "Compare Digits - Hits on A & C side", 600, 600);
  grHitsAzr->Draw("ap");
  grHitsAzr->GetXaxis()->SetLimits(-250, 250);
  grHitsAzr->SetMinimum(80);
  grHitsAzr->SetMaximum(250);
  grDigitsAzr->Draw("p");
  grHitsCzr->Draw("p");
  grHitsCzr->GetXaxis()->SetLimits(-250, 250);
  grHitsCzr->SetMinimum(-250);
  grHitsCzr->SetMaximum(250);
  grDigitsCzr->Draw("p"); 
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
