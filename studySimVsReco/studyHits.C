#if !defined(__CLING__) || defined(__ROOTCLING__)

#include <vector>
#include <fstream>
#include <iostream>

#include "TSystem.h"
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
#include "TChain.h"
#include "TPCSimulation/Point.h"
#include "TPCBase/Digit.h"
#include "TPCSimulation/Digitizer.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "ReconstructionDataFormats/Track.h"
#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#endif

using namespace o2::TPC;

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

void studyHits(TString FileList, int sector=0, const char* Filename="studyHits.root",
               const char *OutputPath="/scratch2/tklemenz/SimAna/Results/newO2/compareSimvsReco/hits", int iEv=-1) {

  int shiftedSector = sector + 1;

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetTitleSize(24);

  TH1F *hEnergyLoss       = new TH1F("hEnergyLoss","; E loss; Counts",200,0,200);
  TH1F *hNHitsPerEvent    = new TH1F("hNHitsPerEvent","; Number of Hits in Event; Counts",400,0,400);
  TH1F *hNElectronsPerCM  = new TH1F("NElectronsPerCM","; Electrons/cm; Counts",300,0,600);

  TGraph *gHitsXY = new TGraph();
  gHitsXY->SetTitle("X-Y plane ; x [cm]; y [cm]");

  TGraph *gHitsXZ = new TGraph();
  gHitsXZ->SetTitle("X-Z plane ; x [cm]; z [cm]");

  TGraph *gHitsYZ = new TGraph();
  gHitsYZ->SetTitle("Y-Z plane ; y [cm]; z [cm]");

  std::vector<o2::TPC::HitGroup> *sectorHitsArray[Sector::MAXSECTOR];

/// ======| get files ready |=============================================================

  std::vector<o2::TPC::HitGroup> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress(Form("TPCHitsShiftedSector%i",shiftedSector),&fChainEvent);

  TString allFiles;
  if (FileList.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else if (FileList.EndsWith(".root")){
    allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else {
    allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      file = FileList+"/"+file;
      fChain.Add(file);
    }
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }

///===================================================================================
///======| process hits |=============================================================

  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile) {
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;
    TTree *hitTree = (TTree*)TreeFile->Get("o2sim");

    float nElectrons = .0;
    float oldDistance = .0;
    float refX = .0;
    float refY = .0;
    float refZ = .0;
    int inputGroupCounter = 0;
    int hitCounter = 0;

    for (int s=0;s<Sector::MAXSECTOR;++s) {

      if(shiftedSector != s) continue;

      sectorHitsArray[s] = nullptr;
      std::stringstream sectornamestr;
      sectornamestr << "TPCHitsShiftedSector" << s;
      hitTree->SetBranchAddress(sectornamestr.str().c_str(), &sectorHitsArray[s]);

    }

    for (int event=0; event<hitTree->GetEntriesFast(); ++event) {
      if (iEv != -1 && event != iEv) continue;

      hitTree->GetEntry(event);

      nElectrons = .0;
      oldDistance = .0;
      refX = .0;
      refY = .0;
      refZ = .0;

      for (auto& inputgroup : *sectorHitsArray[shiftedSector]) {
        const int MCTrackID = inputgroup.GetTrackID();

        int nHitsinGroup = 0;

        nHitsinGroup = inputgroup.getSize();
        hNHitsPerEvent->Fill(nHitsinGroup);

        /// old place of 0-setters

        for (size_t hitindex = 0; hitindex<inputgroup.getSize(); ++hitindex) {
          const auto& eh = inputgroup.getHit(hitindex);

          const float posX  = eh.GetX();
          const float posY  = eh.GetY();
          const float posZ  = eh.GetZ();
          const float ELoss = eh.GetEnergyLoss();

          hEnergyLoss->Fill(ELoss);

          gHitsXY->SetPoint(hitCounter, posX, posY);
          gHitsXZ->SetPoint(hitCounter, posX, posZ);
          gHitsYZ->SetPoint(hitCounter, posY, posZ);

          const float distance = TMath::Sqrt((posX-refX) * (posX-refX) +
                                             (posY-refY) * (posY-refY) +
                                             (posZ-refZ) * (posZ-refZ));

          if (std::abs(distance) < 5.f) {                     // hit within current 5cm bin

            nElectrons += ELoss;

          }
          else if (std::abs(distance) > 10.f) {               // next Track

            refX = posX;
            refY = posY;
            refZ = posZ;
            nElectrons = 0;

          }
          else {                                              // start new 5cm bin

            hNElectronsPerCM->Fill(nElectrons/5);			  // take 5 cm and divide by 5 in the end
            refX = posX;
            refY = posY;
            refZ = posZ;
            nElectrons = 0;

          }
          ++hitCounter;
        }
        ++inputGroupCounter;
      }
    }

    //std::cout << "\n\n\n inputGroups: " << inputGroupCounter << std::endl;
    //std::cout << "hits in group: " << nHitsinGroup << std::endl;
    //std::cout << "hits within 1 cm: " << oneCmCounter << std::endl;

  }

  hEnergyLoss->Scale(1/double(hEnergyLoss->GetEntries()));
  hNElectronsPerCM->Scale(1/double(hNElectronsPerCM->GetEntries()));
  hNHitsPerEvent->Scale(1/double(hNHitsPerEvent->GetEntries()));

///==================================================================================

  TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");

  OutFile->WriteObject(hEnergyLoss, "hEnergyLoss");
  OutFile->WriteObject(hNElectronsPerCM, "hNElectronsPerCM");
  OutFile->WriteObject(hNHitsPerEvent, "hNHitsPerEvent");

  TCanvas *cXY = new TCanvas("cXY","Hits in X-Y plane");
  gHitsXY->Draw("ap");
  gHitsXY->GetXaxis()->SetLimits(-250, 250);
  gHitsXY->SetMinimum(-250);
  gHitsXY->SetMaximum(250);
  drawSectorBoundaries();
  cXY->Print(Form("%s/xy_plane.png", OutputPath));

  TCanvas *cXZ = new TCanvas("cXZ","Hits in X-Z plane");
  gHitsXZ->Draw("ap");
  gHitsXZ->GetXaxis()->SetLimits(-250,250);
  gHitsXZ->SetMinimum(-250);
  gHitsXZ->SetMaximum(250);
  cXZ->Print(Form("%s/xz_plane.png", OutputPath));

  TCanvas *cYZ = new TCanvas("cYZ","Hits in Y-Z plane");
  gHitsYZ->Draw("ap");
  gHitsYZ->GetXaxis()->SetLimits(-250,250);
  gHitsYZ->SetMinimum(-250);
  gHitsYZ->SetMaximum(250);
  cYZ->Print(Form("%s/yz_plane.png", OutputPath));

}
