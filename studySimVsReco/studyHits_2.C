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

void drawSectorBoundaries();

void studyHits_2(TString FileList, int sector=0, const char* Filename="studyHits.root",
               const char *OutputPath="/home/gu74yub/SimAna/Results/compareSimvsReco/hits", int iEv=-1) {

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(0.5);
  gStyle->SetTitleSize(24);

  TH1F *hX                = new TH1F("hX","; X [cm]; Counts",1000,50,150);
  TH1F *hY                = new TH1F("hY","; Y [cm]; Counts",1000,0,80);
  TH1F *hZ                = new TH1F("hZ","; Z [cm]; Counts",1000,240,250);
  TH1F *hEnergyLoss       = new TH1F("hEnergyLoss","; E loss; Counts",200,0,200);
  TH1F *hEnergyLossOneCM  = new TH1F("hEnergyLossOneCM","; E loss; Counts",200,0,200);
  TH1F *hNHitsPerEvent    = new TH1F("hNHitsPerEvent","; Number of Hits in Event; Counts",400,0,400);
  TH1F *hNHitsOneCM       = new TH1F("hNHitsOneCM","; Number of Hits within 1 cm; Counts",15,0,15);
  TH1F *hIntegratedELossInOneCM = new TH1F("hNIntElossInOneCM","; E loss within 1 cm; Counts",1000,0,1000);

  std::vector<o2::TPC::HitGroup> *sectorHitsArray[Sector::MAXSECTOR];

/// ======| get files ready |=============================================================

  std::vector<o2::TPC::HitGroup> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress("TPCHitsShiftedSector0",&fChainEvent);

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

    float startX = .0;
    float startY = .0;
    float startZ = .0;
    float ELossSum = .0;
    int inputGroupCounter = 0;
    int oneCmCounter = 0;

    for (int s=0;s<Sector::MAXSECTOR;++s) {

      if(sector != s) continue;

      sectorHitsArray[s] = nullptr;
      std::stringstream sectornamestr;
      sectornamestr << "TPCHitsShiftedSector" << s;
      hitTree->SetBranchAddress(sectornamestr.str().c_str(), &sectorHitsArray[s]);

    }

    for (int event=0; event<hitTree->GetEntriesFast(); ++event) {
      if (iEv != -1 && event != iEv) continue;

      hitTree->GetEntry(event);

      for (auto& inputgroup : *sectorHitsArray[sector]) {
        const int MCTrackID = inputgroup.GetTrackID();

        int nHitsinGroup = 0;

        nHitsinGroup = inputgroup.getSize();
        hNHitsPerEvent->Fill(nHitsinGroup);
        std::cout << "hitgroup size: " << nHitsinGroup << std::endl;

        oneCmCounter = 0;
        startX = .0;
        startY = .0;
        startZ = .0;
        ELossSum = .0;

        for (size_t hitindex = 0; hitindex<inputgroup.getSize(); ++hitindex) {
          const auto& eh = inputgroup.getHit(hitindex);

          hEnergyLoss->Fill(eh.GetEnergyLoss());

          float posX = eh.GetX();
          float posY = eh.GetY();
          float posZ = eh.GetZ();

          hX->Fill(posX);
          hY->Fill(posY);
          hZ->Fill(posZ);

          if (hitindex == 50) {
            startX = posX;
            startY = posY;
            startZ = posZ;
          }

          const float distance = TMath::Sqrt((startX-posX) * (startX-posX) +
                                             (startY-posY) * (startY-posY) +
                                             (startZ-posZ) * (startZ-posZ));

          if (startX != 0 && startY != 0 && startZ != 0) {
            if (std::abs(distance)< 1.f) {
              hEnergyLossOneCM->Fill(eh.GetEnergyLoss());
              ELossSum += eh.GetEnergyLoss();
              ++oneCmCounter;
            }
          }
        }
        hNHitsOneCM->Fill(oneCmCounter);
        hIntegratedELossInOneCM->Fill(ELossSum);
        ++inputGroupCounter;
      }
    }
    //std::cout << "\n\n\n inputGroups: " << inputGroupCounter << std::endl;
    //std::cout << "hits in group: " << nHitsinGroup << std::endl;
    std::cout << "hits within 1 cm: " << oneCmCounter << std::endl;


  }

///==================================================================================

  TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");

  OutFile->WriteObject(hX, "hX");
  OutFile->WriteObject(hY, "hY");
  OutFile->WriteObject(hZ, "hZ");
  OutFile->WriteObject(hEnergyLoss, "hEnergyLoss");
  OutFile->WriteObject(hEnergyLossOneCM, "hEnergyLossOneCM");
  OutFile->WriteObject(hIntegratedELossInOneCM, "hIntegratedELossOneCM");
  OutFile->WriteObject(hNHitsPerEvent, "hNHitsPerEvent");
  OutFile->WriteObject(hNHitsOneCM, "hNHitsOneCM");
}
