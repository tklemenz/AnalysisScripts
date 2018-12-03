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
#include "TPCSimulation/SAMPAProcessing.h"
#include "TPCBase/CDBInterface.h"

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

void studyHitsAndDigits(TString HitFile, TString DigitFile, int sector=0, const char* Filename="studyHitsAndDigits.root",
               const char *OutputPath="/home/gu74yub/SimAna/Results/newO2/showHits", int iEv=-1) {

  int ShiftedSector = sector + 1;

  auto& cdb = CDBInterface::instance();
  cdb.setUseDefaults();

  Mapper &mapper = Mapper::instance();
  SAMPAProcessing &sampa = SAMPAProcessing::instance();
  int digitCounter = 0;

  gStyle->SetMarkerStyle(20);
  gStyle->SetMarkerSize(1);
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

  TGraph *gDigitsXY = new TGraph();
  gDigitsXY->SetTitle("X-Y plane ; x [cm]; y [cm]");

  TGraph *gDigitsXZ = new TGraph();
  gDigitsXZ->SetTitle("X-Z plane ; x [cm]; z [cm]");

  TGraph *gDigitsYZ = new TGraph();
  gDigitsYZ->SetTitle("Y-Z plane ; y [cm]; z [cm]");

  TGraph *gHitsRowZ = new TGraph();
  gHitsRowZ->SetTitle("z vs row; pad row; z [cm]");

  TGraph *gDigitsRowZ = new TGraph();
  gDigitsRowZ->SetTitle("z vs row; pad row; z [cm]");

  TGraph *gDigitsRowTimeStamp = new TGraph();
  gDigitsRowTimeStamp->SetTitle("timeStamp vs row; pad row; timestamp");

  TGraph *gHitsRowTimeStamp = new TGraph();
  gHitsRowTimeStamp->SetTitle("timeStamp vs row; pad row; timeStamp");


  std::vector<o2::TPC::HitGroup> *sectorHitsArray[Sector::MAXSECTOR];

/// ======| get files ready |=============================================================

  std::vector<o2::TPC::HitGroup> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress("TPCHitsShiftedSector0",&fChainEvent);

  TString allFiles;
  if (HitFile.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",HitFile.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else if (HitFile.EndsWith(".root")){
    allFiles=gSystem->GetFromPipe(Form("ls %s",HitFile.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  else {
    allFiles=gSystem->GetFromPipe(Form("ls %s",HitFile.Data()));

    TObjArray *arr = allFiles.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      file = HitFile+"/"+file;
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

      if(ShiftedSector != s) continue;

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

      for (auto& inputgroup : *sectorHitsArray[ShiftedSector]) {
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
          const float hitTime = eh.GetTime();

          const GlobalPosition3D hitPosGlob(posX,posY,posZ);
          const DigitPos hitDigitPos =  mapper.findDigitPosFromGlobalPosition(hitPosGlob);
          const PadSecPos hitPadSecPos = hitDigitPos.getPadSecPos();
          const PadPos hitPadPos = hitPadSecPos.getPadPos();
          const float hitRow = hitPadPos.getRow();

          gHitsRowZ->SetPoint(hitCounter, hitRow, posZ);

          hEnergyLoss->Fill(ELoss);

          gHitsXY->SetPoint(hitCounter, posX, posY);
          gHitsXZ->SetPoint(hitCounter, posX, posZ);
          gHitsYZ->SetPoint(hitCounter, posY, posZ);
          gHitsRowTimeStamp->SetPoint(hitCounter, hitRow, hitTime);

          const float distance = TMath::Sqrt((posX-refX) * (posX-refX) +
                                             (posY-refY) * (posY-refY) +
                                             (posZ-refZ) * (posZ-refZ));

          if (std::abs(distance) < 1.f) {                     // hit within current cm bin

            nElectrons += ELoss;

          }
          else if (std::abs(distance) > 10.f) {               // next Track

            refX = posX;
            refY = posY;
            refZ = posZ;
            nElectrons = 0;

          }
          else {                                              // start new cm bin

            hNElectronsPerCM->Fill(nElectrons);
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


///===================================================================================
///======| process digits |===========================================================

  std::vector<Digit> *sectorDigitsArray[Sector::MAXSECTOR];
  int sectorDigits = sector;

  /// ======| get files ready |=============================================================

  std::vector<Digit> *fChainEventDigits=0;

  TChain fChainDigits("o2sim");

  fChainDigits.SetBranchAddress("TPCDigit_0",&fChainEventDigits);

  TString allFilesDigits;
  if (DigitFile.EndsWith(".txt")){
    allFilesDigits=gSystem->GetFromPipe(Form("cat %s",DigitFile.Data()));

    TObjArray *arr = allFilesDigits.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChainDigits.Add(file);
    }
  }
  else if (DigitFile.EndsWith(".root")){
    allFilesDigits=gSystem->GetFromPipe(Form("ls %s",DigitFile.Data()));

    TObjArray *arr = allFilesDigits.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChainDigits.Add(file);
    }
  }
  else {
    allFilesDigits=gSystem->GetFromPipe(Form("ls %s",DigitFile.Data()));

    TObjArray *arr = allFilesDigits.Tokenize("\n");

    for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
      TString file=arr->At(ifile)->GetName();
      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      file = DigitFile+"/"+file;
      fChainDigits.Add(file);
    }
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *ChainEntriesDigits = fChainDigits.GetListOfFiles();
  for (int ifile=0; ifile<ChainEntriesDigits->GetEntriesFast(); ++ifile){
    cout<<ChainEntriesDigits->At(ifile)->GetTitle()<<endl;
  }



  for (int ifile=0; ifile<ChainEntriesDigits->GetEntriesFast(); ++ifile) {
    TFile *TreeFile = new TFile(Form("%s", ChainEntriesDigits->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<ChainEntriesDigits->At(ifile)->GetTitle()<<endl;
    TTree *digitTree = (TTree*)TreeFile->Get("o2sim");

    for (int s=0;s<Sector::MAXSECTOR;++s) {

      if(sectorDigits != s) continue;

      sectorDigitsArray[s] = nullptr;
      std::stringstream sectornamestrDigits;
      sectornamestrDigits << "TPCDigit_" << s;
      digitTree->SetBranchAddress(sectornamestrDigits.str().c_str(), &sectorDigitsArray[s]);
    }


    for (int event=0; event<digitTree->GetEntriesFast(); ++event) {
      if (iEv != -1 && event != iEv) continue;
      digitTree->GetEntry(event);
      for (auto& digit : *sectorDigitsArray[sectorDigits]) {
        const CRU cru(digit.getCRU());

        //const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
        const int row = digit.getRow();   // in CRU coordinates
        const GlobalPadNumber pad = mapper.globalPadNumber(PadPos(row, digit.getPad()));
        const PadCentre& padCentre = mapper.padCentre(pad);
        const float localYfactor = (cru.side() == Side::A) ? -1.f : 1.f;
        const float timeStamp = digit.getTimeStamp();
        float zPosition = sampa.getZfromTimeBin(digit.getTimeStamp(), cru.side()); //stimmt evtl nicht, timeBin /= timeStamp

        LocalPosition3D posLoc(padCentre.X(), localYfactor * padCentre.Y(), zPosition);
        GlobalPosition3D posGlob = mapper.LocalToGlobal(posLoc, cru.sector());

        //const float relPadPos = padCentre - pad;
        //std::cout << std::endl << "pad: " << pad << "\t padCentre: " << padCentre << " \t relPadPos: " << relPadPos << std::endl;

        float digitX = posGlob.X();
        float digitY = posGlob.Y();
        float digitZ = posGlob.Z();

        gDigitsRowZ->SetPoint(digitCounter, row, digitZ);
        gDigitsRowTimeStamp->SetPoint(digitCounter, row, timeStamp);
        gDigitsXY->SetPoint(digitCounter, digitX, digitY);
        gDigitsXZ->SetPoint(digitCounter, digitX, digitZ);
        gDigitsYZ->SetPoint(digitCounter, digitY, digitZ);

        ++digitCounter;
      }
    }

  }
///==================================================================================

  TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");

  OutFile->WriteObject(hEnergyLoss, "hEnergyLoss");
  OutFile->WriteObject(hNElectronsPerCM, "hNElectronsPerCM");
  OutFile->WriteObject(hNHitsPerEvent, "hNHitsPerEvent");

  TCanvas *cXY = new TCanvas("cXY","Hits in X-Y plane");
  gDigitsXY->SetMarkerSize(0.5);
  gDigitsXY->SetMarkerColor(kRed);
  gHitsXY->Draw("ap");
  gDigitsXY->Draw("p,same");
  gHitsXY->GetXaxis()->SetLimits(-250, 250);
  gHitsXY->SetMinimum(-250);
  gHitsXY->SetMaximum(250);
  drawSectorBoundaries();
  cXY->Print(Form("%s/xy_plane.png", OutputPath));

  TCanvas *cXZ = new TCanvas("cXZ","Hits in X-Z plane");
  gDigitsXZ->SetMarkerSize(0.5);
  gDigitsXZ->SetMarkerColor(kRed);
  gHitsXZ->Draw("ap");
  gDigitsXZ->Draw("p,same");
  gHitsXZ->GetXaxis()->SetLimits(-250,250);
  gHitsXZ->SetMinimum(-250);
  gHitsXZ->SetMaximum(250);
  cXZ->Print(Form("%s/xz_plane.png", OutputPath));

  TCanvas *cYZ = new TCanvas("cYZ","Hits in Y-Z plane");
  gDigitsYZ->SetMarkerSize(0.5);
  gDigitsYZ->SetMarkerColor(kRed);
  gHitsYZ->Draw("ap");
  gDigitsYZ->Draw("p,same");
  gHitsYZ->GetXaxis()->SetLimits(-250,250);
  gHitsYZ->SetMinimum(-250);
  gHitsYZ->SetMaximum(250);
  cYZ->Print(Form("%s/yz_plane.png", OutputPath));

  TCanvas *cRowZ = new TCanvas("cRowZ","Z position vs Pad Row");
  gDigitsRowZ->SetMarkerSize(0.5);
  gDigitsRowZ->SetMarkerColor(kRed);
  gHitsRowZ->Draw("ap");
  gDigitsRowZ->Draw("p,same");
  cRowZ->Print(Form("%s/row_z.png",OutputPath));

  TCanvas *cRowTime = new TCanvas("cRowTimeStamp","TimeStamp position vs Pad Row");
  gDigitsRowTimeStamp->SetMarkerSize(0.5);
  gDigitsRowTimeStamp->SetMarkerColor(kRed);
  //gHitsRowTimeStamp->Draw("ap");
  gDigitsRowTimeStamp->Draw("ap");
  cRowTime->Print(Form("%s/row_timeStamp.png",OutputPath));
}
