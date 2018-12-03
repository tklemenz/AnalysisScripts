#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TSystem.h"
#include "TPCBase/Digit.h"
#include "TPCBase/CRU.h"
#include "TPCBase/PadRegionInfo.h"
#include "TPCBase/Mapper.h"
#include "TPCSimulation/SAMPAProcessing.h"

using namespace o2::TPC;
//using namespace o2::dataformats;

void compareSimvsReco_digits(TString FileList="")
{

  const Mapper& mapper = Mapper::instance();

  TH1F *hDigiCharge  = new TH1F("hCharge","; Q_{digit}; Counts",1000,0,1000);
  TH2F *hChargeVsX   = new TH2F("hChargeVsX","; X [cm]; Q_{digit}",1000,60,160,1000,0,1000);
  TH2F *hChargeVsY   = new TH2F("hChargeVsY","; Y [cm]; Q_{digit}",500,0,50,1000,0,1000);
  TH2F *hChargeVsZ   = new TH2F("hChargeVsZ","; Z [cm]; Q_{digit}",150,150,300,1000,0,1000);
  TH2F *hChargeVsTime = new TH2F("hChargeVsTime","; time; Q_{digit}", 1000,0,1000,1000,0,1000);

  //TFile *digitFile = TFile::Open(digiFile.data());

  //if (!digitFile->IsOpen()) {
    //std::cout << std::endl << "no file opened" << std::endl;
    //return;
  //}

  std::vector<Digit> *fChainEvent=0;

  TChain fChain("o2sim");

  fChain.SetBranchAddress("TPCDigit",&fChainEvent);

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



  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;

    TTree *digitTree = (TTree*)TreeFile->Get("o2sim");

    std::vector<Digit> *digitsArray = nullptr;
    digitTree->SetBranchAddress("TPCDigit", &digitsArray);

    for (int iEv = 0; iEv<digitTree->GetEntriesFast(); ++iEv) {
      //std::cout << std::endl << "number of events: " << digitTree->GetEntriesFast() << std::endl;
      digitTree->GetEntry(iEv);

      for (auto digitObject : *digitsArray) {
        //Digit oder DigitMC???

        //std::cout << std::endl << "loop over digits" << std::endl;

        Digit *inputdigit = static_cast<Digit*>(&digitObject);

        const CRU cru(inputdigit->getCRU());
        const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
        const int rowInSector       = inputdigit->getRow() + region.getGlobalRowOffset();
        const GlobalPadNumber pad   = mapper.globalPadNumber(PadPos(rowInSector, inputdigit->getPad()));
        const PadCentre& padCentre  = mapper.padCentre(pad);
        const float localYfactor    = (cru.side()==Side::A)?-1.f:1.f;
        float zPosition             = SAMPAProcessing::getZfromTimeBin(inputdigit->getTimeStamp(), cru.side());

        LocalPosition3D posLoc(padCentre.X(), localYfactor*padCentre.Y(), zPosition);
        GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

        const DigitPos digiPadPos = mapper.findDigitPosFromGlobalPosition(posGlob);

        /// take only sector 0
        if (digiPadPos.getCRU() > 3) continue;
        const PadSecPos digiSecPos = digiPadPos.getPadSecPos();
        const Sector digiSector = digiSecPos.getSector();
        if (digiSector.getSector() != 0) continue;

        const float digitX         = posGlob.X();
        const float digitY         = posGlob.Y();
        const float digitZ         = zPosition;
        const float digitTimeStamp = inputdigit->getTimeStamp();

        const float Q = digitObject.getChargeFloat();

        hDigiCharge->Fill(Q);
        hChargeVsX->Fill(digitX,Q);
        hChargeVsY->Fill(digitY,Q);
        hChargeVsZ->Fill(digitZ,Q);
        hChargeVsTime->Fill(digitTimeStamp,Q);
      }
    }
  }

  TFile *Outfile = new TFile("/home/gu74yub/SimAna/Results/compareSimvsReco/digits/compareSimvsReco_digits.root","recreate");

  Outfile->WriteObject(hDigiCharge, "hDigiCharge");
  Outfile->WriteObject(hChargeVsX, "hChargeVsX");
  Outfile->WriteObject(hChargeVsY, "hChargeVsY");
  Outfile->WriteObject(hChargeVsZ, "hChargeVsZ");
  Outfile->WriteObject(hChargeVsTime, "hChargeVsTimeStamp");

  delete Outfile;
  delete hDigiCharge;
  delete hChargeVsX;
  delete hChargeVsY;
  delete hChargeVsZ;
  delete hChargeVsTime;
}
