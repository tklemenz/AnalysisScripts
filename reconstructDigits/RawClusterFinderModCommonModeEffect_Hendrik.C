// Copyright CERN and copyright holders of ALICE O2. This software is
// distributed under the terms of the GNU General Public License v3 (GPL
// Version 3), copied verbatim in the file "COPYING".
//
// See https://alice-o2.web.cern.ch/ for full licensing information.
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef ALICEO2_TPC_CONVERTRAWCLUSTERS_C_
#define ALICEO2_TPC_CONVERTRAWCLUSTERS_C_

/// \file   RawClusterFinder.h
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de
//Modificated version to investigate the influence of the common mode effect

#include <vector>
#include <memory>
#include <fstream>

#include "Rtypes.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TSystem.h"
#include "TGraph.h"
#include "TLinearFitter.h"

#include "TPCBase/Defs.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/CRU.h"
#include "TPCCalibration/CalibRawBase.h"
#include "TPCSimulation/HwClusterer.h"
#include "TPCSimulation/BoxClusterer.h"
#include "TPCSimulation/ClusterContainer.h"
#include "TPCReconstruction/DigitData.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/TrackTPC.h"
#include "TPCSimulation/Cluster.h"

#include "TCanvas.h"
#include "TStyle.h"


namespace o2
{
namespace TPC
{


/// \brief Raw cluster conversion
///
/// This class is used to produce pad wise pedestal and noise calibration data
///
/// origin: TPC
/// \author Jens Wiechula, Jens.Wiechula@ikf.uni-frankfurt.de

void SaveCanvas(TCanvas* c_in);
double Gamma4(double* x, double* par);

class RawClusterFinder : public CalibRawBase
{
  public:
    void CommonModeEffectSearch(TClonesArray* useArrCluster, std::vector<std::unique_ptr<Digit>>& useDigitVector, TH1D* useHistQTrackPads, TH1D* useHistQTrackExcluded, TH1F* useHistTZeroDist, TH1F* useHistTZeroFineDist, int useRunNr, int useEventNr, CalDet<float> *pedestal, int& useNActive, int& useSignalPads, bool& useAskTZeroFine, int& useEventCounter, const int& useTimeAcceptLow, const int& useTimeAcceptUp, TH2D* useHistCmmMdEffPrfl, TProfile2D* useHistCmmMdEffPrflP, TH1F* useHistTimeDistAll, TH1F* useHistTZeroFineSlopes, int useFine, float useCherenkovValue, TH1D* useEMinusSgnl, TH1D* useEMinusCme, TH1D* usePionSgnl, TH1D* usePionCme, int& useNActiveEMinus, int& useNActivePions, int& useSignalPadsEMinus, int& useSignalPadsPions, int& useEventCounterEMinus, int& useEventCounterPions, std::vector<float>& useVecClusterMax, bool& useAcceptedEvent, std::vector<float>& useVecClusterMaxFull);     // ModCmmnMdEff

    enum class ClustererType : char {
      Box,  ///< use box clusterer
      HW    ///< use HW clusterer
    };


    struct EventHeader
    {
      int run;
      float cherenkovValue;
      int beamMomentum;
      int powerSupply;
      int HVSettings;
      int trigger;
      int dataType;
      int driftFieldStrength;
      int runType;

    };

    using vectorType = std::vector<float>;

    /// default constructor
    RawClusterFinder(PadSubset padSubset = PadSubset::ROC) : CalibRawBase(padSubset), mClustererType(ClustererType::HW), mPedestals(nullptr), mVectorDigits() {;}

    /// default destructor
    virtual ~RawClusterFinder() = default;

    /// set clusterer type
    void setClustererType(ClustererType clustererType) { mClustererType = clustererType; }

    /// not used
    Int_t updateROC(const Int_t sector, const Int_t row, const Int_t pad,
                    const Int_t timeBin, const Float_t signal) final { return 0;}

    /// 
    Int_t updateCRU(const CRU& cru, const Int_t row, const Int_t pad,
                    const Int_t timeBin, const Float_t signal) final;


    void TrackFinder(int nTries, int allowedMisses, int minClusters, bool reverse, TClonesArray *pclusterArray, vector<TrackTPC> &arrTracks);
    static void GetFit(TGraph *trackGraph, float& slope, float& offset, float& slopeError, float& offsetError, float& chi2);

    void setPedestals(CalPad* pedestals) { mPedestals = pedestals; }
    static void processEvents(TString fileInfo, TString pedestalFile, TString outputFileName="clustersModCommonModeEffect.root", Int_t maxEvents=-1, TString cherenkovFile="", ClustererType clustererType=ClustererType::HW);

    std::vector<std::unique_ptr<Digit>>& getDigitVector() { return mVectorDigits; }

    /// Dummy end event
    virtual void endEvent() final {};

    TTree *mTOut;

  private:
    ClustererType     mClustererType;
    CalPad       *mPedestals;
    std::vector<std::unique_ptr<Digit>> mVectorDigits;

    /// dummy reset
    void resetEvent() final { mVectorDigits.clear(); }
};

Int_t RawClusterFinder::updateCRU(const CRU& cru, const Int_t row, const Int_t pad,
                                     const Int_t timeBin, const Float_t signal)
{
  float corrSignal = signal;

  // ===| get pedestal |========================================================
  if (mPedestals) {
    corrSignal -= mPedestals->getValue(cru, row, pad);
  }

  // ===| add new digit |=======================================================
  mVectorDigits.emplace_back(new DigitData(cru, corrSignal, row, pad, timeBin));

  return 1;
}

void RawClusterFinder::processEvents(TString fileInfo, TString pedestalFile, TString outputFileName, Int_t maxEvents, TString cherenkovFile, ClustererType clustererType)
{

  // ===| create raw converter |================================================
  RawClusterFinder converter;
  converter.setupContainers(fileInfo);

  // ===| load pedestals |======================================================
  TFile f(pedestalFile);
  CalDet<float> *pedestal = nullptr;
  if (f.IsOpen() && !f.IsZombie()) {
    f.GetObject("Pedestals", pedestal);
    printf("pedestal: %.2f\n", pedestal->getValue(CRU(0), 0, 0));
  }

  // ===| output file and container |===========================================
  TClonesArray arrCluster("o2::TPC::Cluster");
  float cherenkovValue = 0.;
  int runNumber = 0;
  int beamMomentum = 0;
  int powerSupply = 0;
  int HVSettings = 0;
  int trigger = 0;
  int driftFieldStrength = 0;
  int dataType = 0;
  int runType = 0;


  TFile fout(outputFileName,"recreate");
  TTree t("cbmsim","cbmsim");
  t.Branch("TPCClusterHW",           &arrCluster);
  t.Branch("cherenkovValue",         &cherenkovValue);
  t.Branch("runNumber",              &runNumber);
  t.Branch("beamMomentum",           &beamMomentum);
  t.Branch("powerSupply",             &powerSupply);
  t.Branch("HVSettings",             &HVSettings);
  t.Branch("trigger",                &trigger);
  t.Branch("dataType",               &dataType);
  t.Branch("driftFieldStrength",     &driftFieldStrength);                                                                                           
  t.Branch("runType",                &runType);


  // ===| output track file and container |======================================
  //TFile foutTracks("tracksLinear_ModCommonModeEffect.root", "recreate");
  //TTree tOut("events","events");

  //vector<TrackTPC> arrTracks;
  //EventHeader eventHeader;


  //tOut.Branch("header", &eventHeader, "run/I:cherenkovValue/F");
  //tOut.Branch("Tracks", &arrTracks);


  // ===| input cherenkov file |================================================
  ifstream istr(cherenkovFile.Data());

  // ===| fill header |=========================================================
  //eventHeader.run = TString(gSystem->Getenv("RUN_NUMBER")).Atoi();
  //runNumber = eventHeader.run;
  runNumber = TString(gSystem->Getenv("RUN_NUMBER")).Atoi();
  TString logbookData;
  const char* logbook = gSystem->Getenv("LOGBOOK");
  printf("content of variable logbook is %s \n", logbook);

  if (logbook)
  {
    logbookData = gSystem->GetFromPipe(TString::Format("awk '{if($1 == %d) print $3,$4,$5,$6,$7,$8,$9}' %s", runNumber, logbook));
    TObjArray *arrData = logbookData.Tokenize(" ");
    printf("in logbook are %d entries \n", arrData->GetEntriesFast()); 
    if (arrData && arrData->GetEntriesFast()==7)
    {
      beamMomentum = ((TObjString*)arrData->At(0))->String().Atoi();
      printf("the value of beamMomentum is %d \n", beamMomentum);
      powerSupply = ((TObjString*)arrData->At(1))->String().Atoi();
      HVSettings = ((TObjString*)arrData->At(2))->String().Atoi();
      trigger = ((TObjString*)arrData->At(3))->String().Atoi();
      dataType = ((TObjString*)arrData->At(4))->String().Atoi();
      driftFieldStrength = ((TObjString*)arrData->At(5))->String().Atoi();
      runType = ((TObjString*)arrData->At(6))->String().Atoi();
    }
    else 
    {
      printf("Error: Could not get logbook information for run %d", runNumber);
    }
  
    delete arrData;
   }
   else
   {
      printf("Couldn't find the logbook! \n");
   }


  // ===| cluster finder |======================================================
  // HW cluster finder
  std::unique_ptr<Clusterer> cl;
  if (clustererType == ClustererType::HW) {
    HwClusterer *hwCl = new HwClusterer(HwClusterer::Processing::Parallel, 0, 0, 4);
    hwCl->setContinuousReadout(false);
    hwCl->setPedestalObject(pedestal);
    cl = std::unique_ptr<Clusterer>(hwCl);
  }
  else if (clustererType == ClustererType::Box) {
    BoxClusterer *boxCl = new BoxClusterer;
    boxCl->setPedestals(pedestal);
    cl = std::unique_ptr<Clusterer>(boxCl);
  }
  else {
    return;
  }
    
  //cl->setRequirePositiveCharge(false);
  cl->setMinQMax(5.);
  cl->Init();

  // Box cluster finder
  
  // ===| loop over all data |==================================================

  // ModCmmnMdEff; begin 
  int events = 0;
  // ---------------------------------------------------------------------------- for all accepted events-----------------------------------------------------------------
  int NActive = 0;            // introduce this integer to get the number the pads that contribute to the Common Mode Effect Signal
  int SignalPads = 0;
  int eventCounter = 0;
  // ---------------------------------------------------------------------------- for accepted electron events only ------------------------------------------------------
  int NActiveEMinus = 0;
  int SignalPadsEMinus = 0;
  int eventCounterEMinus = 0;
  // ---------------------------------------------------------------------------- for accepted pion events only ----------------------------------------------------------
  int NActivePions = 0;
  int SignalPadsPions = 0;
  int eventCounterPions = 0;

  std::vector<float> vecClusterMax (63, -5000.);   // make a vector to extract the Cluster Max of every used event
  bool acceptedEvent = false;                     // use this bool to avoid adding Cluster Max values to the Cluster Max histogram if an event was not accepted by the applied cuts
  std::vector<float> vecClusterMaxFull (63, -5000.);   


  // Parameter to control the tZero time cut and to adjust the binning of the concerning histograms
  bool askTZeroFine = true;          // change the initial value of this bool to change between a raw and a more fine cut of the time acceptance
  const int timeAcceptLow = -2;      // lowest accepted time bin under tZero respectively tZeroFine, set -50 to set tZero cut off
  const int timeAcceptUp =   5;        // highest accepted time bin over tZero respectively tZeroFine, set 50 to set tZero cut off
  const int fine = 5;               // define how fine the tZero binning of HistQTrackPads and HistQTrackExcluded should be 
  const int bins = (timeAcceptUp - timeAcceptLow +3) * fine; // is the number of bins for HistQTrackPads and HistQTrackExcluded 

  printf("max events %d\n", maxEvents);
  CalibRawBase::ProcessStatus status;
  //while (((status = converter.processEvent()) == CalibRawBase::ProcessStatus::Ok) && (maxEvents>0)?events<maxEvents:1) {
  // skip synch event
  converter.processEvent();

  TH1D* HistQTrackPads = new TH1D("HistQTrackPads", "Track signal in one run", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1D* HistQTrackExcluded = new TH1D("HistQTrackExcluded", "Common Mode signal in  one run ", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1D* HistEMinusSgnl = new TH1D("HistEMinusSgnl", "Track signal by electrons only in one run", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1D* HistEMinusCme = new TH1D("HistEMinusCme", "Common Mode signal by electrons only in one run", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1D* HistPionSgnl = new TH1D("HistPionSgnl", "Track signal by pions only in one run", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1D* HistPionCme = new TH1D("HistPionCme", "Common Mode signal by pions only in one run", bins, timeAcceptLow-1.5, timeAcceptUp+1.5);
  TH1F* HistTZeroDist = new TH1F("HistTZeroDist", Form("Distribution of TZeros from run %d", runNumber), 25, -5, 20);
  TH1F* HistTZeroFineDist = new TH1F("HistTZeroFineDist", Form("Distribution of TZerosFine from run %d", runNumber), 25*fine, -5, 20);
  TH1F* HistTimeDistAll = new TH1F("HistTimeDistAll", Form("Distribution of time of All Digits from run %d", runNumber), 100, 0, 100);
  TH2D* HistCmmMdEffPrfl = new TH2D("HistCmmMdEffPrfl", Form("Profile of the Common Mode Effect Entries, divided by the number of contributing events, run %d", runNumber), 63, 0, 63, 25, 5, 30);
  TProfile2D* HistCmmMdEffPrflP = new TProfile2D("HistCmmMdEffPrflP", Form("Profile of the Common Mode Effect Entries, divided by the number of contributing events, run %d", runNumber), 63, 0, 63, 25, 5, 30);
  TH1F* HistTZeroFineSlopes = new TH1F("HistTZeroFineSlopes", Form("Distribution of Slopes from TZeroFine Fits, run %d", runNumber), 200, -0.1, 0.1);
  TH1F* HistClusterMax = new TH1F("HistClusterMax", Form("Distributon of highest digit charge per row of all used events from run %d", runNumber), 2000, -1000, 1000); 
  TH1F* HistClusterMaxFull = new TH1F("HistClusterMaxFull", Form("Distributon of highest digit charge per row of all events from run %d", runNumber), 2000, -1000, 1000); 

  Mapper &mapper = Mapper::instance();
  // ModCmmnMdEff; end 

  while (converter.processEvent() == CalibRawBase::ProcessStatus::Ok) {
    if (maxEvents>0 && events>=maxEvents) break;             // ModCmmnMdEff; limit the number of used events for test purposes here

    printf("========| Event %4zu %d %d %hhd |========\n", converter.getPresentEventNumber(), events, maxEvents, status);

    auto &arr = converter.getDigitVector();
    if (!arr.size()) {++events; continue;}
    //printf("Converted digits: %zu %f\n", arr.size(), arr.at(0)->getChargeFloat());

    ClusterContainer* clCont = cl->Process(arr);

    clCont->FillOutputContainer(&arrCluster);

    // ---| set cherenkov value|---------------------------------------------
    if (istr.is_open()) {
      float value=0.f;
      istr >> value;
      cherenkovValue = TMath::Abs(value);
    }
    t.Fill();

    //printf("Found clusters: %d\n", arrCluster.GetEntriesFast());

    Int_t nTries = -1; // defaults to number of clusters in the roc
    Int_t allowedMisses = 3;  //number of allowed holes in track before discarding
    Int_t minClusters = 10;
    Bool_t reverse = kFALSE; //start search from pad row opposite to beam entry side
    //TString TrackerOut = "tracksLinear_ModCommonModeEffect.root";
    //converter.TrackFinder(nTries, allowedMisses, minClusters, reverse, &arrCluster, arrTracks);

    //tOut.Fill();

    // reminder: change the initial value of bool askTZeroFine to change between raw and more fine time acceptance 
    converter.CommonModeEffectSearch(&arrCluster, converter.getDigitVector(), HistQTrackPads, HistQTrackExcluded, HistTZeroDist, HistTZeroFineDist, runNumber, converter.getPresentEventNumber(), pedestal, NActive, SignalPads, askTZeroFine, eventCounter, timeAcceptLow, timeAcceptUp, HistCmmMdEffPrfl, HistCmmMdEffPrflP, HistTimeDistAll, HistTZeroFineSlopes, fine, cherenkovValue, HistEMinusSgnl, HistEMinusCme, HistPionSgnl, HistPionCme, NActiveEMinus, NActivePions, SignalPadsEMinus, SignalPadsPions, eventCounterEMinus, eventCounterPions, vecClusterMax, acceptedEvent, vecClusterMaxFull); // ModCmmnMdEff 

    arrCluster.Clear();

    if (acceptedEvent == true){
      for (int  pos = 0; pos < vecClusterMax.size(); pos++){
        HistClusterMax->Fill(vecClusterMax[pos]);
      }
      acceptedEvent = false;
    }
    for (int pos = 0; pos < vecClusterMaxFull.size(); pos++){
      HistClusterMaxFull->Fill(vecClusterMaxFull[pos]);
    }

    ++events;
  }

  //ModCmmnMdEff; begin
  // ---------------------------------------------------------------------------- for all accepted events-----------------------------------------------------------------
  NActive = NActive * (1./8);   // divide the number gained in the for loop that fills the histograms of the signal as well as that of the Common Mode Effect by the number of timebins accepted (number of digits per pad) to get the real number of the red out and considered pads
  SignalPads = SignalPads * (1./8);
  int padsOfIROC = mapper.getPadsInIROC();
  double normFactorSgnl =       1./SignalPads;
  double normFactorCME =        1./NActive;
  double scaleFactor =         -1./padsOfIROC * 1./eventCounter; 
  double scaleFactorNoNorm =   -1./padsOfIROC * NActive * 1./eventCounter;
  bool askIfNormed = false;
  //printf("NActive is after dividing by the number of digits per pad and multiplying with (-1) %d, the number of pads in an IROC %d, the scaling factor %f \n", NActive, padsOfIROC, scaleFactor);

  // ---------------------------------------------------------------------------- for accepted electron events only ------------------------------------------------------
  NActiveEMinus = NActiveEMinus * (1./8);
  SignalPadsEMinus = SignalPadsEMinus * (1./8);
  double normFactorSgnlEMinus =       1./SignalPadsEMinus;
  double normFactorCMEEMinus =        1./NActiveEMinus;
  double scaleFactorEMinus =         -1./ padsOfIROC * 1./eventCounterEMinus;
  double scaleFactorNoNormEMinus =   -1./padsOfIROC * NActiveEMinus * 1./eventCounterEMinus;

  // ---------------------------------------------------------------------------- for accepted pion events only ----------------------------------------------------------
  NActivePions = NActivePions * (1./8);
  SignalPadsPions = SignalPadsPions * (1./8);
  double normFactorSgnlPions =       1./SignalPadsPions;
  double normFactorCMEPions =        1./NActivePions;
  double scaleFactorPions =         -1./ padsOfIROC * 1./eventCounterPions;
  double scaleFactorNoNormPions =   -1./padsOfIROC * NActivePions * 1./eventCounterPions;
  

  TString appendage;
  if (askTZeroFine == true){
    appendage = "TZeroFine_";
    if (askIfNormed == true) appendage = "TZeroFine_Normed_";
    HistQTrackPads->Scale(fine);   // divide the content of each bin by the size of the bin, that is given by the inverse of the int "fine" that defines into how much bins one timebin is divided
    HistQTrackExcluded->Scale(fine);
    HistEMinusSgnl->Scale(fine);
    HistEMinusCme->Scale(fine);
    HistPionSgnl->Scale(fine);
    HistPionCme->Scale(fine);
  }
  else {
    appendage = "";
    if (askIfNormed == true) appendage = "Normed_";
    HistQTrackPads->Rebin(fine, "");
    HistQTrackExcluded->Rebin(fine, "");
    HistEMinusSgnl->Rebin(fine, "");
    HistEMinusCme->Rebin(fine, "");
    HistPionSgnl->Rebin(fine, "");
    HistPionCme->Rebin(fine, "");
  }
  if (timeAcceptLow == -50 && timeAcceptUp == 50) appendage = "no_tZero_limits_";

  //HistQTrackPads->Scale(normFactorSgnl);
  //HistQTrackExcluded->Scale(normFactorCME);
  TH1D* HistQTrackPadsNegativeAndScaled =        new TH1D(*HistQTrackPads);
  TH1D* HistQTrackPadsNegativeAndScaledEMinus =  new TH1D(*HistEMinusSgnl);
  TH1D* HistQTrackPadsNegativeAndScaledPions =   new TH1D(*HistPionSgnl);

  if (askIfNormed == true){
    HistQTrackPads->Scale(normFactorSgnl);
    HistQTrackExcluded->Scale(normFactorCME);
    HistQTrackPadsNegativeAndScaled->Scale(scaleFactor); 
    
    HistEMinusSgnl->Scale(normFactorSgnlEMinus); 
    HistEMinusCme->Scale(normFactorCMEEMinus);
    HistQTrackPadsNegativeAndScaledEMinus->Scale(scaleFactorEMinus);

    HistPionSgnl->Scale(normFactorSgnlPions); 
    HistPionCme->Scale(normFactorCMEPions);
    HistQTrackPadsNegativeAndScaledPions->Scale(scaleFactorPions);
  }
  else {
    HistQTrackPadsNegativeAndScaled->Scale(scaleFactorNoNorm); 
    HistQTrackPadsNegativeAndScaledEMinus->Scale(scaleFactorNoNormEMinus);
    HistQTrackPadsNegativeAndScaledPions->Scale(scaleFactorNoNormPions);
  }

  TF1* fGamma4 = new TF1("fGamma4", Gamma4, -10, 20, 3);
  fGamma4->SetParNames("constant", "position", "width");
  // 160. is the peaking time [ns], 200. the sampling with [ns]
  // 160./200. is the peaking time in terms of time bins
  fGamma4->SetParameters(HistQTrackPads->GetMaximum(), 0, 160./200.);
  fGamma4->SetNpx(2000);

  TF1* fGamma4CpyEMinus = new TF1(*fGamma4);
  fGamma4CpyEMinus->SetParameters(HistEMinusSgnl->GetMaximum(), 0, 160./200);

  TF1* fGamma4CPyPions = new TF1(*fGamma4);
  fGamma4CPyPions->SetParameters(HistPionSgnl->GetMaximum(), 0, 160./200);
  

  // ---------------------------------------------------------------------------- for all accepted events-----------------------------------------------------------------
  TCanvas* c_cmeAndTack = new TCanvas(Form("CME_cmeAndTrack_%d_%s%der", runNumber, appendage.Data(), maxEvents), " ");
  TPad* Pad1 = new TPad("Pad1", "Pad1", 0.,0.,0.5,1.);
  TPad* Pad2 = new TPad("Pad2", "Pad2", 0.5,0.,1.,1.);
  c_cmeAndTack->cd();
  Pad1->Draw();
  Pad2->Draw("same");
  Pad1->SetLeftMargin(0.2);
  Pad2->SetLeftMargin(0.15);
 
  Pad1->cd();
  HistQTrackPads->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistQTrackPads->GetYaxis()->SetTitle("charge per pad per event");
  else HistQTrackPads->GetYaxis()->SetTitle("charge");
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistQTrackPads->GetXaxis()->SetRangeUser(-7.5,7.5);
  //HistQTrackPads->SetStats(0);
  HistQTrackPads->Draw("HIST");
  HistQTrackPads->Fit(fGamma4);
  fGamma4->Draw("same");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1001);

  Pad2->cd();
  //HistQTrackExcluded->GetXaxis()->SetTitle("timebins around t_{0}");
  //HistQTrackExcluded->GetYaxis()->SetTitle("charge");
  //HistQTrackExcluded->SetStats(0);
  HistQTrackPadsNegativeAndScaled->SetStats(0);
  HistQTrackPadsNegativeAndScaled->SetTitle("Charge of pads in one run excluding the trace's pads");
  HistQTrackPadsNegativeAndScaled->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistQTrackPadsNegativeAndScaled->GetYaxis()->SetTitle("charge per pad per event");
  else HistQTrackPadsNegativeAndScaled->GetYaxis()->SetTitle("charge");
  HistQTrackPadsNegativeAndScaled->SetLineColor(kRed);
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistQTrackPadsNegativeAndScaled->GetXaxis()->SetRangeUser(-7.5,7.5);
  HistQTrackPadsNegativeAndScaled->Draw("HIST");
  HistQTrackExcluded->Draw("HIST same");
  

  // ---------------------------------------------------------------------------- for accepted electron events only ------------------------------------------------------
  TCanvas* c_cmeAndTackEMinus = new TCanvas(Form("CME_cmeAndTrackEMinus_%d_%s%der", runNumber, appendage.Data(), maxEvents), " ");
  TPad* Pad1EMinus = new TPad("Pad1EMinus", "Pad1EMinus", 0.,0.,0.5,1.);
  TPad* Pad2EMinus = new TPad("Pad2EMinus", "Pad2EMinus", 0.5,0.,1.,1.);
  c_cmeAndTackEMinus->cd();
  Pad1EMinus->Draw();
  Pad2EMinus->Draw("same");
  Pad1EMinus->SetLeftMargin(0.2);
  Pad2EMinus->SetLeftMargin(0.15);
 
  Pad1EMinus->cd();
  HistEMinusSgnl->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistEMinusSgnl->GetYaxis()->SetTitle("charge per pad per event");
  else HistEMinusSgnl->GetYaxis()->SetTitle("charge");
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistEMinusSgnl->GetXaxis()->SetRangeUser(-7.5,7.5);
  //HistEMinusSgnl->SetStats(0);
  HistEMinusSgnl->Draw("HIST");
  HistEMinusSgnl->Fit(fGamma4CpyEMinus);
  fGamma4CpyEMinus->Draw("same");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1001);

  Pad2EMinus->cd();
  HistQTrackPadsNegativeAndScaledEMinus->SetStats(0);
  HistQTrackPadsNegativeAndScaledEMinus->SetTitle("Common Mode signal by electrons only in one run");
  HistQTrackPadsNegativeAndScaledEMinus->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistQTrackPadsNegativeAndScaledEMinus->GetYaxis()->SetTitle("charge per pad per event");
  else HistQTrackPadsNegativeAndScaledEMinus->GetYaxis()->SetTitle("charge");
  HistQTrackPadsNegativeAndScaledEMinus->SetLineColor(kRed);
  if (askIfNormed == true) HistQTrackPadsNegativeAndScaledEMinus->GetYaxis()->SetRangeUser(-0.9,0.1);
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistQTrackPadsNegativeAndScaledEMinus->GetXaxis()->SetRangeUser(-7.5,7.5);
  HistQTrackPadsNegativeAndScaledEMinus->Draw("HIST");
  HistEMinusCme->Draw("HIST same");


  // ---------------------------------------------------------------------------- for accepted pion events only ------------------------------------------------------
  TCanvas* c_cmeAndTackPions = new TCanvas(Form("CME_cmeAndTrackPions_%d_%s%der", runNumber, appendage.Data(), maxEvents), " ");
  //TCanvas* c_cmeAndTack = new TCanvas(Form("CME_cmeAndTrack_%d_%s%der", runNumber, appendage.Data(), maxEvents), " ");
  TPad* Pad1Pions = new TPad("Pad1Pions", "Pad1Pions", 0.,0.,0.5,1.);
  TPad* Pad2Pions = new TPad("Pad2Pions", "Pad2Pions", 0.5,0.,1.,1.);
  c_cmeAndTackPions->cd();
  Pad1Pions->Draw();
  Pad2Pions->Draw("same");
  Pad1Pions->SetLeftMargin(0.2);
  Pad2Pions->SetLeftMargin(0.15);
 
  Pad1Pions->cd();
  HistPionSgnl->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistPionSgnl->GetYaxis()->SetTitle("charge per pad per event");
  else HistPionSgnl->GetYaxis()->SetTitle("charge");
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistPionSgnl->GetXaxis()->SetRangeUser(-7.5,7.5);
  HistPionSgnl->SetStats(0);
  HistPionSgnl->Draw("HIST");
  HistPionSgnl->Fit(fGamma4CPyPions);
  fGamma4CPyPions->Draw("same");
  gStyle->SetOptStat(11);
  gStyle->SetOptFit(1001);

  Pad2Pions->cd();
  HistQTrackPadsNegativeAndScaledPions->SetStats(0);
  HistQTrackPadsNegativeAndScaledPions->SetTitle("Common Mode signal by electrons only in one run");
  HistQTrackPadsNegativeAndScaledPions->GetXaxis()->SetTitle("timebins around t_{0}");
  if (askIfNormed == true) HistQTrackPadsNegativeAndScaledPions->GetYaxis()->SetTitle("charge per pad per event");
  else HistQTrackPadsNegativeAndScaledPions->GetYaxis()->SetTitle("charge");
  HistQTrackPadsNegativeAndScaledPions->SetLineColor(kRed);
  if (timeAcceptLow == -50 && timeAcceptUp == 50) HistQTrackPadsNegativeAndScaledPions->GetXaxis()->SetRangeUser(-7.5,7.5);
  HistQTrackPadsNegativeAndScaledPions->Draw("HIST");
  HistPionCme->Draw("HIST same");


  TCanvas* c_TZeroDistribution = new TCanvas(Form("CME_TZeroDistribution_%d_%der", runNumber, maxEvents), " ");
  TPad* Pad1TZDist = new TPad("Pad1TZDist", "Pad1TZDist", 0.,0.,0.5,1.);
  TPad* Pad2TZDist = new TPad("Pad2TZDist", "Pad2TZDist", 0.5,0.,1.,1.);
  c_TZeroDistribution->cd();
  Pad1TZDist->Draw();
  Pad2TZDist->Draw("same");
  Pad1TZDist->SetLeftMargin(0.2);
  Pad2TZDist->SetLeftMargin(0.15);
 
  Pad1TZDist->cd();
  HistTZeroDist->SetStats(0);
  HistTZeroDist->GetXaxis()->SetTitle("t_{0}");
  HistTZeroDist->GetYaxis()->SetTitle("entries");
  HistTZeroDist->Draw("HIST");

  Pad2TZDist->cd();
  HistTZeroFineDist->SetStats(0);
  HistTZeroFineDist->GetXaxis()->SetTitle("t_{0}fine");
  HistTZeroFineDist->GetYaxis()->SetTitle("entries");
  HistTZeroFineDist->Draw("HIST");


  TCanvas* c_CmmnMdEffPrfl = new TCanvas(Form("CME_CmmnMdEffPrfl_%d_%s%der", runNumber, appendage.Data(), maxEvents), " ");
  //TPad* Pad1Prfl = new TPad("Pad1Prfl", "Pad1Prfl", 0.,0.,0.5,1.);
  //TPad* Pad2Prfl = new TPad("Pad2Prfl", "Pad2Prfl", 0.5,0.,1.,1.);
  //c_CmmnMdEffPrfl->cd();
  //Pad1Prfl->SetLeftMargin(0.2);
  //Pad2Prfl->SetLeftMargin(0.2);
  //Pad1Prfl->Draw();
  //Pad2Prfl->Draw("same");

  //Pad1Prfl->cd();
  //HistCmmMdEffPrfl->SetStats(0);
  //HistCmmMdEffPrfl->GetXaxis()->SetTitle("row");
  //HistCmmMdEffPrfl->GetYaxis()->SetTitle("pad");
  //HistCmmMdEffPrfl->GetZaxis()->SetTitle("charge from Common Mode Effect");
  //HistCmmMdEffPrfl->GetXaxis()->SetTitleOffset(1.7);
  //HistCmmMdEffPrfl->GetYaxis()->SetTitleOffset(1.7);
  //HistCmmMdEffPrfl->GetZaxis()->SetTitleOffset(1.7);
  //HistCmmMdEffPrfl->Draw("lego");

  //Pad2Prfl->cd();
  HistCmmMdEffPrflP->SetStats(0);
  HistCmmMdEffPrflP->GetXaxis()->SetTitle("row");
  HistCmmMdEffPrflP->GetYaxis()->SetTitle("pad");
  HistCmmMdEffPrflP->GetZaxis()->SetTitle("charge");
  HistCmmMdEffPrflP->Draw("colz");

  TCanvas* c_DigitTimeDist = new TCanvas(Form("CME_DigitTimeDist_%d_%der", runNumber, maxEvents), " ");
  HistTimeDistAll->SetStats(0);
  HistTimeDistAll->GetXaxis()->SetTitle("timebins");
  HistTimeDistAll->GetYaxis()->SetTitle("entries");
  HistTimeDistAll->Draw("colz");

  TCanvas* c_TZeroFitSlopeDist = new TCanvas(Form("CME_TZeroFitSlopeDist_%d_%der", runNumber, maxEvents), " ");
  HistTZeroFineSlopes->GetXaxis()->SetTitle("Slope");
  HistTZeroFineSlopes->GetYaxis()->SetTitle("entries");
  HistTZeroFineSlopes->Draw();

  TCanvas* c_ClusterMax = new TCanvas(Form("CME_ClusterMax_%d_%der", runNumber, maxEvents), " ");
  //gPad->SetLeftMargin(0.12);
  //gPad->SetRightMargin(0.12);
  c_ClusterMax->Divide(2);
  
  c_ClusterMax->cd(1);
  HistClusterMax->GetXaxis()->SetTitle("charge");
  HistClusterMax->GetYaxis()->SetTitle("entries");
  //gPad->SetLeftMargin(0.12);
  HistClusterMax->Draw(); 

  c_ClusterMax->cd(2);
  HistClusterMaxFull->GetXaxis()->SetTitle("charge");
  HistClusterMaxFull->GetYaxis()->SetTitle("entries");
  HistClusterMaxFull->Draw();

  SaveCanvas(c_cmeAndTack);
  SaveCanvas(c_cmeAndTackEMinus);
  SaveCanvas(c_cmeAndTackPions);
  SaveCanvas(c_TZeroDistribution);
  SaveCanvas(c_CmmnMdEffPrfl);
  SaveCanvas(c_DigitTimeDist);
  SaveCanvas(c_TZeroFitSlopeDist);
  SaveCanvas(c_ClusterMax);
  //SaveCanvas(c_test2);
 // ModCmmnMdEff; end

  fout.Write();
  fout.Close();

  TFile* plotFile = TFile::Open("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/recoTestCmmnMd/File_PlotsFromRuns.root", "UPDATE");
  plotFile->cd();

  TDirectory* EMinusDir = NULL;
  if (!(plotFile->GetListOfKeys()->Contains("EMinusPlotsDirectory"))){
        EMinusDir = plotFile->mkdir("EMinusPlotsDirectory");
  }

  TDirectory* PionDir = NULL;
  if (!(plotFile->GetListOfKeys()->Contains("PionPlotsDirectory"))){
        PionDir = plotFile->mkdir("PionPlotsDirectory");
  }

  TDirectory* AllParticlesDir = NULL;
  if (!(plotFile->GetListOfKeys()->Contains("AllParticleEventsPlotsDirectory"))){
        AllParticlesDir = plotFile->mkdir("AllParticleEventsPlotsDirectory");
  }

  TDirectory* CMEMapDir = NULL;
  if (!(plotFile->GetListOfKeys()->Contains("CommonModeEffectMapDirectory"))){
        CMEMapDir = plotFile->mkdir("CommonModeEffectMapDirectory");
  }

  gDirectory->cd("EMinusPlotsDirectory");
  HistEMinusSgnl->Write(Form("%s_Run%d", HistEMinusSgnl->GetName(), runNumber));
  HistEMinusCme->Write(Form("%s_Run%d", HistEMinusCme->GetName(), runNumber));
  HistQTrackPadsNegativeAndScaledEMinus->Write(Form("%s_Run%d", HistQTrackPadsNegativeAndScaledEMinus->GetName(), runNumber));

  gDirectory->cd("PionPlotsDirectory");
  HistPionSgnl->Write(Form("%s_Run%d", HistPionSgnl->GetName(), runNumber));
  HistPionCme->Write(Form("%s_Run%d", HistPionCme->GetName(), runNumber));
  HistQTrackPadsNegativeAndScaledPions->Write(Form("%s_Run%d", HistQTrackPadsNegativeAndScaledPions->GetName(), runNumber));
  
  gDirectory->cd("AllParticleEventsPlotsDirectory");
  HistQTrackPads->Write(Form("%s_Run%d", HistQTrackPads->GetName(), runNumber));
  HistQTrackExcluded->Write(Form("%s_Run%d", HistQTrackExcluded->GetName(), runNumber));
  HistQTrackPadsNegativeAndScaled->Write(Form("%s_Run%d", HistQTrackPadsNegativeAndScaled->GetName(), runNumber));
  
  gDirectory->cd("CommonModeEffectMapDirectory");
  HistCmmMdEffPrflP->Write(Form("%s_Run%d", HistCmmMdEffPrflP->GetName(), runNumber));
  
  plotFile->Close();

  //foutTracks.Write();
  //foutTracks.Close();
}

void RawClusterFinder::TrackFinder(int nTries, int allowedMisses, int minClusters, bool reverse, TClonesArray *pclusterArray, vector<TrackTPC> &arrTracks)
{

int loopcounter = 0;

  float zPosition = 0.;
  float fWindowPad = 4.;
  float fWindowTime = 4.;
 

  Mapper &mapper = Mapper::instance();
  TClonesArray &clusterArray = *pclusterArray;
  std::vector<int> trackNumbers;
  trackNumbers.resize(clusterArray.GetEntriesFast());

  arrTracks.clear();

  if (nTries < 0) nTries = clusterArray.GetEntriesFast();

  int nTracks = 0;

  if(clusterArray.GetEntriesFast()>=1000) {
    printf("Very large event, no tracking will be performed!!");
    return;
  }

  uint nclustersUsed = 0;

  //loop over the n (given bu nTries) clusters and use iFirst as starting point for tracking
  for(Int_t iFirst = 0; (iFirst<nTries&&iFirst<clusterArray.GetEntriesFast()); iFirst++){

    Int_t nClusters = clusterArray.GetEntriesFast();
    Int_t searchIndexArray[1000]; //array of cluster indices to be used in the tracking
    for(Int_t i = 0; i<1000;i++){
      searchIndexArray[i] = -1;
    }
    Int_t freeClusters = 0; //no of clusters not already associated to a track

    for(Int_t iCluster = iFirst; iCluster<nClusters;iCluster++) {
      Cluster *cluster = 0;
      if (reverse) cluster =  (Cluster*) clusterArray.At(clusterArray.GetEntriesFast()-1-iCluster);
      else cluster =  (Cluster*) clusterArray.At(iCluster);
      //if(cluster->getTrackNr() ==0) {
      if(trackNumbers[iCluster] ==0) {
	if (reverse) searchIndexArray[freeClusters] =clusterArray.GetEntriesFast()-1- iCluster;
	else searchIndexArray[freeClusters] = iCluster;
	freeClusters++;
      }
    }
    if(freeClusters < minClusters) continue;

    Int_t iStartCluster = -1;
    Int_t trackIndexArray[100];
    Int_t nClustersInTrackCand = 0;

    TGraph *trackCandidateGraph = new TGraph();
    TGraph *trackCandidateGraphYZ = new TGraph();
    TGraph *zGraph = new TGraph();


    for(Int_t i = 0; i<100;i++) trackIndexArray[i] = -1;

    iStartCluster = searchIndexArray[0];
    Cluster *startCluster = (Cluster*) clusterArray.At(iStartCluster);

    // Cluster *test = 0;
    //if(reverse) test =(Cluster*) clusterArray.At(clusterArray.GetEntriesFast()-1-iFirst);
    //else test = (Cluster*) clusterArray.At(iFirst);
    //if(startCluster->getTrackNr() != 0) continue;
    if(trackNumbers[iStartCluster] != 0) continue;

//     Double_t oldX = startCluster->GetRow();
    DigitPos pos(startCluster->getCRU(), PadPos(startCluster->getRow(), startCluster->getPadMean()));
    float oldRow = pos.getPadSecPos().getPadPos().getRow();
    float oldY = pos.getPadSecPos().getPadPos().getPad();
    //Double_t oldY = startCluster->GetPad();
    float oldZ = startCluster->getTimeMean();
    //Double_t oldZ = startCluster->GetTimeBinW();
    //Int_t oldRow = Int_t(startCluster->GetRow());
    const CRU cru(startCluster->getCRU());

    const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
    const int rowInSector       = startCluster->getRow() + region.getGlobalRowOffset();
    const GlobalPadNumber pad   = mapper.globalPadNumber(PadPos(rowInSector, startCluster->getPadMean()));
    const PadCentre& padCentre  = mapper.padCentre(pad);
    const float localYfactor    = (cru.side()==Side::A)?-1.f:1.f;

    LocalPosition3D clusLoc(padCentre.X(), localYfactor*padCentre.Y(), zPosition);
    trackIndexArray[0] = iStartCluster;
    nClustersInTrackCand = 1;
    Int_t totalMisses = 0;
    Int_t currentSearchRow = oldRow + 1;
    trackCandidateGraph->SetPoint(0,clusLoc.X(),clusLoc.Y());
    trackCandidateGraphYZ->SetPoint(0,clusLoc.Y(),startCluster->getTimeMean());
    zGraph->SetPoint(0,clusLoc.X(),startCluster->getTimeMean());
    Bool_t candidateFound = kFALSE;
    const Double_t searchWindow = fWindowPad;
    const Double_t timeWindow = fWindowTime;


    for(Int_t iCluster = 1; iCluster < freeClusters; iCluster++) {


      Cluster *currentCluster = (Cluster*) clusterArray.At(searchIndexArray[iCluster]);
      //DigitPos currentpos(currentCluster->getCRU(), PadPos(currentCluster->getRow(), currentCluster->getPadMean()));
      //Int_t clusterRow = currentpos.getPadSecPos().getPadPos().getRow();
      const PadRegionInfo& region = mapper.getPadRegionInfo(currentCluster->getCRU());
      Int_t clusterRow = currentCluster->getRow() + region.getGlobalRowOffset();
      if(clusterRow == oldRow) continue;
      //if(clusterRow > currentSearchRow) {
      if(clusterRow != currentSearchRow) {
        totalMisses=totalMisses+(TMath::Abs(clusterRow-oldRow)-1);
        currentSearchRow = clusterRow;
      }
      Double_t closest = 1000.;
      Double_t zClosest = 1000.;
      Int_t iClosest = -1;

      while (clusterRow == currentSearchRow && iCluster < freeClusters ) { //search currentSearchRow to find the cluster closest to the last accepted cluster


        Double_t currentZ = currentCluster -> getTimeMean();
        Double_t currentY = currentCluster -> getPadMean(); //currentpos.getPadSecPos().getPadPos().getPad();

        Double_t distance;
        distance = (TMath::Abs(oldY - currentY));

        Double_t zDistance = (TMath::Abs(oldZ - currentZ));

        if(distance < closest && zDistance<=timeWindow) {

          closest = distance;
          zClosest = zDistance;
          iClosest = iCluster;
        }
        if (searchIndexArray[iCluster+1]==-1) break;
        currentCluster = (Cluster*) clusterArray.At(searchIndexArray[iCluster+1]);
        DigitPos newcurrentpos(currentCluster->getCRU(), PadPos(currentCluster->getRow(), currentCluster->getPadMean()));

        clusterRow = Int_t(newcurrentpos.getPadSecPos().getPadPos().getRow());
        if (clusterRow == currentSearchRow) iCluster++;
      }

      Int_t miss = TMath::Abs(currentSearchRow-oldRow)-1;
      if(miss>0) totalMisses++;
      if(nClustersInTrackCand<6 && miss>1) break;//seems to lower risk for bad seed cluster being accepted
      if(miss>allowedMisses)  break;


      //Now check if cluster is close enough
      if( closest <=searchWindow && zClosest <=timeWindow){
        trackIndexArray[nClustersInTrackCand] = searchIndexArray[iClosest];
        nClustersInTrackCand++;

        Cluster *acceptedCluster = (Cluster*) clusterArray.At(searchIndexArray[iClosest]);
        DigitPos acceptedpos(currentCluster->getCRU(), PadPos(currentCluster->getRow(), currentCluster->getPadMean()));

        // 	oldX = acceptedCluster->GetRow();
        oldY = acceptedpos.getPadSecPos().getPadPos().getPad();
        oldZ = acceptedCluster->getTimeMean();
        oldRow = Int_t(acceptedpos.getPadSecPos().getPadPos().getRow());
        const CRU acceptedcru(acceptedCluster->getCRU());

        const PadRegionInfo& acceptedregion = mapper.getPadRegionInfo(acceptedcru.region());
        const int rowInSector       = acceptedCluster->getRow() + acceptedregion.getGlobalRowOffset();
        const GlobalPadNumber acceptedpad   = mapper.globalPadNumber(PadPos(rowInSector, acceptedCluster->getPadMean()));
        const PadCentre& acceptedpadCentre  = mapper.padCentre(acceptedpad);
        const float acceptedlocalYfactor    = (acceptedcru.side()==Side::A)?-1.f:1.f;

        LocalPosition3D acceptedclusLoc(acceptedpadCentre.X(), acceptedlocalYfactor*acceptedpadCentre.Y(), zPosition);
        GlobalPosition3D acceptedclusGlob = Mapper::LocalToGlobal(acceptedclusLoc, acceptedcru.sector());



        trackCandidateGraph->SetPoint(trackCandidateGraph->GetN(),acceptedclusLoc.X(),acceptedclusLoc.Y());
        trackCandidateGraphYZ->SetPoint(trackCandidateGraph->GetN(),acceptedclusLoc.Y(), acceptedCluster->getTimeMean());
        zGraph->SetPoint(zGraph->GetN(),acceptedclusLoc.X(),acceptedCluster->getTimeMean());
      }

      if(reverse) currentSearchRow--;
      else currentSearchRow++;

      if(nClustersInTrackCand>=minClusters) candidateFound = kTRUE;
    }

    if(candidateFound) {

      float slope=0., offset=0., chi2=0., slopeError=0., offsetError=0.;
      float slopeZ=0., offsetZ=0., chi2Z=0., slopeErrorZ=0., offsetErrorZ=0.;
      GetFit(trackCandidateGraph, slope, offset, slopeError, offsetError, chi2);
      GetFit(trackCandidateGraph, slope, offset, slopeError, offsetError, chi2);
      GetFit(zGraph, slopeZ, offsetZ, slopeErrorZ, offsetErrorZ, chi2Z);

     if(loopcounter == 0) {
    	TCanvas *c1 = new TCanvas(); 	       
	trackCandidateGraph->Draw("APE");
	c1->Update();
        c1->Print("test_cand_1_xy_ModCommonModeEffect.root");
        TCanvas *c2 = new TCanvas(); 	       
	trackCandidateGraphYZ->Draw("APE");
	c2->Update();
        c2->Print("test_cand_1_yz_ModCommonModeEffect.root");
        TCanvas *c3 = new TCanvas(); 	       
	zGraph->Draw("APE");
	c3->Update();
        c3->Print("test_cand_1_xz_ModCommonModeEffect.root");
      }
      if(loopcounter == 1) {
    	TCanvas *c4 = new TCanvas(); 	       
	trackCandidateGraph->Draw("APE");
	c4->Update();
        c4->Print("test_cand_2_xy_ModCommonModeEffect.root");
        TCanvas *c5 = new TCanvas(); 	       
	trackCandidateGraphYZ->Draw("APE");
	c5->Update();
        c5->Print("test_cand_2_yz_ModCommonModeEffect.root");
        TCanvas *c6 = new TCanvas(); 	       
	zGraph->Draw("APE");
	c6->Update();
        c6->Print("test_cand_2_xz_ModCommonModeEffect.root");
      }
      if(loopcounter == 2) {
        TCanvas *c7 = new TCanvas(); 	       
	trackCandidateGraph->Draw("APE");
	c7->Update();
        c7->Print("test_cand_3_xy_ModCommonModeEffect.root");
        TCanvas *c8 = new TCanvas(); 	       
	trackCandidateGraphYZ->Draw("APE");
	c8->Update();
        c8->Print("test_cand_3_yz_ModCommonModeEffect.root");
        TCanvas *c9 = new TCanvas(); 	       
	zGraph->Draw("APE");
	c9->Update();
        c9->Print("test_cand_3_xz_ModCommonModeEffect.root");
      }

      ++loopcounter;
      TrackTPC trackTPC(offset, slope, {offsetError, slopeError, chi2, offsetZ, slopeZ}, {0, 0});
      arrTracks.push_back(trackTPC);
      TrackTPC & storedTrack = arrTracks.back();

      for(Int_t i = 0; i<nClustersInTrackCand;i++) {


        Cluster *chosen = (Cluster*) clusterArray.At(trackIndexArray[i]);
        trackNumbers[trackIndexArray[i]]=nTracks+1;
        //chosen->setTrackNr(nTracks+1);


//        if(fRunResUBiased){

//          TGraph tempGraphXY(*trackCandidateGraph);
//          TGraph tempGraphXZ(*zGraph);
//          tempGraphXY.RemovePoint(i);
//          tempGraphXZ.RemovePoint(i);

//          Double_t slopeTemp=0., offsetTemp=0., chi2Temp=0., slopeErrorTemp=0., offsetErrorTemp=0.;
//          Double_t slopeZTemp=0., offsetZTemp=0., chi2ZTemp=0., slopeErrorZTemp=0., offsetErrorZTemp=0.;

//          GetFit(&tempGraphXY, slopeTemp, offsetTemp, slopeErrorTemp, offsetErrorTemp, chi2Temp);
//          GetFit(&tempGraphXZ, slopeZTemp, offsetZTemp, slopeErrorZTemp, offsetErrorZTemp, chi2ZTemp);

//          chosen->SetResidualYUnBiased(slopeTemp*chosen->GetX() + offsetTemp - chosen->GetY());

//        }

        storedTrack.addCluster(*chosen);

        ++nclustersUsed;
      }

//       trackCandidateGraph->Delete();
      nTracks++;
    }
    delete zGraph;
    delete trackCandidateGraph;
    delete trackCandidateGraphYZ;
  }
  std::cout << arrTracks.size() << " tracks found \n";
}


void RawClusterFinder::GetFit(TGraph *trackGraph, float& slope, float& offset, float& slopeError, float& offsetError, float& chi2) {
   //printf("begin of get fit");
   TLinearFitter fitter(1,"pol1","D");
   //printf("A");
   fitter.AssignData(trackGraph->GetN(), 1, trackGraph->GetX(), trackGraph->GetY());
   //printf("B");
   fitter . Eval();//Robust(0.8); //at least 80% good points, return 0 if fit is ok
   //printf("C");
   fitter . Chisquare();
   chi2 = fitter.GetChisquare();
   offset = fitter.GetParameter(0);
   slope = fitter.GetParameter(1);
   slopeError = fitter.GetParError(1);
   offsetError = fitter.GetParError(0);
   //printf("D");
   //if(degree > 1) track->SetP2(fitter.GetParameter(2));
   //if(degree > 2) track->SetP3(fitter.GetParameter(3));

//    fitter.Delete();
    //printf("end of get fit");
}

// ModCmmnMdEff; begin
void RawClusterFinder::CommonModeEffectSearch(TClonesArray* useArrCluster, std::vector<std::unique_ptr<Digit>>& useDigitVector, TH1D* useHistQTrackPads, TH1D* useHistQTrackExcluded, TH1F* useHistTZeroDist, TH1F* useHistTZeroFineDist, int useRunNr, int useEventNr, CalDet<float> *pedestal, int& useNActive, int& useSignalPads, bool& useAskTZeroFine, int& useEventCounter, const int& useTimeAcceptLow, const int& useTimeAcceptUp, TH2D* useHistCmmMdEffPrfl, TProfile2D* useHistCmmMdEffPrflP,  TH1F* useHistTimeDistAll, TH1F* useHistTZeroFineSlopes, int useFine, float useCherenkovValue, TH1D* useEMinusSgnl, TH1D* useEMinusCme, TH1D* usePionSgnl, TH1D* usePionCme, int& useNActiveEMinus, int& useNActivePions, int& useSignalPadsEMinus, int& useSignalPadsPions, int& useEventCounterEMinus, int& useEventCounterPions, std::vector<float>& useVecClusterMax, bool& useAcceptedEvent, std::vector<float>& useVecClusterMaxFull) {
  const int firstSelection = useArrCluster->GetEntriesFast();
  bool debugOutput = false;

  //if (useEventNr != 65) return;  // option for testing one definite event
  //else debugOutput = true;
  
  int arrLowerEdgeOfRows[63];  
  int arrUpperEdgeOfRows[63];

  for (int pos = 0; pos < 63; pos++){
    arrLowerEdgeOfRows[pos] = +1000;
    arrUpperEdgeOfRows[pos] = -1000;
  }

  if (useEventCounter == 0 || debugOutput == true) printf("cluster found: %d, runnumber_%d, event_%d \n", firstSelection, useRunNr, useEventNr);

  Mapper &mapper = Mapper::instance();

  std::ofstream checkChargeFile("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/recoTestCmmnMd/run000090/CheckChargeOfClusterMax.txt");

  for (int i=0; i<useDigitVector.size(); i++){
    const int   regionNr = useDigitVector.at(i)->getCRU();
    const       PadRegionInfo& region = mapper.getPadRegionInfo(regionNr);
    const int   rowOffset = region.getGlobalRowOffset();
    const int   rowInRegion = useDigitVector.at(i)->getRow();
    const int   row = rowInRegion+rowOffset; // row in sector
    const int   pad = useDigitVector.at(i)->getPad();
    const int   padForHist = useDigitVector.at(i)->getPad()-(mapper.getNumberOfPadsInRowRegion(regionNr, rowInRegion))/2;
    const float chargeF = useDigitVector.at(i)->getChargeFloat()-pedestal->getCalArray(0).getValue(row, pad);
    if (row == 4 && padForHist == 8) continue;
    if (useVecClusterMaxFull[row] < chargeF) useVecClusterMaxFull[row] = chargeF;

    //if (chargeF > 960.) {
    //TFile* checkChargeFile = TFile::Open("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/recoTestCmmnMd/run000090/CheckChargeOfClusterMax.txt", "UPDATE");
    //checkChargeFile.open("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/recoTestCmmnMd/run000090/CheckChargeOfClusterMax.txt", std::ofstream::app);
    //char* writeToTxt = (char *) &chargeF;
    //checkChargeFile << writeToTxt << "\n"; 
    //checkChargeFile.close();
    //}
  }

  //Step 1
  if (!(firstSelection > 55 && firstSelection < 70)) return;

  //Mapper &mapper = Mapper::instance();
  int tZero = -1;
  int thresholdOne = 5;
  bool saveTheCanvas = false;
  TH2F* Hist2DDigitMap = new TH2F("Hist2DDigitMap", "2D Historgram of pads with selection of only events with 55 < #Cluster < 70", 63, 0, 63, 25, 5, 30);
  TH2F* Hist2DTrackMap = new TH2F("Hist2DTrackMap", Form("2D Historgram of pads after selections, 55 < #Cluster < 70, charge > %d, Run %d, Event %d", thresholdOne,  useRunNr, useEventNr), 63, 0, 63,  25, 5, 30);
  TH1F* HistTimeDistribution = new TH1F("HistTimeDistribution", Form("Time distribution of digits, 55 < #Cluster < 70, charge > %d, Run %d, Event %d", thresholdOne, useRunNr, useEventNr), 25, 0, 25);
  TProfile* HistFindTrace = new TProfile("HistFindTrace", Form("Profile of the charge distribution, events: 55 < #Cluster < 70, Chage > 1, Run %d, Event %d", useRunNr, useEventNr), 63, 0, 63, 5, 30, ""); 
  TProfile* HistTZeroFine = new TProfile("HistTZeroFine", Form("Profile of the time distribution, events: 55 < #Cluster < 70, charge > 1, Run %d, Event %d", useRunNr, useEventNr), 63, 0, 63, 5, 30, ""); 
  TH2F* Hist2DCheck = new TH2F("Hist2DCheck", "These digits (summed up to pads) passed all cuts", 63, 0, 63, 25, 5, 30);
  TH1F* HistQTrackPadsPerEvent = new TH1F("HistQTrackPadsPerEvent", Form("Track Pads of event %d from run %d", useEventNr, useRunNr), 10*useFine, -3, 7);
  TH1F* HistQTrackExcludedPerEvent = new TH1F("HistQTrackExcludedPerEvent", Form("Pads of event %d from run %d, track pads are excluded", useEventNr, useRunNr), 10*useFine, -3, 7);
  TH2F* Hist2DTMinusTZero = new TH2F("Hist2DTMinusTZero", Form("3D charge distribution, divided by 8 as the number of digits per pad, run %d event %d", useRunNr, useEventNr), 63, 0, 63, 10*useFine, -3, 6);
  TH2F* Hist2DT           = new TH2F("Hist2DT", Form("3D charge distribution, divided by 8 as the number of digits per pad, run %d event %d", useRunNr, useEventNr), 63, 0, 63, 20, 0, 20);
  TH1F* HistSignalDigitsTimeDist = new TH1F("HistSignalDigitsTimeDist", Form("Distribution of time of digits from the signal of event %d, run %d", useEventNr, useRunNr), 20, 0, 20);
  TH1F* HistSignalTZeroFineDist = new TH1F("HistSignalTZeroFineDist", Form("Distribution of TZeroFines from the signal of event %d, run %d", useEventNr, useRunNr), 20*useFine, 0, 20);
  TH1F* HistSignalChargeVSDigitTime = new TH1F("HistSignalChargeVSDigitTime", Form("The digits charge versus the digits time, signal digits only, event %d, run %d", useEventNr, useRunNr), 20, 0, 20);

  for (int i=0; i<useDigitVector.size(); i++){
    const int   regionNr = useDigitVector.at(i)->getCRU();
    const       PadRegionInfo& region = mapper.getPadRegionInfo(regionNr);
    const int   rowOffset = region.getGlobalRowOffset();
    const int   rowInRegion = useDigitVector.at(i)->getRow();
    const int   row = rowInRegion+rowOffset; // row in sector
    const int   pad = useDigitVector.at(i)->getPad();
    const int   padForHist = useDigitVector.at(i)->getPad()-(mapper.getNumberOfPadsInRowRegion(regionNr, rowInRegion))/2;
    const float chargeF = useDigitVector.at(i)->getChargeFloat()-pedestal->getCalArray(0).getValue(row, pad);
    const int   time = useDigitVector.at(i)->getTimeStamp();

    arrLowerEdgeOfRows[row] = TMath::Min(padForHist, arrLowerEdgeOfRows[row]);
    arrUpperEdgeOfRows[row] = TMath::Max(padForHist, arrUpperEdgeOfRows[row]);

    useHistTimeDistAll->Fill(time);

    if (row == 4 && padForHist == 8) continue;

    //const Int_t bin = Hist2DDigitMap->FindBin(row, padForHist);
    //if (Hist2DDigitMap->GetBinContent(bin)) printf("Problem with bin content (%d, %d): %f\n", row, pad, Hist2DDigitMap->GetBinContent(bin));
    Hist2DDigitMap->Fill(row, padForHist, chargeF);
    
    if (time <= 4 || time >= 14) continue; 
    if (row <= 6 || row >= 56) continue;

    //Step 2
    if (chargeF > thresholdOne){
      HistTimeDistribution->Fill(time);
      Hist2DTrackMap->Fill(row, padForHist, chargeF);
      //HistFindTrace->Fill(row, padForHist, chargeF);
      //HistTZeroFine->Fill(row, time, chargeF);
      //saveTheCanvas = true;
    }
    if (chargeF > 1){
      HistFindTrace->Fill(row, padForHist, chargeF);
      HistTZeroFine->Fill(row, time, chargeF);
    }
  }

  //for (int pos = 0; pos < 62; pos++){
    //printf("row %d: lower edge is %d, upper edge is %d \n", pos, arrLowerEdgeOfRows[pos], arrUpperEdgeOfRows[pos]);
  //}
    //return;

  //Step 2a
  if (HistTimeDistribution->GetBinContent(HistTimeDistribution->GetMaximumBin()) > 20) {
    tZero = HistTimeDistribution->GetMaximumBin()-1;
  }
  else {
    delete Hist2DDigitMap;
    delete Hist2DTrackMap;
    delete HistTimeDistribution;
    delete HistFindTrace;
    delete HistTZeroFine;
    delete Hist2DCheck;
    delete Hist2DTMinusTZero;
    delete Hist2DT;
    delete HistQTrackPadsPerEvent;
    delete HistQTrackExcludedPerEvent;
    delete HistSignalDigitsTimeDist;
    delete HistSignalTZeroFineDist;
    delete HistSignalChargeVSDigitTime;
    //printf("Check mark 1\n");
    return;
  }

  TF1* FitTrace = new TF1("FitTrace", "[0]*x+[1]", Hist2DTrackMap->GetMinimum(), Hist2DTrackMap->GetMaximum()); 
  FitTrace->SetParName(0, "slope");
  FitTrace->SetParName(1, "constant");
  FitTrace->SetParameter(0,2);
  FitTrace->SetParameter(1,0);
  HistFindTrace->Fit(FitTrace);
  
  TF1* FitTZeroFine = new TF1("FitTZeroFine", "[0]*x+[1]", HistTZeroFine->GetMinimum(), HistTZeroFine->GetMaximum());   // linearer Fit richtig????
  FitTZeroFine->SetParName(0, "slopeT");
  FitTZeroFine->SetParName(1, "constantT");
  FitTZeroFine->SetParameter(0,1);
  FitTZeroFine->SetParameter(1,0);
  HistTZeroFine->Fit(FitTZeroFine);

  //printf(" Uebersicht des cut fuer Ein- und Austritt: %d low 0, %d up 0, %d low 62, %d up 62 \n", arrLowerEdgeOfRows[0], arrUpperEdgeOfRows[0], arrLowerEdgeOfRows[62], arrUpperEdgeOfRows[62]);
  //printf(" aus dem Fit: Stelle O: %.2f, Stelle 62: %.2f \n", FitTrace->Eval(0), FitTrace->Eval(62));
  if (FitTrace->Eval(0) < arrLowerEdgeOfRows[0]+2 || FitTrace->Eval(0) > arrUpperEdgeOfRows[0]-2 || FitTrace->Eval(62) < arrLowerEdgeOfRows[62]+2 || FitTrace->Eval(62) > arrUpperEdgeOfRows[62]-2) {                                    // cut to exclude events with tracks that enter and leave the observed pad area out of the acceptance borders
    delete Hist2DDigitMap;
    delete Hist2DTrackMap;
    delete HistFindTrace;
    delete HistTimeDistribution;
    delete HistTZeroFine;
    delete Hist2DCheck;
    delete Hist2DTMinusTZero;
    delete Hist2DT;
    delete HistQTrackPadsPerEvent;
    delete HistQTrackExcludedPerEvent;
    delete HistSignalDigitsTimeDist;
    delete HistSignalTZeroFineDist;
    delete HistSignalChargeVSDigitTime;
    //printf("Check mark 2\n");
    return;
  }

  for (int i = 0; i<63; i++){
    if ((TMath::Abs(FitTZeroFine->Eval(i) - HistTZeroFine->GetBinContent(i+1)) > 1.3) && HistTZeroFine->GetBinContent(i+1) > 0) {  // cut to exclude events where the time the mean charge in a row was detected does not match to the tZeroFine fit
      delete Hist2DDigitMap;
      delete HistTimeDistribution;
      delete Hist2DTrackMap;
      delete HistFindTrace;
      delete HistTZeroFine;
      delete Hist2DCheck;
      delete Hist2DTMinusTZero;
      delete Hist2DT;
      delete HistQTrackPadsPerEvent;
      delete HistQTrackExcludedPerEvent;
      delete HistSignalDigitsTimeDist;
      delete HistSignalTZeroFineDist;
      delete HistSignalChargeVSDigitTime;
      //printf("Check mark 3\n");
      return;
    }
    if (HistTZeroFine->GetBinError(i+1) > 1.65) {      // cut to exclude events where the error per row is above one
      delete Hist2DDigitMap;
      delete HistTimeDistribution;
      delete Hist2DTrackMap;
      delete HistFindTrace;
      delete HistTZeroFine;
      delete Hist2DCheck;
      delete Hist2DTMinusTZero;
      delete Hist2DT;
      delete HistQTrackPadsPerEvent;
      delete HistQTrackExcludedPerEvent;
      delete HistSignalDigitsTimeDist;
      delete HistSignalTZeroFineDist;
      delete HistSignalChargeVSDigitTime;
      //printf("Check mark 4\n");
      return;
    }
    useHistTZeroFineDist->Fill(FitTZeroFine->Eval(i));
  }
  useHistTZeroDist->Fill(tZero);
  useHistTZeroFineSlopes->Fill(FitTZeroFine->GetParameter(0));

  //Step 3
  for (int i=0; i<useDigitVector.size(); i++){
    const int   regionNr = useDigitVector.at(i)->getCRU();
    const       PadRegionInfo& region = mapper.getPadRegionInfo(regionNr);
    const int   rowOffset = region.getGlobalRowOffset();
    const int   rowInRegion = useDigitVector.at(i)->getRow();
    const int   row = rowInRegion+rowOffset; // row in sector
    const int   pad = useDigitVector.at(i)->getPad();
    const int   padForHist = useDigitVector.at(i)->getPad()-(mapper.getNumberOfPadsInRowRegion(regionNr, rowInRegion))/2;
    const float chargeF = useDigitVector.at(i)->getChargeFloat()-pedestal->getCalArray(0).getValue(row, pad);
    const int   padFromFit = FitTrace->Eval(row);
    const int   time = useDigitVector.at(i)->getTimeStamp();
    const int   deltaTime = time - tZero;
    bool        bouncer = false;
    const float tZeroFine = FitTZeroFine->Eval(row);
    const float deltaTimeTwo = time - tZeroFine;

    // define the area of accepted pads
    const int lowerAcceptanceBorder = padForHist - arrLowerEdgeOfRows[row];
    const int upperAcceptanceBorder = padForHist - arrUpperEdgeOfRows[row];
    if (lowerAcceptanceBorder < 2 || upperAcceptanceBorder > -2) bouncer = true;

    //if (bouncer == true) printf("bouncer is true \n");
    //else printf("bouncer is false \n");

    //useHistTZeroFineDist->Fill(tZeroFine);

    if (time <= 4 || time >= 14) continue; 
    if (row <= 6 || row >= 56) continue;
    if (bouncer == true){  // cut to exclude edge pads (pads out of accpetance area) to avoid bias in results caused by high noise of edge pads
      bouncer = false;
      continue;
    }

    if (row == 4 && padForHist == 8) continue;                          // excluding pad with continously to high signal values
    if (useAskTZeroFine == false && (deltaTime < useTimeAcceptLow || deltaTime > useTimeAcceptUp) && (useTimeAcceptLow != -50 || useTimeAcceptUp != 50)) continue; // cut at an accepted period of time from which the pad signals are considered only
    if (useAskTZeroFine == true && (deltaTimeTwo < useTimeAcceptLow || deltaTimeTwo > useTimeAcceptUp) && (useTimeAcceptLow != -50 || useTimeAcceptUp != 50)) continue; // cut at an accepted period of time from which the pad signals are considered only

    if(useVecClusterMax[row] < chargeF) useVecClusterMax[row] = chargeF;

    Hist2DCheck->Fill(row, padForHist, chargeF);

    if (TMath::Abs(padFromFit-padForHist) <= 2 ) {
      if (useAskTZeroFine == false) {
        useHistQTrackPads->Fill(deltaTime, chargeF);
        if (useCherenkovValue < 0.01) {
          usePionSgnl->Fill(deltaTime, chargeF);
          useSignalPadsPions += 1;
        }
        if (useCherenkovValue > 0.014) {
          useEMinusSgnl->Fill(deltaTime, chargeF);
          useSignalPadsEMinus += 1;
        }
      }
      if (useAskTZeroFine == true) {
        useHistQTrackPads->Fill(deltaTimeTwo, chargeF);
        if (useCherenkovValue < 0.01) {
          usePionSgnl->Fill(deltaTimeTwo, chargeF);
          useSignalPadsPions += 1;
        }
        if (useCherenkovValue > 0.014) {
          useEMinusSgnl->Fill(deltaTimeTwo, chargeF);
          useSignalPadsEMinus += 1;
        }
      }
      if (useAskTZeroFine == true) Hist2DTMinusTZero->Fill(row, deltaTimeTwo, chargeF);
      else Hist2DTMinusTZero->Fill(row, deltaTime, chargeF);
      Hist2DT->Fill(row, time, chargeF);
      HistQTrackPadsPerEvent->Fill(deltaTimeTwo, chargeF);
      HistSignalDigitsTimeDist->Fill(time);
      HistSignalTZeroFineDist->Fill(tZeroFine);
      HistSignalChargeVSDigitTime->Fill(time, chargeF);
      useSignalPads += 1;
    }
    else {
      if (chargeF > thresholdOne) continue;
      if (useAskTZeroFine == false) {
        useHistQTrackExcluded->Fill(deltaTime, chargeF);
        if (useCherenkovValue < 0.01) {
          usePionCme->Fill(deltaTime, chargeF);
          useNActivePions += 1;
        }
        if (useCherenkovValue > 0.014) {
          useEMinusCme->Fill(deltaTime, chargeF);
          useNActiveEMinus += 1;
        }
      }
      if (useAskTZeroFine == true) {
        useHistQTrackExcluded->Fill(deltaTimeTwo, chargeF);
        if (useCherenkovValue < 0.01) {
          usePionCme->Fill(deltaTimeTwo, chargeF);
          useNActivePions += 1;
        }
        if (useCherenkovValue > 0.014) {
          useEMinusCme->Fill(deltaTimeTwo, chargeF);
          useNActiveEMinus += 1;
        }
      }
      HistQTrackExcludedPerEvent->Fill(deltaTimeTwo, chargeF);
      //if (useEventNr == 5) Hist2DCheck->Fill(row, padForHist);
      useHistCmmMdEffPrfl  ->Fill(row, padForHist, chargeF);
      useHistCmmMdEffPrflP ->Fill(row, padForHist, chargeF);
      useNActive += 1;
    }
  }

  //printf("NActive is %d at the end of the loop for one event (event %d)\n", useNActive, useEventNr);

  bool askEventPlots = false;
  if (askEventPlots == true || useEventNr <= 0) {
    TCanvas* c_DigitMap = new TCanvas(Form("CME_DigitMap_%d_event_%d", useRunNr, useEventNr), "2D view of charge in ADCs");
    gPad->SetRightMargin(0.15);
    Hist2DDigitMap->SetStats(0);
    Hist2DDigitMap->GetZaxis()->SetTitleOffset(1.3);
    Hist2DDigitMap->GetYaxis()->SetTitle("padNr");
    Hist2DDigitMap->GetXaxis()->SetTitle("row");
    Hist2DDigitMap->GetZaxis()->SetTitle("charge");
    Hist2DDigitMap->Draw("colz");

    TCanvas* c_overViewOfTrack = new TCanvas(Form("CME_overViewOfTrack_%d_event_%d", useRunNr, useEventNr), " ");
    //gPad->SetRightMargin(0.9);
    //TPad* pad1 = new TPad("pad1","pad1",0.,0.,0.33333333,1.);
    //TPad* pad2 = new TPad("pad2","pad2",0.3333333,0.,0.66666666,1.);
    //TPad* pad3 = new TPad("pad3","pad3",0.66666,0.,1.,1.);
    TPad* pad1Trace = new TPad("pad1Trace", "pad1Trace", 0.,0.,0.5,1.);
    TPad* pad2Trace = new TPad("pad2Trace", "pad2Trace", 0.5,0.,1.,1.);
    c_overViewOfTrack->cd();
    pad1Trace->Draw();
    pad2Trace->Draw("same");
    //pad3->Draw("same");
    pad1Trace->SetRightMargin(0.2);

    //pad3->cd();
    //HistTimeDistribution->GetXaxis()->SetTitle("timebins");
    //HistTimeDistribution->GetYaxis()->SetTitle("entries");
    //HistTimeDistribution->Draw();

    pad1Trace->cd();
    Hist2DTrackMap->SetStats(0);
    Hist2DTrackMap->GetZaxis()->SetTitleOffset(1.7);
    Hist2DTrackMap->GetYaxis()->SetTitle("padNr");
    Hist2DTrackMap->GetXaxis()->SetTitle("row");
    Hist2DTrackMap->GetZaxis()->SetTitle("charge");
    Hist2DTrackMap->Draw("colz");

    pad2Trace->cd();
    HistFindTrace->SetStats(0);
    HistFindTrace->GetXaxis()->SetTitle("row");
    HistFindTrace->GetYaxis()->SetTitle("padNr");
    HistFindTrace->Draw();

    TCanvas* c_overViewOfTime = new TCanvas(Form("CME_overViewOfTime_%d_event_%d", useRunNr, useEventNr), " ");
    TPad* pad1TZeroFine = new TPad("pad1TZeroFine", "pad1TZeroFine", 0.,0.,0.5,1.);
    TPad* pad2TZeroFine = new TPad("pad2TZeroFine", "pad2TZeroFine", 0.5,0.,1.,1.);
    c_overViewOfTime->cd();
    pad1TZeroFine->Draw();
    pad2TZeroFine->Draw("same");

    pad1TZeroFine->cd();
    HistTimeDistribution->GetXaxis()->SetTitle("timebins");
    HistTimeDistribution->GetYaxis()->SetTitle("entries");
    HistTimeDistribution->SetStats(0);
    HistTimeDistribution->Draw();

    pad2TZeroFine->cd();
    HistTZeroFine->GetXaxis()->SetTitle("row");
    HistTZeroFine->GetYaxis()->SetTitle("timebins");
    HistTZeroFine->GetYaxis()->SetTitleOffset(1.5);
    HistTZeroFine->SetStats(0);
    HistTZeroFine->Draw();


    //if (useEventNr <= 1000) {
    TCanvas* c_passed = new TCanvas(Form("CME_passed_%d_event_%d", useRunNr, useEventNr), " ");
    TPad* pad1Passed = new TPad("pad1Passed", "pad1Passed", 0.,0.,0.333333,1.);
    TPad* pad2Passed = new TPad("pad2Passed", "pad2Passed", 0.333333,0.,0.666666,1.);
    TPad* pad3Passed = new TPad("pad3Passed", "pad3Passed", 0.666666,0.,1.,1.);
    c_passed->cd();
    pad1Passed->SetRightMargin(0.2);
    pad2Passed->SetLeftMargin(0.2);
    pad1Passed->SetLeftMargin(0.2);
    pad1Passed->Draw();
    pad2Passed->Draw("same");
    pad3Passed->Draw("same");

    pad1Passed->cd();
    Hist2DCheck->GetXaxis()->SetTitle("row");
    Hist2DCheck->GetYaxis()->SetTitle("padNr");
    Hist2DCheck->GetZaxis()->SetTitle("charge");
    Hist2DCheck->GetYaxis()->SetTitleOffset(1.7);
    Hist2DCheck->GetZaxis()->SetTitleOffset(1.7);
    Hist2DCheck->SetStats(0);
    Hist2DCheck->Draw("colz");

    pad2Passed->cd();
    Hist2DTMinusTZero->GetXaxis()->SetTitle("row");
    Hist2DTMinusTZero->GetXaxis()->SetTitleOffset(1.7);
    Hist2DTMinusTZero->GetYaxis()->SetTitle("t-t_{0}");
    Hist2DTMinusTZero->GetYaxis()->SetTitleOffset(1.7);
    Hist2DTMinusTZero->GetZaxis()->SetTitle("charge");
    Hist2DTMinusTZero->GetZaxis()->SetTitleOffset(1.7);
    Hist2DTMinusTZero->SetStats(0);
    if (useAskTZeroFine == false) {Hist2DTMinusTZero->Rebin(useFine, "");}
    Hist2DTMinusTZero->Draw("lego");

    pad3Passed->cd();
    Hist2DT->GetXaxis()->SetTitle("row");
    Hist2DT->GetXaxis()->SetTitleOffset(1.7);
    Hist2DT->GetYaxis()->SetTitle("timebins");
    Hist2DT->GetYaxis()->SetTitleOffset(1.7);
    Hist2DT->GetZaxis()->SetTitle("charge");
    Hist2DT->GetZaxis()->SetTitleOffset(1.7);
    Hist2DT->SetStats(0);
    Hist2DT->Draw("lego");


    TCanvas* c_cmeAndTrackPerEvent = new TCanvas(Form("CME_cmeAndTrackPerEvent_%d_event_%d", useRunNr, useEventNr), " ");
    TPad* pad1CTPE = new TPad("pad1CTPE", "pad1CTPE", 0.,0.,0.5,1.);
    TPad* pad2CTPE = new TPad("pad2CTPE", "pad2CTPE", 0.5,0.,1.,1.);
    c_cmeAndTrackPerEvent->cd();
    pad1CTPE->SetLeftMargin(0.2);
    pad2CTPE->SetLeftMargin(0.2);
    pad1CTPE->Draw();
    pad2CTPE->Draw("same");

    pad1CTPE->cd();
    HistQTrackPadsPerEvent->GetXaxis()->SetTitle("timebins around t_{0}");
    HistQTrackPadsPerEvent->GetYaxis()->SetTitle("charge");
    HistQTrackPadsPerEvent->SetStats(0);
    HistQTrackPadsPerEvent->Draw("HIST same");

    pad2CTPE->cd();
    HistQTrackExcludedPerEvent->GetXaxis()->SetTitle("timebins around t_{0}");
    HistQTrackExcludedPerEvent->GetYaxis()->SetTitle("charge");
    HistQTrackExcludedPerEvent->SetStats(0);
    HistQTrackExcludedPerEvent->Draw("HIST same");
    SaveCanvas(c_cmeAndTrackPerEvent);
      //delete Hist2DCheck;
    //}

    TCanvas* c_investigate = new TCanvas(Form("CME_investigate_%d_event_%d", useRunNr, useEventNr), " ");
    c_investigate->Divide(3);

    c_investigate->cd(1);
    HistSignalDigitsTimeDist->GetXaxis()->SetTitle("Timebins");
    HistSignalDigitsTimeDist->GetYaxis()->SetTitle("entries");
    HistSignalDigitsTimeDist->SetStats(0);
    HistSignalDigitsTimeDist->Draw();

    c_investigate->cd(2);
    HistSignalTZeroFineDist->GetXaxis()->SetTitle("t_{0}_fine");
    HistSignalTZeroFineDist->GetYaxis()->SetTitle("entries");
    HistSignalTZeroFineDist->SetStats(0);
    HistSignalTZeroFineDist->Draw();

    c_investigate->cd(3);
    HistSignalChargeVSDigitTime->GetXaxis()->SetTitle("Timebins");
    HistSignalChargeVSDigitTime->GetYaxis()->SetTitle("charge");
    HistSignalChargeVSDigitTime->SetStats(0);
    HistSignalChargeVSDigitTime->Draw("HIST");

    SaveCanvas(c_DigitMap);
    SaveCanvas(c_overViewOfTrack);
    SaveCanvas(c_overViewOfTime);
    SaveCanvas(c_passed);
    SaveCanvas(c_investigate);
    delete c_DigitMap;
    delete c_overViewOfTrack;
    delete c_overViewOfTime;
    delete c_passed;
    delete c_cmeAndTrackPerEvent;
    delete c_investigate;
  }

  delete Hist2DDigitMap;
  delete HistTimeDistribution;
  delete Hist2DTrackMap;
  delete HistFindTrace;
  delete HistTZeroFine;
  delete Hist2DCheck;
  delete Hist2DTMinusTZero;
  delete Hist2DT;
  delete HistQTrackPadsPerEvent;
  delete HistQTrackExcludedPerEvent;
  delete HistSignalDigitsTimeDist;
  delete HistSignalTZeroFineDist;
  delete HistSignalChargeVSDigitTime;
  //delete c_DigitMap;
  //delete c_overViewOfTrack;
  //delete c_overViewOfTime;
  //delete c_passed;

  useEventCounter += 1;
  if (useCherenkovValue < 0.01) useEventCounterPions += 1;
  if (useCherenkovValue > 0.014) useEventCounterEMinus += 1;
  useAcceptedEvent = true;
}
// ModCmmnMdEff; end

void SaveCanvas(TCanvas* c_in){
  c_in->SaveAs(Form("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/CmmnMdEffGraficsDirectory/%s.pdf", c_in->GetName()));
  c_in->SaveAs(Form("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/CmmnMdEffGraficsDirectory/%s.root", c_in->GetName()));
  c_in->SaveAs(Form("/lustre/nyx/alice/users/schulte/data/GEM_test_1705/CmmnMdEffGraficsDirectory/%s.png", c_in->GetName()));
}  

double Gamma4(double* x, double* par)
{
  constexpr double norm = 54.598150; //std::exp(4.);
  const double val = (x[0]-par[1])/par[2];
  const double val2 = val * val;
  double ret  = 0.0; 
  if(x[0] > par[1])
    ret = par[0] * norm * std::exp(-4. * val) * val2 * val2;
  return ret;
}

} // namespace TPC

} // namespace o2
#endif

void RawClusterFinderModCommonModeEffect(TString fileInfo, TString pedestalFile, TString outputFileName="clustersModCommonModeEffect.root", Int_t maxEvents=-1, TString cherenkovFile="cherenkov.txt", o2::TPC::RawClusterFinder::ClustererType clustererType=o2::TPC::RawClusterFinder::ClustererType::HW)
{
  using namespace o2::TPC;
  RawClusterFinder::processEvents(fileInfo, pedestalFile, outputFileName, maxEvents, cherenkovFile, clustererType);
}
