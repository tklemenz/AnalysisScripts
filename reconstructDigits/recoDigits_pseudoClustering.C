#ifndef ALICEO2_TPC_DIGITDUMPER_C_
#define ALICEO2_TPC_DIGITDUMPER_C_

#include <vector>
#include <fstream>

#include "TPCCalibration/CalibRawBase.h"
#include "TPCBase/Digit.h"
#include "TPCBase/Mapper.h"
#include "TPCBase/CRU.h"

#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TGraph.h"
#include "TF1.h"


namespace o2
{
namespace TPC
{

class DigitDumper : public CalibRawBase
{
  public:

    using vectorType = std::vector<float>;

    enum class EventStatus : char {
      Ok,         ///< one-track event with more than 30 rows with signal
      twoTracks,  ///< more than one track found
      noTrack,    ///< no track found
      shortTrack, ///< track with signal in less than 30 rows
      cosmic,     ///< cosmic event
      unknown     ///< unknown particle
    };

    enum class ParticleType : char {
      pion,       ///< particle is a pion
      electron,   ///< particle is an electron
      unknown     ///< unkown particle
    };

    /// default constructor
    DigitDumper(PadSubset padSubset = PadSubset::ROC)
      : CalibRawBase(padSubset),
        mPedestals(nullptr),
        mVectorDigits()
    {;}

    /// default destructor
    virtual ~DigitDumper() = default;

    /// fill mADCdata
    /// updateROC is not used
    Int_t updateROC(const Int_t sector, const Int_t row, const Int_t pad,
                    const Int_t timeBin, const Float_t signal) final { return 0;}

    Int_t updateCRU(const CRU& cru, const Int_t row, const Int_t pad,
                    const Int_t timeBin, const Float_t signal) final;

    void setPedestals(CalPad* pedestals) { mPedestals = pedestals; }
    static void processEvents(TString fileInfo, TString pedestalFile, TString cherenkovFile, TString outputFileName, Int_t maxEvents);
    EventStatus analyseTimeStructure(std::vector<Digit*>& DigitVector, CalDet<float> *pedestal, DigitDumper::ParticleType ParticleType,
                                     TH2F *SignalCounterHisto, TH2F *EventOccupancy, TH2F *T0Search, std::vector<TH2F*> &SinglePadTimeStructure, Int_t maxEvents);
    void fillCherenkovDist(float cherenkovValue) { mCherenkovDist->Fill(cherenkovValue); }
    std::vector<Digit*>& getDigitVector() { return mVectorDigits; }

    double Gamma4Fit(double *x, double *par);

    TGraph* getChargeGraphPi() { return mCharge_Pi; }
    TGraph* getChargeGraphEle() { return mCharge_Ele; }
    TGraph* getChargeGraphSinglePadPi() { return mChargeSinglePad_Pi; }
    TGraph* getChargeGraphSinglePadEle() { return mChargeSinglePad_Ele; }

    /// Dummy end event
    virtual void endEvent() final {}

  private:
    CalPad *mPedestals;
    std::vector<Digit*> mVectorDigits;
    std::vector<std::vector<float>> mChargeVectorVector;
    std::vector<float> mPadMappingVector;
    std::vector<std::vector<float>> mPadMappingVectorRowPad;
    std::vector<float> mPadT0Vector;

    /////////////////////////////////////////////
    /// pion histos

    TH2F *mPadOccupancy_Pi             = new TH2F("mPadOccupancy_Pi", ";row;pad;counts", 63,0,63,100,-50,50);
    TGraph *mCharge_Pi                 = new TGraph();
    TGraph *mChargeSinglePad_Pi        = new TGraph();

    TH1F *mTimeStampDist_Pi            = new TH1F("mTimeStampDist_Pi",";TimeStamp_Digit_Pi;counts",20,0,20);
    TH2F *mChargeVSTime2D_Pi           = new TH2F("mChargeVSTime2D_Pi",";#Deltat_Digit_Pi;charge;counts",40,-20,20,1100,0,1100);
    TH1F *mChargeVSTime1D_Pi           = new TH1F("mChargeVSTime1D_Pi",";#Deltat_Digit_Pi;charge",40,-20,20);

    TH1F *mTimeStampDistSinglePad_Pi   = new TH1F("mTimeStampDistSinglePad_Pi",";TimeStamp_Digit_Pi;counts",20,0,20);
    TH2F *mChargeVSTime2DSinglePad_Pi  = new TH2F("mChargeVSTime2DSinglePad_Pi",";#Deltat_Digit_Pi;charge;counts",40,-20,20,1100,0,1100);
    TH1F *mChargeVSTime1DSinglePad_Pi  = new TH1F("mChargeVSTime1DSinglePad_Pi",";#Deltat_Digit_Pi;charge",40,-20,20);

    TH1F *mChargeDistAfterCuts_Pi      = new TH1F("mChargeDistAfterCuts_Pi",";charge;counts",1100,0,1100);

    TH1F *mDeltaTDist_Pi               = new TH1F("mDeltaTDist_Pi",";#Deltat_Digit_Pi;counts",40,-20,20);
    TH1F *mT0Dist_Pi                   = new TH1F("mt0Dist_Pi",";t0_Pi;counts",20,0,20);

    TH1F *mTimeStampDistAfterCuts_Pi   = new TH1F("mTimeStampDistAfterCuts_Pi",";TimeStamp_Digit_Pi;counts",20,0,20);

    TH2F *mChargeVSTimePadWithTrack_Pi = new TH2F("mChargeVSTimePadWithTrack_Pi",";#Deltat_Digit_Pi;charge",40,-20,20,200,-20,180);
    TH2F *mChargeVSTimePadNoTrack_Pi   = new TH2F("mChargeVSTimePadNoTrack_Pi",";#Deltat_Digit_Pi;charge",40,-20,20,200,-20,180);

    TH2F *mHighT0Location_Pi           = new TH2F("mHighT0Location_Pi",";row;pad;counts",63,0,63,100,-50,50);
    TH2F *mCheck                       = new TH2F("mCheck",";row;pad;counts",63,0,63,100,-50,50);

    TH2F *mPadMapping                  = new TH2F("mPadMapping",";row;pad;counts",63,0,63,100,0,100);

    /////////////////////////////////////////////
    /// electron histos

    TH2F *mPadOccupancy_Ele             = new TH2F("mPadOccupancy_Ele", ";row;pad;counts", 63,0,63,100,-50,50);
    TGraph *mCharge_Ele                 = new TGraph();
    TGraph *mChargeSinglePad_Ele        = new TGraph();

    TH1F *mTimeStampDist_Ele            = new TH1F("mTimeStampDist_Ele",";TimeStamp_Digit_Ele;counts",20,0,20);
    TH2F *mChargeVSTime2D_Ele           = new TH2F("mChargeVSTime2D_Ele",";#Deltat_Digit_Ele;charge;counts",40,-20,20,1100,0,1100);
    TH1F *mChargeVSTime1D_Ele           = new TH1F("mChargeVSTime1D_Ele",";#Deltat_Digit_Ele;charge",40,-20,20);

    TH1F *mTimeStampDistSinglePad_Ele   = new TH1F("mTimeStampDistSinglePad_Ele",";TimeStamp_Digit_Ele;counts",20,0,20);
    TH2F *mChargeVSTime2DSinglePad_Ele  = new TH2F("mChargeVSTime2DSinglePad_Ele",";#Deltat_Digit_Ele;charge;counts",40,-20,20,1100,0,1100);
    TH1F *mChargeVSTime1DSinglePad_Ele  = new TH1F("mChargeVSTime1DSinglePad_Ele",";#Deltat_Digit_Ele;charge",40,-20,20);

    TH1F *mChargeDistAfterCuts_Ele      = new TH1F("mChargeDistAfterCuts_Ele",";charge;counts",1100,0,1100);

    TH1F *mDeltaTDist_Ele               = new TH1F("mDeltaTDist_Ele",";#Deltat_Digit_Ele;counts",40,-20,20);
    TH1F *mT0Dist_Ele                   = new TH1F("mt0Dist_Ele",";t0_Ele;counts",20,0,20);

    TH1F *mTimeStampDistAfterCuts_Ele   = new TH1F("mTimeStampDistAfterCuts_Ele",";TimeStamp_Digit_Ele;counts",20,0,20);

    TH2F *mChargeVSTimePadWithTrack_Ele = new TH2F("mChargeVSTimePadWithTrack_Ele",";#Deltat_Digit_Ele;charge",40,-20,20,200,-20,180);
    TH2F *mChargeVSTimePadNoTrack_Ele   = new TH2F("mChargeVSTimePadNoTrack_Ele",";#Deltat_Digit_Ele;charge",40,-20,20,200,-20,180);

    TH2F *mHighT0Location_Ele           = new TH2F("mHighT0Location_Ele",";row;pad;counts",63,0,63,100,-50,50);

    ///////////////////////////////////////////
    /// general


    TH2F *mCountedSignals           = new TH2F("mCountedSignals",";row;pad;counted_Signals",63,0,63,100,-50,50);
    TH2F *mHitMap                   = new TH2F("mHitMap",";row;pad",63,0,63,100,-50,50);
    TH1F *mCherenkovDist            = new TH1F("mCherenkovDist", "; Cherenkov counter signal; counts", 100,0,0.06);
    TH2F *mChargeVSTimeBinNoCut     = new TH2F("mChargeVSTimeBinNoCut_all_particles",";TimeStamp_Digit;charge;counts",20,0,20,1100,0,1100);

    int mDigitCounter           = 0;
    int mDigitCounterSinglePad  = 0;


    /// dummy reset
    void resetEvent() final { mVectorDigits.clear(); }
};

//double DigitDumper::Gamma4Fit(double *x, double *par)
//{
//  double fitval = par[0]*exp(par[1]*x)*std::pow(x,4);
//  return fitval;
//}

Int_t DigitDumper::updateCRU(const CRU& cru, const Int_t row, const Int_t pad,
                                     const Int_t timeBin, const Float_t signal)
{
  float corrSignal = signal;

  /// ===| get pedestal |========================================================
  if (mPedestals) {
    corrSignal -= mPedestals->getValue(cru, row, pad);
  }

  /// ===| get correct pad and row | ============================================
  DigitPos pos(cru, PadPos(row, pad));
  const float glrow = pos.getPadSecPos().getPadPos().getRow();
  const float glpad = pos.getPadSecPos().getPadPos().getPad();

  const float cpad  = glpad - mMapper.getNumberOfPadsInRowSector(glrow)/2;

  /// ===| add new digit |=======================================================
  mVectorDigits.emplace_back(new Digit(cru, corrSignal, row, pad, timeBin));

  /// ===| fill charge vs timeBin histo |========================================
  mChargeVSTimeBinNoCut->Fill(timeBin,corrSignal);
  //mPadOccupancy->Fill(glrow,cpad);

  /// ===| fill charge graph |===================================================
  //mCharge->SetPoint(mDigitCounter,mDigitCounter,corrSignal);
  //++mDigitCounter;

  return 1;
}

void DigitDumper::processEvents(TString fileInfo, TString pedestalFile, TString cherenkovFile, TString outputFileName, Int_t maxEvents)
{
  int twoTrackEvents = 0;
  int noTrackEvents = 0;
  int shortTrackEvents = 0;
  int goodEvents = 0;
  int cosmicEvents = 0;
  int unknownEvents = 0;
  float cherenkovValue = 0.;


  TH2F *hSignalCounterHisto[maxEvents];
  TH2F *hEventOccupancy[maxEvents];
  TH2F *hT0Search[maxEvents];
  std::vector<TH2F*> SinglePadTimeStructure;

  /// output file and container
  TFile fout(outputFileName, "recreate");
  TTree tout("events","events");

  std::vector<Digit*> vecDigits;

  tout.Branch("Digits", &vecDigits);
  tout.Branch("noTrackEvents", &noTrackEvents);
  tout.Branch("shortTrackEvents", &shortTrackEvents);
  tout.Branch("goodEvents",&goodEvents);
  tout.Branch("twoTrackEvents", &twoTrackEvents);
  tout.Branch("cosmicEvents", &cosmicEvents);

  /// input cherenkov file
  ifstream istr(cherenkovFile.Data());

  /// create converter
  DigitDumper converter;
  converter.setupContainers(fileInfo);

  /// load pedestals
  TFile f(pedestalFile);
  CalDet<float> *pedestal = nullptr;
  if (f.IsOpen() && !f.IsZombie()) {
    f.GetObject("Pedestals", pedestal);
    //converter.setPedestals(pedestal);         /// use pedestals?
  }

  /// loop over data
  int events = 0;

  DigitDumper::EventStatus eventStatus;

  std::cout << "\n\nmax events " << maxEvents << std::endl;
  CalibRawBase::ProcessStatus status;
  converter.processEvent();
  while (converter.processEvent() == CalibRawBase::ProcessStatus::Ok) {
    if (maxEvents > 0 && events >= maxEvents) break;

    DigitDumper::ParticleType particleType = DigitDumper::ParticleType::unknown;

    hSignalCounterHisto[events] = new TH2F(Form("hSignalCounter_event_%i",events),";row;pad;counted_Signals",63,0,63,100,-50,50);
    hEventOccupancy[events] = new TH2F(Form("hEventOccupancy_event_%i",events),";row;pad;accumulated_charge",63,0,63,100,-50,50);
    hT0Search[events] = new TH2F(Form("hT0Search_event_%i",events),";row;time;charge",63,0,63,1000,0,1000);

    if (istr.is_open()) {
      float value=0.f;
      istr >> value;
      cherenkovValue = TMath::Abs(value);
    }

    converter.fillCherenkovDist(cherenkovValue);

    if (cherenkovValue < 0.009) { particleType = DigitDumper::ParticleType::pion; }
    else if (cherenkovValue > 0.011) { particleType = DigitDumper::ParticleType::electron; }

    printf("========| Event %4zu %d %d %hhd |========\n", converter.getPresentEventNumber(), events, maxEvents, status);

    vecDigits = converter.getDigitVector();

    /// analyse time structure
    eventStatus = converter.analyseTimeStructure(vecDigits,pedestal,particleType,hSignalCounterHisto[events],hEventOccupancy[events],hT0Search[events],SinglePadTimeStructure,maxEvents);

    if (eventStatus == DigitDumper::EventStatus::Ok) {
      ++goodEvents;
      if (particleType == DigitDumper::ParticleType::pion) {
        hSignalCounterHisto[events]->SetTitle("good,pion");
        hEventOccupancy[events]->SetTitle("good,pion");
        hT0Search[events]->SetTitle("good,pion");
      }
      else if (particleType == DigitDumper::ParticleType::electron) {
        hSignalCounterHisto[events]->SetTitle("good,electron");
        hEventOccupancy[events]->SetTitle("good,electron");
        hT0Search[events]->SetTitle("good,electron");
      }
    }
    else if (eventStatus == DigitDumper::EventStatus::twoTracks) {
      ++twoTrackEvents;
      if (particleType == DigitDumper::ParticleType::pion) {
        hSignalCounterHisto[events]->SetTitle("twoTracks,pion");
        hEventOccupancy[events]->SetTitle("twoTracks,pion");
        hT0Search[events]->SetTitle("twoTracks,pion");
      }
      else if (particleType == DigitDumper::ParticleType::electron) {
        hSignalCounterHisto[events]->SetTitle("twoTracks,electron");
        hEventOccupancy[events]->SetTitle("twoTracks,electron");
        hT0Search[events]->SetTitle("twoTracks,electron");
      }
    }
    else if (eventStatus == DigitDumper::EventStatus::shortTrack) {
      ++shortTrackEvents;
      if (particleType == DigitDumper::ParticleType::pion) {
        hSignalCounterHisto[events]->SetTitle("shortTrack,pion");
        hEventOccupancy[events]->SetTitle("shortTrack,pion");
        hT0Search[events]->SetTitle("shortTrack,pion");
      }
      else if (particleType == DigitDumper::ParticleType::electron) {
        hSignalCounterHisto[events]->SetTitle("shortTrack,electron");
        hEventOccupancy[events]->SetTitle("shortTrack,electron");
        hT0Search[events]->SetTitle("shortTrack,electron");
      }
    }
    else if (eventStatus == DigitDumper::EventStatus::cosmic) {
      ++cosmicEvents;
      if (particleType == DigitDumper::ParticleType::pion) {
        hSignalCounterHisto[events]->SetTitle("cosmic,pion");
        hEventOccupancy[events]->SetTitle("cosmic,pion");
        hT0Search[events]->SetTitle("cosmic,pion");
      }
      else if (particleType == DigitDumper::ParticleType::electron) {
        hSignalCounterHisto[events]->SetTitle("cosmic,electron");
        hEventOccupancy[events]->SetTitle("cosmic,electron");
        hT0Search[events]->SetTitle("cosmic,electron");
      }
    }
    else if (eventStatus == DigitDumper::EventStatus::unknown) {
      ++unknownEvents;
      hSignalCounterHisto[events]->SetTitle("unknown");
      hEventOccupancy[events]->SetTitle("unknown");
      hT0Search[events]->SetTitle("unknown");
    }
    else {
      ++noTrackEvents;
      hSignalCounterHisto[events]->SetTitle("noEvent");
      hEventOccupancy[events]->SetTitle("noEvent");
      hT0Search[events]->SetTitle("noEvent");
    }


    ++events;
    converter.resetEvent();

  }

  std::cout << "noTrack: " << noTrackEvents << "\t" << "shortTrack: " << shortTrackEvents
            << "\t" << "twoTrack: " << twoTrackEvents << "\t" << "cosmic: " << cosmicEvents
            << "\t" << "good: " << goodEvents << "\t" << "unknown particle: " << unknownEvents << std::endl;
  fout.Write();
  fout.WriteObject(converter.getChargeGraphPi(), "mChargeAllPads_Pi");
  fout.WriteObject(converter.getChargeGraphSinglePadPi(), "mChargeSinglePad_Pi");
  fout.WriteObject(converter.getChargeGraphEle(), "mChargeAllPads_Ele");
  fout.WriteObject(converter.getChargeGraphSinglePadEle(), "mChargeSinglePad_Ele");
  ///if (maxEvents == 1) {
  ///  for (int i=0; i<SinglePadTimeStructure.size(); ++i) {
  ///    fout.WriteObject(SinglePadTimeStructure[i], Form("hSinglePadTimeStructure_%i",i));
  ///  }
  ///}
  for (int i=0; i<maxEvents; ++i) {
    fout.WriteObject(hSignalCounterHisto[i], Form("hSignalCounter_event_%i",i));
    fout.WriteObject(hEventOccupancy[i], Form("hEventOccupancy_event_%i",i));
    fout.WriteObject(hT0Search[i], Form("hhT0Search_event_%i",i));
  }

  fout.Close();

  //std::cout << std::endl << SinglePadTimeStructure.size() << std::endl;

}

DigitDumper::EventStatus DigitDumper::analyseTimeStructure(std::vector<Digit*>& DigitVector, CalDet<float> *pedestal, DigitDumper::ParticleType particle, TH2F *SignalCounterHisto, TH2F *EventOccupancy, TH2F *T0Search, std::vector<TH2F*> &SinglePadTimeStructure, Int_t maxEvents)
{
  /// TODO: pro pad charge summieren und als occupancy zeichnen (done)
  /// Ziel: nur events mit einem Track nehmen (done)
  ///       pro event (fuer jedes pad?) das richtige t0 finden (more or less done)
  /// -> Gamma4 Signal ist auf jedem Pad, das getroffen wird (funktioniert, wenn man macro ueber ein event laufen laesst)
  /// -> je ein Pad anschauen, das getroffen wird und eins, das nicht getroffen wird (done)
  /// TODO: Gamma4 funktion an charge vs time fitten
  /// TODO: charge vs time als Profilfunktion in 2D als Funktion von cpad und als Funktion von row
  /// TODO: cpad vs row mit cut Delta t > 3 um zu sehen, wo die hohen Delta t herkommen

  EventStatus status = EventStatus::Ok;

  if (particle == ParticleType::unknown) {
    status = EventStatus::unknown;
    return status;
  }

  Mapper &mapper = Mapper::instance();

  /// fit function for t0 fit and "tracking"
  /// _____________________________________________________________________________________________________

  TF1* FitT0 = new TF1("FitT0", "[0]*x+[1]",0,63);
  FitT0->SetParName(0,"slopeT0");
  FitT0->SetParName(1,"constantT0");
  FitT0->SetParameter(0,1);
  FitT0->SetParameter(1,0);

  TF1* FitPadOcc = new TF1("FitPadOcc", "[0]*x+[1]",0,63);
  FitPadOcc->SetParName(0,"slopeTrack");
  FitPadOcc->SetParName(1,"constantTrack");
  FitPadOcc->SetParameter(0,1);
  FitPadOcc->SetParameter(1,0);

  /// _____________________________________________________________________________________________________

  int chargeThreshold = 15;
  int timeLow = 4;
  int timeHigh = 14;
  int rowLeft = 6;
  int rowRight = 57;

  int singleRow = 12;
  int singlePad = 5;


  /// prepare mChargeVectorVector (find out how many pads there are) and exclude left and right rows
  /// _____________________________________________________________________________
  int nPads = 0;
  for (int i=0; i<DigitVector.size(); ++i) {

    CRU cruObj = DigitVector.at(i)->getCRU();
    const int crurow = DigitVector.at(i)->getRow();
    const int crupad = DigitVector.at(i)->getPad();
    DigitPos pos(cruObj, PadPos(crurow, crupad));
    const float row = pos.getPadSecPos().getPadPos().getRow();
    const float pad = pos.getPadSecPos().getPadPos().getPad();
    const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;

    const int timeStamp = DigitVector.at(i)->getTimeStamp();

    mHitMap->Fill(row,cpad);    // for edge checking

    if (row <= rowLeft || row >= rowRight) { continue; }

    if ( timeStamp == 0 ) { ++nPads; }
  }

  mChargeVectorVector.resize(nPads);
  mPadMappingVectorRowPad.resize(nPads);
  mHitMap->Fill(4,-1);
  mHitMap->Fill(58,-14);
  mHitMap->Fill(59,-15);

  ///TH2F *hSinglePadTimeStructure[nPads];


  /// fill mChargeVectorVector with raw data in the format
  /// [pad0(timeBin0, timeBin1,timeBin3,...),pad1(...),pad2(...),pad3(...),...]
  /// exclude left and right rows
  /// _____________________________________________________________________________
  int padCounter = -1;
  for (int i=0; i<DigitVector.size(); ++i) {

    CRU cruObj = DigitVector.at(i)->getCRU();
    const int crurow = DigitVector.at(i)->getRow();
    const int crupad = DigitVector.at(i)->getPad();
    DigitPos pos(cruObj, PadPos(crurow, crupad));
    const float row = pos.getPadSecPos().getPadPos().getRow();
    const float pad = pos.getPadSecPos().getPadPos().getPad();
    const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;

    const int timeStamp = DigitVector.at(i)->getTimeStamp();
    const float charge = DigitVector.at(i)->getChargeFloat() - pedestal->getCalArray(0).getValue(row,pad);

    const int roc = cruObj.roc();
    const GlobalPadNumber padInROC = mMapper.getPadNumberInROC(PadROCPos(roc, row, pad));

    if (row <= rowLeft || row >= rowRight) { continue; }
    if ( timeStamp == 0 ) {
      mPadMappingVector.push_back(padInROC);          // i-ter Eintrag enthaelt globale Padnummer
      ++padCounter;
    }

    mChargeVectorVector[padCounter].push_back(charge);
    mPadMappingVectorRowPad[padCounter].push_back(row);
    mPadMappingVectorRowPad[padCounter].push_back(pad);

    if (mPadMapping->GetBinContent(int(row),int(pad)) == 0) {             // sehr fraglich!!!
      mPadMapping->Fill(row,pad,i);                                       // sollte zeigen, in welcher Reihenfolge die pads im DigitVector durchgegangen werden
    }


    /// check how the data is filled into the DigitVector and mChargeVectorVector
    //std::cout << "pad: " << padInROC << "\t " << "timeStamp: " << timeStamp << "\t" << "charge: " << charge
    //          << "      \t" << "padCounter: " << padCounter << "\t" << "nPads: " << nPads << "\t" << "timeStamp: "
    //          << mChargeVectorVector[padCounter].size() << std::endl;

  }



  /// check EventStatus
  /// _____________________________________________________________________________
  int signalPerRow[63];
  int signalPerRowMultTrack[63];
  float signalPad = -999;
  //int signalThreshold = 30;             /// wieder einfuegen!!!
  //int signalThreshold = 512;
  for (int i=0; i<63; ++i) {
    signalPerRow[i] = 0;
    signalPerRowMultTrack[i] = 0;
  }
  for (int i=0; i<DigitVector.size(); ++i) {
    CRU cruObj = DigitVector.at(i)->getCRU();
    const int crurow = DigitVector.at(i)->getRow();
    const int crupad = DigitVector.at(i)->getPad();
    DigitPos pos(cruObj, PadPos(crurow, crupad));
    const int row = pos.getPadSecPos().getPadPos().getRow();
    const float pad = pos.getPadSecPos().getPadPos().getPad();
    const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;

    const int timeStamp = DigitVector.at(i)->getTimeStamp();
    const float charge = DigitVector.at(i)->getChargeFloat() - pedestal->getCalArray(0).getValue(row,pad);
    const int roc = cruObj.roc();                                                         // not needed right now
    const GlobalPadNumber padInROC = mMapper.getPadNumberInROC(PadROCPos(roc, row, pad)); // not needed right now


    //if (charge < 500 && charge > 0) {
      if (particle == ParticleType::pion) { mPadOccupancy_Pi->Fill(row,cpad,charge); }
      else if (particle == ParticleType::electron) { mPadOccupancy_Ele->Fill(row,cpad,charge); }
      EventOccupancy->Fill(row,cpad,charge);
   // }

    if (charge >= chargeThreshold) {
      if (padInROC != signalPad && (row != 4 && cpad != -1)) {
        ++signalPerRowMultTrack[row];      // number of pads with signal per row   for multiple-track events where one track is at the edge
      }
      bool edge = false;
      if (mHitMap->GetBinContent(row,cpad+51) == 0 || mHitMap->GetBinContent(row,cpad+51-1) == 0 || mHitMap->GetBinContent(row,cpad+51-2) == 0
          || mHitMap->GetBinContent(row,cpad+51+1) == 0 || mHitMap->GetBinContent(row,cpad+51+2) == 0) {
        edge = true;
      }
      if (!edge) {
        if (padInROC != signalPad && (row != 4 && cpad != -1)) {
          ++signalPerRow[row];             // number of pads with signal per row
          mCountedSignals->Fill(row,cpad);
          SignalCounterHisto->Fill(row,cpad);
          signalPad = padInROC;
        }

        T0Search->Fill(row,timeStamp,charge);

      }
    }
  }
  //int multiSignalCounter = 0;
  int multiSignalCounterMultTrack = 0;
  int signalCounter = 0;
  for (int i=0; i<63; ++i) {
    if (signalPerRowMultTrack[i] > 1 && signalPerRowMultTrack[i] < 3) { ++multiSignalCounterMultTrack; }
    //if (signalPerRow[i] > 1 && signalPerRow[i] < 3) { ++multiSignalCounter; }
    /*else*/ if (signalPerRow[i] == 1) { ++signalCounter; }
    else if (signalPerRow[i] >= 3) {
      status = EventStatus::cosmic;
      for (int i=0; i<nPads; ++i) { mChargeVectorVector[i].clear(); }
      return status;
    }
  }

  if (multiSignalCounterMultTrack >= 10) {
    status = EventStatus::twoTracks;
    for (int i=0; i<nPads; ++i) { mChargeVectorVector[i].clear(); }
    return status;
  }

  if (signalCounter < 31 && signalCounter > 5) {
    status = EventStatus::shortTrack;
    for (int i=0; i<nPads; ++i) { mChargeVectorVector[i].clear(); }
    return status;
  }
  else if (signalCounter <= 5) {
    status = EventStatus::noTrack;
    for (int i=0; i<nPads; ++i) { mChargeVectorVector[i].clear(); }
    return status;
  }






  T0Search->Fit(FitT0);

  for (int i=0; i<nPads; ++i) {
    TString histoname;
    histoname = "";
    histoname += i;
    ///hSinglePadTimeStructure[i] = new TH2F(histoname.Data(),";timeStamp;charge",40,-20,20,500,-20,480);
  }

  ///int t0Pad = -1;                //with one t0 per event
  ///int checkForFirstTimeBin = -1;
  /// loop over data and fill histos for electrons and pions seperately
  /// _____________________________________________________________________________
  TString hName;
  TString addPadNumber;
  TString addNumberInROC;
  padCounter = -1;
  for (int i=0; i<DigitVector.size(); ++i) {
    addPadNumber = "";
    addNumberInROC= "";
    CRU cruObj = DigitVector.at(i)->getCRU();
    const int crurow = DigitVector.at(i)->getRow();
    const int crupad = DigitVector.at(i)->getPad();
    DigitPos pos(cruObj, PadPos(crurow, crupad));
    const float row = pos.getPadSecPos().getPadPos().getRow();
    if (row <= rowLeft || row >= rowRight) { continue; }
    const float pad = pos.getPadSecPos().getPadPos().getPad();
    const float cpad  = pad - mapper.getNumberOfPadsInRowSector(row)/2;
    const float chargeUncorrected = DigitVector.at(i)->getChargeFloat();
    const float charge = chargeUncorrected - pedestal->getCalArray(0).getValue(row,pad);
    const int timeStamp = DigitVector.at(i)->getTimeStamp();
    const int roc = cruObj.roc();
    const GlobalPadNumber padInROC = mMapper.getPadNumberInROC(PadROCPos(roc, row, pad));
    addNumberInROC += padInROC;
    if (timeStamp == 0) {
      ++padCounter;
      addPadNumber += padCounter;
        hName = "hSinglePadTimeStructure_"+addPadNumber+"_"+addNumberInROC;
    }

    const float deltaT = timeStamp - FitT0->Eval(row);


    //_________________________________________________________________________________________________________________
    if (particle == ParticleType::pion) {

      mTimeStampDist_Pi->Fill(timeStamp);
      //if (charge < 500 && charge > 0) { mPadOccupancy->Fill(row,cpad,charge); }

      if (charge >= chargeThreshold && (timeStamp >= timeLow && timeStamp <= timeHigh) && (row >= rowLeft && row <= rowRight)) {
        mChargeVSTime2D_Pi->Fill(deltaT,charge);
        mChargeVSTime1D_Pi->Fill(deltaT,charge);
        mChargeDistAfterCuts_Pi->Fill(charge);
        mCharge_Pi->SetPoint(mDigitCounter,mDigitCounter,charge);
        mTimeStampDistAfterCuts_Pi->Fill(timeStamp);
        mDeltaTDist_Pi->Fill(deltaT);
        ++mDigitCounter;
        if (row == singleRow && cpad == singlePad) {
          mTimeStampDistSinglePad_Pi->Fill(timeStamp);
          mChargeVSTime2DSinglePad_Pi->Fill(deltaT,charge);
          mChargeVSTime1DSinglePad_Pi->Fill(deltaT,charge);
          mChargeSinglePad_Pi->SetPoint(mDigitCounterSinglePad,mDigitCounterSinglePad,charge);
          ++mDigitCounterSinglePad;
        }
      }
      if (row == 26 && cpad == -4) { mChargeVSTimePadWithTrack_Pi->Fill(deltaT,charge); }
      if (row == 26 && cpad == -3) { mChargeVSTimePadNoTrack_Pi->Fill(deltaT,charge); }
    }

    //_________________________________________________________________________________________________________________
    else if (particle == ParticleType::electron) {

      ///if (padInROC == t0Pad && checkForFirstTimeBin == 0) { mHighT0Location_Ele->Fill(row,cpad); }                //with one t0 per event
      ///else { mCheck->Fill(row,cpad); }                //with one t0 per event
      ///++checkForFirstTimeBin;                //with one t0 per event

      mTimeStampDist_Ele->Fill(timeStamp);
      //if (charge < 500 && charge > 0) { mPadOccupancy->Fill(row,cpad,charge); }

      if (charge >= chargeThreshold && (timeStamp >= timeLow && timeStamp <= timeHigh) && (row >= rowLeft && row <= rowRight)) {
        mChargeVSTime2D_Ele->Fill(deltaT,charge);
        mChargeVSTime1D_Ele->Fill(deltaT,charge);
        mChargeDistAfterCuts_Ele->Fill(charge);
        mCharge_Ele->SetPoint(mDigitCounter,mDigitCounter,charge);
        mTimeStampDistAfterCuts_Ele->Fill(timeStamp);
        mDeltaTDist_Ele->Fill(deltaT);
        ++mDigitCounter;
        if (row == singleRow && cpad == singlePad) {
          mTimeStampDistSinglePad_Ele->Fill(timeStamp);
          mChargeVSTime2DSinglePad_Ele->Fill(deltaT,charge);
          mChargeVSTime1DSinglePad_Ele->Fill(deltaT,charge);
          mChargeSinglePad_Ele->SetPoint(mDigitCounterSinglePad,mDigitCounterSinglePad,charge);
          ++mDigitCounterSinglePad;
        }
      }
      if (row == 26 && cpad == -4) { mChargeVSTimePadWithTrack_Ele->Fill(deltaT,charge); }
      if (row == 26 && cpad == -3) { mChargeVSTimePadNoTrack_Ele->Fill(deltaT,charge); }
    }
  }

  /// write SinglePadTimeStructure histograms to the vector
  /// _____________________________________________________________________________
  ///for (int i=0; i<nPads; ++i) {
  ///  SinglePadTimeStructure.push_back(hSinglePadTimeStructure[i]);
  ///}

  /// clear mChargeVectorVector
  /// _____________________________________________________________________________
  for (int i=0; i<nPads; ++i) { mChargeVectorVector[i].clear(); }
  mPadMappingVector.clear();
  mPadMappingVectorRowPad.clear();
  mPadT0Vector.clear();

//  TF1 *func = new TF1("Gamma4Fit","[0]*exp([1]*x)*std::pow(x,4)",-10,20);
//  func->SetParNames("Constant","Exponent");
//  mChargeVSTime1D_Pi->Fit("Gamma4Fit","r");


  return status;

}

} /// end namespace TPC
} /// end namespace o2
#endif

void recoDigits_pseudoClustering(TString fileInfo, TString pedestalFile, TString cherenkovFile, TString outputFileName, Int_t maxEvents)
{
  o2::TPC::DigitDumper::processEvents(fileInfo, pedestalFile, cherenkovFile, outputFileName, maxEvents);
}
