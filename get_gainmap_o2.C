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

#include <boost/lambda/lambda.hpp>

using namespace std;
using namespace o2::TPC;

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



void GetBinMinMax(const TH1 *hist, const float frac, int &bin1, int &bin2)
{
  const int binMax=hist->GetMaximumBin();
  const double contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<frac*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<frac*contMax) break;
}



void TruncatedMean(TH1F *his, std::vector<float> *param, float down, float up)
{
  int nbins = his->GetNbinsX();
  float nentries = his-> GetEntries();
  float sum = 0;
  float mean = 0;
  float sigma2 = 0;
  float ncumul = 0;
  for (int ibin = 1; ibin<nbins; ibin++){
    ncumul += his->GetBinContent(ibin);
    float fraction = float(ncumul)/float(nentries);
    if (fraction>down && fraction<up){
      sum+=his->GetBinContent(ibin);
      mean+=his->GetBinCenter(ibin)*his->GetBinContent(ibin);
      sigma2+=his->GetBinCenter(ibin)*his->GetBinCenter(ibin)*his->GetBinContent(ibin);
    }
  }
  mean /= sum;
  sigma2 = TMath::Sqrt(TMath::Abs(sigma2/sum-mean*mean));
  if (param){
    (*param)[0] = his->GetMaximum();
    (*param)[1] = mean;
    (*param)[2] = sigma2;
  }
}

float unbiasedTruncatedMean(TH2D *histo, int runNr, o2::TPC::TrackTPC Track, float low=.0, float high=.7)
{
  std::vector<Cluster> Clusters;
  Track.getClusterVector(Clusters);
  std::vector<float> values;
  Mapper &mapper = Mapper::instance();

  for (auto &clusterObject : Clusters) {
    DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
    float row = pos.getPadSecPos().getPadPos().getRow();
    float pad = pos.getPadSecPos().getPadPos().getPad();
    if (row == 31 || row == 32) {
      histo->Fill(row, pad, 1);
      continue;
    }
    if (runNr <= 243) {
      if (row >= 57 || row < 2) {
        histo->Fill(row, pad, 1);
        continue;
      }
    }
    else if (runNr > 255) {
      if (row == 62 || row == 0) {
        histo->Fill(row, pad, 1);
        continue;
      }
      else if (row > 0 && row < 62) {
        if (pad == mapper.getNumberOfPadsInRowSector(row)/2 - 1) {
          histo->Fill(row,pad, 1);
          continue;
        }
      }
    }
    values.push_back(clusterObject.getQ());
  }
  transform(values.begin(), values.end(), values.begin(), boost::lambda::_1 * cos(atan(Track.getTgl())) * cos(asin(Track.getSnp())));
  std::sort(values.begin(), values.end());

  float dEdx = 0.f;
  int nClustersTrunc = 0;
  int nClustersUsed = static_cast<int>(values.size());

  for (int icl=0; icl<nClustersUsed; ++icl) {
    if (icl<std::round(low*nClustersUsed)) continue;
    if (icl>std::round(high*nClustersUsed)) break;

    dEdx+=values[icl];
    ++nClustersTrunc;
  }

  if (nClustersTrunc>0){
    dEdx/=nClustersTrunc;
  }
  return dEdx;
}



//int get_gainmap_o2(int argc, char **argv){
  int get_gainmap_o2(TString FileList, const char *OutputPath="/home/tom") {
  //const char *sarg = argv[1];
  //const char *spath = argv[2];

  /* ============================================================
   * ====================== Set up chain ========================
   * ============================================================ */

  EventHeader Header;
  Mapper &mapper = Mapper::instance();

  //TString FileList(sarg);

  std::vector<TrackTPC> *fChainEvent=0;

  TChain fChain("events");

  fChain.SetBranchAddress("Tracks",&fChainEvent);
  fChain.SetBranchAddress("header",&Header);

  TString allFiles;
  if (FileList.EndsWith(".txt")){
    allFiles=gSystem->GetFromPipe(Form("cat %s",FileList.Data()));
  }
  else {
  allFiles=gSystem->GetFromPipe(Form("ls %s",FileList.Data()));
  }
  TObjArray *arr = allFiles.Tokenize("\n");

  for (int ifile=0; ifile<arr->GetEntriesFast(); ++ifile){
  TString file=arr->At(ifile)->GetName();
  //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
  fChain.Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
  cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }


  /* ============================================================
   * ====================== Cut parameter =======================
   * ============================================================ */

  int TrPerEv = 1;
  int nclCut = 32;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMeanLow = 5;
  int timeMeanHigh =20;


  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/

  TH1D *hChisquarePi        = new TH1D("hChisquarePi", "; Chisquare Landau; # counts", 400,0,200);
  TH2D *hGainMapPi          = new TH2D("hGainMapPi", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncPi     = new TH2D("hGainMapTruncPi", "; Row; Pad", 63,0,63,100,0,100);
  TH1D *hEntriesPi          = new TH1D("EntriesPi", "; Entries; # counts",300,1,300000);
  TH1D *hChisquareEle       = new TH1D("hChisquareEle", "; Chisquare Landau; # counts", 400,0,200);
  TH2D *hGainMapEle         = new TH2D("hGainMapEle", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncEle    = new TH2D("hGainMapTruncEle", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapHisto       = new TH2D("GainMap", "; Row; Pad", 63,0,63,100,-50,50);
  TH1D *hEntriesEle         = new TH1D("EntriesEle", "; Entries; # counts",300,1,300000);
  TH1D *hGainDist           = new TH1D("GainDist", "; normalized gain; ",100,0,2);
  TH2D *hCheckExclude       = new TH2D("CheckExclude", "; Row; Pad", 63,0,63,100,0,100);


  /* initialize array of histograms used to get the gainmap (for pions and electrons seperately)*/
  const int nrows = 63;
  const int npads = 100;
  TH1F *hQtotPadsPi[nrows][npads];
  TH1F *hQtotPadsTruncPi[nrows][npads];

  for(int irow=0; irow<nrows; ++irow) {
    for (int ipad=0; ipad<npads; ++ipad) {
    hQtotPadsPi[irow][ipad] = new TH1F(Form("hQtotPadPi_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
    hQtotPadsTruncPi[irow][ipad] = new TH1F(Form("hQtotPadPi_trunc_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
    }
  }

  TH1F *hQtotPadsEle[nrows][npads];
  TH1F *hQtotPadsTruncEle[nrows][npads];

  for(int irow=0; irow<nrows; ++irow) {
    for (int ipad=0; ipad<npads; ++ipad) {
    hQtotPadsEle[irow][ipad] = new TH1F(Form("hQtotPadEle_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
    hQtotPadsTruncEle[irow][ipad] = new TH1F(Form("hQtotPadEle_trunc_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; #counts", 500, 0, 50);
    }
  }

  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  int pions = 0;
  int electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  float CherenkovValue = 0;
  int runNr = 0;
  bool isok = true;
  std::vector<TrackTPC> *vecEvent = 0;

  //TObjArray *chainEntries = fChain->GetListOfFiles();


  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;

  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
  TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
  cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");
    //TTree *tree = (TTree*)chainEntries->At(ifile)->Get("events");


    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);
    //cout<<"KOMME BIS HIER!"<<endl;

    CherenkovValue = 0;
    runNr = 0;


    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);

      int NTracks = vecEvent->size();
      CherenkovValue = Header.cherenkovValue;
      runNr = Header.run;

      Tracks += NTracks;



  /*========================================== CUT ================================================*/
      if (NTracks != TrPerEv) continue;																// only one-track events


      ++OneTrackEvents;
      //cout<<"DEBUG POINT 4"<<endl;

    for (auto trackObject : *vecEvent) {
      std::vector<Cluster> clCont;
      trackObject.getClusterVector(clCont);


      int ncl = clCont.size();
      isok = true;

      onetrcl += ncl;


      //cout<<"DEBUG POINT 5"<<endl;


  /*========================================== CUT ================================================*/
      if (ncl < nclCut) continue;																	// cut on number of clusters per track

  /*========================================== CUT ================================================*/
      if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov


      for (auto &clusterObject : clCont) {															// make cuts on cluster properties
        ushort timeMean = clusterObject.getTimeMean();
  /*========================================== CUT ================================================*/
        if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle
      }

      if (isok == true){																			// track accepted
        usedcl += ncl;
        //float dEdx = trackObject.getTruncatedMean(0.,.7,1,true);
        float dEdx = unbiasedTruncatedMean(hCheckExclude, runNr, trackObject);

        if (CherenkovValue < CherCutLow){
              for (auto &clusterObject : clCont) {															// loop over clusters
                DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
                int row = pos.getPadSecPos().getPadPos().getRow();
                int pad = pos.getPadSecPos().getPadPos().getPad();
                float QTot = clusterObject.getQ();

                hQtotPadsPi[row][pad]->Fill(QTot/dEdx);
                hQtotPadsTruncPi[row][pad]->Fill(QTot/dEdx);
              }
              ++pions;

        }
        if (CherenkovValue > CherCutHigh){
              for (auto &clusterObject : clCont) {															// loop over clusters
                DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
                int row = pos.getPadSecPos().getPadPos().getRow();
                int pad = pos.getPadSecPos().getPadPos().getPad();
                float QTot = clusterObject.getQ();

                hQtotPadsEle[row][pad]->Fill(QTot/dEdx);
                hQtotPadsTruncEle[row][pad]->Fill(QTot/dEdx);
              }
              ++electrons;

        }

      }
    } ///break; ///if (isok==true){break;} //use only one track
    }
    delete TreeFile;
    cout<<"Tree deleted"<<endl<<endl;

  }

  /* ======================================================================
   * ========== Fit every histogram in array for MPV -> GainMap ===========
   * ====================================================================== */

  /*TFile *fPi = TFile::Open(Form("%s/GainMapPi_Landau.root",spath),"RECREATE");
  TFile *fEle = TFile::Open(Form("%s/GainMapEle_Landau.root",spath),"RECREATE");
  TFile *chiPi = TFile::Open(Form("%s/ChiSquarePi.root",spath),"RECREATE");
  TFile *chiEle = TFile::Open(Form("%s/ChiSquareEle.root",spath),"RECREATE");
  TFile *entPi = TFile::Open(Form("%s/EntriesPi.root",spath),"RECREATE");
  TFile *entEle = TFile::Open(Form("%s/EntriesEle.root",spath),"RECREATE");
  TFile *gTrPi = TFile::Open(Form("%s/GainMapPi_Trunc.root",spath),"RECREATE");
  TFile *gTrEle = TFile::Open(Form("%s/GainMapEle_Trunc.root",spath),"RECREATE");*/


//  TFile *fPi = TFile::Open("GainMapPi_Landau.root","RECREATE");
//  TFile *fEle = TFile::Open("GainMapEle_Landau.root","RECREATE");
//  TFile *chiPi = TFile::Open("ChiSquarePi.root","RECREATE");
//  TFile *chiEle = TFile::Open("ChiSquareEle.root","RECREATE");
//  TFile *entPi = TFile::Open("EntriesPi.root","RECREATE");
//  TFile *entEle = TFile::Open("EntriesEle.root","RECREATE");
//  TFile *gTrPi = TFile::Open("GainMapPi_Histo.root","RECREATE");
//  TFile *gTrEle = TFile::Open("GainMapEle_Histo.root","RECREATE");
//  TFile *gGM = TFile::Open("GainMap_Histo.root", "recreate");

  ROC roc(0);
  CalDet<float> *gainmapEle = new CalDet<float>(PadSubset::ROC);
  gainmapEle->setName("GainMapEle");
  CalArray<float>& GainMapEle = gainmapEle->getCalArray(roc);

  CalDet<float> *gainmapPi = new CalDet<float>(PadSubset::ROC);
  gainmapPi->setName("GainMapPi");
  CalArray<float>& GainMapPi = gainmapPi->getCalArray(roc);

  TF1 *MPVfit = new TF1("MPV","landau");
  int means = 0;
  int Landau = 0;
  const Float_t frac=0.2;
  Int_t bin1=0,bin2=0;
  //TVectorD param(4);
  std::vector<float> param(4);
  float lowTr = 0.;
  float upTr = .75;
  for (int irow=0; irow<nrows; ++irow){
    for (int ipad=0; ipad<npads; ++ipad){
    /*if (irow == 62 && ipad == 66){
    hQtotPadsTruncPi[irow][ipad]->SetLineWidth(2);
      hQtotPadsTruncPi[irow][ipad]->SetLineColor(kBlue);
      hQtotPadsTruncPi[irow][ipad]->Draw();
    }*/
    int entriesPi = hQtotPadsPi[irow][ipad]->GetEntries();
    hEntriesPi->Fill(entriesPi);
    int entriesEle = hQtotPadsEle[irow][ipad]->GetEntries();
    hEntriesEle->Fill(entriesEle);
    if (entriesPi > 40){
    GetBinMinMax(hQtotPadsPi[irow][ipad],frac,bin1,bin2);
        hQtotPadsPi[irow][ipad]->Fit("MPV","NQ0","", hQtotPadsPi[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPadsPi[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
        float mpv = MPVfit->GetParameter(1);
        float chisquare = MPVfit->GetChisquare()/MPVfit->GetNDF();
        hChisquarePi->Fill(chisquare);
        hGainMapPi->Fill(irow,ipad,mpv);
        ++Landau;


            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTruncPi[irow][ipad],&param,lowTr,upTr);
            hGainMapTruncPi->Fill(irow,ipad,param[1]);
            GainMapPi.setValue(irow,ipad,param[1]);

    }
    else if (entriesPi <= 40 && entriesPi >= 10){
      hGainMapPi->Fill(irow,ipad,hQtotPadsPi[irow][ipad]->GetMean());
        ++means;

            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTruncPi[irow][ipad],&param,lowTr,upTr);
            hGainMapTruncPi->Fill(irow,ipad,param[1]);
            GainMapPi.setValue(irow,ipad,param[1]);
    }


    if (entriesEle > 40){
    GetBinMinMax(hQtotPadsEle[irow][ipad],frac,bin1,bin2);
        hQtotPadsEle[irow][ipad]->Fit("MPV","NQ0","", hQtotPadsEle[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPadsEle[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
        float mpv = MPVfit->GetParameter(1);
        float chisquare = MPVfit->GetChisquare()/MPVfit->GetNDF();
        hChisquareEle->Fill(chisquare);
        hGainMapEle->Fill(irow,ipad,mpv);
        ++Landau;


            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTruncEle[irow][ipad],&param,lowTr,upTr);
            hGainMapTruncEle->Fill(irow,ipad,param[1]);
            GainMapEle.setValue(irow,ipad,param[1]);

    }
    else if (entriesEle <= 40 && entriesEle >= 10){
      hGainMapEle->Fill(irow,ipad,hQtotPadsEle[irow][ipad]->GetMean());
        ++means;

            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTruncEle[irow][ipad],&param,lowTr,upTr);
            hGainMapTruncEle->Fill(irow,ipad,param[1]);
            GainMapEle.setValue(irow,ipad,param[1]);
    }
         /* else if (entries > 0) {
            param[0]=0;
            param[1]=0;
            param[2]=0;
            param[3]=0;
            TruncatedMean(hQtotPadsTrunc[irow][ipad],&param,lowTr,upTr);
            hGainMapTrunc->Fill(irow,ipad,param[1]);
          }*/

    }
  }


//  hGainMapPi->GetZaxis()->SetRangeUser(.5,1.5);
  hGainMapTruncPi->GetZaxis()->SetRangeUser(.5,1.5);
//  hGainMapEle->GetZaxis()->SetRangeUser(.5,1.5);
  hGainMapTruncEle->GetZaxis()->SetRangeUser(.5,1.5);
//  fPi->cd();
//  hGainMapPi->Write("GainMap");
//  delete fPi;
//  fEle->cd();
//  hGainMapEle->Write("GainMap");
//  delete fEle;
//  chiPi->cd();
//  hChisquarePi->Write("ChiSquare");
//  delete chiPi;
//  chiEle->cd();
//  hChisquareEle->Write("ChiSquare");
//  delete chiEle;
//  entPi->cd();
//  hEntriesPi->Write("Entries");
//  delete entPi;
//  entEle->cd();
//  hEntriesEle->Write("Entries");
//  delete entEle;
//  gTrPi->cd();
//  hGainMapTruncPi->Write("GainMapPi");
//  delete gTrPi;
//  gTrEle->cd();
//  hGainMapTruncEle->Write("GainMapEle");
//  delete gTrEle;


  auto gmplotPi = Painter::getHistogram2D(GainMapPi);
  auto gmplotEle = Painter::getHistogram2D(GainMapEle);

  auto cgmplotPi = new TCanvas("cgmPi","GainMap Pions");
  gmplotPi->Draw("colz");

  auto cgmplotEle = new TCanvas("cgmEle","GainMap Electrons");
  gmplotEle->Draw("colz");

  CalDet<float> *gainmap = new CalDet<float>(PadSubset::ROC);
  gainmap->setName("GainMap");
  CalArray<float>& GainMap = gainmap->getCalArray(roc);

  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<100; ++ipad){
      float pival = gainmapPi->getValue(roc, irow, ipad);
      float eleval = gainmapEle->getValue(roc, irow, ipad);
      if (pival != 0 || eleval != 0) {
        float mean = (pival+eleval)/2;
        float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;
        GainMap.setValue(irow, ipad, mean);
        hGainMapHisto->Fill(irow, cpad, mean);
        if (mean > .7) {
          hGainDist->Fill(mean);
        }
      }
    }
  }

//  auto gmplot = Painter::getHistogram2D(GainMap);
//  auto cgmplot = new TCanvas("cgm","GainMap");
//  gmplot->Draw("colz");

  TFile *f = TFile::Open(Form("%s/GainMap.root", OutputPath), "recreate");
  f->WriteObject(gainmap, "GainMap");
  delete f;

//  gGM->cd();
//  hGainMapHisto->Write("GainMapHisto");
//  delete gGM;

  auto gaindist = new TCanvas("gaindist", "Gain Distribution");
  float hMean = hGainDist->GetMean();
  float hRMS = hGainDist->GetRMS();
  TPaveText *pave = new TPaveText(0.6,.7,.9,.9,"NDC");
  pave->SetBorderSize(1);
  pave->AddText(Form("mean: %.2f", hMean));
  pave->AddText(Form("RMS: %.2f", hRMS));
  hGainDist->Draw();
  pave->Draw("same");
  gaindist->Update();
  gaindist->Print(Form("%s/GainDistribution.png", OutputPath));


  TFile *g = TFile::Open(Form("%s/GainMap_Histos.root", OutputPath), "recreate");
  g->WriteObject(hGainMapHisto, "GainMap");
  g->WriteObject(hGainMapPi, "GainMapPi");
  g->WriteObject(hGainMapEle, "GainMapEle");
  g->WriteObject(hGainDist, "GainDistribution");
  g->WriteObject(hCheckExclude, "CheckExclude");
  delete g;

  cout<<endl<<endl<<endl<<"*** =============   done   ============= ***"<<endl<<endl<<endl;
  return 0;
}
