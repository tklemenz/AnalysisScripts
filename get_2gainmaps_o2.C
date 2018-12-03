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

float GetCPad(Cluster cl)
{
  Mapper &mapper = Mapper::instance();
  DigitPos pos(cl.getCRU(), PadPos(cl.getRow(), cl.getPadMean()));
  float row = pos.getPadSecPos().getPadPos().getRow();
  float pad = pos.getPadSecPos().getPadPos().getPad();
  return pad - mapper.getNumberOfPadsInRowSector(row)/2;
}

void TruncatedMean(TH1F *his, std::vector<float> *param, float down, float up, TH1F *hTrunc=nullptr)
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
      if (hTrunc != nullptr) {
        hTrunc->SetBinContent(ibin,his->GetBinContent(ibin));
      }
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
    (*param)[3] = nentries;
  }
}

bool InsideEdge(float row, float cpad, int FECSetting=0, int cut=2)
{
  if (FECSetting == 0) {

    if (row < 5) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 19) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 18) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 17) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 5) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 19) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 18) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 17) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 5 && row < 10) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 20) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 19) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 18) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 10) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 20) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 19) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 18) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 10 && row < 16) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 21) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 20) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 19) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 15) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 8 && cpad < 21) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 9 && cpad < 20) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 10 && cpad < 19) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 16) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 7 && cpad < 20) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 8 && cpad < 19) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 9 && cpad < 18) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 15 && row <= 22) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 8 && cpad < 22) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 9 && cpad < 21) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 10 && cpad < 20) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 22 && row <= 27) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 8 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 9 && cpad < 22) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 10 && cpad < 21) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 27 && row <= 30) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 8 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 9 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 10 && cpad < 22) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 31) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 9 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 10 && cpad < 22) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 11 && cpad < 21) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 32) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 22) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 21) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 32 && row <= 37) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 22) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 21) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 37 && row <= 43) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 22) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 44 || row == 45) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 25) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 23) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 46) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 11 && cpad < 25) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 12 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 13 && cpad < 23) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 47 || row == 48) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 11 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 12 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 13 && cpad < 22) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 49) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 23) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 22) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 50 || row == 51) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 10 && cpad < 25) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 11 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 12 && cpad < 23) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 52 && row <= 55) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 11 && cpad < 25) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 12 && cpad < 24) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 13 && cpad < 23) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row > 55 && row <= 60) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 11 && cpad < 26) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 12 && cpad < 25) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 13 && cpad < 24) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 61 || row == 62) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > 11 && cpad < 27) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > 12 && cpad < 26) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 13 && cpad < 25) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }
  }
  else if (FECSetting == 1) {
    if (row < 6) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -2 && cpad < 10) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -1 && cpad < 9) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 0 && cpad < 8) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 6 && row <=10) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -2 && cpad < 11) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -1 && cpad < 10) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 0 && cpad < 9) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 11) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -2 && cpad < 10) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -1 && cpad < 9) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > 0 && cpad < 8) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 12 && row <= 15) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -3 && cpad < 10) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -2 && cpad < 9) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -1 && cpad < 8) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 16 || row == 17) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -3 && cpad < 7) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -2 && cpad < 6) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -1 && cpad < 5) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 18 && row <= 29) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -8 && cpad < 7) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -7 && cpad < 6) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -6 && cpad < 5) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 30) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -8 && cpad < 6) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -7 && cpad < 5) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -6 && cpad < 4) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 31 || row == 32) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -8 && cpad < 4) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -7 && cpad < 3) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -6 && cpad < 2) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 33) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -9 && cpad < 4) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -8 && cpad < 3) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -7 && cpad < 2) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 34 && row <= 46) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -10 && cpad < 4) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -9 && cpad < 3) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -8 && cpad < 2) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 47) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -10 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -9 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -8 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 48) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -11 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -10 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -9 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 49) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -13 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -12 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -11 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 50) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -14 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -13 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -12 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 51 || row == 52) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -15 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -14 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -13 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row >= 53 && row <= 59) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -16 && cpad < -1) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -15 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -14 && cpad < -3) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 60) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -16 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -15 && cpad < -3) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -14 && cpad < -4) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

    else if (row == 61 || row == 62) {
      if (cut == 0) {return true;}
      else if (cut == 1) {
        if (cpad > -17 && cpad < -2) {return true;}
        else {return false;}
      }
      else if (cut == 2) {
        if (cpad > -16 && cpad < -3) {return true;}
        else {return false;}
      }
      else if (cut == 3) {
        if (cpad > -15 && cpad < -4) {return true;}
        else {return false;}
      }
      else {
        std::cout<<std::endl<<"funtion cannot cut more than three pads on each side. choose 0, 1, 2 or 3"<<std::endl;
        return true;
      }
    }

  }
  else {
    return false;
  }
}

bool removedRow(float row, float cpad, int FECSetting=0)
{
  if (FECSetting == 0) {
    if (row == 31 || row == 32 || row == 9 || row == 17 || row == 25 || row == 39 || row == 45 || row == 51) {
      return true;
    }
    else if (row > 56 || row < 2) {
      return true;
    }
    else {
      return false;
    }
  }

  else if (FECSetting == 1) {
    if (row == 31 || row == 32 || row == 9 || row == 17 || row == 25 || row == 39 || row == 45 || row == 51) {
      return true;
    }
    else if (row > 60 || row < 2) {
      return true;
    }
    else if (row > 1 && row < 61) {
      if (cpad == -1 || cpad == 0 || cpad == -2) {
        return true;
      }
    }
    else {
      return false;
    }
  }

  else {
    std::cout<<"no valid setting: choose 0 or 1."<<std::endl;
    return true;
  }
}

//float unbiasedTruncatedMean(TH2D *histo, int runNr, o2::TPC::TrackTPC Track, CalDet<int> *HitMapS0, CalDet<int> *HitMapS1, float low=.0, float high=.7)
//{
//  ROC roc(0);
//  std::vector<Cluster> Clusters;
//  Track.getClusterVector(Clusters);
//  std::vector<float> values;
//  Mapper &mapper = Mapper::instance();

//  for (auto &clusterObject : Clusters) {
//    DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
//    float row = pos.getPadSecPos().getPadPos().getRow();
//    float pad = pos.getPadSecPos().getPadPos().getPad();
//    float cpad = GetCPad(clusterObject);
//    if (row == 31 || row == 32) {
//      histo->Fill(row, cpad, 1);
//      continue;
//    }
//    if (runNr <= 243) {
//      if (row >= 57 || row < 2) {
//        histo->Fill(row, cpad, 1);
//        continue;
//      }

//      int nexthit = HitMapS0->getValue(roc, row, pad+1);
//      int previoushit = HitMapS0->getValue(roc, row, pad-1);
//      if (nexthit != 1 || previoushit != 1) {
//        histo->Fill(row, cpad, 1);
//        continue;
//      }
//    }
//    else if (runNr > 255) {
//      if (row == 62 || row < 2) {
//        histo->Fill(row, cpad, 1);
//        continue;
//      }
//      else if (row > 0 && row < 62) {
//        if (pad == mapper.getNumberOfPadsInRowSector(row)/2 - 1) {
//          histo->Fill(row, cpad, 1);
//          continue;
//        }
//      }
//      int nexthit = HitMapS1->getValue(roc, row, pad+1);
//      int previoushit = HitMapS1->getValue(roc, row, pad-1);
//      if (nexthit != 1 || previoushit != 1) {
//        histo->Fill(row, cpad, 1);
//        continue;
//      }
//    }
//    values.push_back(clusterObject.getQ());
//  }
//  transform(values.begin(), values.end(), values.begin(), boost::lambda::_1 * cos(atan(Track.getTgl())) * cos(asin(Track.getSnp())));
//  std::sort(values.begin(), values.end());

//  float dEdx = 0.f;
//  int nClustersTrunc = 0;
//  int nClustersUsed = static_cast<int>(values.size());

//  for (int icl=0; icl<nClustersUsed; ++icl) {
//    if (icl<std::round(low*nClustersUsed)) continue;
//    if (icl>std::round(high*nClustersUsed)) break;

//    dEdx+=values[icl];
//    ++nClustersTrunc;
//  }

//  if (nClustersTrunc>0){
//    dEdx/=nClustersTrunc;
//  }
//  return dEdx;
//}



//int get_gainmap_o2(int argc, char **argv){
  //int get_2gainmaps_o2(TString FileList, TString HitMapFile, const char *OutputPath="/home/tom") {
int get_2gainmaps_o2(TString FileList, TString GainMap="", const char *OutputPath="/home/tom/myana/Results") {
  //gROOT->ProcessLine(".x /lustre/nyx/alice/users/tklemenz/rootlogon.C");

  //const char *sarg = argv[1];
  //const char *spath = argv[2];
  /* ============================================================
   * ====================== Set up HitMap =======================
   * ============================================================ */

//    TFile h(HitMapFile);
//    gROOT->cd();

//    CalDet<int> *hitmapS0 = nullptr, *hitmapS1 = nullptr;
//    h.GetObject("HitMapS0", hitmapS0);
//    h.GetObject("HitMapS1", hitmapS1);

    bool usegainmap = true;
    if (GainMap == "") {
      usegainmap = false;
    }

    ROC roc(0);
    CalDet<float> *gainmapcorrS0 = nullptr, *gainmapcorrS1 = nullptr;

    if (usegainmap == true) {
      TFile f(GainMap);
      gROOT->cd();

      f.GetObject("GainMapS0", gainmapcorrS0);
      f.GetObject("GainMapS1", gainmapcorrS1);

//      CalArray<float>& test = gainmapcorrS0->getCalArray(roc);

//      auto gainmapplot = Painter::getHistogram2D(test);
//      auto cgainmapplot = new TCanvas("cgainmapplot", "input gainmap");
//      gainmapplot->Draw("colz");
    }

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
  int nclCut = 48;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMeanLow = 3;
  int timeMeanHigh =20;


  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/

  //TH1D *hChisquarePi            = new TH1D("hChisquarePi", "; Chisquare Landau; # counts", 400,0,200);
  //TH2D *hGainMapPi              = new TH2D("hGainMapPi", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncPiS0       = new TH2D("hGainMapTruncPiS0", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncPiS1       = new TH2D("hGainMapTruncPiS1", "; Row; Pad", 63,0,63,100,0,100);
  //TH1D *hEntriesPi              = new TH1D("EntriesPi", "; Entries; # counts",300,1,300000);
  //TH1D *hChisquareEle           = new TH1D("hChisquareEle", "; Chisquare Landau; # counts", 400,0,200);
  //TH2D *hGainMapEle             = new TH2D("hGainMapEle", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncEleS0      = new TH2D("hGainMapTruncEleS0", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapTruncEleS1      = new TH2D("hGainMapTruncEleS1", "; Row; Pad", 63,0,63,100,0,100);
  TH2D *hGainMapHistoS0         = new TH2D("GainMapS0", "; Row; Pad; Normalized gain", 63,0,63,100,-50,50);
  TH2D *hGainMapHistoS1         = new TH2D("GainMapS1", "; Row; Pad; Normalized gain", 63,0,63,100,-50,50);
  //TH1D *hEntriesEle             = new TH1D("EntriesEle", "; Entries; # counts",300,1,300000);
  TH1D *hGainDistS0             = new TH1D("NGDS1", "; Normalized gain; Counts",100,0,2);
  TH1D *hGainDistS1             = new TH1D("NGDS2", "; Normalized gain; Counts",100,0,2);
  TH2D *hCheckExclude           = new TH2D("CheckExclude", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hOverlapRatio           = new TH2D("OverlapRatio","; Row; Pad; Ratio",63,0,63,100,0,100);
  TH1D *hOverlapRatioDist       = new TH1D("OverlapRatioDist", "; Ratio in overlap region; Counts",100,0.5,2.5);
  TH2D *hStDevGainMapPiS0       = new TH2D("StDevGainMapPiS0", "; Row; Pad; StdDev", 63,0,63,100,-50,50);
  TH2D *hStDevGainMapEleS0      = new TH2D("StDevGainMapEleS0", "; Row; Pad; StdDev", 63,0,63,100,-50,50);
  TH2D *hStDevGainMapPiS1       = new TH2D("StDevGainMapPiS1", "; Row; Pad; StdDev", 63,0,63,100,-50,50);
  TH2D *hStDevGainMapEleS1      = new TH2D("StDevGainMapEleS1", "; Row; Pad; StdDev", 63,0,63,100,-50,50);
  TH1D *hStDevGainMapPiS0Dist   = new TH1D("StDevGainMapPiS0Dist", "; StdDev normalized gain; Counts", 1000,0,10);
  TH1D *hStDevGainMapEleS0Dist  = new TH1D("StDevGainMapEleS0Dist", "; StdDev normalized gain; Counts", 1000,0,10);
  TH1D *hStDevGainMapPiS1Dist   = new TH1D("StDevGainMapPiS1Dist", "; StdDev normalized gain; Counts", 1000,0,10);
  TH1D *hStDevGainMapEleS1Dist  = new TH1D("StDevGainMapEleS1Dist", "; StdDev normalized gain; Counts", 1000,0,10);
  TH1D *hnClTrMean              = new TH1D("NClTrMean",";NClTrMean",70,0,70);




  /* initialize array of histograms used to get the gainmap (for pions and electrons seperately)*/
  const int nrows = 63;
  const int npads = 100;
  TH1F *hQtotPadsPi[nrows][npads];
  TH1F *hQtotPadsTruncPiS0[nrows][npads];
  TH1F *hQtotPadsTruncPiS1[nrows][npads];
  TH1F *hQtotTruncMeanPiS0[nrows][npads];
  TH1F *hQtotTruncMeanPiS1[nrows][npads];

  for(int irow=0; irow<nrows; ++irow) {
    for (int ipad=0; ipad<npads; ++ipad) {
    hQtotPadsPi[irow][ipad] = new TH1F(Form("hQtotPadPi_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotPadsTruncPiS0[irow][ipad] = new TH1F(Form("hQtotPadPi_truncS0_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotPadsTruncPiS1[irow][ipad] = new TH1F(Form("hQtotPadPi_truncS1_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotTruncMeanPiS0[irow][ipad] = new TH1F(Form("hQtotTruncMeanPi_S0_%i_%i", irow, ipad), "; Truncated mean",500,0,10);
    hQtotTruncMeanPiS1[irow][ipad] = new TH1F(Form("hQtotTruncMeanPi_S1_%i_%i", irow, ipad), "; Truncated mean",500,0,10);
    }
  }

  TH1F *hQtotPadsEle[nrows][npads];
  TH1F *hQtotPadsTruncEleS0[nrows][npads];
  TH1F *hQtotPadsTruncEleS1[nrows][npads];
  TH1F *hQtotTruncMeanEleS0[nrows][npads];
  TH1F *hQtotTruncMeanEleS1[nrows][npads];


  for(int irow=0; irow<nrows; ++irow) {
    for (int ipad=0; ipad<npads; ++ipad) {
    hQtotPadsEle[irow][ipad] = new TH1F(Form("hQtotPadEle_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotPadsTruncEleS0[irow][ipad] = new TH1F(Form("hQtotPadEle_truncS0_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotPadsTruncEleS1[irow][ipad] = new TH1F(Form("hQtotPadEle_truncS1_%i_%i", irow, ipad), "; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts", 500, 0, 50);
    hQtotTruncMeanEleS0[irow][ipad] = new TH1F(Form("hQtotTruncMeanEle_S0_%i_%i", irow, ipad), "; Truncated mean",500,0,10);
    hQtotTruncMeanEleS1[irow][ipad] = new TH1F(Form("hQtotTruncMeanEle_S1_%i_%i", irow, ipad), "; Truncated mean",500,0,10);
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
  int setting = 0;
  int loopcounter = 0;
  int loopcounter2 = 0;
  bool setgainmapS1 = false;
  int nClTrMean = 0;

  std::vector<TrackTPC> *vecEvent = 0;


  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;

  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
  TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
  cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");

    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);
    //cout<<"KOMME BIS HIER!"<<endl;

    CherenkovValue = 0;
    runNr = 0;


    for (int iev=0; iev<tree->GetEntriesFast(); ++iev) {
      tree->GetEntry(iev);

      int NTracks = vecEvent->size();
/*========================================== CUT ================================================*/
      if (NTracks != TrPerEv) continue;																// only one-track events

      CherenkovValue = Header.cherenkovValue;
      runNr = Header.run;

      Tracks += NTracks;

      ++OneTrackEvents;
      //cout<<"DEBUG POINT 4"<<endl;

    for (auto trackObject : *vecEvent) {
      std::vector<Cluster> clCont;
      trackObject.getClusterVector(clCont);

      setting = 0;
      if (runNr >= 255) {setting = 1;}
      if (setting == 1) {
        setgainmapS1 = true;
        ++loopcounter2;
      }
      trackObject.getClusterVector(clCont);
      if ((loopcounter == 0) && usegainmap) {
        trackObject.setGainMap(GainMap, setting);
        std::cout<<std::endl<<"set gainmap for setting 0"<<std::endl;
      }
      else if (loopcounter2 == 1 && usegainmap) {
        trackObject.setGainMap(GainMap, setting);
        std::cout<<std::endl<<"set gainmap for setting 1"<<std::endl;
      }
      ++loopcounter;

      int ncl = clCont.size();
      isok = true;

      onetrcl += ncl;


      //cout<<"DEBUG POINT 5"<<endl;


  /*========================================== CUT ================================================*/
      if (ncl < nclCut) continue;																	// cut on number of clusters per track

  /*========================================== CUT ================================================*/
      if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov


      for (auto &clusterObject : clCont) {															// make cuts on cluster properties
        float timeMean = clusterObject.getTimeMean();
  /*========================================== CUT ================================================*/
        if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle
      }

      if (isok == true){																			// track accepted
        usedcl += ncl;
        float dEdx = trackObject.getTruncatedMean(runNr, .0, .7, 1, false, true, true, hCheckExclude, 1, 2, &nClTrMean);
        hnClTrMean->Fill(nClTrMean);
        //float dEdx = unbiasedTruncatedMean(hCheckExclude, runNr, trackObject, hitmapS0, hitmapS1);
        if (CherenkovValue < CherCutLow){
              for (auto &clusterObject : clCont) {															// loop over clusters
                DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
                int row = pos.getPadSecPos().getPadPos().getRow();
                int pad = pos.getPadSecPos().getPadPos().getPad();
                float cpad = GetCPad(clusterObject);
                float QTot = clusterObject.getQ();
                if (usegainmap == true) {
                  if (runNr <= 243) {
                    float corr = gainmapcorrS0->getValue(roc, row, pad);
                    QTot /= corr;
                  }
                  else if (runNr >255) {
                    float corr = gainmapcorrS1->getValue(roc, row, pad);
                    QTot /= corr;
                  }
                }

                if (!InsideEdge(row, cpad, setting, 2)) continue;
                //if (removedRow(row, cpad, setting)) continue;

                if (runNr <= 243) {
                  hQtotPadsPi[row][pad]->Fill(QTot/dEdx);
                  hQtotPadsTruncPiS0[row][pad]->Fill(QTot/dEdx);
                }
                else if (runNr > 255) {
                  hQtotPadsPi[row][pad]->Fill(QTot/dEdx);
                  hQtotPadsTruncPiS1[row][pad]->Fill(QTot/dEdx);
                }
              }
              ++pions;

        }
        if (CherenkovValue > CherCutHigh){
              for (auto &clusterObject : clCont) {															// loop over clusters
                DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
                int row = pos.getPadSecPos().getPadPos().getRow();
                int pad = pos.getPadSecPos().getPadPos().getPad();
                float cpad = GetCPad(clusterObject);
                float QTot = clusterObject.getQ();
                if (usegainmap == true) {
                  if (runNr <= 243) {
                    float corr = gainmapcorrS0->getValue(roc, row, pad);
                    QTot /= corr;
                  }
                  else if (runNr >255) {
                    float corr = gainmapcorrS1->getValue(roc, row, pad);
                    QTot /= corr;
                  }
                }

                if (!InsideEdge(row, cpad, setting, 2)) continue;
                //if (removedRow(row, cpad, setting)) continue;

                if (runNr <= 243) {
                  hQtotPadsEle[row][pad]->Fill(QTot/dEdx);
                  hQtotPadsTruncEleS0[row][pad]->Fill(QTot/dEdx);
                }
                else if (runNr > 255) {
                  hQtotPadsEle[row][pad]->Fill(QTot/dEdx);
                  hQtotPadsTruncEleS1[row][pad]->Fill(QTot/dEdx);
                }
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

  CalDet<float> *gainmapEleS0 = new CalDet<float>(PadSubset::ROC);
  gainmapEleS0->setName("GainMapEle");
  CalArray<float>& GainMapEleS0 = gainmapEleS0->getCalArray(roc);

  CalDet<float> *gainmapPiS0 = new CalDet<float>(PadSubset::ROC);
  gainmapPiS0->setName("GainMapPi");
  CalArray<float>& GainMapPiS0 = gainmapPiS0->getCalArray(roc);

  CalDet<float> *gainmapEleS1 = new CalDet<float>(PadSubset::ROC);
  gainmapEleS1->setName("GainMapEle");
  CalArray<float>& GainMapEleS1 = gainmapEleS1->getCalArray(roc);

  CalDet<float> *gainmapPiS1 = new CalDet<float>(PadSubset::ROC);
  gainmapPiS1->setName("GainMapPi");
  CalArray<float>& GainMapPiS1 = gainmapPiS1->getCalArray(roc);

//  TF1 *MPVfit = new TF1("MPV","landau");
//  int means = 0;
//  int Landau = 0;
//  const Float_t frac=0.2;
//  Int_t bin1=0,bin2=0;
  //TVectorD param(4);
  std::vector<float> param(4);
  float lowTr = 0.;
  float upTr = .7;


  for (int irow=0; irow<nrows; ++irow){
    for (int ipad=0; ipad<npads; ++ipad){
      float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;
      //std::cout<<"row: "<<irow<<"\t"<<"pad: "<<ipad<<"\t"<<"cpad: "<<cpad<<std::endl;

//    if (irow == 62 && ipad == 66){
//      hQtotPadsTruncPi[irow][ipad]->SetLineWidth(2);
//      hQtotPadsTruncPi[irow][ipad]->SetLineColor(kBlue);
//      hQtotPadsTruncPi[irow][ipad]->Draw();
//    }
      //int entriesPi = hQtotPadsPi[irow][ipad]->GetEntries();
      int entriesPiS0 = hQtotPadsTruncPiS0[irow][ipad]->GetEntries();
      //hEntriesPi->Fill(entriesPi);
      int entriesEleS0 = hQtotPadsTruncEleS0[irow][ipad]->GetEntries();
      //hEntriesEleS0->Fill(entriesEle);
      int entriesPiS1 = hQtotPadsTruncPiS1[irow][ipad]->GetEntries();
      //hEntriesPi->Fill(entriesPi);
      int entriesEleS1 = hQtotPadsTruncEleS1[irow][ipad]->GetEntries();
      //hEntriesEleS1->Fill(entriesEle);
      if (entriesPiS0 > 0){
//        GetBinMinMax(hQtotPadsPi[irow][ipad],frac,bin1,bin2);
//        hQtotPadsPi[irow][ipad]->Fit("MPV","NQ0","", hQtotPadsPi[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPadsPi[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
//        float mpv = MPVfit->GetParameter(1);
//        float chisquare = MPVfit->GetChisquare()/MPVfit->GetNDF();
//        hChisquarePi->Fill(chisquare);
//        hGainMapPi->Fill(irow,ipad,mpv);
//        ++Landau;


        param[0]=0;
        param[1]=0;
        param[2]=0;
        param[3]=0;
        TruncatedMean(hQtotPadsTruncPiS0[irow][ipad],&param,lowTr,upTr);
        hGainMapTruncPiS0->Fill(irow,ipad,param[1]);
        GainMapPiS0.setValue(irow,ipad,param[1]);
        hQtotTruncMeanPiS0[irow][ipad]->Fill(param[1]);
        if (!(TMath::IsNaN(param[2]))) {
          //std::cout<<"row: "<<irow<<"\t"<<"cpad: "<<cpad<<"\t"<<"stdev: "<<param[2]<<"\t"<<"entries: "<<param[3]<<std::endl;
          hStDevGainMapPiS0->Fill(irow,cpad,TMath::Sqrt(param[2]));
          hStDevGainMapPiS0Dist->Fill(TMath::Sqrt(param[2]));
        }
      }
      if (entriesPiS1 > 0) {
        param[0]=0;
        param[1]=0;
        param[2]=0;
        param[3]=0;
        TruncatedMean(hQtotPadsTruncPiS1[irow][ipad],&param,lowTr,upTr);
        hGainMapTruncPiS1->Fill(irow,ipad,param[1]);
        GainMapPiS1.setValue(irow,ipad,param[1]);
        hQtotTruncMeanPiS1[irow][ipad]->Fill(param[1]);
        if (!(TMath::IsNaN(param[2]))) {
          hStDevGainMapPiS1->Fill(irow,cpad,TMath::Sqrt(param[2]));
          hStDevGainMapPiS1Dist->Fill(TMath::Sqrt(param[2]));
        }

      }
//      else if (entriesPi <= 40 && entriesPi > 0){
//        hGainMapPi->Fill(irow,ipad,hQtotPadsPi[irow][ipad]->GetMean());
//        ++means;

//        param[0]=0;
//        param[1]=0;
//        param[2]=0;
//        param[3]=0;
//        TruncatedMean(hQtotPadsTruncPiS0[irow][ipad],&param,lowTr,upTr);
//        hGainMapTruncPiS0->Fill(irow,ipad,param[1]);
//        GainMapPiS0.setValue(irow,ipad,param[1]);
//        if (!(TMath::IsNaN(param[2]))) {
//          hStDevGainMapPiS0->Fill(irow,cpad,param[2]);
//          hStDevGainMapPiS0Dist->Fill(param[2]);
//        }

//        param[0]=0;
//        param[1]=0;
//        param[2]=0;
//        param[3]=0;
//        TruncatedMean(hQtotPadsTruncPiS1[irow][ipad],&param,lowTr,upTr);
//        hGainMapTruncPiS1->Fill(irow,ipad,param[1]);
//        GainMapPiS1.setValue(irow,ipad,param[1]);
//        if (!(TMath::IsNaN(param[2]))) {
//          hStDevGainMapPiS1->Fill(irow,cpad,param[2]);
//          hStDevGainMapPiS1Dist->Fill(param[2]);
//        }
//      }


      if (entriesEleS0 > 0){
//        GetBinMinMax(hQtotPadsEle[irow][ipad],frac,bin1,bin2);
//        hQtotPadsEle[irow][ipad]->Fit("MPV","NQ0","", hQtotPadsEle[irow][ipad]->GetXaxis()->GetBinLowEdge(bin1),hQtotPadsEle[irow][ipad]->GetXaxis()->GetBinUpEdge(bin2));
//        float mpv = MPVfit->GetParameter(1);
//        float chisquare = MPVfit->GetChisquare()/MPVfit->GetNDF();
//        hChisquareEle->Fill(chisquare);
//        hGainMapEle->Fill(irow,ipad,mpv);
//        ++Landau;


        param[0]=0;
        param[1]=0;
        param[2]=0;
        param[3]=0;
        TruncatedMean(hQtotPadsTruncEleS0[irow][ipad],&param,lowTr,upTr);
        hGainMapTruncEleS0->Fill(irow,ipad,param[1]);
        GainMapEleS0.setValue(irow,ipad,param[1]);
        hQtotTruncMeanEleS0[irow][ipad]->Fill(param[1]);
        if (!(TMath::IsNaN(param[2]))) {
          hStDevGainMapEleS0->Fill(irow,cpad,TMath::Sqrt(param[2]));
          hStDevGainMapEleS0Dist->Fill(TMath::Sqrt(param[2]));
        }
      }
      if (entriesEleS1 > 0) {
        param[0]=0;
        param[1]=0;
        param[2]=0;
        param[3]=0;
        TruncatedMean(hQtotPadsTruncEleS1[irow][ipad],&param,lowTr,upTr);
        hGainMapTruncEleS1->Fill(irow,ipad,param[1]);
        GainMapEleS1.setValue(irow,ipad,param[1]);
        hQtotTruncMeanEleS1[irow][ipad]->Fill(param[1]);
        if (!(TMath::IsNaN(param[2]))) {
          hStDevGainMapEleS1->Fill(irow,cpad,TMath::Sqrt(param[2]));
          hStDevGainMapEleS1Dist->Fill(TMath::Sqrt(param[2]));
        }

    }
//      else if (entriesEle <= 40 && entriesEle > 0){
//        hGainMapEle->Fill(irow,ipad,hQtotPadsEle[irow][ipad]->GetMean());
//        ++means;

//        param[0]=0;
//        param[1]=0;
//        param[2]=0;
//        param[3]=0;
//        TruncatedMean(hQtotPadsTruncEleS0[irow][ipad],&param,lowTr,upTr);
//        hGainMapTruncEleS0->Fill(irow,ipad,param[1]);
//        GainMapEleS0.setValue(irow,ipad,param[1]);
//        if (!(TMath::IsNaN(param[2]))) {
//          hStDevGainMapEleS0->Fill(irow,cpad,param[2]);
//          hStDevGainMapEleS0Dist->Fill(param[2]);
//        }

//        param[0]=0;
//        param[1]=0;
//        param[2]=0;
//        param[3]=0;
//        TruncatedMean(hQtotPadsTruncEleS1[irow][ipad],&param,lowTr,upTr);
//        hGainMapTruncEleS1->Fill(irow,ipad,param[1]);
//        GainMapEleS1.setValue(irow,ipad,param[1]);
//        if (!(TMath::IsNaN(param[2]))) {
//          hStDevGainMapEleS1->Fill(irow,cpad,param[2]);
//          hStDevGainMapEleS1Dist->Fill(param[2]);
//        }

//    }
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
//  hGainMapTruncPiS0->GetZaxis()->SetRangeUser(.5,1.5);
//  hGainMapTruncPiS1->GetZaxis()->SetRangeUser(.5,1.5);

//  hGainMapEle->GetZaxis()->SetRangeUser(.5,1.5);
//  hGainMapTruncEleS0->GetZaxis()->SetRangeUser(.5,1.5);
//  hGainMapTruncEleS1->GetZaxis()->SetRangeUser(.5,1.5);

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


  auto gmplotPi = Painter::getHistogram2D(GainMapPiS0);
  auto gmplotEle = Painter::getHistogram2D(GainMapEleS0);

  auto cgmplotPi = new TCanvas("cgmPi","GainMap Pions");
  gmplotPi->Draw("colz");

  auto cgmplotEle = new TCanvas("cgmEle","GainMap Electrons");
  gmplotEle->Draw("colz");

  CalDet<float> *gainmapS0 = new CalDet<float>(PadSubset::ROC);
  gainmapS0->setName("GainMap");
  CalArray<float>& GainMapS0 = gainmapS0->getCalArray(roc);

  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<100; ++ipad){
      if (TMath::IsNaN(gainmapPiS0->getValue(roc, irow, ipad))) {
        GainMapPiS0.setValue(irow, ipad, 0);
      }
      if (TMath::IsNaN(gainmapPiS1->getValue(roc, irow, ipad))) {
        GainMapPiS1.setValue(irow, ipad, 0);
      }
      if (TMath::IsNaN(gainmapEleS0->getValue(roc, irow, ipad))) {
        GainMapEleS0.setValue(irow, ipad, 0);
      }
      if (TMath::IsNaN(gainmapEleS1->getValue(roc, irow, ipad))) {
        GainMapEleS1.setValue(irow, ipad, 0);
      }
    }
  }

/*================================ fill Gain map S0 ======================================*/
  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<100; ++ipad){
      float pival = gainmapPiS0->getValue(roc, irow, ipad);
      float eleval = gainmapEleS0->getValue(roc, irow, ipad);
      if (pival != 0 || eleval != 0) {
        float mean = (pival+eleval)/2;
        float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;
        GainMapS0.setValue(irow, ipad, mean);
        hGainMapHistoS0->Fill(irow, cpad, mean);
        //if (mean > .7) {
          hGainDistS0->Fill(mean);
        //}
      }
    }
  }

  /*================================ fill Gain map S1 ======================================*/
  CalDet<float> *gainmapS1 = new CalDet<float>(PadSubset::ROC);
  gainmapS1->setName("GainMap");
  CalArray<float>& GainMapS1 = gainmapS1->getCalArray(roc);

  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<100; ++ipad){
      float pival = gainmapPiS1->getValue(roc, irow, ipad);
      float eleval = gainmapEleS1->getValue(roc, irow, ipad);
      if (pival != 0 || eleval != 0) {
        float mean = (pival+eleval)/2;
        float cpad = ipad - mapper.getNumberOfPadsInRowSector(irow)/2;
        GainMapS1.setValue(irow, ipad, mean);
        hGainMapHistoS1->Fill(irow, cpad, mean);
        //if (mean > .7) {
          hGainDistS1->Fill(mean);
        //}
      }
    }
  }


  for (int irow = 0; irow<63; ++irow) {
    for (int ipad = 0; ipad<100; ++ipad){
      float S0val = gainmapS0->getValue(roc, irow, ipad);
      float S1val = gainmapS1->getValue(roc, irow, ipad);
      if (S0val != 0 && S1val != 0) {
        hOverlapRatio->Fill(irow, ipad, S0val/S1val);
        hOverlapRatioDist->Fill(S0val/S1val);
      }
    }
  }

//  auto gmplot = Painter::getHistogram2D(GainMap);
//  auto cgmplot = new TCanvas("cgm","GainMap");
//  gmplot->Draw("colz");

  TFile *f = TFile::Open(Form("%s/GainMap2.root", OutputPath), "recreate");
  f->WriteObject(gainmapS0, "GainMapS0");
  f->WriteObject(gainmapS1, "GainMapS1");
  delete f;

//  gGM->cd();
//  hGainMapHisto->Write("GainMapHisto");
//  delete gGM;

  auto gaindistS0 = new TCanvas("gaindistS0", "Gain Distribution S0");
  float hMean0 = hGainDistS0->GetMean();
  float hRMS0 = hGainDistS0->GetRMS();
  TPaveText *pave0 = new TPaveText(0.6,.7,.9,.9,"NDC");
  pave0->SetBorderSize(1);
  pave0->AddText(Form("mean: %.2f", hMean0));
  pave0->AddText(Form("RMS: %.2f", hRMS0));
  hGainDistS0->Draw();
  pave0->Draw("same");
  gaindistS0->Update();
  gaindistS0->Print(Form("%s/GainDistributionS0.png", OutputPath));

  auto gaindistS1 = new TCanvas("gaindistS1", "Gain Distribution S1");
  float hMean1 = hGainDistS1->GetMean();
  float hRMS1 = hGainDistS1->GetRMS();
  TPaveText *pave1 = new TPaveText(0.6,.7,.9,.9,"NDC");
  pave1->SetBorderSize(1);
  pave1->AddText(Form("mean: %.2f", hMean1));
  pave1->AddText(Form("RMS: %.2f", hRMS1));
  hGainDistS1->Draw();
  pave1->Draw("same");
  gaindistS1->Update();
  gaindistS1->Print(Form("%s/GainDistributionS1.png", OutputPath));

  TH1F *hTruncDist = new TH1F("hTruncDist","; Q_{tot}/<dE/dx>_{tr} [a.u.]; Counts",500,0,50);
  TruncatedMean(hQtotPadsTruncPiS1[14][38],&param,lowTr,upTr,hTruncDist);



  TFile *g = TFile::Open(Form("%s/GainMap_Histos2.root", OutputPath), "recreate");
  g->WriteObject(hGainMapHistoS0, "GainMapS0");
  g->WriteObject(hGainMapHistoS1, "GainMapS1");
  g->WriteObject(hGainMapTruncPiS0, "GainMapPiS0");
  g->WriteObject(hGainMapTruncPiS1, "GainMapPiS1");
  g->WriteObject(hGainMapTruncEleS0, "GainMapEleS0");
  g->WriteObject(hGainMapTruncEleS1, "GainMapEleS1");
  g->WriteObject(hGainDistS0, "GainDistributionS0");
  g->WriteObject(hGainDistS1, "GainDistributionS1");
  g->WriteObject(hCheckExclude, "CheckExclude");
  g->WriteObject(hOverlapRatio, "OverlapRatio");
  g->WriteObject(hOverlapRatioDist, "OverlapRatioDist");
  g->WriteObject(hStDevGainMapPiS0, "StDevGainMapPiS0");
  g->WriteObject(hStDevGainMapEleS0, "StDevGainMapEleS0");
  g->WriteObject(hStDevGainMapPiS1, "StDevGainMapPiS1");
  g->WriteObject(hStDevGainMapEleS1, "StDevGainMapEleS1");
  g->WriteObject(hStDevGainMapPiS0Dist, "StDevGainMapPiS0Dist");
  g->WriteObject(hStDevGainMapEleS0Dist, "StDevGainMapEleS0Dist");
  g->WriteObject(hStDevGainMapPiS1Dist, "StDevGainMapPiS1Dist");
  g->WriteObject(hStDevGainMapEleS1Dist, "StDevGainMapEleS1Dist");
  g->WriteObject(hnClTrMean, "nClTrMean");
  g->WriteObject(hTruncDist, "hTruncDist_14_38");
//  for (int irow=0; irow<nrows; ++irow){
//    for (int ipad=0; ipad<npads; ++ipad){
//      if (hQtotPadsTruncPiS0[irow][ipad]->GetEntries() > 0) {
//        g->WriteObject(hQtotPadsTruncPiS0[irow][ipad], Form("hQtotPadPi_truncS0_%i_%i", irow, ipad));
//        g->WriteObject(hQtotTruncMeanPiS0[irow][ipad], Form("hTruncMeanPiS0_%i_%i", irow, ipad));
//      }
//      if (hQtotPadsTruncEleS0[irow][ipad]->GetEntries() > 0) {
//        g->WriteObject(hQtotPadsTruncEleS0[irow][ipad], Form("hQtotPadEle_truncS0_%i_%i", irow, ipad));
//        g->WriteObject(hQtotTruncMeanEleS0[irow][ipad], Form("hTruncMeanEleS0_%i_%i", irow, ipad));
//      }
//      if (hQtotPadsTruncPiS1[irow][ipad]->GetEntries() > 0) {
//        g->WriteObject(hQtotPadsTruncPiS1[irow][ipad], Form("hQtotPadPi_truncS1_%i_%i", irow, ipad));
//        g->WriteObject(hQtotTruncMeanPiS1[irow][ipad], Form("hTruncMeanPiS1_%i_%i", irow, ipad));
//      }
//      if (hQtotPadsTruncEleS1[irow][ipad]->GetEntries() > 0) {
//        g->WriteObject(hQtotPadsTruncEleS1[irow][ipad], Form("hQtotPadEle_truncS1_%i_%i", irow, ipad));
//        g->WriteObject(hQtotTruncMeanEleS1[irow][ipad], Form("hTruncMeanEleS1_%i_%i", irow, ipad));
//      }
//    }
//  }
  delete g;

  cout<<endl<<endl<<endl<<"*** =============   done   ============= ***"<<endl<<endl<<endl;
  return 0;
}
