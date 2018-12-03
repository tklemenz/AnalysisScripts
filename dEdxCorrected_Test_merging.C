#include "TPCSimulation/Cluster.h"
#include "TPCBase/Mapper.h"
#include "TPCReconstruction/TrackTPC.h"
#include "TPCBase/CRU.h"
#include "TPCBase/ROC.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/Painter.h"


#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TSystem.h"
#include "TString.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TChain.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <math.h>
#include <boost/lambda/lambda.hpp>



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


void GetBinMinMax(const TH1 *hist, const float fracLow, const float fracHigh, int &bin1, int &bin2)
{
  const int binMax=hist->GetMaximumBin();
  const float contMax=hist->GetBinContent(binMax);
  bin1=binMax;
  bin2=binMax;
  while ( (bin1--)>binMax/3. ) if (hist->GetBinContent(bin1)<fracLow*contMax) break;
  while ( (bin2++)<binMax*3. ) if (hist->GetBinContent(bin2)<fracHigh*contMax) break;
}

float GetCPad(Cluster cl)
{
  Mapper &mapper = Mapper::instance();
  DigitPos pos(cl.getCRU(), PadPos(cl.getRow(), cl.getPadMean()));
  float row = pos.getPadSecPos().getPadPos().getRow();
  float pad = pos.getPadSecPos().getPadPos().getPad();
  return pad - mapper.getNumberOfPadsInRowSector(row)/2;
}

bool InsideEdge(float row, float cpad, int FECSetting=0, int cut=1)
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

void dEdxCorrected_Test_merging(TString FileList, TString GainMap="", const char *OutputPath="/home/tom/myana/Results")
{
  //gROOT->ProcessLine(".x /lustre/nyx/alice/users/tklemenz/rootlogon.C");
  //gStyle->SetOptStat(0);

  EventHeader Header;

  bool usegainmap = true;
  if (GainMap == "") {
    usegainmap = false;
  }

  CalDet<float> *gainmapS0 = nullptr, *gainmapS1 = nullptr;
  ROC roc(0);
  if (GainMap != "") {
    TFile f(GainMap);
    gROOT->cd();

    f.GetObject("GainMapS0", gainmapS0);
    f.GetObject("GainMapS1", gainmapS1);
  }

  float electronmeanTot, electronsigmaTot, pionmeanTot, pionsigmaTot, electronmeanerr, electronsigmaerr, pionmeanerr, pionsigmaerr, electronmeanMax,
      electronsigmaMax, pionmeanMax, pionsigmaMax, electronmeanerrMax, electronsigmaerrMax, pionmeanerrMax, pionsigmaerrMax, pionres, pionresMax, electronres, electronresMax, separationpower, separationpowerMax, nclRun;
  int pions, electrons, pioncl, elecl, nclTrack, nclEvent;

  electronmeanTot = electronsigmaTot = pionmeanTot = pionsigmaTot = electronmeanerr = electronsigmaerr = pionmeanerr = pionsigmaerr = electronmeanMax =
      electronsigmaMax = pionmeanMax = pionsigmaMax = electronmeanerrMax = electronsigmaerrMax = pionmeanerrMax = pionsigmaerrMax = pionres =
      pionresMax = electronres = electronresMax = separationpower = separationpowerMax = nclRun = 0.;
  pions = electrons = pioncl = elecl = nclTrack = nclEvent = 0;

  TFile *g = TFile::Open(Form("%s/Check_Merging.root", OutputPath), "recreate");
  TTree *TrackAna = new TTree("TrackAna","dEdxAna");

  TrackAna->Branch("PidEdxResQ", &pionres, "pionres/F");
  TrackAna->Branch("PidEdxResQmax", &pionresMax, "pionresMax/F");
  TrackAna->Branch("EledEdxResQ", &electronres, "electronres/F");
  TrackAna->Branch("EledEdxResQmax", &electronresMax, "electronresMax/F");
  TrackAna->Branch("SepPowerQ", &separationpower, "separationpower/F");
  TrackAna->Branch("SepPowerQmax", &separationpowerMax);
  TrackAna->Branch("posEle", &electronmeanTot, "electronmeanTot/F");
  TrackAna->Branch("resEle", &electronsigmaTot, "electronsigmaTot/F");
  TrackAna->Branch("posPi", &pionmeanTot, "pionmeanTot/F");
  TrackAna->Branch("resPi", &pionsigmaTot, "pionsigmaTot/F");
  TrackAna->Branch("posEleErr", &electronmeanerr, "electronmeanerr/F");
  TrackAna->Branch("resEleErr", &electronsigmaerr, "electronsigmaerr/F");
  TrackAna->Branch("posPiErr", &pionmeanerr, "pionmeanerr/F");
  TrackAna->Branch("resPiErr", &pionsigmaerr, "pionsigmaerr/F");
  TrackAna->Branch("posEleQmax", &electronmeanMax, "electronmeanMax/F");
  TrackAna->Branch("resEleQmax", &electronsigmaMax, "electronsigmaMax/F");
  TrackAna->Branch("posPiQmax", &pionmeanMax, "pionmeanMax/F");
  TrackAna->Branch("resPiQmax", &pionsigmaMax, "pionsigmaMax/F");
  TrackAna->Branch("posEleQmaxErr", &electronmeanerrMax, "electronmeanerrMax/F");
  TrackAna->Branch("resEleQmaxErr", &electronsigmaerrMax, "electronsigmaerrMax/F");
  TrackAna->Branch("posPiQmaxErr", &pionmeanerrMax, "pionmeanerrMax/F");
  TrackAna->Branch("resPiQmaxErr", &pionsigmaerrMax, "pionsigmaerrMax/F");
  TrackAna->Branch("nPions", &pions, "Pions/I");
  TrackAna->Branch("nElectrons", &electrons, "Electrons/I");
  TrackAna->Branch("nPionCl", &pioncl, "PionClusters/I");
  TrackAna->Branch("nEleCl", &elecl, "ElectronClusters/I");

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

  int setting = 0;
  int nclCut = 47;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMeanLow = 6;
  int timeMeanHigh = 20;
  int nEdgeCut = 2;
  float fracEdgeCl = .3;
  int excludeEdge = 1;

  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/

  /*
  const int nRuns = 350;
  TH1F *hQmaxPi[nRuns];

  const int nHVSettings = 11;
  TH1F *hdEdxEleTotAll[nHVSettings];
  TH1F *hdEdxPionTotAll[nHVSettings];

  TH2D *hQmaxPiRuns         = new TH2D("hQmaxPiRuns","; Run; Q_{max}Pi",350,0,350,300,0,300);
  TH2D *hMPVRuns            = new TH2D("MPVRuns", "; Run; MPV_Qmax",350,0,350,1000,0,200);
  TH2D *hChisquareMPVfit    = new TH2D("hChisqaureMPVfit","; Run; Chisquare",350,0,350,150,0,15);
  TH2D *hExcludeHisto       = new TH2D ("hExcludeHisto","; Pad row; Pad", 63,0,63,100,-50,50);

  for (int i = 0; i < nHVSettings; ++i) {
    hdEdxEleTotAll[i]      = new TH1F(Form("hdEdxEleTotAll_HVSetting_%i", i-1), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 15);
    hdEdxPionTotAll[i]     = new TH1F(Form("hdEdxTotAll_HVSetting_%i", i-1), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 15);
  }

  TGraph *gSepVsHVSet = new TGraph();
  gSepVsHVSet->SetMarkerStyle(20);
  gSepVsHVSet->GetXaxis()->SetTitle("HV Setting");
  gSepVsHVSet->GetYaxis()->SetTitle("Separation power");

  TGraph *gSepVsFeSigma = new TGraph();
  gSepVsFeSigma->SetMarkerStyle(20);
  gSepVsFeSigma->GetXaxis()->SetTitle("#sigma Fe55 [%]");
  gSepVsFeSigma->GetYaxis()->SetTitle("Separation power");

  TGraph *gSepVsIBF = new TGraph();
  gSepVsIBF->SetMarkerStyle(20);
  gSepVsIBF->GetXaxis()->SetTitle("#it{IB} [%]");
  gSepVsIBF->GetYaxis()->SetTitle("Separation power");

  TGraph *gdEdxPiVsFeSigma = new TGraph();
  gdEdxPiVsFeSigma->SetMarkerStyle(20);
  gdEdxPiVsFeSigma->GetXaxis()->SetTitle("#sigma Fe55 [%]");
  gdEdxPiVsFeSigma->GetYaxis()->SetTitle("d#it{E}/d#it{x} resolution #pi [%]");

  TGraph *gdEdxEleVsFeSigma = new TGraph();
  gdEdxEleVsFeSigma->SetMarkerStyle(20);
  gdEdxEleVsFeSigma->GetXaxis()->SetTitle("#sigma Fe55 [%]");
  gdEdxEleVsFeSigma->GetYaxis()->SetTitle("d#it{E}/d#it{x} resolution e [%]");
*/

  const int nBunches = 1000;
  TH1F *hQmaxBunch[nBunches];

  TH1F *hdEdxPiBunch[nBunches];

  TF1 *MPVfit[nBunches];
  float MPVBunch [nBunches];
  float ChisquareBunch [nBunches];

  TH2D *hExcludeHisto = new TH2D ("hExcludeHisto","; Pad row; Pad", 63,0,63,100,-50,50);


  for (int i = 0; i < nBunches; ++i) {
    hQmaxBunch[i] = new TH1F(Form("hQmaxBunch_%i",i),Form("; Q_{max} Bunch %i; Counts",i),600,0,600);
    hdEdxPiBunch[i] = new TH1F(Form("hdEdxPi_Bunch_%i", i), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 15);
    MPVfit[i] = new TF1(Form("MPVfit_Bunch_%i",i),"landau");
  }

  TGraphErrors *gdEdxBunch = new TGraphErrors();
  gdEdxBunch->SetMarkerStyle(20);
  gdEdxBunch->GetXaxis()->SetTitle("bunch number");
  gdEdxBunch->GetYaxis()->SetTitle("d#it{E}/d#it{x} #pi^{-} (a.u.)");

  TGraph *gRunNrVsBunchNr = new TGraph();
  gRunNrVsBunchNr->SetMarkerStyle(20);
  gRunNrVsBunchNr->GetXaxis()->SetTitle("bunch number");
  gRunNrVsBunchNr->GetYaxis()->SetTitle("run number");

  /* =============================================================================
   * ================== Loop over events and apply all cuts ======================
   * ======================= meanwhile fill histograms ===========================
   * ============================================================================= */

  int Tracks = 0;
  float CherenkovValue = 0.;
  int runNr = 0;
  bool newBunch = false;
  int HVSettings = 0;
  int powerSupply = 0;
  int loopcounter = 0;
  int loopcounter2 = 0;
  bool setgainmapS1 = false;
  int ncl = 0;
  int ncledge = 0;
  int pioncounter = 0;
  int bunchcounter = 0;
  int realNBunches = 0;
  int bunchsize = 4000;
//  int runcounter = 0;

  std::vector<TrackTPC> *vecEvent = 0;

  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;


/*=============================================== Loop over runs =========================================================*/
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");

    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);

    CherenkovValue = 0;
    runNr = 0;
    nclRun = 0;


    TString CurrentFileName = chainEntries->At(ifile)->GetTitle();
    TObjArray *SplitFileName = CurrentFileName.Tokenize("/");
    TString runTString = SplitFileName->At(5)->GetName();

    std::string runString = "";
    runString.append(runTString);
    char *run = new char[runString.length() + 1];
    std::strcpy(run, runString.c_str());

//    hQmaxPi[runcounter] = new TH1F(Form("hQmaxPi_%s",run),"; Q_{max}Pi; Counts",600,0,600);

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        if (vecEvent->size() == 1) {
          nclRun += clCont.size();
        }
      }
    }
/*=============================================== Loop over events =========================================================*/

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);
      int NTracks = vecEvent->size();
      if (NTracks != 1) continue;
      CherenkovValue = Header.cherenkovValue;
      if (CherenkovValue >= CherCutLow) {continue;}
      runNr = Header.run;
      gRunNrVsBunchNr->SetPoint(bunchcounter, bunchcounter, runNr);
      powerSupply = Header.powerSupply;
      HVSettings = Header.HVSettings + 1;
      if (HVSettings == 11) {
        HVSettings = 1;
      }



      Tracks += NTracks;

      setting = 0;
      if (runNr >= 255) {setting = 1;}
      if (setting == 1) {
        setgainmapS1 = true;
        ++loopcounter2;
      }
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        if ((loopcounter == 0) && usegainmap) {
          trackObject.setGainMap(GainMap, setting);
          std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
        }
        else if (loopcounter2 == 1 && usegainmap) {
          trackObject.setGainMap(GainMap, setting);
          std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
        }
        ++loopcounter;
        ncl = clCont.size();
        bool isok = true;
        ncledge = 0;

        if (ncl < nclCut) continue;
        if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;	 // PID via Cherenkov

        for (auto &clusterObject : clCont) {
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
          float timeMean = clusterObject.getTimeMean();
          if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;} // cut on time
          if (!InsideEdge(row, GetCPad(clusterObject), setting, nEdgeCut)){++ncledge;}
        } // end of cluster loop
        if (((float(ncledge)/float(ncl)) > fracEdgeCl)) continue; // cut on detector geometry
        if (isok == true){

          if (newBunch == true) {
            ++bunchcounter;
            pioncounter = 0;
            newBunch = false;
          }

          if (CherenkovValue < CherCutLow) {
            for (auto &clusterObject : clCont) {
              float Qmax = clusterObject.getQmax();
              DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
              float row = pos.getPadSecPos().getPadPos().getRow();
              float pad = pos.getPadSecPos().getPadPos().getPad();
              float cpad = GetCPad(clusterObject);
              if (setting == 0) {
                if (row < 2 || row > 57){continue;}
                else if (row == 31 || row == 32) {continue;}
                if (!InsideEdge(row, GetCPad(clusterObject), setting, nEdgeCut)){continue;}
                float norm = gainmapS0->getValue(roc, row, pad);
                float QmaxCorr = Qmax/norm;
                hQmaxBunch[bunchcounter]->Fill(QmaxCorr);
//                hQmaxPiRuns->Fill(runNr,QmaxCorr);
//                hQmaxPi[runcounter]->Fill(QmaxCorr);
              }
              else if (setting == 1) {
                if (row < 2 || row > 60){continue;}
                else if (row > 1 && row < 61) {
                  if (cpad == 0 || cpad == -1 || cpad == -2) {continue;}
                }
                if (!InsideEdge(row, GetCPad(clusterObject), setting, nEdgeCut)){continue;}
                float norm = gainmapS1->getValue(roc, row, pad);
                float QmaxCorr = Qmax/norm;
                hQmaxBunch[bunchcounter]->Fill(QmaxCorr);
//                hQmaxPiRuns->Fill(runNr,QmaxCorr);
//                hQmaxPi[runcounter]->Fill(QmaxCorr);
              }
            }
          }
          ++pioncounter;
          if (pioncounter == bunchsize) {newBunch = true;} //{gRunNrVsBunchNr->SetPoint(bunchcounter, bunchcounter, runNr); ++bunchcounter; pioncounter = 0;}
        } // end of good track loop

      } // end of track loop
    } // end of event loop

    if (newBunch == false) {newBunch = true;};

  } // end of run loop
  const float fracLow=0.75;
  const float fracHigh=0.25;
  int bin1=0,bin2=0;

  for (int i = 0; i <= bunchcounter; ++i) {
    GetBinMinMax(hQmaxBunch[i],fracLow,fracHigh,bin1,bin2);
    hQmaxBunch[i]->Fit(Form("MPVfit_Bunch_%i",i),"","",hQmaxBunch[i]->GetXaxis()->GetBinLowEdge(bin1),hQmaxBunch[i]->GetXaxis()->GetBinUpEdge(bin2));
    MPVBunch[i] = MPVfit[i]->GetParameter(1);
    ChisquareBunch[i] = MPVfit[i]->GetChisquare()/MPVfit[i]->GetNDF();
  }

  realNBunches = bunchcounter + 1;

  std::cout<<std::endl<<"number of bunches: "<< realNBunches <<std::endl;

  bunchcounter = 0;
  pioncounter = 0;

  loopcounter = 0;
  loopcounter2 = 0;
  setting = 0;

  newBunch = false;

  /*=============================================== Loop over runs =========================================================*/
    for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
      TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
      cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" again: "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
      TTree *tree = (TTree*)TreeFile->Get("events");

      vecEvent=0;
      tree->SetBranchAddress("Tracks", &vecEvent);
      tree->SetBranchAddress("header", &Header);

      CherenkovValue = 0;
      runNr = 0;
      nclRun = 0;

      TString CurrentFileName = chainEntries->At(ifile)->GetTitle();
      TObjArray *SplitFileName = CurrentFileName.Tokenize("/");
      TString runTString = SplitFileName->At(5)->GetName();

      std::string runString = "";
      runString.append(runTString);
      char *run = new char[runString.length() + 1];
      std::strcpy(run, runString.c_str());

  //    hQmaxPi[runcounter] = new TH1F(Form("hQmaxPi_%s",run),"; Q_{max}Pi; Counts",600,0,600);

      for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
        tree->GetEntry(iev);
        for (auto trackObject : *vecEvent) {
          std::vector<Cluster> clCont;
          trackObject.getClusterVector(clCont);
          if (vecEvent->size() == 1) {
            nclRun += clCont.size();
          }
        }
      }

      int eventcounter = 0;
  /*=============================================== Loop over events =========================================================*/

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);
      int NTracks = vecEvent->size();
      if (NTracks != 1) continue;
      CherenkovValue = Header.cherenkovValue;
      if (CherenkovValue >= CherCutLow) {continue;}
      runNr = Header.run;
      powerSupply = Header.powerSupply;
      HVSettings = Header.HVSettings + 1;
      if (HVSettings == 11 && powerSupply == 1) {
        HVSettings = 1;
      }

      Tracks += NTracks;

//      std::cout<<std::endl<<"vor setting gainmap"<< setting <<std::endl;

      setting = 0;
      if (runNr >= 255) {setting = 1;}
      if (setting == 1) {
        setgainmapS1 = true;
        ++loopcounter2;
      }
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        if ((loopcounter == 0) && usegainmap) {
          trackObject.setGainMap(GainMap, setting);
          std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
        }
        else if (loopcounter2 == 1 && usegainmap) {
          trackObject.setGainMap(GainMap, setting);
          std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
        }
        ++loopcounter;
//        std::cout<<std::endl<<"loopcounter hoch" <<std::endl;

        ncl = clCont.size();
//        std::cout<<std::endl<<"ncl: " << ncl << std::endl;
//        std::cout<<std::endl<<"cherenkovValue: " << CherenkovValue << std::endl;


        bool isok = true;
        ncledge = 0;

        if (ncl < nclCut) continue;
        if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;	 // PID via Cherenkov

        for (auto &clusterObject : clCont) {
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
          float timeMean = clusterObject.getTimeMean();
          if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;} // cut on time
          if (!InsideEdge(row, GetCPad(clusterObject), setting, nEdgeCut)){++ncledge;}
        } // end of cluster loop
        if (((float(ncledge)/float(ncl)) > fracEdgeCl)) continue; // cut on detector geometry

//        std::cout<<std::endl<<"ende track cuts, setting: "<< setting <<std::endl;

        if (isok == true){
//          std::cout<<std::endl<<"anfang good track loop, setting:"<< setting <<std::endl;

          if (newBunch == true) {
            ++bunchcounter;
            pioncounter = 0;
            newBunch = false;
          }

          float dEdxTot = 0;
          float dEdxMax = 0;

//          std::cout<<std::endl<<"vor dEdx Berechnung"<< std::endl;

          if (setting == 0) {
            dEdxTot = trackObject.getTruncatedMean(runNr,0,.7,1,true,true,true,hExcludeHisto,MPVBunch[bunchcounter]);
            dEdxMax = trackObject.getTruncatedMean(runNr,0,.7,0,true,true,true,hExcludeHisto,MPVBunch[bunchcounter]);
          }

          else if (setting == 1) {
            dEdxTot = trackObject.getTruncatedMean(runNr,0,.7,1,true,true,true,hExcludeHisto,MPVBunch[bunchcounter]);
            dEdxMax = trackObject.getTruncatedMean(runNr,0,.7,0,true,true,true,hExcludeHisto,MPVBunch[bunchcounter]);
          }
//          std::cout<<std::endl<<"nach dEdx Berechnung"<< std::endl;

          if (CherenkovValue < CherCutLow){
            hdEdxPiBunch[bunchcounter]->Fill(dEdxTot);
//            hdEdxPionTotAll[HVSettings]->Fill(dEdxTot);
          }
          else if (CherenkovValue > CherCutHigh){
            std::cout << "\n\n\n Hier sollten keine Elektronen sein!!! \n\n\n";
//            hdEdxEleTotAll[HVSettings]->Fill(dEdxTot);
          }

          ++pioncounter;
          if (pioncounter == bunchsize) {newBunch = true;}
        }
      }
    }

    std::cout << "\nevent " << eventcounter << " fertig \n";
    ++eventcounter;

    if (newBunch == false) {newBunch = true;}




//    g->WriteObject(hQmaxPi[runcounter], Form("hQmaxPi_%s",run));
//    ++runcounter;
  } // end of run loop

//  TCanvas *dEdx[realNBunches];
  TF1 *pionfit[realNBunches];
//  TF1 *electronfit[realNBunches];
//  TPaveText *pave[realNBunches];


  for (int i=0; i<realNBunches; ++i) {
//    dEdx[i] = new TCanvas(Form("dEdx_HV_%i", i-1), Form("dEdx Qtot HV %i", i-1));
    pionfit[i] = new TF1(Form("pionfit_Bunch_%i",i),"gaus"); //,hdEdxPiBunch[i]->GetXaxis()->GetXmin(),hdEdxPiBunch[i]->GetXaxis()->GetXmax());
//    electronfit[i] = new TF1(Form("electronfit_HV_%i",i),"gaus",hdEdxEleTotAll[i]->GetXaxis()->GetXmin(),hdEdxEleTotAll[i]->GetXaxis()->GetXmax());
//    pave[i] = new TPaveText(0.6,.7,.9,.9,"NDC");
//    pave[i]->SetBorderSize(1);
//    pave[i]->SetFillColor(10);
//    hdEdxEleTotAll[i]->SetLineColor(kRed);
//    hdEdxPionTotAll[i]->SetLineColor(kBlue);
    pionfit[i]->SetLineColor(kBlue);
//    electronfit[i]->SetLineColor(kRed);

    const float fracLow=0.2;
    const float fracHigh=0.2;
    int bin1=0,bin2=0;

    GetBinMinMax(hdEdxPiBunch[i],fracLow,fracHigh,bin1,bin2);
    hdEdxPiBunch[i]->Fit(Form("pionfit_Bunch_%i", i),"","",hdEdxPiBunch[i]->GetXaxis()->GetBinLowEdge(bin1),hdEdxPiBunch[i]->GetXaxis()->GetBinUpEdge(bin2));

//    GetBinMinMax(hdEdxEleTotAll[i],frac,bin1,bin2);
//    hdEdxEleTotAll[i]->Fit(Form("electronfit_HV_%i",i),"","",hdEdxEleTotAll[i]->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTotAll[i]->GetXaxis()->GetBinUpEdge(bin2));

    //alternative fit
    //hdEdxPionTotAll[i]->Fit(Form("pionfit_HV_%i",i));
    //hdEdxEleTotAll[i]->Fit(Form("electronfit_HV_%i",i));


//    hdEdxEleTotAll[i]->GetFunction(Form("electronfit_HV_%i",i))->SetLineColor(kRed);
//    hdEdxPiBunch[i]->GetFunction(Form("pionfit_Bunch_%i",i))->SetLineColor(kBlue);

//    hdEdxPionTotAll[i]->Draw();
//    hdEdxEleTotAll[i]->Draw("same");

    pionmeanTot = pionsigmaTot = electronmeanTot = electronsigmaTot = pionres = electronres
        = pionsigmaerr = pionmeanerr = electronsigmaerr = electronmeanerr = separationpower = 0;
    pionmeanTot = pionfit[i]->GetParameter(1);
    pionsigmaTot = pionfit[i]->GetParameter(2);
//    electronmeanTot = electronfit[i]->GetParameter(1);
//    electronsigmaTot = electronfit[i]->GetParameter(2);

    if (pionmeanTot != 0) {
      pionres = pionsigmaTot/pionmeanTot;
    }
//    if (electronmeanTot != 0) {
//      electronres = electronsigmaTot/electronmeanTot;
//    }
    pionsigmaerr = pionfit[i]->GetParError(2);
    pionmeanerr = pionfit[i]->GetParError(1);
//    electronsigmaerr = electronfit[i]->GetParError(2);
//    electronmeanerr = electronfit[i]->GetParError(1);

//    if (pionsigmaTot != 0 || electronsigmaTot != 0) {
//      separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);
//    }

    gdEdxBunch->SetPoint(i,i,pionmeanTot);
    gdEdxBunch->SetPointError(i,0,pionmeanerr);


/*
    float pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
    float electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));

    float separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
        pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));

    float elechisquareTot = electronfit[i]->GetChisquare();
    float eleNDFTot = electronfit[i]->GetNDF();
    float pionchisquareTot = pionfit[i]->GetChisquare();
    float pionNDFTot = pionfit[i]->GetNDF();
*/

//    pave[i]->AddText(Form("HV setting %i", i-1));
//    pave[i]->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanTot,electronsigmaTot, electronsigmaTot/electronmeanTot*100));
//    pave[i]->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanTot,pionsigmaTot,pionsigmaTot/pionmeanTot*100));
//    pave[i]->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanTot-pionmeanTot)/((electronsigmaTot+pionsigmaTot)/2.)));
//    pave[i]->Draw("same");

//    dEdx[i]->Update();
//    dEdx[i]->Print(Form("%s/dEdxQtot_merged_HV_%i.png", OutputPath, i-1));
  }

  //TrackAna->Fill();
  g->WriteObject(gdEdxBunch, "dEdxMeanVsBunchNr");
  g->WriteObject(gRunNrVsBunchNr, "gRunNrVsBunchNr");
  for (int i=0; i<realNBunches; ++i) {
    g->WriteObject(hQmaxBunch[i], Form("hQmaxBunch_%i",i));
    g->WriteObject(hdEdxPiBunch[i], Form("hdEdxPiBunch_%i",i));
  }
//  g->WriteObject(gSepVsHVSet, "gSepVsHVSet");
//  g->WriteObject(gSepVsFeSigma, "gSepVsFeSigma");
//  g->WriteObject(gSepVsIBF, "gSepVsIBF");
//  g->WriteObject(gdEdxPiVsFeSigma, "gdEdxPiVsFeSigma");
//  g->WriteObject(gdEdxEleVsFeSigma, "gdEdxEleVsFeSigma");
//  g->WriteObject(hQmaxPiRuns, "hQmaxPiRuns");
//  g->WriteObject(hMPVRuns, "hMPVRuns");
//  g->WriteObject(hChisquareMPVfit, "hChisquareMPVfit");
//  for (int i=0; i<realNBunches; ++i) {
//    g->WriteObject(hdEdxPionTotAll[i], Form("hdEdxPiTot_HV_%i", i-1));
//    g->WriteObject(hdEdxEleTotAll[i], Form("hdEdxEleTot_HV_%i", i-1));
//  }
  delete g;
}

