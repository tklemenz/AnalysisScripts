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

#include <math.h>


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
  const float contMax=hist->GetBinContent(binMax);
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

int checkCuts(TString FileList, TString GainMap="", const char *OutputPath="/home/tom/myana/Results")
{
  EventHeader Header;

  bool usegainmap = true;
  if (GainMap == "") {
    usegainmap = false;
  }


//  std::vector<int> nClCut { 16,32,47,54 };
//  std::vector<int> nEdgeCut { 0,1,2,3 };
//  std::vector<float> fracEdgeCl { .2,.3,.4 };
//  std::vector<int> timeMeanLow { 3,4,5,6 };
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMeanHigh = 20;
  std::vector<int> nClCut { 16,32,47,54 };
  std::vector<int> nEdgeCut { 0,1,2,3 };
  std::vector<float> fracEdgeCl { .2,.3,.4 };
  std::vector<int> timeMeanLow { 0,1,2,3,4,5,6 };
  std::vector<float> truncHigh { .6,.65,.7,.75,.8 };
  std::vector<float> truncLow { .0,.05,.1 };

//  std::vector<int> nClCut { 54 };
//  std::vector<int> nEdgeCut { 2 };
//  std::vector<float> fracEdgeCl { .3 };
//  std::vector<int> timeMeanLow { 0,1,2,3,4,5,6 };
//  std::vector<float> truncHigh { .7 };
//  std::vector<float> truncLow { .0 };

  int permutations = nClCut.size() * nEdgeCut.size() * fracEdgeCl.size() * timeMeanLow.size() * truncLow.size() * truncHigh.size();


  float electronmeanTot, electronsigmaTot, pionmeanTot, pionsigmaTot, electronmeanerr, electronsigmaerr, pionmeanerr, pionsigmaerr, electronmeanMax,
      electronsigmaMax, pionmeanMax, pionsigmaMax, electronmeanerrMax, electronsigmaerrMax, pionmeanerrMax, pionsigmaerrMax, pionres, pionresMax, electronres,
      electronresMax, separationpower, separationpowerMax, nclRun, CherenkovValue;
  int pions, electrons, pioncl, elecl, nclTrack, nclEvent, runNr,
      loopcounter, loopcounter2, setting, permutationcounter;
  bool setgainmapS1;

  electronmeanTot = electronsigmaTot = pionmeanTot = pionsigmaTot = electronmeanerr = electronsigmaerr = pionmeanerr = pionsigmaerr = electronmeanMax =
      electronsigmaMax = pionmeanMax = pionsigmaMax = electronmeanerrMax = electronsigmaerrMax = pionmeanerrMax = pionsigmaerrMax = pionres =
      pionresMax = electronres = electronresMax = separationpower
      = separationpowerMax = nclRun = CherenkovValue = 0.;
  pions = electrons = pioncl = elecl = nclTrack = nclEvent = runNr =
      loopcounter = loopcounter2 = setting = permutationcounter = 0;
  setgainmapS1 = false;

  TH2D *hSeperation    = new TH2D ("Separationpower",";run number; seperation power",350,0,350,100,0,8);
  TH2D *hExcludeHisto  = new TH2D ("hExcludeHisto","; Pad row; Pad", 63,0,63,100,-50,50);

  TFile *OutFile = new TFile(Form("%s/checkCutsOut.root", OutputPath), "recreate");

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
    fChain.Add(file);
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }

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
    permutationcounter = 0;

    TString CurrentFileName = chainEntries->At(ifile)->GetTitle();
    TObjArray *SplitFileName = CurrentFileName.Tokenize("/");
    TString runTString = SplitFileName->At(5)->GetName();

    std::string runString = "";
    runString.append(runTString);
    char *run = new char[runString.length() + 1];
    std::strcpy(run, runString.c_str());

    TH2F *dEdxPiTotCut   = new TH2F (Form("dEdxPiTotCut_%s", run),"; cut; dE/dx res Pi tot",permutations,0,permutations, 400,0,0.25);
    TH2F *dEdxPiMaxCut   = new TH2F (Form("dEdxPiMaxCut_%s", run),"; cut; dE/dx res Pi max",permutations,0,permutations, 400,0,0.25);
    TH2F *dEdxEleTotCut  = new TH2F (Form("dEdxEleTotCut_%s", run),"; cut; dE/dx res Ele tot",permutations,0,permutations, 400,0,0.25);
    TH2F *dEdxEleMaxCut  = new TH2F (Form("dEdxEleMaxCut_%s", run),"; cut; dE/dx res Ele max",permutations,0,permutations, 400,0,0.25);
    TH2F *SepPowerTotCut = new TH2F (Form("SepPowerTotCut_%s", run),"; cut; separationpower tot",permutations,0,permutations, 100,0,5);
    TH2F *SepPowerMaxCut = new TH2F (Form("SepPowerMaxCut_%s", run),"; cut; separationpower max",permutations,0,permutations, 100,0,5);


    /// create a TTree for each run with the runNr as ID
//    TTree *TrackAna = new TTree(Form("TrackAna_%s",run), Form("dEdxAna_%s", run));
    /// loop over nclCut vector
    for (int i = 0; i < nClCut.size(); ++i) {
      for (int j = 0; j < nEdgeCut.size(); ++j) {
        for (int k = 0; k < fracEdgeCl.size(); ++k) {
          for (int l = 0; l < timeMeanLow.size(); ++l) {
            for (int m = 0; m < truncLow.size(); ++m) {
              for (int n = 0; n < truncHigh.size(); ++n) {
//                TrackAna->Branch(Form("PidEdxResQ_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &pionres, "pionres/F");
//                TrackAna->Branch(Form("PidEdxResQmax_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &pionresMax, "pionresMax/F");
//                TrackAna->Branch(Form("EledEdxResQ_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &electronres, "electronres/F");
//                TrackAna->Branch(Form("EledEdxResQmax_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &electronresMax, "electronresMax/F");
//                TrackAna->Branch(Form("SepPowerQ_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &separationpower, "separationpower/F");
//                TrackAna->Branch(Form("SepPowerQmax_ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]), &separationpowerMax, "separationpowerMax/F");


                TH1F *hdEdxEleTot             = new TH1F("hdEdxEleTot", "; d#it{E}/d#it{x} (a.u.); # counts", 400, 0, 800);
                TH1F *hdEdxPionTot            = new TH1F("hdEdxPionTot", "; d#it{E}/d#it{x} (a.u.); # counts", 400, 0, 800);

                TH1F *hdEdxEleMax             = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); # counts", 300, 0, 360);
                TH1F *hdEdxPionMax            = new TH1F("hdEdxPionMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); # counts", 300, 0, 360);

                pionmeanTot = 0;
                pionsigmaTot = 0;
                electronmeanTot = 0;
                electronsigmaTot = 0;
                pionres = 0;
                electronres = 0;
                separationpower = 0;
                pionmeanMax = 0;
                pionsigmaMax = 0;
                electronmeanMax = 0;
                electronsigmaMax = 0;
                pionresMax = 0;
                electronresMax = 0;
                separationpowerMax = 0;

                for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
                  tree->GetEntry(iev);
                  CherenkovValue = Header.cherenkovValue;
                  runNr = Header.run;
                  int NTracks = vecEvent->size();
                  if (NTracks != 1) continue;
                  setting = 0;
                  if (runNr >= 250) {setting = 1;}
                  if (setting == 1) {
                    setgainmapS1 = true;
                    ++loopcounter2;
                  }
                  for (auto trackObject : *vecEvent) {
                    std::vector<Cluster> clCont;
                    trackObject.getClusterVector(clCont);
                    /// ______________________________________________
                    /// set proper gainmap
                    if ((loopcounter == 0) && usegainmap) {
                      trackObject.setGainMap(GainMap, setting);
                      std::cout<<std::endl<<"set gainmap for setting "<<setting<<std::endl;
                    }
                    else if (loopcounter2 == 1 && usegainmap) {
                      trackObject.setGainMap(GainMap, setting);
                      std::cout<<std::endl<<"set gainmap for setting "<<setting<<std::endl;
                    }
                    ++loopcounter;
                    /// ______________________________________________
                    int ncl = clCont.size();
                    bool isok = true;
                    int ncledge = 0;
                    if (ncl < nClCut[i]) continue;
                    if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;	 // PID via Cherenkov
                    for (auto &clusterObject : clCont) {															             // make cuts on cluster properties
                      DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
                      float row = pos.getPadSecPos().getPadPos().getRow();
                      float timeMean = clusterObject.getTimeMean();
                      if (timeMean < timeMeanLow[l] || timeMean > timeMeanHigh){isok = false; break;} // cut on time
                      if (!InsideEdge(row, GetCPad(clusterObject), setting, nEdgeCut[j])){++ncledge;}
                    } // end of loop over clusters
                    if (((float(ncledge)/float(ncl)) > fracEdgeCl[k])) continue; // cut on detector geometry
                    if (isok == true){
                      //std::cout<<std::endl<<"good track found"<<std::endl;
                      float dEdxTot = 0;
                      float dEdxMax = 0;

                      if (setting == 0) {
                        dEdxTot = trackObject.getTruncatedMean(runNr,truncLow[m],truncHigh[n],1,false,false,true,hExcludeHisto);
                        dEdxMax = trackObject.getTruncatedMean(runNr,truncLow[m],truncHigh[n],0,false,false,true,hExcludeHisto);
                      }

                      else if (setting == 1) {
                        dEdxTot = trackObject.getTruncatedMean(runNr,truncLow[m],truncHigh[n],1,false,false,true,hExcludeHisto);
                        dEdxMax = trackObject.getTruncatedMean(runNr,truncLow[m],truncHigh[n],0,false,false,true,hExcludeHisto);
                      }

                      if (CherenkovValue < CherCutLow){
                        hdEdxPionTot->Fill(dEdxTot);
                        hdEdxPionMax->Fill(dEdxMax);
                      } // end of pion loop
                      else if (CherenkovValue > CherCutHigh){
                        hdEdxEleTot->Fill(dEdxTot);
                        hdEdxEleMax->Fill(dEdxMax);
                      } // end of electron loop
                    } // end of 'track is ok'

                  } // end of loops over tracks
                } // end of loop over events
                TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot->GetXaxis()->GetXmin(),hdEdxPionTot->GetXaxis()->GetXmax());
                TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot->GetXaxis()->GetXmin(),hdEdxEleTot->GetXaxis()->GetXmax());
                const Float_t frac=0.2;
                Int_t bin1=0,bin2=0;

                GetBinMinMax(hdEdxPionTot,frac,bin1,bin2);
                hdEdxPionTot->Fit("pionfit","","",hdEdxPionTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot->GetXaxis()->GetBinUpEdge(bin2));
                GetBinMinMax(hdEdxEleTot,frac,bin1,bin2);
                hdEdxEleTot->Fit("electronfit","","",hdEdxEleTot->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTot->GetXaxis()->GetBinUpEdge(bin2));

                pionmeanTot = pionfit->GetParameter(1);
                pionsigmaTot = pionfit->GetParameter(2);
                electronmeanTot = electronfit->GetParameter(1);
                electronsigmaTot = electronfit->GetParameter(2);

                pionres = pionsigmaTot/pionmeanTot;
                electronres = electronsigmaTot/electronmeanTot;
                separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);

                TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
                TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());

                GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
                hdEdxPionMax->Fit("pionfitMax","","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
                GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
                hdEdxEleMax->Fit("electronfitMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

                pionmeanMax = pionfitMax->GetParameter(1);
                pionsigmaMax = pionfitMax->GetParameter(2);
                electronmeanMax = electronfitMax->GetParameter(1);
                electronsigmaMax = electronfitMax->GetParameter(2);

                pionresMax = pionsigmaMax/pionmeanMax;
                electronresMax = electronsigmaMax/electronmeanMax;
                separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);


                dEdxPiTotCut->Fill(permutationcounter, pionres);
                dEdxPiMaxCut->Fill(permutationcounter, pionresMax);
                dEdxEleTotCut->Fill(permutationcounter, electronres);
                dEdxEleMaxCut->Fill(permutationcounter, electronresMax);
                SepPowerTotCut->Fill(permutationcounter, separationpower);
                SepPowerMaxCut->Fill(permutationcounter, separationpowerMax);

                dEdxPiTotCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));
                dEdxPiMaxCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));
                dEdxEleTotCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));
                dEdxEleMaxCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));
                SepPowerTotCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));
                SepPowerMaxCut->GetXaxis()->SetBinLabel(permutationcounter+1, Form("ncl%i_edge%i_frac%.1f_time%i_tl%.2f_th%.2f", nClCut[i], nEdgeCut[j], fracEdgeCl[k], timeMeanLow[l], truncLow[m], truncHigh[n]));


//                TrackAna->Fill();
                OutFile->cd();
//                TrackAna->Write();
                hSeperation->Fill(runNr,separationpower);
                delete pionfitMax;
                delete electronfitMax;
                delete pionfit;
                delete electronfit;
                delete hdEdxEleMax;
                delete hdEdxEleTot;
                delete hdEdxPionMax;
                delete hdEdxPionTot;

                ++permutationcounter;

              } //end of truncHigh loop
            } // end of truncLow loop
          } // end of timeMeanLow loop
        } // end of fracEdgeCut loop
      } // end of nEdgeCut loop
    } // end of nClCut loop
//    OutFile->cd();
//    TrackAna->Write();
    dEdxPiTotCut->Write();
    dEdxPiMaxCut->Write();
    dEdxEleTotCut->Write();
    dEdxEleMaxCut->Write();
    SepPowerTotCut->Write();
    SepPowerMaxCut->Write();
    delete dEdxPiTotCut;
    delete dEdxPiMaxCut;
    delete dEdxEleTotCut;
    delete dEdxEleMaxCut;
    delete SepPowerTotCut;
    delete SepPowerMaxCut;
//    delete TrackAna;
  } // end of loop over individual runs

  hSeperation->Write();
  OutFile->Write();
  OutFile->Close();

  return 0;
}
