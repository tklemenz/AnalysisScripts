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

void AnalyseRuns(TString FileList, TString GainMap="", const char *OutputPath="/home/tom/myana/Results")
{
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
      electronsigmaMax, pionmeanMax, pionsigmaMax, electronmeanerrMax, electronsigmaerrMax, pionmeanerrMax, pionsigmaerrMax, pionres, pionresMax, electronres,
      electronresMax, separationpower, separationpowerMax, nclRun, pionreserror, electronreserror, separationpowererr, elechisquareTot, eleNDFTot, pionchisquareTot,
      pionNDFTot, pionreserrorMax, electronreserrorMax, separationpowererrMax, elechisquareMax, eleNDFMax, pionchisquareMax, pionNDFMax, CherenkovValue, dEdxTot, dEdxMax,
      gain, IB, sigmaFe, IBErr, sigmaFeErr;
  int nPions, nElectrons, pioncl, elecl, nclTrack, nclEvent, runNr, beamMomentum, powerSupply, HVSetting, trigger, dataType, driftFieldStrenght, runType;

  electronmeanTot = electronsigmaTot = pionmeanTot = pionsigmaTot = electronmeanerr = electronsigmaerr = pionmeanerr = pionsigmaerr = electronmeanMax =
      electronsigmaMax = pionmeanMax = pionsigmaMax = electronmeanerrMax = electronsigmaerrMax = pionmeanerrMax = pionsigmaerrMax = pionres =
      pionresMax = electronres = electronresMax = separationpower = separationpowerMax = nclRun = pionreserror = electronreserror = separationpowererr =
      elechisquareTot = eleNDFTot = pionchisquareTot = pionNDFTot = pionreserrorMax = electronreserrorMax = separationpowererrMax = elechisquareMax =
      eleNDFMax = pionchisquareMax = pionNDFMax = CherenkovValue = dEdxTot = dEdxMax = gain = IB = sigmaFe = IBErr = sigmaFeErr = 0.;
  nPions = nElectrons = pioncl = elecl = nclTrack = nclEvent = beamMomentum = powerSupply = HVSetting = trigger = dataType = driftFieldStrenght = runType = 0;

  TFile *OutFile = new TFile(Form("%s/anaOut_tree.root", OutputPath), "recreate");
  TTree *TrackAna = new TTree("TrackAna","dEdxAna");

  TrackAna->Branch("PidEdxResQtot", &pionres, "pionres/F");
  TrackAna->Branch("PidEdxResQmax", &pionresMax, "pionresMax/F");
  TrackAna->Branch("EledEdxResQtot", &electronres, "electronres/F");
  TrackAna->Branch("EledEdxResQmax", &electronresMax, "electronresMax/F");
  TrackAna->Branch("SepPowerQtot", &separationpower, "separationpower/F");
  TrackAna->Branch("SepPowerQmax", &separationpowerMax);
  TrackAna->Branch("posEleQtot", &electronmeanTot, "electronmeanTot/F");
  TrackAna->Branch("sigmaEleQtot", &electronsigmaTot, "electronsigmaTot/F");
  TrackAna->Branch("posPiQtot", &pionmeanTot, "pionmeanTot/F");
  TrackAna->Branch("sigmaPiQtot", &pionsigmaTot, "pionsigmaTot/F");
  TrackAna->Branch("posEleQtotErr", &electronmeanerr, "electronmeanerr/F");
  TrackAna->Branch("sigmaEleQtotErr", &electronsigmaerr, "electronsigmaerr/F");
  TrackAna->Branch("posPiQtotErr", &pionmeanerr, "pionmeanerr/F");
  TrackAna->Branch("sigmaPiQtotErr", &pionsigmaerr, "pionsigmaerr/F");
  TrackAna->Branch("posEleQmax", &electronmeanMax, "electronmeanMax/F");
  TrackAna->Branch("sigmaEleQmax", &electronsigmaMax, "electronsigmaMax/F");
  TrackAna->Branch("posPiQmax", &pionmeanMax, "pionmeanMax/F");
  TrackAna->Branch("sigmaPiQmax", &pionsigmaMax, "pionsigmaMax/F");
  TrackAna->Branch("posEleQmaxErr", &electronmeanerrMax, "electronmeanerrMax/F");
  TrackAna->Branch("sigmaEleQmaxErr", &electronsigmaerrMax, "electronsigmaerrMax/F");
  TrackAna->Branch("posPiQmaxErr", &pionmeanerrMax, "pionmeanerrMax/F");
  TrackAna->Branch("sigmaPiQmaxErr", &pionsigmaerrMax, "pionsigmaerrMax/F");
  TrackAna->Branch("PidEdxResQtotErr", &pionreserror, "pionreserror/F");
  TrackAna->Branch("PidEdxResQmaxErr", &pionreserrorMax, "pionreserrorMax/F");
  TrackAna->Branch("EledEdxResQtotErr", &electronreserror, "electronreserror/F");
  TrackAna->Branch("EledEdxResQmaxErr", &electronreserrorMax, "electronreserrorMax/F");
  TrackAna->Branch("SepPowerQtotErr", &separationpowererr, "separationpowererr/F");
  TrackAna->Branch("SepPowerQmaxErr", &separationpowererrMax, "separationpowererrMax/F");
  TrackAna->Branch("EleChisquareQtot", &elechisquareTot, "elechisquareTot/F");
  TrackAna->Branch("EleChisquareQmax", &elechisquareMax, "elechisquareMax/F");
  TrackAna->Branch("EleNDFQtot", &eleNDFTot, "eleNDFTot/F");
  TrackAna->Branch("EleNDFQmax", &eleNDFMax, "eleNDFMax/F");
  TrackAna->Branch("PiChisquareQtot", &pionchisquareTot, "pionchisquareTot/F");
  TrackAna->Branch("PiChisquareQmax", &pionchisquareMax, "pionchisquareMax/F");
  TrackAna->Branch("PiNDFQtot", &pionNDFTot, "pionNDFTot/F");
  TrackAna->Branch("PiNDFQmax", &pionNDFMax, "pionNDFMax/F");
  TrackAna->Branch("runNr", &runNr, "runNr/I");
  TrackAna->Branch("CherenkovValue", &CherenkovValue, "CherenkovValue/F");
  TrackAna->Branch("beamMomentum", &beamMomentum, "beamMomentum/I");
  TrackAna->Branch("powerSupply", &powerSupply, "powerSupply/I");
  TrackAna->Branch("HVSetting", &HVSetting, "HVSetting/I");
  TrackAna->Branch("trigger", &trigger, "trigger/I");
  TrackAna->Branch("dataType", &dataType, "dataType/I");
  TrackAna->Branch("driftFieldStrenght", &driftFieldStrenght, "driftFieldStrenght/I");
  TrackAna->Branch("runType", &runType, "runType/I");
  TrackAna->Branch("gain", &gain, "gain/F");
  TrackAna->Branch("IB", &IB, "IB/F");
  TrackAna->Branch("IBErr", &IBErr, "IBErr/F");
  TrackAna->Branch("sigmaFe", &sigmaFe, "sigmaFe/F");
  TrackAna->Branch("sigmaFeErr", &sigmaFeErr, "sigmaFeErr/F");
  TrackAna->Branch("nPions", &nPions, "nPions/I");
  TrackAna->Branch("nElectrons", &nElectrons, "nElectrons/I");


/*
  TrackAna->Branch("dEdxPiQtot",&dEdxTot, "dEdxTot/F");
  TrackAna->Branch("dEdxEleQtot",&dEdxTot, "dEdxTot/F");
  TrackAna->Branch("dEdxPiQmax",&dEdxMax, "dEdxMax/F");
  TrackAna->Branch("dEdxEleQmax",&dEdxMax, "dEdxMax/F");
*/

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
  int TrPerEv = 1;
  int nclCut = 47;
  float CherCutLow = 0.009;
  float CherCutHigh = 0.011;
  int timeMeanLow = 6;
  int timeMeanHigh = 20;
  float nclFracoutofCPad = .3;
  int excludeEdge = 1;

  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/

  TH2D *hCheckExclude = new TH2D("CheckExclude", "; Row; Pad", 63,0,63,100,-50,50);

  const int nRuns = 350;
  TH1F *hdEdxEleTot[nRuns];
  TH1F *hdEdxPionTot[nRuns];
  TH1F *hdEdxEleMax[nRuns];
  TH1F *hdEdxPionMax[nRuns];

  for (int i = 0; i < nRuns; ++i) {
    hdEdxEleTot[i]   = new TH1F(Form("hdEdxEleTot_%i", i), "; d#it{E}/d#it{x}_{Q_{tot}} (a.u.); Counts", 200, 0, 400);
    hdEdxPionTot[i]  = new TH1F(Form("hdEdxPionTot_%i", i), "; d#it{E}/d#it{x}_{Q_{tot}} (a.u.); Counts", 200, 0, 400);
    hdEdxEleMax[i]   = new TH1F(Form("hdEdxEleMax_%i", i), "; d#it{E}/d#it{x}_{Q_{max}} (a.u.); Counts", 100, 0, 120);
    hdEdxPionMax[i]  = new TH1F(Form("hdEdxPionMax_%i", i), "; d#it{E}/d#it{x}_{Q_{max}} (a.u.); Counts", 100, 0, 120);
  }

  /* =============================================================================
   * ================== Loop over events and apply all cuts ======================
   * ======================= meanwhile fill histograms ===========================
   * ============================================================================= */


    int loopcounter = 0;
    int loopcounter2 = 0;
    bool setgainmapS1 = false;
    int ncl = 0;
    int ncledge = 0;
    int runcounter = 0;


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
      beamMomentum = 0;
      powerSupply = 0;
      HVSetting = 0;
      trigger = 0;
      dataType = 0;
      runType = 0;
      driftFieldStrenght = 0;

      nPions = 0;
      nElectrons = 0;

  /*=============================================== Loop over events =========================================================*/

      for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
        tree->GetEntry(iev);

        int NTracks = vecEvent->size();

        if (NTracks != TrPerEv) continue;																// only one-track events
        CherenkovValue = Header.cherenkovValue;
        runNr = Header.run;
        beamMomentum = Header.beamMomentum;
        powerSupply = Header.powerSupply;
        HVSetting = Header.HVSettings;
        trigger = Header.trigger;
        dataType = Header.dataType;
        runType = Header.runType;
        driftFieldStrenght = Header.driftFieldStrength;

        if (HVSetting == -1) {
          gain = 2000;
          IB = 0.65;
          sigmaFe = 12.7;
        }
        else if (HVSetting == 0) {
          gain = 2000;
          IB = 1.12;
          sigmaFe = 10.8;
        }
        else if (HVSetting == 1) {
          gain = 2000;
          IB = 0.98;
          sigmaFe = 12.2;
        }
        else if (HVSetting == 2) {
          gain = 1500;
          IB = 1.12;
          sigmaFe = 12.1;
        }
        else if (HVSetting == 3) {
          gain = 1000;
          IB = 1.6;
          sigmaFe = 12.3;
        }
        else if (HVSetting == 4) {
          gain = 3000;
          IB = 0.92;
          sigmaFe = 12.0;
        }
        else if (HVSetting == 5) {
          gain = 2000;
          IB = 0.34;
          sigmaFe = 17.0;
        }
        else if (HVSetting == 6) {
          gain = 2000;
          IB = 0.51;
          sigmaFe = 13.8;
        }
        else if (HVSetting == 7) {
          gain = 2000;
          IB = 0.65;
          sigmaFe = 12.1;
        }
        else if (HVSetting == 8) {
          gain = 2000;
          IB = 0.89;
          sigmaFe = 10.4;
        }
        else if (HVSetting == 9) {
          gain = 2000;
          IB = 2.05;
          sigmaFe = 9.1;
        }
        else if (HVSetting == 10) {
          gain = 2000;
          IB = 1.12;
          sigmaFe = 10.8;
        }

        IBErr = 0.05*IB;
        sigmaFeErr = 0.05*sigmaFe;

        /// Set proper Gain Map
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
          if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;


          for (auto &clusterObject : clCont) {															// make cuts on cluster properties
            DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
            float row = pos.getPadSecPos().getPadPos().getRow();
  //          if (setting == 0) {
  //            if (row < 2 || row >= 57) continue;
  //          }
  //          else if (setting == 1) {
  //           if (row <= 2 || row >= 61) continue;
  //          }
            float pad = pos.getPadSecPos().getPadPos().getPad();
            float timeMean = clusterObject.getTimeMean();
            float cpad = GetCPad(clusterObject);

            if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle

            if (!InsideEdge(row, cpad, setting, 3)){++ncledge;}
          }

          if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;

          if (isok == true){

            dEdxTot = 0;
            dEdxMax = 0;

            dEdxTot = trackObject.getTruncatedMean(runNr, .0, .7, 1, false, false, true, hCheckExclude, 1, 2);
            dEdxMax = trackObject.getTruncatedMean(runNr, .0, .7, 0, false, false, true, hCheckExclude, 1, 2);

            if (CherenkovValue < CherCutLow) {
              hdEdxPionTot[runcounter]->Fill(dEdxTot);
              hdEdxPionMax[runcounter]->Fill(dEdxMax);
              ++nPions;
            }

            else if (CherenkovValue > CherCutHigh) {
              hdEdxEleTot[runcounter]->Fill(dEdxTot);
              hdEdxEleMax[runcounter]->Fill(dEdxMax);
              ++nElectrons;
            }

          } /// end of 'if true'
        } /// end of loop over tracks
      } /// end of loops over events

      if (hdEdxPionTot[runcounter]->GetEntries() > 0 && hdEdxEleTot[runcounter]->GetEntries() > 0){
        TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot[runcounter]->GetXaxis()->GetXmin(),hdEdxPionTot[runcounter]->GetXaxis()->GetXmax());
        TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot[runcounter]->GetXaxis()->GetXmin(),hdEdxEleTot[runcounter]->GetXaxis()->GetXmax());

        const Float_t frac=0.2;
        Int_t bin1=0,bin2=0;


        GetBinMinMax(hdEdxPionTot[runcounter],frac,bin1,bin2);
        hdEdxPionTot[runcounter]->Fit("pionfit","","",hdEdxPionTot[runcounter]->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot[runcounter]->GetXaxis()->GetBinUpEdge(bin2));
        GetBinMinMax(hdEdxEleTot[runcounter],frac,bin1,bin2);
        hdEdxEleTot[runcounter]->Fit("electronfit");

        pionmeanTot = pionfit->GetParameter(1);
        pionsigmaTot = pionfit->GetParameter(2);
        electronmeanTot = electronfit->GetParameter(1);
        electronsigmaTot = electronfit->GetParameter(2);

        pionres = pionsigmaTot/pionmeanTot;
        electronres = electronsigmaTot/electronmeanTot;
        pionsigmaerr = pionfit->GetParError(2);
        pionmeanerr = pionfit->GetParError(1);
        electronsigmaerr = electronfit->GetParError(2);
        electronmeanerr = electronfit->GetParError(1);


        separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);


        pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
        electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));

        separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
              pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));

        elechisquareTot = electronfit->GetChisquare();
        eleNDFTot = electronfit->GetNDF();
        pionchisquareTot = pionfit->GetChisquare();
        pionNDFTot = pionfit->GetNDF();

        TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax[runcounter]->GetXaxis()->GetXmin(),hdEdxPionMax[runcounter]->GetXaxis()->GetXmax());
        TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax[runcounter]->GetXaxis()->GetXmin(),hdEdxEleMax[runcounter]->GetXaxis()->GetXmax());


        GetBinMinMax(hdEdxPionMax[runcounter],frac,bin1,bin2);
        hdEdxPionMax[runcounter]->Fit("pionfitMax","","r",hdEdxPionMax[runcounter]->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax[runcounter]->GetXaxis()->GetBinUpEdge(bin2));
        GetBinMinMax(hdEdxEleMax[runcounter],frac,bin1,bin2);
        hdEdxEleMax[runcounter]->Fit("electronfitMax");

        pionmeanMax = pionfitMax->GetParameter(1);
        pionsigmaMax = pionfitMax->GetParameter(2);
        electronmeanMax = electronfitMax->GetParameter(1);
        electronsigmaMax = electronfitMax->GetParameter(2);

        pionresMax = pionsigmaMax/pionmeanMax;
        electronresMax = electronsigmaMax/electronmeanMax;
        pionsigmaerrMax = pionfitMax->GetParError(2);
        pionmeanerrMax = pionfitMax->GetParError(1);
        electronsigmaerrMax = electronfitMax->GetParError(2);
        electronmeanerrMax = electronfitMax->GetParError(1);


        separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);


        pionreserrorMax = sqrt(pow((pionsigmaerrMax/pionmeanMax),2)+pow(((pionsigmaMax*pionmeanerrMax)/pow(pionmeanMax,2)),2));
        electronreserrorMax = sqrt(pow((electronsigmaerrMax/electronmeanMax),2)+pow(((electronsigmaMax*electronmeanerrMax)/(electronmeanMax,2)),2));

        separationpowererrMax = sqrt(pow((2*electronmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+pow((2*pionmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+
              pow(((2*electronsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2)+pow(((2*pionsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2));

        elechisquareMax = electronfitMax->GetChisquare();
        eleNDFMax = electronfitMax->GetNDF();
        pionchisquareMax = pionfitMax->GetChisquare();
        pionNDFMax = pionfitMax->GetNDF();

        TrackAna->Fill();

        delete pionfit;
        delete electronfit;
        delete pionfitMax;
        delete electronfitMax;
      }
      delete hdEdxPionTot[runcounter];
      delete hdEdxPionMax[runcounter];
      delete hdEdxEleTot[runcounter];
      delete hdEdxEleMax[runcounter];

      ++runcounter;

    } /// end of loop over runs

    for (int i = runcounter+1; i < nRuns; ++i) {
      delete hdEdxPionTot[i];
      delete hdEdxPionMax[i];
      delete hdEdxEleTot[i];
      delete hdEdxEleMax[i];
    }


    OutFile->WriteObject(TrackAna, "TrackAna");


}
