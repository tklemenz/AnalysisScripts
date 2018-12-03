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

float unbiasedTruncatedMean(TH2D *histo, int runNr, o2::TPC::TrackTPC Track, CalDet<int> *HitMapS0, CalDet<int> *HitMapS1, int type=1, float low=.0, float high=.7)
{
  ROC roc(0);
  float norm = 1.;
  std::vector<Cluster> Clusters;
  Track.getClusterVector(Clusters);
  std::vector<float> values;
  Mapper &mapper = Mapper::instance();

  for (auto &clusterObject : Clusters) {
    DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
    float row = pos.getPadSecPos().getPadPos().getRow();
    float pad = pos.getPadSecPos().getPadPos().getPad();
    float cpad = GetCPad(clusterObject);
    if (row == 31 || row == 32 || row == 9 || row == 17 || row == 25 || row == 39 || row == 45 || row == 51) {
      histo->Fill(row, cpad, 1);
      continue;
    }
    if (runNr <= 243) {
      if (row >= 57 || row < 2) {
        histo->Fill(row, cpad, 1);
        continue;
      }

      int nexthit = HitMapS0->getValue(roc, row, pad+1);
      int previoushit = HitMapS0->getValue(roc, row, pad-1);
      if (nexthit != 1 || previoushit != 1) {
        histo->Fill(row, cpad, 1);
        continue;
      }
    }
    else if (runNr > 255) {
      if (row >= 61 || row <= 1) {
        histo->Fill(row, cpad, 1);
        continue;
      }
      else if (row > 1 && row < 61) {
        if (pad == mapper.getNumberOfPadsInRowSector(row)/2 - 1) {
          histo->Fill(row, cpad, 1);
          continue;
        }
      }
      int nexthit = HitMapS1->getValue(roc, row, pad+1);
      int previoushit = HitMapS1->getValue(roc, row, pad-1);
      if (nexthit != 1 || previoushit != 1) {
        histo->Fill(row, cpad, 1);
        continue;
      }
    }
    values.push_back((type == 0)?clusterObject.getQmax()/norm:clusterObject.getQ()/norm);
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

int testingWithParse(TString FileList, TString GainMap="", const char *OutputPath="/home/tom/myana/Results")
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

  TFile *OutFile = new TFile(Form("%s/testingWithParse.root", OutputPath), "recreate");
  TTree *TrackAna = new TTree("TrackAna","dEdxAna");

  TrackAna->Branch("PidEdxResQ", &pionres, "pionres/F");
  TrackAna->Branch("PidEdxResQmax", &pionresMax, "pionresMax/F");
  TrackAna->Branch("EledEdxResQ", &electronres, "electronres/F");
  TrackAna->Branch("EledEdxResQmax", &electronresMax, "electronresMax/F");
  TrackAna->Branch("SepPowerQ", &separationpower, "separationpower/F");
  TrackAna->Branch("SepPowerQmax", &separationpowerMax);
  TrackAna->Branch("posEle", &electronmeanTot, "electronmeanTot/F");
  TrackAna->Branch("sigmaEle", &electronsigmaTot, "electronsigmaTot/F");
  TrackAna->Branch("posPi", &pionmeanTot, "pionmeanTot/F");
  TrackAna->Branch("sigmaPi", &pionsigmaTot, "pionsigmaTot/F");
  TrackAna->Branch("posEleErr", &electronmeanerr, "electronmeanerr/F");
  TrackAna->Branch("sigmaEleErr", &electronsigmaerr, "electronsigmaerr/F");
  TrackAna->Branch("posPiErr", &pionmeanerr, "pionmeanerr/F");
  TrackAna->Branch("sigmaPiErr", &pionsigmaerr, "pionsigmaerr/F");
  TrackAna->Branch("posEleQmax", &electronmeanMax, "electronmeanMax/F");
  TrackAna->Branch("sigmaEleQmax", &electronsigmaMax, "electronsigmaMax/F");
  TrackAna->Branch("posPiQmax", &pionmeanMax, "pionmeanMax/F");
  TrackAna->Branch("sigmaPiQmax", &pionsigmaMax, "pionsigmaMax/F");
  TrackAna->Branch("posEleQmaxErr", &electronmeanerrMax, "electronmeanerrMax/F");
  TrackAna->Branch("sigmaEleQmaxErr", &electronsigmaerrMax, "electronsigmaerrMax/F");
  TrackAna->Branch("posPiQmaxErr", &pionmeanerrMax, "pionmeanerrMax/F");
  TrackAna->Branch("sigmaPiQmaxErr", &pionsigmaerrMax, "pionsigmaerrMax/F");
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

      file = FileList+"/"+file;

      //cout<<endl<<endl<<arr->At(ifile)->GetName()<<endl<<endl;
      fChain.Add(file);
    }
  }
  cout<<endl<<endl<<endl<<"Chain ready!!"<<endl<<endl;
  TObjArray *chainEntries = fChain.GetListOfFiles();
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    cout<<chainEntries->At(ifile)->GetTitle()<<endl;
  }

  int numberofruns = chainEntries->GetEntriesFast();



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


  //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
  TH1D *hNclustersUsed          = new TH1D("hNclustersUsed", ";Number of clusters per track; Counts", 300,0,300);

  TH1D *hNIROCtracks            = new TH1D("hNtracksIROC","; Number of tracks; Counts",20,0,20);
  TH1D *hNIROCtracksUsed        = new TH1D("hNtracksIROCUsed","; Number of tracks per event; Counts",50,0,50);

  TH1F *hdEdxEleTotAll          = new TH1F("hdEdxEleTotAll", "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
  TH1F *hdEdxPionTotAll         = new TH1F("hdEdxTotAll", "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);

  TH1F *hdEdxEleMax             = new TH1F("hdEdxEleMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 100, 0, 120);
  TH1F *hdEdxPionMax            = new TH1F("hdEdxMax", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 100, 0, 120);

  TH1D *hQ                      = new TH1D("hQ", "; Q_{tot} [ADC counts]; Counts", 600,0,600);
  TH1D *hQUsed                  = new TH1D("hQUsed", "; Q_{tot} [ADC counts]; Counts", 600,0,600);

  TH1D *hCherenkov              = new TH1D("hCherenkov", "; Cherenkov counter signal; Counts", 100,0,0.06);
  TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsed", "; Cherenkov counter signal; Counts", 100,0,0.06);
  TH1D *hCherenkovUsedEle       = new TH1D("hCherenkovUsedEle", "; Cherenkov counter signal; Counts", 100,0,0.06);

  TH2D *hPadOccupancy           = new TH2D("hOcc", "; Row; Pad; Counts", 63,0,63,100,0,100);
  TH2D *hPadOccupancyCPad       = new TH2D("hOccCPad", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsedCPad   = new TH2D("hOccUsed", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsed       = new TH2D("hOccUsedCPad", "; Row; Pad; Counts", 63,0,63,100,0,100);
  TH2D *hPadOccupancyUsedPions  = new TH2D("hOccUsedPions", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsedEle    = new TH2D("hOccUsedEle", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancySetting0   = new TH2D("hOccS0","; Row; Pad; Counts",63,0,63,100,-50,50);
  TH2D *hPadOccupancySetting1   = new TH2D("hOccS1","; Row; Pad; Counts",63,0,63,100,-50,50);

  TH2D *hQTimeMax               = new TH2D("hQTime","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);
  TH2D *hQTimeMaxUsed           = new TH2D("hQTimeUsed","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);

  TH2D *hQmaxRelPadPosUsed      = new TH2D("hQmaxRelPadPos", "; relative position of the cluster; Q_{max}; Counts", 100,-.2,1.2,600,0,1200);

  TH2D *hCherRuns               = new TH2D("hCherRuns", "; Run; Cherenkov value; Normalized counts",350,0,350,300,0,0.03);
  TH2D *hCherRunsUsed           = new TH2D("hCherRunsUsed", "; Run; Cherenkov value; Counts",350,0,350,300,0,0.03);

  TH2D *hClustersRuns           = new TH2D("hClustersRuns", ";Run; #clusters", 350,0,350,200,0,200);
  TH2D *hClustersRunsUsed       = new TH2D("hClustersRunsUsed", ";Run; #clusters", 350,0,350,200,0,200);

  TH2D *hNrTracksRuns           = new TH2D("hNrTracksRuns", "; Run; #tracks", 350,0,350,20,0,20);
  TH2D *hNrTracksRunsUsed       = new TH2D("hNrTracksRunsUsed", "; Run; #tracks", 350,0,350,20,0,20);

  TH2D *hTimeMeanRuns           = new TH2D("hTimeMeanRuns", "; Run; time mean", 350,0,350,50,0,50);
  TH2D *hTimeMeanRunsUsed       = new TH2D("hTimeMeanRunsUsed", "; Run; time mean", 350,0,350,50,0,50);

  TH2D *hdEdxPiRuns             = new TH2D("hdEdxPiRuns", "; Run; dE/dx max pions", 350,0,350,100,0,120);
  TH2D *hdEdxEleRuns            = new TH2D("hdEdxEleRuns", "; Run; dE/dx max electrons", 350,0,350,100,0,120);

  TH2D *hQtotPiRuns             = new TH2D("hQtotPiRuns", "; Run; Q_{tot}", 350,0,350,300,0,300);
  TH2D *hQtotEleRuns            = new TH2D("hQtotEleRuns", "; Run; Q_{tot}", 350,0,350,300,0,300);

//  TH2D *hQmaxPiRuns             = new TH2D("hQmaxPiRuns", "; Run; Q_{max}", 350,0,350,167,0,200);
//  TH2D *hQmaxEleRuns            = new TH2D("hQmaxEleRuns", "; Run; Q_{max}", 350,0,350,167,0,200);
  TH2D *hQmaxPiRuns             = new TH2D("hQmaxPiRuns", "; Run; Q_{max}; Normalized counts", 350,0,350,200,0,200);
  TH2D *hQmaxEleRuns            = new TH2D("hQmaxEleRuns", "; Run; Q_{max}; Normalized counts", 350,0,350,200,0,200);

  TH2D *hCheckExclude           = new TH2D("CheckExclude", "; Row; Pad", 63,0,63,100,-50,50);

  TH1D *hnClTrMeanTot           = new TH1D("NClTrMeanTot",";NClTrMean",70,0,70);
  TH1D *hnClTrMeanMax           = new TH1D("NClTrMeanMax",";NClTrMean",70,0,70);

  TH2D *hQmaxRowS0              = new TH2D("hQmaxRowS0","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS1              = new TH2D("hQmaxRowS1","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS0Corr          = new TH2D("hQmaxRowS0Corr","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS1Corr          = new TH2D("hQmaxRowS1Corr","; Row; Q_{max}; Counts",63,0,63,100,0,400);

  TH2D *hnTruncVsnClTot         = new TH2D("nTruncVsnClTot","; nCl; nClTruncUsedTot",70,0,70,70,0,70);
  TH2D *hnTruncVsnClMax         = new TH2D("nTruncVsnClMax","; nCl; nClTruncUsedMax",70,0,70,70,0,70);

  TH2D *hdEdxNCl2D              = new TH2D("hdEdxNCl2D","; nCl; dE/dx",63,0,63,800,0,200);

  TH1F *hClustersVsRowS0All     = new TH1F("hClustersVsRowS0All","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1All     = new TH1F("hClustersVsRowS1All","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS0UsedPiTracks     = new TH1F("hClustersVsRowS0UsedPiTracks","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1UsedPiTracks     = new TH1F("hClustersVsRowS1UsedPiTracks","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS0UsedEleTracks     = new TH1F("hClustersVsRowS0UsedEleTracks","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1UsedEleTracks     = new TH1F("hClustersVsRowS1UsedEleTracks","; Pad row; Clusters",63,0,63);


  TGraph *gPionMean = new TGraph();
  TGraph *gEleMean  = new TGraph();
  TGraph *gPionSigma = new TGraph();
  TGraph *gEleSigma = new TGraph();
  TGraph *gSeparationPower = new TGraph();
  TGraph *gPiRes = new TGraph();
  TGraph *gEleRes = new TGraph();

  TGraph *gNClPiRes = new TGraph();
  TGraph *gNClEleRes = new TGraph();
  TGraphErrors *gNClSeparation = new TGraphErrors();

  TH1F *hChargeDistTrack[10100];
  TH1F *hTruncTrack[10100];

  for (int iTrack = 0; iTrack <10100; ++iTrack) {
    hChargeDistTrack[iTrack] = new TH1F(Form("hChargeDistTrack_%i", iTrack),";Q_{tot} [ADC counts];Counts",1200,0,1200);
    hTruncTrack[iTrack] = new TH1F(Form("hTruncTrack_%i", iTrack),";Q_{tot} [ADC counts];Counts",1200,0,1200);
  }

  TH1F *hdEdxEleTot[350];
  TH1F *hdEdxPionTot[350];
  TH1F *hdEdxNClPi[150];
  TH1F *hdEdxNClEle[150];

  for (int iNCl = 0; iNCl < 150; ++iNCl) {
    hdEdxNClPi[iNCl] = new TH1F(Form("hdEdxNClPi_%i", iNCl),"; dEdx; Counts",500,0,250);
    hdEdxNClEle[iNCl] = new TH1F(Form("hdEdxNClEle_%i", iNCl),"; dEdx; Counts",500,0,250);
  }

//  TProfile *hQmaxRowS0Prof      = new TProfile();
//  TProfile *hQmaxRowS1Prof      = new TProfile();
//  TProfile *hQmaxRowS0CorrProf  = new TProfile();
//  TProfile *hQmaxRowS1CorrProf  = new TProfile();

/*  TH1F *hQmaxPiSingleRun[numberofruns];
  TH1F *hQmaxEleSingleRun[numberofruns];
  TH1F *hQtotPiSingleRun[numberofruns];
  TH1F *hQtotEleSingleRun[numberofruns];

  TF1 *MPVfitQmaxPiSingleRun[numberofruns];
  TF1 *MPVfitQmaxEleSingleRun[numberofruns];
  TF1 *MPVfitQtotPiSingleRun[numberofruns];
  TF1 *MPVfitQtotEleSingleRun[numberofruns];
*/

/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */


  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  // float pions = 0;
  // float electrons = 0;
  int usedcl = 0;
  int onetrcl = 0;
  // float pioncl = 0;
  // float elecl = 0;
  float CherenkovValue = 0.;
  int runNr = 0;
  int loopcounter = 0;
  int loopcounter2 = 0;
  bool setgainmapS1 = false;
  int nClTruncTot = 0;
  int nClTruncMax = 0;
  int ncl = 0;
  int ncledge = 0;
  int runcounter = 0;
  int eventcounter = 0;


  std::vector<TrackTPC> *vecEvent = 0;

  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;


/*=============================================== Loop over runs =========================================================*/
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
    TTree *tree = (TTree*)TreeFile->Get("events");
    //TTree *tree = (TTree*)chainEntries->At(ifile)->Get("events");

/*
    TString CurrentFileName = chainEntries->At(ifile)->GetTitle();
    TObjArray *SplitFileName = CurrentFileName.Tokenize("/");
    TString runTString = SplitFileName->At(5)->GetName();

    std::string runString = "";
    runString.append(runTString);
    char *run = new char[runString.length() + 1];
    std::strcpy(run, runString.c_str());

    hQmaxPiSingleRun[runcounter] = new TH1F(Form("hQmaxPiSingleRun_%s",run),"; Q_{max} Pi: Counts", 300,0,300);
    hQtotPiSingleRun[runcounter] = new TH1F(Form("hQtotPiSingleRun_%s",run),"; Q_{tot} Pi: Counts", 300,0,300);
    hQmaxEleSingleRun[runcounter] = new TH1F(Form("hQmaxEleSingleRun_%s",run),"; Q_{max} Ele: Counts", 300,0,300);
    hQtotEleSingleRun[runcounter] = new TH1F(Form("hQtotEleSingleRun_%s",run),"; Q_{tot} Ele: Counts", 300,0,300);
*/
    hdEdxEleTot[runcounter] = new TH1F(Form("hdEdxEleTot_%i", runcounter), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
    hdEdxPionTot[runcounter]= new TH1F(Form("hdEdxTot_%i", runcounter), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
    vecEvent=0;
    tree->SetBranchAddress("Tracks", &vecEvent);
    tree->SetBranchAddress("header", &Header);

    CherenkovValue = 0;
    runNr = 0;
    nclRun = 0;

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

    std::cout<<endl<<nclRun<<endl;
/*=============================================== Loop over events =========================================================*/

    for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
      tree->GetEntry(iev);

      int NTracks = vecEvent->size();
      hNrTracksRuns->Fill(runNr,NTracks, 1./nclRun);
/*========================================== CUT ================================================*/
      if (NTracks != TrPerEv) continue;																// only one-track events
      CherenkovValue = Header.cherenkovValue;
      runNr = Header.run;

      hCherRuns->Fill(runNr,CherenkovValue, 1./nclRun);


      Tracks += NTracks;

      hNIROCtracks->Fill(NTracks);

      nclEvent = 0;
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        nclTrack = clCont.size();
        nclEvent += nclTrack;
        //nclRun += nclTrack;
        hClustersRuns->Fill(runNr,nclEvent, 1./nclRun);
//        for (auto &clusterObject : clCont) {
//          float timeMean = clusterObject.getTimeMean();
//          hTimeMeanRuns->Fill(runNr,timeMean);
//        }
      }


      ++OneTrackEvents;
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
/*========================================== CUT ================================================*/
        ///if (track->GetROC() != 0) continue;															// only ROC 0


        ncl = clCont.size();
        bool isok = true;
        ncledge = 0;

        onetrcl += ncl;

        for (auto &clusterObject : clCont) {                     // loop over the clusters of all one-track-events
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
          float pad = pos.getPadSecPos().getPadPos().getPad();
          float cpad = GetCPad(clusterObject);
//        if (setting == 0) {
//          if (row <=2 || row >= 57) continue;
//        }
          //const int cpad = GetCPad(clusterObject);
          float timeMean = clusterObject.getTimeMean();
          double Q = clusterObject.getQ();
          double Qmax = clusterObject.getQmax();
          hTimeMeanRuns->Fill(runNr,timeMean, 1./nclRun);
          if (CherenkovValue < CherCutLow) {
            hQmaxPiRuns->Fill(runNr, Qmax, 1.);
            hQtotPiRuns->Fill(runNr, Q, 1.);
///            hQmaxPiSingleRun[runcounter]->Fill(Qmax);
///            hQtotPiSingleRun[runcounter]->Fill(Q);
          }
          else if (CherenkovValue > CherCutHigh) {
            hQmaxEleRuns->Fill(runNr, Qmax, 1.);
            hQtotEleRuns->Fill(runNr, Q, 1.);
///            hQmaxEleSingleRun[runcounter]->Fill(Qmax);
///            hQtotEleSingleRun[runcounter]->Fill(Q);
          }

          hQTimeMax->Fill(timeMean,Q);
          hQ->Fill(Q);
          hPadOccupancy->Fill(row,pad);
          hPadOccupancyCPad->Fill(row,cpad);
          if (runNr <= 243) {
            hClustersVsRowS0All->Fill(row);
            hPadOccupancySetting0->Fill(row,cpad);
            hQmaxRowS0->Fill(row,Qmax);
            //hQmaxRowS0Prof->Fill(row,Qmax);
            if (GainMap != "") {
              float QmaxCorr = Qmax/gainmapS0->getValue(roc,row,pad);
              hQmaxRowS0Corr->Fill(row,QmaxCorr);
              //hQmaxRowS0CorrProf->Fill(row,QmaxCorr);
            }
          }
          else if (runNr > 255) {
            hClustersVsRowS1All->Fill(row);
            hPadOccupancySetting1->Fill(row,cpad);
            hQmaxRowS1->Fill(row,Qmax);
            //hQmaxRowS1Prof->Fill(row,Qmax);
            if (GainMap != "") {
              float QmaxCorr = Qmax/gainmapS1->getValue(roc,row,pad);
              hQmaxRowS1Corr->Fill(row,QmaxCorr);
              //hQmaxRowS1CorrProf->Fill(row,QmaxCorr);
            }
          }
        }




/*========================================== CUT ================================================*/
        if (ncl < nclCut) continue;																	// cut on number of clusters per track

        hCherenkov->Fill(CherenkovValue);
/*========================================== CUT ================================================*/
        if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov
        if (CherenkovValue < CherCutLow) {
          hCherenkovUsedPions->Fill(CherenkovValue);
        }
        else if (CherenkovValue > CherCutHigh){
          hCherenkovUsedEle->Fill(CherenkovValue);
        }

      //double lambda = atan(trackObject.getTgl());
      //double phi = asin(trackObject.getSnp());

        for (auto &clusterObject : clCont) {															// make cuts on cluster properties
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
//          if (setting == 0) {
//            if (row < 2 || row >= 57) continue;
//          }
          float pad = pos.getPadSecPos().getPadPos().getPad();
          float timeMean = clusterObject.getTimeMean();
          float cpad = GetCPad(clusterObject);

        //cout<<endl<<"Row: "<<double(row)<<"\t"<<"Pad: "<<double(pad)<<"\t"<<"CPad: "<<double(cpad)<<endl;


/*========================================== CUT ================================================*/
          if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle

          if (!InsideEdge(row, cpad, setting, 3)){++ncledge;}
        }

/*========================================== CUT ================================================*/
        if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry


        if (isok == true){																			// track accepted
          usedcl += ncl;

          float dEdxTot = 0;
          float dEdxMax = 0;
          nClTruncTot = 0;
          nClTruncMax = 0;

          if (setting == 0) {
            dEdxTot = trackObject.getTruncatedMean(runNr, .0, .7, 1, false, false, true, hCheckExclude, 1, 2, &nClTruncTot, hTruncTrack[eventcounter], hChargeDistTrack[eventcounter]);
            dEdxMax = trackObject.getTruncatedMean(runNr, .0, .7, 0, false, false, true, hCheckExclude, 1, 2, &nClTruncMax, hTruncTrack[eventcounter], hChargeDistTrack[eventcounter]);
            if (nClTruncTot < 35) continue;
            hnClTrMeanTot->Fill(nClTruncTot);
            hnClTrMeanMax->Fill(nClTruncMax);
          }
          else if (setting == 1) {
            dEdxTot = trackObject.getTruncatedMean(runNr, .0, .7, 1, false, false, true, hCheckExclude, 1, 2, &nClTruncTot, hTruncTrack[eventcounter], hChargeDistTrack[eventcounter]);
            dEdxMax = trackObject.getTruncatedMean(runNr, .0, .7, 0, false, false, true, hCheckExclude, 1, 2, &nClTruncMax, hTruncTrack[eventcounter], hChargeDistTrack[eventcounter]);
            if (nClTruncTot < 35) continue;
            hnClTrMeanTot->Fill(nClTruncTot);
            hnClTrMeanMax->Fill(nClTruncMax);
          }
          hnTruncVsnClTot->Fill(ncl,nClTruncTot);
          hnTruncVsnClMax->Fill(ncl,nClTruncMax);


          if (CherenkovValue < CherCutLow){
            hdEdxPionTot[runcounter]->Fill(dEdxTot);
            hdEdxPionTotAll->Fill(dEdxTot);
            hdEdxNClPi[nClTruncTot]->Fill(dEdxTot);
            hdEdxNCl2D->Fill(nClTruncTot,dEdxTot);
            hdEdxPionMax->Fill(dEdxMax);
            hdEdxPiRuns->Fill(runNr,dEdxMax, 1./nclRun);
            for (auto &clusterObject : clCont) {															// loop over clusters
              DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
              float row = pos.getPadSecPos().getPadPos().getRow();
              if (setting == 0) {
                hClustersVsRowS0UsedPiTracks->Fill(row);
                if (row <=2 || row >= 57) continue;
              }
              else if (setting == 1) {
                hClustersVsRowS1UsedPiTracks->Fill(row);
                if (row <= 2 || row >= 61) continue;
              }
              float pad = pos.getPadSecPos().getPadPos().getPad();
              float cpad = GetCPad(clusterObject);
              float Qmax = clusterObject.getQmax();
              const float padmean = clusterObject.getPadMean();
              float relPadPos = padmean-pad;

              hPadOccupancyUsedPions->Fill(row,cpad);
              hQmaxRelPadPosUsed->Fill(relPadPos,Qmax);
            }
            ++pions;
            pioncl += ncl;
          }
          else if (CherenkovValue > CherCutHigh){
            hdEdxEleTot[runcounter]->Fill(dEdxTot);
            hdEdxEleTotAll->Fill(dEdxTot);
            hdEdxNClEle[nClTruncTot]->Fill(dEdxTot);
            hdEdxNCl2D->Fill(nClTruncTot,dEdxTot);
            hdEdxEleMax->Fill(dEdxMax);
            hdEdxEleRuns->Fill(runNr,dEdxMax, 1./nclRun);
            for (auto &clusterObject : clCont) {															// loop over clusters
              DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
              float row = pos.getPadSecPos().getPadPos().getRow();
              if (setting == 0) {
                hClustersVsRowS0UsedEleTracks->Fill(row);
                if (row <=2 || row >= 57) continue;
              }
              else if (setting == 1) {
                hClustersVsRowS1UsedEleTracks->Fill(row);
                if (row <= 2 || row >= 61) continue;
              }
              float pad = pos.getPadSecPos().getPadPos().getPad();
              float cpad = GetCPad(clusterObject);
              float Qmax = clusterObject.getQmax();
              const float padmean = clusterObject.getPadMean();
              float relPadPos = padmean-pad;

              hPadOccupancyUsedEle->Fill(row,cpad);
              hQmaxRelPadPosUsed->Fill(relPadPos,Qmax);
            }
            ++electrons;
            elecl += ncl;
          }
          ++usedTracks;
          hNIROCtracksUsed->Fill(NTracks);
          hNclustersUsed->Fill(ncl);
          for (auto &clusterObject : clCont) {															// loop over clusters
            DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
            float row = pos.getPadSecPos().getPadPos().getRow();
            if (setting == 0) {
              if (row <=2 || row >= 57) continue;
            }
            else if (setting == 1) {
              if (row <= 2 || row >= 61) continue;
            }
            float pad = pos.getPadSecPos().getPadPos().getPad();
            const int cpad = GetCPad(clusterObject);
            float timeMean = clusterObject.getTimeMean();
            double Q = clusterObject.getQ();

            hQTimeMaxUsed->Fill(timeMean,Q);
            hQUsed->Fill(Q);
            hPadOccupancyUsed->Fill(row,pad);
            hPadOccupancyUsedCPad->Fill(row,cpad);
            hCherRunsUsed->Fill(runNr,CherenkovValue, 1./nclRun);
            hClustersRunsUsed->Fill(runNr,ncl, 1./nclRun);
            hNrTracksRunsUsed->Fill(runNr,NTracks, 1./nclRun);
            hTimeMeanRunsUsed->Fill(runNr,timeMean, 1./nclRun);
          }

        } // end of 'if true'
      } // end of loop over tracks ///break; ///if (isok==true){break;} //use only one track
      ++eventcounter;
    } // end of loop over events
//    for (int q = 0; q < 300; ++q) {
//      float content = hQtotPiRuns->GetBinContent(runNr,q);
//      hQtotPiRuns->SetBinContent(runNr,q,content/nclRun);
//    }

     if (hdEdxPionTot[runcounter]->GetEntries() > 0 && hdEdxEleTot[runcounter]->GetEntries() > 0){
    //  TCanvas *c15 = new TCanvas("c15", "dEdx Qtot");
      TF1 *pionfit = new TF1("pionfit","gaus",hdEdxPionTot[runcounter]->GetXaxis()->GetXmin(),hdEdxPionTot[runcounter]->GetXaxis()->GetXmax());
      TF1 *electronfit = new TF1("electronfit","gaus",hdEdxEleTot[runcounter]->GetXaxis()->GetXmin(),hdEdxEleTot[runcounter]->GetXaxis()->GetXmax());
      hdEdxEleTot[runcounter]->SetLineColor(kRed);
      hdEdxPionTot[runcounter]->SetLineColor(kBlue);

      const Float_t frac=0.2;
      Int_t bin1=0,bin2=0;

      GetBinMinMax(hdEdxPionTot[runcounter],frac,bin1,bin2);
      hdEdxPionTot[runcounter]->Fit("pionfit","","",hdEdxPionTot[runcounter]->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTot[runcounter]->GetXaxis()->GetBinUpEdge(bin2));
      GetBinMinMax(hdEdxEleTot[runcounter],frac,bin1,bin2);
      //hdEdxEleTot[runcounter]->Fit("electronfit","","",hdEdxEleTot[runcounter]->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTot[runcounter]->GetXaxis()->GetBinUpEdge(bin2));

      //alternative fit
      //hdEdxPionTot[runcounter]->Fit("pionfit");
      hdEdxEleTot[runcounter]->Fit("electronfit");


      hdEdxEleTot[runcounter]->GetFunction("electronfit")->SetLineColor(kRed);
      hdEdxPionTot[runcounter]->GetFunction("pionfit")->SetLineColor(kBlue);
      hdEdxPionTot[runcounter]->Draw();
      hdEdxEleTot[runcounter]->Draw("same");
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


      float pionreserror = sqrt(pow((pionsigmaerr/pionmeanTot),2)+pow(((pionsigmaTot*pionmeanerr)/pow(pionmeanTot,2)),2));
      float electronreserror = sqrt(pow((electronsigmaerr/electronmeanTot),2)+pow(((electronsigmaTot*electronmeanerr)/pow(electronmeanTot,2)),2));

      float separationpowererr = sqrt(pow((2*electronmeanerr/(electronsigmaTot+pionsigmaTot)),2)+pow((2*pionmeanerr/(electronsigmaTot+pionsigmaTot)),2)+
            pow(((2*electronsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2)+pow(((2*pionsigmaerr*(pionmeanTot-electronmeanTot))/(pow((electronsigmaTot+pionsigmaTot),2))),2));

      float elechisquareTot = electronfit->GetChisquare();
      float eleNDFTot = electronfit->GetNDF();
      float pionchisquareTot = pionfit->GetChisquare();
      float pionNDFTot = pionfit->GetNDF();

    //  TPaveText *pave1=new TPaveText(0.6,.7,.9,.9,"NDC");
    //  pave1->SetBorderSize(1);
    //  pave1->SetFillColor(10);
    //  pave1->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanTot,electronsigmaTot, electronsigmaTot/electronmeanTot*100));
    //  pave1->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanTot,pionsigmaTot,pionsigmaTot/pionmeanTot*100));
    //  pave1->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanTot-pionmeanTot)/((electronsigmaTot+pionsigmaTot)/2.)));
    //  pave1->Draw("same");

    //  c15->Update();
    //  c15->Print(Form("%s/dEdxQtot.png", OutputPath));

    //  TCanvas *c16 = new TCanvas("c16", "dEdx Qmax");
      TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
      TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
      hdEdxEleMax->SetLineColor(kRed);
      hdEdxPionMax->SetLineColor(kBlue);

      GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
      hdEdxPionMax->Fit("pionfitMax","","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
      GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
      //hdEdxEleMax->Fit("electronfitMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

      hdEdxEleMax->Fit("electronfitMax");


      hdEdxEleMax->GetFunction("electronfitMax")->SetLineColor(kRed);
      hdEdxPionMax->GetFunction("pionfitMax")->SetLineColor(kBlue);
      hdEdxPionMax->Draw();
      hdEdxEleMax->Draw("same");
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


      float pionreserrorMax = sqrt(pow((pionsigmaerrMax/pionmeanMax),2)+pow(((pionsigmaMax*pionmeanerrMax)/pow(pionmeanMax,2)),2));
      float electronreserrorMax = sqrt(pow((electronsigmaerrMax/electronmeanMax),2)+pow(((electronsigmaMax*electronmeanerrMax)/(electronmeanMax,2)),2));

      float separationpowererrMax = sqrt(pow((2*electronmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+pow((2*pionmeanerrMax/(electronsigmaMax+pionsigmaMax)),2)+
            pow(((2*electronsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2)+pow(((2*pionsigmaerrMax*(pionmeanMax-electronmeanMax))/(pow((electronsigmaMax+pionsigmaMax),2))),2));

      float elechisquareMax = electronfitMax->GetChisquare();
      float eleNDFMax = electronfitMax->GetNDF();
      float pionchisquareMax = pionfitMax->GetChisquare();
      float pionNDFMax = pionfitMax->GetNDF();

    //  TPaveText *pave2=new TPaveText(0.6,.7,.9,.9,"NDC");
    //  pave2->SetBorderSize(1);
    //  pave2->SetFillColor(10);
    //  pave2->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanMax,electronsigmaMax, electronsigmaMax/electronmeanMax*100));
    //  pave2->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanMax,pionsigmaMax,pionsigmaMax/pionmeanMax*100));
    //  pave2->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanMax-pionmeanMax)/((electronsigmaMax+pionsigmaMax)/2.)));
    //  pave2->Draw("same");

    //  c16->Update();
    //  c16->Print(Form("%s/dEdxQmax.png", OutputPath));

    gPionMean->SetPoint(runcounter,runNr,pionmeanTot);
    gEleMean->SetPoint(runcounter,runNr,electronmeanTot);
    gPionSigma->SetPoint(runcounter,runNr,pionsigmaTot);
    gEleSigma->SetPoint(runcounter,runNr,electronsigmaTot);
    gSeparationPower->SetPoint(runcounter,runNr,separationpower);
    gPiRes->SetPoint(runcounter,runNr,pionres);
    gEleRes->SetPoint(runcounter,runNr,electronres);
     } // end of if getentries > 0
    delete TreeFile;
    cout<<"Tree deleted"<<endl;
    ++runcounter;

//    hdEdxEleTot[runcounter]->Reset();
//    hdEdxPionTot[runcounter]->Reset();


    TF1 *hdEdxNClPiFit[63];
    TF1 *hdEdxNClEleFit[63];



    for (int iNCl = 0; iNCl < 63; ++iNCl) {
      hdEdxNClPiFit[iNCl] = new TF1(Form("hdEdxNClPiFit[%i]", iNCl),"gaus",0,500);
      hdEdxNClEleFit[iNCl] = new TF1(Form("hdEdxNClEleFit[%i]", iNCl),"gaus",0,500);
      if (hdEdxNClPi[iNCl]->GetEntries() > 0) {
        hdEdxNClPi[iNCl]->Fit(Form("hdEdxNClPiFit[%i]",iNCl));
      }
      if (hdEdxNClEle[iNCl]->GetEntries() > 0) {
        hdEdxNClEle[iNCl]->Fit(Form("hdEdxNClEleFit[%i]",iNCl));
      }

      float pimean = hdEdxNClPiFit[iNCl]->GetParameter(1);
      float pisigma = hdEdxNClPiFit[iNCl]->GetParameter(2);
      float pires = pisigma/pimean;

      float elemean = hdEdxNClEleFit[iNCl]->GetParameter(1);
      float elesigma = hdEdxNClEleFit[iNCl]->GetParameter(2);
      float eleres = elesigma/elemean;

      float pionsigmaerrNCl = hdEdxNClPiFit[iNCl]->GetParError(2);
      float pionmeanerrNCl = hdEdxNClPiFit[iNCl]->GetParError(1);
      float electronsigmaerrNCl = hdEdxNClEleFit[iNCl]->GetParError(2);
      float electronmeanerrNCl = hdEdxNClEleFit[iNCl]->GetParError(1);

      float separationpowerNCl = 2*(elemean-pimean)/(pisigma+elesigma);

      float separationpowerNClerr = sqrt(pow((2*electronmeanerrNCl/(elesigma+pisigma)),2)+pow((2*pionmeanerrNCl/(elesigma+pisigma)),2)+
            pow(((2*electronsigmaerrNCl*(pimean-elemean))/(pow((elesigma+pisigma),2))),2)+pow(((2*pionsigmaerrNCl*(pimean-elemean))/(pow((elesigma+pisigma),2))),2));

      gNClPiRes->SetPoint(iNCl,iNCl,pires);
      gNClEleRes->SetPoint(iNCl,iNCl,eleres);
      gNClSeparation->SetPoint(iNCl,iNCl,separationpowerNCl);
      gNClSeparation->SetPointError(iNCl,0,separationpowerNClerr);
    }
/*
    MPVfitQmaxEleSingleRun[runcounter] = new TF1(Form("MPVfitQmaxEleSingleRun_%s", run),"landau");
    MPVfitQtotEleSingleRun[runcounter] = new TF1(Form("MPVfitQtotEleSingleRun_%s", run),"landau");
    MPVfitQmaxPiSingleRun[runcounter] = new TF1(Form("MPVfitQmaxPiSingleRun_%s", run),"landau");
    MPVfitQtotPiSingleRun[runcounter] = new TF1(Form("MPVfitQtotPiSingleRun_%s", run),"landau");

    hQmaxEleSingleRun[runcounter]->Fit(Form("MPVfitQmaxEleSingleRun_%s", run),"NQ0");
    float QmaxEleMPV = MPVfitQmaxEleSingleRun[runcounter]->GetParameter(1);
    hQtotEleSingleRun[runcounter]->Fit(Form("MPVfitQtotEleSingleRun_%s", run),"NQ0");
    float QtotEleMPV = MPVfitQtotEleSingleRun[runcounter]->GetParameter(1);
    hQmaxPiSingleRun[runcounter]->Fit(Form("MPVfitQmaxPiSingleRun_%s", run),"NQ0");
    float QmaxPiMPV = MPVfitQmaxPiSingleRun[runcounter]->GetParameter(1);
    hQtotPiSingleRun[runcounter]->Fit(Form("MPVfitQtotPiSingleRun_%s", run),"NQ0");
    float QtotPiMPV = MPVfitQtotPiSingleRun[runcounter]->GetParameter(1);

    std::cout << "\n\n\n\n" << "Electrons:" << "\n\n\n\n\n";
    std::cout << "run number: " << runNr << "\t Qmax: " << QmaxEleMPV << "\t Qtot: " <<  QtotEleMPV;
    std::cout << "\n\n" << "Pions:" << "\n\n";
    std::cout << "run number: " << runNr << "\t Qmax: " << QmaxPiMPV << "\t Qtot: " <<  QtotPiMPV;
*/



  } // end of loop over runs


    TCanvas *c15 = new TCanvas("c15", "dEdx Qtot");
    TF1 *pionfitEnd = new TF1("pionfitEnd","gaus",hdEdxPionTotAll->GetXaxis()->GetXmin(),hdEdxPionTotAll->GetXaxis()->GetXmax());
    TF1 *electronfitEnd = new TF1("electronfitEnd","gaus",hdEdxEleTotAll->GetXaxis()->GetXmin(),hdEdxEleTotAll->GetXaxis()->GetXmax());
    hdEdxEleTotAll->SetLineColor(kRed);
    hdEdxPionTotAll->SetLineColor(kBlue);

    const Float_t frac=0.2;
    Int_t bin1=0,bin2=0;

    GetBinMinMax(hdEdxPionTotAll,frac,bin1,bin2);
    hdEdxPionTotAll->Fit("pionfitEnd","","",hdEdxPionTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTotAll->GetXaxis()->GetBinUpEdge(bin2));
    GetBinMinMax(hdEdxEleTotAll,frac,bin1,bin2);
    //hdEdxEleTotAll->Fit("electronfitEnd","","",hdEdxEleTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTotAll->GetXaxis()->GetBinUpEdge(bin2));

    //alternative fit
    //hdEdxPionTotAll->Fit("pionfitEnd");
    hdEdxEleTotAll->Fit("electronfitEnd");


    hdEdxEleTotAll->GetFunction("electronfitEnd")->SetLineColor(kRed);
    hdEdxPionTotAll->GetFunction("pionfitEnd")->SetLineColor(kBlue);
    hdEdxPionTotAll->Draw();
    hdEdxEleTotAll->Draw("same");
    pionmeanTot = pionfitEnd->GetParameter(1);
    pionsigmaTot = pionfitEnd->GetParameter(2);
    electronmeanTot = electronfitEnd->GetParameter(1);
    electronsigmaTot = electronfitEnd->GetParameter(2);

    pionres = pionsigmaTot/pionmeanTot;
    electronres = electronsigmaTot/electronmeanTot;
    pionsigmaerr = pionfitEnd->GetParError(2);
    pionmeanerr = pionfitEnd->GetParError(1);
    electronsigmaerr = electronfitEnd->GetParError(2);
    electronmeanerr = electronfitEnd->GetParError(1);


    separationpower = 2*(electronmeanTot-pionmeanTot)/(pionsigmaTot+electronsigmaTot);

    TPaveText *pave1=new TPaveText(0.6,.7,.9,.9,"NDC");
    pave1->SetBorderSize(1);
    pave1->SetFillColor(10);
    pave1->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanTot,electronsigmaTot, electronsigmaTot/electronmeanTot*100));
    pave1->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanTot,pionsigmaTot,pionsigmaTot/pionmeanTot*100));
    pave1->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanTot-pionmeanTot)/((electronsigmaTot+pionsigmaTot)/2.)));
    pave1->Draw("same");

    c15->Update();
    c15->Print(Form("%s/dEdxQtot.png", OutputPath));

    TCanvas *c16 = new TCanvas("c16", "dEdx Qmax");
    TF1 *pionfitMaxEnd = new TF1("pionfitMaxEnd","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
    TF1 *electronfitEndMax = new TF1("electronfitEndMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
    hdEdxEleMax->SetLineColor(kRed);
    hdEdxPionMax->SetLineColor(kBlue);

    GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
    hdEdxPionMax->Fit("pionfitMaxEnd","","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
    GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
    //hdEdxEleMax->Fit("electronfitEndMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

    hdEdxEleMax->Fit("electronfitEndMax");


    hdEdxEleMax->GetFunction("electronfitEndMax")->SetLineColor(kRed);
    hdEdxPionMax->GetFunction("pionfitMaxEnd")->SetLineColor(kBlue);
    hdEdxPionMax->Draw();
    hdEdxEleMax->Draw("same");
    pionmeanMax = pionfitMaxEnd->GetParameter(1);
    pionsigmaMax = pionfitMaxEnd->GetParameter(2);
    electronmeanMax = electronfitEndMax->GetParameter(1);
    electronsigmaMax = electronfitEndMax->GetParameter(2);

    pionresMax = pionsigmaMax/pionmeanMax;
    electronresMax = electronsigmaMax/electronmeanMax;
    pionsigmaerrMax = pionfitMaxEnd->GetParError(2);
    pionmeanerrMax = pionfitMaxEnd->GetParError(1);
    electronsigmaerrMax = electronfitEndMax->GetParError(2);
    electronmeanerrMax = electronfitEndMax->GetParError(1);


    separationpowerMax = 2*(electronmeanMax-pionmeanMax)/(pionsigmaMax+electronsigmaMax);

    TPaveText *pave2=new TPaveText(0.6,.7,.9,.9,"NDC");
    pave2->SetBorderSize(1);
    pave2->SetFillColor(10);
    pave2->AddText(Form("e: %.2f #pm %.2f (%.2f%%)",electronmeanMax,electronsigmaMax, electronsigmaMax/electronmeanMax*100));
    pave2->AddText(Form("#pi: %.2f #pm %.2f (%.2f%%)",pionmeanMax,pionsigmaMax,pionsigmaMax/pionmeanMax*100));
    pave2->AddText(Form("Separation: %.2f#sigma", TMath::Abs(electronmeanMax-pionmeanMax)/((electronsigmaMax+pionsigmaMax)/2.)));
    pave2->Draw("same");

    c16->Update();
    c16->Print(Form("%s/dEdxQmax.png", OutputPath));

  TCanvas *mean = new TCanvas();
  gPionMean->GetYaxis()->SetRangeUser(0,400);
  gPionMean->GetYaxis()->SetTitle("dE/dx (a.u.)");
  gPionMean->GetXaxis()->SetTitle("Run");
  gPionMean->SetMarkerStyle(kPlus);
  gPionMean->SetMarkerSize(2);
  gPionMean->SetMarkerColor(kBlue);
  gEleMean->SetMarkerStyle(kPlus);
  gEleMean->SetMarkerSize(2);
  gEleMean->SetMarkerColor(kRed);
  gPionMean->Draw("AP");
  gEleMean->Draw("P,same");
  mean->Print(Form("%s/dEdxMeanRun.png", OutputPath));

  TCanvas *sigma = new TCanvas();
  gPionSigma->GetYaxis()->SetRangeUser(0,80);
  gPionSigma->GetYaxis()->SetTitle("Width dE/dx (a.u.)");
  gPionSigma->GetXaxis()->SetTitle("Run");
  gPionSigma->SetMarkerStyle(kPlus);
  gPionSigma->SetMarkerSize(2);
  gPionSigma->SetMarkerColor(kBlue);
  gEleSigma->SetMarkerStyle(kPlus);
  gEleSigma->SetMarkerSize(2);
  gEleSigma->SetMarkerColor(kRed);
  gPionSigma->Draw("AP");
  gEleSigma->Draw("P,same");
  sigma->Print(Form("%s/dEdxSigmaRun.png", OutputPath));

  TCanvas *separation = new TCanvas();
  gSeparationPower->GetYaxis()->SetRangeUser(0,5);
  gSeparationPower->GetYaxis()->SetTitle("Separation power [#sigma]");
  gSeparationPower->GetXaxis()->SetTitle("Run");
  gSeparationPower->SetMarkerStyle(kPlus);
  gSeparationPower->SetMarkerSize(2);
  gSeparationPower->SetMarkerColor(kBlue+2);
  gSeparationPower->Draw("AP");
  separation->Print(Form("%s/SeparationRun.png", OutputPath));

  TCanvas *piores = new TCanvas();
  gPiRes->GetYaxis()->SetRangeUser(0,0.5);
  gPiRes->GetYaxis()->SetTitle("dE/dx resolution");
  gPiRes->GetXaxis()->SetTitle("Run");
  gPiRes->SetMarkerStyle(kPlus);
  gPiRes->SetMarkerSize(2);
  gPiRes->SetMarkerColor(kBlue);
  gEleRes->SetMarkerStyle(kPlus);
  gEleRes->SetMarkerSize(2);
  gEleRes->SetMarkerColor(kRed);
  gPiRes->Draw("AP");
  gEleRes->Draw("P,same");
  piores->Print(Form("%s/dEdxResolution.png", OutputPath));

  TCanvas *pioresNCl = new TCanvas();
  gNClPiRes->GetYaxis()->SetRangeUser(0,0.5);
  gNClPiRes->GetYaxis()->SetTitle("dE/dx resolution NCl");
  gNClPiRes->GetXaxis()->SetTitle("NCl");
  gNClPiRes->SetMarkerStyle(kPlus);
  gNClPiRes->SetMarkerSize(2);
  gNClPiRes->SetMarkerColor(kBlue);
  gNClEleRes->SetMarkerStyle(kPlus);
  gNClEleRes->SetMarkerSize(2);
  gNClEleRes->SetMarkerColor(kRed);
  gNClPiRes->Draw("AP");
  gNClEleRes->Draw("P,same");
  pioresNCl->Print(Form("%s/dEdxResolutionNCl.png", OutputPath));

  TCanvas *separationNCl = new TCanvas();
  gNClSeparation->GetYaxis()->SetRangeUser(0,5);
  gNClSeparation->GetYaxis()->SetTitle("Separation power [#sigma]");
  gNClSeparation->GetXaxis()->SetTitle("NCl");
  gNClSeparation->SetMarkerStyle(kPlus);
  gNClSeparation->SetMarkerSize(2);
  gNClSeparation->SetMarkerColor(kBlue+2);
  gNClSeparation->Draw("AP");
  separationNCl->Print(Form("%s/SeparationRunNCl.png", OutputPath));


  for (int irun = 0; irun < 350; ++irun) {
    float normpi = hQmaxPiRuns->Integral(irun, irun, 1, 400);
    float normele = hQmaxEleRuns->Integral(irun, irun, 1, 400);
    float normcher = hCherRuns->Integral(irun, irun, 1, 300);
    if (normpi > 0) {
      for (float biny = 1; biny <= 400; ++biny) {
        float newvaluepi = hQmaxPiRuns->GetBinContent(irun,biny)/normpi;
        hQmaxPiRuns->SetBinContent(irun,biny,newvaluepi);
      }
    }
    if (normele > 0) {
      for (float biny = 1; biny <= 400; ++biny) {
        float newvalueele = hQmaxEleRuns->GetBinContent(irun,biny)/normele;
        hQmaxEleRuns->SetBinContent(irun,biny,newvalueele);
      }
    }
    if (normcher > 0) {
      for (float biny = 1; biny <= 400; ++biny) {
        float newvaluecher = hCherRuns->GetBinContent(irun,biny)/normcher;
        hCherRuns->SetBinContent(irun,biny,newvaluecher);
      }
    }
  }

//  TProfile *PQmaxRowS0     = (TProfile*)hQmaxRowS0->ProfileX();
//  TProfile *PQmaxRowS0Corr = (TProfile*)hQmaxRowS0Corr->ProfileX();
//  TProfile *PQmaxRowS1     = (TProfile*)hQmaxRowS1->ProfileX();
//  TProfile *PQmaxRowS1Corr = (TProfile*)hQmaxRowS1Corr->ProfileX();

//  TCanvas *CanQmaxRowS0 = new TCanvas("1","QmaxRowS0");
//  hQmaxRowS0->Draw("colz");
//  PQmaxRowS0->Draw("same");
//  CanQmaxRowS0->Update();
//  CanQmaxRowS0->Print(Form("%s/QmaxRowS0.png", OutputPath));
//  TCanvas *CanQmaxRowS0Corr = new TCanvas("2","QmaxRowS0Corr");
//  hQmaxRowS0Corr->Draw("colz");
//  PQmaxRowS0Corr->Draw("same");
//  CanQmaxRowS0Corr->Update();
//  CanQmaxRowS0Corr->Print(Form("%s/QmaxRowS0Corr.png", OutputPath));
//  TCanvas *CanQmaxRowS1 = new TCanvas("3","QmaxRowS1");
//  hQmaxRowS1->Draw("colz");
//  PQmaxRowS1->Draw("same");
//  CanQmaxRowS1->Update();
//  CanQmaxRowS1->Print(Form("%s/QmaxRowS1.png", OutputPath));
//  TCanvas *CanQmaxRowS1Corr = new TCanvas("4","QmaxRowS1Corr");
//  hQmaxRowS1Corr->Draw("colz");
//  PQmaxRowS1Corr->Draw("same");
//  CanQmaxRowS1Corr->Update();
//  CanQmaxRowS1Corr->Print(Form("%s/QmaxRowS1Corr.png", OutputPath));





  /* ======================================================================
   * ===================== Plot everything ================================
   * ====================================================================== */

//  hQ->SetLineColor(kBlue+2);

//  TCanvas *c1 = new TCanvas("c1", "Number of tracks: used events");
//  hNIROCtracksUsed->Draw();

//  TCanvas *c2 = new TCanvas("c2", "Number of tracks: all events");
//  hNIROCtracks->Draw();

//  TCanvas *c3 = new TCanvas("c3", "Number of clusters of used tracks");
//  hNclustersUsed->Draw();

//  TCanvas *c4 = new TCanvas("c4", "Qtot vs. TimeMean: all tracks");
//  ///hQTimeMax->GetZaxis()->SetRangeUser(0,820);
//  hQTimeMax->Draw("colz");

//  TCanvas *c5 = new TCanvas("c5", "Qtot vs. TimeMean: used tracks");
//  ///hQTimeMaxUsed->GetZaxis()->SetRangeUser(0,820);
//  hQTimeMaxUsed->Draw("colz");

//  TCanvas *c6 = new TCanvas("c6", "Cluster charge distribution: all tracks");
//  hQ->Draw();

//  TCanvas *c7 = new TCanvas("c7", "Cluster charge distribution: used tracks");
//  hQUsed->Draw();

//  TCanvas *c8 = new TCanvas("c8", "Number of hits per pad: all particles, all tracks");
//  ///hPadOccupancy->GetZaxis()->SetRangeUser(0,3700);
//  hPadOccupancy->Draw("colz");

//  TCanvas *c9 = new TCanvas("c9", "Number of hits per pad: all particles, used tracks");
//  ///hPadOccupancyUsed->GetZaxis()->SetRangeUser(0,3700);
//  hPadOccupancyUsed->Draw("colz");

//  TCanvas *OccUsedPions = new TCanvas("OccUsedPions", "Number of hits per pad: used pions");
//  hPadOccupancyUsedPions->Draw("colz");

//  TCanvas *OccUsedEle = new TCanvas("OccUsedEle", "Number of hits per pad: used electrons");
//  hPadOccupancyUsedEle->Draw("colz");

//  TCanvas *OccCPad = new TCanvas("OccCPad", "Number of hits per pad: all particles, all tracks, CPad coordinates");
//  hPadOccupancyCPad->Draw("colz");

//  TCanvas *OccUsedCPad = new TCanvas("OccUsedCPad","Number of hits per pad: all particles, used tracks, CPad coordinates");
//  hPadOccupancyUsedCPad->Draw("colz");

//  TCanvas *c10 = new TCanvas("c10", "Cherenkov value: all tracks");
//  c10->SetLogy();
  hCherenkov->SetLineColor(kBlue);
//  hCherenkov->SetLineWidth(2);
//  hCherenkov->Draw();

//  TCanvas *c11 = new TCanvas("c11", "Cherenkov value: used tracks");
//  c11->SetLogy();
//  hCherenkovUsedPions->GetYaxis()->SetRangeUser(0.5,21000);
  hCherenkovUsedPions->SetLineColor(kBlue);
  hCherenkovUsedEle->SetLineColor(kRed);
//  hCherenkovUsedPions->SetLineWidth(2);
//  hCherenkovUsedEle->SetLineWidth(2);

//  hCherenkovUsedPions->Draw();
//  hCherenkovUsedEle->Draw("same");

//  TCanvas *cCherRun = new TCanvas("cCherRun","Cherenkov values per run: all tracks");
//  hCherRuns->Draw("colz");

//  TCanvas *cCherRunUsed = new TCanvas("cCherRunUsed","Cherenkov values per run: used tracks");
//  hCherRunsUsed->Draw("colz");

//  TCanvas *cClustersRuns = new TCanvas("cClustersRuns", "Number of clusters: all events, all tracks");
//  hClustersRuns->Draw("colz");

//  TCanvas *cClustersRunsUsed = new TCanvas("cClustersRunsUsed", "Number of clusters: events with only one track, used tracks");
//  hClustersRunsUsed->Draw("colz");

//  TCanvas *cNTracksRuns = new TCanvas("cNTracksRuns", "Number of tracks per event: all events");
//  hNrTracksRuns->Draw("colz");

//  TCanvas *cNTracksRunsUsed = new TCanvas("cNTracksRunsUsed", "Number of tracks per event: used events");
//  hNrTracksRunsUsed->Draw("colz");

//  TCanvas *cTimeMean = new TCanvas("cTimeMean", "time mean: all events, all tracks, all clusters");
//  hTimeMeanRuns->Draw("colz");

//  TCanvas *cTimeMeanUsed = new TCanvas("cTimeMeanUsed", "time mean: used clusters");
//  hTimeMeanRunsUsed->Draw("colz");

//  TCanvas *cdEdxMaxPiRuns = new TCanvas("cdEdxMaxPiRuns", "dE/dx max pions, all runs");
//  hdEdxPiRuns->Draw("colz");

//  TCanvas *cdEdxMaxEleRuns = new TCanvas("cdEdxMaxEleRuns", "dE/dx max electrons, all runs");
//  hdEdxEleRuns->Draw("colz");



  TrackAna->Fill();
  OutFile->Write();
  OutFile->Close();

//  TCanvas *QmaxRowS1 = new TCanvas();
//  hQmaxRowS1->Draw("colz");
//  hQmaxRowS1Prof->Draw("same");

//  TFile *f = TFile::Open("QA_plots.root","recreate");
//  f->WriteObject(hNclustersUsed, "NClustersTrack");
//  f->WriteObject(hNIROCtracks, "NTracksEvent");
//  f->WriteObject(hNIROCtracksUsed, "NTracksEventUsed");
//  f->WriteObject(hQ, "ClusterChargeDist");
//  f->WriteObject(hQUsed, "ClusterChargeDistUsed");
//  f->WriteObject(hCherenkov, "CherenkovDist");
//  f->WriteObject(hPadOccupancy, "PadHits");
//  f->WriteObject(hPadOccupancyCPad, "PadHitsCPad");
//  f->WriteObject(hPadOccupancyUsed, "PadHitsUsed");
//  f->WriteObject(hPadOccupancyUsedCPad, "PadHitsUsedCpad");
//  f->WriteObject(hPadOccupancyUsedPions, "PadHitsPiUsedCpad");
//  f->WriteObject(hPadOccupancyUsedEle, "PadHitsEleUsedCPad");
//  f->WriteObject(hQTimeMax, "QTime");
//  f->WriteObject(hQTimeMaxUsed, "QTimeUsed");
//  f->WriteObject(hCherRuns, "CherenkovRuns");
//  f->WriteObject(hCherRunsUsed, "CherenkovRunsUsed");
//  f->WriteObject(hClustersRuns, "NClustersRuns");
//  f->WriteObject(hClustersRunsUsed, "NClustersRunsUsed");
//  f->WriteObject(hNrTracksRuns, "NTracksRuns");
//  f->WriteObject(hNrTracksRunsUsed, "NTracksRunsUsed");
//  f->WriteObject(hTimeMeanRuns, "TimeMeanRuns");
//  f->WriteObject(hTimeMeanRunsUsed, "TimeMeanRunsUsed");
//  f->WriteObject(hdEdxPiRuns, "dEdxPiRuns");
//  f->WriteObject(hdEdxEleRuns, "dEdxEleRuns");
//  delete f;

  cout<<endl;cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl;
  cout<<"============================= Total Charge ======================================="<<endl;
//  cout<<"Pion dE/dx resolution:          ("<<pionres*100<<" +- "<<pionreserror*100<<") %"<<endl;
//  cout<<"Electron dE/dx resolution:      ("<<electronres*100<<" +- "<<electronreserror*100<<") %"<<endl;
//  cout<<"Separation power:               ("<<separationpower<<" +- "<<separationpowererr<<") sigma"<<endl;
//  cout<<"Pion Chisquare:                  "<<pionchisquareTot/pionNDFTot<<endl;
//  cout<<"Electron Chisquare:              "<<elechisquareTot/eleNDFTot<<endl;
//  cout<<"Pion mean:                       "<<pionmeanTot<<" +- "<<pionsigmaTot<<endl;
//  cout<<"Electron mean:                   "<<electronmeanTot<<" +- "<<electronsigmaTot<<endl;
//  cout<<endl;
//  cout<<"============================== Max Charge ========================================"<<endl;
//  cout<<"Pion dE/dx resolution:          ("<<pionresMax*100<<" +- "<<pionreserrorMax*100<<") %"<<endl;
//  cout<<"Electron dE/dx resolution:      ("<<electronresMax*100<<" +- "<<electronreserrorMax*100<<") %"<<endl;
//  cout<<"Separation power:               ("<<separationpowerMax<<" +- "<<separationpowererrMax<<") sigma"<<endl;
//  cout<<"Pion Chisquare:                  "<<pionchisquareMax/pionNDFMax<<endl;
//  cout<<"Electron Chisquare:              "<<elechisquareMax/eleNDFMax<<endl;
//  cout<<"Pion mean:                       "<<pionmeanMax<<" +- "<<pionsigmaMax<<endl;
//  cout<<"Electron mean:                   "<<electronmeanMax<<" +- "<<electronsigmaMax<<endl;
//  cout<<endl;
//  cout<<"============================== Statistics ========================================"<<endl;
//  cout<<"Number of tracks:                "<<Tracks<<endl;
//  cout<<"Number of one-track-events:      "<<OneTrackEvents<<endl;
//  cout<<"Number of used tracks:           "<<usedTracks<<endl;
//  cout<<"Number of pions:                 "<<pions<<endl;
//  cout<<"Number of electrons:             "<<electrons<<endl;
//  cout<<endl;
//  /*cout<<"Number of one-track-clusters:    "<<onetrcl<<endl;
//  cout<<"Number of used clusters:         "<<usedcl<<endl;
//  cout<<"Number of electron clusters:     "<<elecl<<endl;
//  cout<<"Number of pion clusters:         "<<pioncl<<endl;*/
//  cout<<"Average clusters per pion:       "<<pioncl/pions<<endl;
//  cout<<"Average clusters per electron:   "<<elecl/electrons<<endl;
//  cout<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<"=================================================================================="<<endl;
  cout<<endl<<endl;

  //delete TreeFile;
  //if (fGain1) delete fGain1;
  //if (fGain2) delete fGain2;
  /*delete hNIROCtracks;
  delete hNIROCtracksUsed;
  delete hdEdxEleTot;
  delete hdEdxPionTot;
  delete hdEdxEleMax;
  delete hdEdxPionMax;
  delete hQ;
  delete hQUsed;
  delete hCherenkov;
  delete hCherenkovUsedPions;
  delete hCherenkovUsedEle;
  delete hPadOccupancy;
  delete hPadOccupancyUsed;
  delete hChisquare;
  delete hGainMap;
  delete hQTimeMax;
  delete hQTimeMaxUsed;
  delete c1;
  delete c2;
  delete c3;
  delete c4;
  delete c5;
  delete c6;
  delete c7;
  delete c8;
  delete c9;
  delete c10;
  delete c11;
  delete c15;
  delete c16;
  delete pionfit;
  delete electronfit; */
  return 0;
}


