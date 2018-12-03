#include "DataFormatsTPC/Cluster.h"
#include "TPCBase/Mapper.h"
#include "DataFormatsTPC/TrackTPC.h"
#include "TPCBase/CRU.h"
#include "TPCBase/ROC.h"
#include "TPCBase/CalDet.h"
#include "TPCBase/Painter.h"
#include "TPCSimulation/SAMPAProcessing.h"


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
#include "TFile.h"
#include "TPCBase/CDBInterface.h"

#include <math.h>
#include <boost/lambda/lambda.hpp>
#include <iostream>
#include <cstring>
#include <stdlib.h>



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
/*
float getTrackOffset(std::vector<TrackTPC> track)
{
  float x = track.getX();
  float y = track.getY();
  float p = track.getSnp();

  return y - tan(asin(p))*x;
}

float getTrackSlope(std::vector<TrackTPC> track)
{
  return tan(asin(track.getSnp()));
}
*/
float evalX(const TrackTPC &track, float pos, TH1F *TrackAngleLocal=nullptr, TH1F *TrackAngleGlobal=nullptr, TH1F *Slope=nullptr, TH1F *Offset=nullptr) // x-pos [cm]
{
  float x = track.getX();
  float y = track.getY();
  float p = asin(track.getSnp());
  float pGlob = p+10*TMath::DegToRad();
  float pDeg = p*TMath::RadToDeg();
  float pGlobDeg = pGlob*TMath::RadToDeg();

  float xGlob = (x*cos(350*TMath::DegToRad()) + y*sin(350*TMath::DegToRad()));
  float yGlob = (-x*sin(350*TMath::DegToRad()) + y*cos(350*TMath::DegToRad()));

  float slope = tan(pGlob);
  float offset = yGlob - slope * xGlob;

  if (TrackAngleLocal!=nullptr) {
    TrackAngleLocal->Fill(pDeg);
  }
  if (TrackAngleGlobal!=nullptr) {
    TrackAngleGlobal->Fill(pGlobDeg);
  }
  if (Slope!=nullptr) {
    Slope->Fill(slope);
  }
  if (Offset!=nullptr) {
    Offset->Fill(offset);
  }

  return offset + slope * pos;
}

int identify_particle(TString inputString)
{
  const char *pion = "p";
  const char *electron = "e";

  TObjArray *split = inputString.Tokenize("/");

  if (strncmp(split->At(split->GetEntriesFast()-1)->GetName(), pion, 1) == 0) {
    for (int i = 0; i<split->GetEntriesFast(); ++i) {
      //std::cout<<std::endl<<split->At(i)->GetName()<<std::endl;
    }
    //std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<split->At(split->GetEntriesFast()-1)->GetName()<<std::endl<<std::endl<<std::endl<<std::endl;
    return 0;
  }

  else if (strncmp(split->At(split->GetEntriesFast()-1)->GetName(), electron, 1) == 0) {
    for (int i = 0; i<split->GetEntriesFast(); ++i) {
      //std::cout<<std::endl<<split->At(i)->GetName()<<std::endl;
    }
    //std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<split->At(split->GetEntriesFast()-1)->GetName()<<std::endl<<std::endl<<std::endl<<std::endl;
    return 1;
  }

  else{
    return 2;
  }
}

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


int compareSimvsReco(TString FileList, TString DataType="sim", TString GainMap="", const char *Filename="compareSimvsReco.root", const char *OutputPath="/scratch2/tklemenz/SimAna/Results/newO2/compareSimvsReco/tracks")
{
  //gROOT->ProcessLine(".x /home/gu74yub/rootlogon.C");
  //gStyle->SetOptStat(0);
  auto& cdb = CDBInterface::instance();
  cdb.setUseDefaults();
  EventHeader Header;
  Mapper &mapper = Mapper::instance();
  SAMPAProcessing &sampa = SAMPAProcessing::instance();


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

  TFile *OutFile = new TFile(Form("%s/%s", OutputPath,Filename), "recreate");
/*  TTree *TrackAna = new TTree("TrackAna","dEdxAna");

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
*/




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


  if (DataType == "sim") {
	  
	  std::vector<TrackTPC> *fChainEvent=0;

	  TChain fChain("o2sim");

	  fChain.SetBranchAddress("TPCTracks",&fChainEvent);
	  //fChain.SetBranchAddress("header",&Header);

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


  /* ============================================================
   * ============= Define histograms and graphs =================
   * ============================================================*/

  TH1D *hNclustersUsed          = new TH1D("hNclustersUsedSim", ";Number of clusters per track; Counts", 300,0,300);

  TH1D *hNIROCtracks            = new TH1D("hNtracksIROCSim","; Number of tracks; Counts",20,0,20);
  TH1D *hNIROCtracksUsed        = new TH1D("hNtracksIROCUsedSim","; Number of tracks per event; Counts",50,0,50);

  TH1F *hdEdxEleTotAll          = new TH1F("hdEdxEleTotAllSim", "; d#it{E}/d#it{x} (a.u.); Counts", 2000, 0, 4000);
  TH1F *hdEdxPionTotAll         = new TH1F("hdEdxTotAllSim", "; d#it{E}/d#it{x} (a.u.); Counts", 2000, 0, 4000);

  TH1F *hdEdxEleMax             = new TH1F("hdEdxEleMaxSim", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 1000, 0, 1200);
  TH1F *hdEdxPionMax            = new TH1F("hdEdxMaxSim", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 1000, 0, 1200);

  TH1D *hQ                      = new TH1D("hQtotSim", "; Q_{tot} [ADC counts]; Counts", 600,0,600);
  TH1D *hQmax                   = new TH1D("hQmaxSim", "; Q_{max} [ADC counts]; Counts", 600,0,600);
  TH1D *hQUsed                  = new TH1D("hQtotUsedSim", "; Q_{tot} [ADC counts]; Counts", 600,0,600);
  TH1D *hQmaxUsed               = new TH1D("hQmaxUsedSim", "; Q_{max} [ADC counts]; Counts", 600,0,600);

  TH1D *hCRU                    = new TH1D("hCRUSim","; CRU; Counts", 10,0,10);
  TH1D *hRow                    = new TH1D("hRowSim","; Row; Counts",75,-5,70);
  TH1D *hPadMean                = new TH1D("hPadMeanSim","; PadMean; Counts",800,20,100);
  TH1D *hPadSigma               = new TH1D("hPadSigmaSim","; PadSigma; Counts",400,-2,2);
  TH1D *hTimeMean               = new TH1D("hTimeMeanSim","; TimeMean; Counts",200,0,20);
  TH1D *hTimeSigma              = new TH1D("hTimeSigmaSim","; TimeSigma; Counts",800,-4,4);

  TH2D *hQtotVsTimeMean         = new TH2D("hQtotVsTimeMeanSim","; TimeMean; Qtot",100,0,20,800,0,800);
  TH2D *hQmaxVsTimeMean         = new TH2D("hQmaxVsTimeMeanSim","; TimeMean; Qmax",100,0,20,800,0,800);

  TH2D *hQmaxVsQtot             = new TH2D("hQmaxVsQtotSim","; Qmax [ADC]; Qtot [ADC]",600,0,600,600,0,600);
  TH2D *hQNeighboringRowLeft    = new TH2D("hQNeighboringRowLeftSim","; Qtot [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
  TH2D *hQNeighboringRowRight   = new TH2D("hQNeighboringRowRightSim","; Qtot [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
  TH2D *hQmaxNeighboringRowLeft  = new TH2D("hQmaxNeighboringRowLeftSim","; Qmax [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
  TH2D *hQmaxNeighboringRowRight = new TH2D("hQmaxNeighboringRowRightSim","; Qmax [ADC]; Q_{neighbor}", 600,0,600,600,0,600);

  TH1D *hTrackResiduals         = new TH1D("hTrackResidualsSim", "; Y [cm]; Counts", 150,-1.5,1.5);
  TH1D *hTrackResidualsTime     = new TH1D("hTrackResidualsTimeReco","; Track residuals t-dir; Counts", 1000,-50,50);
  TH1D *hEvalX                  = new TH1D("hEvalXSim", "; Y [cm]; Counts", 150,0,15);
  TH1D *hyClus                  = new TH1D("hyClusSim", "; Y [cm]; Counts", 200,0,50);
  TH1D *hxClus                  = new TH1D("hxClusSim", "; X [cm]; Counts", 100,80,130);
  TH1D *hCpadMean               = new TH1D("hCpadMeanSim", "; cPad; Counts", 120,0,30);
  TH1D *hTgl                    = new TH1D("hLambdaTimeDirSim","; #lambda [deg]; Counts",100,-2,2);

  TH1F *hTrackAngleLocal        = new TH1F("hTrackAngleLocalSim","; local angle; Counts",100,0,10);
  TH1F *hTrackAngleGlobal       = new TH1F("hTrackAngleGlobalSim","; global angle; Counts",200,0,20);
  TH1F *hSlope                  = new TH1F("hSlopeSim","; slope; Counts",200,-0.4,0.4);
  TH1F *hOffset                 = new TH1F("hOffsetSim","; offset; Counts",100,-10,10);

  TH2D *hTrackResVsxCluster     = new TH2D("hTrackResVsxClusterSim","; x_{cluster} [cm]; Track residual [cm]", 125,80,130,100,-1,1);
  TH2D *hTrackResVsSigCluster   = new TH2D("hTrackResVsSigClusterSim","; #sigma_{cluster} [pad?]; Track residual [cm]", 100,0,2,100,-1,1);
  TH2D *hTrackResVsRelPadPos    = new TH2D("hTrackResVsRelPadPosSim","; relative pad position; Track residual [cm]", 100,0,1,100,-1,1);
  TH2D *hTrackResVsTrackAngle   = new TH2D("hTrackResVsTrackAngleSim","; Track angle [deg]; Track residual [cm]", 500,0,5,100,-1,1);
  TH2D *hSigClusterVsxCluster   = new TH2D("hSigClusterVsxClusterSim","; x_{cluster} [cm]; #sigma_{cluster} [pad?]", 125,80,130,100,0,2);
  TH2D *hSigClusterVsyCentre    = new TH2D("hSigClusterVsyCentreSim","; y_{pad center} [cm]; #sigma_{cluster} [pad?]", 100,15,35,100,-0.5,2);
  TH2D *hSigClusterVsyCluster   = new TH2D("hSigClusterVsyClusterSim","; y_{cluster} [cm]; #sigma_{cluster} [pad?]", 100,15,35,100,-0.5,2);
  TH2D *hSigClusterVsRelPadPos  = new TH2D("hSigClusterVsRelPadPosSim","; relative pad position; #sigma_{cluster} [pad?]", 100,0,1,100,-0.5,2);

  TH1F *hTrackResOnePadClus     = new TH1F("hTrackResOnePadClusterSim",";Track residual (only one pad clusters); Counts",150,-1.5,1.5);
  TH1F *hTrackResNoOnePadClus   = new TH1F("hTrackResNoOnePadClusterSim",";Track residual (no one pad clusters); Counts",150,-1.5,1.5);

  TH1F *hETRes                  = new TH1F("hETResSim",";Track residual;Counts;",100,-1,1);
  TH1F *hETPadMean              = new TH1F("hETPadMeanSim",";pad mean;Counts;",1000,0,100);
  TH1F *hETPad                  = new TH1F("hETPadSim",";pad;Counts;",100,0,100);
  TH1F *hETRelPad               = new TH1F("hETRelPadPositionSim",";Relative pad position;Counts;",100,0,1);

  TH1F *hRelPadPos              = new TH1F("hRelPadPosSim",";relative pad position;Counts",400,-2,2);

  TH2D *hPadOccupancy           = new TH2D("hOccSim", "; Row; Pad; Counts", 63,0,63,100,0,100);
  TH2D *hPadOccupancyCPad       = new TH2D("hOccCPadSim", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsedCPad   = new TH2D("hOccUsedSim", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsed       = new TH2D("hOccUsedCPadSim", "; Row; Pad; Counts", 63,0,63,100,0,100);
  TH2D *hPadOccupancyUsedPions  = new TH2D("hOccUsedPionsSim", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancyUsedEle    = new TH2D("hOccUsedEleSim", "; Row; Pad; Counts", 63,0,63,100,-50,50);
  TH2D *hPadOccupancySetting0   = new TH2D("hOccS0Sim","; Row; Pad; Counts",63,0,63,100,-50,50);
  TH2D *hPadOccupancySetting1   = new TH2D("hOccS1Sim","; Row; Pad; Counts",63,0,63,100,-50,50);

  TH2D *hQTimeMax               = new TH2D("hQTimeSim","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);
  TH2D *hQTimeMaxUsed           = new TH2D("hQTimeUsedSim","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);

  TH2D *hQmaxRelPadPosUsed      = new TH2D("hQmaxRelPadPosSim", "; relative position of the cluster; Q_{max}; Counts", 100,-.2,1.2,600,0,1200);

  TH2D *hQmaxRowS0              = new TH2D("hQmaxRowS0Sim","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS1              = new TH2D("hQmaxRowS1Sim","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS0Corr          = new TH2D("hQmaxRowS0CorrSim","; Row; Q_{max}; Counts",63,0,63,100,0,400);
  TH2D *hQmaxRowS1Corr          = new TH2D("hQmaxRowS1CorrSim","; Row; Q_{max}; Counts",63,0,63,100,0,400);

  TH1F *hClustersVsRowS0All     = new TH1F("hClustersVsRowS0AllSim","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1All     = new TH1F("hClustersVsRowS1AllSim","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS0UsedPiTracks     = new TH1F("hClustersVsRowS0UsedPiTracksSim","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1UsedPiTracks     = new TH1F("hClustersVsRowS1UsedPiTracksSim","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS0UsedEleTracks     = new TH1F("hClustersVsRowS0UsedEleTracksSim","; Pad row; Clusters",63,0,63);
  TH1F *hClustersVsRowS1UsedEleTracks     = new TH1F("hClustersVsRowS1UsedEleTracksSim","; Pad row; Clusters",63,0,63);

  TGraph *gTrack = new TGraph();
  TGraph *gCluster = new TGraph();

  TGraph *gETRes = new TGraph();
  TGraph *gETPadMean = new TGraph();
  TGraph *gETPad = new TGraph();
  TGraph *gETRelPad = new TGraph();

  TH1F *hdEdxEleTot[350];
  TH1F *hdEdxPionTot[350];

/* =============================================================================
 * ================== Loop over events and apply all cuts ======================
 * ======================= meanwhile fill histograms ===========================
 * ============================================================================= */


  int usedTracks = 0;
  int Tracks = 0;
  int OneTrackEvents = 0;
  int usedcl = 0;
  int onetrcl = 0;
  int runNr = 0;
  int loopcounter = 0;
  int loopcounter2 = 0;
  bool setgainmapS1 = false;
  int ncl = 0;
  int ncledge = 0;
  int runcounter = 0;
  int eventcounter = 0;
  int firsttrack = 0;

  bool pion = false;
  bool electron = false;


  std::vector<TrackTPC> *vecEvent = 0;

  cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;


/*=============================================== Loop over runs =========================================================*/
  for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
    TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
    cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;
    TTree *tree = (TTree*)TreeFile->Get("o2sim");
    //TTree *tree = (TTree*)chainEntries->At(ifile)->Get("events");


    pion = false;
    electron = false;

    if (identify_particle(chainEntries->At(ifile)->GetTitle()) == 0) {     /// GetName() ist falsch, sollte GetTitle sein!!!
      pion = true;
     // std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<"pions"<<std::endl<<std::endl<<std::endl<<std::endl;
    }
    else if (identify_particle(chainEntries->At(ifile)->GetTitle()) == 1) {
      electron = true;
      //std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<"electrons"<<std::endl<<std::endl<<std::endl<<std::endl;
    }
    else {
      //std::cout<<std::endl<<std::endl<<std::endl<<std::endl<<"ERROR: Particles could not be identified"<<std::endl<<std::endl<<std::endl<<std::endl;
    }





    hdEdxEleTot[runcounter] = new TH1F(Form("hdEdxEleTot_%i", runcounter), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
    hdEdxPionTot[runcounter]= new TH1F(Form("hdEdxTot_%i", runcounter), "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
    vecEvent=0;
    nclRun = 0;
    tree->SetBranchAddress("TPCTracks", &vecEvent);


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
/*========================================== CUT ================================================*/
      if (NTracks != TrPerEv) continue;																// only one-track events



      Tracks += NTracks;

      hNIROCtracks->Fill(NTracks);

      nclEvent = 0;
      for (auto trackObject : *vecEvent) {
        std::vector<Cluster> clCont;
        trackObject.getClusterVector(clCont);
        nclTrack = clCont.size();
        nclEvent += nclTrack;
        //nclRun += nclTrack;
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
          //trackObject.setGainMap(GainMap, setting);
          std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
        }
        else if (loopcounter2 == 1 && usegainmap) {
          //trackObject.setGainMap(GainMap, setting);
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
          float timeMean = clusterObject.getTimeMean();
          float Q = clusterObject.getQ();
          float Qmax = clusterObject.getQmax();
/*          if (pion && !electron) {
            hQmaxPiRuns->Fill(runNr, Qmax, 1.);
            hQtotPiRuns->Fill(runNr, Q, 1.);
          }
          else if (!pion && electron) {
            hQmaxEleRuns->Fill(runNr, Qmax, 1.);
            hQtotEleRuns->Fill(runNr, Q, 1.);
          }
*/
          hQTimeMax->Fill(timeMean,Q);
          hQ->Fill(Q);
          hQmax->Fill(Qmax);
          hPadOccupancy->Fill(row,pad);
          hPadOccupancyCPad->Fill(row,cpad);
          if (runNr <= 243) {
            hClustersVsRowS0All->Fill(row);
            hPadOccupancySetting0->Fill(row,cpad);
            hQmaxRowS0->Fill(row,Qmax);
            if (GainMap != "") {
              float QmaxCorr = Qmax/gainmapS0->getValue(roc,row,pad);
              hQmaxRowS0Corr->Fill(row,QmaxCorr);
            }
          }
          else if (runNr > 255) {
            hClustersVsRowS1All->Fill(row);
            hPadOccupancySetting1->Fill(row,cpad);
            hQmaxRowS1->Fill(row,Qmax);
            if (GainMap != "") {
              float QmaxCorr = Qmax/gainmapS1->getValue(roc,row,pad);
              hQmaxRowS1Corr->Fill(row,QmaxCorr);
            }
          }
        }




/*========================================== CUT ================================================*/
        if (ncl < nclCut) continue;																	// cut on number of clusters per track

/*========================================== CUT ================================================*/

        for (auto &clusterObject : clCont) {															// make cuts on cluster properties
          DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
          float row = pos.getPadSecPos().getPadPos().getRow();
//          if (setting == 0) {
//            if (row < 2 || row >= 57) continue;
//          }
          float pad = pos.getPadSecPos().getPadPos().getPad();
          float timeMean = clusterObject.getTimeMean();
          float cpad = GetCPad(clusterObject);


/*========================================== CUT ================================================*/
//          if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle

//          if (!InsideEdge(row, cpad, setting, 3)){++ncledge;}
        }

/*========================================== CUT ================================================*/
        if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry


        if (isok == true){																			// track accepted
          usedcl += ncl;

          float dEdxTot = 0;
          float dEdxMax = 0;

          float angle = asin(trackObject.getSnp())*TMath::RadToDeg();

          if (setting == 0) {
            dEdxTot = trackObject.getTruncatedMean(0,0.7,1);
            dEdxMax = trackObject.getTruncatedMean(0,0.7,0);
          }
          else if (setting == 1) {
            dEdxTot = trackObject.getTruncatedMean(0,0.7,1);
            dEdxMax = trackObject.getTruncatedMean(0,0.7,0);
          }

          if (pion && !electron){
            hdEdxPionTot[runcounter]->Fill(dEdxTot);
            hdEdxPionTotAll->Fill(dEdxTot);
            hdEdxPionMax->Fill(dEdxMax);
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
          else if (!pion && electron){
            hdEdxEleTot[runcounter]->Fill(dEdxTot);
            hdEdxEleTotAll->Fill(dEdxTot);
            hdEdxEleMax->Fill(dEdxMax);
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

          int clustercounter = 0;
          for (auto &clusterObject : clCont) {            // loop over clusters
            const CRU cru(clusterObject.getCRU());
            const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
            const int rowInSector = clusterObject.getRow(); // + region.getGlobalRowOffset();
            const GlobalPadNumber pad = mapper.globalPadNumber(PadPos(rowInSector, clusterObject.getPadMean()));
            const PadCentre& padCentre = mapper.padCentre(pad);
            const float localYfactor = (cru.side() == Side::A) ? -1.f : 1.f;
            float zPosition = sampa.getZfromTimeBin(clusterObject.getTimeMean(), cru.side());

            LocalPosition3D posLoc(padCentre.X(), localYfactor * padCentre.Y(), zPosition);
            GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

            DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
            float row = pos.getPadSecPos().getPadPos().getRow();
            if (setting == 0) {
              if (row <= 2 || row >= 57) continue;
            }
            else if (setting == 1) {
              if (row <= 2 || row >= 61) continue;
            }
            //float pad = pos.getPadSecPos().getPadPos().getPad();
            const int cpad = GetCPad(clusterObject);
            float timeMean = clusterObject.getTimeMean();
            float Q = clusterObject.getQ();
            float Qmax = clusterObject.getQmax();
            float padMean = clusterObject.getPadMean();
            float padSigma = clusterObject.getPadSigma();
            float padLoc = pos.getPadSecPos().getPadPos().getPad();

            int clusCRU = clusterObject.getCRU();
            float timeMeanSigma = clusterObject.getTimeSigma();

            for (auto &clusterObject2 : clCont) {
              DigitPos pos(clusterObject2.getCRU(), PadPos(clusterObject2.getRow(), clusterObject2.getPadMean()));
              float row2 = pos.getPadSecPos().getPadPos().getRow();
              if (setting == 0) {
                if (row2-1 <=2 || row2+1 >= 57) continue;
              }
              else if (setting == 1) {
                if (row2-1 <= 2 || row2+1 >= 61) continue;
              }
              if (row2 == row-1) {
                hQNeighboringRowLeft->Fill(Q,clusterObject2.getQ());
                hQmaxNeighboringRowLeft->Fill(Qmax,clusterObject2.getQmax());
              }
              else if (row2 == row+1) {
                hQNeighboringRowRight->Fill(Q,clusterObject2.getQ());
                hQmaxNeighboringRowRight->Fill(Qmax,clusterObject2.getQmax());
              }
            }

            float relPadPos = padMean - padLoc;

            float xClus = posGlob.X();
            float yClus = posGlob.Y();//-0.2;//+relPadPos*0.4;
            float yCentreGlob = posGlob.Y();

            float CpadMean = cpad + relPadPos;

            float yTrack = evalX(trackObject, xClus);

            float residual = yClus - yTrack;

            hTrackResiduals->Fill(residual);

            hCRU->Fill(clusCRU);
            hRow->Fill(row);
            hTimeMean->Fill(timeMean);
            hTimeSigma->Fill(timeMeanSigma);
            hPadMean->Fill(padMean);
            hPadSigma->Fill(padSigma);
            hQtotVsTimeMean->Fill(timeMean,Q);
            hQmaxVsTimeMean->Fill(timeMean,Qmax);

            hTrackResVsxCluster->Fill(xClus,residual);
            hTrackResVsSigCluster->Fill(padSigma,residual);
            hTrackResVsRelPadPos->Fill(relPadPos,residual);
            hTrackResVsTrackAngle->Fill(angle,residual);
            hSigClusterVsxCluster->Fill(xClus,padSigma);
            hSigClusterVsyCentre->Fill(yCentreGlob,padSigma);
            hSigClusterVsyCluster->Fill(yClus,padSigma);
            hSigClusterVsRelPadPos->Fill(relPadPos,padSigma);

            hEvalX->Fill(evalX(trackObject, xClus));
            hCpadMean->Fill(CpadMean);
            hxClus->Fill(xClus);
            hyClus->Fill(yClus);

            hTgl->Fill(atan(trackObject.getTgl())*TMath::RadToDeg());

            hRelPadPos->Fill(relPadPos);

            hQTimeMaxUsed->Fill(timeMean,Q);
            hQUsed->Fill(Q);
            hQmaxUsed->Fill(Qmax);
            hPadOccupancyUsed->Fill(row,pad);
            hPadOccupancyUsedCPad->Fill(row,cpad);
            hQmaxVsQtot->Fill(Qmax,Q);

            if (padSigma == 0) {
              hTrackResOnePadClus->Fill(residual);
            }
            else {
              hTrackResNoOnePadClus->Fill(residual);
            }

            if (Qmax > Q) {
              hQmaxVsQtot->Fill(Qmax,Q,100000000);
            }

            if (firsttrack == 1) {
              gCluster->SetPoint(clustercounter,xClus,yClus);
              gETRes->SetPoint(clustercounter,xClus,residual);
              gETPadMean->SetPoint(clustercounter,xClus,padMean);
              gETPad->SetPoint(clustercounter,xClus,padLoc);
              gETRelPad->SetPoint(clustercounter,xClus,relPadPos);
              hETRes->Fill(residual);
              hETPadMean->Fill(padMean);
              hETPad->Fill(padLoc);
              hETRelPad->Fill(relPadPos);
            }

            ++clustercounter;

          }

          if (firsttrack == 1) {
            for (int step = 0; step < 500; ++step) {
              gTrack->SetPoint(step,step*0.375,evalX(trackObject, step*0.375));
            }
          }

          evalX(trackObject, 0, hTrackAngleLocal, hTrackAngleGlobal, hSlope, hOffset);
          ++firsttrack;
        } // end of 'if true'

      } // end of loop over tracks ///break; ///if (isok==true){break;} //use only one track
      ++eventcounter;
    } // end of loop over events

     if (hdEdxPionTot[runcounter]->GetEntries() > 0 && hdEdxEleTot[runcounter]->GetEntries() > 0){
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



     } // end of if getentries > 0
    delete TreeFile;
    cout<<"Tree deleted"<<endl;
    ++runcounter;




  } // end of loop over runs

    TF1 *pionfitEnd = new TF1("pionfitEnd","gaus",hdEdxPionTotAll->GetXaxis()->GetXmin(),hdEdxPionTotAll->GetXaxis()->GetXmax());
    TF1 *electronfitEnd = new TF1("electronfitEnd","gaus",hdEdxEleTotAll->GetXaxis()->GetXmin(),hdEdxEleTotAll->GetXaxis()->GetXmax());
    hdEdxEleTotAll->SetLineColor(kRed);
    hdEdxPionTotAll->SetLineColor(kBlue);

    const Float_t frac=0.2;
    Int_t bin1=0,bin2=0;

    GetBinMinMax(hdEdxPionTotAll,frac,bin1,bin2);
    hdEdxPionTotAll->Fit("pionfitEnd","0","",hdEdxPionTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTotAll->GetXaxis()->GetBinUpEdge(bin2));
    GetBinMinMax(hdEdxEleTotAll,frac,bin1,bin2);
    //hdEdxEleTotAll->Fit("electronfitEnd","","",hdEdxEleTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTotAll->GetXaxis()->GetBinUpEdge(bin2));

    //alternative fit
    //hdEdxPionTotAll->Fit("pionfitEnd");
    hdEdxEleTotAll->Fit("electronfitEnd","0");


    electronfitEnd->SetLineColor(kRed);
    pionfitEnd->SetLineColor(kBlue);
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


    TF1 *pionfitMaxEnd = new TF1("pionfitMaxEnd","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
    TF1 *electronfitEndMax = new TF1("electronfitEndMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
    hdEdxEleMax->SetLineColor(kRed);
    hdEdxPionMax->SetLineColor(kBlue);

    GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
    hdEdxPionMax->Fit("pionfitMaxEnd","","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
    GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
    //hdEdxEleMax->Fit("electronfitEndMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

    hdEdxEleMax->Fit("electronfitEndMax");


    electronfitEndMax->SetLineColor(kRed);
    pionfitMaxEnd->SetLineColor(kBlue);
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

    hQmaxUsed->Scale(1/hQmaxUsed->GetEntries());
    hQUsed->Scale(1/hQUsed->GetEntries());
    hTrackResiduals->Scale(1/hTrackResiduals->GetEntries());
    hdEdxPionTotAll->Scale(1/hdEdxPionTotAll->GetEntries());
    hdEdxEleTotAll->Scale(1/hdEdxEleTotAll->GetEntries());
    hdEdxPionMax->Scale(1/hdEdxPionMax->GetEntries());
    hdEdxEleMax->Scale(1/hdEdxEleMax->GetEntries());
    hTrackResOnePadClus->Scale(1/hTrackResiduals->GetEntries());
    hTrackResNoOnePadClus->Scale(1/hTrackResiduals->GetEntries());
    hCRU->Scale(1/hCRU->GetEntries());
    hRow->Scale(1/hRow->GetEntries());
    hPadMean->Scale(1/hPadMean->GetEntries());
    hPadSigma->Scale(1/hPadSigma->GetEntries());
    hTimeMean->Scale(1/hTimeMean->GetEntries());
    hTimeSigma->Scale(1/hTimeSigma->GetEntries());

    hTrackResiduals->SetLineColor(kCyan+2);
    hTrackResOnePadClus->SetLineColor(kMagenta+2);
    hTrackResNoOnePadClus->SetLineColor(kYellow+2);

    gTrack->SetMarkerStyle(20);
    gTrack->SetMarkerColor(kGreen+2);
    gCluster->SetMarkerStyle(20);
    gCluster->SetMarkerColor(kBlue+2);
    gTrack->GetXaxis()->SetTitle("x [cm]");
    gTrack->GetYaxis()->SetTitle("y [cm]");
    gCluster->GetXaxis()->SetTitle("x [cm]");
    gCluster->GetYaxis()->SetTitle("y [cm]");

    gETRes->SetMarkerStyle(20);
    gETPadMean->SetMarkerStyle(20);
    gETPad->SetMarkerStyle(20);
    gETRelPad->SetMarkerStyle(20);
    gETRes->SetMarkerColor(kRed+2);
    gETPadMean->SetMarkerColor(kCyan+2);
    gETPad->SetMarkerColor(kMagenta+2);
    gETRelPad->SetMarkerColor(kYellow+2);
    gETRes->GetXaxis()->SetTitle("x [cm]");
    gETPadMean->GetXaxis()->SetTitle("x [cm]");
    gETPad->GetXaxis()->SetTitle("x [cm]");
    gETRelPad->GetXaxis()->SetTitle("x [cm]");
    gETRes->GetYaxis()->SetTitle("Track residual");
    gETPadMean->GetYaxis()->SetTitle("Pad mean");
    gETPad->GetYaxis()->SetTitle("Pad");
    gETRelPad->GetYaxis()->SetTitle("Relative pad position");

    //TProfile *pQNeighboringRowRight   = (TProfile*)hQNeighboringRowRight->ProfileX();

    hQNeighboringRowLeft->Add(hQNeighboringRowRight);
    TProfile *pQNeighboringRowSim    = (TProfile*)hQNeighboringRowLeft->ProfileX();
    hQmaxNeighboringRowLeft->Add(hQmaxNeighboringRowRight);
    TProfile *pQmaxNeighboringRowSim    = (TProfile*)hQmaxNeighboringRowLeft->ProfileX();



  //TrackAna->Fill();
  OutFile->Write();
  OutFile->WriteObject(gTrack, "gTrack");
  OutFile->WriteObject(gCluster, "gCluster");
  OutFile->WriteObject(gETRes, "gETResidual");
  OutFile->WriteObject(gETPadMean, "gETPadMean");
  OutFile->WriteObject(gETPad, "gETPad");
  OutFile->WriteObject(gETRelPad, "gETRelPadPosition");
  OutFile->WriteObject(pQNeighboringRowSim, "pQNeighboringRowSim");
  OutFile->WriteObject(pQmaxNeighboringRowSim, "pQmaxNeighboringRowSim");
  //OutFile->WriteObject(pQNeighboringRowRight, "pQNeighboringRowRightSim");
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


  return 0;
  }
































































 else if (DataType == "reco") {
	 
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
	  
	  
    //* ============================================================
    //* ============= Define histograms and graphs =================
    //* ============================================================

    //TH1D *hNclusters   = new TH1D("hNclusters", ";Number of clusters; Counts", 65,0,65);
    TH1D *hNclustersUsed          = new TH1D("hNclustersUsedReco", ";Number of clusters per track; Counts", 300,0,300);

    TH1D *hNIROCtracks            = new TH1D("hNtracksIROCReco","; Number of tracks; Counts",20,0,20);
    TH1D *hNIROCtracksUsed        = new TH1D("hNtracksIROCUsedReco","; Number of tracks per event; Counts",50,0,50);

    TH1F *hdEdxEleTotAll          = new TH1F("hdEdxEleTotAllReco", "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);
    TH1F *hdEdxPionTotAll         = new TH1F("hdEdxTotAllReco", "; d#it{E}/d#it{x} (a.u.); Counts", 200, 0, 400);

    TH1F *hdEdxEleMax             = new TH1F("hdEdxEleMaxReco", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 100, 0, 120);
    TH1F *hdEdxPionMax            = new TH1F("hdEdxMaxReco", "; d#it{E}/d#it{x} Q_{max} (a.u.); Counts", 100, 0, 120);

    TH1D *hQ                      = new TH1D("hQtotReco", "; Q_{tot} [ADC counts]; Counts", 600,0,600);
    TH1D *hQUsed                  = new TH1D("hQtotUsedReco", "; Q_{tot} [ADC counts]; Counts", 600,0,600);
    TH1D *hQmax                   = new TH1D("hQmaxReco", "; Q_{max} [ADC counts]; Counts", 600,0,600);
    TH1D *hQmaxUsed               = new TH1D("hQmaxUsedReco", "; Q_{max} [ADC counts]; Counts", 600,0,600);

    TH1D *hCRU                    = new TH1D("hCRUReco","; CRU; Counts", 10,0,10);
    TH1D *hRow                    = new TH1D("hRowReco","; Row; Counts",75,-5,70);
    TH1D *hPadMean                = new TH1D("hPadMeanReco","; PadMean; Counts",800,20,100);
    TH1D *hPadSigma               = new TH1D("hPadSigmaReco","; PadSigma; Counts",400,-2,2);
    TH1D *hTimeMean               = new TH1D("hTimeMeanReco","; TimeMean; Counts",200,0,20);
    TH1D *hTimeSigma              = new TH1D("hTimeSigmaReco","; TimeSigma; Counts",800,-4,4);

    TH2D *hQtotVsTimeMean         = new TH2D("hQtotVsTimeMeanReco","; TimeMean; Qtot",100,0,20,800,0,800);
    TH2D *hQmaxVsTimeMean         = new TH2D("hQmaxVsTimeMeanReco","; TimeMean; Qmax",100,0,20,800,0,800);

    TH2D *hQmaxVsQtot             = new TH2D("hQmaxVsQtotReco","; Qmax [ADC]; Qtot [ADC]",600,0,600,600,0,600);
    TH2D *hQNeighboringRowLeft    = new TH2D("hQNeighboringRowLeftReco","; Qtot [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
    TH2D *hQNeighboringRowRight   = new TH2D("hQNeighboringRowRightReco","; Qtot [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
    TH2D *hQmaxNeighboringRowLeft  = new TH2D("hQmaxNeighboringRowLeftReco","; Qmax [ADC]; Q_{neighbor}", 600,0,600,600,0,600);
    TH2D *hQmaxNeighboringRowRight = new TH2D("hQmaxNeighboringRowRightReco","; Qmax [ADC]; Q_{neighbor}", 600,0,600,600,0,600);

    TH1D *hTrackResiduals         = new TH1D("hTrackResidualsReco", "; Y [cm]; Counts", 150,-1.5,1.5);
    TH1D *hTrackResidualsTime     = new TH1D("hTrackResidualsTimeReco","; Track residuals t-dir; Counts", 1000,-50,50);
    TH1D *hEvalX                  = new TH1D("hEvalXReco", "; Y [cm]; Counts", 150,0,15);
    TH1D *hyClus                  = new TH1D("hyClusReco", "; Y [cm]; Counts", 150,0,15);
    TH1D *hxClus                  = new TH1D("hxClusReco", "; X [cm]; Counts", 100,80,130);
    TH1D *hCpadMean               = new TH1D("hCpadMeanReco", "; cPad; Counts", 120,0,30);
    TH1D *hTgl                    = new TH1D("hLambdaTimeDirReco","; #lambda [deg]; Counts",100,-2,2);

    TH1F *hTrackAngleLocal        = new TH1F("hTrackAngleLocalReco","; local angle; Counts",100,0,10);
    TH1F *hTrackAngleGlobal       = new TH1F("hTrackAngleGlobalReco","; global angle; Counts",200,0,20);
    TH1F *hSlope                  = new TH1F("hSlopeReco","; slope; Counts",200,-0.4,0.4);
    TH1F *hOffset                 = new TH1F("hOffsetReco","; offset; Counts",100,-10,10);

    TH2D *hTrackResVsxCluster     = new TH2D("hTrackResVsxClusterReco","; x_{cluster} [cm]; Track residual [cm]", 125,80,130,100,-1,1);
    TH2D *hTrackResVsSigCluster   = new TH2D("hTrackResVsSigClusterReco","; #sigma_{cluster} [pad?]; Track residual [cm]", 100,0,2,100,-1,1);
    TH2D *hTrackResVsRelPadPos    = new TH2D("hTrackResVsRelPadPosReco","; relative pad position; Track residual [cm]", 100,0,1,100,-1,1);
    TH2D *hTrackResVsTrackAngle   = new TH2D("hTrackResVsTrackAngleReco","; Track angle [deg]; Track residual [cm]", 500,0,5,100,-1,1);
    TH2D *hSigClusterVsxCluster   = new TH2D("hSigClusterVsxClusterReco","; x_{cluster} [cm]; #sigma_{cluster} [pad?]", 125,80,130,100,0,2);
    TH2D *hSigClusterVsyCentre    = new TH2D("hSigClusterVsyCentreReco","; y_{pad center} [cm]; #sigma_{cluster} [pad?]", 100,15,35,100,-0.5,2);
    TH2D *hSigClusterVsyCluster   = new TH2D("hSigClusterVsyClusterReco","; y_{cluster} [cm]; #sigma_{cluster} [pad?]", 100,15,35,100,-0.5,2);
    TH2D *hSigClusterVsRelPadPos  = new TH2D("hSigClusterVsRelPadPosReco","; relative pad position; #sigma_{cluster} [pad?]", 100,0,1,100,-0.5,2);

    TH1F *hTrackResOnePadClus     = new TH1F("hTrackResOnePadClusterReco",";Track residual (only one pad clusters); Counts",150,-1.5,1.5);
    TH1F *hTrackResNoOnePadClus   = new TH1F("hTrackResNoOnePadClusterReco",";Track residual (no one pad clusters); Counts",150,-1.5,1.5);

    TH1F *hETRes                  = new TH1F("hETResReco",";Track residual;Counts;",100,-1,1);
    TH1F *hETPadMean              = new TH1F("hETPadMeanReco",";pad mean;Counts;",1000,0,100);
    TH1F *hETPad                  = new TH1F("hETPadReco",";pad;Counts;",100,0,100);
    TH1F *hETRelPad               = new TH1F("hETRelPadPositionReco",";Relative pad position;Counts;",100,0,1);

    TH1F *hRelPadPos              = new TH1F("hRelPadPosReco",";relative pad position;Counts",400,-2,2);

    TH1D *hCherenkov              = new TH1D("hCherenkovReco", "; Cherenkov counter signal; Counts", 100,0,0.06);
    TH1D *hCherenkovUsedPions     = new TH1D("hCherenkovUsedReco", "; Cherenkov counter signal; Counts", 100,0,0.06);
    TH1D *hCherenkovUsedEle       = new TH1D("hCherenkovUsedEleReco", "; Cherenkov counter signal; Counts", 100,0,0.06);

    TH2D *hPadOccupancy           = new TH2D("hOccReco", "; Row; Pad; Counts", 63,0,63,100,0,100);
    TH2D *hPadOccupancyCPad       = new TH2D("hOccCPadReco", "; Row; Pad; Counts", 63,0,63,100,-50,50);
    TH2D *hPadOccupancyUsedCPad   = new TH2D("hOccUsedReco", "; Row; Pad; Counts", 63,0,63,100,-50,50);
    TH2D *hPadOccupancyUsed       = new TH2D("hOccUsedCPadReco", "; Row; Pad; Counts", 63,0,63,100,0,100);
    TH2D *hPadOccupancyUsedPions  = new TH2D("hOccUsedPionsReco", "; Row; Pad; Counts", 63,0,63,100,-50,50);
    TH2D *hPadOccupancyUsedEle    = new TH2D("hOccUsedEleReco", "; Row; Pad; Counts", 63,0,63,100,-50,50);
    TH2D *hPadOccupancySetting0   = new TH2D("hOccS0Reco","; Row; Pad; Counts",63,0,63,100,-50,50);
    TH2D *hPadOccupancySetting1   = new TH2D("hOccS1Reco","; Row; Pad; Counts",63,0,63,100,-50,50);

    TH2D *hQTimeMax               = new TH2D("hQTimeReco","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);
    TH2D *hQTimeMaxUsed           = new TH2D("hQTimeUsedReco","; TimeBin; Q_{tot}; Counts",100,0,20,300,0,1000);

    TH2D *hQmaxRelPadPosUsed      = new TH2D("hQmaxRelPadPosReco", "; relative position of the cluster; Q_{max}; Counts", 100,-.2,1.2,600,0,1200);

    TH2D *hQtotPiRuns             = new TH2D("hQtotPiRunsReco", "; Run; Q_{tot}", 350,0,350,300,0,300);
    TH2D *hQtotEleRuns            = new TH2D("hQtotEleRunsReco", "; Run; Q_{tot}", 350,0,350,300,0,300);

  //  TH2D *hQmaxPiRuns             = new TH2D("hQmaxPiRuns", "; Run; Q_{max}", 350,0,350,167,0,200);
  //  TH2D *hQmaxEleRuns            = new TH2D("hQmaxEleRuns", "; Run; Q_{max}", 350,0,350,167,0,200);
    TH2D *hQmaxPiRuns             = new TH2D("hQmaxPiRunsReco", "; Run; Q_{max}; Normalized counts", 350,0,350,200,0,200);
    TH2D *hQmaxEleRuns            = new TH2D("hQmaxEleRunsReco", "; Run; Q_{max}; Normalized counts", 350,0,350,200,0,200);

    TH2D *hQmaxRowS0              = new TH2D("hQmaxRowS0Reco","; Row; Q_{max}; Counts",63,0,63,100,0,400);
    TH2D *hQmaxRowS1              = new TH2D("hQmaxRowS1Reco","; Row; Q_{max}; Counts",63,0,63,100,0,400);
    TH2D *hQmaxRowS0Corr          = new TH2D("hQmaxRowS0CorrReco","; Row; Q_{max}; Counts",63,0,63,100,0,400);
    TH2D *hQmaxRowS1Corr          = new TH2D("hQmaxRowS1CorrReco","; Row; Q_{max}; Counts",63,0,63,100,0,400);

    TH1F *hClustersVsRowS0All     = new TH1F("hClustersVsRowS0AllReco","; Pad row; Clusters",63,0,63);
    TH1F *hClustersVsRowS1All     = new TH1F("hClustersVsRowS1AllReco","; Pad row; Clusters",63,0,63);
    TH1F *hClustersVsRowS0UsedPiTracks     = new TH1F("hClustersVsRowS0UsedPiTracksReco","; Pad row; Clusters",63,0,63);
    TH1F *hClustersVsRowS1UsedPiTracks     = new TH1F("hClustersVsRowS1UsedPiTracksReco","; Pad row; Clusters",63,0,63);
    TH1F *hClustersVsRowS0UsedEleTracks     = new TH1F("hClustersVsRowS0UsedEleTracksReco","; Pad row; Clusters",63,0,63);
    TH1F *hClustersVsRowS1UsedEleTracks     = new TH1F("hClustersVsRowS1UsedEleTracksReco","; Pad row; Clusters",63,0,63);

    TGraph *gTrack = new TGraph();
    TGraph *gCluster = new TGraph();

    TGraph *gETRes = new TGraph();
    TGraph *gETPadMean = new TGraph();
    TGraph *gETPad = new TGraph();
    TGraph *gETRelPad = new TGraph();

    TH1F *hdEdxEleTot[350];
    TH1F *hdEdxPionTot[350];


  //* =============================================================================
  //* ================== Loop over events and apply all cuts ======================
  //* ======================= meanwhile fill histograms ===========================
  //* =============================================================================


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
    int ncl = 0;
    int ncledge = 0;
    int runcounter = 0;
    int eventcounter = 0;
    int firsttrack = 0;




    std::vector<TrackTPC> *vecEvent = 0;

    cout<<endl<<endl<<endl<<"Number of files to process: "<<chainEntries->GetEntriesFast()<<endl<<endl<<endl;


  //*=============================================== Loop over runs =========================================================
    for (int ifile=0; ifile<chainEntries->GetEntriesFast(); ++ifile){
      TFile *TreeFile = new TFile(Form("%s", chainEntries->At(ifile)->GetTitle()));
      cout<<endl<<endl<<"processing file Nr. "<<ifile+1<<" : "<<chainEntries->At(ifile)->GetTitle()<<endl;//<<endl;
      TTree *tree = (TTree*)TreeFile->Get("events");
      //TTree *tree = (TTree*)chainEntries->At(ifile)->Get("events");

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
  //*=============================================== Loop over events =========================================================

      for (int iev=0; iev<tree->GetEntriesFast(); ++iev){
        tree->GetEntry(iev);

        int NTracks = vecEvent->size();
  //*========================================== CUT ================================================
        if (NTracks != TrPerEv) continue;																// only one-track events
        CherenkovValue = Header.cherenkovValue;
        runNr = Header.run;



        Tracks += NTracks;

        hNIROCtracks->Fill(NTracks);

        nclEvent = 0;
        for (auto trackObject : *vecEvent) {
          std::vector<Cluster> clCont;
          trackObject.getClusterVector(clCont);
          nclTrack = clCont.size();
          nclEvent += nclTrack;
          //nclRun += nclTrack;
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
            //trackObject.setGainMap(GainMap, setting);
            std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
          }
          else if (loopcounter2 == 1 && usegainmap) {
            //trackObject.setGainMap(GainMap, setting);
            std::cout<<std::endl<<"set gainmap for setting "<< setting <<std::endl;
          }
          ++loopcounter;
  //*========================================== CUT ================================================
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
            float timeMean = clusterObject.getTimeMean();
            float Q = clusterObject.getQ();
            float Qmax = clusterObject.getQmax();
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
            hQmax->Fill(Qmax);
            hPadOccupancy->Fill(row,pad);
            hPadOccupancyCPad->Fill(row,cpad);
            if (runNr <= 243) {
              hClustersVsRowS0All->Fill(row);
              hPadOccupancySetting0->Fill(row,cpad);
              hQmaxRowS0->Fill(row,Qmax);
              if (GainMap != "") {
                float QmaxCorr = Qmax/gainmapS0->getValue(roc,row,pad);
                hQmaxRowS0Corr->Fill(row,QmaxCorr);
              }
            }
            else if (runNr > 255) {
              hClustersVsRowS1All->Fill(row);
              hPadOccupancySetting1->Fill(row,cpad);
              hQmaxRowS1->Fill(row,Qmax);
              if (GainMap != "") {
                float QmaxCorr = Qmax/gainmapS1->getValue(roc,row,pad);
                hQmaxRowS1Corr->Fill(row,QmaxCorr);
              }
            }
          }




  //*========================================== CUT ================================================
          if (ncl < nclCut) continue;																	// cut on number of clusters per track

          hCherenkov->Fill(CherenkovValue);
  //*========================================== CUT ================================================
          if (CherenkovValue >= CherCutLow && CherenkovValue <= CherCutHigh) continue;					// PID via Cherenkov
          if (CherenkovValue < CherCutLow) {
            hCherenkovUsedPions->Fill(CherenkovValue);
          }
          else if (CherenkovValue > CherCutHigh){
            hCherenkovUsedEle->Fill(CherenkovValue);
          }

          for (auto &clusterObject : clCont) {															// make cuts on cluster properties
            DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
            float row = pos.getPadSecPos().getPadPos().getRow();
  //          if (setting == 0) {
  //            if (row < 2 || row >= 57) continue;
  //          }
            float pad = pos.getPadSecPos().getPadPos().getPad();
            float timeMean = clusterObject.getTimeMean();
            float cpad = GetCPad(clusterObject);


  //*========================================== CUT ================================================
            if (timeMean < timeMeanLow || timeMean > timeMeanHigh){isok = false; break;}					// cut on time max, only clusters from a certain range on z-axis can come from a particle

            if (!InsideEdge(row, cpad, setting, 3)){++ncledge;}
          }

  //*========================================== CUT ================================================
          if (((float(ncledge)/float(ncl)) > nclFracoutofCPad) && excludeEdge) continue;					// cut on detector geometry


          if (isok == true){																			// track accepted
            usedcl += ncl;

            float dEdxTot = 0;
            float dEdxMax = 0;

            float angle = asin(trackObject.getSnp())*TMath::RadToDeg();

            if (setting == 0) {
              dEdxTot = trackObject.getTruncatedMean(.0,.7,1);
              dEdxMax = trackObject.getTruncatedMean(.0,.7,0);
            }
            else if (setting == 1) {
              dEdxTot = trackObject.getTruncatedMean(.0,.7,1);
              dEdxMax = trackObject.getTruncatedMean(.0,.7,0);
            }


            if (CherenkovValue < CherCutLow){
              hdEdxPionTot[runcounter]->Fill(dEdxTot);
              hdEdxPionTotAll->Fill(dEdxTot);
              hdEdxPionMax->Fill(dEdxMax);
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
              hdEdxEleMax->Fill(dEdxMax);
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

            int clustercounter = 0;
            for (auto &clusterObject : clCont) {															// loop over clusters
              const CRU cru(clusterObject.getCRU());

              const PadRegionInfo& region = mapper.getPadRegionInfo(cru.region());
              const int rowInSector = clusterObject.getRow() + region.getGlobalRowOffset();
              const GlobalPadNumber pad = mapper.globalPadNumber(PadPos(rowInSector, clusterObject.getPadMean()));
              const PadCentre& padCentre = mapper.padCentre(pad);
              const float localYfactor = (cru.side() == Side::A) ? -1.f : 1.f;
              float zPosition = sampa.getZfromTimeBin(clusterObject.getTimeMean(), cru.side());

              LocalPosition3D posLoc(padCentre.X(), localYfactor * padCentre.Y(), zPosition);
              GlobalPosition3D posGlob = Mapper::LocalToGlobal(posLoc, cru.sector());

              DigitPos pos(clusterObject.getCRU(), PadPos(clusterObject.getRow(), clusterObject.getPadMean()));
              float row = pos.getPadSecPos().getPadPos().getRow();

              if (setting == 0) {
                if (row <=2 || row >= 57) continue;
              }
              else if (setting == 1) {
                if (row <= 2 || row >= 61) continue;
              }
              //float pad = pos.getPadSecPos().getPadPos().getPad();
              const int cpad = GetCPad(clusterObject);
              float timeMean = clusterObject.getTimeMean();
              float Q = clusterObject.getQ();
              float Qmax = clusterObject.getQmax();
              float padMean = clusterObject.getPadMean();
              float padSigma = clusterObject.getPadSigma();
              float padLoc = pos.getPadSecPos().getPadPos().getPad();

              int clusCRU = clusterObject.getCRU();
              float timeMeanSigma = clusterObject.getTimeSigma();

              for (auto &clusterObject2 : clCont) {
                DigitPos pos(clusterObject2.getCRU(), PadPos(clusterObject2.getRow(), clusterObject2.getPadMean()));
                float row2 = pos.getPadSecPos().getPadPos().getRow();
                if (setting == 0) {
                  if (row2-1 <=2 || row2+1 >= 57) continue;
                }
                else if (setting == 1) {
                  if (row2-1 <= 2 || row2+1 >= 61) continue;
                }
                if (row2 == row-1) {
                  hQNeighboringRowLeft->Fill(Q,clusterObject2.getQ());
                  hQmaxNeighboringRowLeft->Fill(Qmax,clusterObject2.getQmax());
                }
                else if (row2 == row+1) {
                  hQNeighboringRowRight->Fill(Q,clusterObject2.getQ());
                  hQmaxNeighboringRowRight->Fill(Qmax,clusterObject2.getQmax());
                }
              }

              float relPadPos = padMean - padLoc;

              float xClus = posGlob.X();
              float yClus = posGlob.Y()-0.2;//+relPadPos*0.4;
              float yCentreGlob = posGlob.Y();

              float CpadMean = cpad + relPadPos;

              //float yClus = 0.2 + 0.4 * CpadMean;
              //float xClus = 84.1 + 0.375 + 0.75 * row;
              float yTrack = evalX(trackObject, xClus);

              float residual = yClus - yTrack;

              hTrackResiduals->Fill(residual);

              hCRU->Fill(clusCRU);
              hRow->Fill(row);
              hTimeMean->Fill(timeMean);
              hTimeSigma->Fill(timeMeanSigma);
              hPadMean->Fill(padMean);
              hPadSigma->Fill(padSigma);
              hQtotVsTimeMean->Fill(timeMean,Q);
              hQmaxVsTimeMean->Fill(timeMean,Qmax);

              hTrackResVsxCluster->Fill(xClus,residual);
              hTrackResVsSigCluster->Fill(padSigma,residual);
              hTrackResVsRelPadPos->Fill(relPadPos,residual);
              hTrackResVsTrackAngle->Fill(angle,residual);
              hSigClusterVsxCluster->Fill(xClus,padSigma);
              hSigClusterVsyCentre->Fill(yCentreGlob,padSigma);
              hSigClusterVsyCluster->Fill(yClus,padSigma);
              hSigClusterVsRelPadPos->Fill(relPadPos,padSigma);

              hEvalX->Fill(yTrack);
              hCpadMean->Fill(CpadMean);
              hxClus->Fill(xClus);
              hyClus->Fill(yClus);

              hTgl->Fill(atan(trackObject.getTgl())*TMath::RadToDeg());

              hRelPadPos->Fill(relPadPos);

              hQTimeMaxUsed->Fill(timeMean,Q);
              hQUsed->Fill(Q);
              hQmaxUsed->Fill(Qmax);
              hPadOccupancyUsed->Fill(row,pad);
              hPadOccupancyUsedCPad->Fill(row,cpad);
              hQmaxVsQtot->Fill(Qmax,Q);

              if (padSigma == 0) {
                hTrackResOnePadClus->Fill(residual);
              }
              else {
                hTrackResNoOnePadClus->Fill(residual);
              }

              //if (Qmax > Q) {
              //  hQmaxVsQtot->Fill(Qmax,Q,100000000);
              //}

              if (firsttrack == 1) {
                gCluster->SetPoint(clustercounter,xClus,yClus);
                gETRes->SetPoint(clustercounter,xClus,residual);
                gETPadMean->SetPoint(clustercounter,xClus,padMean);
                gETPad->SetPoint(clustercounter,xClus,padLoc);
                gETRelPad->SetPoint(clustercounter,xClus,relPadPos);
                hETRes->Fill(residual);
                hETPadMean->Fill(padMean);
                hETPad->Fill(padLoc);
                hETRelPad->Fill(relPadPos);
              }

              ++clustercounter;
            }

            if (firsttrack == 1) {
              for (int step = 0; step < 500; ++step) {
                gTrack->SetPoint(step,step*0.375,evalX(trackObject, step*0.375));
              }
            }

          } // end of 'if true'
          ++firsttrack;
        } // end of loop over tracks ///break; ///if (isok==true){break;} //use only one track
        ++eventcounter;
      } // end of loop over events

       if (hdEdxPionTot[runcounter]->GetEntries() > 0 && hdEdxEleTot[runcounter]->GetEntries() > 0){
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


        TF1 *pionfitMax = new TF1("pionfitMax","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
        TF1 *electronfitMax = new TF1("electronfitMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
        hdEdxEleMax->SetLineColor(kRed);
        hdEdxPionMax->SetLineColor(kBlue);

        GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
        hdEdxPionMax->Fit("pionfitMax","","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
        GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
        //hdEdxEleMax->Fit("electronfitMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

        hdEdxEleMax->Fit("electronfitMax");


        hdEdxEleMax->GetFunction("electronfitMax")->SetLineColor(kRed);Mapper.h
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

       } // end of if getentries > 0
      delete TreeFile;
      cout<<"Tree deleted"<<endl;
      ++runcounter;


    } // end of loop over runs


      TF1 *pionfitEnd = new TF1("pionfitEnd","gaus",hdEdxPionTotAll->GetXaxis()->GetXmin(),hdEdxPionTotAll->GetXaxis()->GetXmax());
      TF1 *electronfitEnd = new TF1("electronfitEnd","gaus",hdEdxEleTotAll->GetXaxis()->GetXmin(),hdEdxEleTotAll->GetXaxis()->GetXmax());
      hdEdxEleTotAll->SetLineColor(kRed);
      hdEdxPionTotAll->SetLineColor(kBlue);

      const Float_t frac=0.2;
      Int_t bin1=0,bin2=0;

      GetBinMinMax(hdEdxPionTotAll,frac,bin1,bin2);
      hdEdxPionTotAll->Fit("pionfitEnd","0","",hdEdxPionTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionTotAll->GetXaxis()->GetBinUpEdge(bin2));
      GetBinMinMax(hdEdxEleTotAll,frac,bin1,bin2);
      //hdEdxEleTotAll->Fit("electronfitEnd","","",hdEdxEleTotAll->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleTotAll->GetXaxis()->GetBinUpEdge(bin2));

      //alternative fit
      //hdEdxPionTotAll->Fit("pionfitEnd");
      hdEdxEleTotAll->Fit("electronfitEnd","0");


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



      TF1 *pionfitMaxEnd = new TF1("pionfitMaxEnd","gaus",hdEdxPionMax->GetXaxis()->GetXmin(),hdEdxPionMax->GetXaxis()->GetXmax());
      TF1 *electronfitEndMax = new TF1("electronfitEndMax","gaus",hdEdxEleMax->GetXaxis()->GetXmin(),hdEdxEleMax->GetXaxis()->GetXmax());
      hdEdxEleMax->SetLineColor(kRed);
      hdEdxPionMax->SetLineColor(kBlue);

      GetBinMinMax(hdEdxPionMax,frac,bin1,bin2);
      hdEdxPionMax->Fit("pionfitMaxEnd","0","r",hdEdxPionMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxPionMax->GetXaxis()->GetBinUpEdge(bin2));
      GetBinMinMax(hdEdxEleMax,frac,bin1,bin2);
      //hdEdxEleMax->Fit("electronfitEndMax","","r",hdEdxEleMax->GetXaxis()->GetBinLowEdge(bin1),hdEdxEleMax->GetXaxis()->GetBinUpEdge(bin2));

      hdEdxEleMax->Fit("electronfitEndMax","0");


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




    hCherenkov->SetLineColor(kBlue);

    hCherenkovUsedPions->SetLineColor(kBlue);
    hCherenkovUsedEle->SetLineColor(kRed);

    hQmaxUsed->Scale(1/hQmaxUsed->GetEntries());
    hQUsed->Scale(1/hQUsed->GetEntries());
    hTrackResiduals->Scale(1/hTrackResiduals->GetEntries());
    hdEdxPionTotAll->Scale(1/hdEdxPionTotAll->GetEntries());
    hdEdxEleTotAll->Scale(1/hdEdxEleTotAll->GetEntries());
    hdEdxPionMax->Scale(1/hdEdxPionMax->GetEntries());
    hdEdxEleMax->Scale(1/hdEdxEleMax->GetEntries());
    hTrackResOnePadClus->Scale(1/hTrackResiduals->GetEntries());
    hTrackResNoOnePadClus->Scale(1/hTrackResiduals->GetEntries());
    hCRU->Scale(1/hCRU->GetEntries());
    hRow->Scale(1/hRow->GetEntries());
    hPadMean->Scale(1/hPadMean->GetEntries());
    hPadSigma->Scale(1/hPadSigma->GetEntries());
    hTimeMean->Scale(1/hTimeMean->GetEntries());
    hTimeSigma->Scale(1/hTimeSigma->GetEntries());

    hTrackResiduals->SetLineColor(kCyan+2);
    hTrackResOnePadClus->SetLineColor(kMagenta+2);
    hTrackResNoOnePadClus->SetLineColor(kYellow+2);

    gTrack->SetMarkerStyle(20);
    gTrack->SetMarkerColor(kGreen+2);
    gCluster->SetMarkerStyle(20);
    gCluster->SetMarkerColor(kBlue+2);
    gTrack->GetXaxis()->SetTitle("x [cm]");
    gTrack->GetYaxis()->SetTitle("y [cm]");
    gCluster->GetXaxis()->SetTitle("x [cm]");
    gCluster->GetYaxis()->SetTitle("y [cm]");

    gETRes->SetMarkerStyle(20);
    gETPadMean->SetMarkerStyle(20);
    gETPad->SetMarkerStyle(20);
    gETRelPad->SetMarkerStyle(20);
    gETRes->SetMarkerColor(kRed+2);
    gETPadMean->SetMarkerColor(kCyan+2);
    gETPad->SetMarkerColor(kMagenta+2);
    gETRelPad->SetMarkerColor(kYellow+2);
    gETRes->GetXaxis()->SetTitle("x [cm]");
    gETPadMean->GetXaxis()->SetTitle("x [cm]");
    gETPad->GetXaxis()->SetTitle("x [cm]");
    gETRelPad->GetXaxis()->SetTitle("x [cm]");
    gETRes->GetYaxis()->SetTitle("Track residual");
    gETPadMean->GetYaxis()->SetTitle("Pad mean");
    gETPad->GetYaxis()->SetTitle("Pad");
    gETRelPad->GetYaxis()->SetTitle("Relative pad position");

    //gTrack->GetXaxis()->SetRangeUser(0,50);
    //gTrack->GetYaxis()->SetRangeUser(0,10);

    //gCluster->GetXaxis()->SetRangeUser(0,50);
    //gCluster->GetYaxis()->SetRangeUser(0,10);

    hQNeighboringRowLeft->Add(hQNeighboringRowRight);
    TProfile *pQNeighboringRowReco    = (TProfile*)hQNeighboringRowLeft->ProfileX();
    hQmaxNeighboringRowLeft->Add(hQmaxNeighboringRowRight);
    TProfile *pQmaxNeighboringRowReco    = (TProfile*)hQmaxNeighboringRowRight->ProfileX();

    //TrackAna->Fill();
    OutFile->Write();
    OutFile->WriteObject(gTrack, "gTrack");
    OutFile->WriteObject(gCluster, "gCluster");
    OutFile->WriteObject(gETRes, "gETResidual");
    OutFile->WriteObject(gETPadMean, "gETPadMean");
    OutFile->WriteObject(gETPad, "gETPad");
    OutFile->WriteObject(gETRelPad, "gETRelPadPosition");
    OutFile->WriteObject(pQNeighboringRowReco, "pQNeighboringRowReco");
    OutFile->WriteObject(pQmaxNeighboringRowReco, "pQmaxNeighboringRowReco");
    //OutFile->WriteObject(pQNeighboringRowRight, "pQNeighboringRowRightReco");

    OutFile->Close();

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

    return 0;
  }
}

