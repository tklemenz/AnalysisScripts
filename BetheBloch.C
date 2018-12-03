#include "TPCSimulation/Detector.h"
#include "TPCBase/ParameterGas.h"

#include <iostream>

using std::cout;
using std::endl;

using namespace o2::TPC;

float BetheBlochAleph(float bg, float kp1, float kp2, float kp3, float kp4, float kp5)
{
  float beta = bg/std::sqrt(1.+ bg*bg);
  float aa = std::pow(beta,kp4);
  float bb = std::pow(1./bg,kp5);
  bb=std::log(kp3+bb);

  return (kp2-aa-bb)*kp1/aa;
}

float BetheBloch(std::string Particle = "", int BeamMomentum = 1, float gain = 2000)
{
  const static ParameterGas &gasParam = ParameterGas::defaultInstance();

  float Mpi = 139.57018*1e6;               // eV/c^2
  float Mele = 0.5109989461*1e6;           // eV/c^2
  float Mproton = 938.2720813*1e6;         // eV/c^2
  float beamMomentum = BeamMomentum*1e9;   // eV/c

  float NElectronsPad = 0;
  float m = 0;
  float betaGamma = 0;
  float ConversionGain = 2*1e13;           // Volt/Coulomb
//  float ConversionGain = 1.2*1e13;
  float ChargeTot = 0;
  float ChargeEle = 1.6021766208*1e-19;    // Coulomb
  float ADCChannels = 1024;
  float ADCRange = 2.2;                    // Volt
//  float ADCRange = 2.;
  float ADCValue = 0;
  float NeCON = 36.1;                      // electrons per cm produced by MIP

  if (Particle == "pi") {
    m = Mpi;
  }
  else if (Particle == "e") {
    m = Mele;
  }
  else if (Particle == "proton") {
    m = Mproton;
  }
  else {
    cout << endl << endl << "Choose either pi or e as particle" << endl << endl;
    return 0;
  }

  betaGamma = beamMomentum/m;
//  NElectronsPad = gain * 0.7513 * gasParam.getNprim() * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));
//  NElectronsPad = gain * 20;



//  NElectronsPad = gain * 0.7513 * NeCON * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));

  NElectronsPad = gain * NeCON * 0.75 * BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4));

  ChargeTot = NElectronsPad * ChargeEle;
  ADCValue = ConversionGain * ChargeTot * 0.64/(ADCRange/ADCChannels);          ///0.64 is the ratio of Qmax/Qtot coming from the sampling for 5 MHz
//  ADCValue = ConversionGain * ChargeTot * 0.3/(ADCRange/ADCChannels);
  cout << endl << endl << endl << "====================================================" << endl << endl;
  cout << "particle: " << Particle << ", beam momentum: " << BeamMomentum << ", gain: " << gain << endl << endl;
  cout << "expected ADC value: " << ADCValue << endl << endl;
  cout << "betaGamma: " << betaGamma << endl << endl;
  cout << "electrons at pad: " << NElectronsPad << endl << endl;
  cout << "total charge at pad: " << ChargeTot << endl << endl;
  cout << "signal per ADC channel: " << ADCRange/ADCChannels << endl << endl;
  cout << "conversion voltage: " << ConversionGain*ChargeTot << endl << endl;
  cout << "Bethe Bloch correction: " << BetheBlochAleph(betaGamma, gasParam.getBetheBlochParam(0), gasParam.getBetheBlochParam(1), gasParam.getBetheBlochParam(2), gasParam.getBetheBlochParam(3), gasParam.getBetheBlochParam(4))<< endl << endl;

  cout << "====================================================" << endl << endl << endl;

  return ADCValue;
}
