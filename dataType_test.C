#include <iostream>
using namespace std;

void dataType_test() {
  double qTot = 0;
  double charge = 0;
  double meanP = 0;
  double meanT = 0;
  double sigmaP = 0;
  double sigmaT = 0;

  short deltaT, deltaP;
  short t, p;
  for (t = 0; t < 5; ++t) {
    deltaT = t - 5/2;
    cout << endl << "time bin: " << t << "\t" << "delta t: " << deltaT;
    for (p = 0; p < 5; ++p) {
      deltaP = p - 5/2;
      cout << endl << "p: " << p << "\t" << "delta p: " << deltaP;
      qTot += charge;

      meanP += charge * deltaP;
      meanT += charge * deltaT;

      sigmaP += charge * deltaP*deltaP;
      sigmaT += charge * deltaT*deltaT;
    }
  }
}
