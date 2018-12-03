#include "TH1F.h"

double convert(double x)
{
  return ((-67./4035.)+(0.1/(470.5-67)*x));
}

void draw_nele_plot_Andi() {
  TH1F *nele = new TH1F("hnele",";Total number of electrons per cm;Normalized counts",75,0,150);

  nele->Fill(0.,convert(67));
  nele->Fill(2.,convert(67.84));
  nele->Fill(4.,convert(68.66));
  nele->Fill(6.,convert(73.6));
  nele->Fill(8.,convert(90));
  nele->Fill(10.,convert(118));
  nele->Fill(12.,convert(166));
  nele->Fill(14.,convert(200));
  nele->Fill(16.,convert(252));
  nele->Fill(18.,convert(288));
  nele->Fill(20.,convert(292));
  nele->Fill(22.,convert(311));
  nele->Fill(24.,convert(313));
  nele->Fill(26.,convert(291));
  nele->Fill(28.,convert(271));
  nele->Fill(30.,convert(257));
  nele->Fill(32.,convert(232.5));
  nele->Fill(34.,convert(206.5));
  nele->Fill(36.,convert(196));
  nele->Fill(38.,convert(187.5));
  nele->Fill(40.,convert(167));
  nele->Fill(42.,convert(165));
  nele->Fill(44.,convert(148));
  nele->Fill(46.,convert(133));
  nele->Fill(48.,convert(130));
  nele->Fill(50.,convert(123));
  nele->Fill(52.,convert(113));
  nele->Fill(54.,convert(110));
  nele->Fill(56.,convert(106));
  nele->Fill(58.,convert(101.5));
  nele->Fill(60.,convert(97.8));
  nele->Fill(62.,convert(98.1));
  nele->Fill(64.,convert(91.5));
  nele->Fill(66.,convert(91.4));
  nele->Fill(68.,convert(90.2));
  nele->Fill(70.,convert(90));
  nele->Fill(72.,convert(83.5));
  nele->Fill(74.,convert(86.7));
  nele->Fill(76.,convert(81.6));
  nele->Fill(78.,convert(82.6));
  nele->Fill(80.,convert(80.9));
  nele->Fill(82.,convert(79.6));
  nele->Fill(84.,convert(79.7));
  nele->Fill(86.,convert(79.4));
  nele->Fill(88.,convert(79.6));
  nele->Fill(90.,convert(78.7));
  nele->Fill(92.,convert(76.8));
  nele->Fill(94.,convert(77.7));
  nele->Fill(96.,convert(77.2));
  nele->Fill(98.,convert(76.9));
  nele->Fill(100.,convert(76.4));
  nele->Fill(102.,convert(75.6));
  nele->Fill(104.,convert(76.3));
  nele->Fill(106.,convert(75.1));
  nele->Fill(108.,convert(74.8));
  nele->Fill(110.,convert(74.4));
  nele->Fill(112.,convert(73.2));
  nele->Fill(114.,convert(74.1));
  nele->Fill(116.,convert(72.9));
  nele->Fill(118.,convert(71.8));
  nele->Fill(120.,convert(72.8));
  nele->Fill(122.,convert(73));
  nele->Fill(124.,convert(72));
  nele->Fill(126.,convert(72.8));
  nele->Fill(128.,convert(71.6));
  nele->Fill(130.,convert(73));
  nele->Fill(132.,convert(71.2));
  nele->Fill(134.,convert(72));
  nele->Fill(136.,convert(72.2));
  nele->Fill(138.,convert(72.2));
  nele->Fill(140.,convert(72));
  nele->Fill(142.,convert(70.6));
  nele->Fill(144.,convert(71));
  nele->Fill(146.,convert(71));
  nele->Fill(148.,convert(72.5));


  TFile *OutFile = new TFile("/home/tom/myana/studySimVsReco/Results_studyHits/Nele_AliRoot_Andi.root", "recreate");

  OutFile->WriteObject(nele, "NElePerCM");

  OutFile->Close();


}
