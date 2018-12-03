#include "Riostream.h"
#include "TROOT.h"
#include "TSystem.h"


void test(){
float x,y,z;
x = 0;
y = 1;
z = 2;

TFile *f = new TFile("basic.root","RECREATE");
TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","x:y:z");

h1->Fill(x);
h1->Fill(y);
h1->Fill(z);

ntuple->Fill(x,y,z);
ntuple->Fill(x,y,z);

f->Write();
f->Close();

 std::cout<< std::endl << "test done" << std::endl;
}

