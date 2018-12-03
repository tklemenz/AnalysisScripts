static UChar_t GetNPads(UChar_t row) { return (row==0) ? 68 : 2 *UChar_t(Double_t(row)/3. +33.67); }

void CompareGainMaps (const char *Map1 /*to be tested*/, const char *Map2 /* reference map*/){

gSystem->Load("~/GEMtest/code/libGEMEvent.so");

TFile *mapFile1 = new TFile(Form("%s", Map1));
GainMap1 = (TH2D*)mapFile1->Get("GainMap");

TFile *mapFile2 = new TFile(Form("%s", Map2));
GainMap2 = (AliTPCCalROC*)mapFile2->Get("gainMap");



const int nrows = 74;
const int npads = 119;
double frac = 0;
double TestValue = 0;
double RefValue = 0;
int counter = 0;
double sum = 0;
double mean = 0;


TH2D *hGainMapTest = new TH2D("hGainMapTest", "; Row; Pad", 63,0,62,35,-5,30);

cout<<endl<<endl<<endl<<endl<<"========================== Test value / Reference value =========================="<<endl<<endl;
for(int irow=0; irow<nrows; ++irow){
  for(int ipad=0; ipad<npads; ++ipad){
	TestValue = GainMap1->GetBinContent(irow+6,ipad+60);
	RefValue = GainMap2->GetValue(irow,ipad+GetNPads(irow)/2);
	if(TestValue != 0 && RefValue != 0){
	  counter++;  
	  frac = TestValue/RefValue; 
	  hGainMapTest->Fill(irow,ipad,frac);
	  sum += frac;
	  cout<<frac<<"\t";
	  if (counter%9 == 0){
	    cout<<endl;
	  }
	  
    }
    else if(TestValue == 0 && RefValue == 0){
	  hGainMapTest->Fill(irow,ipad,0);
	}
	else if(TestValue == 0 && RefValue != 0){
	  hGainMapTest->Fill(irow,ipad,0.1);
	}
	else if(TestValue != 0 && RefValue == 0){
	  hGainMapTest->Fill(irow,ipad,0.2);
	}
  }
}
mean = sum/counter;
cout<<endl<<endl<<"__________________________________________________________________________________"<<endl<<endl;
cout<<"Mean:\t"<<mean;
cout<<endl<<endl<<"=================================================================================="<<endl<<endl;

TCanvas *c1 = new TCanvas("Test");
GainMap1->Draw("colz");

TCanvas *c2 = new TCanvas("Reference");
GainMap2->Draw("colz");

hGainMapTest->GetZaxis()->SetRangeUser(0.5,1.8);
TCanvas *c3 = new TCanvas("Frac");
hGainMapTest->Draw("colz");
}
