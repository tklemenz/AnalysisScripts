void CompareGainMaps (const char *Map1 /*to be tested*/, const char *Map2 /* Jens' map*/){

gSystem->Load("~/GEMtest/code/libGEMEvent.so");

TFile *mapFile1 = new TFile(Form("%s", Map1));
GainMap1 = (TH2D*)mapFile1->Get("GainMap");

TFile *mapFile2 = new TFile(Form("%s", Map2));
GainMap2 = (TH2D*)mapFile2->Get("gainMap");


const int nrows = 63;
const int npads = 30;
double frac = 0;
double arr[nrows][npads];

cout<<GainMap2->GetBinContent(10,15)<<endl;

for (int ipad = 0; ipad < npads; ++ipad){
  int counter = 0;
  for (int irow = 0; irow < nrows; ++irow, ++counter){
	arr[irow][ipad] = -1;
	/*if (counter%62 == 0 && counter != 0 && ipad != 0){cout<<arr[irow][ipad]<<endl;}			//print formatted array to terminal with -1 entries
    else if (counter%62 == 0 && ipad == 0 && irow ==62){cout<<arr[irow][ipad]<<endl;}
	else if (irow == 62 && ipad == 29){cout<<endl;}    
	else {cout<<arr[irow][ipad];}*/ 	
  }
}

for (int irow = 0; irow < nrows; ++irow){
  for (int ipad = 0; ipad < npads; ++ipad){ 
    arr[irow][ipad] = GainMap1->GetBinContent(irow,ipad);
  }
}
/*cout<<setprecision(2);	
for (int ipad = 0; ipad < npads; ++ipad){														//print real entries to terminal
  int counter = 0;
  for (int irow = 0; irow < nrows; ++irow, ++counter){
	if (counter%62 == 0 && counter != 0 && ipad != 0){cout<<arr[irow][ipad]<<endl;}
    else if (counter%62 == 0 && ipad == 0 && irow ==62){cout<<arr[irow][ipad]<<endl;}
	else if (irow == 62 && ipad == 29){cout<<endl;}    
	else {cout<<arr[irow][ipad]<<"\t";} 	
  }
}*/



/*for (int ipad = 0; ipad < npads; ++ipad){
  for (int irow = 0; irow < nrows; ++irow){
	  cout<<"Row: "<<irow<<"  Pad: "<<ipad<<endl;
  }
}*/



TCanvas *c1 = new TCanvas("mine");
GainMap1->Draw("colz");

TCanvas *c2 = new TCanvas("Jens");
GainMap2->Draw("colz");
}
