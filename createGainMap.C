void testTChain(TString FileList){
	gSystem->Load("/home/tklemenz/GEMtest/code/libGEMEvent.so");
	
	GEMEvent *fEvent=0x0;
    fEvent=0x0;
	
	fChain = new TChain("GEM");
	//fChain = 0x0;
	fChain->SetBranchAddress("Events",&fEvent);
	
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
		fChain->Add(file);
	}
	
																						//pointer rausfinden
	cout<<arr->At(0)->GetName()<<endl;
	cout<<endl;
	cout<<arr->At(1)->GetName()<<endl;
	cout<<endl;
	cout<<allFiles<<endl;
	cout<<arr->GetEntriesFast()<<endl;
	cout<<fChain->GetListOfFiles()<<endl;

}
