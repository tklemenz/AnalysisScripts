.X /home/tklemenz/GEMtest/code/macros/loadlibs.C
gSystem->Load("~/GEMtest/code/libGEMtest.so")
gSystem->Load("~/GEMtest/code/libGEMEvent.so")
.include /home/tklemenz/GEMtest/code/Event
.include /home/tklemenz/AliSoftware/aliroot/master/inst/include
.L code/Analysis/FillHistos.C
TFile *TreeFile = new TFile("~/ROC_Raw_data/data1499.root")
TTree *tree = (TTree*)TreeFile->Get("GEM")
FillHistos(tree, "/home/tklemenz/")





// open root in ~/GEMtest and just copy the lines (in AliEnv)
