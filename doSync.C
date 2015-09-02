{

  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");

  //gROOT->ProcessLine(".L GetCutFlow.C+"); 
  // gROOT->ProcessLine(".L SyncCutflow.C+");
  const unsigned int chainsize = 2;
  TChain *ch[chainsize];
  string dataset[chainsize];


  //original "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_2jskim/
  //intermed "/nfs-7/userdata/stop2015/babies_4May2015/"
  // /nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1lep_2jskim/
  //string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1lep_2jskim__fullLepSelection__20150527/";
  string babylocation = "/nfs-7/userdata/mliu/";
  string helper;
  dataset[0] = "ttbarskim15";
  ch[0] = new TChain("t"); 
  helper = babylocation+"ttbar_skim15.root";      ch[0]->Add(helper.c_str());
  dataset[1] = "ttbarskim30";
  ch[1] = new TChain("t"); 
  helper = babylocation+"ttbar_skim30.root";      ch[1]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
    ScanChain(ch[i],true,-1,dataset[i]); 
  }

}
