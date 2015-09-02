{
  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");

  gROOT->ProcessLine(".L GetHistoW_selection.C+");
  
  const unsigned int chainsize = 3;
  TChain *ch[chainsize];
  string dataset[chainsize];
  string babylocation = "/nfs-7/userdata/mliu/Comissioning/";
  string helper;
  dataset[0] = "TTbar";
  ch[0] = new TChain("t"); 
  helper = babylocation+"ttbar.root"; ch[0]->Add(helper.c_str());
  dataset[1] = "WJets";
  ch[1] = new TChain("t");
  helper = babylocation+"t_sch.root";    ch[1]->Add(helper.c_str());
  dataset[2] = "ZJets";
  ch[2] = new TChain("t");
  helper = babylocation+"t_sch.root";    ch[2]->Add(helper.c_str());
  dataset[3] = "data";
  ch[3] = new TChain("t");
  helper = babylocation+"t_sch.root";    ch[3]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
    ScanChain(ch[i],true,-1,dataset[i]); 
  }
  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
  //gROOT->ProcessLine(".x makePlotOneCut.C+");
  //gROOT->ProcessLine(".x makePlotCutFlow.C+");
  //gROOT->ProcessLine(".x GetSRTable.C+");
  //gROOT->ProcessLine(".x GetCFNumbers.C+");
  //gROOT->ProcessLine(".x makePlotIP.C+");
}
