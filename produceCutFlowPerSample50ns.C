{
  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");

  gROOT->ProcessLine(".L GetCutFlow.C+");
  const unsigned int chainsize = 5;
  TChain *ch[chainsize];
  string dataset[chainsize];
//  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/";
//  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V728_met30_ge1l_ge2j_Phys14/";
  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_50ns/";
  string helper;
  dataset[0] = "T_tbarW_5f_powheg50ns";
  ch[0] = new TChain("t"); 
  helper = babylocation+"T_tbarW_5f_powheg50ns.root"; ch[0]->Add(helper.c_str());

  dataset[1] = "T_tch_4f_amc50ns";
  ch[1] = new TChain("t");
  helper = babylocation+"T_tch_4f_amc50ns.root";    ch[1]->Add(helper.c_str());

  dataset[2] = "T_tch_4f_powheg50ns";
  ch[2] = new TChain("t");
  helper = babylocation+"T_tch_4f_powheg50ns.root"; ch[2]->Add(helper.c_str());

  dataset[3] = "Tbar_tch_4f_powheg50ns";
  ch[3] = new TChain("t");
  helper = babylocation+"Tbar_tch_4f_powheg50ns.root";    ch[3]->Add(helper.c_str());

  dataset[4] = "ttbar_amc50ns";
  ch[4] = new TChain("t");
  helper = babylocation+"ttbar_amc50ns.root"; ch[4]->Add(helper.c_str());

 // helper = babylocation+"stop_650_325.root";      ch[6]->Add(helper.c_str());
 // dataset[7] = "Stop_850_100";
 // ch[7] = new TChain("t"); 
 // helper = babylocation+"stop_850_100.root";      ch[7]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
 //  if(i==3||i==4) continue; 
   ScanChain(ch[i],true,-1,dataset[i]); 
  }
//  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
 // gROOT->ProcessLine(".x makePlotPlot.C+");
  //gROOT->ProcessLine(".x makePlotCutFlow.C+");
  gROOT->ProcessLine(".x GetSRTable.C+");
  //gROOT->ProcessLine(".x GetCFNumbers.C+");
  //gROOT->ProcessLine(".x makePlotIP.C+");
}
