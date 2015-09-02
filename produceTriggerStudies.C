{
  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");
  gROOT->ProcessLine(".L GetHistosTriggerEff.C+");
  const unsigned int chainsize = 11;
  TChain *ch[chainsize];
  string dataset[chainsize];
  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/";
  string helper;

  dataset[0] = "T_sch_4f_amc25ns";
  ch[0] = new TChain("t"); 
  helper = babylocation+"T_sch_4f_amc25ns.root"; ch[0]->Add(helper.c_str());

  dataset[1] = "T_tbarW_5f_powheg25ns";
  ch[1] = new TChain("t");
  helper = babylocation+"T_tbarW_5f_powheg25ns.root";    ch[1]->Add(helper.c_str());

  dataset[2] = "T_tch_4f_amc25ns";
  ch[2] = new TChain("t");
  helper = babylocation+"T_tch_4f_amc25ns.root"; ch[2]->Add(helper.c_str());

  dataset[3] = "T_tch_5f_amc25ns";
  ch[3] = new TChain("t");
  helper = babylocation+"T_tch_5f_amc25ns.root";    ch[3]->Add(helper.c_str());

  dataset[4] = "WJets";
  ch[4] = new TChain("t");
  helper = babylocation+"WJetsToLNu_HT400To600_25ns.root"; ch[4]->Add(helper.c_str());
  helper = babylocation+"WJetsToLNu_HT100to200_25ns.root";     ch[4]->Add(helper.c_str());
  helper = babylocation+"WJetsToLNu_HT200To400_25ns.root";     ch[4]->Add(helper.c_str());
  helper = babylocation+"WJetsToLNu_HT600ToInf_25ns.root";     ch[4]->Add(helper.c_str());

  dataset[5] = "WW_25ns";
  ch[5] = new TChain("t");
  helper = babylocation+"WW_25ns.root";  ch[5]->Add(helper.c_str());

  dataset[6] = "WZ_25ns";
  ch[6] = new TChain("t");
  helper = babylocation+"WZ_25ns.root"; ch[6]->Add(helper.c_str());

  dataset[7] = "ZZ_25ns";
  ch[7] = new TChain("t");
  helper = babylocation+"ZZ_25ns.root"; ch[7]->Add(helper.c_str());

  ch[8] = new TChain("t");
  dataset[8] = "ttbar_madgraph_25ns";
  helper = babylocation+"ttbar_madgraph_25ns.root";  ch[8]->Add(helper.c_str());

  dataset[9] = "ttbar_powheg25ns";
  ch[9] = new TChain("t");
  helper = babylocation+"ttbar_powheg25ns.root";      ch[9]->Add(helper.c_str());

  dataset[10] = "SingleEl";
  ch[10] = new TChain("t");
//  helper = "/home/users/mliu/CMSSW_7_4_1_patch1/src/StopAnalysis/StopBabyMaker/data_single_muon_2015B_50ns.root";      ch[10]->Add(helper.c_str());
  helper = "/home/users/mliu/CMSSW_7_4_1_patch1/src/StopAnalysis/StopBabyMaker/data_single_electron_2015B_50ns.root";      ch[10]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
  //for(int i = 0; i<1; ++i){
   if(!(i==10)) continue; 
   //if(!(i==11)) continue; 
   ScanChain(ch[i],true,-1,dataset[i]); 
  }
//  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
}
