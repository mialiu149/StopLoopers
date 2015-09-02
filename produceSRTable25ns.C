{
  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");

  gROOT->ProcessLine(".L GetCutFlow.C+");
  const unsigned int chainsize = 9;
  TChain *ch[chainsize];
  string dataset[chainsize];
//  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/";
  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V728_met30_ge1l_ge2j_Phys14/";
//  string babylocation = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_50ns/";
  string helper;
  dataset[0] = "TTbar";
  ch[0] = new TChain("t"); 
  //helper = babylocation+"ttbar.root"; ch[0]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/ttbar_powheg25ns.root";ch[0]->Add(helper.c_str());
  dataset[1] = "SingleT";
  ch[1] = new TChain("t");
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/T_sch_4f_amc25ns.root";ch[1]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/T_tbarW_5f_powheg25ns.root";ch[1]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/T_tch_4f_amc25ns.root";ch[1]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/T_tch_5f_amc25ns.root";ch[1]->Add(helper.c_str());
  /*dataset[2] = "Rare";
  ch[2] = new TChain("t");
  helper = babylocation+"ttwjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"ttzjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"wzjets.root";  ch[2]->Add(helper.c_str());
  helper = babylocation+"zz.root";      ch[2]->Add(helper.c_str());
  */
  dataset[2] = "TTV";
  ch[2] = new TChain("t");
  helper = babylocation+"ttwjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"ttzjets.root"; ch[2]->Add(helper.c_str());
  ch[8] = new TChain("t");
  dataset[8] = "VV";
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/WW_25ns.root";ch[8]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/WZ_25ns.root";ch[8]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V74x_met30_ge1l_ge2j_25ns/ZZ_25ns.root";ch[8]->Add(helper.c_str());

  dataset[3] = "WJets";
  ch[3] = new TChain("t"); 
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1l_ge2j_tempsync_20150618/wjetsHT100.root";      ch[3]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1l_ge2j_tempsync_20150618/wjetsHT200.root";      ch[3]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1l_ge2j_tempsync_20150618/wjetsHT400.root";      ch[3]->Add(helper.c_str());
  helper = "/nfs-7/userdata/stopRun2/StopBabies_V07_02_08_met30_ge1l_ge2j_tempsync_20150618/wjetsHT600.root";      ch[3]->Add(helper.c_str());
  dataset[4] = "Stop_425_325";
  ch[4] = new TChain("t"); 
  helper = babylocation+"stop_425_325.root";      ch[4]->Add(helper.c_str());
  dataset[5] = "Stop_500_325";
  ch[5] = new TChain("t"); 
  helper = babylocation+"stop_500_325.root";      ch[5]->Add(helper.c_str());
  dataset[6] = "Stop_650_325";
  ch[6] = new TChain("t"); 
  helper = babylocation+"stop_650_325.root";      ch[6]->Add(helper.c_str());
  dataset[7] = "Stop_850_100";
  ch[7] = new TChain("t"); 
  helper = babylocation+"stop_850_100.root";      ch[7]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
   if(i==1) continue; 
   ScanChain(ch[i],true,-1,dataset[i]); 
  }
//  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
 // gROOT->ProcessLine(".x makePlotPlot.C+");
  //gROOT->ProcessLine(".x makePlotCutFlow.C+");
  gROOT->ProcessLine(".x GetSRTable.C+");
  //gROOT->ProcessLine(".x GetCFNumbers.C+");
  //gROOT->ProcessLine(".x makePlotIP.C+");
}
