{
  gSystem->Load("libMinuit2.so");
  gROOT->ProcessLine(".L stop_variables/JetUtil.cc+");
  gROOT->ProcessLine(".L stop_variables/topness.cc+");
  gROOT->ProcessLine(".L stop_variables/Davismt2.cc+");
  gROOT->ProcessLine(".L stop_variables/MT2_implementations.cc+");
  gROOT->ProcessLine(".L stop_variables/mt2w_bisect.cpp+");
  gROOT->ProcessLine(".L stop_variables/mt2w.cc+");
  gROOT->ProcessLine(".L GetHistoFatJets.C+");
  
  const unsigned int chainsize = 7;
  TChain *ch[chainsize];
  string dataset[chainsize];
  string babylocation = "/nfs-7/userdata/mliu/fatjet/";
  string helper;
  dataset[0] = "TTbar_fatjet";
  ch[0] = new TChain("t"); 
  helper = babylocation+"ttbar_fatjet.root"; ch[0]->Add(helper.c_str());
/*  dataset[1] = "SingleT";
  ch[1] = new TChain("t");
  helper = babylocation+"t_sch.root";    ch[1]->Add(helper.c_str());
  helper = babylocation+"tbar_sch.root"; ch[1]->Add(helper.c_str());
  helper = babylocation+"t_tch.root";    ch[1]->Add(helper.c_str());
  helper = babylocation+"tbar_tch.root"; ch[1]->Add(helper.c_str());
  helper = babylocation+"t_tW.root";     ch[1]->Add(helper.c_str());
  helper = babylocation+"tbar_tW.root";  ch[1]->Add(helper.c_str());
*/ 
 /*dataset[2] = "Rare";
  ch[2] = new TChain("t");
  helper = babylocation+"ttwjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"ttzjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"wzjets.root";  ch[2]->Add(helper.c_str());
  helper = babylocation+"zz.root";      ch[2]->Add(helper.c_str());
  */
  /*dataset[2] = "TTV";
  ch[2] = new TChain("t");
  helper = babylocation+"ttwjets.root"; ch[2]->Add(helper.c_str());
  helper = babylocation+"ttzjets.root"; ch[2]->Add(helper.c_str());
  ch[8] = new TChain("t");
  dataset[8] = "VV";
  helper = babylocation+"wzjets.root";  ch[8]->Add(helper.c_str());
  helper = babylocation+"zz.root";      ch[8]->Add(helper.c_str());

  dataset[3] = "WJets";
  ch[3] = new TChain("t"); 
  helper = babylocation+"wjets.root";      ch[3]->Add(helper.c_str());
  dataset[4] = "Stop_425_325";
  ch[4] = new TChain("t"); 
  helper = babylocation+"stop_425_325.root";      ch[4]->Add(helper.c_str());
  dataset[5] = "Stop_500_325";
  ch[5] = new TChain("t"); 
  helper = babylocation+"stop_500_325.root";      ch[5]->Add(helper.c_str());
  dataset[6] = "Stop_650_325";
  ch[6] = new TChain("t"); 
  helper = babylocation+"stop_650_325.root";      ch[6]->Add(helper.c_str());
 */
  dataset[1] = "Stop_850_100_fatjet";
  ch[1] = new TChain("t"); 
  helper = babylocation+"stop_850_100_fatjet.root";      ch[1]->Add(helper.c_str());

  dataset[2] = "ttw_fatjet";
  ch[2] = new TChain("t"); 
  helper = babylocation+"ttw_fatjet.root";      ch[2]->Add(helper.c_str());

  dataset[3] = "ttz_fatjet";
  ch[3] = new TChain("t"); 
  helper = babylocation+"ttz_fatjet.root";      ch[3]->Add(helper.c_str());

  dataset[4] = "wjets_fatjet";
  ch[4] = new TChain("t"); 
  helper = babylocation+"wjets_fatjet_100.root";      ch[4]->Add(helper.c_str());
  helper = babylocation+"wjets_fatjet_200.root";      ch[4]->Add(helper.c_str());
  helper = babylocation+"wjets_fatjet_400.root";      ch[4]->Add(helper.c_str());
  helper = babylocation+"wjets_fatjet_600.root";      ch[4]->Add(helper.c_str());

  dataset[5] = "sch_tch_fatjet";
  ch[5] = new TChain("t"); 
  helper = babylocation+"tbar_s_fatjet.root";       ch[5]->Add(helper.c_str());
  helper = babylocation+"t_s_fatjet.root";          ch[5]->Add(helper.c_str());
  helper = babylocation+"tbar_t_fatjet.root";       ch[5]->Add(helper.c_str());
  helper = babylocation+"t_t_fatjet.root";          ch[5]->Add(helper.c_str());

  dataset[6] = "tW_fatjet";
  ch[6] = new TChain("t"); 
  helper = babylocation+"t_tw_fatjet.root";         ch[6]->Add(helper.c_str());
  helper = babylocation+"tbar_tw_fatjet.root";      ch[6]->Add(helper.c_str());

  for(int i = 0; i<chainsize; ++i){
    //if(i<1) continue;
    ScanChain(ch[i],true,-1,dataset[i]); 
  }
  gROOT->ProcessLine(".L dataMCplotMaker.cc+");
  gROOT->ProcessLine(".x makePlotFatJet.C+");
  //gROOT->ProcessLine(".x makePlotCutFlow.C+");
  gROOT->ProcessLine(".x GetSRTableFatJet.C+");
  //gROOT->ProcessLine(".x GetCFNumbers.C+");
  //gROOT->ProcessLine(".x makePlotIP.C+");
}
