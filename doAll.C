{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("t"); 
  ch->Add("/nfs-7/userdata/mliu/ttz/ttg_25ns_test.root");
  ScanChain(ch); 
}