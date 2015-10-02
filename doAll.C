{

  gROOT->ProcessLine(".L ScanChain.C+");

  TChain *ch = new TChain("t"); 
  ch->Add("/nfs-7/userdata/mliu/ttz/scale_unc/ttz_zllnunu_25ns.root");
  ScanChain(ch); 
}