// Usage:
// > root -b doAll.C

// C++
#include <iostream>
#include <vector>
#include <map>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"

// CMS3
#include "StopBabies06232015.cc"

//MT2 variants
#include "stop_variables/Davismt2.h"
#include "stop_variables/topness.h"
#include "stop_variables/MT2_implementations.h"
#include "stop_variables/JetUtil.h"
#include "stop_variables/mt2w.h"
#include "stop_variables/mt2w_bisect.h"

using namespace std;
using namespace tas;

inline float getMT(LorentzVector lep,LorentzVector met){
  // From cmssw reco::deltaPhi()
  return TMath::Sqrt(2*met.Et()*lep.Et()*(1-TMath::Cos(JetUtil::deltaPhi(lep,met) ) ) );
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  map<string, TH1F*> histos; //massive
  vector<string> histonames; histonames.clear();
  vector<int> histobinn; histobinn.clear();
  vector<double> histobinl; histobinl.clear();
  vector<double> histobinu; histobinu.clear();
  map<string, float> value;

  histonames.push_back("MT2W");               histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("M3b");                histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("MTqq");               histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("Topness");            histobinn.push_back(30); histobinl.push_back(-15.); histobinu.push_back(15.);
  histonames.push_back("MT");                 histobinn.push_back(40); histobinl.push_back(0.); histobinu.push_back(400.);
  histonames.push_back("MET");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("HT");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("METoverHT");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1.25);
  histonames.push_back("METoverSqrtHT");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(40.);
  histonames.push_back("METoverSqrtHTplusleppt");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(40.);
  histonames.push_back("HTratio");            histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1.);
  histonames.push_back("dRLepBJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("dRLepBMax");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("dRLepBMin");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("dRLepFatJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("ak8_Over_ak10");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(3.);
  histonames.push_back("dRJet1FatJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("dRJet2FatJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("dRJet3FatJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("dRJet4FatJet");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("dR_ak8_ak10");          histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1.);
  histonames.push_back("dRbb");               histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("chi2");               histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(30.);
  histonames.push_back("NBJets");             histobinn.push_back( 5); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("NJets");              histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt40");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt50");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt60");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt70");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt80");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt90");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt100");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt125");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt150");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt175");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("NJets_pt200");         histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);
  histonames.push_back("minDPhi");            histobinn.push_back(32); histobinl.push_back(0.); histobinu.push_back(3.2);
  histonames.push_back("minDPhiJ3");          histobinn.push_back(32); histobinl.push_back(0.); histobinu.push_back(3.2);
  histonames.push_back("minDPhiB");           histobinn.push_back(32); histobinl.push_back(0.); histobinu.push_back(3.2);
  histonames.push_back("NBJetsOverNJets");    histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1.);
  histonames.push_back("pTl");                histobinn.push_back(27); histobinl.push_back(0.); histobinu.push_back(405.);
  histonames.push_back("pT_diff_jet1_FatJet"); histobinn.push_back(27); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("pTlb");               histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(600.);
  histonames.push_back("pTlbb");              histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(600.);
  histonames.push_back("pTleadj");            histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("pTleadb");            histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("pTtrailb");           histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("pTl4j");              histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("pTljets");            histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("pTl4joverMET");       histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);

  histonames.push_back("ak4_leadpt");         histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak4_subleadpt");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak4_leadmass");         histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak4_subleadmass");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);

  histonames.push_back("ak8_leadpt");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak8mass");         histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak8prunedmass");   histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak8trimmedmass");  histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak8filteredmass"); histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);

  histonames.push_back("ak10_leadpt");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak10_subleadpt");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("ak10prunedmass");   histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10mass");         histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10trimmedmass");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10filteredmass");     histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10genjets_mass");     histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10mass_Corrected");   histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);

  histonames.push_back("ak10mass_sublead");         histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10prunedmass_sublead");   histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10trimmedmass_sublead");      histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10filteredmass_sublead");     histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak10genjets_mass_sublead");     histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);

  histonames.push_back("ak12prunedmass");   histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak12trimmedmass");  histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("ak12filteredmass"); histobinn.push_back(20); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("DeltaPhiWl");       histobinn.push_back(32); histobinl.push_back(0.); histobinu.push_back(3.2);
  histonames.push_back("DeltaPhiSubFatJet_MET");histobinn.push_back(32); histobinl.push_back(0.); histobinu.push_back(3.2);

// book Cutflow hists
    map<string, TH1D*> histos_cutflow;
    vector<string> histonames_cutflow; histonames_cutflow.clear();
    histonames_cutflow.push_back("NEventsPerSignalRegion_nominal");
    histonames_cutflow.push_back("NEventsPerSignalRegion_nominal_mod");//nominal with jet pt requirement
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost1"); //only requirement on boosted object
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost2"); //both requirement
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost3"); //pt lead pt fatjet  > 250
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost4"); //pt lead pt fatjet  > 275
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost5"); //pt lead pt fatjet  > 325
    histonames_cutflow.push_back("NEventsPerSignalRegion_boost6"); //pt lead pt fatjet  > 350

  const unsigned int poststringsize =10 ;
  string poststring[poststringsize] = {"SR","SR250","SR250MT2W","SR_boost_200","SR_boost_250",
                                       "SR_boost_300","SR_boost_300_mass75","SR_boost_250_mass75","SR_boost_375_mass75","SR_boost_375_mass75_ak10only"};

  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix.find("TTbar")==std::string::npos&&b>0) continue;
    if(skimFilePrefix.find("TTbar")!=std::string::npos&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix.find("TTbar")!=std::string::npos&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix.find("TTbar")!=std::string::npos&&b==2) samplename = "TTbarH";
    for(unsigned int c = 0; c<poststringsize; ++c){
      for(unsigned int i = 0; i<histonames.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames[i] + "_"+samplename;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1F(mapname.c_str(), "", histobinn[i], histobinl[i], histobinu[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	if(c==0) value[histonames[i] ] = -1;
      }
    } 
   for(unsigned int i = 0;i<histonames_cutflow.size();i++){
    string histoname = histonames_cutflow.at(i)+"_"+samplename;
    histos_cutflow[histoname] = new TH1D(histoname.c_str(),"",5,0,5);
    histos_cutflow[histoname]->Sumw2(); histos_cutflow[histoname]->SetDirectory(rootdir);
    histos_cutflow[histoname]->GetYaxis()->SetTitle("Events");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(1,"MET>300,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(2,"MET>325,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(3,"MET>350,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(4,"MET>375,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(5,"MET>400,hDM");
    }
  }

  int c1vtx(0), c1l(0), cno2l(0), cno2track(0), cnotau(0), c2j(0), c1b(0), cmet(0);
  int cmt(0), cmdphi(0), cchi(0);
  // Loop over events to Analyze
  unsigned int nEventsTotal = 0;
  unsigned int nEventsChain = chain->GetEntries();
  if( nEvents >= 0 ) nEventsChain = nEvents;
  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;

  // File Loop
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    // Get File Content
    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    if(fast) TTreeCache::SetLearnEntries(10);
    if(fast) tree->SetCacheSize(128*1024*1024);
    cms3.Init(tree);
    
    // Loop over Events in current file
    if( nEventsTotal >= nEventsChain ) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for( unsigned int event = 0; event < nEventsTree; ++event) {
 
      // Get Event Content
      if( nEventsTotal >= nEventsChain ) continue;
      if(fast) tree->LoadTree(event);
      cms3.GetEntry(event);
      ++nEventsTotal;
   
      // Progress
      StopBabies06232015::progress( nEventsTotal, nEventsChain );

      string samplename = skimFilePrefix;
      if(skimFilePrefix.find("TTbar")!=std::string::npos){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
        else samplename = "TTbarH";
      }
      // Analysis Code
      //unsigned int run = cms3.run();
      //unsigned int ls = cms3.ls();
      //unsigned int evt = cms3.evt();
      float weight = cms3.scale1fb()*5.;
      int NLeps = cms3.ngoodleps();
      string ds = cms3.dataset();
      float MET = cms3.pfmet();
      float METPhi = cms3.pfmet_phi();
      float METx = MET*TMath::Cos(METPhi);
      float METy = MET*TMath::Sin(METPhi);
      float MT2W = cms3.MT2W();
      float MT = cms3.mt_met_lep();
      float dRLepBJet = cms3.dR_lep_leadb();
      float chi2 = cms3.hadronic_top_chi2();
      float HT = cms3.ak4_HT();
      float HTratio = cms3.ak4_htratiom();
      int nvtxs = cms3.nvtxs();
      float boostjet_ak8_pt = 0;
      float boostjet_ak10_pt = 0;
      float boostjet_ak10_pt_sublead = 0;
      float boostjet_ak12_pt = 0;
      if(cms3.ak8pfjets_p4().size())  boostjet_ak8_pt = cms3.ak8pfjets_p4()[0].pt();
      if(cms3.ak10pfjets_p4().size())  boostjet_ak10_pt = cms3.ak10pfjets_p4()[0].pt();
      if(cms3.ak10pfjets_p4().size()>1)  boostjet_ak10_pt_sublead = cms3.ak10pfjets_p4()[1].pt();
      if(cms3.ak12pfjets_p4().size())  boostjet_ak12_pt = cms3.ak12pfjets_p4()[0].pt();
      float minDPhi = cms3.mindphi_met_j1_j2();
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetslv = cms3.ak4pfjets_p4();
      vector<float> jetsbtag = cms3.ak4pfjets_CSV();
      vector<bool> jetsID = cms3.ak4pfjets_loose_pfid();
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep1lv = cms3.lep1_p4();
      float lep1pt = cms3.lep1_pt();
      float lep1eta = cms3.lep1_eta();
      float lep1dr03isoDB = cms3.lep1_relIso03DB();
      bool lep1eIDl = cms3.lep1_is_eleid_loose();
      bool lep1eIDm = cms3.lep1_is_eleid_medium();
      bool lep1mIDt = cms3.lep1_is_muoid_tight();
      bool lep1ismu = cms3.lep1_is_mu();
      bool lep1isel = cms3.lep1_is_el();
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep2lv = cms3.lep1_p4();
      float lep2pt = cms3.lep2_pt();
      float lep2eta = cms3.lep2_eta();
      float lep2dr03isoDB = cms3.lep2_relIso03DB();
      bool lep2eIDl = cms3.lep2_is_eleid_loose();
      bool lep2eIDm = cms3.lep2_is_eleid_medium();
      bool lep2mIDt = cms3.lep2_is_muoid_tight();
      bool lep2ismu = cms3.lep2_is_mu();
      bool lep2isel = cms3.lep2_is_el();
      
      bool trackveto = cms3.PassTrackVeto();
      bool tauveto = cms3.PassTauVeto();

      int NGLeps = 0;
      int NSLeps = 0;
      int NGJets = 0;
      int NGBJets = 0;

      int NGJets_pt40 = 0;
      int NGJets_pt50 = 0;
      int NGJets_pt60 = 0;
      int NGJets_pt70 = 0;
      int NGJets_pt80 = 0;
      int NGJets_pt90 = 0;
      int NGJets_pt100 = 0;
      int NGJets_pt125 = 0;
      int NGJets_pt150 = 0;
      int NGJets_pt175 = 0;
      int NGJets_pt200 = 0;

      int l1=-1;
      if(lep1ismu){
	if(lep1pt>30&&fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1isel){
	if(lep1pt>40&&fabs(lep1eta)<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      /*
      if(lep2ismu){
	//if(lep2pt>20&&fabs(lep2eta)<99&&lep2mIDt) ++NGLeps;
	//if(lep2pt>25&&fabs(lep2eta)<2.1&&lep2mIDt&&lep2dr03isoDB*lep2pt<TMath::Min(5.,0.15*lep2pt)) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
	if(lep2pt>30&&fabs(lep2eta)<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; if(l1!=1) l1 = 2; else l1=-2; }
      } else if (lep2isel){
	//if(lep2pt>20&&fabs(lep2eta)<99&&lep2eIDl) ++NGLeps;
	//if(lep2pt>30&&fabs(lep2eta)<1.442&&lep2eIDm&&lep2dr03isoDB*lep2pt<TMath::Min(5.,0.15*lep2pt)) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
	if(lep2pt>40&&fabs(lep2eta)<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
      }
      */
      NGLeps = NLeps;

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) leplv = lep1lv;
      else if(l1==2) leplv = lep2lv;

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      metlv.SetPxPyPzE(METx,METy,0.,MET);

      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > boostjetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > bjetlv;
      vector<float> btag;

      int jj1(-1), jj2(-1), jj3(-1);// j1(-1), j2(-1), jm(-1);
      float jj1d(-1), jj2d(-1), jj3d(-1);// jjmm(-1);
      for(unsigned int nj = 0; nj<jetsbtag.size(); ++nj){
	if(jetslv[nj].Pt()<30) continue;
	if(fabs(jetslv[nj].Eta())>2.4) continue;
	if(jetsID[nj]==false) continue;
	jind.push_back(nj);
	jetlv.push_back(jetslv[nj]);
	if(jetslv[nj].Pt()>250) boostjetlv.push_back(jetslv[nj]);
	btag.push_back(jetsbtag[nj]);
	++NGJets;

        if(jetslv[nj].Pt()>40) ++NGJets_pt40;
        if(jetslv[nj].Pt()>50) ++NGJets_pt50;
        if(jetslv[nj].Pt()>60) ++NGJets_pt60;
        if(jetslv[nj].Pt()>70) ++NGJets_pt70;
        if(jetslv[nj].Pt()>80) ++NGJets_pt80;
        if(jetslv[nj].Pt()>90) ++NGJets_pt90;
        if(jetslv[nj].Pt()>100) ++NGJets_pt100;
        if(jetslv[nj].Pt()>125) ++NGJets_pt125;
        if(jetslv[nj].Pt()>150) ++NGJets_pt150;
        if(jetslv[nj].Pt()>175) ++NGJets_pt175;
        if(jetslv[nj].Pt()>200) ++NGJets_pt200;

	if(jetsbtag[nj]>0.814) {++NGBJets; bind.push_back(nj); bjetlv.push_back(jetslv[nj]);}
	if(jetslv[nj].Pt()>jj1d){
	  jj3d = jj2d; jj2d = jj1d; jj1d = jetslv[nj].Pt();
	  jj3 = jj2; jj2 = jj1; jj1 = nj;
	} else if(jetslv[nj].Pt()>jj2d){
	  jj3d = jj2d; jj2d = jetslv[nj].Pt();
	  jj3 = jj2; jj2 = nj;
	} else if(jetslv[nj].Pt()>jj3d){
	  jj3d = jetslv[nj].Pt();
	  jj3 = nj;
	}
      }

      if(nvtxs<0) continue; ++c1vtx;
      if(NGLeps!=1) continue; ++c1l;
      if(NSLeps!=1) continue; ++cno2l;
      if(!trackveto) continue; ++cno2track;
      if(!tauveto) continue; ++cnotau;
      if(NGJets<3) continue; ++c2j;
      if(NGBJets<1) continue; ++c1b;
      if(MET<200) continue; ++cmet;
      if(MT<150) continue;
//      if(!boostjet_ak10_pt) continue;

      bool SR = false; bool SR250(false),SR250MT2W(false),SR_boost_200(false),SR_boost_250(false),SR_boost_300(false),SR_boost_300_mass75(false),SR_boost_250_mass75(false),SR_boost_375_mass75(false),SR_boost_375_mass75_ak10only(false);

      //if(MT>150 && MET>250 && MT2W>-1 && boostjet_ak10_pt>250) SR = true;
      if(MT>150 && MET>300 && MT2W>-1 ) SR = true;
      if(MT>150 && MET>300 && MT2W>-1 ) SR250 = true;
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>200) SR_boost_200 = true;
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>250) SR_boost_250 = true;
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300) SR_boost_300 = true;
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>250 && cms3.ak10pfjets_trimmed_mass()[0]>75 && cms3.ak8pfjets_trimmed_mass()[0]>50) SR_boost_250_mass75 = true;
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }

      vector<LorentzVector > btaggedjets = JetUtil::BJetSelector(jetlv,btag,0.814,2,2,2);

      vector<LorentzVector > dummybjets;dummybjets.clear();
      LorentzVector dummybjet; dummybjet.SetPxPyPzE(0.,0.,0.,0.); dummybjets.push_back(dummybjet);dummybjets.push_back(dummybjet);
      if(MT>150 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75) SR_boost_300_mass75 = true;
      //if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75 && Gettopness_(MET,METPhi,leplv,btaggedjets,1)>7) SR_boost_300_mass75 = true;
//      if(MT>150 && MET>300 && MT2W>200 && minDPhi>0.8) SR250MT2W = true;
      if(MT>150 && MET>300 && MT2W>200 && chi2<10 && minDPhi>0.8) SR250MT2W = true;
      if(MT>150 && MET>375 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75 && cms3.ak8pfjets_trimmed_mass()[0]>75) SR_boost_375_mass75 = true;
      //if(MT>150 && MET>375 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75&&Gettopness_(MET,METPhi,leplv,btaggedjets,1)>5 && cms3.ak8pfjets_trimmed_mass()[0]>75) SR_boost_375_mass75 = true;
      if(MT>150 && MET>375 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75&&Gettopness_(MET,METPhi,leplv,btaggedjets,1)>5 ) SR_boost_375_mass75_ak10only = true;
      
      value["MT2W"] = MT2W;
      //value["Topness"] = Gettopness_(MET,METPhi,leplv,btaggedjets,1);
      value["Topness"] = topnessMod();
      value["ak8_leadpt"] = boostjet_ak8_pt;
      value["ak10_leadpt"] = boostjet_ak10_pt;
      value["ak10_subleadpt"] = boostjet_ak10_pt_sublead;
      value["ak4_leadpt"] = jetlv[0].pt();
      value["ak4_subleadpt"] = jetlv[1].pt();
      value["ak4_leadmass"] = jetlv[0].mass();
      value["ak4_subleadmass"] = jetlv[1].mass();

      if(ak10genjets_p4().size()) value["ak10genjets_mass"] = ak10genjets_p4()[0].mass(); 

      int leadb(-1), trailb(-1);
      for(unsigned int n = 0; n<btaggedjets.size();++n){
	if(leadb<0) leadb = n;
	else if(trailb<0) trailb = n;
	else if(btaggedjets[n].Pt() > btaggedjets[leadb].Pt()){ trailb = leadb; leadb = n; }
	else if(btaggedjets[n].Pt() > btaggedjets[trailb].Pt()){ trailb = n; }
	LorentzVector temp = btaggedjets[n]+leplv;
	if(temp.Pt()<value["pTlb"]) value["pTlb"] = temp.Pt();  else if(value["pTlb"]<0) value["pTlb"] = temp.Pt();
	float tmp = getMT(btaggedjets[n],metlv);
	tmp = JetUtil::deltaR(btaggedjets[n],leplv);
	if(tmp<value["dRLepBMin"]) value["dRLepBMin"] = tmp;  else if(value["dRLepBMin"]<0) value["dRLepBMin"] = tmp;
	if(tmp>value["dRLepBMax"]) value["dRLepBMax"] = tmp;  else if(value["dRLepBMax"]<0) value["dRLepBMax"] = tmp;
	tmp = JetUtil::deltaPhi(btaggedjets[n],metlv);
	if(tmp<value["minDPhiB"]) value["minDPhiB"] = tmp;  else if(value["minDPhiB"]<0) value["minDPhiB"] = tmp;
	for(unsigned int m = n+1; m<btaggedjets.size();++m){
	  temp = btaggedjets[n]+btaggedjets[m]+leplv;
	  if(temp.Pt()<value["pTlbb"]) value["pTlbb"] = temp.Pt();  else if(value["pTlbb"]<0) value["pTlbb"] = temp.Pt();
	  tmp = JetUtil::deltaR(btaggedjets[n],btaggedjets[m]);
	  if(tmp<value["dRbb"]) value["dRbb"] = tmp;  else if(value["dRbb"]<0) value["dRbb"] = tmp;
	}
      }
      float myleadjpt = -1.;
      float myleadjpt300 = -1.;
      float myleadjptno300 = -1.;
      int bj1(-1), bj2(-1),bj3(-1);
      for(unsigned int n = 0; n<jetlv.size();++n){
	float tmp = getMT(jetlv[n],metlv);
	for(unsigned int m = n+1; m<jetlv.size();++m){
	  tmp = getMT(jetlv[n]+jetlv[m],metlv);
	  if(tmp<value["MTqq"]) value["MTqq"] = tmp;  else if(value["MTqq"]<0) value["MTqq"] = tmp;
	}
	if(n<3){
	  tmp = JetUtil::deltaPhi(jetlv[n],metlv);
	  if(tmp<value["minDPhiJ3"]) value["minDPhiJ3"] = tmp;  else if(value["minDPhiJ3"]<0) value["minDPhiJ3"] = tmp;
	}
	float dP1l(-1), dP2l(-1), dP3l(-1);
	if(bj1>0) dP1l = JetUtil::deltaPhi(jetlv[bj1],leplv);
	if(bj2>0) dP2l = JetUtil::deltaPhi(jetlv[bj2],leplv);
	if(bj3>0) dP3l = JetUtil::deltaPhi(jetlv[bj3],leplv);
	if(JetUtil::deltaPhi(jetlv[n],leplv)>dP1l){ bj3 = bj2; bj2 = bj1; bj1 = n;}
	else if(JetUtil::deltaPhi(jetlv[n],leplv)>dP2l){ bj3 = bj2; bj2 = n;}
	else if(JetUtil::deltaPhi(jetlv[n],leplv)>dP3l){ bj3 = n;}
      }
      
      LorentzVector bjsumlep = jetlv[bj1] + jetlv[bj2] + jetlv[bj3];
      value["M3b"] = bjsumlep.M();


      if(cms3.ak8pfjets_p4().size()) value["DeltaPhiSubFatJet_MET"]   = JetUtil::deltaPhi(cms3.ak8pfjets_p4()[0],metlv);
      if(cms3.ak8pfjets_p4().size()&&cms3.ak10pfjets_p4().size()) value["ak8_Over_ak10"]   = cms3.ak8pfjets_trimmed_mass()[0]/cms3.ak10pfjets_trimmed_mass()[0];
      if(cms3.ak10pfjets_p4().size()) value["dRLepFatJet"]   = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],leplv);
      if(cms3.ak10pfjets_p4().size()) value["dRJet1FatJet"]   = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],jetlv[0]);
      if(cms3.ak10pfjets_p4().size()) value["pT_diff_jet1_FatJet"]   = fabs(cms3.ak10pfjets_p4()[0].pt()-jetlv[0].pt());
      if(cms3.ak10pfjets_p4().size()) value["dRJet2FatJet"]   = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],jetlv[1]);
      if(cms3.ak10pfjets_p4().size()) value["dRJet3FatJet"]   = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],jetlv[2]);
      if(cms3.ak10pfjets_p4().size()) value["dRJet4FatJet"]   = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],jetlv[3]);
      if(cms3.ak8pfjets_p4().size()&&cms3.ak10pfjets_p4().size()) value["dR_ak8_ak10"] = JetUtil::deltaR(cms3.ak10pfjets_p4()[0],cms3.ak8pfjets_p4()[0]); 
      if(cms3.ak8pfjets_p4().size()) value["ak8mass"]   = cms3.ak8pfjets_p4()[0].mass();
      if(cms3.ak8pfjets_trimmed_mass().size()) value["ak8prunedmass"]   = cms3.ak8pfjets_pruned_mass()[0];
      if(cms3.ak8pfjets_pruned_mass().size()) value["ak8trimmedmass"]  = cms3.ak8pfjets_trimmed_mass()[0];
      if(cms3.ak8pfjets_filtered_mass().size())value["ak8filteredmass"] = cms3.ak8pfjets_filtered_mass()[0];

      if(cms3.ak10pfjets_p4().size()) value["ak10mass"]   = cms3.ak10pfjets_p4()[0].mass();
      if(cms3.ak10pfjets_trimmed_mass().size()) value["ak10prunedmass"]   = cms3.ak10pfjets_pruned_mass()[0];
      if(cms3.ak10pfjets_pruned_mass().size()) value["ak10trimmedmass"]  = cms3.ak10pfjets_trimmed_mass()[0];
      if(cms3.ak10pfjets_filtered_mass().size())value["ak10filteredmass"] = cms3.ak10pfjets_filtered_mass()[0];

      if(cms3.ak10pfjets_p4().size()>1) value["ak10mass_sublead"]   = cms3.ak10pfjets_p4()[1].mass();
      if(cms3.ak10pfjets_trimmed_mass().size()>1) value["ak10prunedmass_sublead"]   = cms3.ak10pfjets_pruned_mass()[1];
      if(cms3.ak10pfjets_pruned_mass().size()>1) value["ak10trimmedmass_sublead"]  = cms3.ak10pfjets_trimmed_mass()[1];
      if(cms3.ak10pfjets_filtered_mass().size()>1)value["ak10filteredmass_sublead"] = cms3.ak10pfjets_filtered_mass()[1];

      if(cms3.ak12pfjets_trimmed_mass().size()) value["ak12prunedmass"]   = cms3.ak12pfjets_pruned_mass()[0];
      if(cms3.ak12pfjets_pruned_mass().size()) value["ak12trimmedmass"]  = cms3.ak12pfjets_trimmed_mass()[0];
      if(cms3.ak12pfjets_filtered_mass().size())value["ak12filteredmass"] = cms3.ak12pfjets_filtered_mass()[0];

      value["MT"] = MT;
      value["MET"] = MET;
      value["HT"] = HT;
      if(HT>0) {
	value["METoverHT"] = MET/HT;
	value["METoverSqrtHT"] = MET/TMath::Sqrt(HT);
	value["METoverSqrtHTplusleppt"] = MET/TMath::Sqrt(HT+leplv.Pt());
      }

      value["HTratio"] = HTratio;
      value["dRLepBJet"] = dRLepBJet;
      value["chi2"] = chi2;
      value["NBJets"] = NGBJets;
      value["NJets"] = NGJets;
      value["NJets_pt40"] = NGJets_pt40;
      value["NJets_pt50"] = NGJets_pt50;
      value["NJets_pt60"] = NGJets_pt60;
      value["NJets_pt70"] = NGJets_pt70;
      value["NJets_pt80"] = NGJets_pt80;
      value["NJets_pt90"] = NGJets_pt90;
      value["NJets_pt100"] = NGJets_pt100;
      value["NJets_pt125"] = NGJets_pt125;
      value["NJets_pt150"] = NGJets_pt150;
      value["NJets_pt175"] = NGJets_pt175;
      value["NJets_pt200"] = NGJets_pt200;
      value["NBJetsOverNJets"] = float(NGBJets)/float(NGJets);
      value["minDPhi"] = minDPhi;
      value["pTl"] = leplv.Pt();
      value["pTl4j"] = (leplv+jetslv[0]+jetslv[1]+jetslv[2]+jetslv[2]).Pt();
      value["pTl4joverMET"] = (leplv+jetslv[0]+jetslv[1]+jetslv[2]+jetslv[2]).Pt()/metlv.Pt();
      value["pTleadj"] = myleadjpt;
      value["pTleadb"] = btaggedjets[leadb].Pt();
      value["pTtrailb"] = btaggedjets[trailb].Pt();
      value["DeltaPhiWl"] = JetUtil::deltaPhi(leplv,metlv+leplv);

      for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
	if(SR){
	  mname = "SR"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR250){
	  mname = "SR250"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR250MT2W){
	  mname = "SR250MT2W"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_200){
	  mname = "SR_boost_200"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_250){
	  mname = "SR_boost_250"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_300){
	  mname = "SR_boost_300"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_300_mass75){
	  mname = "SR_boost_300_mass75"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_375_mass75){
	  mname = "SR_boost_375_mass75"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_375_mass75_ak10only){
	  mname = "SR_boost_375_mass75_ak10only"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR_boost_250_mass75){
	  mname = "SR_boost_250_mass75"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
      }

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
      //boosted regions
//      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && fabs(ak10pfjets_p4()[0].eta()) < 2.4 &&boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75&&Gettopness_(MET,METPhi,leplv,btaggedjets,1)>7){
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && boostjet_ak10_pt>300  && fabs(ak10pfjets_p4()[0].eta()) < 2.4 && cms3.ak10pfjets_trimmed_mass()[0]>75){
      if(MET>300&&MET<350) histos_cutflow["NEventsPerSignalRegion_boost1_"+samplename]  ->Fill(0.5,weight);
      if(MET>350&&MET<400) histos_cutflow["NEventsPerSignalRegion_boost1_"+samplename]  ->Fill(1.5,weight);
      if(MET>400&&MET<450) histos_cutflow["NEventsPerSignalRegion_boost1_"+samplename]  ->Fill(2.5,weight);
      if(MET>450&&MET<500) histos_cutflow["NEventsPerSignalRegion_boost1_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost1_"+samplename]  ->Fill(4.5,weight);
      }
      //if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && fabs(ak10pfjets_p4()[0].eta()) < 2.4 &&boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75&&Gettopness_(MET,METPhi,leplv,btaggedjets,1)>7){
      if(MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && fabs(ak10pfjets_p4()[0].eta()) < 2.4 &&boostjet_ak10_pt>300 && cms3.ak10pfjets_trimmed_mass()[0]>75&&  topnessMod() > 0 ){
//      if(SR_boost_300_mass75 && jetlv[0].pt()>70 && jetlv[1].pt()>70&&jetlv[2].pt()>50){
      if(MET>300&&MET<350) histos_cutflow["NEventsPerSignalRegion_boost2_"+samplename]  ->Fill(0.5,weight);
      if(MET>350&&MET<400) histos_cutflow["NEventsPerSignalRegion_boost2_"+samplename]  ->Fill(1.5,weight);
      if(MET>400&&MET<450) histos_cutflow["NEventsPerSignalRegion_boost2_"+samplename]  ->Fill(2.5,weight);
      if(MET>450&&MET<500) histos_cutflow["NEventsPerSignalRegion_boost2_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost2_"+samplename]  ->Fill(4.5,weight);
      }    
      //if(boostjet_ak10_pt>250 && MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && cms3.ak10pfjets_trimmed_mass()[0]>75 && jetlv[0].pt()>70 && jetlv[1].pt()>70&&jetlv[2].pt()>50){
      if(boostjet_ak10_pt>250 && MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && cms3.ak10pfjets_trimmed_mass()[0]>75){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_boost3_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_boost3_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_boost3_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_boost3_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost3_"+samplename]  ->Fill(4.5,weight);
      }
      if(boostjet_ak10_pt>275 && MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && cms3.ak10pfjets_trimmed_mass()[0]>75){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_boost4_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_boost4_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_boost4_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_boost4_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost4_"+samplename]  ->Fill(4.5,weight);
      }
      if(boostjet_ak10_pt>325 && MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && cms3.ak10pfjets_trimmed_mass()[0]>75){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_boost5_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_boost5_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_boost5_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_boost5_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost5_"+samplename]  ->Fill(4.5,weight);
      }
      if(boostjet_ak10_pt>350 && MT>150 && MET>300 && minDPhi>0.8 && MT2W>200 && cms3.ak10pfjets_trimmed_mass()[0]>75){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_boost6_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_boost6_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_boost6_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_boost6_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_boost6_"+samplename]  ->Fill(4.5,weight);
      }
      //nominal
      //if(SR250MT2W && NGJets>=4 && Gettopness_(MET,METPhi,leplv,btaggedjets,1)>7){
      if(SR250MT2W && NGJets>=4){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(4.5,weight);
      //modified signal region
      //if(SR250MT2W && NGJets>=3 && jetlv[0].pt()>150 && jetlv[1].pt()>80 && jetlv[2].pt()>50){
      if(SR250MT2W && NGJets>=3 && topnessMod() > 0){
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_nominal_mod_"+samplename]  ->Fill(0.5,weight);
      if(MET>350) histos_cutflow["NEventsPerSignalRegion_nominal_mod_"+samplename]  ->Fill(1.5,weight);
      if(MET>400) histos_cutflow["NEventsPerSignalRegion_nominal_mod_"+samplename]  ->Fill(2.5,weight);
      if(MET>450) histos_cutflow["NEventsPerSignalRegion_nominal_mod_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_nominal_mod_"+samplename]  ->Fill(4.5,weight);
      }
     }
      if(MT<120) continue; ++cmt;
      if(minDPhi<0.8) continue; ++cmdphi;
      if(chi2>5) continue; ++cchi;
    }
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  
  // Example Histograms
  // samplehisto->Draw();

  for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
  }

  string filename = "rootfiles/CutHistos/FatJetPlots/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
  for(map<string,TH1D*>::iterator h=    histos_cutflow.begin(); h!=    histos_cutflow.end();++h) h->second->Write();
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;
  // return
  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed" << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:	" << Form( "%.01f", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:	" << Form( "%.01f", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;
  return 0;
}
