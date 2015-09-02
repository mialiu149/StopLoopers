// Usage:
// > root -b doAll.C

// C++
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <iomanip>
#include <algorithm>

// ROOT
#include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"

// StopBabies20150603
//#include "StopBabies20150603_old20150505.cc"
#include "StopBabies20150603.cc"

//MT2 variants
#include "stop_variables/Davismt2.h"
#include "stop_variables/topness.h"
#include "stop_variables/MT2_implementations.h"
#include "stop_variables/JetUtil.h"
#include "stop_variables/mt2w.h"
#include "stop_variables/mt2w_bisect.h"

using namespace std;
using namespace tas;

inline double DeltaPhi(double phi1, double phi2){
  // From cmssw reco::deltaPhi()
  double result = phi1 - phi2;
  while( result >   TMath::Pi() ) result -= TMath::TwoPi();
  while( result <= -TMath::Pi() ) result += TMath::TwoPi();
  return TMath::Abs(result);
}

struct PrintOut{
  int Run;
  int LS;
  unsigned int Evt;
  float leppt;
  int leppdgid;
  float met;
  float mt;
  int njets;
  int nbjets;
  int ngenleps;
  float DPhiWlep;
  float htssm;
  float Mlb_cb;
  float Mjjj;
  float mt2w;
  float Topness;
  float chi2;
};
bool SortPrintOut(PrintOut const& lhs, PrintOut const& rhs) {
    if (lhs.Run != rhs.Run)
        return lhs.Run < rhs.Run;
    if (lhs.LS != rhs.LS)
        return lhs.LS < rhs.LS;
    return lhs.Evt < rhs.Evt;
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  /*
  map<string, TH1F*> histos; //massive
  vector<string> histonames; histonames.clear();
  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix!="TTbar"&&b>0) continue;
    if(skimFilePrefix=="TTbar"&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix=="TTbar"&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix=="TTbar"&&b==2) samplename = "TTbarH";
    for(unsigned int i = 0; i<histonames.size(); ++i){
      string mapname;
      mapname = histonames[i] + "_"+samplename;
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1F(mapname.c_str(), "", 50, 0, 750);
      histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
    }
  }
  */

  vector <PrintOut> POel; POel.clear();
  vector <PrintOut> POmu; POmu.clear();
  bool fWriteToFile = true;
  bool fAppend = false;
  TString outputdir = "/home/users/mliu/CMSSW_7_2_0/src/stop_2015/AnalysisCode_newbabyformat/";
  std::ostringstream* fLogStreamEl     = 0;
  fLogStreamEl = new std::ostringstream();
  std::ostringstream* fLogStreamMu     = 0;
  fLogStreamMu = new std::ostringstream();

  
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
      StopBabies20150603::progress( nEventsTotal, nEventsChain );

      string samplename = skimFilePrefix;
      if(skimFilePrefix=="TTbar"){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
	else samplename = "TTbarH";
	//if(cms3.gen_nmus_() + cms3.gen_nels_() + cms3.gen_ntaus_() ==2) samplename = "TTbar2l";
	//else if(cms3.gen_nmus_() + cms3.gen_nels_() + cms3.gen_ntaus_() ==1) samplename = "TTbar1l";
	//else samplename = "TTbarH";
      }

      // Analysis Code
      bool print = false;

      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetslv = cms3.ak4pfjets_p4();
      vector<float> jetsbtag = cms3.ak4pfjets_CSV();
      vector<bool> jetsID = cms3.ak4pfjets_loose_pfid();


      int NGLeps = 0;
      int NSLeps = 0;
      int NGJets = 0;
      int NGBJets = 0;

      int l1=-1; bool lele = false; bool lmuo = false; int lpdg = -1; float lpt = -1;
     /* if(lep1_is_mu()){
	if(lep1_pt()>30&&fabs(lep1_eta())<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1_is_el()){
	if(lep1_pt()>40&&fabs(lep1_eta())<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      if(lep2_is_mu()){
	if(lep2_pt()>30&&fabs(lep2_eta())<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; if(l1!=1) l1 = 2; else l1=-2; }
      } else if (lep2_is_el()){
	if(lep2_pt()>40&&fabs(lep2_eta())<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
      }*/
      bool lep1ismu = cms3.lep1_is_mu();
      bool lep1isel = cms3.lep1_is_el();
      bool lep2ismu = cms3.lep2_is_mu();
      bool lep2isel = cms3.lep2_is_el();
      float lep2pt = cms3.lep2_pt();
      float lep2eta = cms3.lep2_eta();
      float lep1pt = cms3.lep1_pt();
      float lep1eta = cms3.lep1_eta();
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep1lv = cms3.lep1_p4();
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep2lv = cms3.lep2_p4();
      if(lep1ismu){
	if(lep1pt>30&&fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1isel){
	if(lep1pt>40&&fabs(lep1eta)<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      if(lep2ismu){
	if(lep2pt>30&&fabs(lep2eta)<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 2; }
      } else if (lep2isel){
	if(lep2pt>40&&fabs(lep2eta)<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 2; }
      }
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(NSLeps==1) {if(l1==1) leplv = lep1lv;
      else if(l1==2) leplv = lep2lv;}
      if(NSLeps>1&&l1>0) leplv = lep1lv;
      NGLeps = ngoodleps();

/*
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) { leplv = lep1lv; lele  = lep1_is_el(); lmuo = lep1_is_mu(); lpdg = lep1_pdgid(); lpt = lep1_pt(); }
      else if(l1==2) { leplv = lep2lv; lele  = lep2_is_el(); lmuo = lep2_is_mu(); lpdg = lep2_pdgid(); lpt = lep2_pt(); }
*/
      //     if(NGLeps != NLeps) cout << "NGLeps = " << NGLeps << " NLeps = " << NLeps << endl;
      //NGLeps = NLeps;
      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > bjetlv;
      vector<float> btag;
      for(unsigned int nj = 0; nj<ak4pfjets_CSV().size(); ++nj){
	if(ak4pfjets_p4()[nj].Pt()<30) continue;
	if(fabs(ak4pfjets_p4()[nj].Eta())>2.4) continue;
	//if(ak4pfjets_loose_pfid()[nj]==false) continue;
	jind.push_back(nj);
	jetlv.push_back(ak4pfjets_p4()[nj]);
	btag.push_back(ak4pfjets_CSV()[nj]);
	++NGJets;
	if(ak4pfjets_CSV()[nj]>0.814) {++NGBJets; bind.push_back(nj); bjetlv.push_back(ak4pfjets_p4()[nj]);}
      }

      int a = ls(); int b = evt();
      if(
	 //(a==1&&b==42)||(a==1&&b==87)||(a==235&&b==23442)||(a==236&&b==23524)||(a==229&&b==22838)||(a==27&&b==2603)||(a==27&&b==2686)||
	 //(a==1&&b==78)||(a==2&&b==181)||(a==4&&b==381)||(a==236&&b==23527)||(a==305&&b==30452)||(a==305&&b==30493)||(a==67&&b==6661)||
	 //(a==28&&b==2747)||(a==226&&b==22541)||(a==227&&b==22689)||(a==228&&b==22791)||
	 //(a==347&&b==34634)||(a==8&&b==701)||(a==933&&b==93279)||(a==935&&b==93442)||(a==936&&b==93512)
	 (a==2&&b==136)||(a==1&&b==94)||(a==3&&b==210)
         ){
      	print = true;
      }

      if(print){
	cout << "My Event dump: run:ls:evt " << run() << ":" << ls() << ":" << evt() << endl;
	cout << "lepton 1: pt " << lep1_pt() << " eta " << lep1_eta() << " phi " << lep1_phi() << " mass " << lep1_mass() << endl;
	cout << "          pdgID " << lep1_pdgid() << " passVeto " << lep1_passVeto() << " passMedium " << lep1_passMediumID() << " (note true="<<true << ")" << endl;
	cout << "          dz " << lep1_dz() << " d0 " << lep1_d0() << " MiniIso " << lep1_MiniIso() << endl;
	if(lep2_pt()>0){
	  cout << "lepton 2: pt " << lep2_pt() << " eta " << lep2_eta() << " phi " << lep2_phi() << " mass " << lep2_mass() << endl;
	  cout << "          pdgID " << lep2_pdgid() << " passVeto " << lep2_passVeto() << " passMedium " << lep2_passMediumID() << " (note true="<<true << ")" << endl;
	  cout << "          dz " << lep2_dz() << " d0 " << lep2_d0() << " MiniIso " << lep2_MiniIso() << endl;
	}
	for(unsigned int j = 0; j<ak4pfjets_pt().size();++j){
	  cout << "jet " << j+1 << ": pt " << ak4pfjets_pt()[j] << " eta " << ak4pfjets_eta()[j] << " phi " << ak4pfjets_phi()[j] << " mass " << ak4pfjets_mass()[j] << endl;
	  cout << "       CSV " << ak4pfjets_CSV()[j] << " ID " << ak4pfjets_loose_pfid()[j] << " PUid " << ak4pfjets_puid()[j] << " metDPhi " << dphi_ak4pfjet_met()[j] << endl;
	}
	cout << "PassTrackVeto " << PassTrackVeto() << " PassTauVeto " << PassTauVeto() << " rho " << ak4pfjets_rho() << " MET " << pfmet() << " METphi " << pfmet_phi() << endl;
	cout << "MT " << mt_met_lep() << " NJets " << ngoodjets() << " NBJets " << ngoodbtags() << " NLeps " << ngoodleps() << " NVetoLeps " << nvetoleps() << " NGenLeps " << genlepsfromtop() << endl;
	cout << "MT2W " << MT2W() << " chi2 " << hadronic_top_chi2() << " topness " << topness() << " dphi_Wlep " << dphi_Wlep() << " HT " << ak4_HT() << " HTssm " << ak4_htssm() << " HTosm " << ak4_htosm() << endl;
	cout << "Mlb_closestb " << Mlb_closestb() << " Mlb_lead_bdiscr " << Mlb_lead_bdiscr() << " M3b " << M3b() << " dR_lep_leadb " << dR_lep_leadb() << endl;
	cout << "HLT_SingleMu " << HLT_SingleMu() << " HLT_SingleEl " << HLT_SingleEl() << " scale1fb " << scale1fb() << " mass_stop " << mass_stop() << " mass_lsp " << mass_lsp() << " mass_chargino " << mass_chargino() << endl;
	cout << endl;
      }
      /*if(print){
	vector<LorentzVector > mybjets = JetUtil::BJetSelector(ak4pfjets_p4(),ak4pfjets_CSV(),0.814,2,3,2);
	cout << "My input bjets are";
	for(unsigned int j = 0; j<mybjets.size(); ++j) { cout << " b" << j+1 << " " << mybjets[j].Pt(); }
	cout << endl;
	float mytopness=Gettopness_(pfmet(),pfmet_phi(),lep1_p4(),mybjets,0);
	}*/

      if(nvtxs()<0)        continue; ++c1vtx;
      if(NGLeps!=1)        continue; ++c1l;
      if(NSLeps!=1)        continue; ++cno2l;
      //if(!PassTrackVeto()) continue; ++cno2track;
      //if(!PassTauVeto())   continue; ++cnotau;
      if(NGJets<4)         continue; ++c2j;
      if(NGBJets<1)        continue; ++c1b;
      if(pfmet()<80)       continue; ++cmet;

      if(lele!=lep1_is_el()) cout << "lep1_is_el "   << lep1_is_el() << " lele " << lele << " lep2_is_el "   << lep2_is_el() << endl;
      if(lmuo!=lep1_is_mu()) cout << "lep1_is_mu "   << lep1_is_mu() << " lmuo " << lmuo << " lep2_is_mu "   << lep2_is_mu() << endl;
      if(lpdg!=lep1_pdgid()) cout << "lep1_pdgid() " << lep1_pdgid() << " lpdg " << lpdg << " lep2_pdgid() " << lep2_pdgid() << endl;
      if(lpt!=lep1_pt())     cout << "lep1_pt() "    << lep1_pt()    << " lpt "  << lpt  << " lep2_pt() "    << lep2_pt()    << endl;
      if(lep1_miniRelIsoDB()!=lep1_MiniIso() ) cout << "lep1_miniRelIsoDB() "    << lep1_miniRelIsoDB()    << " lep1_MiniIso() "  << lep1_MiniIso()  << endl;

      PrintOut PO;
      PO.Run = run(); PO.LS = ls(); PO.Evt = evt(); PO.leppt = lep1_pt(); PO.leppdgid = lep1_pdgid(); PO.met = pfmet(); PO.mt = mt_met_lep(); PO.njets = ngoodjets(); PO.nbjets = ngoodbtags(); PO.ngenleps = genlepsfromtop();
      PO.DPhiWlep = dphi_Wlep(); PO.htssm = ak4_htssm(); PO.Mlb_cb = Mlb_closestb(); PO.Mjjj = M3b(); PO.mt2w = MT2W(); PO.Topness = topness(); PO.chi2 = hadronic_top_chi2();
      //Print out (4 digits after comma - as agreed in meeting) run:ls:evt:lep1_pt:lep1_pdgid:pfmet:mt_met_lep:ngoodjets:ngoodbtags:genlepsfromtop:dphi_Wlep:ak4_htssm:Mlb_leadb:M3b:MT2W:topness:hadronic_top_chi2 - v2
      if(lele){
	POel.push_back(PO);
	//	cout << "TEST Ele " << fixed << setprecision(4) << run() << " " << evt() << " " << lpt << " " << lpdg << endl;
	//*fLogStreamEl << " " << fixed << setprecision(4) << run() << " " << ls() << " " << evt() << " " << lep1_pt() << " " << lep1_pdgid() << " " << pfmet() << " " << mt_met_lep() << " " << ngoodjets() << " " << ngoodbtags() 
	//	      << " " << genlepsfromtop() << " " << dphi_Wlep() << " " << ak4_htssm() << " " << Mlb_closestb() << " " << M3b() << " " << MT2W() << " " << topness() << " " << hadronic_top_chi2() << endl;
	//*fLogStreamEl << " " << fixed << setprecision(4) << run() << " " << setw(4) << ls() << " " << setw(7) << evt() << " " << setw(8) << lep1_pt() << " " << setw(3) << lep1_pdgid() << " " << setw(9)
	//	      << pfmet() << " " << setw(9) << mt_met_lep() << " " << setw(2) << ngoodjets() << " " << setw(1) << ngoodbtags() << " " << setw(1) << genlepsfromtop() << " " << setw(6)
	//	      << dphi_Wlep() << " " << setw(9) << ak4_htssm() << " " << setw(9) << Mlb_closestb() << " " << setw(9) << M3b() << " " << setw(9) << MT2W() << " " << setw(9) << topness() << " " << setw(9) << hadronic_top_chi2() << endl;
      }
      if(lmuo){
	POmu.push_back(PO);
	//	cout << "TEST Muo " << fixed << setprecision(4) << run() << " " << evt() << " " << lpt << " " << lpdg << endl;
	//*fLogStreamMu << " " << fixed << setprecision(4) << run() << " " << ls() << " " << evt() << " " << lep1_pt() << " " << lep1_pdgid() << " " << pfmet() << " " << mt_met_lep() << " " << ngoodjets() << " " << ngoodbtags() 
	//	      << " " << genlepsfromtop() << " " << dphi_Wlep() << " " << ak4_htssm() << " " << Mlb_closestb() << " " << M3b() << " " << MT2W() << " " << topness() << " " << hadronic_top_chi2() << endl;
	//*fLogStreamMu << " " << fixed << setprecision(4) << run() << " " << setw(4) << ls() << " " << setw(7) << evt() << " " << setw(8) << lep1_pt() << " " << setw(3) << lep1_pdgid() << " " << setw(9)
	//	      << pfmet() << " " << setw(9) << mt_met_lep() << " " << setw(2) << ngoodjets() << " " << setw(1) << ngoodbtags() << " " << setw(1) << genlepsfromtop() << " " << setw(6)
	//	      << dphi_Wlep() << " " << setw(9) << ak4_htssm() << " " << setw(9) << Mlb_closestb() << " " << setw(9) << M3b() << " " << setw(9) << MT2W() << " " << setw(9) << topness() << " " << setw(9) << hadronic_top_chi2() << endl;
      }
      if(mt_met_lep()<120) continue; ++cmt;
      if(mindphi_met_j1_j2()<0.8) continue; ++cmdphi;
      if(hadronic_top_chi2()>5) continue; ++cchi;


    }
  
    // Clean Up
    delete tree;
    file->Close();
    delete file;
  }
  if ( nEventsChain != nEventsTotal ) {
    cout << Form( "ERROR: number of events from files (%d) is not equal to total number of events (%d)", nEventsChain, nEventsTotal ) << endl;
  }
  //std::sort(POel.begin(), POel.end(), &SortPrintOut);
  //std::sort(POmu.begin(), POmu.end(), &SortPrintOut);

  for(unsigned int n=0; n<POel.size();++n){
    PrintOut p = POel[n];
    *fLogStreamEl << " " << fixed << setprecision(4) << p.Run << " " << setw(4) << p.LS << " " << setw(6) << p.Evt << " " << setw(8) << p.leppt << " " << setw(3) << p.leppdgid << " " << setw(9)
		  << p.met << " " << setw(9) << p.mt << " " << setw(2) << p.njets << " " << setw(1) << p.nbjets << " " << setw(1) << p.ngenleps << " " << setw(6)
		  << p.DPhiWlep << " " << setw(9) << p.htssm << " " << setw(9) << p.Mlb_cb << " " << setw(9) << p.Mjjj << " " << setw(9) << p.mt2w << " " << setw(9) << p.Topness << " " << setw(9) << p.chi2 << endl;
  }
    for(unsigned int n=0; n<POmu.size();++n){
    PrintOut p = POmu[n];
    *fLogStreamMu << " " << fixed << setprecision(4) << p.Run << " " << setw(4) << p.LS << " " << setw(6) << p.Evt << " " << setw(8) << p.leppt << " " << setw(3) << p.leppdgid << " " << setw(9)
		  << p.met << " " << setw(9) << p.mt << " " << setw(2) << p.njets << " " << setw(1) << p.nbjets << " " << setw(1) << p.ngenleps << " " << setw(6)
		  << p.DPhiWlep << " " << setw(9) << p.htssm << " " << setw(9) << p.Mlb_cb << " " << setw(9) << p.Mjjj << " " << setw(9) << p.mt2w << " " << setw(9) << p.Topness << " " << setw(9) << p.chi2 << endl;
  }
    
  //for(unsigned int n=0; n<10;++n){
  //  PrintOut p = POel[n];
  //  cout << " " << fixed << setprecision(4) << p.Run << " " << setw(4) << p.LS << " " << setw(7) << p.Evt << " " << setw(8) << p.leppt << " " << setw(3) << p.leppdgid << " " << setw(9)
  //	 << p.met << " " << setw(9) << p.mt << " " << setw(2) << p.njets << " " << setw(1) << p.nbjets << " " << setw(1) << p.ngenleps << " " << setw(6)
  //	 << p.DPhiWlep << " " << setw(9) << p.htssm << " " << setw(9) << p.Mlb_cb << " " << setw(9) << p.Mjjj << " " << setw(9) << p.mt2w << " " << setw(9) << p.Topness << " " << setw(9) << p.chi2 << endl;
  //}

  if(fWriteToFile && fAppend){
    TString lognameEl =outputdir + skimFilePrefix+"_el_SnT.txt"; 
    ofstream f_logEl (lognameEl.Data(), ios::app);
    f_logEl << fLogStreamEl->str();
    cout << "wrote results into  " << lognameEl << " (appended at the end of old file)" << endl;
    TString lognameMu =outputdir +skimFilePrefix+"_mu_SnT.txt"; 
    ofstream f_logMu (lognameMu.Data(), ios::app);
    f_logEl << fLogStreamMu->str();
    cout << "wrote results into  " << lognameMu << " (appended at the end of old file)" << endl;
  }else if(fWriteToFile){
    TString lognameEl =outputdir + skimFilePrefix+"_el_SnT.txt"; 
    ofstream f_logEl (lognameEl.Data(), ios::trunc);
    f_logEl << fLogStreamEl->str();
    cout << "wrote results into  " << lognameEl <<  " (old file replaced)" << endl;
    TString lognameMu =outputdir + skimFilePrefix+"_mu_SnT.txt"; 
    ofstream f_logMu (lognameMu.Data(), ios::trunc);
    f_logMu << fLogStreamMu->str();
    cout << "wrote results into  " << lognameMu <<  " (old file replaced)" << endl;
  } else{
    cout << "Electron events:" << endl;
    cout << fLogStreamEl->str();
    cout << endl << "Muon events:" << endl;
    cout << fLogStreamMu->str();
  }
  delete fLogStreamEl;
  delete fLogStreamMu;
  
  // Example Histograms
  // samplehisto->Draw();
  /*
  for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
  }
  string filename = "rootfiles/MT2MTMassStudies/TestMT2input/HistosFine_"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
  f->Close();
  cout << "Saved histos in " << f->GetName() << endl;
*/
  cout << "For " << skimFilePrefix << ":" << endl;
  cout << "Events passing NVtx>=1: " << c1vtx << endl;
  cout << "Events with one good lepton: " << c1l << endl;
  cout << "Events with one selection lepton: " << cno2l << endl;
  cout << "Events passing track veto: " << cno2track << endl;
  cout << "Events passing tau veto: " << cnotau << endl;
  cout << "Events with at least 4 jets: " << c2j << endl;
  cout << "Events with at least 1 b-jet: " << c1b << endl;
  cout << "Events with MET>80: " << cmet << endl;
  cout << "Events with MT>120: " << cmt << endl;
  cout << "Events with minDPhi>0.8: " << cmdphi << endl;
  cout << "Events with chi2<5: " << cchi << endl;
  
  //cout << endl << "Cutflow" << endl;
  //cout << "Events " << cevtcf << endl << "NVtx>0 " << cvtxcf << endl << "MET>30 " << cmetcf << endl << "NLeps>0 " << cnlepcf << endl << "NJets>1 " << cnjetcf << endl;

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
