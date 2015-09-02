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
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TTreeCache.h"
#include "TLorentzVector.h"

// CMS3
//#include "CMS3.cc"
#include "StopBabiesRun2_07142015.cc"
//#include "/home/users/mliu/Tools/goodrun.h"
//#include "/home/users/mliu/Tools/goodrun.cc"
//MT2 variants
//#include "stop_variables/Davismt2.h"
//#include "stop_variables/topness.h"
//#include "stop_variables/MT2_implementations.h"
#include "stop_variables/JetUtil.h"
//#include "stop_variables/mt2w.h"
//#include "stop_variables/mt2w_bisect.h"

using namespace std;
using namespace tas;

inline float getMT(LorentzVector lep,LorentzVector met){
  // From cmssw reco::deltaPhi()
  return TMath::Sqrt(2*met.Et()*lep.Et()*(1-TMath::Cos(JetUtil::deltaPhi(lep,met) ) ) );
}

inline float dot(LorentzVector a,LorentzVector b){
  // From cmssw reco::deltaPhi()
  return a.Px()*b.Px()+a.Py()*b.Py()+a.Pz()*b.Pz();
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  map<string, TH1F*> histos1d; //massive
  vector<string> histonames; histonames.clear();
  vector<int> histobinn; histobinn.clear();
  vector<double> histobinl; histobinl.clear();
  vector<double> histobinu; histobinu.clear();
  map<string, float> value;
  map<string, TH2F*> histos2d; //massive
  vector<string> histonames2d; histonames2d.clear();
  vector<int> histobinnx; histobinnx.clear();
  vector<double> histobinlx; histobinlx.clear();
  vector<double> histobinux; histobinux.clear();
  vector<int> histobinny; histobinny.clear();
  vector<double> histobinly; histobinly.clear();
  vector<double> histobinuy; histobinuy.clear();

  histonames.push_back("MET_pre_met170");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met_phi");                histobinn.push_back(30); histobinl.push_back(-3.4); histobinu.push_back(3.4);
  histonames.push_back("MET_after_met170");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_el");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_el");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_el_offline20");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_el_offline20");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_el_offline40");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_el_offline40");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_mu");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_mu");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_mu_offline20");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_mu_offline20");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_met170_mu_offline30");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170_mu_offline30");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);

  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix!="TTbar"&&b>0) continue;
    if(skimFilePrefix=="TTbar"&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix=="TTbar"&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix=="TTbar"&&b==2) samplename = "TTbarH";
    for(unsigned int i = 0; i<histonames.size(); ++i){
      string mapname;
      mapname = histonames[i] + "_"+samplename;
      if(histos1d.count(mapname) == 0 ) histos1d[mapname] = new TH1F(mapname.c_str(), "", histobinn[i], histobinl[i], histobinu[i]);
      histos1d[mapname]->Sumw2(); histos1d[mapname]->SetDirectory(rootdir);
      value[histonames[i] ] = -1;
    }
    for(unsigned int i = 0; i<histonames2d.size(); ++i){
      string mapname;
     mapname = histonames2d[i] + "_"+samplename;
     //cout << mapname << endl;
      if(histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2F(mapname.c_str(), "", histobinnx[i], histobinlx[i], histobinux[i], histobinny[i], histobinly[i], histobinuy[i]);
      histos2d[mapname]->Sumw2(); histos2d[mapname]->SetDirectory(rootdir);
      value[histonames2d[i] ] = -1;
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

  //load json file
//  const char* json_file = "json_DCSONLY_Run2015B_snt.txt";
  //set_goodrun_file(json_file);
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
      //CMS3::progress( nEventsTotal, nEventsChain );
      StopBabiesRun2_07142015::progress( nEventsTotal, nEventsChain );
//      if( !goodrun(run(), ls()) ) continue;
      //if(run()!=251244) continue;   
      //if(filt_hbhenoise()) continue; 
      //double weight = scale1fb()*10.;
      double weight = 1;

      string samplename = skimFilePrefix;
      if(skimFilePrefix=="TTbar"){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
	else samplename = "TTbarH";
      }
      // Analysis Code
      int NGLeps = 0;
      int NSLeps = 0;
      int NGJets = 0;
      int NGBJets = 0;

      int l1=-1;
      if(lep1_is_mu()){
	if(lep1_pt()>20&&fabs(lep1_eta())<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1_is_el()){
	if(lep1_pt()>20&&fabs(lep1_eta())<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      if(lep2_is_mu()){
	if(lep2_pt()>20&&fabs(lep2_eta())<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; if(l1!=1) l1 = 2; else l1=-2; }
      } else if (lep2_is_el()){
	if(lep2_pt()>20&&fabs(lep2_eta())<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
      }
      NGLeps = ngoodleps();

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) leplv = lep1_p4();
      else if(l1==2) leplv = lep2_p4();
      bool isMu = false;
      bool isEl = false;
      if(l1==1&&lep1_is_mu()||l1==2&&lep2_is_mu()) isMu=true;
      if(l1==1&&lep1_is_el()||l1==2&&lep2_is_el()) isEl=true;
     //select jets.
      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > bjetlv;
      vector<float> btag;
      for(unsigned int nj = 0; nj<ak4pfjets_CSV().size(); ++nj){
	if(ak4pfjets_p4()[nj].Pt()<30) continue;
	if(fabs(ak4pfjets_p4()[nj].Eta())>2.4) continue;
	if(ak4pfjets_loose_pfid()[nj]==false) continue;
	jind.push_back(nj);
	jetlv.push_back(ak4pfjets_p4()[nj]);
	btag.push_back(ak4pfjets_CSV()[nj]);
	++NGJets;
	if(ak4pfjets_CSV()[nj]>0.814) {++NGBJets; bind.push_back(nj); bjetlv.push_back(ak4pfjets_p4()[nj]);}
      }
      //met
      float METx = pfmet()*TMath::Cos(pfmet_phi());
      float METy = pfmet()*TMath::Sin(pfmet_phi());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      metlv.SetPxPyPzE(METx,METy,0.,pfmet());
      //if(!(ak4_HT()>10)) continue;
      if(pfmet()<300) continue; ++c2j;
      histos1d["MET_pre_met170_"+samplename]->Fill(metlv.pt(),weight); 
     if(NGJets==0)  histos1d["MET_pre_met_phi_"+samplename]->Fill(metlv.phi(),weight); 

      if(HLT_MET170()) histos1d["MET_after_met170_"+samplename]->Fill(metlv.pt(),weight);

      if(isEl&&leplv.pt()>20){
         histos1d["MET_pre_met170_el_offline20_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_el_offline20_"+samplename]->Fill(metlv.pt(),weight);
      }
      if(isEl&&leplv.pt()>40){
         histos1d["MET_pre_met170_el_offline40_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_el_offline40_"+samplename]->Fill(metlv.pt(),weight);
      }
      if(isMu&&leplv.pt()>20){
         histos1d["MET_pre_met170_mu_offline20_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_mu_offline20_"+samplename]->Fill(metlv.pt(),weight);
      }
      if(isMu&&leplv.pt()>30){
         histos1d["MET_pre_met170_mu_offline30_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_mu_offline30_"+samplename]->Fill(metlv.pt(),weight);
      }
      if(HLT_SingleEl()){
         histos1d["MET_pre_met170_el_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_el_"+samplename]->Fill(metlv.pt(),weight);
      }
      if(HLT_SingleMu()){
         histos1d["MET_pre_met170_mu_"+samplename]->Fill(metlv.pt(),weight);
         if(HLT_MET170()) histos1d["MET_after_met170_mu_"+samplename]->Fill(metlv.pt(),weight);
      }


      if(nvtxs()<0) continue; ++c1vtx;
      if(NGLeps!=1) continue; ++c1l;
      if(NSLeps!=1) continue; ++cno2l;
      if(!PassTrackVeto()) continue; ++cno2track;
      if(!PassTauVeto()) continue; ++cnotau;
      //if(NGJets<4) continue; ++c2j;
      //if(NGBJets<1) continue; ++c1b;
      //if(pfmet()<100) continue; ++cmet;
      //if(mt_met_lep()<120) continue; ++cmt;

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
	    
      float minDR(99), mingenjDR(99), mingenqDR(99);
      float bminDR(99), bmingenjDR(99), bmingenqDR(99);
      //closest jet to lepton
      //closest gen-quark to lepton
      //2d plot of lepton pt vs DR
      //2d plot of jet pt vs DR
      //2d plot of lepton+jet pt vs DR
  
      //for(unsigned int i = 0; i<histonames.size(); ++i){
      //string d = "_";
      //string mname;
      //mname= histonames[i]+d+samplename;
      //histos[mname]->Fill(value[histonames[i] ],weight);
      //}
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
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
  
  // Example Histograms
  // samplehisto->Draw();

  for(map<string,TH1F*>::iterator h=histos1d.begin(); h!=histos1d.end();++h){
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
  }

  for(map<string,TH2F*>::iterator h=histos2d.begin(); h!=histos2d.end();++h){
    for(int by = 1; by<=h->second->GetNbinsY();++by){
      h->second->SetBinContent(h->second->GetNbinsX(), by, h->second->GetBinContent(h->second->GetNbinsX(),by)+ h->second->GetBinContent(h->second->GetNbinsX()+1,by) );
      h->second->SetBinError(h->second->GetNbinsX(), by, sqrt(pow(h->second->GetBinError(h->second->GetNbinsX(),by),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1,by),2) ) );
    }
    for(int bx = 1; bx<=h->second->GetNbinsX();++bx){
      h->second->SetBinContent(bx, h->second->GetNbinsY(), h->second->GetBinContent(bx, h->second->GetNbinsY() )+ h->second->GetBinContent(bx, h->second->GetNbinsY()+1) );
      h->second->SetBinError(bx, h->second->GetNbinsY(), sqrt(pow(h->second->GetBinError(bx, h->second->GetNbinsY() ),2)+pow(h->second->GetBinError(bx, h->second->GetNbinsY()+1),2) ) );
    }
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			     h->second->GetBinContent(h->second->GetNbinsX(),h->second->GetNbinsY() )+ h->second->GetBinContent(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			   sqrt(pow(h->second->GetBinError(h->second->GetNbinsX(),h->second->GetNbinsY() ),2)+
				pow(h->second->GetBinError(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1),2) ) );
  }

  string filename = "rootfiles/TriggerEffs/Histos_"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos1d.begin(); h!=    histos1d.end();++h) h->second->Write();
  for(map<string,TH2F*>::iterator h=    histos2d.begin(); h!=    histos2d.end();++h) h->second->Write();
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
