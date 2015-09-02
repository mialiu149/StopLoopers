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
#include "StopBabies20150603.cc"
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
  histonames.push_back("DR");     histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("GenjDR"); histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("GenqDR"); histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("boost500DR");     histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("boost500GenjDR"); histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("boost500GenqDR"); histobinn.push_back(64); histobinl.push_back(0.); histobinu.push_back(6.4);
  histonames.push_back("MET_pre_met170");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_met170");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_pre_ht350met120");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MET_after_ht350met120");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("HT_pre_ht350met120");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("HT_after_ht350met120");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
//  histonames.push_back("HT_MET_pre_ht350met120");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
 // histonames.push_back("MET_after_ht350met120");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
 // histonames.push_back("HT_pre_ht350met120");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
 // histonames.push_back("HT_after_ht350met120");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames2d.push_back("DPtVsDR");    histobinnx.push_back(60); histobinlx.push_back(0.); histobinux.push_back(3.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("LepptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("JetptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("LepjetptVsDR"); histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(750.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenjDPtVsDR");    histobinnx.push_back(60); histobinlx.push_back(0.); histobinux.push_back(3.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenjLepptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenjJetptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenjLepjetptVsDR"); histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(750.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenqDPtVsDR");    histobinnx.push_back(60); histobinlx.push_back(0.); histobinux.push_back(3.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenqLepptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenqJetptVsDR");    histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(500.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  histonames2d.push_back("GenqLepjetptVsDR"); histobinnx.push_back(50); histobinlx.push_back(0.); histobinux.push_back(750.); histobinny.push_back(64); histobinly.push_back(0.); histobinuy.push_back(6.4);
  
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
      StopBabies20150603::progress( nEventsTotal, nEventsChain );

      double weight = scale1fb()*10.;

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
	if(lep1_pt()>30&&fabs(lep1_eta())<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1_is_el()){
	if(lep1_pt()>40&&fabs(lep1_eta())<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      if(lep2_is_mu()){
	if(lep2_pt()>30&&fabs(lep2_eta())<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; if(l1!=1) l1 = 2; else l1=-2; }
      } else if (lep2_is_el()){
	if(lep2_pt()>40&&fabs(lep2_eta())<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps; if(l1!=1) l1 = 2; else l1=-2;}
      }
      NGLeps = ngoodleps();

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) leplv = lep1_p4();
      else if(l1==2) leplv = lep2_p4();
      float METx = pfmet()*TMath::Cos(pfmet_phi());
      float METy = pfmet()*TMath::Sin(pfmet_phi());
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      metlv.SetPxPyPzE(METx,METy,0.,pfmet());
      histos1d["MET_pre_met170_"+samplename]->Fill(metlv.pt(),weight); 
      if(HLT_MET170()) histos1d["MET_after_met170_"+samplename]->Fill(metlv.pt(),weight);
      //if(ak4_HT()>400) 
      histos1d["MET_pre_ht350met120_"+samplename]->Fill(metlv.pt(),weight); 
      //if(metlv.pt()>200) 
      histos1d["HT_pre_ht350met120_"+samplename]->Fill(ak4_HT(),weight); 
      if(HLT_ht350met120()) 
      {
       //if(ak4_HT()>400) 
       histos1d["MET_after_ht350met120_"+samplename]->Fill(metlv.pt(),weight);
       //if(metlv.pt()>200)
       histos1d["HT_after_ht350met120_"+samplename]->Fill(ak4_HT(),weight);
      }
      //if(NGLeps != NLeps) cout << "NGLeps = " << NGLeps << " NLeps = " << NLeps << endl;
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
	if(ak4pfjets_loose_pfid()[nj]==false) continue;
	jind.push_back(nj);
	jetlv.push_back(ak4pfjets_p4()[nj]);
	btag.push_back(ak4pfjets_CSV()[nj]);
	++NGJets;
	if(ak4pfjets_CSV()[nj]>0.814) {++NGBJets; bind.push_back(nj); bjetlv.push_back(ak4pfjets_p4()[nj]);}
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
      for(unsigned int j =0; j<ak4pfjets_p4().size();++j){
	float DR = JetUtil::deltaR(ak4pfjets_p4()[j],leplv);
	float DPtPt = TMath::Abs(ak4pfjets_p4()[j].Pt()-leplv.Pt())/leplv.Pt();
	if(DR<0.01 && DPtPt<0.01) continue;
	histos2d["DPtVsDR_"+samplename]->Fill(DPtPt,DR,weight);
	histos2d["LepptVsDR_"+samplename]->Fill(leplv.Pt(),DR,weight);
	histos2d["JetptVsDR_"+samplename]->Fill(ak4pfjets_p4()[j].Pt(),DR,weight);
	histos2d["LepjetptVsDR_"+samplename]->Fill((ak4pfjets_p4()[j]+leplv).Pt(),DR,weight);
	if(DR<minDR) minDR = DR;
	if(ak4pfjets_p4()[j].Pt()>500 && DR<bminDR) bminDR = DR;
      }
      histos1d["DR_"+samplename]->Fill(minDR,weight);
      if(bminDR<90.) histos1d["boost500DR_"+samplename]->Fill(bminDR,weight);
      for(unsigned int j =0; j<ak4genjets_p4().size();++j){
	float DR = JetUtil::deltaR(ak4genjets_p4()[j],leplv);
	float DPtPt = TMath::Abs(ak4genjets_p4()[j].Pt()-leplv.Pt())/leplv.Pt();
	if(DR<0.01 && DPtPt<0.01) continue;
	histos2d["GenjDPtVsDR_"+samplename]->Fill(DPtPt,DR,weight);
	histos2d["GenjLepptVsDR_"+samplename]->Fill(leplv.Pt(),DR,weight);
	histos2d["GenjJetptVsDR_"+samplename]->Fill(ak4genjets_p4()[j].Pt(),DR,weight);
	histos2d["GenjLepjetptVsDR_"+samplename]->Fill((ak4genjets_p4()[j]+leplv).Pt(),DR,weight);
	if(DR<mingenjDR) mingenjDR = DR;
	if(ak4genjets_p4()[j].Pt()>500 && DR<bmingenjDR) bmingenjDR = DR;
      }
      histos1d["GenjDR_"+samplename]->Fill(mingenjDR,weight);
      if(bmingenjDR<90.) histos1d["boost500GenjDR_"+samplename]->Fill(bmingenjDR,weight);
      for(unsigned int j=0; j<genbs_p4().size();++j){
	if(TMath::Abs(genbs_motherid()[j])<26||(TMath::Abs(genbs_motherid()[j])>1000000&&TMath::Abs(genbs_motherid()[j])<2000016) ){
	  float DR = JetUtil::deltaR(genbs_p4()[j],leplv);
	  float DPtPt = TMath::Abs(genbs_p4()[j].Pt()-leplv.Pt())/leplv.Pt();
	  if(DR<0.01 && DPtPt<0.01) continue;
	  histos2d["GenqDPtVsDR_"+samplename]->Fill(DPtPt,DR,weight);
	  histos2d["GenqLepptVsDR_"+samplename]->Fill(leplv.Pt(),DR,weight);
	  histos2d["GenqJetptVsDR_"+samplename]->Fill(genbs_p4()[j].Pt(),DR,weight);
	  histos2d["GenqLepjetptVsDR_"+samplename]->Fill((genbs_p4()[j]+leplv).Pt(),DR,weight);
	  if(DR<mingenqDR) mingenqDR = DR;  
	  if(genbs_p4()[j].Pt()>500 && DR<bmingenqDR) bmingenqDR = DR;
	  }
      }
      for(unsigned int j=0; j<genqs_p4().size();++j){
	if(TMath::Abs(genqs_motherid()[j])<26||(TMath::Abs(genqs_motherid()[j])>1000000&&TMath::Abs(genqs_motherid()[j])<2000016) ){
	  float DR = JetUtil::deltaR(genqs_p4()[j],leplv);
	  float DPtPt = TMath::Abs(genqs_p4()[j].Pt()-leplv.Pt())/leplv.Pt();
	  if(DR<0.01 && DPtPt<0.01) continue;
	  histos2d["GenqDPtVsDR_"+samplename]->Fill(DPtPt,DR,weight);
	  histos2d["GenqLepptVsDR_"+samplename]->Fill(leplv.Pt(),DR,weight);
	  histos2d["GenqJetptVsDR_"+samplename]->Fill(genqs_p4()[j].Pt(),DR,weight);
	  histos2d["GenqLepjetptVsDR_"+samplename]->Fill((genqs_p4()[j]+leplv).Pt(),DR,weight);
	  if(DR<mingenqDR) mingenqDR = DR;  
	  if(genqs_p4()[j].Pt()>500 && DR<bmingenqDR) bmingenqDR = DR;
	  }
      }
      histos1d["GenqDR_"+samplename]->Fill(mingenqDR,weight);
      if(bmingenqDR<90.) histos1d["boost500GenqDR_"+samplename]->Fill(bmingenqDR,weight);

      for(unsigned int i=0; i<ak4genjets_p4().size(); ++i){
	if(ak4genjets_p4()[i].Pt()<30) continue;
	if(fabs(ak4genjets_p4()[i].Eta())>2.5) continue;
	  LorentzVector reduce = ak4genjets_p4()[i]-leplv;
	  LorentzVector genjet = ak4genjets_p4()[i];
	  if(JetUtil::deltaR(leplv,genjet)>0.5) continue;
	  bool matchedtoq = false; bool matchedtoqp = false;
	  bool matchedtol = false; bool matchedtolp = false;
	  for(unsigned int j =0; j<genbs_p4().size(); ++j){
	    if( (JetUtil::deltaR(genbs_p4()[j],genjet)<0.1) && (TMath::Abs(genbs_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtoq = true;
	    }
	    if( (JetUtil::deltaR(genbs_p4()[j],reduce)<0.1) && (TMath::Abs(genbs_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtoqp = true;
	    }
	  }
	  for(unsigned int j =0; j<genqs_p4().size(); ++j){
	    if( (JetUtil::deltaR(genqs_p4()[j],genjet)<0.1) && (TMath::Abs(genqs_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtoq = true;
	    }
	    if( (JetUtil::deltaR(genqs_p4()[j],reduce)<0.1) && (TMath::Abs(genqs_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtoqp = true;
	    }
	  }
	  for(unsigned int j =0; j<genels_p4().size(); ++j){
	    if( (JetUtil::deltaR(genels_p4()[j],genjet)<0.1) && (TMath::Abs(genels_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtol = true;
	    }
	    if( (JetUtil::deltaR(genels_p4()[j],reduce)<0.1) && (TMath::Abs(genels_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtolp = true;
	    }
	  }
	  for(unsigned int j =0; j<genmus_p4().size(); ++j){
	    if( (JetUtil::deltaR(genmus_p4()[j],genjet)<0.1) && (TMath::Abs(genmus_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtol = true;
	    }
	    if( (JetUtil::deltaR(genmus_p4()[j],reduce)<0.1) && (TMath::Abs(genmus_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtolp = true;
	    }
	  }
	  for(unsigned int j =0; j<genleptau_els_p4().size(); ++j){
	    if( (JetUtil::deltaR(genleptau_els_p4()[j],genjet)<0.1) && (TMath::Abs(genleptau_els_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtol = true;
	    }
	    if( (JetUtil::deltaR(genleptau_els_p4()[j],reduce)<0.1) && (TMath::Abs(genleptau_els_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtolp = true;
	    }
	  }
	  for(unsigned int j =0; j<genleptau_mus_p4().size(); ++j){
	    if( (JetUtil::deltaR(genleptau_mus_p4()[j],genjet)<0.1) && (TMath::Abs(genleptau_mus_p4()[j].Pt()-genjet.Pt())<0.1*genjet.Pt() ) ){
	      matchedtol = true;
	    }
	    if( (JetUtil::deltaR(genleptau_mus_p4()[j],reduce)<0.1) && (TMath::Abs(genleptau_mus_p4()[j].Pt()-reduce.Pt())<0.1*reduce.Pt() ) ){
	      matchedtolp = true;
	    }
	  }
	  float fillindex = -1;
	  if(!matchedtoq && !matchedtoqp && !matchedtol && !matchedtolp) fillindex =  1;
	  else if(!matchedtoq && !matchedtoqp && !matchedtol &&  matchedtolp) fillindex =  2;
	  else if(!matchedtoq && !matchedtoqp &&  matchedtol && !matchedtolp) fillindex =  3;
	  else if(!matchedtoq &&  matchedtoqp && !matchedtol && !matchedtolp) fillindex =  4;
	  else if( matchedtoq && !matchedtoqp && !matchedtol && !matchedtolp) fillindex =  5;
	  else if(!matchedtoq && !matchedtoqp &&  matchedtol &&  matchedtolp) fillindex =  6;
	  else if(!matchedtoq &&  matchedtoqp && !matchedtol &&  matchedtolp) fillindex =  7;
	  else if(!matchedtoq &&  matchedtoqp &&  matchedtol && !matchedtolp) fillindex =  8;
	  else if( matchedtoq && !matchedtoqp && !matchedtol &&  matchedtolp) fillindex =  9;
	  else if( matchedtoq && !matchedtoqp &&  matchedtol && !matchedtolp) fillindex = 10;
	  else if( matchedtoq &&  matchedtoqp && !matchedtol && !matchedtolp) fillindex = 11;
	  else if(!matchedtoq &&  matchedtoqp &&  matchedtol &&  matchedtolp) fillindex = 12;
	  else if( matchedtoq && !matchedtoqp &&  matchedtol &&  matchedtolp) fillindex = 13;
	  else if( matchedtoq &&  matchedtoqp && !matchedtol &&  matchedtolp) fillindex = 14;
	  else if( matchedtoq &&  matchedtoqp &&  matchedtol && !matchedtolp) fillindex = 15;
	  else if( matchedtoq &&  matchedtoqp &&  matchedtol &&  matchedtolp) fillindex = 16;
	  else  fillindex = 19;
//	  histos1d["GenJetMatch_"+samplename]->Fill(fillindex,weight);
      }
  
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

  string filename = "rootfiles/DRHistos/Histos_"+skimFilePrefix+".root";
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
