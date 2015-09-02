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
  histonames.push_back("lep_pt");             histobinn.push_back(27); histobinl.push_back(0.); histobinu.push_back(405.);
  histonames.push_back("lep_eta");            histobinn.push_back(27); histobinl.push_back(0.); histobinu.push_back(405.);
  histonames.push_back("lep_phi");            histobinn.push_back(27); histobinl.push_back(0.); histobinu.push_back(405.);
  histonames.push_back("MET");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);
  histonames.push_back("MT");                 histobinn.push_back(40); histobinl.push_back(0.); histobinu.push_back(400.);
  histonames.push_back("HT");                 histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(1000.);
  histonames.push_back("NJets");              histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(10.);

  const unsigned int poststringsize = 2;
  string poststring[poststringsize] = {"noMET","MET50"};

  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix!="TTbar"&&b>0) continue;
    if(skimFilePrefix=="TTbar"&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix=="TTbar"&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix=="TTbar"&&b==2) samplename = "TTbarH";
    //cout << b << " " << samplename << endl;
    for(unsigned int c = 0; c<poststringsize; ++c){
      for(unsigned int i = 0; i<histonames.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames[i] + "_"+samplename;
	//cout << mapname << endl;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1F(mapname.c_str(), "", histobinn[i], histobinl[i], histobinu[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	if(c==0) value[histonames[i] ] = -1;
      }
    }
  }

  int c1vtx(0), c1l(0), cno2l(0), cno2track(0), cnotau(0), c2j(0), c1b(0), cmet(0);
  int cmt(0), cmdphi(0), cchi(0);
  //  unsigned int cevtcf(0), cvtxcf(0), cmetcf(0), cnlepcf(0), cnjetcf(0);

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
      //unsigned int run = cms3.run();
      //unsigned int ls = cms3.ls();
      //unsigned int evt = cms3.evt();
      float weight = cms3.scale1fb()*10.;
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
      //float genmet = cms3.genmet();
      //int NJets = cms3.ak4GoodPFJets();
      //int NBJets = cms3.ak4_nBTags_Med();
      float HT = cms3.ak4_HT();
      float HTratio = cms3.ak4_htratiom();
      int nvtxs = cms3.nvtxs();
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
//      if(NSLeps!=1) continue; ++cno2l;
//      if(!trackveto) continue; ++cno2track;
//      if(!tauveto) continue; ++cnotau;
//      if(NGJets<4) continue; ++c2j;
//      if(NGBJets<1) continue; ++c1b;
//      if(MET<0) continue; ++cmet;
      bool noMET = false; bool MET50 = false; 
      if(MET>200){
	noMT = true;
	if(MT>80) MT80 = true;
	if(MT>100) MT100 = true;
	if(MT>125) MT125 = true;
	if(MT>150) MT150 = true;
      }

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
      //vector<LorentzVector > btaggedjets = JetUtil::BJetSelector(jetlv,btag,0.814,2,3,1);
      vector<LorentzVector > btaggedjets = JetUtil::BJetSelector(jetlv,btag,0.814,2,2,2);

      vector<LorentzVector > dummybjets;dummybjets.clear();
      LorentzVector dummybjet; dummybjet.SetPxPyPzE(0.,0.,0.,0.); dummybjets.push_back(dummybjet);dummybjets.push_back(dummybjet);

      int leadb(-1), trailb(-1);
      for(unsigned int n = 0; n<btaggedjets.size();++n){
	if(leadb<0) leadb = n;
	else if(trailb<0) trailb = n;
	else if(btaggedjets[n].Pt() > btaggedjets[leadb].Pt()){ trailb = leadb; leadb = n; }
	else if(btaggedjets[n].Pt() > btaggedjets[trailb].Pt()){ trailb = n; }
      }
      float myleadjpt = -1.;
      int bj1(-1), bj2(-1),bj3(-1);
      for(unsigned int n = 0; n<jetlv.size();++n){
	float tmp = getMT(jetlv[n],metlv);
  	if(tmp>value["MTqmax"]) value["MTqmax"] = tmp;  else if(value["MTqmax"]<0) value["MTqmax"] = tmp;
	if(jetlv[n].Pt()>250.){ if(tmp<value["MTq_boost"]) value["MTq_boost"] = tmp;  else if(value["MTq_boost"]<0) value["MTq_boost"] = tmp; }
	if(jetlv[n].Pt()>250.){ if(tmp>value["MTq_boostmax"]) value["MTq_boostmax"] = tmp;  else if(value["MTq_boostmax"]<0) value["MTq_boostmax"] = tmp; }
	if(jetlv[n].Pt()>200.){ if(tmp<value["MTq_boost200"]) value["MTq_boost200"] = tmp;  else if(value["MTq_boost200"]<0) value["MTq_boost200"] = tmp; }
	if(jetlv[n].Pt()>300.){ if(tmp<value["MTq_boost300"]) value["MTq_boost300"] = tmp;  else if(value["MTq_boost300"]<0) value["MTq_boost300"] = tmp; }
	if(jetlv[n].Pt()>400.){ if(tmp<value["MTq_boost400"]) value["MTq_boost400"] = tmp;  else if(value["MTq_boost400"]<0) value["MTq_boost400"] = tmp; }
	if(jetlv[n].Pt()>450.){ if(tmp<value["MTq_boost450"]) value["MTq_boost450"] = tmp;  else if(value["MTq_boost450"]<0) value["MTq_boost450"] = tmp; }
	if(jetlv[n].Pt()>500.){ if(tmp<value["MTq_boost500"]) value["MTq_boost500"] = tmp;  else if(value["MTq_boost500"]<0) value["MTq_boost500"] = tmp; }
	if(jetlv[n].Pt()>myleadjpt){ value["MTq_boostLeadJet"] = tmp; myleadjpt = jetlv[n].Pt(); }// else if(value["MTq_boostLeadJet"]<0) value["MTq_boostLeadJet"] = tmp;
	if(jetlv[n].Pt()>300. && jetlv[n].Pt()>myleadjpt300){ value["MTq_boostLeadJet300"] = tmp; myleadjpt300 = jetlv[n].Pt(); }// else if(value["MTq_boostLeadJet300"]<0) value["MTq_boostLeadJet300"] = tmp;
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
      
      value["MT"] = MT;
      value["MET"] = MET;
      value["HT"] = HT;
      value["NJets"] = NGJets;
      value["lep_pt"] = leplv.Pt();
      value["lep_eta"] = leplv.Eta();
      value["lep_phi"] = leplv.Phi();
//      value["pTleadj"] = myleadjpt;
//      value["pTleadb"] = btaggedjets[leadb].Pt();
//      value["pTtrailb"] = btaggedjets[trailb].Pt();
      //cout << __LINE__<<endl;
      for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
	if(noMET){
	  mname= "noMET"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(MET50){
	  mname= "MET50"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
      }

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
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

  string filename = "rootfiles/CutHistos/StackPlots/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
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
