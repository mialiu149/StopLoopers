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
//#include "CMS3.cc"
#include "StopBabiesPhys14_06252015.cc"

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

  map<string, TH1D*> histos;

  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix!="TTbar"&&b>0) continue;
    if(skimFilePrefix=="TTbar"&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix=="TTbar"&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix=="TTbar"&&b==2) samplename = "TTbarH";

    string histoname = "NEventsPerSignalRegion";
    histos[histoname] = new TH1D(histoname.c_str(),"",8,0,8);
    histos[histoname]->Sumw2(); histos[histoname]->SetDirectory(rootdir);
    histos[histoname]->GetYaxis()->SetTitle("Events");
    histos[histoname]->GetXaxis()->SetBinLabel(1,"MET>150,lDM");
    histos[histoname]->GetXaxis()->SetBinLabel(2,"MET>200,lDM");
    histos[histoname]->GetXaxis()->SetBinLabel(3,"MET>250,lDM");
    histos[histoname]->GetXaxis()->SetBinLabel(4,"MET>300,lDM");
    histos[histoname]->GetXaxis()->SetBinLabel(5,"MET>150,hDM");
    histos[histoname]->GetXaxis()->SetBinLabel(6,"MET>200,hDM");
    histos[histoname]->GetXaxis()->SetBinLabel(7,"MET>250,hDM");
    histos[histoname]->GetXaxis()->SetBinLabel(8,"MET>300,hDM");
    histoname = "CutflowLowDM";
    histos[histoname] = new TH1D(histoname.c_str(),"",15,0,15);
    histos[histoname]->Sumw2(); histos[histoname]->SetDirectory(rootdir);
    histos[histoname]->GetYaxis()->SetTitle("Events");
    histos[histoname]->GetXaxis()->SetBinLabel(1,"NVtx");
    histos[histoname]->GetXaxis()->SetBinLabel(2,"GLeps");
    histos[histoname]->GetXaxis()->SetBinLabel(3,"SLeps");
    histos[histoname]->GetXaxis()->SetBinLabel(4,"track");
    histos[histoname]->GetXaxis()->SetBinLabel(5,"tau");
    histos[histoname]->GetXaxis()->SetBinLabel(6,"NJets");
    histos[histoname]->GetXaxis()->SetBinLabel(7,"NBJets");
    //histos[histoname]->GetXaxis()->SetBinLabel(1,"Preselection");
    histos[histoname]->GetXaxis()->SetBinLabel(8,"$M_{T}>120$");
    histos[histoname]->GetXaxis()->SetBinLabel(9,"$E_{T}^{miss}>100$");
    histos[histoname]->GetXaxis()->SetBinLabel(10,"min#Delta#phi>0.8");
    histos[histoname]->GetXaxis()->SetBinLabel(11,"#chi^{2}<5");
    histos[histoname]->GetXaxis()->SetBinLabel(12,"E_{T}^{miss}>150");
    histos[histoname]->GetXaxis()->SetBinLabel(13,"E_{T}^{miss}>200");
    histos[histoname]->GetXaxis()->SetBinLabel(14,"E_{T}^{miss}>250");
    histos[histoname]->GetXaxis()->SetBinLabel(15,"E_{T}^{miss}>300");
    histoname = "CutflowHighDM";
    histos[histoname] = new TH1D(histoname.c_str(),"",16,0,16);
    histos[histoname]->Sumw2(); histos[histoname]->SetDirectory(rootdir);
    histos[histoname]->GetYaxis()->SetTitle("Events");
    histos[histoname]->GetXaxis()->SetBinLabel(1,"NVtx");
    histos[histoname]->GetXaxis()->SetBinLabel(2,"GLeps");
    histos[histoname]->GetXaxis()->SetBinLabel(3,"SLeps");
    histos[histoname]->GetXaxis()->SetBinLabel(4,"track");
    histos[histoname]->GetXaxis()->SetBinLabel(5,"tau");
    histos[histoname]->GetXaxis()->SetBinLabel(6,"NJets");
    histos[histoname]->GetXaxis()->SetBinLabel(7,"NBJets");
    //histos[histoname]->GetXaxis()->SetBinLabel(1,"Preselection");
    histos[histoname]->GetXaxis()->SetBinLabel(8,"$M_{T}>120$");
    histos[histoname]->GetXaxis()->SetBinLabel(9,"$E_{T}^{miss}>100$");
    histos[histoname]->GetXaxis()->SetBinLabel(10,"min#Delta#phi>0.8");
    histos[histoname]->GetXaxis()->SetBinLabel(11,"#chi^{2}<5");
    histos[histoname]->GetXaxis()->SetBinLabel(12,"M_{T2}^{W}>200");
    histos[histoname]->GetXaxis()->SetBinLabel(13,"E_{T}^{miss}>150");
    histos[histoname]->GetXaxis()->SetBinLabel(14,"E_{T}^{miss}>200");
    histos[histoname]->GetXaxis()->SetBinLabel(15,"E_{T}^{miss}>250");
    histos[histoname]->GetXaxis()->SetBinLabel(16,"E_{T}^{miss}>300");
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
      StopBabiesPhys14_06252015::progress( nEventsTotal, nEventsChain );

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
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > lep2lv = cms3.lep2_p4();
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

      int NGLeps = cms3.ngoodleps();
      int NSLeps = 0;
      int NGJets = 0;
      int NGBJets = 0;

      int l1=-1; int muon_count; int el_count;
      if(lep1ismu){
	if(lep1pt>30&&fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps;++muon_count ;l1 = 1; }
      } else if (lep1isel){
	if(lep1pt>40&&fabs(lep1eta)<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps;++el_count;l1 = 1; }
      }
      if(lep2ismu){
	if(lep2pt>30&&fabs(lep2eta)<2.1&&fabs(cms3.lep2_d0())<0.02&&fabs(cms3.lep2_dz())<0.1&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps;++muon_count ;l1 = 2; }
      } else if (lep2isel){
	if(lep2pt>40&&fabs(lep2eta)<2.1&&cms3.lep2_is_phys14_medium_noIso()&&cms3.lep2_miniRelIsoDB()<0.1) {++NSLeps;++el_count;l1 = 2; }
      }

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(NSLeps==1) {if(l1==1) leplv = lep1lv;
      else if(l1==2) leplv = lep2lv;}
      if(NSLeps>1&&l1==1) leplv = lep1lv;
      else leplv = lep2lv;

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      metlv.SetPxPyPzE(METx,METy,0.,MET);

      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > boostjetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > bjetlv;
      vector<float> btag;

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
      }
//      if(HT<400) continue;
      double weighttemp = 1;
  //    if(!HLT_MET170()&&!HLT_ht350met120()) continue;
      if(nvtxs>=1){
	histos["CutflowLowDM" ]->Fill(0.5,weight);
	histos["CutflowHighDM"]->Fill(0.5,weight);
	if(NGLeps==1){
	  histos["CutflowLowDM" ]->Fill(1.5,weight);
	  histos["CutflowHighDM"]->Fill(1.5,weight);
	  if(NSLeps==1){
	    histos["CutflowLowDM" ]->Fill(2.5,weight);
	    histos["CutflowHighDM"]->Fill(2.5,weight);
	  if(trackveto){
	    histos["CutflowLowDM" ]->Fill(3.5,weight);
	    histos["CutflowHighDM"]->Fill(3.5,weight);
	    if(tauveto){
	      histos["CutflowLowDM" ]->Fill(4.5,weight);
	      histos["CutflowHighDM"]->Fill(4.5,weight);
	      if(NGJets>=4){
		histos["CutflowLowDM" ]->Fill(5.5,weight);
		histos["CutflowHighDM"]->Fill(5.5,weight);
		if(NGBJets>=1){
		  histos["CutflowLowDM" ]->Fill(6.5,weight);
		  histos["CutflowHighDM"]->Fill(6.5,weight);
		  if(MT>120){
		    histos["CutflowLowDM" ]->Fill(7.5,weight);
		    histos["CutflowHighDM"]->Fill(7.5,weight);
		    if(MET>100){
		      histos["CutflowLowDM" ]->Fill(8.5,weight);
		      histos["CutflowHighDM"]->Fill(8.5,weight);
		      if(minDPhi>0.8){
			histos["CutflowLowDM" ]->Fill(9.5,weight);
			histos["CutflowHighDM"]->Fill(9.5,weight);
			if(chi2<5){
			  histos["CutflowLowDM" ]->Fill(10.5,weight);
			  histos["CutflowHighDM"]->Fill(10.5,weight);
			  if(MET>150) histos["CutflowLowDM"] ->Fill(11.5,weight);
			  if(MET>200) histos["CutflowLowDM"] ->Fill(12.5,weight);
			  if(MET>250) histos["CutflowLowDM"] ->Fill(13.5,weight);
			  if(MET>300) histos["CutflowLowDM"] ->Fill(14.5,weight);
			  if(MT2W>200){
			    histos["CutflowHighDM"]->Fill(11.5,weight);
			    if(MET>150) histos["CutflowHighDM"] ->Fill(12.5,weight);
			    if(MET>200) histos["CutflowHighDM"] ->Fill(13.5,weight);
			    if(MET>250) histos["CutflowHighDM"] ->Fill(14.5,weight);
			    if(MET>300) histos["CutflowHighDM"] ->Fill(15.5,weight);
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      }

      //if(nvtxs<0) continue; ++c1vtx;
      if(firstGoodVtxIdx() != 0) continue; ++c1vtx;
      if(NGLeps!=1) continue; ++c1l;
      if(NSLeps!=1) continue; ++cno2l;
      if(!trackveto) continue; ++cno2track;
      if(!tauveto) continue; ++cnotau;
      if(NGJets<4) continue; ++c2j;
      if(NGBJets<1) continue; ++c1b;
      if(MT<120) continue; ++cmt;
      if(minDPhi<0.8) continue; ++cmdphi;
      if(chi2>5) continue; ++cchi;
      if(MET<150) continue;

      if(MET>150) histos["NEventsPerSignalRegion"] ->Fill(0.5,weight);
      if(MET>200) histos["NEventsPerSignalRegion"]  ->Fill(1.5,weight);
      if(MET>250) histos["NEventsPerSignalRegion"]  ->Fill(2.5,weight);
      if(MET>300) histos["NEventsPerSignalRegion"]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(MET>150) histos["NEventsPerSignalRegion"]  ->Fill(4.5,weight);
	if(MET>200) histos["NEventsPerSignalRegion"]  ->Fill(5.5,weight);
	if(MET>250) histos["NEventsPerSignalRegion"]  ->Fill(6.5,weight);
	if(MET>300) histos["NEventsPerSignalRegion"]  ->Fill(7.5,weight);
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
  


  string filename = "rootfiles/CutflowExtra/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1D*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
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
