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
#include "TH2F.h"

// CMS3
#include "StopBabies10012015.cc"

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
  return TMath::Sqrt(2*met.Et()*lep.Et()*(1-TMath::Cos(JetUtil::deltaPhi(lep,met) ) ) );
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  // Example Histograms
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");

  //--------------------//
  // book 1d histograms //
  //--------------------//
 
  map<string, TH1F*> histos; //massive
  vector<string> histonames; histonames.clear();
  vector<int> histobinn; histobinn.clear();
  vector<double> histobinl; histobinl.clear();
  vector<double> histobinu; histobinu.clear();
  map<string, float> value;

  histonames.push_back("MT2W");               histobinn.push_back(25); histobinl.push_back(0.); histobinu.push_back(500.);
  histonames.push_back("MT");                 histobinn.push_back(40); histobinl.push_back(0.); histobinu.push_back(400.);
  histonames.push_back("MET");                histobinn.push_back(30); histobinl.push_back(0.); histobinu.push_back(750.);

  //--------------------//
  // book 2d histograms //
  //--------------------//

  //2d hists for distributions
  map<string, TH2F*> histos2d;
  vector<string> histonames2d; histonames2d.clear();
  //x
  vector<int> histobinnX; histobinnX.clear();
  vector<float> histobinlX; histobinlX.clear();
  vector<float> histobinuX; histobinuX.clear();
  //y
  vector<int> histobinnY; histobinnY.clear();
  vector<float> histobinlY; histobinlY.clear();
  vector<float> histobinuY; histobinuY.clear();

  map<string, vector<float>> value2d;
  histonames2d.push_back("boostjet_ak10_mass_vs_pt"); 
  histobinnX.push_back(100); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(60); histobinlY.push_back(0.); histobinuY.push_back(300);

  //--------------------// 
  // book Cutflow hists //
  //--------------------//
  map<string, TH1D*> histos_cutflow;
  vector<string> histonames_cutflow; histonames_cutflow.clear();
  histonames_cutflow.push_back("NEventsPerSignalRegion_nominal");
  histonames_cutflow.push_back("NEventsPerSignalRegion_ScaleUp");
  histonames_cutflow.push_back("NEventsPerSignalRegion_ScaleDown");
  histonames_cutflow.push_back("NEventsPerSignalRegion_pdfUp");
  histonames_cutflow.push_back("NEventsPerSignalRegion_pdfDown"); 
  
 //---------------------------------------------------------//
 // initialize histograms here for all regions of interests //
 //---------------------------------------------------------//
 
  //regions
  const unsigned int poststringsize =1 ;
  string poststring[poststringsize] = {"SR"};

  //1d hists
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
 
  //cutflow histograms
   for(unsigned int i = 0;i<histonames_cutflow.size();i++){
    string histoname = histonames_cutflow.at(i)+"_"+samplename;
    histos_cutflow[histoname] = new TH1D(histoname.c_str(),"",10,0,10);
    histos_cutflow[histoname]->Sumw2(); histos_cutflow[histoname]->SetDirectory(rootdir);
    histos_cutflow[histoname]->GetYaxis()->SetTitle("Events");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(1,"MET>250,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(2,"MET>300,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(3,"MET>350,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(4,"MET>400,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(5,"MET>450,lDM");

    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(7,"MET>250,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(8,"MET>300,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(9,"MET>350,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(10,"MET>400,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(11,"MET>450,hDM");
    }

  //2D hists
   for(unsigned int c = 0; c<poststringsize; ++c){
     for(unsigned int i = 0; i<histonames2d.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames2d[i] + "_"+samplename;
	if(histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2F(mapname.c_str(), "", histobinnX[i], histobinlX[i], histobinuX[i],histobinnY[i], histobinlY[i], histobinuY[i]);
	histos2d[mapname]->Sumw2(); histos2d[mapname]->SetDirectory(rootdir);
      }
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

    StopBabies10012015::progress( nEventsTotal, nEventsChain );
    //--------------------------------------------// 
    //       specify processes                    //
    //--------------------------------------------//

    string samplename = skimFilePrefix;
    if(skimFilePrefix.find("TTbar")!=std::string::npos){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
        else samplename = "TTbarH";
      }

     //------------------------------------//
     //  read information from stopbabies  //
     //------------------------------------//
    
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
      float chi2 = cms3.hadronic_top_chi2();
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

     //----------------------------//
     //       lepton selection     //
     //----------------------------//
      int l1=-1;
      if(lep1ismu){
	if(lep1pt>35&&fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1isel){
	if(lep1pt>40&&fabs(lep1eta)<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      NGLeps = NLeps;

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) leplv = lep1lv;
      else if(l1==2) leplv = lep2lv;

     //----------------------------//
     //              met           //
     //----------------------------//

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      metlv.SetPxPyPzE(METx,METy,0.,MET);

     //----------------------------//
     //       jet selection        //
     //----------------------------//
      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
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

     //----------------------------//
     //   apply selection cuts     //
     //----------------------------//
    
      if(nvtxs<0) continue; ++c1vtx;
      if(NGLeps!=1) continue; ++c1l;
      if(NSLeps!=1) continue; ++cno2l;
      if(!trackveto) continue; ++cno2track;
      if(!tauveto) continue; ++cnotau;
      if(NGJets<2) continue; ++c2j;
      if(NGBJets<1) continue; ++c1b;
      if(MET<250) continue; ++cmet;
      if(MT<150) continue;
      if(minDPhi<0.8) continue; ++cmdphi;
      //if(chi2>5) continue; ++cchi;

     //------------------------------//
     // Check which region we are in //
     //------------------------------//

      bool SR = false;

      if(MT>150 && MET>250 && MT2W>-1 ) SR = true;

     //----------------------------//
     //   fill in 1d histograms    //
     //----------------------------//
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
      
      value["MT2W"] = MT2W;
      value["MT"] = MT;
      value["MET"] = MET;

      for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
	if(SR){
	  mname = "SR"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
      }

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }

     //----------------------------//
     //   fill in 2d histograms    //
     //----------------------------//
      
     string d = "_";
     if(SR){
	string mname ="";
	string histname ="";
        //fill each histogram
	histname ="boostjet_ak10_mass_vs_pt";
	mname= "SR"+d+histname+d+samplename;
      //  histos2d[mname]->Fill(boostjet_ak10_pt,boostjet_ak10_mass);
      }

     //----------------------------//
     // Fill in cutflow histograms //
     //----------------------------//
      if(SR){
      if(MET>300&&MET<350) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(0.5,weight);
      if(MET>350&&MET<400) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(1.5,weight);
      if(MET>400&&MET<450) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(2.5,weight);
      if(MET>450&&MET<500) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(3.5,weight);
      if(MET>500) histos_cutflow["NEventsPerSignalRegion_nominal_"+samplename]  ->Fill(4.5,weight);
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
  
  for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
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
  } 
  string filename = "/home/users/mliu/public_html/rootfiles/CutHistos/TTVScaleVar/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
  for(map<string,TH1D*>::iterator h=    histos_cutflow.begin(); h!=    histos_cutflow.end();++h) h->second->Write();
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
