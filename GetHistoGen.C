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
//#include "TH1F.h"
#include "TH2F.h"

// CMS3
#include "StopBabies08052015.cc"

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
// Delare histonames and bins.
  map<string, TH1F*> histos; //1d hists for distributions
  vector<string> histonames; histonames.clear();
  vector<int> histobinn; histobinn.clear();
  vector<double> histobinl; histobinl.clear();
  vector<double> histobinu; histobinu.clear();
  map<string, float> value;

  map<string, TH1F*> histos_vec; //1d hists for distributions
  vector<string> histonames_vec; histonames_vec.clear();
  vector<int> histobinn_vec; histobinn_vec.clear();
  vector<double> histobinl_vec; histobinl_vec.clear();
  vector<double> histobinu_vec; histobinu_vec.clear();
  map<string, vector<float>> value_vec;

  histonames.push_back("genb_n");     histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("gentop_n");   histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("genstop_n");  histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("genel_n");    histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("genmu_n");    histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("gentau_n");   histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);
  histonames.push_back("genqs_n");    histobinn.push_back(10); histobinl.push_back(-0.5); histobinu.push_back(10.);

  histonames_vec.push_back("gen_dRtt");     histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(5.);
  histonames_vec.push_back("gen_dRlt");     histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(5.);
  histonames_vec.push_back("gen_dRlw");     histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(5.);

  histonames.push_back("gentop_pt");  histobinn.push_back(85); histobinl.push_back(0.); histobinu.push_back(1700.);
  histonames.push_back("gentop_eta"); histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("gen_wb");    histobinn.push_back(100); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("gentop_phi"); histobinn.push_back(10); histobinl.push_back(-5); histobinu.push_back(5.);

  map<string, TH2F*> histos2d; //2d hists for distributions
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
  histonames2d.push_back("gentop_pt_vs_dRWlep_b"); 
  histobinnX.push_back(100); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(100); histobinlY.push_back(0.); histobinuY.push_back(3.);
  histonames2d.push_back("gentop_pt_vs_dRWhad_b"); 
  histobinnX.push_back(100); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(100); histobinlY.push_back(0.); histobinuY.push_back(3.);
  histonames2d.push_back("genw_pt_vs_dRqq");
  histobinnX.push_back(100); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(100); histobinlY.push_back(0.); histobinuY.push_back(3.);

  const unsigned int poststringsize = 2;
  string poststring[poststringsize] = {"noMET","MET300"};

  //initialize histos

  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix!="TTbar"&&b>0) continue;
    if(skimFilePrefix=="TTbar"&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix=="TTbar"&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix=="TTbar"&&b==2) samplename = "TTbarH";
  //book 1D hists
    for(unsigned int c = 0; c<poststringsize; ++c){
      for(unsigned int i = 0; i<histonames.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames[i] + "_"+samplename;
	if(histos.count(mapname) == 0 ) histos[mapname] = new TH1F(mapname.c_str(), "", histobinn[i], histobinl[i], histobinu[i]);
	histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
	if(c==0) value[histonames[i] ] = -1;
      }
    }
     for(unsigned int c = 0; c<poststringsize; ++c){
      for(unsigned int i = 0; i<histonames_vec.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames_vec[i] + "_"+samplename;
	if(histos_vec.count(mapname) == 0 ) histos_vec[mapname] = new TH1F(mapname.c_str(), "", histobinn_vec[i], histobinl_vec[i], histobinu_vec[i]);
	histos_vec[mapname]->Sumw2(); histos_vec[mapname]->SetDirectory(rootdir);
	//if(c==0) value_vec[histonames[i] ] = -1;
      }
    }
  //book 2D hists
   for(unsigned int c = 0; c<poststringsize; ++c){
     for(unsigned int i = 0; i<histonames2d.size(); ++i){
	string mapname;
	mapname = poststring[c] + "_" + histonames2d[i] + "_"+samplename;
	if(histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2F(mapname.c_str(), "", histobinnX[i], histobinlX[i], histobinuX[i],histobinnY[i], histobinlY[i], histobinuY[i]);
	histos2d[mapname]->Sumw2(); histos2d[mapname]->SetDirectory(rootdir);
//	if(c==0) value2d[histonames2d[i] ].clear();
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
   
      // Progress
      StopBabies08052015::progress( nEventsTotal, nEventsChain );

      string samplename = skimFilePrefix;
      if(skimFilePrefix=="TTbar"){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
	else samplename = "TTbarH";
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

      vector<int> bs;
      vector<int> ts;
      vector<int> qs_fromWplus;
      vector<int> qs_fromWminus;
      vector<int> els_fromWplus;
      vector<int> els_fromWminus;
      vector<int> mus_fromWplus;
      vector<int> mus_fromWminus;
      vector<int> taus_fromWplus;
      vector<int> taus_fromWminus;
      vector<int> nuels_fromWplus;
      vector<int> nuels_fromWminus;
      vector<int> numus_fromWplus;
      vector<int> numus_fromWminus;
      vector<int> nutaus_fromWplus;
      vector<int> nutaus_fromWminus;


      vector<LorentzVector> ws_lep;
      vector<LorentzVector> ws_had;
      vector<int> gen_els;
      vector<bool> fromTop;
      LorentzVector bquark;
      LorentzVector bbarquark;

      for (int i=0;i<cms3.genbs_p4().size();i++) {
        if(fabs(cms3.genbs_motherid().at(i))==6){
         bs.push_back(i);
       }
      }
    //find b and bbar
      if(bs.size()<2) continue;//skipping events with less than 2 b quarks.
      if(genbs_charge()[bs.at(0)]>0) {bquark = genbs_p4().at(bs.at(0)); bbarquark = genbs_p4().at(bs.at(1));}
      else {bquark = genbs_p4().at(bs.at(1)); bbarquark = genbs_p4().at(bs.at(0));}

    //find electrons coming from W plus
      for (int i=0;i<cms3.genels_p4().size();i++) {
        if(genels_isfromt()[i] && genels_motherid()[i] == 24)   els_fromWplus.push_back(i);         
        if(genels_isfromt()[i] && genels_motherid()[i] == -24)  els_fromWminus.push_back(i);         
      }
    //find muons coming from W plus
      for (int i=0;i<cms3.genmus_p4().size();i++) {
         if(genmus_isfromt()[i] && genmus_motherid()[i] == 24)  mus_fromWplus.push_back(i);
         if(genmus_isfromt()[i] && genmus_motherid()[i] == -24) mus_fromWminus.push_back(i);
      }
    //find taus coming from W plus
      for (int i=0;i<cms3.gentaus_p4().size();i++) {
         if(gentaus_isfromt()[i] && gentaus_motherid()[i] == 24)     taus_fromWplus.push_back(i);
         if(gentaus_isfromt()[i] && gentaus_motherid()[i] == -24)  taus_fromWminus.push_back(i);
      }
    //find neutrinos.
     for(unsigned int n = 0; n< gennus_id().size(); ++n){
        if(abs(gennus_id()[n])!=12&&abs(gennus_id()[n])!=14&&abs(gennus_id()[n])!=16) continue;
        if(abs(gennus_motherid()[n])!=24) continue;
        if(abs(gennus_id()[n])==12&&gennus_motherid()[n]==24) nuels_fromWplus.push_back(n); 
        if(abs(gennus_id()[n])==12&&gennus_motherid()[n]==-24) nuels_fromWminus.push_back(n); 
        if(abs(gennus_id()[n])==14&&gennus_motherid()[n]==24) numus_fromWplus.push_back(n); 
        if(abs(gennus_id()[n])==14&&gennus_motherid()[n]==-24) numus_fromWminus.push_back(n); 
        if(abs(gennus_id()[n])==16&&gennus_motherid()[n]==24) nutaus_fromWplus.push_back(n); 
        if(abs(gennus_id()[n])==16&&gennus_motherid()[n]==-24) nutaus_fromWminus.push_back(n); 
        ws_lep.push_back(gennus_motherp4()[n]);
//        bool Wtop = false;
//        for(unsigned int m = 0; m< gents_id().size(); ++m){
//          //t -> b W+, tbar -> bbar W-, i.e. 6->5 24, -6 -> -5 -24
//        if(gennus_motherid()[n]==24&&gents_id()[m]==6) Wtop = true;
//        if(gennus_motherid()[n]==-24&&gents_id()[m]==-6) Wtop = true;
//         }
//        fromTop.push_back(Wtop);
     }
      for (int i=0;i<cms3.genqs_p4().size();i++) {
        if(genqs_motherid()[i]==24) qs_fromWplus.push_back(i);
        if(genqs_motherid()[i]==-24) qs_fromWminus.push_back(i);
      }
//      if(qs_fromWminus.size()!=2&&qs_fromWplus.size()!=2) cout<<"qs from w plus :"<<qs_fromWplus.size()<< "qs from w minus"<<qs_fromWminus.size()<<endl;
      for (int i=0;i<cms3.gents_p4().size();i++) {
        if(fabs(cms3.gents_motherid().at(i)) == 1000003 || fabs(cms3.gents_motherid().at(i)) == 2000003){
          ts.push_back(i);
        }
      }

      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }

      //find what decay mode is this.
      bool isHad(false); bool is1l(false); bool is2l(false);
      if((qs_fromWplus.size()+qs_fromWminus.size())==4) isHad = true;
      if((qs_fromWplus.size()==2 || qs_fromWminus.size()==2) && (els_fromWplus.size()+els_fromWminus.size()+mus_fromWplus.size()+mus_fromWminus.size()+taus_fromWplus.size()+taus_fromWminus.size())==1) is1l = true;
      if((els_fromWplus.size()+els_fromWminus.size()+mus_fromWplus.size()+mus_fromWminus.size()+taus_fromWplus.size()+taus_fromWminus.size())==2) is2l = true;

      if(!isHad&&!is1l&&!is2l) cout<<"how did the top quark decay?!" << endl;

    //find pairs of t's
      LorentzVector bfrom_Whad;
      LorentzVector bfrom_Wlep;
      LorentzVector Whad;
      LorentzVector Wlep;
      LorentzVector lep_fromW;
      vector<LorentzVector> qs_fromW;
      if(is1l) {
       if(qs_fromWplus.size()==2) { Whad = (genqs_p4().at(qs_fromWplus.at(0))+genqs_p4().at(qs_fromWplus.at(1))); bfrom_Whad = bbarquark; qs_fromW.push_back(genqs_p4().at(qs_fromWplus.at(0)));qs_fromW.push_back(genqs_p4().at(qs_fromWplus.at(1)));}
       if(qs_fromWminus.size()==2) { Whad = (genqs_p4().at(qs_fromWminus.at(0))+genqs_p4().at(qs_fromWminus.at(1))); bfrom_Whad = bquark; qs_fromW.push_back(genqs_p4().at(qs_fromWminus.at(0)));qs_fromW.push_back(genqs_p4().at(qs_fromWminus.at(1)));}
      //wplus
       if(els_fromWplus.size() && nuels_fromWplus.size()) {Wlep = (genels_p4().at(els_fromWplus.at(0))+genels_p4().at(nuels_fromWplus.at(0)));bfrom_Wlep = bbarquark;lep_fromW = genels_p4().at(els_fromWplus.at(0));}     
       if(mus_fromWplus.size()&& numus_fromWplus.size()) {Wlep = (genmus_p4().at(mus_fromWplus.at(0))+genmus_p4().at(numus_fromWplus.at(0)));bfrom_Wlep = bbarquark;lep_fromW = genmus_p4().at(mus_fromWplus.at(0));}     
       if(taus_fromWplus.size()&& nutaus_fromWplus.size()) {Wlep = (gentaus_p4().at(taus_fromWplus.at(0))+gentaus_p4().at(nutaus_fromWplus.at(0)));bfrom_Wlep = bbarquark;lep_fromW = gentaus_p4().at(taus_fromWplus.at(0));}     
      //wminus
       if(els_fromWminus.size() && nuels_fromWminus.size()) {Wlep = (genels_p4().at(els_fromWminus.at(0))+gennus_p4().at(nuels_fromWminus.at(0)));bfrom_Wlep = bquark;lep_fromW = genels_p4().at(els_fromWminus.at(0));}     
       if(mus_fromWminus.size() && numus_fromWminus.size()) {Wlep = (genmus_p4().at(mus_fromWminus.at(0))+gennus_p4().at(numus_fromWminus.at(0)));bfrom_Wlep = bquark;lep_fromW = genmus_p4().at(mus_fromWminus.at(0));}
       if(taus_fromWminus.size() && nutaus_fromWminus.size()) {Wlep = (gentaus_p4().at(taus_fromWminus.at(0))+gennus_p4().at(nutaus_fromWminus.at(0)));bfrom_Wlep = bquark;lep_fromW = gentaus_p4().at(taus_fromWminus.at(0));}     
      } 
     
      // gen level 
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gentslv = cms3.gents_p4();
      value["genel_n"]=cms3.genels_p4().size();  
      value["genmu_n"]=cms3.genmus_p4().size();  
      value["gentau_n"]=cms3.gentaus_p4().size();  
      value["genb_n"]=cms3.genbs_p4().size();  
      value["genqs_n"]=cms3.genqs_p4().size();  
      value["genstop_n"]=cms3.genstop_p4().size();  
      value["gentop_n"]=ts.size();  

      if(is1l && (Wlep+bfrom_Wlep).pt()>10) {
      value["gen_wb"] = JetUtil::deltaR(Wlep,bfrom_Wlep);
      value["gentop_pt"]=(Wlep+bfrom_Wlep).pt();  
      value["gen_dRtt"] = JetUtil::deltaR(Wlep+bfrom_Wlep,Whad+bfrom_Whad);
      value["gen_dRlt"] = JetUtil::deltaR(lep_fromW,Whad+bfrom_Whad);
      value["gen_dRlw"] = JetUtil::deltaR(lep_fromW,Whad);
      }

    bool noMET(true),MET100(false),MET300(false);

    if(pfmet()>300) MET300 = true;
  
    for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
	if(noMET){
	  mname= "noMET"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ]);
        }
	if(MET300){
	  mname= "MET300"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ]);
        }
      }

     for(unsigned int i = 0; i<histonames_vec.size(); ++i){
	string d = "_";
	string mname;
	if(noMET){
	  mname= "noMET"+d+histonames_vec[i]+d+samplename;
	  for (int j=0;j < value_vec[histonames_vec[i]].size();j++){
          histos_vec[mname]->Fill(value_vec[histonames_vec[i]].at(j));
         }
       }
     }

    // fill in 2D hists
	string d = "_";
	if(noMET && is1l){
	  string mname ="";
	  string histname ="";
          //fill each histogram
	  histname ="gentop_pt_vs_dRWlep_b";
	  mname= "noMET"+d+histname+d+samplename;
          histos2d[mname]->Fill((Wlep+bfrom_Wlep).pt(),JetUtil::deltaR(Wlep,bfrom_Wlep));

	  histname ="gentop_pt_vs_dRWhad_b";
	  mname= "noMET"+d+histname+d+samplename;
          histos2d[mname]->Fill((Whad+bfrom_Whad).pt(),JetUtil::deltaR(Whad,bfrom_Whad));

	  histname ="genw_pt_vs_dRqq";
	  mname= "noMET"+d+histname+d+samplename;
          histos2d[mname]->Fill(Whad.pt(),JetUtil::deltaR(qs_fromW.at(0),qs_fromW.at(1)));
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
  //set the overflow bin
  for(map<string,TH1F*>::iterator h=histos.begin(); h!=histos.end();++h){
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetBinContent(h->second->GetNbinsX() )+ h->second->GetBinContent(h->second->GetNbinsX()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), sqrt(pow(h->second->GetBinError(h->second->GetNbinsX() ),2)+pow(h->second->GetBinError(h->second->GetNbinsX()+1),2) ) );
   }
  for(map<string,TH1F*>::iterator h=histos_vec.begin(); h!=histos_vec.end();++h){
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
    /*h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			     h->second->GetBinContent(h->second->GetNbinsX(),h->second->GetNbinsY() )+ h->second->GetBinContent(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			   sqrt(pow(h->second->GetBinError(h->second->GetNbinsX(),h->second->GetNbinsY() ),2)+
				pow(h->second->GetBinError(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1),2) ) );
    */
 }
  string filename = "rootfiles/CutHistos/FatJetPlots/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
  for(map<string,TH1F*>::iterator h=    histos_vec.begin(); h!=    histos_vec.end();++h) h->second->Write();
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
