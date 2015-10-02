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
#include "StopBabies09042015.cc"

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

  histonames_vec.push_back("genb_status");     histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("gentop_status");   histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genstop_status");  histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genel_status");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genmu_status");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("gentau_status");   histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genqs_status");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);

  histonames_vec.push_back("genb_motherid");     histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("gentop_motherid");   histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genstop_motherid");  histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genel_motherid");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genmu_motherid");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("gentau_motherid");   histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);
  histonames_vec.push_back("genqs_motherid");    histobinn_vec.push_back(1000); histobinl_vec.push_back(0.); histobinu_vec.push_back(1000.);

  histonames.push_back("gentop_pt");  histobinn.push_back(85); histobinl.push_back(0.); histobinu.push_back(1700.);
  histonames.push_back("genz_pt");    histobinn.push_back(25); histobinl.push_back(100.); histobinu.push_back(600.);
  histonames.push_back("geng_pt");    histobinn.push_back(85); histobinl.push_back(0.); histobinu.push_back(1700.);
  histonames.push_back("gamma_pt");   histobinn.push_back(25); histobinl.push_back(100.); histobinu.push_back(600.);
  histonames.push_back("gamma_genmatched_pt");   histobinn.push_back(25); histobinl.push_back(100.); histobinu.push_back(600.);
  histonames.push_back("genz_mass");  histobinn.push_back(85); histobinl.push_back(0.); histobinu.push_back(1700.);
  histonames.push_back("gentop_eta"); histobinn.push_back(10); histobinl.push_back(0.); histobinu.push_back(5.);
  histonames.push_back("gen_wb");     histobinn.push_back(100); histobinl.push_back(0.); histobinu.push_back(5.);
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

  histonames2d.push_back("gentop_pt_vs_dRWb"); 
  histobinnX.push_back(100); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(100); histobinlY.push_back(0.); histobinuY.push_back(3.);

  //const unsigned int poststringsize = 3;
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
      StopBabies09042015::progress( nEventsTotal, nEventsChain );

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
      float weight = cms3.scale1fb()*5.;
      int NLeps = cms3.ngoodleps();
      string ds = cms3.dataset();
      float MET = cms3.pfmet();
      float METPhi = cms3.pfmet_phi();
      float METx = MET*TMath::Cos(METPhi);
      float METy = MET*TMath::Sin(METPhi);
      float genMET = cms3.genmet();
      float genMETPhi = cms3.genmet_phi();
      float genMETx = genMET*TMath::Cos(genMETPhi);
      float genMETy = genMET*TMath::Sin(genMETPhi);
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

      vector<LorentzVector> nus_fromZ;
      vector<LorentzVector> gs_fromT;
      vector<LorentzVector> els_fromZ;
      vector<LorentzVector> mus_fromZ;
      vector<LorentzVector> taus_fromZ;
      vector<LorentzVector> reco_gamma;
      vector<LorentzVector> reco_gamma_genmatched;

      vector<LorentzVector> ws_lep;
      vector<LorentzVector> ws_had;
      vector<int> gen_els;
      vector<bool> fromTop;
      LorentzVector bquark;
      LorentzVector bbarquark;
      LorentzVector zboson;

      for (int i=0;i<cms3.genbs_p4().size();i++) {
        if(fabs(cms3.genbs_motherid().at(i))==6 && cms3.genbs_p4()[i].pt() > 30){
         bs.push_back(i);
       }
      }
      for (int i=0;i<cms3.genphs_p4().size();i++) {
        if (!(cms3.genphs_status().at(i)==23 || cms3.genphs_status().at(i)==1))  continue;
        if (fabs(cms3.genphs_simplemotherid().at(i)) > 22  && cms3.genphs_simplemotherid().at(i)!=2212) continue;  
        if (fabs(genphs_p4().at(i).eta())>2.5) continue;
        gs_fromT.push_back(genphs_p4().at(i));
      }
     for (int i=0;i<ph_p4().size();i++) {
        if(ph_chiso().at(i) < 2.5 && ph_idCutBased().at(i)) reco_gamma.push_back(ph_p4().at(i));
        if(ph_chiso().at(i) < 2.5 && ph_idCutBased().at(i) && ph_mcMatchId().at(i) > 0) {
         reco_gamma_genmatched.push_back(ph_p4().at(i));
        }
      }
    //find b and bbar
//      if(bs.size()<2) continue;//skipping events with less than 2 b quarks.
//      if(genbs_charge()[bs.at(0)]>0) {bquark = genbs_p4().at(bs.at(0)); bbarquark = genbs_p4().at(bs.at(1));}
//      else {bquark = genbs_p4().at(bs.at(1)); bbarquark = genbs_p4().at(bs.at(0));}
    //find electrons coming from W and Z
      for (int i=0;i<cms3.genels_p4().size();i++) {
        if(genels_isfromt()[i] && genels_motherid()[i] == 24)   els_fromWplus.push_back(i);         
        if(genels_isfromt()[i] && genels_motherid()[i] == -24)  els_fromWminus.push_back(i);         
        if(genels_motherid()[i] == 23 && genels_p4().at(i).pt() > 20)  els_fromZ.push_back(genels_p4().at(i));         
      }
    //find muons coming from W and Z
      for (int i=0;i<cms3.genmus_p4().size();i++) {
         if(genmus_isfromt()[i] && genmus_motherid()[i] == 24)  mus_fromWplus.push_back(i);
         if(genmus_isfromt()[i] && genmus_motherid()[i] == -24) mus_fromWminus.push_back(i);
         if(genmus_motherid()[i] == 23 && genmus_p4().at(i).pt() > 20)  mus_fromZ.push_back(genmus_p4().at(i));         
      }
    //find taus coming from W and Z
      for (int i=0;i<cms3.gentaus_p4().size();i++) {
         if(gentaus_isfromt()[i] && gentaus_motherid()[i] == 24)   taus_fromWplus.push_back(i);
         if(gentaus_isfromt()[i] && gentaus_motherid()[i] == -24)  taus_fromWminus.push_back(i);
         if(gentaus_motherid()[i] == 23 && gentaus_p4().at(i).pt() > 20 )  taus_fromZ.push_back(gentaus_p4().at(i));         
      }
//   find the neutrinos from Z
      for(unsigned int n = 0; n< gennus_id().size(); ++n){
        if(abs(gennus_id()[n])!=12&&abs(gennus_id()[n])!=14&&abs(gennus_id()[n])!=16) continue;
        if(abs(gennus_motherid()[n])==23) nus_fromZ.push_back(gennus_p4().at(n));
     }
    //if(nus_fromZ.size()<2) continue;
    //cout<<"n neutrinos"<< nus_fromZ.size() <<endl;
    //find neutrinos from w decay.
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
     }
      for (int i=0;i<cms3.genqs_p4().size();i++) {
        if(genqs_motherid()[i]==24 && cms3.genqs_p4()[i].pt()>30) qs_fromWplus.push_back(i);
        if(genqs_motherid()[i]==-24&& cms3.genqs_p4()[i].pt()>30) qs_fromWminus.push_back(i);
      }

      for (int i=0;i<cms3.gents_p4().size();i++) {
        if(fabs(cms3.gents_motherid().at(i)) == 1000003 || fabs(cms3.gents_motherid().at(i)) == 2000003){
          ts.push_back(i);
        }
      }
     // fill in genMET vector
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv_mod;
      genmetlv.SetPxPyPzE(genMETx,genMETy,0.,genMET);
      if(gs_fromT.size()) genmetlv_mod = gs_fromT.at(0)+genmetlv; 
 
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
      //find what decay mode is this.
      bool isHad; bool is1l; bool is2l;
      if((qs_fromWplus.size()+qs_fromWminus.size())==4) isHad = true;
      if((qs_fromWplus.size()+qs_fromWminus.size())==2 && (els_fromWplus.size()+els_fromWminus.size()+mus_fromWplus.size()+mus_fromWminus.size()+taus_fromWplus.size()+taus_fromWminus.size())==1) is1l = true;
      if((els_fromWplus.size()+els_fromWminus.size()+mus_fromWplus.size()+mus_fromWminus.size()+taus_fromWplus.size()+taus_fromWminus.size())==2) is2l = true;
      if(!isHad&&!is1l&&!is2l) cout<<"how did the top quark decay?!" << endl;
    //find pairs of t's
      LorentzVector bfrom_Whad;
      LorentzVector bfrom_Wlep;
      LorentzVector Whad;
      LorentzVector Wlep;

      if(is1l){
       if(qs_fromWplus.size()==2)  {Whad = (genqs_p4().at(qs_fromWplus.at(0))+genqs_p4().at(qs_fromWplus.at(1))); bfrom_Whad = bbarquark;}
       if(qs_fromWminus.size()==2) {Whad = (genqs_p4().at(qs_fromWminus.at(0))+genqs_p4().at(qs_fromWminus.at(1))); bfrom_Whad = bquark;}
      //wplus
       if(els_fromWplus.size() && nuels_fromWplus.size()) {Wlep = (genels_p4().at(els_fromWplus.at(0))+gennus_p4().at(nuels_fromWplus.at(0)));bfrom_Wlep = bbarquark;}     
       if(mus_fromWplus.size()&& numus_fromWplus.size()) {Wlep = (genmus_p4().at(mus_fromWplus.at(0))+gennus_p4().at(numus_fromWplus.at(0)));bfrom_Wlep = bbarquark;}     
       if(taus_fromWplus.size()&& nutaus_fromWplus.size()) {Wlep = (gentaus_p4().at(taus_fromWplus.at(0))+gennus_p4().at(nutaus_fromWplus.at(0)));bfrom_Wlep = bbarquark;}     
      //wminus
       if(els_fromWminus.size() && nuels_fromWminus.size()) {Wlep = (genels_p4().at(els_fromWminus.at(0))+gennus_p4().at(nuels_fromWminus.at(0)));bfrom_Wlep = bquark;}     
       if(mus_fromWminus.size() && numus_fromWminus.size()) {Wlep = (genmus_p4().at(mus_fromWminus.at(0))+gennus_p4().at(numus_fromWminus.at(0)));bfrom_Wlep = bquark;}     
       if(taus_fromWminus.size()&& nutaus_fromWminus.size()) {Wlep = (gentaus_p4().at(taus_fromWminus.at(0))+gennus_p4().at(nutaus_fromWminus.at(0)));bfrom_Wlep = bquark;}     
      } 
      // gen level 
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > gentslv = cms3.gents_p4();
      value["genel_n"]=cms3.genels_p4().size();  
      //value["genmu_n"]=cms3.genmus_p4().size();  
      //value["gentau_n"]=cms3.gentaus_p4().size();  
      //value["genb_n"]=cms3.genbs_p4().size();  
      //value["genqs_n"]=cms3.genqs_p4().size();  
      value["genstop_n"]=cms3.genstop_p4().size();  
      value["gentop_n"]=ts.size();  
      if(is1l && (Wlep+bfrom_Wlep).pt()>10) {
      value["gen_wb"] = JetUtil::deltaR(Wlep,bfrom_Wlep);
      value["gentop_pt"]=(Wlep+bfrom_Wlep).pt();  
      value2d["gentop_pt_vs_dRWb"].push_back((Wlep+bfrom_Wlep).pt());
      value2d["gentop_pt_vs_dRWb"].push_back(JetUtil::deltaR(Wlep,bfrom_Wlep));
      }

    if(els_fromZ.size()==2) zboson = els_fromZ.at(0)+els_fromZ.at(1); 
    if(mus_fromZ.size()==2) zboson = mus_fromZ.at(0)+mus_fromZ.at(1); 
    if(taus_fromZ.size()==2) zboson = taus_fromZ.at(0)+taus_fromZ.at(1); 
    if(nus_fromZ.size()==2) zboson = nus_fromZ.at(0)+nus_fromZ.at(1); 
    //if(zboson.mass() > 40) value["genz_pt"] = zboson.pt();
    if(is1l){
    if(gs_fromT.size()) { if(gs_fromT.at(0).pt()>100) value["genz_pt"] = gs_fromT.at(0).pt(); }
    if(zboson.mass() > 40) value["genz_mass"] = zboson.mass();
    if(gs_fromT.size())  value["geng_pt"] = gs_fromT.at(0).pt(); 
    if(reco_gamma.size()) { if(reco_gamma.at(0).pt()>100) value["gamma_pt"] = reco_gamma.at(0).pt();};
    if(reco_gamma_genmatched.size()) {if(reco_gamma_genmatched.at(0).pt()>100) value["gamma_genmatched_pt"] = reco_gamma_genmatched.at(0).pt();}
    }
    bool noMET(true),MET100(false),MET300(false);

    if(genmetlv_mod.pt() > 300) MET300 = true;
     
    for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
	if(noMET){
	  mname= "noMET"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
        }
	if(MET300){
	  mname= "MET300"+d+histonames[i]+d+samplename;
	  histos[mname]->Fill(value[histonames[i] ],weight);
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

     for(unsigned int i = 0; i<histonames2d.size(); ++i){
	string d = "_";
	string mname;
	if(noMET){
	  mname= "noMET"+d+histonames2d[i]+d+samplename;
          //histos2d[mname]->Fill(value2d[histonames2d[i]].at(0),value2d[histonames2d[i]].at(1));
          //histos2d[mname]->Fill(value2d[histonames2d[i]].at(1),value2d[histonames2d[i]].at(0));
          histos2d[mname]->Fill((Wlep+bfrom_Wlep).pt(),JetUtil::deltaR(Wlep,bfrom_Wlep));
       }
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
  string filename = "../rootfiles/CutHistos/FatJetPlots/"+skimFilePrefix+".root";
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
