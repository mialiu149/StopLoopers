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
#include "TRandom.h"

// CMS3
#include "StopBabiesRun2_06252015.cc"

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

inline float calculateMt(const LorentzVector p4, double met, double met_phi)
{  
  float phi1 = p4.Phi();
  float phi2 = met_phi;
  float Et1  = p4.Et();
  float Et2  = met;
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

int ScanChain( TChain* chain, bool fast = true, int nEvents = -1, string skimFilePrefix = "test") {
// Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");
// Example Histograms
  TRandom myRandom(12345);
  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
//Delare histonames and bins.

  map<string, TH1F*> histos; 
  vector<string> histonames; histonames.clear();
  vector<int> histobinn; histobinn.clear();
  vector<double> histobinl; histobinl.clear();
  vector<double> histobinu; histobinu.clear();
  map<string, float> value;

  map<string, TH2F*> histos2d;
  vector<string> histonames2d; histonames2d.clear();
  vector<int> histobinnX; histobinnX.clear();
  vector<double> histobinlX; histobinlX.clear();
  vector<double> histobinuX; histobinuX.clear();
  vector<int> histobinnY; histobinnY.clear();
  vector<double> histobinlY; histobinlY.clear();
  vector<double> histobinuY; histobinuY.clear();
  map<string, float> value2d;

// fill in histogram names and binning
  histonames.push_back("MT_RecoLep_pfMET");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_pfMET2Xworse");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_pfMET2Xworse_genMET");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_genMET");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_genMET_smeared");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_genMET_smeared2X");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_genMET_BWsmeared");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_RecoLep_genMET_BWsmeared2X");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_genLep_pfMET");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("MT_genLep_genMET");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("pfMETMinusGenMET_x");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("pfMETMinusGenMET_y");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("smearedGenMETMinusGenMET_x");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("smearedGenMETMinusGenMET_y");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("BWsmearedGenMETMinusGenMET_x");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("BWsmearedGenMETMinusGenMET_y");   histobinn.push_back(80); histobinl.push_back(-800); histobinu.push_back(800.);
  histonames.push_back("pfMETMinusGenMET_x_overGenMET");   histobinn.push_back(80); histobinl.push_back(-3); histobinu.push_back(3.);
  histonames.push_back("pfMETMinusGenMET_y_overGenMET");   histobinn.push_back(80); histobinl.push_back(-3); histobinu.push_back(3);
  histonames.push_back("myMT");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);
  histonames.push_back("myMTn");   histobinn.push_back(80); histobinl.push_back(0.); histobinu.push_back(800.);

  histonames2d.push_back("gentop_pt_vs_dRWb"); 
  histobinnX.push_back(10); histobinlX.push_back(0.); histobinuX.push_back(1000.);
  histobinnY.push_back(10); histobinlY.push_back(0.); histobinuY.push_back(3.);
  const unsigned int poststringsize = 5;
  //string poststring[poststringsize] = {"W_selection","Preselection","MET300_noMT","MET100_MT120","noMT","MT80","MT100","MT125","MT150","SR","SR200","SR200MT2W"};
  //string poststring[poststringsize] = {"W_selection","Wenu_selection","Wmunu_selection","W_selection_notau","Preselection","MET300_noMT","SR","SR200","SR200MT2W"};
  string poststring[poststringsize] = {"W_selection","Wenu_selection","Wmunu_selection","W_selection_notau","Preselection"};

// book Cutflow hists
    map<string, TH1D*> histos_cutflow;
    vector<string> histonames_cutflow; histonames_cutflow.clear();
    histonames_cutflow.push_back("NEventsPerSignalRegion");
    histonames_cutflow.push_back("NEventsPerSignalRegion_2Xworse");
    histonames_cutflow.push_back("NEventsPerSignalRegion_2Xworse_genMET");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT150");
    histonames_cutflow.push_back("MT150_2Xworse_wrt_nominal");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT175");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT150_2Xworse");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT150_10Xworse");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT150_2Xworse_genMET");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT150_10Xworse_genMET");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT175_2Xworse_genMET");
    histonames_cutflow.push_back("NEventsPerSignalRegionMT175_10Xworse_genMET");

//initialize histos for ttbar
  for(unsigned int b = 0; b<3; ++b){
    string samplename = skimFilePrefix;
    if(skimFilePrefix.find("ttbar")==std::string::npos&&b>0) continue;
    if(skimFilePrefix.find("ttbar")!=std::string::npos&&b==0) samplename = "TTbar1l";
    if(skimFilePrefix.find("ttbar")!=std::string::npos&&b==1) samplename = "TTbar2l";
    if(skimFilePrefix.find("ttbar")!=std::string::npos&&b==2) samplename = "TTbarH";

    for(unsigned int i = 0;i<histonames_cutflow.size();i++){
    string histoname = histonames_cutflow.at(i)+"_"+samplename;
    histos_cutflow[histoname] = new TH1D(histoname.c_str(),"",8,0,8);
    histos_cutflow[histoname]->Sumw2(); histos_cutflow[histoname]->SetDirectory(rootdir);
    histos_cutflow[histoname]->GetYaxis()->SetTitle("Events");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(1,"MET>150,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(2,"MET>200,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(3,"MET>250,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(4,"MET>300,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(5,"MET>150,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(6,"MET>200,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(7,"MET>250,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(8,"MET>300,hDM");
    }
  }
//initialize histos for wjets
  for(unsigned int b = 0; b<2; ++b){
    string samplename = skimFilePrefix;

    if(skimFilePrefix.find("WJets")==std::string::npos&&b>0) continue;
    if(skimFilePrefix.find("WJets")!=std::string::npos&&b==0) samplename = "WJetsLight";
    if(skimFilePrefix.find("WJets")!=std::string::npos&&b==1) samplename = "WJetsHeavy";

    for(unsigned int i = 0;i<histonames_cutflow.size();i++){
    string histoname = histonames_cutflow.at(i)+"_"+samplename;
    histos_cutflow[histoname] = new TH1D(histoname.c_str(),"",8,0,8);
    histos_cutflow[histoname]->Sumw2(); histos_cutflow[histoname]->SetDirectory(rootdir);
    histos_cutflow[histoname]->GetYaxis()->SetTitle("Events");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(1,"MET>150,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(2,"MET>200,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(3,"MET>250,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(4,"MET>300,lDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(5,"MET>150,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(6,"MET>200,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(7,"MET>250,hDM");
    histos_cutflow[histoname]->GetXaxis()->SetBinLabel(8,"MET>300,hDM");
    }
}
//book 1D hists
for(unsigned int c = 0; c<poststringsize; ++c){
   for(unsigned int i = 0; i<histonames.size(); ++i){
      string mapname;
      mapname = poststring[c] + "_" + histonames[i];
      if(histos.count(mapname) == 0 ) histos[mapname] = new TH1F(mapname.c_str(), "", histobinn[i], histobinl[i], histobinu[i]);
       histos[mapname]->Sumw2(); histos[mapname]->SetDirectory(rootdir);
       if(c==0) value[histonames[i] ] = -1;
      }
  }
//book 2D hists
 for(unsigned int c = 0; c<poststringsize; ++c){
    for(unsigned int i = 0; i<histonames2d.size(); ++i){ 
      string mapname;
      mapname = poststring[c] + "_" + histonames2d[i];
      if(histos2d.count(mapname) == 0 ) histos2d[mapname] = new TH2F(mapname.c_str(), "", histobinnX[i], histobinlX[i], histobinuX[i],histobinnY[i], histobinlY[i], histobinuY[i]);
      histos2d[mapname]->Sumw2(); histos2d[mapname]->SetDirectory(rootdir);
   }
}	//if(c==0) value[histonames2D[i] ] = -1;
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
      StopBabiesRun2_06252015::progress( nEventsTotal, nEventsChain );

      string samplename = skimFilePrefix;
      if(skimFilePrefix.find("ttbar")!=std::string::npos){
	if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==2) samplename = "TTbar2l";
	else if(cms3.gen_nfromtmus_() + cms3.gen_nfromtels_() + cms3.gen_nfromttaus_() ==1) samplename = "TTbar1l";
	else samplename = "TTbarH";
      }

     int nbs = 0;
     int ncs = 0;
     if(skimFilePrefix.find("WJets")!=std::string::npos){
     for(unsigned int j=0; j<genbs_p4().size();++j){
        if(TMath::Abs(genbs_motherid()[j])<26||(TMath::Abs(genbs_motherid()[j])>1000000&&TMath::Abs(genbs_motherid()[j])<2000016) ) nbs++;
      }
      for(unsigned int j=0; j<genqs_p4().size();++j){
        if(TMath::Abs(genqs_motherid()[j])<26||(TMath::Abs(genqs_motherid()[j])>1000000&&TMath::Abs(genqs_motherid()[j])<2000016) && fabs(genqs_id()[j])==4) ncs++;}
     if(nbs!=0) samplename = "WJetsHeavy";
     else samplename = "WJetsLight";
//	cout << "f " << __LINE__ << endl;
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

      int NGLeps = 0;
      int NSLeps = 0;
      int NGJets = 0;
      int NGBJets = 0;

      //find one good lepton.
      int l1=-1;
      if(lep1ismu){
	if(lep1pt>20&&fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
	//if(fabs(lep1eta)<2.1&&fabs(cms3.lep1_d0())<0.02&&fabs(cms3.lep1_dz())<0.1&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      } else if (lep1isel){
	if(lep1pt>20&&fabs(lep1eta)<2.1&&cms3.lep1_is_phys14_medium_noIso()&&cms3.lep1_miniRelIsoDB()<0.1) {++NSLeps; l1 = 1; }
      }
      NGLeps = NLeps;
      if(firstGoodVtxIdx() != 0) continue; ++c1vtx; 
      if(NGLeps!=1) continue; ++c1l;

//	cout << "f " << __LINE__ << endl;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > leplv;
      if(l1==1) leplv = lep1lv;
      else if(l1==2) leplv = lep2lv;
      //count jets
      vector<int> jind;
      vector<int> bind;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > jsumlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > jetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > boostjetlv;
      vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > bjetlv;
      vector<float> btag;
//	cout << "f " << __LINE__ << endl;

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
      // gen level 
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genleplv; 
      // find gen lepton
//      if(skimFilePrefix=="WJetsToLNu_HT400To600_25ns"||skimFilePrefix=="WJetsToLNu_HT100to200_25ns"||skimFilePrefix=="WJetsToLNu_HT200To400_25ns"||skimFilePrefix=="WJetsToLNu_HT600ToInf_25ns")
      if(skimFilePrefix=="WJets")  
     {
      for  (int i=0;i<genels_p4().size();i++){
      if(abs(genels_motherid().at(i))==24) genleplv = genels_p4().at(i);
      }
      for  (int i=0;i<genmus_p4().size();i++){
      if(abs(genmus_motherid().at(i))==24) genleplv = genmus_p4().at(i);
      }
      } 
      else
      {
      for (int i=0;i<genels_p4().size();i++){
      if (genels_isfromt().at(i)) genleplv = genels_p4().at(i);
      }          
      for (int i=0;i<genmus_p4().size();i++){
      if (genmus_isfromt().at(i)) genleplv = genmus_p4().at(i);
      }         
      }
      bool isOneLepEvent=true;
      if(samplename == "TTbar2l") isOneLepEvent = false;
     //	cout << "f " << __LINE__ << endl;
     // find good w's from top
      int ng =0;
      vector<LorentzVector> W;
      vector<int> Wid;
      vector<bool> fromTop;
      vector<int> tid;
      for(unsigned int n = 0; n< gents_id().size(); ++n){
        tid.push_back(gents_id()[n]);
      }
      for(unsigned int n = 0; n< gennus_id().size(); ++n){
        if(abs(gennus_id()[n])!=12&&abs(gennus_id()[n])!=14&&abs(gennus_id()[n])!=16) continue;
        if(abs(gennus_motherid()[n])!=24) continue;
        W.push_back(gennus_motherp4()[n]);
        Wid.push_back(gennus_motherid()[n]);
        bool Wtop = false;
        for(unsigned int m = 0; m< gents_id().size(); ++m){
          //t -> b W+, tbar -> bbar W-, i.e. 6->5 24, -6 -> -5 -24
        if(gennus_motherid()[n]==24&&gents_id()[m]==6) Wtop = true;
        if(gennus_motherid()[n]==-24&&gents_id()[m]==-6) Wtop = true;
         }
        fromTop.push_back(Wtop);
       ++ng;
      }
//	cout << "f " << __LINE__ << endl;
      float thegenlepDR = 9999.;
      LorentzVector thegenlep;
      LorentzVector thegennu;
      LorentzVector thegentaunu;
      int thegenlepid;
      bool thegenlepfromtau = false; 

       for(unsigned int n = 0; n< genmus_p4().size(); ++n) {
        if(genmus_p4()[n].Pt()<1) continue;
        float DR = JetUtil::deltaR(leplv,genmus_p4()[n]);
        if(DR<thegenlepDR&&abs(genmus_motherid().at(n))==24){
          thegenlepDR = DR;
          thegenlep = genmus_p4()[n];
          thegenlepid = genmus_id()[n];
        }
      }
      for(unsigned int n = 0; n< genels_p4().size(); ++n) {
        if(genels_p4()[n].Pt()<1) continue;
        
        float DR = JetUtil::deltaR(leplv,genels_p4()[n]);
        if(DR<thegenlepDR&&abs(genels_motherid().at(n))==24){
          thegenlepDR = DR;
          thegenlep = genels_p4()[n];
          thegenlepid = genels_id()[n];
        }
      }
       for(unsigned int n = 0; n< genleptau_mus_p4().size(); ++n) {
        if(genleptau_mus_p4()[n].Pt()<1) continue;
        float DR = JetUtil::deltaR(leplv,genleptau_mus_p4()[n]);
        if(DR<thegenlepDR){
          thegenlepDR = DR;
          thegenlep = genleptau_mus_p4()[n];
          thegenlepid = genleptau_mus_id()[n];
          thegenlepfromtau = true;
        }
      }
      for(unsigned int n = 0; n< genleptau_els_p4().size(); ++n) {
        if(genleptau_els_p4()[n].Pt()<1) continue;
        float DR = JetUtil::deltaR(leplv,genleptau_els_p4()[n]);
        if(DR<thegenlepDR){
          thegenlepDR = DR;
          thegenlep = genleptau_els_p4()[n];
          thegenlepid = genleptau_els_id()[n];
          thegenlepfromtau = true;
        }
      }

      for(unsigned int n = 0; n< gennus_id().size(); ++n){
        if(abs(gennus_id()[n])!=12&&abs(gennus_id()[n])!=14&&abs(gennus_id()[n])!=16) continue;
        if(thegenlepid==11&&gennus_id()[n]==-12) thegennu = gennus_p4()[n];
        if(thegenlepid==-11&&gennus_id()[n]==12) thegennu = gennus_p4()[n];
        if(thegenlepid==13&&gennus_id()[n]==-14) thegennu = gennus_p4()[n];
        if(thegenlepid==-13&&gennus_id()[n]==14) thegennu = gennus_p4()[n];
        if(thegenlepfromtau){
        if(thegenlepid==11&&gennus_id()[n]==16) thegentaunu = gennus_p4()[n];
        if(thegenlepid==-11&&gennus_id()[n]==-16) thegentaunu = gennus_p4()[n];
        if(thegenlepid==13&&gennus_id()[n]==16) thegentaunu = gennus_p4()[n];
        if(thegenlepid==-13&&gennus_id()[n]==-16) thegentaunu = gennus_p4()[n];
        }
      }
      if(thegenlepfromtau) thegennu = thegennu+thegentaunu;

      float myMT = sqrt(2.*genmet()*thegenlep.Et()*(1.-TMath::Cos(JetUtil::DeltaPhi(genmet_phi(),thegenlep.Phi() ) ) ) );
      float myMTn = sqrt(2.*thegennu.Et()*thegenlep.Et()*(1.-TMath::Cos(JetUtil::DeltaPhi(thegennu.Phi(),thegenlep.Phi() ) ) ) );

      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv2Xworse;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv2Xworse_genMET;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > metlv10Xworse_genMET;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv_smeared;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv_smeared2X;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv_BWsmeared;
      ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > genmetlv_BWsmeared2X;
      float thegenMETx = thegennu.Et()*TMath::Cos(thegennu.Phi());
      float thegenMETy = thegennu.Et()*TMath::Sin(thegennu.Phi());
      metlv.SetPxPyPzE(METx,METy,0.,MET);
      float pfMET2Xworse_METx = thegenMETx+2*(METx-thegenMETx);
      float pfMET2Xworse_METy = thegenMETy+2*(METy-thegenMETy);
      float pfMET2Xworse = sqrt((thegenMETx+2*(METx-thegenMETx))*(thegenMETx+2*(METx-thegenMETx))
                            +(thegenMETy+2*(METy-thegenMETy))*(thegenMETy+2*(METy-thegenMETy)));
      metlv2Xworse.SetPxPyPzE(thegenMETx+2*(METx-thegenMETx),thegenMETy+2*(METy-thegenMETy),0,pfMET2Xworse);
      float pfMET2Xworse_genMET_METx = genMETx+2*(METx-genMETx);
      float pfMET2Xworse_genMET_METy = genMETy+2*(METy-genMETy);
      float pfMET2Xworse_genMET = sqrt((genMETx+2*(METx-genMETx))*(genMETx+2*(METx-genMETx))
                            +(genMETy+2*(METy-genMETy))*(genMETy+2*(METy-genMETy)));
      metlv2Xworse_genMET.SetPxPyPzE(genMETx+2*(METx-genMETx),genMETy+2*(METy-genMETy),0,pfMET2Xworse_genMET);

      float pfMET10Xworse_genMET_METx = genMETx+10*(METx-genMETx);
      float pfMET10Xworse_genMET_METy = genMETy+10*(METy-genMETy);
      float pfMET10Xworse_genMET = sqrt((genMETx+10*(METx-genMETx))*(genMETx+10*(METx-genMETx))
                            +(genMETy+10*(METy-genMETy))*(genMETy+10*(METy-genMETy)));
      metlv10Xworse_genMET.SetPxPyPzE(genMETx+10*(METx-genMETx),genMETy+10*(METy-genMETy),0,pfMET10Xworse_genMET);
      genmetlv.SetPxPyPzE(thegenMETx,thegenMETy,0.,thegennu.Et());
     
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
      float MET_smeared = sqrt((thegenMETx+myRandom.Gaus(0,32))*(thegenMETx+myRandom.Gaus(0,32))+(thegenMETy+myRandom.Gaus(0,32))*(thegenMETy+myRandom.Gaus(0,32)));
      float MET_smeared2X = sqrt((thegenMETx+myRandom.Gaus(0,32*2))*(thegenMETx+myRandom.Gaus(0,32*2))+(thegenMETy+myRandom.Gaus(0,32*2))*(thegenMETy+myRandom.Gaus(0,32*2)));
      float MET_BWsmeared = sqrt((thegenMETx+myRandom.BreitWigner(0,32))*(thegenMETx+myRandom.BreitWigner(0,32))+(thegenMETy+myRandom.BreitWigner(0,32))*(thegenMETy+myRandom.BreitWigner(0,32)));
      float MET_BWsmeared2X = sqrt((thegenMETx+myRandom.BreitWigner(0,32*2))*(thegenMETx+myRandom.BreitWigner(0,32*2))+(thegenMETy+myRandom.BreitWigner(0,32*2))*(thegenMETy+myRandom.BreitWigner(0,32*2)));
      genmetlv_smeared2X.SetPxPyPzE(thegenMETx+myRandom.Gaus(0,32*2),thegenMETy+myRandom.Gaus(0,32*2),0,MET_smeared);
      genmetlv_smeared.SetPxPyPzE(thegenMETx+myRandom.Gaus(0,32),thegenMETy+myRandom.Gaus(0,32),0,MET_smeared);
      genmetlv_BWsmeared2X.SetPxPyPzE(thegenMETx+myRandom.BreitWigner(0,32*2),thegenMETy+myRandom.BreitWigner(0,32*2),0,MET_BWsmeared);
      genmetlv_BWsmeared.SetPxPyPzE(thegenMETx+myRandom.BreitWigner(0,32),thegenMETy+myRandom.BreitWigner(0,32),0,MET_BWsmeared);
  //only fill in histograms for one-lep events
//	cout << "f " << __LINE__ << endl;
      if(isOneLepEvent){
      value["MT_RecoLep_pfMET"]=getMT(leplv,metlv);
      value["MT_RecoLep_pfMET2Xworse"]=getMT(leplv,metlv2Xworse);
      value["MT_RecoLep_pfMET2Xworse_genMET"]=getMT(leplv,metlv2Xworse_genMET);
      value["MT_RecoLep_genMET"]=getMT(leplv,genmetlv);      
      value["MT_genLep_pfMET"]=getMT(thegenlep,metlv);
      value["MT_genLep_genMET"]=getMT(thegenlep,genmetlv);
      value["pfMETMinusGenMET_x"]=METx-genmetlv.Px(); 
      value["pfMETMinusGenMET_y"]=METy-genmetlv.Py();
      value["pfMETMinusGenMET_x_overGenMET"]=(METx-genmetlv.Px())/thegennu.Et(); 
      value["pfMETMinusGenMET_y_overGenMET"]=(METy-genmetlv.Py())/thegennu.Et();
      value["smearedGenMETMinusGenMET_x"]=myRandom.Gaus(0,32);
      value["smearedGenMETMinusGenMET_y"]=myRandom.Gaus(0,32);
      value["BWsmearedGenMETMinusGenMET_x"]=myRandom.BreitWigner(0,32);
      value["BWsmearedGenMETMinusGenMET_y"]=myRandom.BreitWigner(0,32);
      value["MT_RecoLep_genMET_smeared2X"]=getMT(leplv,genmetlv_smeared2X);
      value["MT_RecoLep_genMET_smeared"]=getMT(leplv,genmetlv_smeared);
      value["MT_RecoLep_genMET_BWsmeared2X"]=getMT(leplv,genmetlv_BWsmeared2X);
      value["MT_RecoLep_genMET_BWsmeared"]=getMT(leplv,genmetlv_BWsmeared);
      value["myMT"] = myMT;
      value["myMTn"] = myMTn;
//      cout<<"get gen mt"<<getMT(genleplv,genmetlv)<<endl; 
      bool MET300_noMT = false; bool MET100_MT120 = false; bool noMT = false; bool MT80 = false; bool MT100 = false; bool MT125 = false; bool MT150 = false;
      bool SR = false; bool SRbutMT(false), SRbutMET(false), SRbutMT2W(false), SRbutchi2(false), SRbutminDPhi(false);
      bool Preselection(false);
      if(NGJets>1&&MET>50) Preselection = true;
      if(MT>120) MET100_MT120 = true;
      if(MET>200){
	noMT = true;
	if(MT>80) MT80 = true;
	if(MT>100) MT100 = true;
	if(MT>125) MT125 = true;
	if(MT>150) MT150 = true;
      }
      if(MET>300) MET300_noMT = true;
      if(MT>120 && MET>100 && MT2W>-1 && minDPhi>0.8 && chi2<5) SR = true;
      if(MT>-.1 && MET>100 && MT2W>-1. && minDPhi>0.8 && chi2<5) SRbutMT = true;
      if(MT>120 && MET>-.1 && MT2W>-1. && minDPhi>0.8 && chi2<5) SRbutMET = true;
      if(MT>120 && MET>100 && MT2W>-.1 && minDPhi>0.8 && chi2<5) SRbutMT2W = true;
      if(MT>120 && MET>100 && MT2W>-1. && minDPhi>0.8 && chi2>-1) SRbutchi2 = true;
      if(MT>120 && MET>100 && MT2W>-1. && minDPhi>-99 && chi2<5) SRbutminDPhi = true;

      for(unsigned int i = 0; i<histonames.size(); ++i){
	string d = "_";
	string mname;
        mname="W_selection"+d+histonames[i];
        histos[mname]->Fill(value[histonames[i] ],weight); 
        if (!thegenlepfromtau) {
        mname="W_selection_notau"+d+histonames[i];
        histos[mname]->Fill(value[histonames[i] ],weight); 
        }
        if (lep1ismu&&!thegenlepfromtau){
        mname="Wmunu_selection"+d+histonames[i];
        histos[mname]->Fill(value[histonames[i] ],weight); 
        }
        if (lep1isel&&!thegenlepfromtau){
        mname="Wenu_selection"+d+histonames[i];
        histos[mname]->Fill(value[histonames[i] ],weight); 
        }
        if(Preselection){
        mname= "Preselection"+d+histonames[i];
        histos[mname]->Fill(value[histonames[i] ],weight);
        }/*
	if(MET300_noMT){
	  mname= "MET300_noMT"+d+histonames[i];
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
        if(SR){
	  mname = "SR"+d+histonames[i];
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR&&(MET>200)){
	  mname = "SR200"+d+histonames[i];
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}
	if(SR&&(MET>200)&&(MT2W>200)){
	  mname = "SR200MT2W"+d+histonames[i];
	  histos[mname]->Fill(value[histonames[i] ],weight);
	}*/
      }
      for(unsigned int i = 0; i<histonames.size(); ++i){
	value[histonames[i] ] = -1;//reset values
      }
     }  //end of check if it is one lepton event 

  //fill numbers in signal regions.
      if(firstGoodVtxIdx() != 0) continue; ++c1vtx;
      if(NGLeps!=1) continue; ++c1l;
      if(NSLeps!=1) continue; ++cno2l;
      if(!trackveto) continue; ++cno2track;
      if(!tauveto) continue; ++cnotau;
      if(NGJets<4) continue; ++c2j;
      if(NGBJets<1) continue; ++c1b;
     // if(MT<120) continue; ++cmt;
      if(minDPhi<0.8) continue; ++cmdphi;
      if(chi2>5) continue; ++cchi;
      //if(MET<150) continue;

      if(MET>150&&MT>120){
      if(MET>150) histos_cutflow["NEventsPerSignalRegion_"+samplename] ->Fill(0.5,weight);
      if(MET>200) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(1.5,weight);
      if(MET>250) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(2.5,weight);
      if(MET>300) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(MET>150) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(4.5,weight);
	if(MET>200) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(5.5,weight);
	if(MET>250) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(6.5,weight);
	if(MET>300) histos_cutflow["NEventsPerSignalRegion_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(pfMET2Xworse>150&&getMT(leplv,metlv2Xworse)>120){
      if(pfMET2Xworse>150) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename] ->Fill(0.5,weight);
      if(pfMET2Xworse>200) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(1.5,weight);
      if(pfMET2Xworse>250) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(2.5,weight);
      if(pfMET2Xworse>300) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET2Xworse>150) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(4.5,weight);
	if(pfMET2Xworse>200) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(5.5,weight);
	if(pfMET2Xworse>250) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(6.5,weight);
	if(pfMET2Xworse>300) histos_cutflow["NEventsPerSignalRegion_2Xworse_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(pfMET2Xworse_genMET>150&&getMT(leplv,metlv2Xworse_genMET)>120){
      if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename] ->Fill(0.5,weight);
      if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(1.5,weight);
      if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(2.5,weight);
      if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(4.5,weight);
	if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(5.5,weight);
	if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(6.5,weight);
	if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegion_2Xworse_genMET_"+samplename]  ->Fill(7.5,weight);
       }
      }

      if(MET>150&&MT>150){
      if(MET>150) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename] ->Fill(0.5,weight);
      if(MET>200) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(1.5,weight);
      if(MET>250) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(2.5,weight);
      if(MET>300) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(MET>150) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(4.5,weight);
	if(MET>200) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(5.5,weight);
	if(MET>250) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(6.5,weight);
	if(MET>300) histos_cutflow["NEventsPerSignalRegionMT150_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(MET>150&&MT>175){
      if(MET>150) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename] ->Fill(0.5,weight);
      if(MET>200) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(1.5,weight);
      if(MET>250) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(2.5,weight);
      if(MET>300) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(MET>150) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(4.5,weight);
	if(MET>200) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(5.5,weight);
	if(MET>250) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(6.5,weight);
	if(MET>300) histos_cutflow["NEventsPerSignalRegionMT175_"+samplename]  ->Fill(7.5,weight);
       }
      }

      if(pfMET2Xworse>150&&getMT(leplv,metlv2Xworse)>150){
      if(pfMET2Xworse>150) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename] ->Fill(0.5,weight);
      if(pfMET2Xworse>200) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(1.5,weight);
      if(pfMET2Xworse>250) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(2.5,weight);
      if(pfMET2Xworse>300) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET2Xworse>150) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(4.5,weight);
	if(pfMET2Xworse>200) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(5.5,weight);
	if(pfMET2Xworse>250) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(6.5,weight);
	if(pfMET2Xworse>300) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_"+samplename]  ->Fill(7.5,weight);
       }
      }

      if(pfMET2Xworse_genMET>150&&getMT(leplv,metlv2Xworse_genMET)>150){
      if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename] ->Fill(0.5,weight);
      if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(1.5,weight);
      if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(2.5,weight);
      if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(4.5,weight);
	if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(5.5,weight);
	if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(6.5,weight);
	if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT150_2Xworse_genMET_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(pfMET10Xworse_genMET>150&&getMT(leplv,metlv10Xworse_genMET)>150){
      if(pfMET10Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename] ->Fill(0.5,weight);
      if(pfMET10Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(1.5,weight);
      if(pfMET10Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(2.5,weight);
      if(pfMET10Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET10Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(4.5,weight);
	if(pfMET10Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(5.5,weight);
	if(pfMET10Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(6.5,weight);
	if(pfMET10Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT150_10Xworse_genMET_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(pfMET2Xworse_genMET>150&&getMT(leplv,metlv2Xworse_genMET)>175){
      if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename] ->Fill(0.5,weight);
      if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(1.5,weight);
      if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(2.5,weight);
      if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET2Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(4.5,weight);
	if(pfMET2Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(5.5,weight);
	if(pfMET2Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(6.5,weight);
	if(pfMET2Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT175_2Xworse_genMET_"+samplename]  ->Fill(7.5,weight);
       }
      }
      if(pfMET10Xworse_genMET>150&&getMT(leplv,metlv10Xworse_genMET)>175){
      if(pfMET10Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename] ->Fill(0.5,weight);
      if(pfMET10Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(1.5,weight);
      if(pfMET10Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(2.5,weight);
      if(pfMET10Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(3.5,weight);
      if(MT2W>200){
	if(pfMET10Xworse_genMET>150) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(4.5,weight);
	if(pfMET10Xworse_genMET>200) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(5.5,weight);
	if(pfMET10Xworse_genMET>250) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(6.5,weight);
	if(pfMET10Xworse_genMET>300) histos_cutflow["NEventsPerSignalRegionMT175_10Xworse_genMET_"+samplename]  ->Fill(7.5,weight);
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
    h->second->SetBinContent(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			     h->second->GetBinContent(h->second->GetNbinsX(),h->second->GetNbinsY() )+ h->second->GetBinContent(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1) );
    h->second->SetBinError(h->second->GetNbinsX(), h->second->GetNbinsY(), 
			   sqrt(pow(h->second->GetBinError(h->second->GetNbinsX(),h->second->GetNbinsY() ),2)+
				pow(h->second->GetBinError(h->second->GetNbinsX()+1,h->second->GetNbinsY()+1),2) ) );
    }
//save histograms
  string filename = "rootfiles/MTStudies/bvetoCR/"+skimFilePrefix+".root";
  TFile *f = new TFile(filename.c_str(),"RECREATE");
  f->cd();
  for(map<string,TH1F*>::iterator h=    histos.begin(); h!=    histos.end();++h) h->second->Write();
  for(map<string,TH2F*>::iterator h=    histos2d.begin(); h!=    histos2d.end();++h) h->second->Write();
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
