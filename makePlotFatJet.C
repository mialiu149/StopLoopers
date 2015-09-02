#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TString.h"

using namespace std;

void makePlotFatJet(){

  //  vector<char*> bgnames, signames;
  vector<string> bgstrings, sigstrings;
  vector<Color_t> colors;
  const unsigned int datasetsize = 8;//12
  const unsigned int bgsetsize = 7;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  string datasets[datasetsize]={"TTbar1l","TTbar2l","ttw_fatjet","ttz_fatjet","wjets_fatjet","sch_tch_fatjet","tW_fatjet","Stop_850_100_fatjet"};
  char* dataset[datasetsize]={"TTbar1l","TTbar2l","ttw_fatjet","ttz_fatjet","wjets_fatjet","sch_tch_fatjet","tW_fatjet","Stop_850_100_fatjet"};
  //string datasets[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //char* dataset[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
//  const unsigned int poststringsize = 6;
 // string poststring[poststringsize] = {"SR","SR250","SR250MT2W","SR_boost_200","SR_boost_250","SR_boost_300"};
  //const Color_t mccolor[datasetsize]={kRed-7,kCyan-3,kOrange-2,kMagenta-5,kGreen+2,kYellow+1,kBlue,kMagenta};
  const unsigned int poststringsize =10 ;
  string poststring[poststringsize] = {"SR","SR250","SR250MT2W","SR_boost_200","SR_boost_250",
                                       "SR_boost_300","SR_boost_300_mass75","SR_boost_250_mass75","SR_boost_375_mass75","SR_boost_375_mass75_ak10only"};
  const Color_t mccolor[datasetsize]={kRed-7,kCyan-3,kOrange-2,kMagenta-5,kGreen+2,kBlue,kYellow+1,kMagenta};
  //background color
  for(unsigned int n=0; n<bgsetsize; ++n) {
    bgstrings.push_back(dataset[n]);
    colors.push_back(mccolor[n]);
  }
 //signal color
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    sigstrings.push_back(dataset[n]);
    colors.push_back(mccolor[n]);
  }

  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar_fatjet";
    else if(n==1) continue;
    //else if(n==2) continue;
    TString x = ("rootfiles/CutHistos/FatJetPlots/"+datasetstring+".root");
      fbg[n] = TFile::Open(x);
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
    string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/CutHistos/FatJetPlots/"+datasetstring+".root");
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "~mliu/public_html/FatJetPlots/";

  vector<string> histonames;

  histonames.push_back("MT2W");               
  histonames.push_back("Topness");            
  histonames.push_back("MT");                 
  histonames.push_back("dRLepFatJet");                 
  histonames.push_back("dR_ak8_ak10");                 
  histonames.push_back("pT_diff_jet1_FatJet");                 
  histonames.push_back("dRJet1FatJet");                 
  histonames.push_back("dRJet2FatJet");                 
  histonames.push_back("dRJet3FatJet");                 
  histonames.push_back("dRJet4FatJet");                 
  histonames.push_back("MET");                
  histonames.push_back("HT");                 
  histonames.push_back("METoverHT");          
  histonames.push_back("METoverSqrtHT");      
  histonames.push_back("METoverSqrtHTplusleppt");      
  histonames.push_back("HTratio");            
  histonames.push_back("chi2");               
  histonames.push_back("NBJets");             
  histonames.push_back("NJets");              
  histonames.push_back("NJets_pt40");              
  histonames.push_back("NJets_pt50");              
  histonames.push_back("NJets_pt60");              
  histonames.push_back("NJets_pt70");              
  histonames.push_back("NJets_pt80");              
  histonames.push_back("NJets_pt90");              
  histonames.push_back("NJets_pt100");              
  histonames.push_back("NJets_pt125");              
  histonames.push_back("NJets_pt150");              
  histonames.push_back("NJets_pt175");              
  histonames.push_back("NJets_pt200");              
  histonames.push_back("minDPhi");            
  histonames.push_back("minDPhiJ3");          
  histonames.push_back("minDPhiB");           
  histonames.push_back("pTl4j");               
  histonames.push_back("pTl4joverMET");               
  histonames.push_back("DeltaPhiWl"); 
  histonames.push_back("ak12prunedmass");
  histonames.push_back("ak12trimmedmass");
  histonames.push_back("ak12filteredmass");
  histonames.push_back("ak4_leadpt");            
  histonames.push_back("ak4_subleadpt");            
  histonames.push_back("ak4_leadmass");            
  histonames.push_back("ak4_subleadmass");            
  histonames.push_back("ak10_leadpt");            
  histonames.push_back("ak10_subleadpt");            
  histonames.push_back("ak10mass");
  histonames.push_back("ak10prunedmass");
  histonames.push_back("ak10trimmedmass");
  histonames.push_back("ak10filteredmass");
  histonames.push_back("ak10mass_Corrected");
  histonames.push_back("ak10genjets_mass");            
  histonames.push_back("ak10mass_sublead");
  histonames.push_back("ak10prunedmass_sublead");
  histonames.push_back("ak10trimmedmass_sublead");
  histonames.push_back("ak10filteredmass_sublead");
  histonames.push_back("ak8_leadpt");            
  histonames.push_back("ak8mass");
  histonames.push_back("ak8_Over_ak10");
  histonames.push_back("ak8prunedmass");
  histonames.push_back("ak8trimmedmass");
  histonames.push_back("ak8filteredmass");
  histonames.push_back("DeltaPhiSubFatJet_MET");

  TH1F* null = new TH1F("","",1,0,1);

    for(unsigned int c = 0; c<poststringsize; ++c){
    string prestring = poststring[c];
    for(unsigned int i = 0; i<histonames.size();++i){
      vector<TH1F*> bghist; bghist.clear();
      vector<TH1F*> sighist; sighist.clear();
      string options = "--outputName " + outputdir + prestring + "/" + histonames[i] + " --xAxisLabel " + histonames[i] + " --energy 13 --lumi 5 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05";
      for(unsigned int n=0; n<bgsetsize; ++n){
	int fileidentifier = n;
	if(n==1) fileidentifier = 0;//TTbar2l
	fbg[fileidentifier]->cd();
	string name = prestring + "_" + histonames[i]+"_"+bgstrings[n];
	TH1F *h = (TH1F*)fbg[fileidentifier]->Get(name.c_str());
	bghist.push_back(h);
      }
      for(unsigned int n=0; n<sigsetsize; ++n){
	fsig[n]->cd();
	string name = prestring + "_" + histonames[i]+"_"+sigstrings[n];
	TH1F *h = (TH1F*)fsig[n]->Get(name.c_str());
	h->Scale(10.);
	sighist.push_back(h);
      }
      dataMCplotMaker(null,bghist,bgstrings,"sig x10","",options,sighist,sigstrings,colors);
    }
  }
} 
