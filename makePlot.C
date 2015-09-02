#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TString.h"

using namespace std;

void makePlot(){

  vector<char*> bgnames, signames;
  const unsigned int datasetsize = 8;//12
  const unsigned int bgsetsize = 4;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  string datasets[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  char* dataset[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  for(unsigned int n=0; n<bgsetsize; ++n) {
	bgnames.push_back(dataset[n]);
  }
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) signames.push_back(dataset[n]);

  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar";
    else if(n==1) continue;
    //else if(n==2) continue;
    TString x = ("rootfiles/MT2MTMassStudies/TestX_"+datasetstring+".root");
    cout << x << endl;
      fbg[n] = TFile::Open(x);
      //  cout << fbg[n].IsZombie();
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
    string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/MT2MTMassStudies/TestX_"+datasetstring+".root");
    //TString x = ("rootfiles/MT2MTMassStudies/TestTopness/ExtraPrintNoCMRecoTop_NoCMRecoTop/Histoadd2b_"+datasetstring+".root");
    //TString x = ("rootfiles/MT2MTMassStudies/TestTopness/TestMinimizer/ChangeFitParameters/Full_"+datasetstring+".root");
    //TString x = ("rootfiles/MT2MTMassStudies/TestMT2input/Histosadd2b_"+datasetstring+".root");
    cout << x << endl;
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "rootfiles/MT2MTMassStudies/plots/TestX/";
  //string outputdir = "rootfiles/MT2MTMassStudies/TestTopness/plots/ExtraPrintNoCMRecoTop_NoCMRecoTop/add2b/";
  //string outputdir = "rootfiles/MT2MTMassStudies/TestTopness/plots/Print/OnlyCM/";
  //string outputdir = "rootfiles/MT2MTMassStudies/TestTopness/Cor/OneTerm/";
  // string outputdir = "rootfiles/MT2MTMassStudies/TestTopness/TestMinimizer/ChangeFitParameters/Full/";
  //string outputdir = "rootfiles/MT2MTMassStudies/TestMT2input/add2b/";

  vector<string> histnames;
  
    histnames.push_back("MT2W");
   /*histnames.push_back("MT2_b_b_MW");
    histnames.push_back("MT2_b_b_MW_mless");
    histnames.push_back("MT2_b_b");
    histnames.push_back("MT2_b_b_mless");
    histnames.push_back("MT2_b_b_lMET_MW");
    histnames.push_back("MT2_b_b_lMET_MW_mless");
    histnames.push_back("MT2_b_b_lMET");
    histnames.push_back("MT2_b_b_lMET_mless");
    histnames.push_back("MT2_l_q");
    histnames.push_back("MT2_l_q_boost");
    histnames.push_back("MT2_l_q_mless");
    histnames.push_back("MT2_l_q_boost_mless");
    histnames.push_back("MT2_l_qq");
    histnames.push_back("MT2_l_qq_mless");
    histnames.push_back("MT2_lb_b");
    histnames.push_back("MT2_lb_b_mless");
    histnames.push_back("MT2_lb_bq");
    histnames.push_back("MT2_lb_bq_boost");
    histnames.push_back("MT2_lb_bq_mless");
    histnames.push_back("MT2_lb_bq_boost_mless");
    histnames.push_back("MT2_lb_bqq");
    histnames.push_back("MT2_lb_bqq_mless");
    histnames.push_back("MTb");
    histnames.push_back("MTq");
    histnames.push_back("MTq_boost");
    histnames.push_back("MTqq");
    histnames.push_back("Mbb");
    histnames.push_back("Mbq");
    histnames.push_back("Mbqq");
    histnames.push_back("Mlb");
    histnames.push_back("Mlbb");
  */ histnames.push_back("Topness");
    //histnames.push_back("GenTopness");
    //histnames.push_back("NewTopness");
    //histnames.push_back("MT2W");
    //histnames.push_back("MET");
    //histnames.push_back("MT");

    histnames.push_back("MT2_b_b");
    histnames.push_back("MT2_lb_b");
    histnames.push_back("MT2_lb_b_mless");
    histnames.push_back("MT2_lb_bq");
    histnames.push_back("MT2_lb_bq_mless");
    histnames.push_back("MT2_lb_bqq");
    histnames.push_back("MT2_lb_bqq_mless");

  TH1F* null = new TH1F("","",1,0,1);

  for(unsigned int c = 0; c<3; ++c){
    string prestring = "";
    if(c==1) prestring = "PreMT_";
    if(c==2) prestring = "PostAll_";
    for(unsigned int i = 0; i<histnames.size();++i){
      vector<TH1F*> bghist; bghist.clear();
      vector<TH1F*> sighist; sighist.clear();
      string options = "--outputName " + outputdir + prestring + histnames[i] + " --xAxisLabel " + prestring + histnames[i] + " --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05";
      // string options = "--outputName " + outputdir + prestring + histnames[i] + " --xAxisLabel " + prestring + histnames[i] + " --noXaxisUnit --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05 --legendRight -0.2";//topness
      for(unsigned int n=0; n<bgsetsize; ++n){
	int fileidentifier = n;
	if(n==1) fileidentifier = 0;//TTbar2l
	//if(n==2) fileidentifier = 0;//TTbarH
	fbg[fileidentifier]->cd();
	string name = prestring + histnames[i]+"_"+bgnames[n];
	//cout << fbg[fileidentifier]->GetName() << " " << name << endl;
	TH1F *h = (TH1F*)fbg[fileidentifier]->Get(name.c_str());
	bghist.push_back(h);
      }
      for(unsigned int n=0; n<sigsetsize; ++n){
	fsig[n]->cd();
	string name = prestring + histnames[i]+"_"+signames[n];
	TH1F *h = (TH1F*)fsig[n]->Get(name.c_str());
	//cout << fsig[n]->GetName() << " " << name << endl;
	if(c==2) h->Scale(10.);
	else if(c==1) h->Scale(100.);
	else if(c==0) h->Scale(15.);
	sighist.push_back(h);
      }
      if(c==2) dataMCplotMaker(null,bghist,bgnames,"sig x10","",options,sighist,signames);
      else if(c==1) dataMCplotMaker(null,bghist,bgnames,"sig x100","",options,sighist,signames);
      else if(c==0) dataMCplotMaker(null,bghist,bgnames,"sig x15","",options,sighist,signames);
    }
  }

} 
