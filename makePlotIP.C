#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TString.h"

using namespace std;

void makePlotIP(){

  //  vector<char*> bgnames, signames;
  vector<string> bgstrings, sigstrings;
  vector<Color_t> colors;
  const unsigned int datasetsize = 8;//12
  const unsigned int bgsetsize = 4;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  //string datasets[datasetsize]={"TTbar","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //char* dataset[datasetsize]={"TTbar","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //string datasets[datasetsize]={"TTbar1l","TTbar2l","TTbarH","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //char* dataset[datasetsize]={"TTbar1l","TTbar2l","TTbarH","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  string datasets[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  char* dataset[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  const unsigned int poststringsize = 6;
  //string poststring[poststringsize] = {"MET300_noMT","MET100_MT120","noMT","MT80","MT100","MT125","MT150","SR","SR200","SR200MT2W"};
  string poststring[poststringsize] = {"fakeMu","fakeE","promptMu","promptE","promptfromTauMu","promptfromTauE"};
  const Color_t mccolor[datasetsize]={kRed-7,kCyan-3,kOrange-2,kMagenta-5,kGreen+2,kYellow+1,kBlue,kMagenta};

  for(unsigned int n=0; n<bgsetsize; ++n) {
    //bgnames.push_back(dataset[n]);
    bgstrings.push_back(dataset[n]);
    colors.push_back(mccolor[n]);
  }
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    //signames.push_back(dataset[n]);
    sigstrings.push_back(dataset[n]);
    colors.push_back(mccolor[n]);
  }
  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar";
    else if(n==1) continue;
    //else if(n==2) continue;
    TString x = ("rootfiles/TestIP/Histo"+datasetstring+".root");
    cout << x << endl;
      fbg[n] = TFile::Open(x);
      //  cout << fbg[n].IsZombie();
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
    string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/TestIP/Histo"+datasetstring+".root");
    cout << x << endl;
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "rootfiles/TestIP/";

  vector<string> histonames;

  histonames.push_back("d0");     
  histonames.push_back("dz");      
  histonames.push_back("ip3d");    
  histonames.push_back("d0err");   
  histonames.push_back("dzerr");   
  histonames.push_back("ip3derr"); 
  histonames.push_back("SSd0");    

  TH1F* null = new TH1F("","",1,0,1);

    for(unsigned int c = 0; c<poststringsize; ++c){
    string prestring = poststring[c];
    for(unsigned int i = 0; i<histonames.size();++i){
      vector<TH1F*> bghist; bghist.clear();
      vector<TH1F*> sighist; sighist.clear();
      string options = "--outputName " + outputdir + prestring + "_" + histonames[i] + " --xAxisLabel " + histonames[i] + " --noXaxisUnit --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05";
      // string options = "--outputName " + outputdir + prestring + "/" + histonames[i] + " --xAxisLabel " + prestring + histonames[i] + " --noXaxisUnit --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05 --legendRight -0.2";//topness
      for(unsigned int n=0; n<bgsetsize; ++n){
	int fileidentifier = n;
	if(n==1) fileidentifier = 0;//TTbar2l
	//if(n==2) fileidentifier = 0;//TTbarH
	fbg[fileidentifier]->cd();
	string name = prestring + "_" + histonames[i]+"_"+bgstrings[n];
	//cout << fbg[fileidentifier]->GetName() << " " << name << endl;
	TH1F *h = (TH1F*)fbg[fileidentifier]->Get(name.c_str());
	bghist.push_back(h);
      }
      for(unsigned int n=0; n<sigsetsize; ++n){
	fsig[n]->cd();
	string name = prestring + "_" + histonames[i]+"_"+sigstrings[n];
	TH1F *h = (TH1F*)fsig[n]->Get(name.c_str());
	//cout << fsig[n]->GetName() << " " << name << endl;
	h->Scale(25.);
	sighist.push_back(h);
      }
      dataMCplotMaker(null,bghist,bgstrings,"sig x25","",options,sighist,sigstrings,colors);
    }
  }

} 
