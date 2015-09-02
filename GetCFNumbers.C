#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TString.h"

using namespace std;

void GetCFNumbers(){

  //  vector<char*> bgnames, signames;
  vector<string> bgstrings, sigstrings;
  const unsigned int datasetsize = 11;//12
  const unsigned int bgsetsize = 7;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;

  string datasets[datasetsize]={"TTbar2l","TTbar1l","TTbarH","SingleT","TTV","WJets","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  string datasetL[datasetsize]={"$\\mathrm{t}\\bar{\\mathrm{t}}\\rightarrow \\ell\\ell$","$1\\ell$ Top","$0\\ell$ Top","single Top","$\\mathrm{t}\\bar{\\mathrm{t}}+\\mathrm{V}$","W + jets","VV","$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (425,325)","$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (500,325)","$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (650,325)","$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (850,100)"};
  char* dataset[datasetsize]={"TTbar2l","TTbar1l","TTbarH","SingleT","TTV","WJets","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  bool skipping[datasetsize]={false,false,false,false,false,false,true,false,false,false,false};
  const unsigned int poststringsize = 7;
  string poststring[poststringsize] = {""};

  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) continue;
    //bgnames.push_back(dataset[n]);
    bgstrings.push_back(dataset[n]);
  }
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    if(skipping[n]) continue;
    //signames.push_back(dataset[n]);
    sigstrings.push_back(dataset[n]);
  }
  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) continue;
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar";
    else if(n==1) continue;
    else if(n==2) continue;
    TString x = ("rootfiles/CutflowExtra"+datasetstring+".root");
    cout << x << endl;
      fbg[n] = TFile::Open(x);
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
     if(skipping[bgsetsize+n]) continue;
   string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/CutflowExtra"+datasetstring+".root");
    cout << x << endl;
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "rootfiles/CutPlots/CutFlow/";

  vector<string> histonames;

  histonames.push_back("CutflowHighDM");               

  TH1D* null = new TH1D("","",1,0,1);
  cout << __LINE__ << endl;
    for(unsigned int i = 0; i<histonames.size();++i){
      vector<TH1D*> bghist; bghist.clear();
      vector<TH1D*> sighist; sighist.clear();
      for(unsigned int n=0; n<bgsetsize; ++n){
	cout << n << " " << __LINE__ << endl;
	if(skipping[n]) continue;
	int fileidentifier = n;
	if(n==1) fileidentifier = 0;//TTbar2l
	if(n==2) fileidentifier = 0;//TTbarH
	fbg[fileidentifier]->cd();
	string name = histonames[i]+"_"+bgstrings[n];
	//cout << fbg[fileidentifier]->GetName() << " " << name << endl;
	TH1D *h = (TH1D*)fbg[fileidentifier]->Get(name.c_str());
	//cout << __LINE__<< " " << h->Integral() << endl;
	bghist.push_back(h);
      }
      for(unsigned int n=0; n<sigsetsize; ++n){
	cout << n << " " << __LINE__ << endl;
	if(skipping[bgsetsize+n]) continue;
	fsig[n]->cd();
	string name = histonames[i]+"_"+sigstrings[n];
	TH1D *h = (TH1D*)fsig[n]->Get(name.c_str());
	h->Scale(1.);
	sighist.push_back(h);
      }
      TH1D *bgsum = (TH1D*)bghist[0]->Clone("bgsum");
      for(unsigned int n=1; n<bgsetsize; ++n){
	cout << n << " " << __LINE__ << endl;
	if(skipping[n]) continue;
	bgsum->Add(bghist[n],1.);
      }
      //HACK!!!
      bghist[4]->Add(bghist[2],1);
      //HACK!!
      //std::ostringstream* fLogStream     = 0;
      std::ostringstream* fLogStream = new std::ostringstream();
      for(unsigned int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	if(n==2) continue;//HACK!!
	for(unsigned int m=1; m<bghist[n]->GetNbinsX();++m){
	  * fLogStream << datasetL[n] << "  " << bghist[n]->GetXaxis()->GetBinLabel(m) << " " << bghist[n]->GetBinContent(m) << endl;
	}
      }
      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	if(n==2) continue;//HACK!!
	for(unsigned int m=1; m<sighist[n]->GetNbinsX();++m){
	  * fLogStream << datasetL[bgsetsize+n] << "  " << sighist[n]->GetXaxis()->GetBinLabel(m) << " " << sighist[n]->GetBinContent(m) << endl;
	}
      }
      cout << fLogStream->str();

    }
} 
