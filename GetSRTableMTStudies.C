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

void GetSRTableMTStudies(){
  vector<string> bgstrings, sigstrings;
  const unsigned int datasetsize = 12;//12
  const unsigned int bgsetsize = 8;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  //string datasets[datasetsize]={"TTbar2l","TTbar1l","TTbarH","SingleT","TTV","WJets","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};string datasets[datasetsize]={"TTbar2l","TTbar1l","ttw","ttz","sch_tch","tW","WJets","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  string datasets[datasetsize]={"TTbar2l","TTbar1l","WJetsHeavy","WJetsLight","ttw","ttz","sch_tch","tW","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  string datasetL[datasetsize]={"$\\mathrm{t}\\bar{\\mathrm{t}}\\rightarrow \\ell\\ell$",
                                "semileptonic $\\mathrm{t}\\bar{\\mathrm{t}}$",
                                "W+b",
                                "W+light jets", 
                                "$\\mathrm{t}\\bar{\\mathrm{t}}+\\mathrm{W}$",
                                "$\\mathrm{t}\\bar{\\mathrm{t}}+\\mathrm{Z}$",
                                "Single top (sch+tch)",
                                "Single top (tW)",
                                 "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (425,325)",
                                 "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (500,325)",
                                 "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (650,325)",
                                 "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (850,100)"};
 // datasets to read and get numbers from 
  string ttbar = "ttbar_powheg_25ns";
  string wjets = "WJets";

  //char* dataset[datasetsize]={"TTbar2l","TTbar1l","TTbarH","SingleT","TTV","WJets","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  char* dataset[datasetsize]={"TTbar2l","TTbar1l","WJetsHeavy","WJetsLight","ttw","ttz","sch_tch","tW","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  bool skipping[datasetsize]={false,false,false,false,false,false,false,false,
                              true,true,true,true};
  const unsigned int poststringsize = 7;
  string poststring[poststringsize] = {""};
//backgrounds
  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) bgstrings.push_back("skip");
    else bgstrings.push_back(dataset[n]);
  }
//signals
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    if(skipping[n]) sigstrings.push_back("skip");
    else sigstrings.push_back(dataset[n]);
  }
//open background files
  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) continue;
    string datasetstring = datasets[n];
    if(n==0) datasetstring = ttbar;
    else if(n==1) datasetstring = ttbar;
    else if(n==2) datasetstring = wjets;
    else if(n==3) datasetstring = wjets;
    TString x = ("rootfiles/MTStudies/MT150/"+datasetstring+".root");
      fbg[n] = TFile::Open(x);
  }
//open signal files
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
    if(skipping[bgsetsize+n]) continue;
    string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/MTStudies/MT150/"+datasetstring+".root");
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "rootfiles/CutPlots/CutFlow/";

  vector<string> histonames;

  histonames.push_back("NEventsPerSignalRegion");               
  histonames.push_back("NEventsPerSignalRegion_2Xworse_genMET");               

  TH1D* null = new TH1D("","",1,0,1);
  vector<TH1D*> bghist; bghist.clear();
  vector<TH1D*> sighist; sighist.clear();

  for(unsigned int n=0; n<bgsetsize; ++n){
    int fileidentifier = n;
    if(skipping[n]) bghist.push_back(null);
    else {
    fbg[fileidentifier]->cd();
    string name_nominal =  histonames[0]+"_"+bgstrings[n];
    string name_2X =  histonames[1]+"_"+bgstrings[n];
    TH1D *h1 = (TH1D*)fbg[fileidentifier]->Get(name_nominal.c_str());
    TH1D *h2 = (TH1D*)fbg[fileidentifier]->Get(name_2X.c_str());
    //h2->Divide(h1);
    bghist.push_back(h2);
    }
  }
      //std::ostringstream* fLogStream     = 0;
      std::ostringstream* fLogStream = new std::ostringstream();
      cout<<"wtf2"<<endl;
      fLogStream->precision(3);
      cout<<"wtf3"<<endl;
      *fLogStream << "*********************************************************************" << endl;
      cout<<"wtf4"<<endl;
      *fLogStream << "\%BEGINLATEX\%"             << endl;
      cout<<"wtf5"<<endl;
      *fLogStream << "\\begin{table}"             << endl
		  << "\\begin{center}"            << endl
		  << "\\small"                    << endl;
      *fLogStream << "\\begin{tabular}{lcccc}"    << endl;	     
      *fLogStream << "\\hline"                    << endl;
      *fLogStream << " Sample & $E_\\mathrm{T}^\\mathrm{miss} > 150$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 200$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 250$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 300$ GeV \\\\" << endl;
      *fLogStream << "\\hline" << endl;
      *fLogStream << "Low $\\Delta M$ selection & & & & \\\\" << endl;
      cout<<"wtf6"<<endl;
      for(unsigned int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
//	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(1) << " \\pm " << bghist[n]->GetBinError(1) 
		    << "$ & $" << bghist[n]->GetBinContent(2) << " \\pm " << bghist[n]->GetBinError(2) 
		    << "$ & $" << bghist[n]->GetBinContent(3) << " \\pm " << bghist[n]->GetBinError(3) 
		    << "$ & $" << bghist[n]->GetBinContent(4) << " \\pm " << bghist[n]->GetBinError(4) 
		    << "$ \\\\" << endl;
      }
/*      *fLogStream << "Total & $"
		  << bgsum->GetBinContent(1) << " \\pm " << bgsum->GetBinError(1) 
		  << "$ & $" << bgsum->GetBinContent(2) << " \\pm " << bgsum->GetBinError(2) 
		  << "$ & $" << bgsum->GetBinContent(3) << " \\pm " << bgsum->GetBinError(3) 
		  << "$ & $" << bgsum->GetBinContent(4) << " \\pm " << bgsum->GetBinError(4) 
		  << "$ \\\\" << endl;
*/      cout<<"wtf8"<<endl;
      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(1) << " \\pm " << sighist[n]->GetBinError(1) 
		    << "$ & $" << sighist[n]->GetBinContent(2) << " \\pm " << sighist[n]->GetBinError(2) 
		    << "$ & $" << sighist[n]->GetBinContent(3) << " \\pm " << sighist[n]->GetBinError(3) 
		    << "$ & $" << sighist[n]->GetBinContent(4) << " \\pm " << sighist[n]->GetBinError(4) 
		    << "$ \\\\" << endl;

      }
      *fLogStream << "\\hline" << endl;
      *fLogStream << "High $\\Delta M$ selection & & & & \\\\" << endl;
      for(unsigned int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
//	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(5) << " \\pm " << bghist[n]->GetBinError(5) 
		    << "$ & $" << bghist[n]->GetBinContent(6) << " \\pm " << bghist[n]->GetBinError(6) 
		    << "$ & $" << bghist[n]->GetBinContent(7) << " \\pm " << bghist[n]->GetBinError(7) 
		    << "$ & $" << bghist[n]->GetBinContent(8) << " \\pm " << bghist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
 /*     *fLogStream << "Total & $"
		  << bgsum->GetBinContent(5) << " \\pm " << bgsum->GetBinError(5) 
		  << "$ & $" << bgsum->GetBinContent(6) << " \\pm " << bgsum->GetBinError(6) 
		  << "$ & $" << bgsum->GetBinContent(7) << " \\pm " << bgsum->GetBinError(7) 
		  << "$ & $" << bgsum->GetBinContent(8) << " \\pm " << bgsum->GetBinError(8) 
		  << "$ \\\\" << endl;
*/      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(5) << " \\pm " << sighist[n]->GetBinError(5) 
		    << "$ & $" << sighist[n]->GetBinContent(6) << " \\pm " << sighist[n]->GetBinError(6) 
		    << "$ & $" << sighist[n]->GetBinContent(7) << " \\pm " << sighist[n]->GetBinError(7) 
		    << "$ & $" << sighist[n]->GetBinContent(8) << " \\pm " << sighist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
      *fLogStream << "\\hline\\hline" << endl
		  << "\\end{tabular}" << endl
		  << "\\end{center}"  << endl
		  << "\\end{table}"   << endl
		  << "\%ENDLATEX\%"   << endl
		  << endl;
      cout << fLogStream->str();
     
//print out
//
  for(unsigned int i = 0; i<histonames.size();++i){
      vector<TH1D*> bghist; bghist.clear();
      vector<TH1D*> sighist; sighist.clear();
      for(unsigned int n=0; n<bgsetsize; ++n){
        int fileidentifier = n;
	if(skipping[n]) bghist.push_back(null);
	else {
	fbg[fileidentifier]->cd();
	string name = histonames[i]+"_"+bgstrings[n];        
	TH1D *h = (TH1D*)fbg[fileidentifier]->Get(name.c_str());
	bghist.push_back(h);
         }
       }
      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	fsig[n]->cd();
	string name = histonames[i]+"_"+sigstrings[n];
	TH1D *h = (TH1D*)fsig[n]->Get(name.c_str());
	h->Scale(1.);
	sighist.push_back(h);
       }
//sim of backgrounds
      TH1D *bgsum = (TH1D*)bghist[0]->Clone("bgsum");
      for(unsigned int n=1; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	bgsum->Add(bghist[n],1.);
      }
      //HACK!!!
  //    bghist[4]->Add(bghist[2],1);
      //HACK!!
      //std::ostringstream* fLogStream     = 0;
      std::ostringstream* fLogStream = new std::ostringstream();
      cout<<"wtf2"<<endl;
      fLogStream->precision(3);
      cout<<"wtf3"<<endl;
      *fLogStream << "*********************************************************************" << endl;
      cout<<"wtf4"<<endl;
      *fLogStream << "\%BEGINLATEX\%"             << endl;
      cout<<"wtf5"<<endl;
      *fLogStream << "\\begin{table}"             << endl
		  << "\\begin{center}"            << endl
		  << "\\small"                    << endl;
      *fLogStream << "\\begin{tabular}{lcccc}"    << endl;	     
      *fLogStream << "\\hline"                    << endl;
      *fLogStream << " Sample & $E_\\mathrm{T}^\\mathrm{miss} > 150$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 200$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 250$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 300$ GeV \\\\" << endl;
      *fLogStream << "\\hline" << endl;
      *fLogStream << "Low $\\Delta M$ selection & & & & \\\\" << endl;
      cout<<"wtf6"<<endl;
      for(unsigned int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
//	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(1) << " \\pm " << bghist[n]->GetBinError(1) 
		    << "$ & $" << bghist[n]->GetBinContent(2) << " \\pm " << bghist[n]->GetBinError(2) 
		    << "$ & $" << bghist[n]->GetBinContent(3) << " \\pm " << bghist[n]->GetBinError(3) 
		    << "$ & $" << bghist[n]->GetBinContent(4) << " \\pm " << bghist[n]->GetBinError(4) 
		    << "$ \\\\" << endl;
      }
/*      *fLogStream << "Total & $"
		  << bgsum->GetBinContent(1) << " \\pm " << bgsum->GetBinError(1) 
		  << "$ & $" << bgsum->GetBinContent(2) << " \\pm " << bgsum->GetBinError(2) 
		  << "$ & $" << bgsum->GetBinContent(3) << " \\pm " << bgsum->GetBinError(3) 
		  << "$ & $" << bgsum->GetBinContent(4) << " \\pm " << bgsum->GetBinError(4) 
		  << "$ \\\\" << endl;
*/      cout<<"wtf8"<<endl;
      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(1) << " \\pm " << sighist[n]->GetBinError(1) 
		    << "$ & $" << sighist[n]->GetBinContent(2) << " \\pm " << sighist[n]->GetBinError(2) 
		    << "$ & $" << sighist[n]->GetBinContent(3) << " \\pm " << sighist[n]->GetBinError(3) 
		    << "$ & $" << sighist[n]->GetBinContent(4) << " \\pm " << sighist[n]->GetBinError(4) 
		    << "$ \\\\" << endl;

      }
      *fLogStream << "\\hline" << endl;
      *fLogStream << "High $\\Delta M$ selection & & & & \\\\" << endl;
      for(unsigned int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
//	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(5) << " \\pm " << bghist[n]->GetBinError(5) 
		    << "$ & $" << bghist[n]->GetBinContent(6) << " \\pm " << bghist[n]->GetBinError(6) 
		    << "$ & $" << bghist[n]->GetBinContent(7) << " \\pm " << bghist[n]->GetBinError(7) 
		    << "$ & $" << bghist[n]->GetBinContent(8) << " \\pm " << bghist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
 /*     *fLogStream << "Total & $"
		  << bgsum->GetBinContent(5) << " \\pm " << bgsum->GetBinError(5) 
		  << "$ & $" << bgsum->GetBinContent(6) << " \\pm " << bgsum->GetBinError(6) 
		  << "$ & $" << bgsum->GetBinContent(7) << " \\pm " << bgsum->GetBinError(7) 
		  << "$ & $" << bgsum->GetBinContent(8) << " \\pm " << bgsum->GetBinError(8) 
		  << "$ \\\\" << endl;
*/      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(5) << " \\pm " << sighist[n]->GetBinError(5) 
		    << "$ & $" << sighist[n]->GetBinContent(6) << " \\pm " << sighist[n]->GetBinError(6) 
		    << "$ & $" << sighist[n]->GetBinContent(7) << " \\pm " << sighist[n]->GetBinError(7) 
		    << "$ & $" << sighist[n]->GetBinContent(8) << " \\pm " << sighist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
      *fLogStream << "\\hline\\hline" << endl
		  << "\\end{tabular}" << endl
		  << "\\end{center}"  << endl
		  << "\\end{table}"   << endl
		  << "\%ENDLATEX\%"   << endl
		  << endl;
      cout << fLogStream->str();
    }
} 
