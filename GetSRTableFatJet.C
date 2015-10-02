#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "TString.h"
#include "TMath.h"

using namespace std;
float Zbi_(float sig, float bg, float bgrelunc=0.3){
     double bgunc = bgrelunc*bg;
     double tau = bg/pow(bgunc,2);//bgunc is absolute
     double n_on = sig+bg;//total yield in SR = sig + bg
     double n_off = tau*bg;
     double P_Bi = TMath::BetaIncomplete(1./(1.+tau),n_on,n_off+1);
     double Z_Bi = sqrt(2.)*TMath::ErfInverse(1.-2.*P_Bi);
     return Z_Bi;
}

void GetSRTableFatJet(){
  vector<string> bgstrings, sigstrings;
  const unsigned int datasetsize = 11;//12
  const unsigned int bgsetsize = 7;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  string datasets[datasetsize]={"TTbar2l","TTbar1l","TTbarH","tW_fatjet","ttz_fatjet","wjets_fatjet","VV","Stop_425_325_fatjet","Stop_500_325","Stop_650_325","Stop_850_100_fatjet"};
  string datasetL[datasetsize]={"$\\mathrm{t}\\bar{\\mathrm{t}}\\rightarrow \\ell\\ell$",
                                "semileptonic $\\mathrm{t}\\bar{\\mathrm{t}}$",
                                "$0\\ell$ Top","single Top","$\\mathrm{t}\\bar{\\mathrm{t}}+\\mathrm{V}$",
                                "W + jets",
                                "VV",
                                "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (425,325)",
                                "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (500,325)",
                                "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (650,325)",
                                "$\\tilde{\\mathrm{t}}\\rightarrow\\mathrm{t}\\tilde{\\chi}_{1}^{0}$ (850,100)"};
  char* dataset[datasetsize]={"TTbar2l","TTbar1l","TTbarH","tW_fatjet","ttz_fatjet","wjets_fatjet","VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100_fatjet"};
  char* datacard[datasetsize]={"TTbar2l","TTbar1l","TTbarH","SingleTop",      "ttV"       ,"WJets",       "VV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  bool skipping[datasetsize]={false,     false,    true,    false,     false,        false,         true,     true,          true,          true,          false};
  const unsigned int poststringsize = 7;
  string poststring[poststringsize] = {""};

  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) bgstrings.push_back("skip");
    else bgstrings.push_back(dataset[n]);
  }
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    if(skipping[n]) sigstrings.push_back("skip");
    else sigstrings.push_back(dataset[n]);
  }

  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    if(skipping[n]) continue;
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar_fatjet";
    else if(n==1) datasetstring = "TTbar_fatjet";
    else if(n==2) datasetstring = "TTbar_fatjet";
    TString x = ("~mliu/public_html/rootfiles/CutHistos/FatJetPlots/"+datasetstring+".root");
    cout << x << endl;
      fbg[n] = TFile::Open(x);
      //  cout << fbg[n].IsZombie();
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
     if(skipping[bgsetsize+n]) continue;
   string datasetstring = datasets[bgsetsize+n];
    TString x = ("~mliu/public_html/rootfiles/CutHistos/FatJetPlots/"+datasetstring+".root");
    cout << x << endl;
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "~mliu/public_html/rootfiles/CutPlots/CutFlow/";

  vector<string> histonames;

  histonames.push_back("NEventsPerSignalRegion_nominal");               
  histonames.push_back("NEventsPerSignalRegion_nominal_mod");               
  histonames.push_back("NEventsPerSignalRegion_boost1");               
  histonames.push_back("NEventsPerSignalRegion_boost2");               
  histonames.push_back("NEventsPerSignalRegion_boost3");               
  histonames.push_back("NEventsPerSignalRegion_boost4");               
  histonames.push_back("NEventsPerSignalRegion_boost5");               
  histonames.push_back("NEventsPerSignalRegion_boost6");               

  TH1D* null = new TH1D("","",1,0,1);
    for(unsigned int i = 0; i<histonames.size();++i){
      vector<TH1D*> bghist; bghist.clear();
      vector<TH1D*> sighist; sighist.clear();
      for(unsigned int n=0; n<bgsetsize; ++n){
        if(skipping[n]) bghist.push_back(null);
	if(skipping[n]) continue;
	int fileidentifier = n;
	if(n==1) fileidentifier = 0;//TTbar2l
	if(n==2) fileidentifier = 0;//TTbarH
	fbg[fileidentifier]->cd();
	string name = histonames[i]+"_"+bgstrings[n];
	TH1D *h = (TH1D*)fbg[fileidentifier]->Get(name.c_str());
	bghist.push_back(h);
      }

      for(unsigned int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]){sighist.push_back(null);continue;}
	fsig[n]->cd();
	string name = histonames[i]+"_"+sigstrings[n];
	TH1D *h = (TH1D*)fsig[n]->Get(name.c_str());
	h->Scale(1.);
	sighist.push_back(h);
      }

      TH1D *bgsum = (TH1D*)bghist[0]->Clone("bgsum");
      for(unsigned int n=1; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	bgsum->Add(bghist[n],1.);
      }

     // calculate zbi

      TH1D *h_zbi = new TH1D("h_zbi","h_zbi",5,0,5);;
      for(unsigned int n=1;n<6;++n){
       float bkg = bgsum->GetBinContent(n); 
       float sig = sighist[3]->GetBinContent(n);
       float zbi = Zbi_(sig,bkg);
       h_zbi->SetBinContent(n,zbi);
      }
//      cout<<Zbi_(17,36,0.3)<<endl;
 //     cout<<Zbi_(17,10,0.3)<<endl;
//      cout<<Zbi_(10,5,0.3)<<endl;
//      cout<<Zbi_(10,7,0.3)<<endl;
//      cout<<Zbi_(10,10,0.3)<<endl;
  
      //HACK!!!
      //bghist[4]->Add(bghist[2],1);
      //HACK!!
      //cout << __LINE__<< "wtf " << bghist[0]->GetBinContent(1) << " " <<bgsum->GetBinContent(1) << " " << sighist[1]->GetBinContent(1) << endl;

      std::ostringstream* fLogStream = new std::ostringstream();
      fLogStream->precision(3);
      *fLogStream << "***************************"<<histonames[i]<<"***********************************" << endl;
      *fLogStream << "********************* input card for limit setting*******************************"<< endl;
       for(int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	if(n==2) continue;//HACK!!
	*fLogStream << datacard[n]<<"SR"<<" "<<"1 "<<bghist[n]->GetBinContent(1)<<" "<<"2 "<<bghist[n]->GetBinContent(2)<<" "<<"3 "<<bghist[n]->GetBinContent(3)<< " "<<"4 "<<bghist[n]->GetBinContent(4)<< " 5 "<<bghist[n]->GetBinContent(5)<< endl;
        *fLogStream << datacard[n]<<"stat1"<<" 1 "<< bghist[n]->GetBinError(1)/bghist[n]->GetBinContent(1)+1<<endl;
        *fLogStream << datacard[n]<<"stat2"<<" 2 "<< bghist[n]->GetBinError(2)/bghist[n]->GetBinContent(2)+1<<endl;
        *fLogStream << datacard[n]<<"stat3"<<" 3 "<< bghist[n]->GetBinError(3)/bghist[n]->GetBinContent(3)+1<<endl;
        *fLogStream << datacard[n]<<"stat4"<<" 4 "<< bghist[n]->GetBinError(4)/bghist[n]->GetBinContent(4)+1<<endl;
        *fLogStream << datacard[n]<<"stat5"<<" 5 "<< bghist[n]->GetBinError(5)/bghist[n]->GetBinContent(5)+1<<endl;
        *fLogStream << datacard[n]<<"syst"<<" 1 1.3 2 1.3 3 1.3 4 1.3 5 1.3"<<endl;
      }

      for(int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
      	*fLogStream << datacard[bgsetsize+n]<<"SR"<<" "<<"1 "<<sighist[n]->GetBinContent(1)<<" "<<"2 "<<sighist[n]->GetBinContent(2)<<" "<<"3 "<<sighist[n]->GetBinContent(3)<< " "<<"4 "<<sighist[n]->GetBinContent(4)<< " 5 "<<sighist[n]->GetBinContent(5)<< endl;
        *fLogStream << datacard[bgsetsize+n]<<"stat1"<<" 1 "<< sighist[n]->GetBinError(1)/sighist[n]->GetBinContent(1)+1<<endl;
        *fLogStream << datacard[bgsetsize+n]<<"stat2"<<" 2 "<< sighist[n]->GetBinError(2)/sighist[n]->GetBinContent(2)+1<<endl;
        *fLogStream << datacard[bgsetsize+n]<<"stat3"<<" 3 "<< sighist[n]->GetBinError(3)/sighist[n]->GetBinContent(3)+1<<endl;
        *fLogStream << datacard[bgsetsize+n]<<"stat4"<<" 4 "<< sighist[n]->GetBinError(4)/sighist[n]->GetBinContent(4)+1<<endl;
        *fLogStream << datacard[bgsetsize+n]<<"stat5"<<" 5 "<< sighist[n]->GetBinError(5)/sighist[n]->GetBinContent(5)+1<<endl;
      }
      
      *fLogStream << "***************************"<<histonames[i]<<"****************************************************" << endl;
      *fLogStream << "\%BEGINLATEX\%"             << endl;
      *fLogStream << "\\begin{table}"             << endl
		  << "\\begin{center}"            << endl
		  << "\\small"                    << endl;
      *fLogStream << "\\begin{tabular}{lccccc}"    << endl;	     
      *fLogStream << "\\hline"                    << endl;
      *fLogStream << " Sample & $E_\\mathrm{T}^\\mathrm{miss} > 300$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 350$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 400$ GeV & $E_\\mathrm{T}^\\mathrm{miss} > 450$ GeV  & $E_\\mathrm{T}^\\mathrm{miss} > 500$ GeV\\\\" << endl;
      *fLogStream << "\\hline" << endl;
      *fLogStream << "High $\\Delta M$ selection & & & & \\\\" << endl;
      for(int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(1) << " \\pm " << bghist[n]->GetBinError(1)/bghist[n]->GetBinContent(1) 
		    << "$ & $" << bghist[n]->GetBinContent(2) << " \\pm " << bghist[n]->GetBinError(2)/bghist[n]->GetBinContent(2)
		    << "$ & $" << bghist[n]->GetBinContent(3) << " \\pm " << bghist[n]->GetBinError(3)/bghist[n]->GetBinContent(3) 
		    << "$ & $" << bghist[n]->GetBinContent(4) << " \\pm " << bghist[n]->GetBinError(4)/bghist[n]->GetBinContent(4)
		    << "$ & $" << bghist[n]->GetBinContent(5) << " \\pm " << bghist[n]->GetBinError(5)/bghist[n]->GetBinContent(5)
		    << "$ \\\\" << endl;
      }
      *fLogStream << "Total & $"
		  << bgsum->GetBinContent(1) << " \\pm " << bgsum->GetBinError(1) 
		  << "$ & $" << bgsum->GetBinContent(2) << " \\pm " << bgsum->GetBinError(2) 
		  << "$ & $" << bgsum->GetBinContent(3) << " \\pm " << bgsum->GetBinError(3) 
		  << "$ & $" << bgsum->GetBinContent(4) << " \\pm " << bgsum->GetBinError(4) 
		  << "$ & $" << bgsum->GetBinContent(5) << " \\pm " << bgsum->GetBinError(5) 
		  << "$ \\\\" << endl;
      for(int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(1) << " \\pm " << sighist[n]->GetBinError(1)/sighist[n]->GetBinContent(1) 
		    << "$ & $" << sighist[n]->GetBinContent(2) << " \\pm " << sighist[n]->GetBinError(2)/sighist[n]->GetBinContent(2) 
		    << "$ & $" << sighist[n]->GetBinContent(3) << " \\pm " << sighist[n]->GetBinError(3)/sighist[n]->GetBinContent(3) 
		    << "$ & $" << sighist[n]->GetBinContent(4) << " \\pm " << sighist[n]->GetBinError(4)/sighist[n]->GetBinContent(4) 
		    << "$ & $" << sighist[n]->GetBinContent(5) << " \\pm " << sighist[n]->GetBinError(5)/sighist[n]->GetBinContent(5) 
		    << "$ \\\\" << endl;
      }
      *fLogStream << "Zbi & $"
		  << h_zbi->GetBinContent(1) 
		  << "$ & $" << h_zbi->GetBinContent(2) 
		  << "$ & $" << h_zbi->GetBinContent(3) 
		  << "$ & $" << h_zbi->GetBinContent(4) 
		  << "$ & $" << h_zbi->GetBinContent(5) 
		  << "$ \\\\" << endl;
       
/*      *fLogStream << "\\hline" << endl;
      *fLogStream << "High $\\Delta M$ selection & & & & \\\\" << endl;
      for(int n=0; n<bgsetsize; ++n){
	if(skipping[n]) continue;
	if(n==2) continue;//HACK!!
	*fLogStream << datasetL[n] << " & $" 
		    << bghist[n]->GetBinContent(5) << " \\pm " << bghist[n]->GetBinError(5) 
		    << "$ & $" << bghist[n]->GetBinContent(6) << " \\pm " << bghist[n]->GetBinError(6) 
		    << "$ & $" << bghist[n]->GetBinContent(7) << " \\pm " << bghist[n]->GetBinError(7) 
		    << "$ & $" << bghist[n]->GetBinContent(8) << " \\pm " << bghist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
      *fLogStream << "Total & $"
		  << bgsum->GetBinContent(5) << " \\pm " << bgsum->GetBinError(5) 
		  << "$ & $" << bgsum->GetBinContent(6) << " \\pm " << bgsum->GetBinError(6) 
		  << "$ & $" << bgsum->GetBinContent(7) << " \\pm " << bgsum->GetBinError(7) 
		  << "$ & $" << bgsum->GetBinContent(8) << " \\pm " << bgsum->GetBinError(8) 
		  << "$ \\\\" << endl;
      for(int n=0; n<sigsetsize; ++n){
	if(skipping[bgsetsize+n]) continue;
	*fLogStream << datasetL[bgsetsize+n] << " & $" 
		    << sighist[n]->GetBinContent(5) << " \\pm " << sighist[n]->GetBinError(5) 
		    << "$ & $" << sighist[n]->GetBinContent(6) << " \\pm " << sighist[n]->GetBinError(6) 
		    << "$ & $" << sighist[n]->GetBinContent(7) << " \\pm " << sighist[n]->GetBinError(7) 
		    << "$ & $" << sighist[n]->GetBinContent(8) << " \\pm " << sighist[n]->GetBinError(8) 
		    << "$ \\\\" << endl;
      }
*/
      *fLogStream << "\\hline\\hline" << endl
		  << "\\end{tabular}" << endl
		  << "\\end{center}"  << endl
		  << "\\end{table}"   << endl
		  << "\%ENDLATEX\%"   << endl
		  << endl;
      cout << fLogStream->str();
    }
} 
