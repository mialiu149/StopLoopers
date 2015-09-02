#include "dataMCplotMaker.h"
#include "TFile.h"
#include "TH1F.h"
#include <iostream>
#include <vector>
#include <string>
#include "TString.h"

using namespace std;

void makePlotPlot(){

  //  vector<char*> bgnames, signames;
  vector<string> bgstrings, sigstrings;
  const unsigned int datasetsize = 8;//12
  const unsigned int bgsetsize = 4;//8
  const unsigned int sigsetsize = datasetsize-bgsetsize;
  //string datasets[datasetsize]={"TTbar","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //char* dataset[datasetsize]={"TTbar","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //string datasets[datasetsize]={"TTbar1l","TTbar2l","TTbarH","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  //char* dataset[datasetsize]={"TTbar1l","TTbar2l","TTbarH","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  string datasets[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  char* dataset[datasetsize]={"TTbar1l","TTbar2l","SingleT","TTV","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
  const unsigned int poststringsize = 7;
  string poststring[poststringsize] = {"MET300_noMT","MET100_MT120","noMT","MT80","MT100","MT125","MT150"};

  for(unsigned int n=0; n<bgsetsize; ++n) {
    //bgnames.push_back(dataset[n]);
    bgstrings.push_back(dataset[n]);
  }
  for(unsigned int n=bgsetsize; n<datasetsize; ++n) {
    //signames.push_back(dataset[n]);
    sigstrings.push_back(dataset[n]);
  }
  TFile *fbg[bgsetsize];
  for(unsigned int n=0; n<bgsetsize; ++n) {
    string datasetstring = datasets[n];
    if(n==0) datasetstring = "TTbar";
    else if(n==1) continue;
    //else if(n==2) continue;
    TString x = ("rootfiles/CorHistos/Cor"+datasetstring+".root");
    cout << x << endl;
      fbg[n] = TFile::Open(x);
      //  cout << fbg[n].IsZombie();
  }
  TFile *fsig[sigsetsize];
  for(unsigned int n=0; n<sigsetsize; ++n) {
    string datasetstring = datasets[bgsetsize+n];
    TString x = ("rootfiles/CorHistos/Cor"+datasetstring+".root");
    cout << x << endl;
    fsig[n] = TFile::Open(x);
  }
  string outputdir = "rootfiles/CorPlots/Cor/";

  vector<string> histonames;

  histonames.push_back("MT2W");               
  histonames.push_back("MT2_b_b");            
  histonames.push_back("MT2_lb_b");           
  histonames.push_back("MT2_lb_bq");          
  //histonames.push_back("MT2_lb_bq_boost");    
  histonames.push_back("MT2_lb_bqq");         
  histonames.push_back("MT2_l_q");            
  histonames.push_back("MT2_lb_b_mless");        
  histonames.push_back("MT2_lb_bq_mless");       
  //histonames.push_back("MT2_lb_bq_boost_mless"); 
  histonames.push_back("MT2_lb_bqq_mless");      
  //histonames.push_back("MT2_l_q_boost");      
  histonames.push_back("MT2_l_qq_mless");     
  histonames.push_back("Mlb");                
  histonames.push_back("Mlbb");               
  histonames.push_back("MTb");                
  histonames.push_back("MTq");                
  histonames.push_back("MTq_boost");          
  histonames.push_back("MTqq");               
  histonames.push_back("Topness");            
  histonames.push_back("MT");                 
  histonames.push_back("MET");                
  histonames.push_back("HT");                 
  histonames.push_back("METoverHT");          
  histonames.push_back("METoverSqrtHT");      
  histonames.push_back("HTratio");            
  histonames.push_back("dRLepBJet");          
  histonames.push_back("dRLepBMax");          
  histonames.push_back("dRLepBMin");          
  histonames.push_back("dRbb");               
  histonames.push_back("chi2");               
  histonames.push_back("NBJets");             
  histonames.push_back("NJets");              
  histonames.push_back("minDPhi");            
  histonames.push_back("minDPhiJ3");          
  histonames.push_back("minDPhiB");           
  histonames.push_back("NBJetsOverNJets");    
  histonames.push_back("pTlb");               
  histonames.push_back("pTlbb");              
  histonames.push_back("pTl");                
  histonames.push_back("pTleadb");            
  histonames.push_back("pTtrailb");           
  //histonames.push_back("sumak8topmass");      
  histonames.push_back("sumak8prunedmass");   
  histonames.push_back("sumak8trimmedmass");  
  histonames.push_back("sumak8filteredmass"); 

  TH1F* null = new TH1F("","",1,0,1);

    for(unsigned int c = 0; c<poststringsize; ++c){
    string prestring = poststring[c];
    for(unsigned int i = 0; i<histonames.size();++i){
      for(unsigned int j = i+1; j<histonames.size();++j){
	vector<TH2F*> bghist; bghist.clear();
	vector<TH2F*> sighist; sighist.clear();
	string options = "--outputName " + outputdir + prestring + "/" + histonames[i] + "VS"  + histonames[j] + " --xAxisLabel " + histonames[i] + "VS"  + histonames[j] + " --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05";
	// string options = "--outputName " + outputdir + prestring + "/" + histonames[i] + " --xAxisLabel " + prestring + histonames[i] + " --noXaxisUnit --energy 13 --lumi 10 --legendTextSize 0.02 --type 1l --preserveBackgroundOrder --legendUp 0.05 --legendRight -0.2";//topness
	for(unsigned int n=0; n<bgsetsize; ++n){
	  int fileidentifier = n;
	  if(n==1) fileidentifier = 0;//TTbar2l
	  //if(n==2) fileidentifier = 0;//TTbarH
	  fbg[fileidentifier]->cd();
	  string name = prestring + "_" + histonames[i]+ "VS"  + histonames[j]+"_"+bgstrings[n];
	  //cout << fbg[fileidentifier]->GetName() << " " << name << endl;
	  TH1F *h = (TH2F*)fbg[fileidentifier]->Get(name.c_str());
	  bghist.push_back(h);
	}
	for(unsigned int n=0; n<sigsetsize; ++n){
	  fsig[n]->cd();
	  string name = prestring + "_" + histonames[i]+ "VS"  + histonames[j]+"_"+sigstrings[n];
	  TH1F *h = (TH2F*)fsig[n]->Get(name.c_str());
	  //cout << fsig[n]->GetName() << " " << name << endl;
	  h->Scale(25.);
	  sighist.push_back(h);
	}
	//dataMCplotMaker(null,bghist,bgstrings,"sig x25","",options,sighist,sigstrings);//get a mass plotter for 2d correlations
      }
    }
  }

} 
