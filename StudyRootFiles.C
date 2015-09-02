#include "TEfficiency.h"
#include "TEventList.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TString.h"
#include "TObject.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TCut.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TMap.h"
#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <string>
#include <cmath>

using namespace std;

//const int samplenamesize = 10;
//string samplename[samplenamesize]={"TTbar","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};
const int samplenamesize = 12;
string samplename[samplenamesize]={"TTbar1l","TTbar2l","TTbarH","SingleT","TTV","WJets","VV","DYJets","Stop_425_325","Stop_500_325","Stop_650_325","Stop_850_100"};

void StudyRootFiles(){

  map<string, TH1F*> histosMT2; //massive
  map<string, TH1F*> histosMT2m;//massless
  map<string, TH1F*> histos; //massive
  vector<string> histonamesMT2; histonamesMT2.clear();
  vector<string> histonames; histonames.clear();
  map<string,float> bgyield;
  map<string,float> sigyield;
  vector<float> histoMT2cut; histoMT2cut.clear();
  vector<float> histocut; histocut.clear();
  vector<float> histoMT2cutmless; histoMT2cutmless.clear();
  
  histonamesMT2.push_back("MT2W");             histoMT2cut.push_back(250.); histoMT2cutmless.push_back(-1.);
  /*  histonamesMT2.push_back("MT2_b_b");          histoMT2cut.push_back(175.); histoMT2cutmless.push_back(250.);
  histonamesMT2.push_back("MT2_b_b_MW");       histoMT2cut.push_back(250.); histoMT2cutmless.push_back(250.);
  histonamesMT2.push_back("MT2_b_b_lMET");     histoMT2cut.push_back(250.); histoMT2cutmless.push_back(125.);
  histonamesMT2.push_back("MT2_b_b_lMET_MW");  histoMT2cut.push_back(275.); histoMT2cutmless.push_back(275.);
  histonamesMT2.push_back("MT2_lb_b");         histoMT2cut.push_back(275.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_lb_bq");        histoMT2cut.push_back(350.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_lb_bq_boost");  histoMT2cut.push_back(350.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_lb_bqq");       histoMT2cut.push_back(275.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_l_q");          histoMT2cut.push_back(100.); histoMT2cutmless.push_back(100.);
  histonamesMT2.push_back("MT2_l_q_boost");    histoMT2cut.push_back(200.); histoMT2cutmless.push_back(200.);
  histonamesMT2.push_back("MT2_l_qq");         histoMT2cut.push_back(225.); histoMT2cutmless.push_back(125.);
  histonames.push_back("Mlb");                 histocut.push_back(180.);
  histonames.push_back("Mbq");                 histocut.push_back(-1.);
  histonames.push_back("Mbqq");                histocut.push_back(-1.);
  histonames.push_back("Mbb");                 histocut.push_back(-1.);
  histonames.push_back("Mlbb");                histocut.push_back(-1.);
  histonames.push_back("MTb");                 histocut.push_back(180.);
  histonames.push_back("MTq");                 histocut.push_back(180.);
  histonames.push_back("MTq_boost");           histocut.push_back(600.);
  histonames.push_back("MTqq");                histocut.push_back(225.);
*/  histonames.push_back("Topness");             histocut.push_back(9.);
  histonamesMT2.push_back("MT2_b_b");          histoMT2cut.push_back(175.); histoMT2cutmless.push_back(250.);
  histonamesMT2.push_back("MT2_lb_b");         histoMT2cut.push_back(275.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_lb_bq");        histoMT2cut.push_back(350.); histoMT2cutmless.push_back(175.);
  histonamesMT2.push_back("MT2_lb_bqq");       histoMT2cut.push_back(275.); histoMT2cutmless.push_back(175.);

  TFile *fin[samplenamesize];
  for(unsigned int j = 0; j<samplenamesize; ++j){
    //string rootname = "rootfiles/MT2MTMassStudies/Histos_"+samplename[j]+".root";
    //if(j<=2)  rootname = "rootfiles/MT2MTMassStudies/Histos_TTbar.root";
    string rootname = "rootfiles/MT2MTMassStudies/TestX_"+samplename[j]+".root";
    if(j<=2)  rootname = "rootfiles/MT2MTMassStudies/TestX_TTbar.root";
    //string rootname = "rootfiles/MT2MTMassStudies/TestTopness/ExtraPrintNoCMRecoTop_NoCMRecoTop/Histoadd2j_"+samplename[j]+".root";
    //if(j<=2)  rootname = "rootfiles/MT2MTMassStudies/TestTopness/ExtraPrintNoCMRecoTop_NoCMRecoTop/Histoadd2j_TTbar.root";
    //string rootname = "rootfiles/MT2MTMassStudies/Histosadd1b_"+samplename[j]+".root";
    //if(j<=2)  rootname = "rootfiles/MT2MTMassStudies/Histosadd1b_TTbar.root";
    fin[j] = TFile::Open(rootname.c_str());
    fin[j]->cd();
    for(unsigned int i = 0; i<histonamesMT2.size(); ++i){
      string mapname;
      mapname = histonamesMT2[i] + "_"+samplename[j];
      if(histosMT2.count(mapname) == 0 ) histosMT2[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
      mapname = "PreMT_" + histonamesMT2[i] + "_"+samplename[j];
      if(histosMT2.count(mapname) == 0 ) histosMT2[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
      mapname = histonamesMT2[i] + "_mless_"+samplename[j];
      if(histosMT2.count(mapname) == 0 ) histosMT2[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
      mapname = "PreMT_"+ histonamesMT2[i] + "_mless_"+samplename[j];
      if(histosMT2.count(mapname) == 0 ) histosMT2[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
    }
    for(unsigned int i = 0; i<histonames.size(); ++i){
      string mapname;
      fin[j]->cd();
      mapname = histonames[i] + "_"+samplename[j];
      if(histos.count(mapname) == 0 ) histos[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
      mapname = "PreMT_" + histonames[i] + "_"+samplename[j];
      if(histos.count(mapname) == 0 ) histos[mapname] = (TH1F*)fin[j]->Get(mapname.c_str());
    }
  }

    for(unsigned int i = 0; i<histonamesMT2.size(); ++i){
      float bgTop(0),bgST(0),bgVV(0),bgTTV(0);
      //float sig425(0),sig500(0),sig650(0),sig850(0);
      string mapname;
      mapname = histonamesMT2[i];
      bgTop = histosMT2[mapname+"_TTbar1l"]->Integral(histosMT2[mapname+"_TTbar1l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbar2l"]->Integral(histosMT2[mapname+"_TTbar2l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbarH" ]->Integral(histosMT2[mapname+"_TTbarH" ]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histosMT2[mapname+"_SingleT"]->Integral(histosMT2[mapname+"_SingleT"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histosMT2[mapname+"_VV"]->Integral(histosMT2[mapname+"_VV"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histosMT2[mapname+"_TTV"]->Integral(histosMT2[mapname+"_TTV"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV; 
      sigyield[mapname+"_Stop_425_325"] = histosMT2[mapname+"_Stop_425_325"]->Integral(histosMT2[mapname+"_Stop_425_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histosMT2[mapname+"_Stop_500_325"]->Integral(histosMT2[mapname+"_Stop_500_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histosMT2[mapname+"_Stop_650_325"]->Integral(histosMT2[mapname+"_Stop_650_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histosMT2[mapname+"_Stop_850_100"]->Integral(histosMT2[mapname+"_Stop_850_100"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_850_100"]->GetNbinsX() );
      mapname = "PreMT_"+histonamesMT2[i];
      bgTop = histosMT2[mapname+"_TTbar1l"]->Integral(histosMT2[mapname+"_TTbar1l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbar2l"]->Integral(histosMT2[mapname+"_TTbar2l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbarH" ]->Integral(histosMT2[mapname+"_TTbarH" ]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histosMT2[mapname+"_SingleT"]->Integral(histosMT2[mapname+"_SingleT"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histosMT2[mapname+"_VV"]->Integral(histosMT2[mapname+"_VV"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histosMT2[mapname+"_TTV"]->Integral(histosMT2[mapname+"_TTV"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV;
      sigyield[mapname+"_Stop_425_325"] = histosMT2[mapname+"_Stop_425_325"]->Integral(histosMT2[mapname+"_Stop_425_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histosMT2[mapname+"_Stop_500_325"]->Integral(histosMT2[mapname+"_Stop_500_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histosMT2[mapname+"_Stop_650_325"]->Integral(histosMT2[mapname+"_Stop_650_325"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histosMT2[mapname+"_Stop_850_100"]->Integral(histosMT2[mapname+"_Stop_850_100"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_Stop_850_100"]->GetNbinsX() );
      if(histoMT2cutmless[i]<0) continue;
      mapname = histonamesMT2[i] + "_mless";
      bgTop = histosMT2[mapname+"_TTbar1l"]->Integral(histosMT2[mapname+"_TTbar1l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbar2l"]->Integral(histosMT2[mapname+"_TTbar2l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbarH" ]->Integral(histosMT2[mapname+"_TTbarH" ]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histosMT2[mapname+"_SingleT"]->Integral(histosMT2[mapname+"_VV"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histosMT2[mapname+"_VV"]->Integral(histosMT2[mapname+"_SingleT"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histosMT2[mapname+"_TTV"]->Integral(histosMT2[mapname+"_TTV"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV;
      sigyield[mapname+"_Stop_425_325"] = histosMT2[mapname+"_Stop_425_325"]->Integral(histosMT2[mapname+"_Stop_425_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histosMT2[mapname+"_Stop_500_325"]->Integral(histosMT2[mapname+"_Stop_500_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histosMT2[mapname+"_Stop_650_325"]->Integral(histosMT2[mapname+"_Stop_650_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histosMT2[mapname+"_Stop_850_100"]->Integral(histosMT2[mapname+"_Stop_850_100"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_850_100"]->GetNbinsX() );
      mapname = "PreMT_"+histonamesMT2[i] + "_mless";
      bgTop = histosMT2[mapname+"_TTbar1l"]->Integral(histosMT2[mapname+"_TTbar1l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbar2l"]->Integral(histosMT2[mapname+"_TTbar2l"]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histosMT2[mapname+"_TTbarH" ]->Integral(histosMT2[mapname+"_TTbarH" ]->FindBin(histoMT2cut[i]),histosMT2[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histosMT2[mapname+"_SingleT"]->Integral(histosMT2[mapname+"_SingleT"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histosMT2[mapname+"_VV"]->Integral(histosMT2[mapname+"_VV"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histosMT2[mapname+"_TTV"]->Integral(histosMT2[mapname+"_TTV"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV;
      sigyield[mapname+"_Stop_425_325"] = histosMT2[mapname+"_Stop_425_325"]->Integral(histosMT2[mapname+"_Stop_425_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histosMT2[mapname+"_Stop_500_325"]->Integral(histosMT2[mapname+"_Stop_500_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histosMT2[mapname+"_Stop_650_325"]->Integral(histosMT2[mapname+"_Stop_650_325"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histosMT2[mapname+"_Stop_850_100"]->Integral(histosMT2[mapname+"_Stop_850_100"]->FindBin(histoMT2cutmless[i]),histosMT2[mapname+"_Stop_850_100"]->GetNbinsX() );
    }
    for(unsigned int i = 0; i<histonames.size(); ++i){
      if(histocut[i]<0) continue;
      float bgTop(0),bgST(0),bgVV(0),bgTTV(0);
      //float sig425(0),sig500(0),sig650(0),sig850(0);
      string mapname;
      mapname = histonames[i];
      bgTop = histos[mapname+"_TTbar1l"]->Integral(histos[mapname+"_TTbar1l"]->FindBin(histocut[i]),histos[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histos[mapname+"_TTbar2l"]->Integral(histos[mapname+"_TTbar2l"]->FindBin(histocut[i]),histos[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histos[mapname+"_TTbarH" ]->Integral(histos[mapname+"_TTbarH" ]->FindBin(histocut[i]),histos[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histos[mapname+"_SingleT"]->Integral(histos[mapname+"_SingleT"]->FindBin(histocut[i]),histos[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histos[mapname+"_VV"]->Integral(histos[mapname+"_VV"]->FindBin(histocut[i]),histos[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histos[mapname+"_TTV"]->Integral(histos[mapname+"_TTV"]->FindBin(histocut[i]),histos[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV;
      sigyield[mapname+"_Stop_425_325"] = histos[mapname+"_Stop_425_325"]->Integral(histos[mapname+"_Stop_425_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histos[mapname+"_Stop_500_325"]->Integral(histos[mapname+"_Stop_500_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histos[mapname+"_Stop_650_325"]->Integral(histos[mapname+"_Stop_650_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histos[mapname+"_Stop_850_100"]->Integral(histos[mapname+"_Stop_850_100"]->FindBin(histocut[i]),histos[mapname+"_Stop_850_100"]->GetNbinsX() );
      mapname = "PreMT_"+histonames[i];
      bgTop = histos[mapname+"_TTbar1l"]->Integral(histos[mapname+"_TTbar1l"]->FindBin(histocut[i]),histos[mapname+"_TTbar1l"]->GetNbinsX() );
      bgTop+= histos[mapname+"_TTbar2l"]->Integral(histos[mapname+"_TTbar2l"]->FindBin(histocut[i]),histos[mapname+"_TTbar2l"]->GetNbinsX() );
      bgTop+= histos[mapname+"_TTbarH" ]->Integral(histos[mapname+"_TTbarH" ]->FindBin(histocut[i]),histos[mapname+"_TTbarH" ]->GetNbinsX() );
      bgST  = histos[mapname+"_SingleT"]->Integral(histos[mapname+"_SingleT"]->FindBin(histocut[i]),histos[mapname+"_SingleT"]->GetNbinsX() );
      bgVV  = histos[mapname+"_VV"]->Integral(histos[mapname+"_VV"]->FindBin(histocut[i]),histos[mapname+"_VV"]->GetNbinsX() );
      bgTTV = histos[mapname+"_TTV"]->Integral(histos[mapname+"_TTV"]->FindBin(histocut[i]),histos[mapname+"_TTV"]->GetNbinsX() );
      bgyield[mapname] = bgTop+bgST+bgVV+bgTTV;
      sigyield[mapname+"_Stop_425_325"] = histos[mapname+"_Stop_425_325"]->Integral(histos[mapname+"_Stop_425_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_425_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_500_325"] = histos[mapname+"_Stop_500_325"]->Integral(histos[mapname+"_Stop_500_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_500_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_650_325"] = histos[mapname+"_Stop_650_325"]->Integral(histos[mapname+"_Stop_650_325"]->FindBin(histocut[i]),histos[mapname+"_Stop_650_325"]->GetNbinsX() );
      sigyield[mapname+"_Stop_850_100"] = histos[mapname+"_Stop_850_100"]->Integral(histos[mapname+"_Stop_850_100"]->FindBin(histocut[i]),histos[mapname+"_Stop_850_100"]->GetNbinsX() );
    }
    cout << "Event yields" << endl;
    for(unsigned int i = 0; i<histonamesMT2.size(); ++i){
      string mapname = "PreMT_"+histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
      mapname = histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
      if(histoMT2cutmless[i]<0) continue;
      mapname = "PreMT_"+histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
      mapname = histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
    }
    for(unsigned int i = 0; i<histonames.size(); ++i){
      if(histocut[i]<0) continue;
      string mapname = "PreMT_"+histonames[i];
      cout << mapname << ", cut at " << histocut[i] <<": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
      mapname = histonames[i];
      cout << mapname << ", cut at " << histocut[i] <<": bg=" << bgyield[mapname ] << ", sig: 425=" << sigyield[mapname+"_Stop_425_325" ]
	   << ", 500=" << sigyield[mapname+"_Stop_500_325" ] << ", 650=" << sigyield[mapname+"_Stop_650_325" ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ] << endl;
    }

    cout << endl << "now S/B" << endl;
    for(unsigned int i = 0; i<histonamesMT2.size(); ++i){
      string mapname = "PreMT_"+histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
      mapname = histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
      if(histoMT2cutmless[i]<0) continue;
      mapname = "PreMT_"+histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
      mapname = histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
    }
    for(unsigned int i = 0; i<histonames.size(); ++i){
      if(histocut[i]<0) continue;
      string mapname = "PreMT_"+histonames[i];
      cout << mapname << ", cut at " << histocut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
      mapname = histonames[i];
      cout << mapname << ", cut at " << histocut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/bgyield[mapname ] << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/bgyield[mapname ] 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/bgyield[mapname ] << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/bgyield[mapname ] << endl;
    }

    cout << endl << "now S/sqrt(B)" << endl;
    for(unsigned int i = 0; i<histonamesMT2.size(); ++i){
      string mapname = "PreMT_"+histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
      mapname = histonamesMT2[i];
      cout << mapname << ", cut at " << histoMT2cut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
      if(histoMT2cutmless[i]<0) continue;
      mapname = "PreMT_"+histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
      mapname = histonamesMT2[i]+"_mless";
      cout << mapname << ", cut at " << histoMT2cutmless[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
    }
    for(unsigned int i = 0; i<histonames.size(); ++i){
      if(histocut[i]<0) continue;
      string mapname = "PreMT_"+histonames[i];
      cout << mapname << ", cut at " << histocut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
      mapname = histonames[i];
      cout << mapname << ", cut at " << histocut[i] << ": sig: 425=" << sigyield[mapname+"_Stop_425_325" ]/sqrt(bgyield[mapname]) << ", 500=" << sigyield[mapname+"_Stop_500_325" ]/sqrt(bgyield[mapname]) 
           << ", 650=" << sigyield[mapname+"_Stop_650_325" ]/sqrt(bgyield[mapname]) << ", 800=" << sigyield[mapname+"_Stop_850_100" ]/sqrt(bgyield[mapname]) << endl;
    }

}
