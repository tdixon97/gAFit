
#ifndef THEORY_H
#define THEORY_H
#include "TF1.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TFile.h"

#include "TH2D.h"
#include "TStyle.h"
#include <chrono>
#include <vector>
#include "TTree.h"
#include <iostream>
#include "TGraph.h"
#include "TString.h"
using namespace std;








class gATheoryHandler
{

public:
  gATheoryHandler(std::string path){SetInput(path); };
  ~gATheoryHandler(){};

  void SetInput(TString path){TFile *f = new TFile(path); fTree = (TTree*)f->Get("spectrum_conv");};
  void bilinear_interpolate(double &f,double x,double y,double x1,double x2,double y1,double y2,double fQ11,double fQ12,double fQ21,double fQ22);
  void GetLowAndHigh(double val,std::vector<double> vals,double &low,double &high,bool verbose=true);
  void CreateArbitaryGraph(double ga,double sNME,int verbose);
  void Normalise();

  void Prepare(TString name);
  TH1D *fHist;
private:
  TTree *fTree;
  std::map<std::tuple<double,double,double>,double> fMapofvalues;     
  double low_ga,high_ga,low_s,high_s;
  vector<double> *s_s;
  vector<double>*E_s;
  vector<double>*ga_s;
};



#endif //THEORY_H
