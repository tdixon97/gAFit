// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__GA_FIT__H
#define __BAT__GA_FIT__H

#include <BAT/BCModel.h>
#include "TH1D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include <string>
#include <vector>
#include <BAT/BCMath.h>
#include "TH3D.h"
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


// This is a gA_fit header file.
// Model source code is located in file gA_fit/gA_fit.cxx

// ---------------------------------------------------------
class gA_fit : public BCModel
{

public:

    // Constructor

  gA_fit(const std::string& name,TString path_theory,TString path_data,TString path_bkg,TString path_eff,bool float_eff,bool float_rb,bool posS);
  void GetEffTH3(TString PathMarkov);

    // Destructor
    ~gA_fit();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars);
  void ReconstructFit(TH1D *&h,TH1D *&hs,double ga,double s,double b,std::vector<double>effpars);
  void UpdateEff(std::vector<double>effpars);

    // Overload LogAprioriProbability if not using built-in 1D priors
    double LogAPrioriProbability(const std::vector<double> & pars);

    // Overload CalculateObservables if using observables
  void CalculateObservables(const std::vector<double> & pars);
  double GetT12(double ga,double s);
  void SetFloatEff(bool f){fFloatEfficiency=f;};
  void SetInputs(TString in, TString in_bkg)
  {
    fFileIn = new TFile(in);
    TH1D * h = (TH1D*)fFileIn->Get("Energy1");
    TF1 *f = (TF1*)fFileIn->Get("funeff");
    TGraph *g = (TGraph*)fFileIn->Get("Graph");
    
    SetData(h);
    SetEff(f);

    fFileInBkg = new TFile(in_bkg);

    TH1D *hb = (TH1D*)fFileInBkg->Get("MC");
    SetBkg(hb);
  };
    
  void SetData(TH1D *h){fData=h;};
  void SetBkg(TH1D *h){fBkg=h;};
  void SetEff(TF1 *f){fEff=f;};
  
  TH1D* GetData(){return fData;};
  TH1D * GetBkg(){return fBkg;};
  TH3D * fEffTH3;
  
  bool verbose;
  int counter=0;
  int fIndexEffp0;
  bool fFloatEfficiency;
  bool fFloatRb;
private:
  TH1D *fData;
  TH1D *fBkg;
  gATheoryHandler *fTheory;
  TF1 * fEff;
  
  double N,t;
  int first;
  int last;
  TFile *fFileIn;
  TFile *fFileInBkg;
  bool fFloatBias;
};
// ---------------------------------------------------------

#endif
