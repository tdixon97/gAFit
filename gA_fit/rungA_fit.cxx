// ***************************************************************
// This file was created using the bat-project script
// for project gA_fit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include "gA_fit.h"
#include "TString.h"
#include "TFile.h"
#include "TH1D.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
#include <utility>
#include "TSpectrum.h"
#include <fstream>
#include "TFile.h"
#include <BAT/BCLog.h>
#include "TRandom3.h"
#include "TF1.h"
#include "TApplication.h"
#include <sstream>


void MakeOutputPlots(gA_fit &m,std::string name)
{

  TCanvas * can = new TCanvas();
  
  std::vector<double> pars = m.GetBestFitParameters();

  
  TH1D * h =(TH1D*)m.GetData();
  h->GetListOfFunctions()->Clear();
  h->SetLineColor(kBlack);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(6);

  h->SetTitle(" ; Energy [keV] ; counts/keV; ");

  TH1D * hb =(TH1D*)m.GetBkg();

  hb->SetLineColor(8);
  double ga = pars[0];
  double s =pars[1];
  double b =pars[2];

  std::vector<double>effpars;
  if (m.fFloatEfficiency)
    {
      effpars[0]=pars[3];
      effpars[1]=pars[4];
      effpars[2]=pars[5];
    }
      
  TH1D *hmodel = (TH1D*)h->Clone();
  hmodel->SetName("hmodel");
  hmodel->Clear(); hmodel->Reset();

  TH1D *hsignal = (TH1D*)h->Clone();
  hsignal->SetName("hsignal");
  hsignal->Clear(); hsignal->Reset();

  m.ReconstructFit(hmodel,hsignal,ga,s,b,effpars);

  hb->Scale(b);
  hmodel->SetLineColor(2);

  hsignal->SetLineColor(6);

  TLegend * l = new TLegend(0.7,0.7,0.9,0.9);
  gStyle->SetOptStat(0);
  l->AddEntry(h,"Data");
  l->AddEntry(hb,"Bkg");
  l->AddEntry(hmodel,"Reconstruction");
  l->AddEntry(hsignal,"^{113}Cd #beta^{-} decay");

  
  h->Draw("EP");
  hmodel->SetLineWidth(1);
  hb->SetLineWidth(1);
  hsignal->SetLineWidth(1);

  hmodel->Draw("sameHIST");
  hb->Draw("sameHIST");
  hmodel->Draw("sameHIST");
  hsignal->Draw("sameHIST");
  l->Draw();
  can->Draw();
  can->SaveAs(Form("output/output_%s.C",name.data()));

}
void Usage()
{

  std::cout << std::endl << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if( std::getenv("USER") != NULL )
        std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
    else
      std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
    std::cout<<"./rungA_fit"<<std::endl;
    std::cout<<"options "<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;
    std::cout<<"-n (or --name)            [fit name]           (default: test_fit)"<<std::endl;
    std::cout<<"-t (or --theory-path)     [path for templates] (default: /home/tdixon/Cd113Shape/output_shell.root) "<<std::endl;
    std::cout<<"-d (or --data-path)       [path for data]      (default: /home/tdixon/BATGraphHistoFit/inputs/cross/data_CWOnat.root) note: shouuld be after running eff fit"<<std::endl;
    std::cout<<"-b (or --bkg-path)        [path for bkg]       (default: /home/tdixon/Cd113Shape/CROSS_CWOnat/histo_bkg_model.root) "<<std::endl;
    std::cout<<"-f (or --eff-path)        [path for eff prior  (default: ~/BATGraphHistoFit/output/cross/eff_output_mcmc.root) "<<std::endl;
    std::cout<<"-p (or --is-positive-m)   [force s-NME>0 ?]    (default: 1)"<<std::endl;
    std::cout<<"-e (or --is-float-eff)    [float eff?]         (default: 1)"<<std::endl;
    std::cout<<"-r (or --is-float-Rb)     [add Rb to fit?]     (default: 0)"<<std::endl;
    std::cout<<"-h (or --help)  "<<std::endl;
}
int main(int argc,char** argv)
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    std::string name="test_fit";
    TString path_theory="/home/tdixon/Cd113Shape/output_shell.root";
    TString path_data="/home/tdixon/BATGraphHistoFit/inputs/cross/data_CWOnat.root";
    TString path_bkg="/home/tdixon/Cd113Shape/CROSS_CWOnat/histo_bkg_model.root";
    TString path_eff="~/BATGraphHistoFit/output/cross/eff_output_mcmc.root";
    bool isPositiveM=1;
    bool isFloatEff=1;
    bool isFloatRb=0;

    {
      
      static struct option long_options[] = {{"name",        required_argument, nullptr, 'n'},
					     {"theory-path", required_argument, nullptr, 't'},
					     {"data-path",   required_argument, nullptr, 'd'},
					     {"bkg-path",    required_argument, nullptr, 'b'},
					     {"eff-path",    required_argument,nullptr,'f'},
					     {"is-positive-m",required_argument,nullptr, 'p'},
					     {"is-float-eff",required_argument, nullptr, 'e'},
					     {"is-float-Rb",required_argument, nullptr,  'r'},
					     { "help",                       no_argument, nullptr, 'h' },
					     {nullptr, 0, nullptr, 0}
      };
      const char* const short_options = "n:t:d:b:p:e:r:f:h";
      int c;
    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
      switch (c) {
      case 'n': {
        name = optarg;
        break;
      }
      case 't':{
        path_theory=optarg;
        break;
      }
      case 'd':{
	path_data=optarg;
	break;
      }
      case 'b':{
	path_bkg=optarg;
	break;
      }
      case 'p':{
	isPositiveM=atoi(optarg);
	break;
      }
      case 'f':{
	path_eff=optarg;
	break;
      }
      case 'e':{
	isFloatEff=atoi(optarg);
	break;
      }
      case 'r':{
	isFloatRb=atoi(optarg);
	break;
      }

      case 'h':{
        Usage();
        return 0;
        break;
      }

      default: {
	exit(1);
      }
      }
    }
  }
                      


    // Create the fitter
    //---------------------------------------------------------------------------------------
    gA_fit m(name,path_theory,path_data,path_bkg,path_eff,isFloatEff,isFloatRb,isPositiveM);


    // Seed initial pars
    //-------------------------------------------------------------------------------------------
    std::vector<double> initialPos ;
    m.FindMode(m.GetBestFitParameters());
    
    std::vector<double> bestFitMinuit = m.GetBestFitParameters();                                                                                                                                       
    m.SetInitialPositions(bestFitMinuit);


    
    // set precision
    m.SetPrecision(BCEngineMCMC::kHigh);
    m.SetNChains(4);
    m.FindMode(m.GetBestFitParameters());

    std::cout<<"Pars: "<<std::endl;
    std::cout<<"gA   : "<<m.GetBestFitParameters()[0]<<std::endl;
    std::cout<<"sNME : "<<m.GetBestFitParameters()[1]<<std::endl;
    std::cout<<"b    : "<<m.GetBestFitParameters()[2]<<std::endl;
    std::cout<<"eff  : "<<m.GetBestFitParameters()[3]<<" , "<<m.GetBestFitParameters()[4]<<" , "<<m.GetBestFitParameters()[5]<<std::endl;
    std::cout<<"T12  " <<m.GetT12(m.GetBestFitParameters()[0],m.GetBestFitParameters()[1])<<std::endl;
    BCLog::OutSummary("Test model created");
    

    //////////////////////////////
    // perform your analysis here
    
    // Write Markov Chain to a ROOT file as a TTree
     m.WriteMarkovChain("output/"+m.GetSafeName() + "_mcmc.root", "RECREATE");

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // run mode finding; by default using Minuit
    m.FindMode(m.GetBestFitParameters());

    // draw all marginalized distributions into a PDF file
    m.PrintAllMarginalized("output/"+m.GetSafeName() + "_plots.pdf");

    // print summary plots
    m.PrintParameterPlot("output/"+m.GetSafeName() + "_parameters.pdf");
    m.PrintCorrelationPlot("output/"+m.GetSafeName() + "_correlation.pdf");
    m.PrintCorrelationMatrix("output/"+m.GetSafeName() + "_correlationMatrix.pdf");
    m.PrintKnowledgeUpdatePlots("output/"+m.GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    m.PrintSummary();
    MakeOutputPlots(m,argv[1]);
    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
