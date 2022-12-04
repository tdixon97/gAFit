// ***************************************************************
// This file was created using the bat-project script
// for project gA_fit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>

#include "gA_fit.h"
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
  
  TH1D *hmodel = (TH1D*)h->Clone();
  hmodel->SetName("hmodel");
  hmodel->Clear(); hmodel->Reset();

  TH1D *hsignal = (TH1D*)h->Clone();
  hsignal->SetName("hsignal");
  hsignal->Clear(); hsignal->Reset();

  m.ReconstructFit(hmodel,hsignal,ga,s,b);

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
int main(int argc,char** argv)
{
    // open log file
    BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);

    if (argc!=4)
      {
	std::cout<<"Error run as ./gA_fit <name> <prefix> <isPositiveM>"<<std::endl;
	return -1;
      }
    // create new gA_fit object
    gA_fit m(argv[1],argv[2],atoi(argv[3]));

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNChains(4);
    m.FindMode(m.GetBestFitParameters());

    std::cout<<"Pars: "<<std::endl;
    std::cout<<"gA   : "<<m.GetBestFitParameters()[0]<<std::endl;
    std::cout<<"sNME : "<<m.GetBestFitParameters()[1]<<std::endl;
    std::cout<<"b    : "<<m.GetBestFitParameters()[2]<<std::endl;
    std::cout<<"T12  " <<m.GetT12(m.GetBestFitParameters()[0],m.GetBestFitParameters()[1])<<std::endl;
    BCLog::OutSummary("Test model created");

    vector<double> initialPos ;
    m.FindMode(initialPos);
    vector<double> bestFitMinuit = m.GetBestFitParameters();                                                                                                                                       
    m.SetInitialPositions(bestFitMinuit);  
    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // m.Normalize();

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
