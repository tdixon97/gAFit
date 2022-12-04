// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "gA_fit.h"

// #include <BAT/BCMath.h>



void gATheoryHandler::bilinear_interpolate(double &f,double x,double y,double x1,double x2,double y1,double y2,double fQ11,double fQ12,double fQ21,double fQ22)
{
  double prefactor =1;
  prefactor *=1/(x2-x1);
  prefactor *=1/(y2-y1);

  double term1=fQ11*(x2-x)*(y2-y);
  double term2=fQ21*(x-x1)*(y2-y);
  double term3=fQ12*(x2-x)*(y-y1);
  double term4=fQ22*(x-x1)*(y-y1);

  f=prefactor*(term1+term2+term3+term4);


}

void gATheoryHandler::GetLowAndHigh(double val,std::vector<double> vals,double &low,double &high,bool verbose)
{
  double dist_low=1e100;
  double dist_high=1e100;

  for (auto & value: vals)
    {
      if (value<val)
        {
          if(fabs(value-val)<dist_low)
            {
              dist_low=fabs(value-val);
              low=value;
            }
        }
      else
        {
          if(fabs(value-val)<dist_high)
            {
              dist_high=fabs(value-val);
              high=value;
            }
        }



    }

  if (verbose)
    {
      std::cout<<"For the list: "<<std::endl;
      for (auto &value:vals)
        {
          std::cout<<value<<" , ";
        }
      std::cout<<" " <<std::endl;

      std::cout<<"low, val, high  = "<<low<<" , "<<val<<" , "<<high<<std::endl;
    }




}

void gATheoryHandler::Prepare(TString name)
{

  gErrorIgnoreLevel=kWarning;
  double one_ga=0.5;
  double one_s=0;
  double one_E=0.5;

  int N=fTree->Draw("ga",Form("energy==%f &&sNME==%f",one_E,one_s),"goff");

  ga_s = new vector<double>(fTree->GetV1(),fTree->GetV1()+N);
  N=fTree->Draw("sNME",Form("energy==%f &&ga==%f",one_E,one_ga),"goff");
  s_s =new vector<double>(fTree->GetV1(),fTree->GetV1()+N);

  N=fTree->Draw("energy",Form("ga==%f &&sNME==%f",one_ga,one_s),"goff");
  E_s=new vector<double>(fTree->GetV1(),fTree->GetV1()+N);

  fTree->SetEstimate(fTree->GetEntries());
  N = fTree->Draw("ga","","goff");
  std::vector<double>ga_all(fTree->GetV1(),fTree->GetV1()+N);

  N = fTree->Draw("energy","","goff");
  std::vector<double>energy_all(fTree->GetV1(),fTree->GetV1()+N);

  N = fTree->Draw("sNME","","goff");
  std::vector<double>s_all(fTree->GetV1(),fTree->GetV1()+N);

  N = fTree->Draw("spectrum","","goff");
  std::vector<double>val_all(fTree->GetV1(),fTree->GetV1()+N);

  for (int i=0;i<N;i++)
    {
      double gA_tmp = ga_all[i];
      double s_tmp = s_all[i];
      double energy_tmp=energy_all[i];
      double val_tmp=val_all[i];

      fMapofvalues[std::make_tuple(gA_tmp,s_tmp,energy_tmp)]=val_tmp;
    }
  double e0=energy_all[0];
  double e1=energy_all[1];
  double bin=e1-e0;
  std::cout<<bin<<std::endl;
  std::cout<<"e0,e1 = "<<e0<<" , "<<e1<<std::endl;
  double low=e0-bin/2.;
  std::cout<<low<<std::endl;

  fHist= new TH1D(name,name,1000,low,1000*bin);

}
void gATheoryHandler::CreateArbitaryGraph(double ga,double sNME,int verbose)
{



  auto start = chrono::steady_clock::now();
  GetLowAndHigh(ga,*ga_s,low_ga,high_ga,verbose);
  GetLowAndHigh(sNME,*s_s,low_s,high_s,verbose);

  for (auto  & E: *E_s)
    {

      // now we have an energy, we have gA and sNME and we need the 4 relevant fQij                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 \
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     

      double fQ11=fMapofvalues[std::make_tuple(low_ga,low_s,E)];
      double fQ21=fMapofvalues[std::make_tuple(high_ga,low_s,E)];
      double fQ12=fMapofvalues[std::make_tuple(low_ga,high_s,E)];
      double fQ22=fMapofvalues[std::make_tuple(high_ga,high_s,E)];


      double val;
      bilinear_interpolate(val,ga,sNME,low_ga,high_ga,low_s,high_s,fQ11,fQ12,fQ21,fQ22);
      fHist->SetBinContent(fHist->FindBin(E),val);

    }
  auto end = chrono::steady_clock::now();
  if (verbose)
    {
  cout << "Elapsed time: "
       << chrono::duration_cast<chrono::microseconds>(end - start).count()
       << " us" << endl;
    }

  Normalise();





}

void gATheoryHandler::Normalise()
{
                                                                                                 


  double eta = 0.1222;
  double Na= 6.022e26;
  double W=360.24;
  double M = 0.433;


  const double N = eta*Na*M/W;
  const double t =(2281526.266/60./60.)/24./365;


  double norm = log(2)*N*t*(1.022)*(0.002/6289)*60*60*24*365;



  this->fHist->Scale(norm);
}

// ---------------------------------------------------------
gA_fit::gA_fit(const std::string& name,TString in_prefix,bool posS)
    : BCModel(name)
{

  SetInputs("/home/tdixon/Cd113Shape/CROSS_CWOnat/data_CWOnat.root","/home/tdixon/Cd113Shape/CROSS_CWOnat/histo_bkg_model.root");

  fBkg->Scale(1/(double)fBkg->Integral());

  // check binnings

  
  
  fTheory= new gATheoryHandler(Form("../output_%s.root",in_prefix.Data()));
  fTheory->Prepare("h");

  //check binning

  int Ndata =fData->GetNbinsX();
  int Nbkg = fBkg->GetNbinsX();
  int Ntheory = fTheory->fHist->GetNbinsX();

  std::cout<<"In data we have "<<Ndata<<" in bkg "<<Nbkg<<" and in theory "<<Ntheory<<" bins "<<std::endl;

  double maxdata = fData->GetBinCenter(fData->FindLastBinAbove(-1));
  double maxbkg = fBkg->GetBinCenter(fBkg->FindLastBinAbove(-1));
  double maxtheory = fTheory->fHist->GetBinCenter(fTheory->fHist->FindLastBinAbove(-1));

  std::cout<<"In data we have  max "<<maxdata<<" in bkg "<<maxbkg<<" and in theory "<<maxtheory<<" of histos "<<std::endl;

  double threshold =10;
  double maxt=1000;
  first= fData->FindBin(threshold);
  last =fData->FindBin(maxt);
  AddParameter("ga",0.5,1.3,"ga","[a.u]");

  if (posS)
    AddParameter("s",0,3,"S","[a.u]");
  else
    AddParameter("s",-3,0,"S","[a.u]");

  AddParameter("b",0,100000,"b","[cts]");
  AddObservable("T",0,2e16,"T_{1/2}","[yrs]");
  GetParameters().SetPriorConstantAll();
  GetParameters().SetNBins(1000);
  GetObservables().SetNBins(1000);

  double eta = 0.1222;
  double Na= 6.022e26;
  double W=360.24;
  double M = 0.433;

  
  N = eta*Na*M/W;
  std::cout<<"N = "<<N<<std::endl;

  t =(2281526.266/60./60.)/24./365;
  std::cout<<"t = "<<t*24*365<<std::endl;
}

// ---------------------------------------------------------
gA_fit::~gA_fit()
{
    // destructor
}


//----------------------------------------------------------

void gA_fit::ReconstructFit(TH1D *&h,TH1D *&hs,double ga,double s,double b)
{
  fTheory->CreateArbitaryGraph(ga,s,0);
  for (int i=fData->FindBin(0);i<last;i++)
    {

      double energy =fData->GetBinCenter(i);
      double sig = fTheory->fHist->GetBinContent(fTheory->fHist->FindBin(energy))*fEff->Eval(energy);
      double b_prob = fBkg->GetBinContent(fBkg->FindBin(energy));

      double lambda = sig+b*b_prob;

      h->SetBinContent(h->FindBin(energy),lambda);
      hs->SetBinContent(hs->FindBin(energy),sig);

    }


}
// ---------------------------------------------------------
double gA_fit::LogLikelihood(const std::vector<double>& pars)
{
  // return the log of the conditional probability p(data|pars).
  // This is where you define your model.
  // BCMath contains many functions you will find helpful.

  if (counter%100000==0)
    {
      verbose=1;
    }
  else
    {
      verbose=0;
    }
    
  double gA = pars[0];
  double s = pars[1];
  double b =pars[2];

  // model the data
  fTheory->CreateArbitaryGraph(gA,s,0);
  double logL=0;
  if (verbose)
    {
      std::cout<<"ga = "<<gA<<std::endl;
      std::cout<<"s  = "<<s<<std::endl;
    }
  for (int i=first;i<last;i++)
    {
      double energy =fData->GetBinCenter(i);
      double N = fData->GetBinContent(i);
      double b_prob = fBkg->GetBinContent(fBkg->FindBin(energy));
      double sig = fTheory->fHist->GetBinContent(fTheory->fHist->FindBin(energy));
      double lambda = sig*fEff->Eval(energy)+b*b_prob;
      if( lambda <= 0. ) continue;
      logL += -lambda + N * log( lambda ) - BCMath::LogFact( N );

      if (verbose &&i%100==0)
	{
	  std::cout<<"energy = "<<energy<<std::endl;
	  std::cout<<"N      = "<<N<<std::endl;
	  std::cout<<"b      = "<<b*b_prob<<std::endl;
	  std::cout<<"sig    = "<<sig<<std::endl;
	  std::cout<<"lambda = "<<lambda<<std::endl;
	  // std::cout<<" "<<std::endl;
	  std::cout<<"logL   = "<<logL<<std::endl;
	  std::cout<<" "<<std::endl;   
	}


    }
  counter++;
  return logL;
  return 0;
}






// ---------------------------------------------------------
// double gA_fit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

double gA_fit::GetT12(double ga,double s)
{
  fTheory->CreateArbitaryGraph(ga,s,0);

   double counts=0;
   for (int i=fData->FindBin(0);i<last;i++)
     {
       double energy=fData->GetBinCenter(i);
       counts+=fTheory->fHist->GetBinContent(fTheory->fHist->FindBin(energy));

     }
   if (verbose)
     {
       std::cout<<"T12 = "<<log(2)*N*t/counts<<std::endl;
       std::cout<<" "<<std::endl;
       std::cout<<" "<<std::endl;

     }
   return log(2)*N*t/counts;
}
// ---------------------------------------------------------
 void gA_fit::CalculateObservables(const std::vector<double>& pars)
 {

   GetObservable(0)= GetT12(pars[0],pars[1]);   
 }
