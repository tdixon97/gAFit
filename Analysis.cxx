#include "TTree.h"
#include "TString.h"
#include "TF1.h"
#include "TFile.h"
#include <map>
#include <tuple>
#include "TGraph.h"
#include <chrono>
using namespace std;

void bilinear_interpolate(double &f,double x,double y,double x1,double x2,double y1,double y2,double fQ11,double fQ12,double fQ21,double fQ22)
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

void GetLowAndHigh(double val,std::vector<double> vals,double &low,double &high,bool verbose=true)
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

void CreateArbitaryGraph(double ga,double sNME,TTree *&T)
{

  // ingredients we need
  // 1) a list of gA
  // 2) a list of sNME
  // 3) the input file

  double one_ga=0.5;
  double one_s=0;
  double one_E=1.022;

  int N=T->Draw("ga",Form("energy==%f &&sNME==%f",one_E,one_s),"goff");
 
  std::vector<double>ga_s(T->GetV1(),T->GetV1()+N);
  N=T->Draw("sNME",Form("energy==%f &&ga==%f",one_E,one_ga),"goff");
  std::vector<double>s_s(T->GetV1(),T->GetV1()+N);

  N=T->Draw("energy",Form("ga==%f &&sNME==%f",one_ga,one_s),"goff");
  std::vector<double> E_s(T->GetV1(),T->GetV1()+N);
  
  // now lets call gA axis x and sNME y

  // we need to save the data as a simple c++ struct
  T->SetEstimate(T->GetEntries());
  N = T->Draw("ga","","goff");
  std::vector<double>ga_all(T->GetV1(),T->GetV1()+N);

  N = T->Draw("energy","","goff");
  std::vector<double>energy_all(T->GetV1(),T->GetV1()+N);

  N = T->Draw("sNME","","goff");
  std::vector<double>s_all(T->GetV1(),T->GetV1()+N);

  N = T->Draw("spectrum","","goff");
  std::vector<double>val_all(T->GetV1(),T->GetV1()+N);

  std::map<std::tuple<double,double,double>,double> mapofvalues;
 
  for (int i=0;i<N;i++)
    {
      double gA_tmp = ga_all[i];
      double s_tmp = s_all[i];
      double energy_tmp=energy_all[i];
      double val_tmp=val_all[i];

      mapofvalues[std::make_tuple(gA_tmp,s_tmp,energy_tmp)]=val_tmp;
    }
      
    
  TGraph * g= new TGraph();
  int counter=0;
  double low_ga,high_ga,low_s,high_s;
  GetLowAndHigh(ga,ga_s,low_ga,high_ga,1);
  GetLowAndHigh(sNME,s_s,low_s,high_s,1);
  auto start = chrono::steady_clock::now();
  for (auto  & E: E_s)
    {

      // now we have an energy, we have gA and sNME and we need the 4 relevant fQij

      double fQ11=mapofvalues[std::make_tuple(low_ga,low_s,E)];
      double fQ21=mapofvalues[std::make_tuple(high_ga,low_s,E)];
      double fQ12=mapofvalues[std::make_tuple(low_ga,high_s,E)];
      double fQ22=mapofvalues[std::make_tuple(high_ga,high_s,E)];

      std::cout<<"energy = "<<E<<std::endl;
      
      // now interpolate
      double val;
      bilinear_interpolate(val,ga,sNME,low_ga,high_ga,low_s,high_s,fQ11,fQ12,fQ21,fQ22);
      std::cout<<"val = "<<val<<std::endl;
      g->SetPoint(counter,E,val);
      counter++;

    }
  auto end = chrono::steady_clock::now();
  cout << "Elapsed time: "
       << chrono::duration_cast<chrono::microseconds>(end - start).count()
       << " us" << endl;
  g->SetTitle(Form("Spectrum for gA = %f, sNME = %f; Energy [keV]; Probability [a.u] ; ",ga,sNME));
  g->Draw("AC");

  





  
}
  
  


    
int run()
{
  TFile *f = new TFile("output.root");
  TTree * T =(TTree*)f->Get("spectrum");
  
  CreateArbitaryGraph(0.76253,1.562,T);
  return 1;
}
