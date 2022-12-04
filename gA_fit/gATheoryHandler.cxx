
#include "gATheoryHandler.hh"


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

      // now we have an energy, we have gA and sNME and we need the 4 relevant fQij                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

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

  // now we need to normalise, we want a spectrum in cts/kg-yr
  // integral of spectrum should be
  //Gamma = S/(eff*NumberOfAtoms*time)
  // ln(2)/T = S/(eff*N*t)
  // we are told by joel that
  // half-life (in seconds) = 6289 / (0.002*sum(spectrum))
  // note we changed the binning
  // or 1/T = sum(spectrum)*0.002/6289
  // norm 0 adjusts for the binning
  // spectrum^{(0)} = spectrum*1.022/1
  // so norm 1 will make the sum of the spectrum gamma
  // spectrum^{(1)} = spectrum^{(0)}*6829/(0.002)
  // now we have 1/T = sum(spectrum^{(1)} but we want
  // sum(spectrum^{2}) = counts
  // so we multiply by eff*N*t*ln(2)
  // spectrum^{2} = eff*N*t*ln(2)*spectrum^{1}
  // in practice efficiency is handled seperatly so our factor is


  double eta = 0.1222;
  double Na= 6.022e26;
  double W=360.24;
  double M = 0.433;


  const double N = eta*Na*M/W;
  const double t = (2281526.266/60./60.)/24./365;
  std::cout<<" N = "<<N<<std::endl;
  std::cout<<" t = "<<t<<std::endl;
  
  
  double norm = 60*60*24*365*N*t*log(2)*(1.022/1)*0.002/6289;

  std::cout<<"norm = "<<norm<<std::endl;
  

  this->fHist->Scale(norm);



  
















}



int run()
{
  gATheoryHandler *gA= new gATheoryHandler("../output_conv.root");

  double s=2.1;
  TCanvas *c = new TCanvas();
  gA->Prepare("h");
  TH2D * h2 = new TH2D("h2","h2",400,0,400,1000,0.5,1.2);
  TH2D *hT = new TH2D("hT","hT",100,0.5,1.2,100,-3,3);
  TH2D *hT2 = new TH2D("hT2","hT2",100,0.5,1.2,100,-3,3);

  for (int sn=0;sn<100;sn++)
    {
      s=-3+6*sn/100.;
  for (int i=0;i<1e2;i++)
    {
      int verbose=0;
      double ga=0.5+0.7*i*1e-2;

      if (i%100==0)
	{
	  verbose=1;
	  std::cout<<"Running for ga = "<<ga<<" "<<i<<" / 1e4 "<<std::endl;
	}
      gA->CreateArbitaryGraph(ga,s,verbose);
      auto start = chrono::steady_clock::now();

      gA->fHist->SetTitle(Form("Spectrum for gA = %.3f, sNME = %3f; Energy [keV]; Probability [a.u] ; ",ga,s));
      gA->fHist->GetXaxis()->SetRangeUser(0,400);
      gStyle->SetOptStat(0);

      gA->fHist->GetYaxis()->SetRangeUser(1e-22,3e-19);
      gA->fHist->Draw("HIST");

      c->SetLogy();
      c->Draw();
      c->SaveAs(Form("plots/plot_%i.pdf",i));
      
      auto end = chrono::steady_clock::now();
       if (verbose)
	 {
	   cout << "Elapsed time: "
		<< chrono::duration_cast<chrono::microseconds>(end - start).count()
		<< " us" << endl;
	 }
       for (int b=1;b<gA->fHist->GetNbinsX();b++)
	 {
	   h2->SetBinContent(b,i,gA->fHist->GetBinContent(b));
	 }
       double T12 = log(2)*8.858e22*(634./24./365.) / (gA->fHist->Integral());
       if (verbose)
	 {
	   std::cout<<" N  = "<<gA->fHist->Integral()<<std::endl;
	   std::cout<<"T12 = "<<T12<<std::endl;
	 }
       
       hT->SetBinContent(i,sn,T12);
       if (fabs(T12-8e15)<1e15)
	 {
	   hT2->SetBinContent(i,sn,T12);
	 }
       else
	 {
	   hT2->SetBinContent(i,sn,0);
	 }
    }
  h2->GetZaxis()->SetRangeUser(1e-22,3e-19);
  h2->SetTitle(Form(" sNME = %3f; Energy [keV]; ga ;  ",s));

  h2->SetContour(1000);
  h2->Draw("colz");
  c->SetLogy(0);
  c->SetLogz();
  c->Draw();
  c->Print("plots2/plot_%i.pdf");
  
    }

  hT->SetContour(1000);
  
  hT->SetTitle("sNME ; gA ;  ");
  hT2->SetContour(1000);

  hT2->SetTitle("sNME ; gA ;  ");

  c->SetLogy(0);
  c->SetLogz();
  c->Draw();
  c->SaveAs("hT.C");
  return 100;
}

