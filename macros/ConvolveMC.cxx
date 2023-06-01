#include "TH1D.h"
#include "TTree.h"
#include "TF1.h"
#include "TChain.h"
#include "TRandom3.h"
#include <getopt.h>
#include "TString.h"
#include <map>
#include <vector>
#include "TFile.h"
#include <fstream>
#include <iostream>
#include "TCanvas.h"
#include "TStyle.h"

void ConvolveMC(TString name,TString mc_path,TString reso_path,TString outpath,TString prefix,TString plot_dir,TString reso_label)
{
  

  // get the MC

  double Elow=0;
  double Ehigh=1022;
  double bin=1.022;

  double Ein_mc; 
  double Eout_mc;

  //  TFile *fMC = new TFile(mc_path);
  TChain *TMC = new TChain("Ntuple0");
  TMC->Add(mc_path);
  TMC->SetBranchAddress("primary_energy",&Ein_mc);
  TMC->SetBranchAddress("crystal_cwo_nat_totalEdep",&Eout_mc);

  TFile *freso = new TFile(reso_path);
  TF1 *fun_res = (TF1*)freso->Get(Form("Resolution_%s",reso_label.Data()));
  
  
  double E1=0; double E2=bin;

  std::map<int,TH1D*> map_of_hist;

  TRandom3 *rand =new TRandom3();
  
  
  
  while (E2<Ehigh)
    {
     
      TH1D * h  =new TH1D(Form("h_%i",(int)(round(E1/bin))),Form("h_%i",(int)(round(E1/bin))),1000,0,1000);
      h->SetTitle(Form("histogram for energy %f, bin %i ; Energy [keV] ; counts/keV; ",E1,(int)(round(E1/bin))));
      map_of_hist[(int)(round(E1/bin))]=h;
      E1+=bin;
      E2+=bin;
      
    }


  for (long int i=0;i<TMC->GetEntries();i++)
    {

      TMC->GetEntry(i);
      double Elow = bin*trunc((Ein_mc-bin/2.)/bin)+bin/2.;
      double Ehigh = Elow+bin;
      double Eoutput = rand->Gaus(Eout_mc,fun_res->Eval(Eout_mc)/2.355);
      map_of_hist[(int)(round((Elow-bin/2.)/bin))]->Fill(Ein_mc);
      if (i%1000000==0)
	{
	  std::cout<<100*(double)i/(double)TMC->GetEntries()<<" %"<<std::endl;
	}
    }

  // now draw
  TCanvas *can = new TCanvas();
  gStyle->SetOptStat(0);
  can->Print(Form("%s/plots_%s.pdf(",plot_dir.Data(),name.Data()),"pdf");
  int c=0;
  for (auto it =map_of_hist.begin();it!=map_of_hist.end();it++)
    {
      if (it->second->Integral()>0)
	it->second->Scale(1/(double)it->second->Integral());
      it->second->GetXaxis()->SetRangeUser(it->second->GetBinCenter(it->second->GetMaximumBin())-4,
					   it->second->GetBinCenter(it->second->GetMaximumBin())+4);
      it->second->Draw("HIST");
      
      can->Draw();

      
      if (c%1000==0)
	{
	  can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");
	}
      c++;
      

    }
  can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");


  
  
  // get the input file

  
  std::vector<std::string> small_matrix{"0",
                                        "-0.5",
                                        "0.5",
                                        "-1.0",
                                        "1.0",
                                        "-1.5",
                                        "1.5",
                                        "1.55",
                                        "-1.60",
                                        "1.60",
                                        "-1.65",
                                        "1.65",
                                        "-1.70",
                                        "1.70",
                                        "-1.75",
                                        "1.75",
                                        "-1.80",
                                        "1.80",
                                        "-1.85",
                                        "1.85",
                                        "-1.90",
                                        "1.90",
                                        "-1.95",
                                        "1.95",
                                        "-2.0",
                                        "2.0",
                                        "-2.05",
                                        "-2.10",
                                        "-2.5",
                                        "2.5",
                                        "-3.0",
                                        "3.0",
                                        "-3.5",
                                        "-4.0",
                                        };


  
  std::vector<double> ga;
  double ga_low=0.5; double ga_high=1.59; double ga_step=0.01;
  double ga_tmp=ga_low;
  while (ga_tmp<=ga_high)
    {
      ga.push_back(ga_tmp);
      ga_tmp+=ga_step;
    }

  TFile *f= new TFile(outpath+Form("/%s_%s_%s.root",name.Data(),reso_label.Data(),prefix.Data()),"recreate");
  TTree * T = new TTree("spectrum","spectrum");
  double g,s,E,v;

  T->Branch("ga",&g);
  T->Branch("sNME",&s);
  T->Branch("energy",&E);
  T->Branch("spectrum",&v);

  TTree * T_conv = new TTree("spectrum_conv","spectrum_conv");
  double g_conv,s_conv,E_conv,v_conv;

  T_conv->Branch("ga",&g_conv);
  T_conv->Branch("sNME",&s_conv);
  T_conv->Branch("energy",&E_conv);
  T_conv->Branch("spectrum",&v_conv);

  TH1D * htmp = new TH1D("htmp","htmp",1000,0,1000);
  std::fstream file;
  int count=0;
  can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");
  for (auto &ga_tmp_i:ga)
    {
      for (auto & s_tmp: small_matrix)
        {
          file.open(Form("Input_%s/%s/%1.2f.dat",prefix.Data(),s_tmp.data(),ga_tmp_i));


	  count++;
	  if (count%100==0)
	    std::cout<<"opening file"<<Form("Input/%s/%1.2f.dat",s_tmp.data(),ga_tmp_i)<<std::endl;
	  

	  // now we have the file

	  
	  htmp->Reset();
	  htmp->Clear();

	  c=0;
          while (file.is_open() &&!file.eof())
            {
	      c++;
	      // for the raw spectrun
              file>>E;
              file>>v;
              g=ga_tmp_i;
              s=std::stod(s_tmp);
              T->Fill();
	      int bin2 = (int)round(E/bin)-1;
	      if (map_of_hist[bin2]->Integral()>0)
		
		map_of_hist[bin2]->Scale(v/(double)map_of_hist[bin2]->Integral());
		  map_of_hist[bin2]->Draw("HIST");

		  if (ga_tmp_i==0.51 &&s_tmp=="0" &&c==100)
		    can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");
		  
		  htmp->Add(map_of_hist[bin2]);
		  htmp->Draw("HIST");
		  if (ga_tmp_i==0.51 &&s_tmp=="0"&&c==100 )
		    can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");
		
	      
	    }


	  
	  for (int b=1;b<htmp->GetNbinsX();b++)
	    {
	      E_conv = htmp->GetXaxis()->GetBinCenter(b);
	      g_conv=g;
	      s_conv=s;
	      v_conv=htmp->GetBinContent(b);
	      T_conv->Fill();
	    }


	  htmp->Draw("HIST");
	  if (ga_tmp_i==0.51 &&s_tmp=="0")
	    can->Print(Form("%s/plots_%s.pdf",plot_dir.Data(),name.Data()),"pdf");

	  
	  file.close();

        }
      
    }
  can->Print(Form("%s/plots_%s.pdf)",plot_dir.Data(),name.Data()),"pdf");

  f->Write();

  
















}


void Usage()
{

  std::cout << std::endl << std::endl;
    std::cout << "---------------------------------------------------------------------------------" << std::endl;
    if( std::getenv("USER") != NULL )
        std::cout << "Hi " << std::getenv("USER") <<"! The usage of this wonderful program is: " << std::endl;
    else
      std::cout << "Hi there! The usage of this wonderful program is:" << std::endl;
    std::cout<<"./Convolve_MC"<<std::endl;
    std::cout<<"This macro convolves the detector response with the beta spectrum templates, the inputs are set on command line. The output is a root file in {OUTPATH}/output_name.root"<<std::endl;
    std::cout<<"options "<<std::endl;
    std::cout<<"------------------------------------------------------------"<<std::endl;
    std::cout<<"-n (or --name             [name ]              (default: def)"<<std::endl;
    std::cout<<"-b (or --reso-label)      [type of reso sqrt or linear) (default: sqrt)"<<std::endl;
    std::cout<<"-m (or --mc-path)         [path to MC]         (default: None)"<<std::endl;
    std::cout<<"-r (or --reso-path)       [path for data]      (default: /home/tdixon/BATGraphHistoFit/inputs/cross/data_CWOnat.root) note: shouuld be after running eff fit"<<std::endl;
    std::cout<<"-o (or --out-path)        [path for out ]      (default: output/ "<<std::endl;
    std::cout<<"-p (or --prefix)          [which NME model]    (default: shell) "<<std::endl;
    std::cout<<"-l (or --plot-dir)        [directory for plot] (default: plots/)"<<std::endl;
    std::cout<<"-h (or --help)  "<<std::endl;
}


int main(int argc,char**argv)
{

  TString mc_path="/home/tdixon/Cd113Shape/data/cwa_nat_uniform/rootFiles/*";
  TString reso_path="/home/tdixon/BATGraphHistoFit/inputs/cross/data_CWOnat.root";
  TString name="def";
  TString out_path="gA_inputs";
  TString prefix="shell";
  TString plot_dir="plots/";
  TString reso_label="sqrt";

  {
    
    static struct option long_options[] =  {{ "name",required_argument,nullptr,'n'},
					    { "reso-label",required_argument,nullptr,'b'},
					    {"mc-path",        required_argument, nullptr, 'm'},
                                             {"reso-path", required_argument, nullptr, 'r'},
                                             {"out-path",   required_argument, nullptr, 'o'},
                                             {"prefix",    required_argument, nullptr, 'p'},
                                             {"plot-dir",    required_argument,nullptr,'l'},
					   { "help",                       no_argument, nullptr, 'h' },
                                             {nullptr, 0, nullptr, 0}
      };
    const char* const short_options = "n:b:m:r:o:p:l:h";
    int c;
    while ((c = getopt_long(argc, argv, short_options, long_options, nullptr)) != -1 ) {
      switch (c) {

      case 'n':{
	name=optarg;
	break;
      }
      case 'b':{
	reso_label=optarg;
	break;
      }
      case 'm': {
        mc_path = optarg;
        break;
      }
      case 'r': {
        reso_path = optarg;
        break;
      }
      case 'o': {
        out_path = optarg;
        break;
      }
      case 'p':{
	prefix=optarg;
	break;
      }
      case 'l': {
        plot_dir = optarg;
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
  

  
  

  ConvolveMC(name,mc_path,reso_path,out_path,prefix,plot_dir,reso_label);

  return 1;






}
