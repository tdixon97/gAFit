

void Run(TString output_name="output.root")
{



  std::vector<std::string> small_matrix{"0",
					"-0.5",
					"0.5",
					"-1.0",
					"1.0",
					"-1.5",
					"1.5",
					"-1.55",
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
					"2.05",					
					"-2.10",
					"2.10",
					"-2.5",
					"2.5",
					"-3.0",
					"3.0",
					"-3.5",
					"3.5",
					"-4.0",
					"4.0"};

  

  std::vector<double> ga;
  double ga_low=0.5; double ga_high=1.59; double ga_step=0.01;
  double ga_tmp=ga_low;
  while (ga_tmp<=ga_high)
    {
      ga.push_back(ga_tmp);
      ga_tmp+=ga_step;
    }


  // make the input ROOT file

  TFile *f= new TFile(output_name,"recreate");
  TTree * T = new TTree("spectrum","spectrum");
  double g,s,E,v;

  T->Branch("ga",&g);
  T->Branch("sNME",&s);
  T->Branch("energy",&E);
  T->Branch("spectrum",&v);
  
  std::fstream file;
  for (auto &ga_tmp_i:ga)
    {
      for (auto & s_tmp: small_matrix)
	{
	  file.open(Form("Input/%s/%1.2f.dat",s_tmp.data(),ga_tmp_i));

	  std::cout<<"opening file"<<Form("Input/%s/%1.2f.dat",s_tmp.data(),ga_tmp_i)<<std::endl;

	  while (file.is_open() &&!file.eof())
	    {
	      file>>E;
	      
	      file>>v;
	      g=ga_tmp_i;
	      s=std::stod(s_tmp);
	      T->Fill();
	    }
	  file.close();
	      
	}
    }

  f->Write();
    
}
	  
