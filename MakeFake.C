void MakeFake(TString out)
{

  TString reso_path = "CROSS_CWOnat/data_CWOnat.root";
  TFile *fReso = new TFile(reso_path);

  TF1 * funcReso = (TF1*)fReso->Get("fit2_resolution");

  TFile * f =new TFile(out,"recreate");

  TTree * T = new TTree("mc_tree","mc_tree");

  double Ein,Eout;
  T->Branch("Ein",&Ein);
  T->Branch("Eout",&Eout);

  

  int N=1e9;
  double Emax=500;
  
  TRandom3 *rand = new TRandom3(0);
  
  for (int i=0;i<N;i++)
    {
      Ein=rand->Rndm()*Emax;

      double reso=funcReso->Eval(Ein)/2.355;

      Eout=rand->Gaus(Ein,reso);

      T->Fill();
    }


  f->Write();



}
