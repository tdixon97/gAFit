TH1D* TTree2Hist(TTree *&T,double gA,double S,TString name)
{

 
 int N = T->Draw("energy:spectrum",Form("ga==%f &&sNME==%f",gA,S),"goff");

 std::vector<double>energy(T->GetV1(),T->GetV1()+N);
 std::vector<double>spectrum(T->GetV2(),T->GetV2()+N);

 double e0=energy[0];
 double e1=energy[1];
 double bin=e1-e0;
 std::cout<<bin<<std::endl;
 std::cout<<"e0,e1 = "<<e0<<" , "<<e1<<std::endl;
 double low=e0-bin/2.;
 std::cout<<low<<std::endl;

 TH1D *h = new TH1D(name,name,10000,low,10000*bin);


 for (int i=0;i<N;i++)
   {
     double e_tmp = energy[i];
     double s_tmp=spectrum[i];

     h->SetBinContent(h->FindBin(e_tmp),s_tmp);
   }

 h->Draw("HIST");


 return h;
}
   
     
