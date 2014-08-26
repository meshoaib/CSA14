#include <TFile.h>
#include <TROOT.h>
#include <TH1.h>
#include "UserCode/TopAnalysis/interface/MiniEvent.h"
 
#include <vector>
#include <iostream>

void ReadTree(TString filename,TString output)
{
  gROOT->Reset();

  TH1F *jetpt = new TH1F("jetpt","jet pt " ,100,0.,300.);
  TH1F *jeteta = new TH1F("jeteta","jet eta" ,100,0.,3.);
//  TH1F *h2 = new TH1F("h2","csv" ,100,0.,300.);
  TH1F *ptofthirdjet = new TH1F("ptofthirdjet","third jet pt" ,100,0.,300.);
  TH1F *h4 = new TH1F("h4","cut flow" ,4,0.,4.);
  
  //read tree from file
  MiniEvent_t ev;
  TFile *f = TFile::Open(filename);
  TTree *t = (TTree*)f->Get("demo/AnaTree");
  attachToMiniEventTree(t,ev);

  //fill histograms, loop over all entries
//uint32_t nCSVMtags(0);
float pt=0.0,eta=0.0, thirdjetpt=0.0, csv=0.0;

  Int_t nentries = (Int_t)t->GetEntriesFast();
  for (Int_t i=0;i<nentries;i++)
    {
      t->GetEntry(i);
    for (int j=0; j<ev.nj;j++)
      {
//check pt and eta of this jet
	 pt = ev.j_pt[j];
         eta = ev.j_eta[j];
       	

	if (pt > 30 && abs(eta) < 2.5)        
	{
	    if (j >= 3)
	    {
	    thirdjetpt = ev.j_pt[j];
            h4->Fill(0);
	    }
            if(j >=4)
  	    h4->Fill(1);

//float csv=j.bDiscriminator("combinedSecondaryVertexBJetTags");
        csv = ev.j_csv[j];
if (csv>0.679) //nCSVMtags++;
//{
if(csv >=1)
h4->Fill(2);
if(csv >=2)
h4->Fill(3);
	}
}
ptofthirdjet->Fill(thirdjetpt);

	jetpt->Fill(pt);
        jeteta->Fill(eta); 
//	h3->Fill(jetpt);
//        h4->Fill(0);
//	h4->Fill(1);
    }
  
  //close file
  f->Close();

  //open output file
  TFile *fOut=TFile::Open(output+"/"+filename,"RECREATE");
  jetpt->Write();
  jeteta->Write();
//  h2->Write();
  ptofthirdjet->Write();
  h4->Write();
  fOut->Close();
}

