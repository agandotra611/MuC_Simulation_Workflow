/*
Authors: A. Gandotra, T. Shemma, P. Das, I. Ojalvo

To run: 
root -l
.L quadCandTree.C++
.x quadCandTree.C("test.root")

or

root -b -l -q 'quadCandTree.C("test.root")'

nohup root -b -l -q quandCandTree.C(\"test.root\") & >log

FINISH ME: read in file to process from the commmand line
 */
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include <iostream>
#include "TString.h"

void quadCandTree(TString newFileName)
{

  //gROOT->SetBatch(1);

  //TFile *f1 = new TFile("/afs/cern.ch/work/t/tshemma/public/signal_bkg_1000/mmTOvvhh/root/b_test_1000.root","UPDATE");
  //TTree *MyLCTuple1 = (TTree*)f1->Get("MyLCTuple");
  // Dangerous to hardcode file names, should be input at the command line, google a tutorial for an example after making sure the rest of the code works
  TFile *fNew = new TFile(newFileName,"UPDATE");
  TTree *quadCandTree = new TTree("quadCandTree","new cand tree");

  Float_t n_jets_ev;
  Float_t m1;
  Float_t pt1;
  Float_t m2;
  Float_t pt2; 
  TLorentzVector h1, h2;
  TLorentzVector m1_leg1, m1_leg2, m2_leg1, m2_leg2;
  TLorentzVector jet1, jet2, jet3, jet4;
  Float_t m1_l1_pt, m1_l1_eta, m1_l1_phi, m1_l1_mass,
    m1_l2_pt, m1_l2_eta, m1_l2_phi, m1_l2_mass;

  Float_t m2l1_pt, m2l1_eta, m2l1_phi, m2l1_mass, m2l2_pt, m2l2_eta, m2l2_phi, m2l2_mass;

  Float_t jet1_pt, jet1_eta, jet1_phi, jet1_mass,  jet2_pt, jet2_eta, jet2_phi, jet2_mass, jet3_pt, jet3_eta, jet3_phi, jet3_mass, jet4_pt, jet4_eta, jet4_phi, jet4_mass;

  TBranch *n_jets_branch = quadCandTree->Branch("n_jets_ev", &n_jets_ev, "n_jets_ev/F");  
  TBranch *h1_mass_branch = quadCandTree->Branch("m1", &m1, "m1/F");
  TBranch *h1_pt_branch   = quadCandTree->Branch("pt1", &pt1, "pt1/F");
  TBranch *m1_leg1_pt_branch  = quadCandTree->Branch( "m1_l1_pt",  &m1_l1_pt,  "m1_l1_pt/F");
  TBranch *m1_leg1_eta_branch = quadCandTree->Branch("m1_l1_eta", &m1_l1_eta, "m1_l1_eta/F");
  TBranch *m1_leg1_phi_branch = quadCandTree->Branch("m1_l1_phi", &m1_l1_phi, "m1_l1_phi/F");
  TBranch *m1_leg2_pt_branch  = quadCandTree->Branch( "m1_l2_pt",  &m1_l2_pt,  "m1_l2_pt/F");
  TBranch *m1_leg2_eta_branch = quadCandTree->Branch("m1_l2_eta", &m1_l2_eta, "m1_l2_eta/F");
  TBranch *m1_leg2_phi_branch = quadCandTree->Branch("m1_l2_phi", &m1_l2_phi, "m1_l2_phi/F");
  

  TBranch *h2_mass_branch = quadCandTree->Branch("m2", &m2, "m2/F");
  TBranch *h2_pt_branch = quadCandTree->Branch("pt2", &pt2, "pt2/F");


  //  Add the rest of the branches, create variables and fill them
  TBranch *m2_leg1_pt_branch  = quadCandTree->Branch( "m2l1_pt",  &m2l1_pt,  "m2l1_pt/F");
  TBranch *m2_leg1_eta_branch = quadCandTree->Branch("m2l1_eta", &m2l1_eta, "m2l1_eta/F");
  TBranch *m2_leg1_phi_branch = quadCandTree->Branch("m2l1_phi", &m2l1_phi, "m2l1_phi/F");
  TBranch *m2_leg2_pt_branch  = quadCandTree->Branch( "m2l2_pt",  &m2l2_pt,  "m2l2_pt/F");
  TBranch *m2_leg2_eta_branch = quadCandTree->Branch("m2l2_eta", &m2l2_eta, "m2l2_eta/F");
  TBranch *m2_leg2_phi_branch = quadCandTree->Branch("m2l2_phi", &m2l2_phi, "m2l2_phi/F");
  
  TBranch *jet1_pt_branch   = quadCandTree->Branch( "jet1_pt",  &jet1_pt,  "jet1_pt/F");
  TBranch *jet1_eta_branch  = quadCandTree->Branch( "jet1_eta", &jet1_eta, "jet1_eta/F");
  TBranch *jet1_phi_branch  = quadCandTree->Branch( "jet1_phi", &jet1_phi, "jet1_phi/F");

  // Add the rest of the jet branches, create variables and fill them
  TBranch *jet2_pt_branch   = quadCandTree->Branch( "jet2_pt",  &jet2_pt,  "jet2_pt/F");
  TBranch *jet2_eta_branch  = quadCandTree->Branch( "jet2_eta", &jet2_eta, "jet2_eta/F");
  TBranch *jet2_phi_branch  = quadCandTree->Branch( "jet2_phi", &jet2_phi, "jet2_phi/F");
  
  TBranch *jet3_pt_branch   = quadCandTree->Branch( "jet3_pt",  &jet3_pt,  "jet3_pt/F");
  TBranch *jet3_eta_branch  = quadCandTree->Branch( "jet3_eta", &jet3_eta, "jet3_eta/F");
  TBranch *jet3_phi_branch  = quadCandTree->Branch( "jet3_phi", &jet3_phi, "jet3_phi/F");
  
  TBranch *jet4_pt_branch   = quadCandTree->Branch( "jet4_pt",  &jet4_pt,  "jet4_pt/F");
  TBranch *jet4_eta_branch  = quadCandTree->Branch( "jet4_eta", &jet4_eta, "jet4_eta/F");
  TBranch *jet4_phi_branch  = quadCandTree->Branch( "jet4_phi", &jet4_phi, "jet4_phi/F");
  

  // set mass of higgs to 125, try to match di-candidate pairs that are close to this mass
  double mh=125.;

  // Declaration of leaf types in tree
  Int_t           nj;
  Float_t         jmox[8], jmoy[8], jmoz[8], jmas[8], jene[8];

  // List of branches
  TBranch        *b_njet, *b_jmox, *b_jmoy, *b_jmoz, *b_jmas, *b_jene;  

  TChain* fChain = new TChain("fChain");
  // Dangerous to hardcode file names, should be input at the command line, google a tutorial for an example after making sure the rest of the code works
  fChain->Add("signal_vvHH_REC_1k_tree.root/MyLCTuple");
   
  fChain->SetBranchAddress("nj", &nj, &b_njet);
  // the following is one vector of jets
  fChain->SetBranchAddress("jmox", jmox, &b_jmox); //px
  fChain->SetBranchAddress("jmoy", jmoy, &b_jmoy); //py
  fChain->SetBranchAddress("jmoz", jmoz, &b_jmoz); //pz
  fChain->SetBranchAddress("jmas", jmas, &b_jmas); //jet mass
  fChain->SetBranchAddress("jene", jene, &b_jene); //jet energy
  
  double inv_mass[1000][1000];
  
  double delta_min=10000000;
  double delta_tot;
  double delta_m1;
  double delta_m2;
  int jet1_num=0;
  int jet2_num=0;
  int jet3_num=0;
  int jet4_num=0;

  //to test try setting nLoops to 10
  unsigned int nLoops =10;
  unsigned int iLoop =0;

  for(unsigned int ientry=0; ientry< fChain->GetEntries(); ++ientry)
    {

      fChain->GetEntry(ientry);
      //zero out all storage variables first      
      m1  = 0.0;
      pt1 = 0.0;
      m2  = 0.0; 
      pt2 = 0.0;
      m1_l1_pt=0, m1_l1_eta=0, m1_l1_phi=0, m1_l2_pt=0, m1_l2_eta=0, m1_l2_phi;
      jet1_pt=0, jet1_eta=0, jet1_phi=0;

      delta_min = 10000000;
      jet1_num=0;
      jet2_num=0;
      jet3_num=0;
      jet4_num=0;

      n_jets_ev=nj;
      
      /* remove for now
      for (int k=0;k<nj; k++)
	{ 
          double jmom=sqrt(jmox[k]*jmox[k]+jmoy[k]*jmoy[k]+jmoz[k]*jmoz[k]);          
	  eta_histo = (-TMath::Log(TMath::Tan(TMath::ACos(jmoz[k]/jmom)/2.)));
 	  pt_histo = (jmom/TMath::CosH(-TMath::Log(TMath::Tan(TMath::ACos(jmoz[k]/jmom)/2.))));
	  theta_histo = (TMath::ACos(jmoz[k]/jmom));
	  phi_histo = (TMath::ATan(jmoy[k]/jmox[k]));
	  E_histo = (jene[k]);
	}
      */

      if(nj<=3){continue;}// Continue without filling anything
      if(nj>3) //selection of events with more than four jets
	{
	  //remove if not a test
	  /*iLoop++;
	  if(iLoop>nLoops)
	  break;*/

	  for (int i=0; i<nj-1; i++) //Calculation of invariant mass 
	    {
	      for (int j=i+1; j<nj;j++)
		{
		  double pij_mod=sqrt(((jmox[i]+jmox[j])*(jmox[i]+jmox[j]))+((jmoy[i]+jmoy[j])*(jmoy[i]+jmoy[j]))+((jmoz[i]+jmoz[j])*(jmoz[i]+jmoz[j])));
		  inv_mass[i][j]=sqrt((jene[i]+jene[j])*(jene[i]+jene[j])-pij_mod*pij_mod);
		}
	    }
	  
	  for (int i=0; i<nj; i++)
	    { 
	      for (int j=i+1; j<nj;j++)
		{
		  for(int k=0; k<nj; k++)
		    {
		      for(int l=k+1; l<nj; l++)
			{   
			  if(k!=i &&k!=j &&l!=i &&l!=j)
			    {
			      delta_m1=mh-inv_mass[i][j]; 
			      delta_m2=mh-inv_mass[k][l]; 
			      delta_tot=(delta_m1*delta_m1)+(delta_m2*delta_m2); 
			      if (delta_tot<delta_min ){delta_min=delta_tot; jet1_num=i;jet2_num=j;jet3_num=k;jet4_num=l;}
			    } 
			} //end l
		    } //end k
		} //end j
	    } //end i

	  //Set the PT ordered jet lorentz vectors
	  jet1.SetPxPyPzE(jmox[0],jmoy[0],jmoz[0],jene[0]);
	  jet2.SetPxPyPzE(jmox[1],jmoy[1],jmoz[1],jene[1]);
	  jet3.SetPxPyPzE(jmox[2],jmoy[2],jmoz[2],jene[2]);
	  jet4.SetPxPyPzE(jmox[3],jmoy[3],jmoz[3],jene[3]);
	  
	  //ordering of H1 and H2
	  double PT_H1=sqrt((jmox[jet1_num]+jmox[jet2_num])*(jmox[jet1_num]+jmox[jet2_num])+(jmoy[jet1_num]+jmoy[jet2_num])*(jmoy[jet1_num]+jmoy[jet2_num]));
	  double PT_H2=sqrt((jmox[jet3_num]+jmox[jet4_num])*(jmox[jet3_num]+jmox[jet4_num])+(jmoy[jet3_num]+jmoy[jet4_num])*(jmoy[jet3_num]+jmoy[jet4_num]));
	  
	  if (PT_H1>=PT_H2){	    
	    m1 = (inv_mass[jet1_num][jet2_num]);
	    m2 = (inv_mass[jet3_num][jet4_num]); 
	    //set the 
	    m1_leg1.SetPxPyPzE(jmox[jet1_num], jmoy[jet1_num], jmoz[jet1_num], jene[jet1_num]);
	    m1_leg2.SetPxPyPzE(jmox[jet2_num], jmoy[jet2_num], jmoz[jet2_num], jene[jet2_num]);
	    m2_leg1.SetPxPyPzE(jmox[jet3_num], jmoy[jet3_num], jmoz[jet3_num], jene[jet3_num]);
	    m2_leg2.SetPxPyPzE(jmox[jet4_num], jmoy[jet4_num], jmoz[jet4_num], jene[jet4_num]);
	    //std::cout<<"m1 "<<m1<<std::endl;
	  }
	  else{  //(PT_H1<PT_H2) - (note from Isobel: removing this as it is superfluous)
	    m1 = (inv_mass[jet3_num][jet4_num]);
	    m2 = (inv_mass[jet1_num][jet2_num]); 

	    m1_leg1.SetPxPyPzE(jmox[jet3_num], jmoy[jet3_num], jmoz[jet3_num], jene[jet3_num]);
	    m1_leg2.SetPxPyPzE(jmox[jet4_num], jmoy[jet4_num], jmoz[jet4_num], jene[jet4_num]);
	    m2_leg1.SetPxPyPzE(jmox[jet1_num], jmoy[jet1_num], jmoz[jet1_num], jene[jet1_num]);
	    m2_leg2.SetPxPyPzE(jmox[jet2_num], jmoy[jet2_num], jmoz[jet2_num], jene[jet2_num]);
	  }
	  //now set the lorentz vector for the higgs candidates	  
	  h1 = (m1_leg1+m1_leg2);
	  h2 = (m2_leg1+m2_leg2);

	} //end nj>3
      //Set all the variables 
      //Higgs Cand 1
      m1 = h1.M();
      pt1 = h1.Pt();
      m1_l1_pt  = m1_leg1.Pt();
      m1_l1_eta = m1_leg1.Eta();
      m1_l1_phi = m1_leg1.Phi();
      m1_l2_pt  = m1_leg2.Pt();
      m1_l2_eta = m1_leg2.Eta();
      m1_l2_phi = m1_leg2.Phi();

      //Higgs Cand 2 
      m2 = h2.M();
      pt2 = h1.Pt();
      m2l1_pt = m2_leg1.Pt();
      m2l1_eta = m2_leg1.Eta();
      m2l1_phi = m2_leg1.Phi();
      m2l2_pt = m2_leg2.Pt();
      m2l2_eta = m2_leg2.Eta();
      m2l2_phi = m2_leg2.Phi();
      //Now the jets 
      //jet 1
      jet1_pt  = jet1.Pt();
      jet1_eta = jet1.Eta();
      jet1_phi = jet1.Phi();
     
      // jet 2
      jet2_pt  = jet2.Pt();
      jet2_eta = jet2.Eta();
      jet2_phi = jet2.Phi();

      // jet 3                                                                                                                                                                                              
      jet3_pt  = jet3.Pt();
      jet3_eta = jet3.Eta();
      jet3_phi = jet3.Phi();

      // jet 4                                                                                                                                                                                              
      jet4_pt  = jet4.Pt();
      jet4_eta = jet4.Eta();
      jet4_phi = jet4.Phi();

      printf("event no: %d ,n_jets: %f, m1: %f, m2: %f\n ", ientry, n_jets_ev, m1, m2);

      // And now fill the tree once per event
      quadCandTree->Fill();
      /*
      m1_branch->Fill();
      m1_leg1_pt_branch->Fill();
      m1_leg1_eta_branch->Fill();
      m1_leg1_phi_branch->Fill();
      m1_leg2_pt_branch->Fill();
      m1_leg2_eta_branch->Fill();
      m1_leg2_phi_branch->Fill();
      
      jet1_pt->Fill();
      */


    } // end of loop on fChain->GetEntries()

  fNew->cd();
  quadCandTree->Write();
  fNew->Write("", TObject::kOverwrite);

}
