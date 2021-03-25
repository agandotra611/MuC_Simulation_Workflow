#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TCut.h>
#include <TMath.h>
#include <TPad.h>

void overlay(){

gROOT->SetBatch(1);

TCanvas *c1 = new TCanvas("c1","",800,800);

//root files with branches from Invariant_mass.C
TFile *f0 = new TFile("signal_vvHH_REC_1k_tree.root");
TFile *f1 = new TFile("bkg4b_REC_1k_test.root");
TFile *f2 = new TFile("bkg2bh_REC_1k_tree.root");
TFile *f3 = new TFile("bkgz_REC_1k_tree.root");

//Trees from the root files
TTree *t0 = (TTree*)f0->Get("quadCandTree");
TTree *t1 = (TTree*)f1->Get("quadCandTree");
TTree *t2 = (TTree*)f3->Get("quadCandTree");
TTree *t3 = (TTree*)f3->Get("quadCandTree");

//Stack plots
THStack *hs1 = new THStack("hs1","Invariant Mass of Jet Pair with Highest P_{T}");
THStack *hs2 = new THStack("hs2","Invariant Mass of Jet Pair with Lowest P_{T}");

//histograms (of signal and each of the bkgs for m1 and m2) to be stacked
TH1F *h0_m1 = new TH1F("h0_m1","signal_mmTOvvHH",25,0,250);
TH1F *h1_m1 = new TH1F("h1_m1","bkg1_mmTOvvbbbb",25,0,250);
TH1F *h2_m1 = new TH1F("h2_m1","bkg2_mmTOvvbbH",25,0,250);
TH1F *h3_m1 = new TH1F("h3_m1","bkg3_mmTOvvbbZ",25,0,250);

TH1F *h0_m2 = new TH1F("h0","signal_mmTOvvHH",25,0,250);
TH1F *h1_m2 = new TH1F("h1","bkg1_mmTOvvbbbb",25,0,250);
TH1F *h2_m2 = new TH1F("h2","bkg2_mmTOvvbbH",25,0,250);
TH1F *h3_m2 = new TH1F("h3","bkg3_mmTOvvbbZ",25,0,250);

//for error bars in ratio plot
h0_m1->Sumw2();
h1_m1->Sumw2();
h2_m1->Sumw2();
h3_m1->Sumw2();
h0_m2->Sumw2();
h1_m2->Sumw2();
h2_m2->Sumw2();
h3_m2->Sumw2();

//declare local variables
Float_t n_jets_ev, m1, pt1, m2, pt2;

//link local variables to corresponding branches (for making cuts and filling the histograms)
t0->SetBranchAddress("n_jets_ev",&n_jets_ev);
t0->SetBranchAddress("m1",&m1);
t0->SetBranchAddress("pt1",&pt1);
t0->SetBranchAddress("m2",&m2);
t0->SetBranchAddress("pt2",&pt2);

t1->SetBranchAddress("n_jets_ev",&n_jets_ev);
t1->SetBranchAddress("m1",&m1);
t1->SetBranchAddress("pt1",&pt1);
t1->SetBranchAddress("m2",&m2);
t1->SetBranchAddress("pt2",&pt2);

t2->SetBranchAddress("n_jets_ev",&n_jets_ev);
t2->SetBranchAddress("m1",&m1);
t2->SetBranchAddress("pt1",&pt1);
t2->SetBranchAddress("m2",&m2);
t2->SetBranchAddress("pt2",&pt2);

t3->SetBranchAddress("n_jets_ev",&n_jets_ev);
t3->SetBranchAddress("m1",&m1);
t3->SetBranchAddress("pt1",&pt1);
t3->SetBranchAddress("m2",&m2);
t3->SetBranchAddress("pt2",&pt2);

//get number of entries from each root file for looping on them
Long64_t nentries_sig = t0->GetEntries();
Long64_t nentries_bkg1 = t1->GetEntries();
Long64_t nentries_bkg2 = t2->GetEntries();
Long64_t nentries_bkg3 = t3->GetEntries(); 

//In each loop: fill the histograms for the corresponding root file, and make cuts

for (Long64_t i=0;i<nentries_sig;i++){
  t0->GetEntry(i);
  if(m1>100 && pt1>10 && pt2>10){
  h0_m1->Fill(m1);
  h0_m2->Fill(m2);
  }
  else{continue;}
}

for (Long64_t i=0;i<nentries_bkg1;i++){
  t1->GetEntry(i);
  if(m1>100 && pt1>10 && pt2>10){
  h1_m1->Fill(m1);
  h1_m2->Fill(m2);
  }
  else{continue;}
}

for (Long64_t i=0;i<nentries_bkg2;i++){
   t2->GetEntry(i);
   if(m1>100 && pt1>10 && pt2>10){
   h2_m1->Fill(m1);
   h2_m2->Fill(m2);
   }
   else{continue;}
}
for (Long64_t i=0;i<nentries_bkg3;i++){
   t3->GetEntry(i);
   if(m1>100 && pt1>10 && pt2>10){
   h3_m1->Fill(m1);
   h3_m2->Fill(m2);
   }
   else{continue;}
}


//xsec normalization

//bkg1: mm->vvbbbb
h1_m1->Scale(0.74);
h1_m2->Scale(0.74);
//bkg2: mm->vvbbH
h2_m1->Scale(3.74);
h2_m2->Scale(3.74);
//bkg3: mm->vvbbZ
h3_m1->Scale(31.25);
h3_m2->Scale(31.25);


//Set histogram style
h0_m1->SetFillStyle(0);
h0_m2->SetFillStyle(0);
h1_m1->SetFillStyle(1001);
h1_m2->SetFillStyle(1001);
h2_m1->SetFillStyle(1001);
h2_m2->SetFillStyle(1001);
h3_m1->SetFillStyle(1001);
h3_m2->SetFillStyle(1001);

h1_m1->SetFillColor(kCyan+3);
h1_m2->SetFillColor(kCyan+3);
h2_m1->SetFillColor(kCyan-6);
h2_m2->SetFillColor(kCyan-6);
h3_m1->SetFillColor(kCyan-9);
h3_m2->SetFillColor(kCyan-9);


//Add histograms to stack plot
hs1->Add(h0_m1);
hs1->Add(h1_m1);
hs1->Add(h2_m1);
hs1->Add(h3_m1);

hs2->Add(h0_m2);
hs2->Add(h1_m2);
hs2->Add(h2_m2);
hs2->Add(h3_m2);



//ratio plot

//clone signal histograms
TH1F *h0_m1_c1 = (TH1F*) h0_m1->Clone();
TH1F *h0_m2_c1 = (TH1F*) h0_m2->Clone();

//clone histograms to fill with bkgs
TH1F *m1_S_B = (TH1F*) h0_m1->Clone();
TH1F *m2_S_B = (TH1F*) h0_m2->Clone();

int n_bins;
float_t bin_content, bin_content_s, bin_content_b;
n_bins = h0_m1->GetNbinsX();

//loop over each of the bins in each of the bksgs; add content bin-by-bin to signal histo
for(int i=0; i<(n_bins+1);i++){

  bin_content = h1_m1->GetBinContent(i);
  m1_S_B->AddBinContent(i,bin_content);
  bin_content = h2_m1->GetBinContent(i);
  m1_S_B->AddBinContent(i,bin_content);
  bin_content = h3_m1->GetBinContent(i);
  m1_S_B->AddBinContent(i,bin_content);

  bin_content = h1_m2->GetBinContent(i);
  m2_S_B->AddBinContent(i,bin_content);
  bin_content = h2_m2->GetBinContent(i);
  m2_S_B->AddBinContent(i,bin_content);
  bin_content = h3_m2->GetBinContent(i);
  m2_S_B->AddBinContent(i,bin_content);

  bin_content=0.;
}

//divide each bin in the signal histo by sqrt(S+bkg1+bkg2+bkg3)
for(int i=0; i<(n_bins+1);i++){
  
  bin_content_s = h0_m1_c1->GetBinContent(i);
  bin_content_b = m1_S_B->GetBinContent(i);

  if(bin_content_s==0 || bin_content_b==0){
    h0_m1_c1->SetBinContent(i,0);
  }  
  else{
    bin_content = bin_content_s/sqrt(bin_content_b);
    h0_m1_c1->SetBinContent(i,bin_content);
  }
  
  bin_content_s = h0_m2_c1->GetBinContent(i);
  bin_content_b = m2_S_B->GetBinContent(i);
  
  if(bin_content_s==0 || bin_content_b==0){
    h0_m2_c1->SetBinContent(i,0);
  }
  else{
    bin_content = bin_content_s/sqrt(bin_content_b);
    h0_m2_c1->SetBinContent(i,bin_content);
  }

  //print ratio of 125 GeV bin
  if(i==12){
  printf("m1 125GeV bin ratio: %f \n", h0_m1_c1->GetBinContent(i));
  printf("m2 125GeV bin ratio: %f \n", h0_m2_c1->GetBinContent(i));
  }

  bin_content=0.; bin_content_s=0.; bin_content_b=0.;
}


//Stack histogram

//m1

//create pad for upper part (stack plot)
TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
pad1->Draw();
pad1->cd();
pad1->SetBottomMargin(0.05);
//draw stack plot
hs1->Draw("hist");
hs1->GetYaxis()->SetTitle("Number of Events");
hs1->GetXaxis()->SetTitle("");
hs1->GetXaxis()->SetLabelSize(0);
hs1->GetXaxis()->SetTitleSize(0);
//legend for m1 stack_histo
auto legend1 = new TLegend(0.6,0.6,0.9,0.9);
legend1->AddEntry(h0_m1,"Signal: #mu^{+} #mu^{-} #rightarrow #nu #nu H H","f");
legend1->AddEntry(h1_m1,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} b #bar{b}","f");
legend1->AddEntry(h2_m1,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} H","f");
legend1->AddEntry(h3_m1,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} Z","f");
legend1->Draw();

//exit pad1
c1->cd();

//create pad for lower part (ratio plot)
TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.3);
pad2->Draw();
pad2->cd();
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.25);
//set style for ratio plot
h0_m1_c1->SetMarkerStyle(kFullCircle);
h0_m1_c1->SetStats(0);
h0_m1_c1->Draw("hist p e1");

h0_m1_c1->GetXaxis()->SetLabelSize(0.1);
h0_m1_c1->GetXaxis()->SetTitleSize(0.1);
h0_m1_c1->GetXaxis()->SetTitle("M_{H_{1}} [GeV]");

h0_m1_c1->GetYaxis()->SetLabelSize(0.1);
h0_m1_c1->GetYaxis()->SetTitleSize(0.1);
h0_m1_c1->GetYaxis()->SetTitle("#frac{S}{#sqrt{S+B}} ratio");
h0_m1_c1->GetYaxis()->SetTitleOffset(0.4);

h0_m1_c1->SetTitle("");

//draw line on ratio plot for visual comparison
float xmin = h0_m1_c1->GetXaxis()->GetXmin();
float xmax = h0_m1_c1->GetXaxis()->GetXmax();
TLine *line = new TLine(xmin, 4.,xmax, 4.);
line->SetLineColor(kMagenta-10);
line->SetLineWidth(1.);
line->Draw();

c1->Update();
c1->SaveAs("m1_m1>100_pt10_ratio_cut.png");

////m2
c1->cd();
pad1->cd();

hs2->Draw("hist");
hs2->GetXaxis()->SetTitle("");
hs2->GetXaxis()->SetLabelSize(0);
hs2->GetXaxis()->SetTitleSize(0);
hs2->GetYaxis()->SetTitle("Number of Events");
//legend for m2 stack_histo
auto legend2 = new TLegend(0.6,0.6,0.9,0.9);
legend2->AddEntry(h0_m2,"Signal: #mu^{+} #mu^{-} #rightarrow #nu #nu H H","f");
legend2->AddEntry(h1_m2,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} b #bar{b}","f");
legend2->AddEntry(h2_m2,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} H","f");
legend2->AddEntry(h3_m2,"bkg: #mu^{+} #mu^{-} #rightarrow #nu #nu b #bar{b} Z","f");
legend2->Draw();

c1->cd();
pad2->cd();

h0_m2_c1->SetMarkerStyle(kFullCircle);
h0_m2_c1->SetStats(0);
h0_m2_c1->Draw("hist p e1");

h0_m2_c1->GetXaxis()->SetLabelSize(0.1);
h0_m2_c1->GetXaxis()->SetTitleSize(0.1);
h0_m2_c1->GetXaxis()->SetTitle("M_{H_{2}} [GeV]");

h0_m2_c1->GetYaxis()->SetLabelSize(0.1);
h0_m2_c1->GetYaxis()->SetTitleSize(0.1);
h0_m2_c1->GetYaxis()->SetTitle("#frac{S}{#sqrt{S+B}} ratio");
h0_m2_c1->GetYaxis()->SetTitleOffset(0.4);

h0_m2_c1->SetTitle("");

float xmin2 = h0_m2_c1->GetXaxis()->GetXmin();
float xmax2 = h0_m2_c1->GetXaxis()->GetXmax();
TLine *line2 = new TLine(xmin, 4.,xmax, 4.);
line2->SetLineColor(kMagenta-10);
line2->SetLineWidth(1.);
line2->Draw();

c1->Update();
c1->SaveAs("m2_m1>100_pt10_ratio_cut.png");

}
