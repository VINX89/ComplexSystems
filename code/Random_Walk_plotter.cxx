#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TLegend.h"

#define N_STEPS 1000000
#define N_WALKERS 100000

void Random_Walk_plotter(){
  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk1d.txt","r");
  FILE *g = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk1d_bias.txt","r");
  FILE *h = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk2d_traj.txt","r");
  FILE *z = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk2d_prob.txt","r");

  int x_temp=0;
  int t_temp=0;
  double x[N_STEPS]={0};
  double t[N_STEPS]={0};
  double xmean[N_STEPS]={0};
  double xmean_err[N_STEPS]={0};
  double x2[N_STEPS]={0};
  double x2mean[N_STEPS]={0};
  double x2mean_err[N_STEPS]={0};
  double xb[N_STEPS]={0};
  double xmeanb[N_STEPS]={0};
  double xmean_errb[N_STEPS]={0};
  double x2b[N_STEPS]={0};
  double x2meanb[N_STEPS]={0};
  double x2mean_errb[N_STEPS]={0};
  double empty[N_STEPS]={0};
  
  int x_1_temp=0;
  int x_2_temp=0;
  
  int i=0;
  while(fscanf(f,"%d %d %lf %lf %lf %lf %lf\n",&t_temp,&x_temp,&xmean[i],&xmean_err[i],&x2[i],&x2mean[i],&x2mean_err[i])!=EOF){
    t[i]=(double)(t_temp);
    x[i]=(double)(x_temp);
    i++;
  }

  TCanvas *c1 = new TCanvas("c1","1d Random Walk: x(t)");
  c1->cd();
  TGraph *x_t = new TGraph(N_STEPS,t,x);
  x_t->Draw("ALP");
  x_t->GetXaxis()->SetTitle("t");
  x_t->GetYaxis()->SetTitle("x(t)");
  x_t->SetTitle("1d Random Walk: x(t)");
  x_t->SetLineColor(kBlue);
  x_t->SetMarkerColor(kBlue);
  c1->SetGridx();
  c1->SetGridy();
  c1->SetTickx();
  c1->SetTicky();

  TCanvas *c2 = new TCanvas("c2","1d Random Walk: <x(t)>");
  c2->cd();
  TGraphErrors *xmean_t = new TGraphErrors(N_STEPS,t,xmean,empty,xmean_err);
  xmean_t->Draw("ALP");
  xmean_t->GetXaxis()->SetTitle("t");
  xmean_t->GetYaxis()->SetTitle("<x(t)>");
  xmean_t->SetTitle("1d Random Walk: <x(t)>");
  xmean_t->SetLineColor(kBlue);
  xmean_t->SetMarkerColor(kBlack);
  c2->SetGridx();
  c2->SetGridy();
  c2->SetTickx();
  c2->SetTicky();


  TCanvas *c3 = new TCanvas("c3","1d Random Walk: x^{2}(t) (log-log)");
  c3->cd();
  TGraph *x2_t = new TGraph(N_STEPS,t,x2);
  x2_t->Draw("ALP");
  x2_t->GetXaxis()->SetTitle("t");
  x2_t->GetYaxis()->SetTitle("x^{2}(t)");
  x2_t->SetTitle("1d Random Walk: x^{2}(t)");
  x2_t->SetMarkerColor(kBlue);
  x2_t->SetLineColor(kBlue);
  x2_t->GetXaxis()->SetRangeUser(1E+0,1E+06);
  x2_t->GetYaxis()->SetRangeUser(1E+0,1E+06);
  c3->SetGridx();
  c3->SetGridy();
  c3->SetTickx();
  c3->SetTicky();
  c3->SetLogx();
  c3->SetLogy();
  TF1 *line = new TF1("line","x",0.0,1E+09);
  line->Draw("SAME");
  line->SetLineColor(kRed);
  TLegend *legend3 = new TLegend(0.6,0.7,0.8,0.8);
  legend3->AddEntry(x2_t,"x^{2}(t) (simulated)","L");
  legend3->AddEntry(line,"x^{2}(t)=t (theoretical)","L");
  legend3->Draw("SAME");

  TCanvas *c4 = new TCanvas("c4","1d Random Walk: <x^{2}(t)> (log-log)");
  c4->cd();
  TGraphErrors *x2mean_t = new TGraphErrors(N_STEPS,t,x2mean,empty,x2mean_err);
  x2mean_t->Draw("ALP");
  x2mean_t->GetXaxis()->SetTitle("t");
  x2mean_t->GetYaxis()->SetTitle("<x^{2}(t)>");
  x2mean_t->SetTitle("1d Random Walk: <x^{2}(t)>");
  x2mean_t->SetMarkerColor(kBlue);
  x2mean_t->SetLineColor(kBlue);
  x2mean_t->GetXaxis()->SetRangeUser(1E+0,1E+06);
  x2mean_t->GetYaxis()->SetRangeUser(1E+0,1E+06);
  c4->SetGridx();
  c4->SetGridy();
  c4->SetTickx();
  c4->SetTicky();
  c4->SetLogx();
  c4->SetLogy();
  TF1 *line2 = new TF1("line2","x",0.0,1E+09);
  line2->Draw("SAME");
  line2->SetLineColor(kRed);
  TLegend *legend4 = new TLegend(0.6,0.7,0.8,0.8);
  legend4->AddEntry(x2mean_t,"<x^{2}(t)> (simulated)","L");
  legend4->AddEntry(line2,"<x^{2}(t)>=t (theoretical)","L");
  legend4->Draw("SAME");

  i=0;
  while(fscanf(g,"%d %d %lf %lf %lf %lf %lf\n",&t_temp,&x_temp,&xmeanb[i],&xmean_errb[i],&x2b[i],&x2meanb[i],&x2mean_errb[i])!=EOF){
    xb[i]=(double)(x_temp);
    i++;
  }

  TCanvas *c5 = new TCanvas("c5","1d Random Walk: x(t) (bias)");
  c5->cd();
  TGraph *x_tb = new TGraph(N_STEPS,t,xb);
  x_tb->Draw("ALP");
  x_tb->GetXaxis()->SetTitle("t");
  x_tb->GetYaxis()->SetTitle("x(t)");
  x_tb->SetTitle("1d Random Walk: x(t) (p_{+}=0.52, p_{-}=0.48)");
  x_tb->SetLineColor(kBlue);
  x_tb->SetMarkerColor(kBlue);
  c5->SetGridx();
  c5->SetGridy();
  c5->SetTickx();
  c5->SetTicky();


  TCanvas *c6 = new TCanvas("c6","1d Random Walk: <x(t)> (log-log, bias)");
  c6->cd();
  TGraphErrors *xmean_tb = new TGraphErrors(N_STEPS,t,xmeanb,empty,xmean_errb);
  xmean_tb->Draw("ALP");
  xmean_tb->GetXaxis()->SetTitle("t");
  xmean_tb->GetYaxis()->SetTitle("<x(t)>");
  xmean_tb->SetTitle("1d Random Walk: <x(t)> (p_{+}=0.52, p_{-}=0.48)");
  xmean_tb->SetMarkerColor(kBlue);
  xmean_tb->SetLineColor(kBlue);
  xmean_tb->GetXaxis()->SetRangeUser(1E+0,1E+06);
  xmean_tb->GetYaxis()->SetRangeUser(1E+0,1E+06);
  c6->SetGridx();
  c6->SetGridy();
  c6->SetTickx();
  c6->SetTicky();
  c6->SetLogx();
  c6->SetLogy();
  TF1 *line3 = new TF1("line2","[0]*x",0.0,1E+09);
  line3->SetParameter(0,0.04);
  line3->Draw("SAME");
  line3->SetLineColor(kRed);
  TLegend *legend6 = new TLegend(0.6,0.7,0.8,0.8);
  legend6->AddEntry(xmean_tb,"<x(t)> (simulated)","L");
  legend6->AddEntry(line3,"<x(t)>=(p_{+}-p_{-})t=0.04t (theoretical)","L");
  legend6->Draw("SAME");


  TH2D *x1x2 = new TH2D("x1x2","2d Random Walk: trajectory",1400,-500.0,900.0,1600,-800.0,800.0); 
  while(fscanf(h,"%d %d %d\n",&t_temp,&x_1_temp,&x_2_temp)!=EOF){
    x1x2->Fill((double)x_1_temp,(double)x_2_temp);
  }
  TCanvas *c7 = new TCanvas("c7","2d Random Walk: trajectory");
  c7->cd();
  x1x2->Draw("SCAT");
  x1x2->GetXaxis()->SetTitle("x");
  x1x2->GetYaxis()->SetTitle("y");
  x1x2->SetMarkerColor(kBlue);
  c7->SetGridx();
  c7->SetGridy();
  c7->SetTickx();
  c7->SetTicky();

  TH2D *gauss2d_4 = new TH2D("gauss2d_4","2d Random Walk: positions distribution @ t=10^{4}",100,-1000.0,1000.0,100,-1000.0,1000.0);
  TH2D *gauss2d_5 = new TH2D("gauss2d_5","2d Random Walk: positions distribution @ t=10^{5}",100,-1000.0,1000.0,100,-1000.0,1000.0);
  TH2D *gauss2d_6 = new TH2D("gauss2d_6","2d Random Walk: positions distribution @ t=10^{6}",100,-1000.0,1000.0,100,-1000.0,1000.0);
  while(fscanf(z,"%d %d %d\n",&t_temp,&x_1_temp,&x_2_temp)!=EOF){
    if(t_temp==10000){
      gauss2d_4->Fill((double)x_1_temp,(double)x_2_temp);
    }
    else if(t_temp==100000){
      gauss2d_5->Fill((double)x_1_temp,(double)x_2_temp);
    }
    else{
      gauss2d_6->Fill((double)x_1_temp,(double)x_2_temp);
    }
  }
  TCanvas *c8 = new TCanvas("c8","2d Random Walk: positions distribution @ t=10^{4}");
  c8->cd();
  gStyle->SetPalette(1);
  gauss2d_4->Scale(1/gauss2d_4->Integral("width"));
  gauss2d_4->Draw("LEGO2");
  gauss2d_4->GetZaxis()->SetTitle("Relative frequency/(0.0025)");
  TCanvas *c9 =new TCanvas("c9","2d Random Walk: positions distribution @ t=10^{5}");
  c9->cd();
  gStyle->SetPalette(1);
  gauss2d_5->Scale(1/gauss2d_5->Integral("width"));
  gauss2d_5->Draw("LEGO2");
  gauss2d_5->GetZaxis()->SetTitle("Relative frequency/(0.0025)");
  TCanvas *c10 =new TCanvas("c10","2d Random Walk: positions distribution @ t=10^{6}");
  c10->cd();
  gStyle->SetPalette(1);
  gauss2d_6->Scale(1/gauss2d_6->Integral("width"));
  gauss2d_6->Draw("LEGO2");
  gauss2d_6->GetZaxis()->SetTitle("Relative frequency/(0.0025)");
  fclose(z);

  fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk2d_prob.txt","r");
  TH1D *gauss1d_4 = new TH1D("gauss1d_4","2d Random Walk: x positions distribution",100,-1000.0,1000.0);
  TH1D *gauss1d_5 = new TH1D("gauss1d_5","2d Random Walk: x positions distribution",100,-1000.0,1000.0);
  TH1D *gauss1d_6 = new TH1D("gauss1d_6","2d Random Walk: x positions distribution",100,-1000.0,1000.0);
  while(fscanf(z,"%d %d %d\n",&t_temp,&x_1_temp,&x_2_temp)!=EOF){
    if(t_temp==10000){
      gauss1d_4->Fill((double)x_1_temp);
    }
    else if(t_temp==100000){
      gauss1d_5->Fill((double)x_1_temp);
    }
    else{
      gauss1d_6->Fill((double)x_1_temp);
    }
  }
  TCanvas *c11 = new TCanvas("c11","2d Random Walk: x positions distributions");
  c11->cd();
  gauss1d_4->Sumw2();
  gauss1d_5->Sumw2();
  gauss1d_6->Sumw2();
  gauss1d_4->Scale(1/gauss1d_4->Integral("width"));
  gauss1d_5->Scale(1/gauss1d_5->Integral("width"));
  gauss1d_6->Scale(1/gauss1d_6->Integral("width"));
  
  TF1 *gauss4 = new TF1("gauss4","gaus(0)");
  gauss4->SetLineColor(kBlack);
  gauss4->SetParameter(1,gauss1d_4->GetMean());
  gauss4->SetParameter(2,gauss1d_4->GetRMS());
  TF1 *gauss5 = new TF1("gauss5","gaus(0)");
  gauss5->SetLineColor(kBlue);
  gauss5->SetParameter(1,gauss1d_5->GetMean());
  gauss5->SetParameter(2,gauss1d_5->GetRMS());
  TF1 *gauss6 = new TF1("gauss6","gaus(0)");
  gauss6->SetLineColor(kRed);
  gauss6->SetParameter(1,gauss1d_6->GetMean());
  gauss6->SetParameter(2,gauss1d_6->GetRMS());
  gauss1d_4->Fit("gauss4");
  gauss1d_5->Fit("gauss5");
  gauss1d_6->Fit("gauss6");

  gauss1d_4->Draw("E1");
  gauss1d_5->Draw("E1SAME");
  gauss1d_6->Draw("E1SAME");
  gauss1d_4->SetLineColor(kBlack);
  gauss1d_5->SetLineColor(kBlue);
  gauss1d_6->SetLineColor(kRed);
  gauss1d_4->GetXaxis()->SetTitle("x");
  gauss1d_4->GetYaxis()->SetTitle("Relative frequency/(0.05)");
  gauss1d_4->GetYaxis()->SetTitleOffset(1.5);

  TLegend *legend = new TLegend(0.6,0.7,0.8,0.8);
  legend->AddEntry(gauss4,"t=10^{4}","L");
  legend->AddEntry(gauss5,"t=10^{5}","L");
  legend->AddEntry(gauss6,"t=10^{6}","L");
  legend->Draw("SAME");

  gStyle->SetOptStat(0);
  

}
int main(){
  Random_Walk_plotter();
}
