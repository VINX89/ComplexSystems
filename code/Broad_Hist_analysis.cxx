#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Float.h>
#include <iostream>
#include <cstdlib>

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
#include "TMath.h"
#include "TRandom.h"

#define threshold_inf 5
#define threshold_sup 1E+10
#define BETA_MIN 0.1
#define BETA_MAX 1.0

void Broad_Hist_analysis(){
  char path[200] = "/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/test_L8.txt";
  FILE *f = fopen(path,"r");
  int L=0;
  double buffer=0.0;
  double buffer_e=0.0;
  double E_MIN=0.0;
  double E_MAX=0.0;
  int i=0,j=0,k=0;
  while(fscanf(f,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&L,&buffer_e,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(i==threshold_inf){
      E_MIN = buffer_e;
      k=1;
    }
    else if(i==threshold_sup){
      E_MAX = buffer_e;
      k=0;
    }
    else{}
    if(k==1){
      j++;
    }
    i++;
  }
  fclose(f);
  int N_BINS = j;
  double BIN_SIZE = (E_MAX-E_MIN)/((double)(N_BINS));
  printf("Lattice size %d, # bins %d, bin size %lf, E MIN %lf, E MAX %lf\n",L,N_BINS,BIN_SIZE,E_MIN,E_MAX);
  TH1D *e_freq = new TH1D("e_freq","Energy Distribution",N_BINS,E_MIN,E_MAX);
  TH1D *N_down_4 = new TH1D("N_down_4","N_{down}(4)",N_BINS,E_MIN,E_MAX);
  TH1D *N_up_4 = new TH1D("N_up_4","N_{up}(4)",N_BINS,E_MIN,E_MAX);
  TH1D *N_down_8 = new TH1D("N_down_8","N_{down}(8)",N_BINS,E_MIN,E_MAX);
  TH1D *N_up_8 = new TH1D("N_up_8","N_{up}(8)",N_BINS,E_MIN,E_MAX);
  double *m;
  m = (double *)calloc(N_BINS,sizeof(double));
  double *m_err;
  m_err = (double *)calloc(N_BINS,sizeof(double));
  double *e;
  e = (double *)calloc(N_BINS,sizeof(double));
  double e_temp,e_freq_temp,m_temp,m_err_temp,nd4,nu4,nd8,nu8;

  fopen(path,"r");
  i=0;
  j=0;
  printf("Data:\n");
  while(fscanf(f,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",&L,&e_temp,&e_freq_temp,&m_temp,&m_err_temp,&nd4,&nu4,&nd8,&nu8)!=EOF){
    if(i>=threshold_inf && i<=threshold_sup){
      e[j] = e_temp;
      e_freq->Fill(e[j],e_freq_temp);
      m[j]=m_temp;
      m_err[j]=m_err_temp;
      N_down_4->Fill(e[j],nd4/(nd4+nu4));
      N_down_4->SetBinError(j+1,N_down_4->GetBinContent(j+1)*(sqrt(nd4)/nd4+(sqrt(nd4)+sqrt(nu4))/(nd4+nu4)));
      N_up_4->Fill(e[j],nu4/(nd4+nu4));
      N_up_4->SetBinError(j+1,N_up_4->GetBinContent(j+1)*(sqrt(nu4)/nd4+(sqrt(nd4)+sqrt(nu4))/(nd4+nu4)));
      N_down_8->Fill(e[j],nd8/(nd8+nu8));
      N_down_8->SetBinError(j+1,N_down_8->GetBinContent(j+1)*(sqrt(nd8)/nd8+(sqrt(nd8)+sqrt(nu8))/(nd8+nu8)));
      N_up_8->Fill(e[j],nu8/(nd8+nu8));
      N_up_8->SetBinError(j+1,N_up_8->GetBinContent(j+1)*(sqrt(nu8)/nd8+(sqrt(nd8)+sqrt(nu8))/(nd8+nu8)));
      printf("e %lf, m %lf +- %lf, N_down(4) %lf +- %lf, N_up(4) %lf +- %lf, N_down(8) %lf +- %lf, N_up(8) %lf +- %lf\n",e[j],m[j],m_err[j],N_down_4->GetBinContent(j+1),N_down_4->GetBinError(j+1),N_up_4->GetBinContent(j+1),N_up_4->GetBinError(j+1),N_down_8->GetBinContent(j+1),N_down_8->GetBinError(j+1),N_up_8->GetBinContent(j+1),N_up_8->GetBinError(j+1));
      j++;
    }
    i++;
  }
  
  //Plot data
  TCanvas *c1 = new TCanvas("c1","N_down_4 and N_up_4");
  c1->cd();
  N_down_4->Draw("E1");
  N_down_4->SetLineColor(kBlue);
  N_down_4->SetMarkerColor(kBlue);
  N_down_4->SetMarkerStyle(21);
  N_down_4->SetTitle("#DeltaE=4");
  N_down_4->GetXaxis()->SetTitle("E");
  N_down_4->GetYaxis()->SetTitle("P(E)");
  N_down_4->GetYaxis()->SetRangeUser(0.0,1.1*N_up_4->GetBinContent(N_up_4->GetMaximumBin()));
  N_up_4->Draw("SAMEE1");
  N_up_4->SetLineColor(kRed);
  N_up_4->SetMarkerColor(kRed);
  N_up_4->SetMarkerStyle(21);
  N_up_4->SetTitle("#DeltaE=4");
  N_up_4->GetXaxis()->SetTitle("E");
  N_up_4->GetYaxis()->SetTitle("P(E)");
  TLegend *leg1 = new TLegend(0.6,0.7,0.8,0.8);
  leg1->AddEntry(N_down_4,"P_{down}(E)=#frac{N_{down}(E)}{N_{down}(E)+N_{up}(E)}","L");
  leg1->AddEntry(N_up_4,"P_{up}(E)=#frac{N_{up}(E)}{N_{down}(E)+N_{up}(E)}","L");
  leg1->Draw("SAME");
  TCanvas *c2 = new TCanvas("c2","N_down_8 and N_up_8");
  c2->cd();
  N_down_8->Draw("E1");
  N_down_8->SetLineColor(kBlue);
  N_down_8->SetMarkerColor(kBlue);
  N_down_8->SetMarkerStyle(21);
  N_down_8->SetTitle("#DeltaE=8");
  N_down_8->GetXaxis()->SetTitle("E");
  N_down_8->GetYaxis()->SetTitle("P(E)");
  N_down_8->GetYaxis()->SetRangeUser(0.0,1.1*N_up_8->GetBinContent(N_up_8->GetMaximumBin()));
  N_up_8->Draw("SAMEE1");
  N_up_8->SetLineColor(kRed);
  N_up_8->SetMarkerColor(kRed);
  N_up_8->SetMarkerStyle(21);
  N_up_8->SetTitle("#DeltaE=8");
  N_up_8->GetXaxis()->SetTitle("E");
  N_up_8->GetYaxis()->SetTitle("P(E)");
  TLegend *leg2 = new TLegend(0.6,0.7,0.8,0.8);
  leg2->AddEntry(N_down_8,"P_{down}(E)=#frac{N_{down}(E)}{N_{down}(E)+N_{up}(E)}","L");
  leg2->AddEntry(N_up_8,"P_{up}(E)=#frac{N_{up}(E)}{N_{down}(E)+N_{up}(E)}","L");
  leg2->Draw("SAME");
  TCanvas *c3 = new TCanvas("c3","Energy distribution");
  c3->cd();
  e_freq->Draw();
  e_freq->GetXaxis()->SetTitle("E");
  e_freq->GetYaxis()->SetTitle("entries/4.0");
  e_freq->SetFillColor(kBlue);

  //Compute log g(E) and g(E)
  double *log_ge4;
  log_ge4 = (double *)calloc(N_BINS-1,sizeof(double));
  double *log_ge4_err;
  log_ge4_err = (double *)calloc(N_BINS-1,sizeof(double));
  double *ge4;
  ge4 = (double *)calloc(N_BINS-1,sizeof(double));
  double *ge4_err;
  ge4_err = (double *)calloc(N_BINS-1,sizeof(double));
  double *log_ge8;
  log_ge8 = (double *)calloc(N_BINS-1,sizeof(double));
  double *log_ge8_err;
  log_ge8_err = (double *)calloc(N_BINS-1,sizeof(double));
  double *ge8;
  ge8 = (double *)calloc(N_BINS-1,sizeof(double));
  double *ge8_err;
  ge8_err= (double *)calloc(N_BINS-1,sizeof(double));
  double *empty;
  empty = (double *)calloc(N_BINS-1,sizeof(double));
  double *E;
  E = (double *)calloc(N_BINS-1,sizeof(double));
  for(i=0; i<N_BINS-1; i++){
    E[i]=e[i];
  }
  TH1D *integrand4 = new TH1D("integrand4","integrand4",N_BINS,E_MIN,E_MAX);
  TH1D *integrand8 = new TH1D("integrand8","integrand8",N_BINS,E_MIN,E_MAX);
  printf("Degeneracy:\n");
  for(i=1; i<=N_BINS-1; i++){
    integrand4->SetBinContent(i,N_up_4->GetBinContent(i)/N_down_4->GetBinContent(i+1));
    integrand4->SetBinError(i,integrand4->GetBinContent(i)*(N_up_4->GetBinError(i)/N_up_4->GetBinContent(i)+N_down_4->GetBinError(i+1)/N_down_4->GetBinContent(i+1)));
    integrand8->SetBinContent(i,N_up_8->GetBinContent(i)/N_down_8->GetBinContent(i+1));
    integrand8->SetBinError(i,integrand8->GetBinContent(i)*(N_up_8->GetBinError(i)/N_up_8->GetBinContent(i)+N_down_8->GetBinError(i+1)/N_down_8->GetBinContent(i+1)));
    log_ge4[i-1] = (1/4.0)*integrand4->IntegralAndError(1,i,log_ge4_err[i-1],"width");
    log_ge8[i-1] = (1/8.0)*integrand8->IntegralAndError(1,i,log_ge8_err[i-1],"width");
    log_ge4_err[i-1] = log_ge4_err[i-1]/4.0;
    log_ge8_err[i-1] = log_ge8_err[i-1]/8.0;
    ge4[i-1] = exp(log_ge4[i-1]);
    ge8[i-1] = exp(log_ge8[i-1]);
    ge4_err[i-1] = exp(log_ge8_err[i-1]);
    ge8_err[i-1] = exp(log_ge8_err[i-1]);
    empty[i-1] = 0.0;
    printf("DeltaE = 4: logg(E) %lf +- %lf, Deltag(E)/g(E) %lf. DeltaE=8: logg(E) %lf +- %lf, Deltag(E)/g(E) %lf\n",log_ge4[i-1],log_ge4_err[i-1],ge4_err[i-1]/ge4[i-1],log_ge8[i-1],log_ge8_err[i-1],ge8_err[i-1]/ge8[i-1]);
  }

  //Plot log g(E)
  TCanvas *c4 = new TCanvas("c4","log[g(E)]");
  c4->cd();
  TGraphErrors *log_ge4_gr = new TGraphErrors(N_BINS-1,E,log_ge4,empty,log_ge4_err);
  log_ge4_gr->Draw("ALP");
  log_ge4_gr->SetMarkerColor(kBlue);
  log_ge4_gr->SetMarkerStyle(21);
  log_ge4_gr->SetLineColor(kBlue);
  log_ge4_gr->SetTitle("log[g(E)]");
  log_ge4_gr->GetXaxis()->SetTitle("E");
  log_ge4_gr->GetYaxis()->SetTitle("log[g(E)]");
  TGraphErrors *log_ge8_gr = new TGraphErrors(N_BINS-1,E,log_ge8,empty,log_ge8_err);
  log_ge8_gr->Draw("SAMELP");
  log_ge8_gr->SetMarkerColor(kRed);
  log_ge8_gr->SetMarkerStyle(21);
  log_ge8_gr->SetLineColor(kRed);
  log_ge8_gr->SetTitle("log[g(E)]");
  log_ge8_gr->GetXaxis()->SetTitle("E");
  log_ge8_gr->GetYaxis()->SetTitle("log[g(E)]");
  TLegend *leg4 = new TLegend(0.6,0.7,0.8,0.8);
  leg4->AddEntry(log_ge4_gr,"log[g(E)] (#DeltaE=4)","LP");
  leg4->AddEntry(log_ge8_gr,"log[g(E)] (#DeltaE=8)","LP");
  leg4->Draw("SAME");

  //Compute m(beta)
  double beta;
  int N_BETA=0;
  for(beta=BETA_MIN; beta<=BETA_MAX; beta = beta+0.01){
    N_BETA++;
  }
  printf("# beta points: %d\n",N_BETA);
  double *m_beta4;
  m_beta4 = (double *)calloc(N_BETA,sizeof(double));
  double *m_beta4_err;
  m_beta4_err = (double *)calloc(N_BETA,sizeof(double));
  double *m_beta8;
  m_beta8 = (double *)calloc(N_BETA,sizeof(double));
  double *m_beta8_err;
  m_beta8_err = (double *)calloc(N_BETA,sizeof(double));
  double *beta_val;
  beta_val = (double *)calloc(N_BETA,sizeof(double));
  double *empty_beta;
  empty_beta = (double *)calloc(N_BETA,sizeof(double));
  if(m_beta4==NULL || m_beta4_err==NULL || m_beta8==NULL || m_beta8_err==NULL || beta_val==NULL || empty_beta==NULL){
    printf("ERROR!\n");
    exit(EXIT_FAILURE);
  }
  double num_m4_err,den_4_err,num_m8_err,den_8_err;
  i=0;
  j=0;
  for(beta=BETA_MIN; beta<=BETA_MAX; beta = beta+0.01){
    TH1D *num_m4 = new TH1D("num_m4","num_m4",N_BINS,E_MIN,E_MAX);
    TH1D *den_4 = new TH1D("den_4","den_4",N_BINS,E_MIN,E_MAX);
    TH1D *num_m8 = new TH1D("num_m8","num_m8",N_BINS,E_MIN,E_MAX);
    TH1D *den_8 = new TH1D("den_8","den_8",N_BINS,E_MIN,E_MAX);
    for(j=1; j<=N_BINS-1; j++){
      num_m4->SetBinContent(j,m[j-1]*exp(log_ge4[j-1]-beta*e[j-1]));
      num_m4->SetBinError(j,sqrt(pow(exp(log_ge4[j-1]-beta*e[j-1])*m_err[j-1],2)+pow(m[j-1]*exp(log_ge4[j-1]-beta*e[j-1])*log_ge4_err[j-1],2)));
      den_4->SetBinContent(j,exp(log_ge4[j-1]-beta*e[j-1]));
      den_4->SetBinError(j,exp(log_ge4[j-1]-beta*e[j-1])*log_ge4_err[j-1]);
      num_m8->SetBinContent(j,m[j-1]*exp(log_ge8[j-1]-beta*e[j-1]));
      num_m8->SetBinError(j,sqrt(pow(exp(log_ge8[j-1]-beta*e[j-1])*m_err[j-1],2)+pow(m[j-1]*exp(log_ge8[j-1]-beta*e[j-1])*log_ge8_err[j-1],2)));
      den_8->SetBinContent(j,exp(log_ge8[j-1]-beta*e[j-1]));
      den_8->SetBinError(j,exp(log_ge8[j-1]-beta*e[j-1])*log_ge8_err[j-1]);
    }
    m_beta4[i] = num_m4->IntegralAndError(1,N_BINS,num_m4_err,"width")/den_4->IntegralAndError(1,N_BINS,den_4_err,"width");
    m_beta8[i] = num_m8->IntegralAndError(1,N_BINS,num_m8_err,"width")/den_8->IntegralAndError(1,N_BINS,den_8_err,"width");
    m_beta4_err[i] = m_beta4[i]*(num_m4_err/num_m4->Integral(1,N_BINS,"width")+den_4_err/den_4->Integral(1,N_BINS,"width"));
    m_beta8_err[i] = m_beta8[i]*(num_m8_err/num_m8->Integral(1,N_BINS,"width")+den_8_err/den_8->Integral(1,N_BINS,"width"));
    beta_val[i] = beta;
    printf("beta %lf, m4 %lf +- %lf, m8 %lf +- %lf\n",beta_val[i],m_beta4[i],m_beta4_err[i],m_beta8[i],m_beta8_err[i]); 
    empty_beta[i] = 0.0;
    num_m4->Delete();
    den_4->Delete();
    num_m8->Delete();
    den_8->Delete();
    i++;
  }

  //Plot m(beta)
  TCanvas *c5 = new TCanvas("c5","m(#beta)");
  c5->cd();
  TGraphErrors *m_beta4_gr = new TGraphErrors(N_BETA,beta_val,m_beta4,empty_beta,m_beta4_err);
  m_beta4_gr->Draw("ALP");
  m_beta4_gr->SetLineColor(kBlue);
  m_beta4_gr->SetMarkerStyle(21);
  m_beta4_gr->SetMarkerColor(kBlue);
  m_beta4_gr->GetXaxis()->SetTitle("#beta");
  m_beta4_gr->GetYaxis()->SetTitle("m(#beta)");
  m_beta4_gr->SetTitle("m(#beta)");
  TGraphErrors *m_beta8_gr = new TGraphErrors(N_BETA,beta_val,m_beta8,empty_beta,m_beta8_err);
  m_beta8_gr->Draw("SAMELP");
  m_beta8_gr->SetLineColor(kRed);
  m_beta8_gr->SetMarkerStyle(21);
  m_beta8_gr->SetMarkerColor(kRed);
  m_beta8_gr->GetXaxis()->SetTitle("#beta");
  m_beta8_gr->GetYaxis()->SetTitle("m(#beta)");
  m_beta8_gr->SetTitle("m(#beta)");
  TLegend *leg5 = new TLegend(0.6,0.7,0.8,0.8);
  leg5->AddEntry(m_beta4_gr,"m(#beta) (#DeltaE=4)","LP");
  leg5->AddEntry(m_beta8_gr,"m(#beta) (#DeltaE=8)","LP");
  leg5->Draw("SAME");

  gStyle->SetOptStat(0);
}
