#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

#define L 8
#define SEED 1057
#define N_STEPS 1000000
#define E_MIN -2.0*L*L
#define E_MAX 0.0
#define BIN_SIZE 4.0
#define BETA_C 0.44

double boltz_prob(int deltaE, double beta);
int Delta_E(int x, int y, int site[][L], int est[], int ovest[], int nord[], int sud[]);
double Energy(int site[][L], int est[], int nord[], int dim);
double Mag(int site[][L], int dim);
void Broad_Hist(){
  srand48(SEED);
  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/test_L8.txt","w");

  int N_BINS=(int)(E_MAX-E_MIN-8.0)/(BIN_SIZE);
  printf("Energy axis: from %lf to %lf (%d bins)\n",E_MIN+8.0,E_MAX,N_BINS);
  printf("Bin size: %lf\n",BIN_SIZE);
  int check = N_STEPS/100;
  printf("Total number of steps: %d (check number each %d steps)\n",N_STEPS,check);

  int i=0, j=0;
  int x=0;
  int y=0;
  int s[L][L];
  
  //First neighbors                                                                                                                                   
  int E[L],O[L],N[L],S[L];
  E[0]=1;
  E[L-1]=0;
  O[0]=L-1;
  O[L-1]=L-2;
  N[0]=1;
  N[L-1]=0;
  S[0]=L-1;
  S[L-1]=L-2;
  for(i=1; i<L-1; i++){
    E[i] = i+1;
    O[i] = i-1;
    N[i] = i+1;
    S[i] = i-1;
  }

  //Initial conditions                                                                                                                                 
  for(x=0; x<L; x++){
    for(y=0; y<L; y++){
      s[x][y]=+1;
    }
  }
  
  //Metropolis thermalization to beta_c 
  double prob[5]={1.0,0.0,boltz_prob(4,BETA_C),0.0,boltz_prob(8,BETA_C)};
  int delta_e=0;
  printf("Metropolis thermalization to critical temperature\n");
  for(i=0; i<100000; i++){
    for(y=0; y<L; y++){
      for(x=0; x<L; x++){
        delta_e = Delta_E(x,y,s,E,O,N,S);
        if(delta_e<=0 || drand48()<=prob[delta_e/2]){
          s[x][y]=-s[x][y];
        }
      }
    }
  }
  printf("Initial energy = %lf (normalized: %lf)\n",Energy(s,E,N,L),Energy(s,E,N,L)/((double)(L*L)));
  printf("Initial magnetization = %lf\n",Mag(s,L));

  //Creating and initializing histograms
  TH1D *N_up_4 = new TH1D("N_up_4","N_{up}(e) (#DeltaE = 4)",N_BINS,E_MIN+8.0,E_MAX);
  TH1D *N_up_8 = new TH1D("N_up_8","N_{up}(e) (#DeltaE = 8)",N_BINS,E_MIN+8.0,E_MAX);
  TH1D *N_down_4 = new TH1D("N_down_4","N_{down}(e) (#DeltaE = 4)",N_BINS,E_MIN+8.0,E_MAX);
  TH1D *N_down_8 = new TH1D("N_down_8","N_{down}(e) (#DeltaE = 8)",N_BINS,E_MIN+8.0,E_MAX);
  TH1D *en = new TH1D("en","Energy distribution",N_BINS,E_MIN+8.0,E_MAX);
  double *m_e;
  m_e = (double *)calloc(N_BINS,sizeof(double));
  double *m_e2;
  m_e2 = (double *)calloc(N_BINS,sizeof(double));
  double *m_e_err;
  m_e_err = (double *)calloc(N_BINS,sizeof(double));
  double e;
  if(m_e==NULL || m_e2==NULL || m_e_err==NULL){
    printf("ERROR!\n");
    exit(EXIT_FAILURE);
  }
  for(e=E_MIN; e<=E_MAX; e = e + BIN_SIZE){
    N_up_4->Fill(e,1);
    N_up_8->Fill(e,1);
    N_down_4->Fill(e,1);
    N_down_8->Fill(e,1);
  }
    
  //E axis random walk
  int m_bin = 0;
  double beta = 0.0;
  int e_bin = 0, e_bin_plus_delta_e = 0;
  double N_up_e = 0.0, N_down_e_plus_delta_e = 0.0;
  delta_e = 0;
  e = 0;
  i=0;
  while(i<N_STEPS){
    e = Energy(s,E,N,L);
    if(e<=E_MAX && e>=E_MIN+8.0){
      en->Fill(e);
      m_bin = en->FindBin(e)-1;
      m_e[m_bin] = m_e[m_bin] + Mag(s,L);
      m_e2[m_bin] = m_e2[m_bin] + Mag(s,L)*Mag(s,L);
      //Statistics improvement
      for(j=0; j<1000; j++){
	for(y=0; y<L; y++){
	  for(x=0; x<L; x++){
	    delta_e = Delta_E(x,y,s,E,O,N,S);
	    if(delta_e == -4){
	      N_down_4->Fill(e);
	    }
	    else if(delta_e == -8){
	      N_down_8->Fill(e);
	    }
	    else if(delta_e == 4 && drand48()<= N_down_4->GetBinContent(N_down_4->FindBin(e))/N_up_4->GetBinContent(N_up_4->FindBin(e))){
	      N_up_4->Fill(e);
	    }
	    else if(delta_e == 8 && drand48()<= N_down_8->GetBinContent(N_down_8->FindBin(e))/N_up_8->GetBinContent(N_up_8->FindBin(e))){
	      N_up_8->Fill(e);
	    }
	    else{}
	  }
	}
      }
      x = drand48()*L;
      y = drand48()*L;
      delta_e = Delta_E(x,y,s,E,O,N,S);
      if(delta_e < 0 || (delta_e == 4 && drand48()<= N_down_4->GetBinContent(N_down_4->FindBin(e))/N_up_4->GetBinContent(N_up_4->FindBin(e))) 
	 || (delta_e == 8 && drand48()<= N_down_8->GetBinContent(N_down_8->FindBin(e))/N_up_8->GetBinContent(N_up_8->FindBin(e)))){
	s[x][y]=-s[x][y];
      }
      //Thermalization to current beta
      e_bin = en->FindBin(e);
      e_bin_plus_delta_e = en->FindBin(e+delta_e);
      beta=0.0;
      if(delta_e==4.0){
	N_up_e = N_up_4->GetBinContent(e_bin)/(N_up_4->GetBinContent(e_bin)+N_down_4->GetBinContent(e_bin));
	N_down_e_plus_delta_e = N_down_4->GetBinContent(e_bin_plus_delta_e)/(N_up_4->GetBinContent(e_bin_plus_delta_e)+N_down_4->GetBinContent(e_bin_plus_delta_e));
	beta = fabs((1.0/delta_e)*TMath::Log(N_up_e/N_down_e_plus_delta_e));
      }
      else if(delta_e==8.0){
	N_up_e = N_up_8->GetBinContent(e_bin)/(N_up_8->GetBinContent(e_bin)+N_down_8->GetBinContent(e_bin));
	N_down_e_plus_delta_e = N_down_8->GetBinContent(e_bin_plus_delta_e)/(N_up_8->GetBinContent(e_bin_plus_delta_e)+N_down_8->GetBinContent(e_bin_plus_delta_e));
	beta = fabs((1.0/delta_e)*TMath::Log(N_up_e/N_down_e_plus_delta_e));
      }
      else{}
      double prob_beta[5]={1.0,0.0,boltz_prob(4,beta),0.0,boltz_prob(8,beta)};
      for(j=0; j<10; j++){
	for(y=0; y<L; y++){
	  for(x=0; x<L; x++){
	    delta_e = Delta_E(x,y,s,E,O,N,S);
	    if(delta_e<=0 || drand48()<=prob_beta[delta_e/2]){
	      s[x][y]=-s[x][y];
	    }
	  }
	}
      }
      if(i%check==0){
	printf("Step number %d, e %lf, m %lf, beta %lf\n",i,e,Mag(s,L),beta);
      }
      i++;
    }
    else{
      //e outside [E_MIN,E_MAX] range: new starting configuration (where visits to E are lowest)
      x = drand48()*L;
      y = drand48()*L;
      double beta_new = 0.2 + drand48()*0.4;
      double prob_beta[5]={1.0,0.0,boltz_prob(4,beta_new),0.0,boltz_prob(8,beta_new)};
      for(j=0; j<10000; j++){
	for(y=0; y<L; y++){
	  for(x=0; x<L; x++){
	    int delta_e_new = Delta_E(x,y,s,E,O,N,S);
	    if(delta_e_new<=0 || drand48()<=prob_beta[delta_e_new/2]){
	      s[x][y]=-s[x][y];
	    }
	  }
	}
      }
    }
  }
  //Show results
  printf("Histogram contents:\n");
  for(i=1; i<=N_BINS; i++){
    printf("energy=%lf, visits=%lf, N_down_4=%lf, N_up_4=%lf, N_down_8=%lf, N_up_8=%lf\n",en->GetBinCenter(i),en->GetBinContent(i),N_down_4->GetBinContent(i),N_up_4->GetBinContent(i),N_down_8->GetBinContent(i),N_up_8->GetBinContent(i));
  }
  printf("Magnetization:\n");
  for(i=0; i<N_BINS; i++){
    m_e[i] = m_e[i]/(en->GetBinContent(i+1));
    m_e2[i] = m_e2[i]/(en->GetBinContent(i+1));
    m_e_err[i] = sqrt((m_e2[i] - m_e[i]*m_e[i])/(en->GetBinContent(i+1)-1));
    printf("m(e=%lf) = %lf +- %lf\n",en->GetBinCenter(i+1),m_e[i],m_e_err[i]);
  }
  //Store results
  for(i=1; i<=N_BINS; i++){
    fprintf(f,"%d %lf %lf %lf %lf %lf %lf %lf %lf\n",L,en->GetBinCenter(i),en->GetBinContent(i),m_e[i-1],m_e_err[i-1],N_down_4->GetBinContent(i),N_up_4->GetBinContent(i),N_down_8->GetBinContent(i),N_up_8->GetBinContent(i));
  }
  TCanvas *c1 = new TCanvas("c1","N_down_4 and N_up_4");
  c1->cd();
  N_down_4->Draw();
  N_down_4->SetLineColor(kBlue);
  N_down_4->SetTitle("#DeltaE=4");
  N_down_4->GetYaxis()->SetRangeUser(0.0,1.1*N_up_4->GetBinContent(N_up_4->GetMaximumBin()));
  N_up_4->Draw("SAME");
  N_up_4->SetLineColor(kRed);
  N_up_4->SetTitle("#DeltaE=4");
  TLegend *leg1 = new TLegend(0.6,0.7,0.8,0.8);
  leg1->AddEntry(N_down_4,"N_{down}(4)","L");
  leg1->AddEntry(N_up_4,"N_{up}(4)","L");
  leg1->Draw("SAME");
  TCanvas *c2 = new TCanvas("c2","N_down_8 and N_up_8");
  c2->cd();
  N_down_8->Draw();
  N_down_8->SetLineColor(kBlue);
  N_down_8->SetTitle("#DeltaE=8");
  N_down_8->GetYaxis()->SetRangeUser(0.0,1.1*N_up_8->GetBinContent(N_up_8->GetMaximumBin()));
  N_up_8->Draw("SAME");
  N_up_8->SetLineColor(kRed);
  N_up_8->SetTitle("#DeltaE=8");
  TLegend *leg2 = new TLegend(0.6,0.7,0.8,0.8);
  leg2->AddEntry(N_down_8,"N_{down}(8)","L");
  leg2->AddEntry(N_up_8,"N_{up}(8)","L");
  leg2->Draw("SAME");
  TCanvas *c3 = new TCanvas("c3","Energy distribution");
  c3->cd();
  en->Draw();

}
double boltz_prob(int deltaE, double beta){
  return exp(-(double)(deltaE)*beta);
}
int Delta_E(int x, int y, int site[][L], int est[], int ovest[], int nord[], int sud[]){
  return 2*site[x][y]*(site[est[x]][y]+site[ovest[x]][y]+site[x][nord[y]]+site[x][sud[y]]);
}
double Energy(int site[][L],int est[], int nord[], int dim){
  int X=0, Y=0;
  double E=0.0;
  for(Y=0; Y<dim; Y++){
    for(X=0; X<dim; X++){
      E = E - (double)(site[X][Y]*(site[est[X]][Y] + site[X][nord[Y]]));
    }
  }
  return E;
}
double Mag(int site[][L], int dim){
  int X=0, Y=0;
  double m=0;
  for(Y=0; Y<dim; Y++){
    for(X=0; X<dim; X++){
      m = m + site[X][Y];
    }
  }
  return fabs(m/((double)(dim*dim)));
}
