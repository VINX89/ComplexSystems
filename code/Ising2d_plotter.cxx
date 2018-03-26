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
#include "TMath.h"

#include "myfunc.h"

void Ising2d_plotter(){
  char path_50[200]="/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Ising2d_results_L50.txt";
  char path_100[200]="/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Ising2d_results_L100.txt";
  char path_150[200]="/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Ising2d_results_L150.txt";

  //Retrieve L=50 data
  FILE *f_50 = fopen(path_50,"r");
  int i=0, j=0, k=0;
  double buffer=0.0;
  double beta=0.0;
  while(fscanf(f_50,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(beta<BETA_C_50){
      j++;
    }
    else if(beta>BETA_C_50){
      k++;
    }
    i++;
  }
  fclose(f_50);
  const int N_50 = i;
  const int scal_low_50 = j;
  const int scal_high_50 = k;

  double empty_50[N_50];
  double beta_50[N_50];
  int L=0;
  double m_50[N_50];
  double m_err_50[N_50];
  double e_50[N_50];
  double e_err_50[N_50];
  double Cv_50[N_50];
  double Cv_err_50[N_50];
  double X_50[N_50];
  double X_err_50[N_50];
  double B_50[N_50];
  double B_err_50[N_50];
  double t_low_50[scal_low_50];
  double t_low_err_50[scal_low_50];
  double t_high_50[scal_high_50];
  double t_high_err_50[scal_high_50];
  double m_high_50[scal_high_50];
  double m_high_err_50[scal_high_50];
  double Cv_low_50[scal_low_50];
  double Cv_low_err_50[scal_low_50];
  double Cv_high_50[scal_high_50];
  double Cv_high_err_50[scal_high_50];
  double X_low_50[scal_low_50];
  double X_low_err_50[scal_low_50];
  double X_high_50[scal_high_50];
  double X_high_err_50[scal_high_50];

  fopen(path_50,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_50,"%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta_50[i],&L,&m_50[i],&m_err_50[i],&e_50[i],&e_err_50[i],&Cv_50[i],&Cv_err_50[i],&X_50[i],&X_err_50[i],&B_50[i],&B_err_50[i])!=EOF){
    if(beta_50[i]<BETA_C_50){
      t_low_50[j]=TMath::Log(fabs((1/beta_50[i] - 1/BETA_C_50)/(1/BETA_C_50)));
      t_low_err_50[j]=0.0;
      Cv_low_50[j]=Cv_50[i];
      Cv_low_err_50[j]=Cv_err_50[i];
      X_low_50[j]=TMath::Log(X_50[i]);
      X_low_err_50[j]=(1.0/X_50[i])*X_err_50[i];
      j++;
    }
    else if(beta_50[i]>BETA_C_50){
      t_high_50[k]=TMath::Log(fabs((1/beta_50[i] - 1/BETA_C_50)/(1/BETA_C_50)));
      t_high_err_50[k]=0.0;
      m_high_50[k]=TMath::Log(m_50[i]);
      m_high_err_50[k]=fabs(1.0/(m_50[i])*m_err_50[i]);
      Cv_high_50[k]=Cv_50[i];
      Cv_high_err_50[k]=Cv_err_50[i];
      X_high_50[k]=TMath::Log(X_50[i]);
      X_high_err_50[k]=(1.0/X_50[i])*X_err_50[i];
      k++;
    }
    empty_50[i]=0.0;
    i++;
  }
  fclose(f_50);
  printf("L=%d data: %d\n",L,N_50);
  printf("Scaling laws; %d points below beta_c, %d points above beta_c\n",scal_low_50,scal_high_50);
  
  //Retrieve L=100 data
  FILE *f_100 = fopen(path_100,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_100,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(beta<BETA_C_100){
      j++;
    }
    else if(beta>BETA_C_100){
      k++;
    }
    i++;
  }
  fclose(f_100);
  const int N_100 = i;
  const int scal_low_100 = j;
  const int scal_high_100 = k;

  double empty_100[N_100];
  double beta_100[N_100];
  double m_100[N_100];
  double m_err_100[N_100];
  double e_100[N_100];
  double e_err_100[N_100];
  double Cv_100[N_100];
  double Cv_err_100[N_100];
  double X_100[N_100];
  double X_err_100[N_100];
  double B_100[N_100];
  double B_err_100[N_100];
  double t_low_100[scal_low_100];
  double t_low_err_100[scal_low_100];
  double t_high_100[scal_high_100];
  double t_high_err_100[scal_high_100];
  double m_high_100[scal_high_100];
  double m_high_err_100[scal_high_100];
  double Cv_low_100[scal_low_100];
  double Cv_low_err_100[scal_low_100];
  double Cv_high_100[scal_high_100];
  double Cv_high_err_100[scal_high_100];
  double X_low_100[scal_low_100];
  double X_low_err_100[scal_low_100];
  double X_high_100[scal_high_100];
  double X_high_err_100[scal_high_100];

  fopen(path_100,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_100,"%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta_100[i],&L,&m_100[i],&m_err_100[i],&e_100[i],&e_err_100[i],&Cv_100[i],&Cv_err_100[i],&X_100[i],&X_err_100[i],&B_100[i],&B_err_100[i])!=EOF){
    if(beta_100[i]<BETA_C_100){
      t_low_100[j]=TMath::Log(fabs((1/beta_100[i] - 1/BETA_C_100)/(1/BETA_C_100)));
      t_low_err_100[j]=0.0;
      Cv_low_100[j]=Cv_100[i];
      Cv_low_err_100[j]=Cv_err_100[i];
      X_low_100[j]=TMath::Log(X_100[i]);
      X_low_err_100[j]=(1.0/X_100[i])*X_err_100[i];
      j++;
    }
    else if(beta_100[i]>BETA_C_100){
      t_high_100[k]=TMath::Log(fabs((1/beta_100[i] - 1/BETA_C_100)/(1/BETA_C_100)));
      t_high_err_100[k]=0.0;
      m_high_100[k]=TMath::Log(m_100[i]);
      m_high_err_100[k]=fabs(1.0/(m_100[i])*m_err_100[i]);
      Cv_high_100[k]=Cv_100[i];
      Cv_high_err_100[k]=Cv_err_100[i];
      X_high_100[k]=TMath::Log(X_100[i]);
      X_high_err_100[k]=(1.0/X_100[i])*X_err_100[i];
      k++;
    }
    empty_100[i]=0.0;
    i++;
  }
  fclose(f_100);
  printf("L=%d data: %d\n",L,N_100);
  printf("Scaling laws; %d points below beta_c, %d points above beta_c\n",scal_low_100,scal_high_100);

  //Retrieve L=150 data 
  FILE *f_150 = fopen(path_150,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_150,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(beta<BETA_C_150){
      j++;
    }
    else if(beta>BETA_C_150){
      k++;
    }
    i++;
  }
  fclose(f_150);
  const int N_150 = i;
  const int scal_low_150 = j;
  const int scal_high_150 = k;

  double empty_150[N_150];
  double beta_150[N_150];
  double m_150[N_150];
  double m_err_150[N_150];
  double e_150[N_150];
  double e_err_150[N_150];
  double Cv_150[N_150];
  double Cv_err_150[N_150];
  double X_150[N_150];
  double X_err_150[N_150];
  double B_150[N_150];
  double B_err_150[N_150];
  double t_low_150[scal_low_150];
  double t_low_err_150[scal_low_150];
  double t_high_150[scal_high_150];
  double t_high_err_150[scal_high_150];
  double m_high_150[scal_high_150];
  double m_high_err_150[scal_high_150];
  double Cv_low_150[scal_low_150];
  double Cv_low_err_150[scal_low_150];
  double Cv_high_150[scal_high_150];
  double Cv_high_err_150[scal_high_150];
  double X_low_150[scal_low_150];
  double X_low_err_150[scal_low_150];
  double X_high_150[scal_high_150];
  double X_high_err_150[scal_high_150];

  fopen(path_150,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_150,"%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta_150[i],&L,&m_150[i],&m_err_150[i],&e_150[i],&e_err_150[i],&Cv_150[i],&Cv_err_150[i],&X_150[i],&X_err_150[i],&B_150[i],&B_err_150[i])!=EOF){
    if(beta_150[i]<BETA_C_150){
      t_low_150[j]=TMath::Log(fabs((1/beta_150[i] - 1/BETA_C_150)/(1/BETA_C_150)));
      t_low_err_150[j]=0.0;
      Cv_low_150[j]=Cv_150[i];
      Cv_low_err_150[j]=Cv_err_150[i];
      X_low_150[j]=TMath::Log(X_150[i]);
      X_low_err_150[j]=(1.0/X_150[i])*X_err_150[i];
      j++;
    }
    else if(beta_150[i]>BETA_C_150){
      t_high_150[k]=TMath::Log(fabs((1/beta_150[i] - 1/BETA_C_150)/(1/BETA_C_150)));
      t_high_err_150[k]=0.0;
      m_high_150[k]=TMath::Log(m_150[i]);
      m_high_err_150[k]=fabs(1.0/(m_150[i])*m_err_150[i]);
      Cv_high_150[k]=Cv_150[i];
      Cv_high_err_150[k]=Cv_err_150[i];
      X_high_150[k]=TMath::Log(X_150[i]);
      X_high_err_150[k]=(1.0/X_150[i])*X_err_150[i];
      k++;
    }
    empty_150[i]=0.0;
    i++;
  }
  fclose(f_150);
  printf("L=%d data: %d\n",L,N_150);
  printf("Scaling laws; %d points below beta_c, %d points above beta_c\n",scal_low_150,scal_high_150);

  //Plots
  TCanvas *c1 = new TCanvas("c1","<|m|> vs beta");
  c1->cd();
  TGraphErrors *m_gr_50 = new TGraphErrors(N_50,beta_50,m_50,empty_50,m_err_50);
  m_gr_50->Draw("ALP*");
  m_gr_50->GetXaxis()->SetTitle("#beta");
  m_gr_50->GetYaxis()->SetTitle("<|m|>");
  m_gr_50->SetLineColor(kBlue);
  m_gr_50->SetMarkerColor(kBlue);
  m_gr_50->SetTitle("<|m|> vs #beta");
  TGraphErrors *m_gr_100 = new TGraphErrors(N_100,beta_100,m_100,empty_100,m_err_100);
  m_gr_100->Draw("SAMELP*");
  m_gr_100->GetXaxis()->SetTitle("#beta");
  m_gr_100->GetYaxis()->SetTitle("<|m|>");
  m_gr_100->SetLineColor(kRed);
  m_gr_100->SetMarkerColor(kRed);
  m_gr_100->SetTitle("<|m|> vs #beta");
  TGraphErrors *m_gr_150 = new TGraphErrors(N_150,beta_150,m_150,empty_150,m_err_150);
  m_gr_150->Draw("SAMELP*");
  m_gr_150->GetXaxis()->SetTitle("#beta");
  m_gr_150->GetYaxis()->SetTitle("<|m|>");
  m_gr_150->SetLineColor(kBlack);
  m_gr_150->SetMarkerColor(kBlack);
  m_gr_150->SetTitle("<|m|> vs #beta");
  TLegend *legend1 = new TLegend(0.6,0.7,0.8,0.8);
  legend1->AddEntry(m_gr_50,"L=50","L");
  legend1->AddEntry(m_gr_100,"L=100","L");
  legend1->AddEntry(m_gr_150,"L=150","L");
  legend1->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2","<e> vs beta");
  c2->cd();
  TGraphErrors *e_gr_50 = new TGraphErrors(N_50,beta_50,e_50,empty_50,e_err_50);
  e_gr_50->Draw("ALP*");
  e_gr_50->GetXaxis()->SetTitle("#beta");
  e_gr_50->GetYaxis()->SetTitle("<e>");
  e_gr_50->SetLineColor(kBlue);
  e_gr_50->SetMarkerColor(kBlue);
  e_gr_50->SetTitle("<e> vs #beta");
  TGraphErrors *e_gr_100 = new TGraphErrors(N_100,beta_100,e_100,empty_100,e_err_100);
  e_gr_100->Draw("SAMELP*");
  e_gr_100->GetXaxis()->SetTitle("#beta");
  e_gr_100->GetYaxis()->SetTitle("<e>");
  e_gr_100->SetLineColor(kRed);
  e_gr_100->SetMarkerColor(kRed);
  e_gr_100->SetTitle("<e> vs #beta");
  TGraphErrors *e_gr_150 = new TGraphErrors(N_150,beta_150,e_150,empty_150,e_err_150);
  e_gr_150->Draw("SAMELP*");
  e_gr_150->GetXaxis()->SetTitle("#beta");
  e_gr_150->GetYaxis()->SetTitle("<e>");
  e_gr_150->SetLineColor(kBlack);
  e_gr_150->SetMarkerColor(kBlack);
  e_gr_150->SetTitle("<e> vs #beta");
  TLegend *legend2 = new TLegend(0.6,0.7,0.8,0.8);
  legend2->AddEntry(e_gr_50,"L=50","L");
  legend2->AddEntry(e_gr_100,"L=100","L");
  legend2->AddEntry(e_gr_150,"L=150","L");
  legend2->Draw("SAME");

  TCanvas *c3 =new TCanvas("c3","Cv vs beta");
  c3->cd();
  TGraphErrors *Cv_gr_50 = new TGraphErrors(N_50,beta_50,Cv_50,empty_50,Cv_err_50);
  Cv_gr_50->Draw("ALP*");
  Cv_gr_50->GetXaxis()->SetTitle("#beta");
  Cv_gr_50->GetYaxis()->SetTitle("C_{V}");
  Cv_gr_50->SetLineColor(kBlue);
  Cv_gr_50->SetMarkerColor(kBlue);
  Cv_gr_50->SetTitle("C_{V} vs #beta");
  TGraphErrors *Cv_gr_100 = new TGraphErrors(N_100,beta_100,Cv_100,empty_100,Cv_err_100);
  Cv_gr_100->Draw("SAMELP*");
  Cv_gr_100->GetXaxis()->SetTitle("#beta");
  Cv_gr_100->GetYaxis()->SetTitle("C_{V}");
  Cv_gr_100->SetLineColor(kRed);
  Cv_gr_100->SetMarkerColor(kRed);
  Cv_gr_100->SetTitle("C_{V} vs #beta");
  TGraphErrors *Cv_gr_150 = new TGraphErrors(N_150,beta_150,Cv_150,empty_150,Cv_err_150);
  Cv_gr_150->Draw("SAMELP*");
  Cv_gr_150->GetXaxis()->SetTitle("#beta");
  Cv_gr_150->GetYaxis()->SetTitle("C_{V}");
  Cv_gr_150->SetLineColor(kBlack);
  Cv_gr_150->SetMarkerColor(kBlack);
  Cv_gr_150->SetTitle("C_{V} vs #beta");
  TLegend *legend3 = new TLegend(0.6,0.7,0.8,0.8);
  legend3->AddEntry(Cv_gr_50,"L=50","L");
  legend3->AddEntry(Cv_gr_100,"L=100","L");
  legend3->AddEntry(Cv_gr_150,"L=150","L");
  legend3->Draw("SAME");

  TCanvas *c4 =new TCanvas("c4","X vs beta");
  c4->cd();
  TGraphErrors *X_gr_50 = new TGraphErrors(N_50,beta_50,X_50,empty_50,X_err_50);
  X_gr_50->Draw("ALP*");
  X_gr_50->GetXaxis()->SetTitle("#beta");
  X_gr_50->GetYaxis()->SetTitle("#chi");
  X_gr_50->SetLineColor(kBlue);
  X_gr_50->SetMarkerColor(kBlue);
  X_gr_50->SetTitle("#chi vs #beta");
  TGraphErrors *X_gr_100 = new TGraphErrors(N_100,beta_100,X_100,empty_100,X_err_100);
  X_gr_100->Draw("SAMELP*");
  X_gr_100->GetXaxis()->SetTitle("#beta");
  X_gr_100->GetYaxis()->SetTitle("#chi");
  X_gr_100->SetLineColor(kRed);
  X_gr_100->SetMarkerColor(kRed);
  X_gr_100->SetTitle("#chi vs #beta");
  TGraphErrors *X_gr_150 = new TGraphErrors(N_150,beta_150,X_150,empty_150,X_err_150);
  X_gr_150->Draw("SAMELP*");
  X_gr_150->GetXaxis()->SetTitle("#beta");
  X_gr_150->GetYaxis()->SetTitle("#chi");
  X_gr_150->SetLineColor(kBlack);
  X_gr_150->SetMarkerColor(kBlack);
  X_gr_150->SetTitle("#chi vs #beta");
  TLegend *legend4 = new TLegend(0.6,0.7,0.8,0.8);
  legend4->AddEntry(X_gr_50,"L=50","L");
  legend4->AddEntry(X_gr_100,"L=100","L");
  legend4->AddEntry(X_gr_150,"L=150","L");
  legend4->Draw("SAME");

  TCanvas *c5 =new TCanvas("c5","B vs beta");
  c5->cd();
  TGraphErrors *B_gr_50 = new TGraphErrors(N_50,beta_50,B_50,empty_50,B_err_50);
  B_gr_50->Draw("ALP*");
  B_gr_50->GetXaxis()->SetTitle("#beta");
  B_gr_50->GetYaxis()->SetTitle("B");
  B_gr_50->SetLineColor(kBlue);
  B_gr_50->SetMarkerColor(kBlue);
  B_gr_50->SetTitle("B vs #beta");
  TGraphErrors *B_gr_100 = new TGraphErrors(N_100,beta_100,B_100,empty_100,B_err_100);
  B_gr_100->Draw("SAMELP*");
  B_gr_100->GetXaxis()->SetTitle("#beta");
  B_gr_100->GetYaxis()->SetTitle("B");
  B_gr_100->SetLineColor(kRed);
  B_gr_100->SetMarkerColor(kRed);
  B_gr_100->SetTitle("B vs #beta");
  TGraphErrors *B_gr_150 = new TGraphErrors(N_150,beta_150,B_150,empty_150,B_err_150);
  B_gr_150->Draw("SAMELP*");
  B_gr_150->GetXaxis()->SetTitle("#beta");
  B_gr_150->GetYaxis()->SetTitle("B");
  B_gr_150->SetLineColor(kBlack);
  B_gr_150->SetMarkerColor(kBlack);
  B_gr_150->SetTitle("B vs #beta");
  TLegend *legend5 = new TLegend(0.6,0.7,0.8,0.8);
  legend5->AddEntry(B_gr_50,"L=50","L");
  legend5->AddEntry(B_gr_100,"L=100","L");
  legend5->AddEntry(B_gr_150,"L=150","L");
  legend5->Draw("SAME");

  //m scaling law
  TCanvas *c8 =new TCanvas("c8","m scaling law (L=150)");
  c8->cd();
  TGraphErrors *m_high_gr_150 = new TGraphErrors(scal_high_150,t_high_150,m_high_150,t_high_err_150,m_high_err_150);
  m_high_gr_150->Draw("AP");
  m_high_gr_150->GetXaxis()->SetTitle("log(-t)");
  m_high_gr_150->GetYaxis()->SetTitle("log(m)");
  m_high_gr_150->SetLineColor(kBlack);
  m_high_gr_150->SetMarkerColor(kBlack);
  m_high_gr_150->SetMarkerStyle(21);
  m_high_gr_150->SetTitle("m scaling law (L=150)");
  TF1 *m_high_fun_150 = new TF1("m_high_fun_150","[0]*x+[1]",-1000,1000);
  m_high_fun_150->SetParName(0,"#beta");
  m_high_gr_150->Fit("m_high_fun_150","","",beta_high_min,beta_high_max);

  //Cv scaling law
  TCanvas *c11 =new TCanvas("c11","Cv scaling law (L=150)");
  c11->cd();
  TGraphErrors *Cv_low_gr_150 = new TGraphErrors(scal_low_150,t_low_150,Cv_low_150,t_low_err_150,Cv_low_err_150);
  Cv_low_gr_150->Draw("AP");
  Cv_low_gr_150->GetXaxis()->SetTitle("log(|t|)");
  Cv_low_gr_150->GetYaxis()->SetTitle("C_{V}");
  Cv_low_gr_150->SetLineColor(kBlue);
  Cv_low_gr_150->SetMarkerColor(kBlue);
  Cv_low_gr_150->SetMarkerStyle(21);
  Cv_low_gr_150->SetTitle("C_{V} scaling law (L=150)");
  TF1 *Cv_low_fun_150 = new TF1("Cv_low_fun_150","[0]*x+[1]",-1000,1000);
  Cv_low_fun_150->SetLineColor(kBlue);
  Cv_low_gr_150->Fit("Cv_low_fun_150","","",alpha_low_min,alpha_low_max);
  TGraphErrors *Cv_high_gr_150 = new TGraphErrors(scal_high_150,t_high_150,Cv_high_150,t_high_err_150,Cv_high_err_150);
  Cv_high_gr_150->Draw("P*SAME");
  Cv_high_gr_150->GetXaxis()->SetTitle("log(|t|)");
  Cv_high_gr_150->GetYaxis()->SetTitle("C_{V}");
  Cv_high_gr_150->SetLineColor(kRed);
  Cv_high_gr_150->SetMarkerColor(kRed);
  Cv_high_gr_150->SetTitle("C_{V} scaling law (L=150)");
  TF1 *Cv_high_fun_150 = new TF1("Cv_high_fun_150","[0]*x+[1]",0.0,1.0);
  Cv_high_fun_150->SetLineColor(kRed);
  Cv_high_gr_150->Fit("Cv_high_fun_150","","",alpha_high_min,alpha_high_max);
  TLegend *legend11 = new TLegend(0.6,0.7,0.8,0.8);                                                                                                        
  legend11->AddEntry(Cv_low_gr_150,"T>T_{C}","LP");                                                                                                         
  legend11->AddEntry(Cv_high_gr_150,"T<T_{C}","LP");                                                                                                     
  legend11->Draw("SAME");

  TCanvas *c12 =new TCanvas("c12","X scaling law (L=150)");
  c12->cd();
  TGraphErrors *X_low_gr_150 = new TGraphErrors(scal_low_150,t_low_150,X_low_150,t_low_err_150,X_low_err_150);
  X_low_gr_150->Draw("AP");
  X_low_gr_150->GetXaxis()->SetTitle("log(|t|)");
  X_low_gr_150->GetYaxis()->SetTitle("log(#chi)");
  X_low_gr_150->SetLineColor(kBlue);
  X_low_gr_150->SetMarkerColor(kBlue);
  X_low_gr_150->SetMarkerStyle(21);
  X_low_gr_150->SetTitle("#chi scaling law (L=150)");
  TF1 *X_low_fun_150 = new TF1("X_low_fun_150","-[0]*x+[1]",-1000,1000);
  X_low_fun_150->SetLineColor(kBlue);
  X_low_fun_150->SetParName(0,"#gamma");
  X_low_gr_150->Fit("X_low_fun_150","","",gamma_low_min,gamma_low_max);
  TGraphErrors *X_high_gr_150 = new TGraphErrors(scal_high_150,t_high_150,X_high_150,t_high_err_150,X_high_err_150);
  X_high_gr_150->Draw("P*SAME");
  X_high_gr_150->GetXaxis()->SetTitle("log(|t|)");
  X_high_gr_150->GetYaxis()->SetTitle("log(#chi)");
  X_high_gr_150->SetLineColor(kRed);
  X_high_gr_150->SetMarkerColor(kRed);
  X_high_gr_150->SetTitle("#chi scaling law (L=150)");
  TF1 *X_high_fun_150 = new TF1("X_high_fun_150","-[0]*x+[1]",0.0,1.0);
  X_high_fun_150->SetLineColor(kRed);
  X_high_fun_150->SetParName(0,"#gamma");
  X_high_gr_150->Fit("X_high_fun_150","","",gamma_high_min,gamma_high_max);
  TLegend *legend12 = new TLegend(0.6,0.7,0.8,0.8);
  legend12->AddEntry(X_low_gr_150,"T>T_{C}","LP");
  legend12->AddEntry(X_high_gr_150,"T<T_{C}","LP");
  legend12->Draw("SAME");
}

