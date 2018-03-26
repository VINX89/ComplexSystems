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

void Heisenberg_plotter(){
  char path_10[200]="/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Heisenberg_results_L10.txt";
  char path_20[200]="/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Heisenberg_results_L20.txt";

  //Retrieve L=10 data                                                                                                                              
  FILE *f_10 = fopen(path_10,"r");
  int i=0, j=0, k=0;
  double buffer=0.0;
  double beta=0.0;
  while(fscanf(f_10,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(beta<BETA_C_10_H){
      j++;
    }
    else if(beta>BETA_C_10_H){
      k++;
    }
    i++;
  }
  fclose(f_10);
  const int N_10 = i;
  const int scal_low_10 = j;
  const int scal_high_10 = k;

  double empty_10[N_10];
  double beta_10[N_10];
  int L=0;
  double m_10[N_10];
  double m_err_10[N_10];
  double e_10[N_10];
  double e_err_10[N_10];
  double Cv_10[N_10];
  double Cv_err_10[N_10];
  double X_10[N_10];
  double X_err_10[N_10];
  double B_10[N_10];
  double B_err_10[N_10];
  double t_low_10[scal_low_10];
  double t_low_err_10[scal_low_10];
  double t_high_10[scal_high_10];
  double t_high_err_10[scal_high_10];
  double m_high_10[scal_high_10];
  double m_high_err_10[scal_high_10];
  double Cv_low_10[scal_low_10];
  double Cv_low_err_10[scal_low_10];
  double Cv_high_10[scal_high_10];
  double Cv_high_err_10[scal_high_10];
  double X_low_10[scal_low_10];
  double X_low_err_10[scal_low_10];
  double X_high_10[scal_high_10];
  double X_high_err_10[scal_high_10];

  fopen(path_10,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_10,"%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta_10[i],&L,&m_10[i],&m_err_10[i],&e_10[i],&e_err_10[i],&Cv_10[i],&Cv_err_10[i],&X_10[i],&X_err_10[i],&B_10[i],&B_err_10[i])!=EOF){
    if(beta_10[i]<BETA_C_10_H){
      t_low_10[j]=TMath::Log(fabs((1/beta_10[i] - 1/BETA_C_10_H)/(1/BETA_C_10_H)));
      t_low_err_10[j]=0.0;
      Cv_low_10[j]=TMath::Log(Cv_10[i]);
      Cv_low_err_10[j]=(1.0/Cv_10[i])*Cv_err_10[i];
      X_low_10[j]=TMath::Log(X_10[i]);
      X_low_err_10[j]=(1.0/X_10[i])*X_err_10[i];
      j++;
    }
    else if(beta_10[i]>BETA_C_10_H){
      t_high_10[k]=TMath::Log(fabs((1/beta_10[i] - 1/BETA_C_10_H)/(1/BETA_C_10_H)));
      t_high_err_10[k]=0.0;
      m_high_10[k]=TMath::Log(m_10[i]);
      m_high_err_10[k]=fabs(1.0/(m_10[i])*m_err_10[i]);
      Cv_high_10[k]=TMath::Log(Cv_10[i]);
      Cv_high_err_10[k]=(1.0/Cv_10[i])*Cv_err_10[i];
      X_high_10[k]=TMath::Log(X_10[i]);
      X_high_err_10[k]=(1.0/X_10[i])*X_err_10[i];
      k++;
    }
    empty_10[i]=0.0;
    i++;
  }
  fclose(f_10);
  printf("L=%d data: %d\n",L,N_10);
  printf("Scaling laws; %d points below beta_c, %d points above beta_c\n",scal_low_10,scal_high_10);

  
  //Retrieve L=20 data                                                                                                                                 
  FILE *f_20 = fopen(path_20,"r");
  i=0, j=0, k=0;
  buffer=0.0;
  beta=0.0;
  while(fscanf(f_20,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer,&buffer)!=EOF){
    if(beta<BETA_C_20_H){
      j++;
    }
    else if(beta>BETA_C_20_H){
      k++;
    }
    i++;
  }
  fclose(f_20);
  const int N_20 = i;
  const int scal_low_20 = j;
  const int scal_high_20 = k;

  double empty_20[N_20];
  double beta_20[N_20];
  double m_20[N_20];
  double m_err_20[N_20];
  double e_20[N_20];
  double e_err_20[N_20];
  double Cv_20[N_20];
  double Cv_err_20[N_20];
  double X_20[N_20];
  double X_err_20[N_20];
  double B_20[N_20];
  double B_err_20[N_20];
  double t_low_20[scal_low_20];
  double t_low_err_20[scal_low_20];
  double t_high_20[scal_high_20];
  double t_high_err_20[scal_high_20];
  double m_high_20[scal_high_20];
  double m_high_err_20[scal_high_20];
  double Cv_low_20[scal_low_20];
  double Cv_low_err_20[scal_low_20];
  double Cv_high_20[scal_high_20];
  double Cv_high_err_20[scal_high_20];
  double X_low_20[scal_low_20];
  double X_low_err_20[scal_low_20];
  double X_high_20[scal_high_20];
  double X_high_err_20[scal_high_20];

  fopen(path_20,"r");
  i=0;
  j=0;
  k=0;
  while(fscanf(f_20,"%lf %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&beta_20[i],&L,&m_20[i],&m_err_20[i],&e_20[i],&e_err_20[i],&Cv_20[i],&Cv_err_20[i],&X_20[i],&X_err_20[i],&B_20[i],&B_err_20[i])!=EOF){
    if(beta_20[i]<BETA_C_20_H){
      t_low_20[j]=TMath::Log(fabs((1/beta_20[i] - 1/BETA_C_20_H)/(1/BETA_C_20_H)));
      t_low_err_20[j]=0.0;
      Cv_low_20[j]=TMath::Log(Cv_20[i]);
      Cv_low_err_20[j]=(1.0/Cv_20[i])*Cv_err_20[i];
      X_low_20[j]=TMath::Log(X_20[i]);
      X_low_err_20[j]=(1.0/X_20[i])*X_err_20[i];
      j++;
    }
    else if(beta_20[i]>BETA_C_20_H){
      t_high_20[k]=TMath::Log(fabs((1/beta_20[i] - 1/BETA_C_20_H)/(1/BETA_C_20_H)));
      t_high_err_20[k]=0.0;
      m_high_20[k]=TMath::Log(m_20[i]);
      m_high_err_20[k]=fabs(1.0/(m_20[i])*m_err_20[i]);
      Cv_high_20[k]=TMath::Log(Cv_20[i]);
      Cv_high_err_20[k]=(1.0/Cv_20[i])*Cv_err_20[i];
      X_high_20[k]=TMath::Log(X_20[i]);
      X_high_err_20[k]=(1.0/X_20[i])*X_err_20[i];
      k++;
    }
    empty_20[i]=0.0;
    i++;
  }
  fclose(f_20);
  printf("L=%d data: %d\n",L,N_20);
  printf("Scaling laws; %d points below beta_c, %d points above beta_c\n",scal_low_20,scal_high_20);
  
  //Plots                                                                                                                                  
  TCanvas *c1 = new TCanvas("c1","<|m|> vs beta");
  c1->cd();
  TGraphErrors *m_gr_10 = new TGraphErrors(N_10,beta_10,m_10,empty_10,m_err_10);
  m_gr_10->Draw("ALP*");
  m_gr_10->GetXaxis()->SetTitle("#beta");
  m_gr_10->GetYaxis()->SetTitle("<|m|>");
  m_gr_10->SetLineColor(kRed);
  m_gr_10->SetMarkerColor(kRed);
  m_gr_10->SetTitle("<|m|> vs #beta");
  TGraphErrors *m_gr_20 = new TGraphErrors(N_20,beta_20,m_20,empty_20,m_err_20);
  m_gr_20->Draw("SAMELP*");
  m_gr_20->GetXaxis()->SetTitle("#beta");
  m_gr_20->GetYaxis()->SetTitle("<|m|>");
  m_gr_20->SetLineColor(kBlue);
  m_gr_20->SetMarkerColor(kBlue);
  m_gr_20->SetTitle("<|m|> vs #beta");
  TLegend *legend1 = new TLegend(0.6,0.7,0.8,0.8);
  legend1->AddEntry(m_gr_10,"L=10","L");
  legend1->AddEntry(m_gr_20,"L=20","L");
  legend1->Draw("SAME");

  TCanvas *c2 = new TCanvas("c2","<e> vs beta");
  c2->cd();
  TGraphErrors *e_gr_10 = new TGraphErrors(N_10,beta_10,e_10,empty_10,e_err_10);
  e_gr_10->Draw("ALP*");
  e_gr_10->GetXaxis()->SetTitle("#beta");
  e_gr_10->GetYaxis()->SetTitle("<e>");
  e_gr_10->SetLineColor(kRed);
  e_gr_10->SetMarkerColor(kRed);
  e_gr_10->SetTitle("<e> vs #beta");
  TGraphErrors *e_gr_20 = new TGraphErrors(N_20,beta_20,e_20,empty_20,e_err_20);
  e_gr_20->Draw("SAMELP*");
  e_gr_20->GetXaxis()->SetTitle("#beta");
  e_gr_20->GetYaxis()->SetTitle("<e>");
  e_gr_20->SetLineColor(kBlue);
  e_gr_20->SetMarkerColor(kBlue);
  e_gr_20->SetTitle("<e> vs #beta");
  TLegend *legend2 = new TLegend(0.6,0.7,0.8,0.8);
  legend2->AddEntry(e_gr_10,"L=10","L");
  legend2->AddEntry(e_gr_20,"L=20","L");
  legend2->Draw("SAME");

  TCanvas *c3 =new TCanvas("c3","Cv vs beta");
  c3->cd();
  TGraphErrors *Cv_gr_10 = new TGraphErrors(N_10,beta_10,Cv_10,empty_10,Cv_err_10);
  Cv_gr_10->Draw("ALP*");
  Cv_gr_10->GetXaxis()->SetTitle("#beta");
  Cv_gr_10->GetYaxis()->SetTitle("C_{V}");
  Cv_gr_10->SetLineColor(kRed);
  Cv_gr_10->SetMarkerColor(kRed);
  Cv_gr_10->SetTitle("C_{V} vs #beta");
  TGraphErrors *Cv_gr_20 = new TGraphErrors(N_20,beta_20,Cv_20,empty_20,Cv_err_20);
  Cv_gr_20->Draw("SAMELP*");
  Cv_gr_20->GetXaxis()->SetTitle("#beta");
  Cv_gr_20->GetYaxis()->SetTitle("C_{V}");
  Cv_gr_20->SetLineColor(kBlue);
  Cv_gr_20->SetMarkerColor(kBlue);
  Cv_gr_20->SetTitle("C_{V} vs #beta");
  TLegend *legend3 = new TLegend(0.6,0.7,0.8,0.8);
  legend3->AddEntry(Cv_gr_10,"L=10","L");
  legend3->AddEntry(Cv_gr_20,"L=20","L");
  legend3->Draw("SAME");

  TCanvas *c4 =new TCanvas("c4","X vs beta");
  c4->cd();
  TGraphErrors *X_gr_10 = new TGraphErrors(N_10,beta_10,X_10,empty_10,X_err_10);
  X_gr_10->Draw("ALP*");
  X_gr_10->GetXaxis()->SetTitle("#beta");
  X_gr_10->GetYaxis()->SetTitle("#chi");
  X_gr_10->SetLineColor(kRed);
  X_gr_10->SetMarkerColor(kRed);
  X_gr_10->SetTitle("#chi vs #beta");
  TGraphErrors *X_gr_20 = new TGraphErrors(N_20,beta_20,X_20,empty_20,X_err_20);
  X_gr_20->Draw("SAMELP*");
  X_gr_20->GetXaxis()->SetTitle("#beta");
  X_gr_20->GetYaxis()->SetTitle("#chi");
  X_gr_20->SetLineColor(kBlue);
  X_gr_20->SetMarkerColor(kBlue);
  X_gr_20->SetTitle("#chi vs #beta");
  TLegend *legend4 = new TLegend(0.6,0.7,0.8,0.8);
  legend4->AddEntry(X_gr_10,"L=10","L");
  legend4->AddEntry(X_gr_20,"L=20","L");
  legend4->Draw("SAME");

  TCanvas *c5 =new TCanvas("c5","B vs beta");
  c5->cd();
  TGraphErrors *B_gr_10 = new TGraphErrors(N_10,beta_10,B_10,empty_10,B_err_10);
  B_gr_10->Draw("ALP*");
  B_gr_10->GetXaxis()->SetTitle("#beta");
  B_gr_10->GetYaxis()->SetTitle("B");
  B_gr_10->SetLineColor(kRed);
  B_gr_10->SetMarkerColor(kRed);
  B_gr_10->SetTitle("B vs #beta");
  TGraphErrors *B_gr_20 = new TGraphErrors(N_20,beta_20,B_20,empty_20,B_err_20);
  B_gr_20->Draw("SAMELP*");
  B_gr_20->GetXaxis()->SetTitle("#beta");
  B_gr_20->GetYaxis()->SetTitle("B");
  B_gr_20->SetLineColor(kBlue);
  B_gr_20->SetMarkerColor(kBlue);
  B_gr_20->SetTitle("B vs #beta");
  TLegend *legend5 = new TLegend(0.6,0.7,0.8,0.8);
  legend5->AddEntry(B_gr_10,"L=10","L");
  legend5->AddEntry(B_gr_20,"L=20","L");
  legend5->Draw("SAME");

  //m scaling law                                                                                                                                    
  TCanvas *c6 =new TCanvas("c6","m scaling law (L=20)");
  c6->cd();
  TGraphErrors *m_high_gr_20 = new TGraphErrors(scal_high_20,t_high_20,m_high_20,t_high_err_20,m_high_err_20);
  m_high_gr_20->Draw("AP");
  m_high_gr_20->GetXaxis()->SetTitle("log(-t)");
  m_high_gr_20->GetYaxis()->SetTitle("log(m)");
  m_high_gr_20->SetLineColor(kBlue);
  m_high_gr_20->SetMarkerColor(kBlue);
  m_high_gr_20->SetMarkerStyle(21);
  m_high_gr_20->SetTitle("m scaling law (L=20)");
  TF1 *m_high_fun_20 = new TF1("m_high_fun_20","[0]*x+[1]",-1000,1000);
  m_high_fun_20->SetParName(0,"#beta");
  m_high_gr_20->Fit("m_high_fun_20","","",beta_high_minH,beta_high_maxH);

  //Cv scaling law                                                                                                                                     
  TCanvas *c7 =new TCanvas("c7","Cv scaling law (L=20)");
  c7->cd();
  TGraphErrors *Cv_low_gr_20 = new TGraphErrors(scal_low_20,t_low_20,Cv_low_20,t_low_err_20,Cv_low_err_20);
  Cv_low_gr_20->Draw("AP");
  Cv_low_gr_20->GetXaxis()->SetTitle("log(|t|)");
  Cv_low_gr_20->GetYaxis()->SetTitle("log(C_{V})");
  Cv_low_gr_20->SetLineColor(kBlue);
  Cv_low_gr_20->SetMarkerColor(kBlue);
  Cv_low_gr_20->SetMarkerStyle(21);
  Cv_low_gr_20->SetTitle("C_{V} scaling law (L=20)");
  TF1 *Cv_low_fun_20 = new TF1("Cv_low_fun_20","-[0]*x+[1]",-1000,1000);
  Cv_low_fun_20->SetParName(0,"#alpha");
  Cv_low_fun_20->SetLineColor(kBlue);
  Cv_low_gr_20->Fit("Cv_low_fun_20","","",alpha_low_minH,alpha_low_maxH);
  TGraphErrors *Cv_high_gr_20 = new TGraphErrors(scal_high_20,t_high_20,Cv_high_20,t_high_err_20,Cv_high_err_20);
  Cv_high_gr_20->Draw("P*SAME");
  Cv_high_gr_20->GetXaxis()->SetTitle("log(|t|)");
  Cv_high_gr_20->GetYaxis()->SetTitle("log(C_{V})");
  Cv_high_gr_20->SetLineColor(kRed);
  Cv_high_gr_20->SetMarkerColor(kRed);
  Cv_high_gr_20->SetTitle("C_{V} scaling law (L=20)");
  TF1 *Cv_high_fun_20 = new TF1("Cv_high_fun_20","-[0]*x+[1]",-1000,1000);
  Cv_high_fun_20->SetParName(0,"#alpha");
  Cv_high_fun_20->SetLineColor(kRed);
  Cv_high_gr_20->Fit("Cv_high_fun_20","","",alpha_high_minH,alpha_high_maxH);
  TLegend *legend6 = new TLegend(0.6,0.7,0.8,0.8);
  legend6->AddEntry(Cv_low_gr_20,"T>T_{C}","LP");
  legend6->AddEntry(Cv_high_gr_20,"T<T_{C}","LP");
  legend6->Draw("SAME");

  TCanvas *c8 =new TCanvas("c8","X scaling law (L=20)");
  c8->cd();
  TGraphErrors *X_low_gr_20 = new TGraphErrors(scal_low_20,t_low_20,X_low_20,t_low_err_20,X_low_err_20);
  X_low_gr_20->Draw("AP");
  X_low_gr_20->GetXaxis()->SetTitle("log(|t|)");
  X_low_gr_20->GetYaxis()->SetTitle("log(#chi)");
  X_low_gr_20->SetLineColor(kBlue);
  X_low_gr_20->SetMarkerColor(kBlue);
  X_low_gr_20->SetMarkerStyle(21);
  X_low_gr_20->SetTitle("#chi scaling law (L=20)");
  TF1 *X_low_fun_20 = new TF1("X_low_fun_20","-[0]*x+[1]",-1000,1000);
  X_low_fun_20->SetLineColor(kBlue);
  X_low_fun_20->SetParName(0,"#gamma");
  X_low_gr_20->Fit("X_low_fun_20","","",gamma_low_minH,gamma_low_maxH);
  TGraphErrors *X_high_gr_20 = new TGraphErrors(scal_high_20,t_high_20,X_high_20,t_high_err_20,X_high_err_20);
  X_high_gr_20->Draw("P*SAME");
  X_high_gr_20->GetXaxis()->SetTitle("log(|t|)");
  X_high_gr_20->GetYaxis()->SetTitle("log(#chi)");
  X_high_gr_20->SetLineColor(kRed);
  X_high_gr_20->SetMarkerColor(kRed);
  X_high_gr_20->SetTitle("#chi scaling law (L=20)");
  TF1 *X_high_fun_20 = new TF1("X_high_fun_20","-[0]*x+[1]",-1000,1000);
  X_high_fun_20->SetLineColor(kRed);
  X_high_fun_20->SetParName(0,"#gamma");
  X_high_gr_20->Fit("X_high_fun_20","","",gamma_high_minH,gamma_high_maxH);
  TLegend *legend7 = new TLegend(0.6,0.7,0.8,0.8);
  legend7->AddEntry(X_low_gr_20,"T>T_{C}","LP");
  legend7->AddEntry(X_high_gr_20,"T<T_{C}","LP");
  legend7->Draw("SAME");


}
