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

#include "myfunc.h"

void Ising2d_analysis(){
  double t[N_DATA];
  int L;
  double beta;
  double m[N_DATA]={0};
  double abs_m[N_DATA]={0};
  double m2[N_DATA]={0};
  double m4[N_DATA]={0};
  double e[N_DATA]={0};
  double e2[N_DATA]={0};

  TH1D *P_beta_m = new TH1D("P_beta_m","P_{#beta}(m)",200,-1.2,1.2);
  TH1D *P_beta_m_low = new TH1D("P_beta_m_low","P_{#beta}(m) (low)",120,-1.2,0.0);
  TH1D *P_beta_m_high = new TH1D("P_beta_m_high","P_{#beta}(m) (high)",120,0.0,1.2);

  int i=0;
  while(fscanf(f,"%lf %lf %d %lf %lf %lf %lf %lf\n",&t[i],&beta,&L,&m[i],&m2[i],&m4[i],&e[i],&e2[i])!=EOF){
    if(t[i]>TMIN){
      P_beta_m->Fill(m[i]);
      if(m[i]>=0){
	P_beta_m_high->Fill(m[i]);
      }
      else{
	P_beta_m_low->Fill(m[i]);
      }
    }
    abs_m[i] = fabs(m[i]);
    i++;
  }
  printf("Beta %lf, lattice size %d\n",beta,L);
  
  //Magnetization and energy plots
  TCanvas *c1 = new TCanvas("c1","m(t) and e(t)");
  c1->cd();
  c1->SetLogx();
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx();
  c1->SetGridy();
  TGraph *m_vs_t = new TGraph(N_DATA,t,m);
  m_vs_t->Draw("ALP");
  m_vs_t->SetTitle("m(t) and e(t)");
  m_vs_t->GetXaxis()->SetTitle("t");
  m_vs_t->SetLineColor(kBlue);
  m_vs_t->SetMarkerColor(kBlue);
  TGraph *e_vs_t = new TGraph(N_DATA,t,e);
  e_vs_t->Draw("LPSAME");
  e_vs_t->SetTitle("m(t) and e(t)");
  e_vs_t->GetXaxis()->SetTitle("t");
  e_vs_t->SetLineColor(kRed);
  e_vs_t->SetMarkerColor(kRed);
  m_vs_t->GetYaxis()->SetRangeUser(-2.0,2.0);
  TLegend *legend1 = new TLegend(0.6,0.7,0.8,0.8);
  legend1->AddEntry(m_vs_t,"m(t)","L");
  legend1->AddEntry(e_vs_t,"e(t)","L");
  legend1->Draw("SAME");

  //Magnetization probability density function
  TCanvas *c2 = new TCanvas("c2","P_beta(m)");
  c2->cd();
  P_beta_m->Draw();
  P_beta_m->Scale(1/(P_beta_m->Integral("width")));
  P_beta_m->SetLineColor(kBlue);
  P_beta_m->GetXaxis()->SetTitle("m");
  P_beta_m->GetYaxis()->SetTitle("entries/(0.012)");

  //Mean values computation  
  double abs_m_mean = mean(abs_m,TMIN,TMAX);
  double e_mean = mean(e,TMIN,TMAX);
  double m2_mean = mean(m2,TMIN,TMAX);
  double m4_mean = mean(m4,TMIN,TMAX);
  double cv_mean = (mean(e2,TMIN,TMAX) - pow(e_mean,2))*pow(beta,2);
  double chi_mean = (mean(m2,TMIN,TMAX) - pow(abs_m_mean,2))*beta;

  //P_m(beta) fitting
  TF1 *gauss_low = new TF1("gauss_low","gaus(0)",-1.0,0.0);
  TF1 *gauss_high = new TF1("gauss_high","gaus(0)",0.0,1.0);
  gauss_low->SetParameter(0,P_beta_m_low->GetBinContent(P_beta_m_low->GetMaximumBin()));
  gauss_low->SetParameter(1,abs_m_mean);
  gauss_low->SetParameter(2,P_beta_m_low->GetRMS());
  gauss_high->SetParameter(0,P_beta_m_high->GetBinContent(P_beta_m_high->GetMaximumBin()));
  gauss_high->SetParameter(1,abs_m_mean);
  gauss_high->SetParameter(2,P_beta_m_high->GetRMS());
  gauss_low->SetParLimits(0,0.0,50.0);
  gauss_low->SetParLimits(1,-1.0,0.05);
  gauss_low->SetParLimits(2,0.0,0.5);
  gauss_high->SetParLimits(0,0.0,50.0);
  gauss_high->SetParLimits(1,-0.05,1.0);
  gauss_high->SetParLimits(2,0.0,0.5);
  P_beta_m_low->Fit("gauss_low","M");
  P_beta_m_high->Fit("gauss_high","M");
  TF1 *gauss = new TF1("gauss","gaus(0)+gaus(3)",-1.0,1.0);
  gauss->SetParameter(0,gauss_low->GetParameter(0));
  gauss->SetParameter(1,gauss_low->GetParameter(1));
  gauss->SetParameter(2,gauss_low->GetParameter(2));
  gauss->SetParameter(3,gauss_high->GetParameter(0));
  gauss->SetParameter(4,gauss_high->GetParameter(1));
  gauss->SetParameter(5,gauss_high->GetParameter(2));
  gauss->SetParLimits(0,0.0,50.0);
  gauss->SetParLimits(1,-1.0,0.0);
  gauss->SetParLimits(2,0.0,0.5);
  gauss->SetParLimits(3,0.0,50.0);
  gauss->SetParLimits(4,0.0,1.0);
  gauss->SetParLimits(5,0.0,0.5);
  P_beta_m->Fit("gauss","M");
  double mu = gauss->GetParameter(4);
  double mu_err = gauss->GetParError(4);

  if(mu+mu_err>=0.0 && mu-mu_err<=0.0){
    printf("Gaussian means are compatible with m=0 -> NOT magnetized system\n");
  }
  else if(mu+mu_err<0.0 || mu-mu_err>0.0){
    printf("Gaussian means are NOT compatible with m=0 -> magnetized system\n");
  }
  if(ABS==0){
    printf("m taken\n");
    abs_m_mean = mean(m,TMIN,TMAX);
  }
  else if(ABS==1){
    printf("|m| taken\n");
  }

  //Abs m plot
  if(ABS==1){
    TCanvas *c3 = new TCanvas("c3","|m(t)|");
    c3->cd();
    c3->SetLogx();
    c3->SetTickx();
    c3->SetTicky();
    c3->SetGridx();
    c3->SetGridy();
    TGraph *absm_vs_t = new TGraph(N_DATA,t,abs_m);
    absm_vs_t->Draw("ALP");
    absm_vs_t->SetTitle("|m(t)|");
    absm_vs_t->GetXaxis()->SetTitle("t");
    absm_vs_t->SetLineColor(kBlue);
    absm_vs_t->SetMarkerColor(kBlue);
  }
  else{
    for(i=0; i<N_DATA; i++){
      abs_m[i] = m[i];
    }
  }

  //Binning counting
  int b=1;
  int n_temp=0;
  int n_meas=TMAX-TMIN;
  while(n_meas>=5){
    n_meas = (TMAX-TMIN)/b;
    n_temp++;
    b = b*2;
  }
  const int n = n_temp;
  double abs_m_err[n];
  double e_err[n];
  double m2_err[n];
  double m4_err[n];
  double b_values[n];
  
  //Error computation with binning method
  b=1;
  int count_min=0, count_max=0;
  int bmin, bmax;
  for(i=0; i<n; i++){
    abs_m_err[i] = binning(abs_m,b,TMIN,TMAX);
    e_err[i] = binning(e,b,TMIN,TMAX);
    m2_err[i] = binning(m2,b,TMIN,TMAX);
    m4_err[i] = binning(m4,b,TMIN,TMAX);
    b_values[i]=b;
    if(b >= BMIN && count_min==0){
      bmin=i;
      count_min++;
    }
    else if(b >= BMAX && count_max==0){
      bmax=i;
      count_max++;
    }
    b = b*2;
  }
  printf("Binning method, size of bins: min %d, max %d\n",bmin,bmax);
  double abs_m_error = mean(abs_m_err,bmin,bmax);
  double e_error = mean(e_err,bmin,bmax);
  double m2_error = mean(m2_err,bmin,bmax);
  double m4_error = mean(m4_err,bmin,bmax);
  printf("<|m|> = %lf +- %lf \n<m2> = %lf +- %lf \n<m4> = %lf +- %lf \n<e> = %lf +- %lf \n",abs_m_mean,abs_m_error,m2_mean,m2_error,m4_mean,m4_error,e_mean,e_error);
  TCanvas *c4 = new TCanvas("c4","Binning");
  c4->cd();
  c4->SetLogx();
  c4->SetLogy();
  c4->SetTickx();
  c4->SetTicky();
  c4->SetGridx();
  c4->SetGridy();
  TGraph *err_abs_m = new TGraph(n,b_values,abs_m_err);
  err_abs_m->Draw("ALP");
  err_abs_m->GetXaxis()->SetTitle("b");
  err_abs_m->SetTitle("<|m|> and <e> errors (binning method)");
  err_abs_m->SetLineColor(kBlue);
  err_abs_m->SetMarkerColor(kBlue);
  err_abs_m->SetMarkerStyle(21);
  TGraph *err_e = new TGraph(n,b_values,e_err);
  err_e->Draw("LPSAME");
  err_e->GetXaxis()->SetTitle("b");
  err_e->SetTitle("<|m|> and <e> errors (binning method)");
  err_e->SetLineColor(kRed);
  err_e->SetMarkerColor(kRed);
  err_e->SetMarkerStyle(21);
  TLegend *legend4 = new TLegend(0.6,0.7,0.8,0.8);
  legend4->AddEntry(err_abs_m,"#sigma_{|m|}","L");
  legend4->AddEntry(err_e,"#sigma_{e}","L");
  legend4->Draw("SAME");
  
  //Canvas edges tuning
  double y_min=abs_m_err[0], y_max = abs_m_err[0];
  for(i=0; i<n; i++){
    if(y_min>abs_m_err[i]){
      y_min=abs_m_err[i];
    }
    if(y_max<abs_m_err[i]){
      y_max=abs_m_err[i];
    }
  }
  for(i=0; i<n;i++){
    if(y_min>e_err[i]){
      y_min=e_err[i];
    }
    if(y_max<e_err[i]){
      y_max=e_err[i];
    }
  }
  err_abs_m->GetYaxis()->SetRangeUser(0.9*y_min,1.1*y_max);

  //Errors with jackknife method
  double binder = 0.5*(3 - mean(m4,TMIN,TMAX)/pow(mean(m2,TMIN,TMAX),2));
  long double binder_err = jackknife_binder(m2,m4,b_FINAL,TMIN,TMAX);
  long double cv_err = jackknife_cv(e2,e,beta,b_FINAL,TMIN,TMAX);
  long double chi_err = jackknife_chi(m2,abs_m,beta,b_FINAL,TMIN,TMAX);
  printf("Cv = %lf +- %Lf \nX = %lf +- %Lf \nB = %lf +- %Lf\n",cv_mean,cv_err,chi_mean,chi_err,binder,binder_err);

  //Results storage
  if(WRITE==1){
    fprintf(out,"%lf %d %lf %lf %lf %lf %lf %Lf %lf %Lf %lf %Lf\n",beta,L,abs_m_mean,abs_m_error,e_mean,e_error,cv_mean,cv_err,chi_mean,chi_err,binder,binder_err);
    printf("Results stored\n");
  }
  else{
    printf("Results not stored\n");
  }
}
