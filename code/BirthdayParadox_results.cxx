#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAxis.h"
#include "TSystem.h"
#include "TLegend.h"

#define N 99

void BirthdayParadox_results(){
  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Segala_prob.txt","r");
  int k_temp=0;
  double ks[N]={0.0};
  double P_k_mc[N]={0.0};
  
  TCanvas *c1 = new TCanvas("Birthday Paradox");
  c1->cd();
  c1->SetTickx();
  c1->SetTicky();
  c1->SetGridx();
  c1->SetGridy();

  int i=0;
  while(fscanf(f,"%d %lf\n",&k_temp,&P_k_mc[i])!=EOF){
    ks[i] = (double)(k_temp);
    i++;
  }
  
  TGraph *mc = new TGraph(N,ks,P_k_mc);
  mc->Draw("AP");
  mc->SetMarkerColor(kBlue);
  mc->SetMarkerStyle(21);
  mc->SetTitle("Birthday Paradox");
  mc->GetXaxis()->SetTitle("# of people");
  mc->GetYaxis()->SetTitle("Probability [%]");

  double P_k_rec_temp[N]={0.0};
  P_k_rec_temp[0]=(1.0)/(365.0);
  double P_k_rec[N]={0.0};
  P_k_rec[0]=P_k_rec_temp[0]*100;
  double k = 2.0;
  printf("%lf people; probability=%lf \n",k,P_k_rec[0]);

  for(i=1; i<N; i++){
    P_k_rec_temp[i] = P_k_rec_temp[i-1]*(1-k/(365.0))+k/(365.0);
    P_k_rec[i] = P_k_rec_temp[i]*100;
    printf("%lf people; probability=%lf \n",k+1,P_k_rec[i]);
    k++;
  }

  TGraph *rec = new TGraph(N,ks,P_k_rec);
  rec->Draw("SAMEL");
  rec->SetLineColor(kRed);
  rec->SetLineWidth(3);
  
  TLegend *legend = new TLegend(0.6,0.7,0.8,0.8);
  legend->AddEntry(mc,"Montecarlo","P");
  legend->AddEntry(rec,"Recursion Formula","L");
  legend->Draw("SAME");

}
int main(){
  BirthdayParadox_results();
}
