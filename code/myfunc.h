#ifndef myfunc_h
#define myfunc_h

#define N_DATA_H 100000
#define N_DATA 500000
#define TMIN 10000
#define TMAX 500000 //Ising2d data end at t=500000                                                                                                     
#define TMIN_H 10000
#define TMAX_H 100000 //Heisenberg data end at t=100000 (L=10), t=50000 (L=30)
#define BMIN 100
#define BMAX 1000 
#define b_FINAL 1000
#define ABS 1 //take |m|, 0 take m (for Ising2d only)
#define WRITE 1 //1 write final results, 0 don't write final results

//Parameters for scaling laws 
#define BETA_C_50 0.4365       
#define BETA_C_100 0.439    
#define BETA_C_150 0.4395

#define BETA_C_10_H 0.6935
#define BETA_C_20_H 0.6935

#define beta_high_min -4.3
#define beta_high_max -1.5
#define alpha_low_min -5.8
#define alpha_low_max -0.5
#define alpha_high_min -5.8
#define alpha_high_max -0.5
#define gamma_low_min -4.7
#define gamma_low_max -0.5
#define gamma_high_min -5
#define gamma_high_max -0.5 

#define beta_high_minH -4.3
#define beta_high_maxH -1.0
#define alpha_low_minH -5.0
#define alpha_low_maxH -3.0
#define alpha_high_minH -5.0
#define alpha_high_maxH -3.0
#define gamma_low_minH -3.5
#define gamma_low_maxH 0.4
#define gamma_high_minH -2.9
#define gamma_high_maxH -0.9


FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Heisenberg_B0694_L20.txt","r");
FILE *out = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Heisenberg_results_L20.txt","a");

double mean(double data[], int tmin, int tmax){
  int i;
  int N=0;
  double m=0;
  for(i=tmin; i<tmax; i++){
    m = m + (double)(data[i]);
    N++;
  }
  return m/((double)(N));
}

double binning(double data[], int b, int tmin, int tmax){
  int j=0, k=0;
  int n_bins = (tmax-tmin)/b;
  double A=0;
  double *A_b;
  A_b = (double *)calloc(n_bins,sizeof(double));
  double mean=0.0;
  for(j=0; j<n_bins; j++){
    A=0;
    for(k=b*j+tmin; k<=b*(j+1)-1+tmin; k++){
      A = A + data[k];
    }
    A_b[j] = A/((double)(b));
    mean = mean + A_b[j];
  }
  mean = mean/((double)(n_bins));
  double variance=0.0;
  for(j=0; j<n_bins; j++){
    variance = variance + pow(A_b[j]-mean,2);
  }
  free(A_b);
  return sqrt((variance/((double)(n_bins*(n_bins-1)))));
}

long double jackknife_binder(double m2[], double m4[], int b, int tmin, int tmax){
  int j=0, k=0;
  int n_bins = (tmax-tmin)/b;
  double m2_tmp=0.0;
  double m4_tmp=0.0;
  double *m2_b;
  m2_b = (double *)calloc(n_bins,sizeof(double));
  double *m4_b;
  m4_b = (double *)calloc(n_bins,sizeof(double));
  for(j=0; j<n_bins; j++){
    m2_tmp=0.0;
    m4_tmp=0.0;
    for(k=b*j+tmin; k<=b*(j+1)-1+tmin; k++){
      m2_tmp = m2_tmp + m2[k];
      m4_tmp = m4_tmp + m4[k];
    }
    m2_b[j] = m2_tmp/((double)(b));
    m4_b[j] = m4_tmp/((double)(b));
  }
  double m2_mean = mean(m2_b,0,n_bins-1);
  double m4_mean = mean(m4_b,0,n_bins-1);
  double m2_j=0.0;
  double m4_j=0.0;
  long double bind=0.0;
  for(j=0; j<n_bins; j++){
    m2_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*m2_mean - m2_b[j]);
    m4_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*m4_mean - m4_b[j]);
    bind = bind + pow(0.5*(3.0 - m4_j/(m2_j*m2_j)) - 0.5*(3.0 - m4_mean/(m2_mean*m2_mean)),2);
  }
  free(m2_b);
  free(m4_b);
  return sqrt(((double)(n_bins - 1))/((double)(n_bins))*bind);
}

long double jackknife_cv(double e2[], double e[], double beta, int b, int tmin, int tmax){
  int j=0, k=0;
  int n_bins = (tmax-tmin)/b;
  double e2_tmp=0.0;
  double e_tmp=0.0;
  double *e2_b;
  e2_b = (double *)calloc(n_bins,sizeof(double));
  double *e_b;
  e_b = (double *)calloc(n_bins,sizeof(double));
  for(j=0; j<n_bins; j++){
    e2_tmp=0.0;
    e_tmp=0.0;
    for(k=b*j+tmin; k<=b*(j+1)-1+tmin; k++){
      e2_tmp = e2_tmp + e2[k];
      e_tmp = e_tmp + e[k];
    }
    e2_b[j] = e2_tmp/((double)(b));
    e_b[j] = e_tmp/((double)(b));
  }
  double e2_mean = mean(e2_b,0,n_bins-1);
  double e_mean = mean(e_b,0,n_bins-1);
  double e2_j=0.0;
  double e_j=0.0;
  long double err=0.0;
  for(j=0; j<n_bins; j++){
    e2_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*e2_mean - e2_b[j]);
    e_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*e_mean - e_b[j]);
    err = err + pow(beta*beta*((e2_j-e_j*e_j)-(e2_mean-e_mean*e_mean)),2);
  }
  free(e2_b);
  free(e_b);
  return sqrt(((double)(n_bins - 1))/((double)(n_bins))*err);
}

long double jackknife_chi(double m2[], double m[], double beta, int b, int tmin, int tmax){
  int j=0, k=0;
  int n_bins = (tmax-tmin)/b;
  double m2_tmp=0.0;
  double m_tmp=0.0;
  double *m2_b;
  m2_b = (double *)calloc(n_bins,sizeof(double));
  double *m_b;
  m_b = (double *)calloc(n_bins,sizeof(double));
  for(j=0; j<n_bins; j++){
    m2_tmp=0.0;
    m_tmp=0.0;
    for(k=b*j+tmin; k<=b*(j+1)-1+tmin; k++){
      m2_tmp = m2_tmp + m2[k];
      m_tmp = m_tmp + m[k];
    }
    m2_b[j] = m2_tmp/((double)(b));
    m_b[j] = m_tmp/((double)(b));
  }
  double m2_mean = mean(m2_b,0,n_bins-1);
  double m_mean = mean(m_b,0,n_bins-1);
  double m2_j=0.0;
  double m_j=0.0;
  long double err=0.0;
  for(j=0; j<n_bins; j++){
    m2_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*m2_mean - m2_b[j]);
    m_j = (1.0/(double)(n_bins-1))*((double)(n_bins)*m_mean - m_b[j]);
    err = err + pow(beta*((m2_j-m_j*m_j)-(m2_mean-m_mean*m_mean)),2);
  }
  free(m2_b);
  free(m_b);
  return sqrt(((double)(n_bins - 1))/((double)(n_bins))*err);
}

#endif
