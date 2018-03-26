#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define BETA 0.4399
#define L 150
#define SEED 1057
#define N_SWEEPS 500000

double boltz_prob(int deltaE, double beta);
int Delta_E(int x, int y, int site[][L], int est[], int ovest[], int nord[], int sud[]);
int Energy(int x, int y, int site[][L], int est[], int nord[]);
main(){
  srand48(SEED);

  printf("Beta=%lf, lattice size=%d\n",BETA,L);

  int i=0;
  int x=0;
  int y=0;
  int s[L][L];
  int sites = L*L;
  char data_path[100] = "/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Ising2d_B04399_L150.txt";
  
  //Initial conditions
  for(x=0; x<L; x++){
    for(y=0; y<L; y++){
      s[x][y]=1;
    }
  }

  //First neighbours
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
  
  //Onsager result
  double M=0;
  if(BETA>=0.4407){
    M = pow(1.0-pow((exp(2*BETA)-exp(-2*BETA))/2.0,-4.0),1.0/8.0);
  }
  printf("Onsager expected result: %lf\n",M);
  
  double prob[5]={1.0,0.0,boltz_prob(4,BETA),0.0,boltz_prob(8,BETA)};
  printf("Boltzmann prob: DeltaE=0: %lf, DeltaE=4: %lf, DeltaE=8: %lf\n",prob[0],prob[2],prob[4]);
  int delta = 0;
  int sum = 0;
  int en = 0;
  double m,m2,m4,En,En2;
  FILE *f = fopen(data_path,"w");

  //Montecarlo sweeps (Metropolis)
  for(i=0; i<N_SWEEPS; i++){
    sum=0;
    en=0;
    for(y=0; y<L; y++){
      for(x=0; x<L; x++){
	delta = Delta_E(x,y,s,E,O,N,S);
	if(delta<=0 || drand48()<=prob[delta/2]){
	  s[x][y]=-s[x][y];
	}
	sum = sum + s[x][y];
	en = en + Energy(x,y,s,E,N);
      }    
    }
    m = (double)(sum)/(double)(sites);
    m2 = m*m;
    m4 = m2*m2;
    En = (double)(en)/(double)(sites);
    En2 = En*En;
    fprintf(f,"%d %lf %d %lf %lf %lf %lf %lf\n",i+1,BETA,L,m,m2,m4,En,En2);
    if(i%100000==0){
      printf("Sweep num %d; m(%d)=%lf; E(%d)=%lf\n",i+1,i+1,m,i+1,En); 
    }
  }
}
double boltz_prob(int deltaE, double beta){
  return exp(-(double)(deltaE)*beta);
}
int Delta_E(int x, int y, int site[][L], int est[], int ovest[], int nord[], int sud[]){
  return 2*site[x][y]*(site[est[x]][y]+site[ovest[x]][y]+site[x][nord[y]]+site[x][sud[y]]);
}
int Energy(int x, int y, int site[][L], int est[], int nord[]){
  return -site[x][y]*(site[est[x]][y] + site[x][nord[y]]);
}

