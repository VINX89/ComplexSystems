#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define SEED 1057
#define N_STEPS 1000000
#define N_WALKERS 10000
#define PROB 0.5

main(){
  srand48(SEED);

  int t,i;
  long int x[N_WALKERS]={0};
  long double x2[N_WALKERS]={0};
  long double xmean,xdevst,x2mean,x2devst;
  double r;
  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk1d.txt","w");

  for(t=0; t<N_STEPS; t++){
    r=0;
    x2mean=0;
    for(i=0; i<N_WALKERS; i++){
      r = drand48();
      if(r<PROB){
	x[i] = x[i]-1;
      }
      else{
	x[i] = x[i]+1;
      }
      xmean = xmean + (long double)(x[i]);
      x2[i] = x[i]*x[i];
      x2mean = x2mean + x2[i];
    }
    xmean = xmean/N_WALKERS;
    x2mean = x2mean/N_WALKERS;
    xdevst=0;
    x2devst=0;
    for(i=0; i<N_WALKERS; i++){
      xdevst = xdevst + pow(x[i]-xmean,2);
      x2devst = x2devst + pow(x2[i]-x2mean,2);
    }
    xdevst = sqrt(xdevst/(N_WALKERS-1))/sqrt(N_WALKERS);
    x2devst = sqrt(x2devst/(N_WALKERS-1))/sqrt(N_WALKERS);
    fprintf(f,"%d %ld %Lf %Lf %Lf %Lf %Lf\n",t+1,x[0],xmean,xdevst,x2[0],x2mean,x2devst);
    if(t%100000==0){
      printf("Time %d, position %ld, mean pos. %Lf, sigma mean pos %Lf, position sq. %Lf, mean pos. sq. %Lf, sigma mean pos. sq. %Lf\n",t+1,x[0],xmean,xdevst,x2[0],x2mean,x2devst);
    }
  }
  
}
