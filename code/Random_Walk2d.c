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
  long int x1[N_WALKERS]={0};
  long int x2[N_WALKERS]={0};
  double r1,r2;
  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk2d_traj.txt","w");
  FILE *g = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Random_Walk2d_prob.txt","w");

  for(t=1; t<=N_STEPS; t++){
    r1=0;
    r2=0;
    for(i=0; i<N_WALKERS; i++){
      r1 = drand48();
      r2 = drand48();
      if(r1<PROB){
        x1[i] = x1[i]-1;
      }
      else{
        x1[i] = x1[i]+1;
      }
      if(r2<PROB){
	x2[i] = x2[i]-1;
      }
      else{
	x2[i] = x2[i]+1;
      }
      if(t==10000 || t==100000 || t==N_STEPS){
	fprintf(g,"%d %ld %ld\n",t,x1[i],x2[i]);
      }
    }
    fprintf(f,"%d %ld %ld\n",t,x1[0],x2[0]);
    if(t%10000==0){
      printf("Time %d, x1=%ld, x2=%ld\n",t,x1[0],x2[0]);
    }
  }
}
