#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SEED 1057

main(){
  srand48(SEED);
  
  int i=0;
  int j=0;
  int k=0;
  int l=0;
  int yes=0;
  double match=0;
  double tot=1000000;
  int N=365;
  double P_k=0;

  FILE *f = fopen("/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Segala_prob.txt","w");
  
  for(k=2; k<=100; k++){
    printf("%d people in the room...\n",k);
    P_k=0;
    yes=0;
    match=0;
    int *day;
    day = (int*) calloc(k,sizeof(int));
    if(day==NULL){
      printf("ALLOCATION OF MEMORY FAILED\n");
      break;
    }
    for(j=0; j<tot; j++){
      for(i=0; i<k; i++){
	day[i] = 0;
	day[i] = (int)(drand48()*N);
	//printf("%d people, simulation %d, person %d, day %d\n",k,j,i,day[i]);
      }
      yes=0;
      for(i=0; i<k; i++){
	for(l=i+1; l<k; l++){
	  if(day[i]==day[l]){
	    yes++;
	  }
	}
      }
      if(yes!=0){
	match++;
      }
      //printf("%lf matches\n",match);
    }
    P_k = (match/tot)*100;
    printf("%d people: probability=%lf (%lf matches on %lf trials)\n",k,P_k,match,tot);
    fprintf(f,"%d %lf\n",k,P_k);
  }
  
}
