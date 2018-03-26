#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.1415
#define BETA 0.694
#define L 10
#define SEED 1057
#define N_SWEEPS 100000

double Random_zn(double r, double beta, double n);
double Rotate_xn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n);
double Rotate_yn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n);
double Rotate_zn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n);
main(){
  srand48(SEED);

  printf("Beta=%lf, lattice size=%d\n",BETA,L);

  int i=0;
  int x=0;
  int y=0;
  int z=0;
  double sx[L][L][L]={0};
  double sy[L][L][L]={0};
  double sz[L][L][L]={0};
  int sites = L*L*L;
  char data_path[100] = "/Users/Vincenzo/Desktop/University/Fisica_dei_Sistemi_Complessi/data/Heisenberg_B0694_L10.txt";

  //First neighbours                                                                                                                                   
  int E[L]={0},O[L]={0},N[L]={0},S[L]={0},U[L]={0},D[L]={0};
  E[0]=1;
  E[L-1]=0;
  O[0]=L-1;
  O[L-1]=L-2;
  N[0]=1;
  N[L-1]=0;
  S[0]=L-1;
  S[L-1]=L-2;
  U[0]=1;
  U[L-1]=0;
  D[0]=L-1;
  D[L-1]=L-2;
  for(i=1; i<L-1; i++){
    E[i] = i+1;
    O[i] = i-1;
    N[i] = i+1;
    S[i] = i-1;
    U[i] = i+1;
    D[i] = i-1;
  }

  //Initial conditions                                                                                                                                 
  double s_tot_x=0.0, s_tot_y=0.0, s_tot_z=0.0;
  double m,e=0.0,ex=0.0,ey=0.0,ez=0.0;
  for(x=0; x<L; x++){
    for(y=0; y<L; y++){
      for(z=0; z<L; z++){
	sz[x][y][z]= 1.0;
	sx[x][y][z] = 0.0;
	sy[x][y][z] = 0.0;
	s_tot_x = s_tot_x + sx[x][y][z];
	s_tot_y = s_tot_y + sy[x][y][z];
	s_tot_z = s_tot_z + sz[x][y][z];
	ex = ex - sx[x][y][z]*(sx[E[x]][y][z] + sx[x][N[y]][z] + sx[x][y][U[z]]);                                                                          
	ey = ey - sy[x][y][z]*(sy[E[x]][y][z] + sy[x][N[y]][z] + sy[x][y][U[z]]);                                                  
	ez = ez - sz[x][y][z]*(sz[E[x]][y][z] + sz[x][N[y]][z] + sz[x][y][U[z]]);                                                                   
	e = ex + ey + ez;                                               
      }
    }
  }
  m = sqrt(s_tot_x*s_tot_x + s_tot_y*s_tot_y + s_tot_z*s_tot_z)/((double)(sites));
  printf("Initial |m|: %lf\n",m);
  printf("Initial s components: sx %lf, sy %lf, sz %lf\n",s_tot_x,s_tot_y,s_tot_z);
  e = e/((double)(sites));
  printf("Initial e: %lf\n",e);

  //Montecarlo sweeps (heat bath)
  double n;
  double n_x, n_y, n_z;
  double x_n, y_n, z_n;
  double phi_n;
  double proj_xy_n;
  double mz,m2,m4,e2;
  FILE *f = fopen(data_path,"w");

  for(i=0; i<N_SWEEPS; i++){
    s_tot_x=0.0;
    s_tot_y=0.0;
    s_tot_z=0.0;
    ex=0.0;
    ey=0.0;
    ez=0.0;
    e=0.0;
    for(x=0; x<L; x++){
      for(y=0; y<L; y++){
	for(z=0; z<L; z++){
	  n_x = sx[E[x]][y][z]+sx[O[x]][y][z]+sx[x][N[y]][z]+sx[x][S[y]][z]+sx[x][y][U[z]]+sx[x][y][D[z]];  //
	  n_y = sy[E[x]][y][z]+sy[O[x]][y][z]+sy[x][N[y]][z]+sy[x][S[y]][z]+sy[x][y][U[z]]+sy[x][y][D[z]];  // Compute first neighbours spin sum in the  
	  n_z = sz[E[x]][y][z]+sz[O[x]][y][z]+sz[x][N[y]][z]+sz[x][S[y]][z]+sz[x][y][U[z]]+sz[x][y][D[z]];  // "absolute" frame
	  n = sqrt(n_x*n_x + n_y*n_y + n_z*n_z);                                                            // 
	  z_n = Random_zn(drand48(),BETA,n);  //
	  phi_n = drand48()*2*PI;             // Compute spin coordinates in n vector frame, in which n is aligned with z_n axis 
	  proj_xy_n = sqrt(1-z_n*z_n);        // 
	  x_n = proj_xy_n*cos(phi_n);         // 
	  y_n = proj_xy_n*sin(phi_n);         //
	  n_x = n_x/n; //
	  n_y = n_y/n; // Normalize n to one (n->k_n unit vector coordinates in absolute frame)
	  n_z = n_z/n; //
	  sx[x][y][z]=Rotate_xn(n_x,n_y,n_z,x_n,y_n,z_n);  //
	  sy[x][y][z]=Rotate_yn(n_x,n_y,n_z,x_n,y_n,z_n);  // Compute spin coordinates in absolute frame
	  sz[x][y][z]=Rotate_zn(n_x,n_y,n_z,x_n,y_n,z_n);  //
	  s_tot_x = s_tot_x + sx[x][y][z];  //
	  s_tot_y = s_tot_y + sy[x][y][z];  // Accumulate total spin coordinates sum
	  s_tot_z = s_tot_z + sz[x][y][z];  //
	  ex = ex - sx[x][y][z]*(sx[E[x]][y][z] + sx[x][N[y]][z] + sx[x][y][U[z]]);  //
	  ey = ey - sy[x][y][z]*(sy[E[x]][y][z] + sy[x][N[y]][z] + sy[x][y][U[z]]);  // Accumulate total energy sum
	  ez = ez - sz[x][y][z]*(sz[E[x]][y][z] + sz[x][N[y]][z] + sz[x][y][U[z]]);  //
	  e = ex + ey + ez;                                                          //
	}
      }
    }
    if(i%5000==0){
      printf("Over-relaxation...\n");
      n_x=0.0;
      n_y=0.0;
      n_z=0.0;
      for(x=0; x<L; x++){
	for(y=0; y<L; y++){
	  for(z=0; z<L; z++){
	    n_x = sx[E[x]][y][z]+sx[O[x]][y][z]+sx[x][N[y]][z]+sx[x][S[y]][z]+sx[x][y][U[z]]+sx[x][y][D[z]];
	    n_y = sy[E[x]][y][z]+sy[O[x]][y][z]+sy[x][N[y]][z]+sy[x][S[y]][z]+sy[x][y][U[z]]+sy[x][y][D[z]];
	    n_z = sz[E[x]][y][z]+sz[O[x]][y][z]+sz[x][N[y]][z]+sz[x][S[y]][z]+sz[x][y][U[z]]+sz[x][y][D[z]];
	    double n2 = n_x*n_x + n_y*n_y + n_z*n_z;
	    double n_dot_s = sx[x][y][z]*n_x + sy[x][y][z]*n_y + sz[x][y][z]*n_z;
	    sx[x][y][z] = (2*(n_dot_s)/(n2))*n_x - sx[x][y][z];
	    sy[x][y][z] = (2*(n_dot_s)/(n2))*n_y - sy[x][y][z];
	    sz[x][y][z] = (2*(n_dot_s)/(n2))*n_z - sz[x][y][z];
	  }
	}
      }
    }
    mz = s_tot_z/((double)(sites));
    m = sqrt(s_tot_x*s_tot_x + s_tot_y*s_tot_y + s_tot_z*s_tot_z)/((double)(sites));
    m2 = m*m;
    m4 = m2*m2;
    e = e/((double)(sites));
    e2 = e*e;
    fprintf(f,"%d %lf %d %lf %lf %lf %lf %lf\n",i+1,BETA,L,mz,m,m2,e,e2);
    if(i%1000==0){
      printf("---Sweep number %d:---\n",i+1);
      printf("|m|(%d)=%lf, m_z(%d)=%lf e(%d)=%lf\n",i+1,m,i+1,mz,i+1,e);
    }
  }
}
double Random_zn(double r, double beta, double n){
  double exp_plus = exp(beta*n);
  double exp_minus = exp(-beta*n);
  double N = 2*PI*(exp_plus - exp_minus);
  double min = (1/(N*beta*n))*exp_minus;
  double max = (1/(N*beta*n))*exp_plus;
  double rho = min + (max-min)*r;
  return (1/(beta*n))*log(N*beta*n*rho);
}
double Rotate_xn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n){
  double i_x_n,i_y_n,i_z_n,j_x_n,j_y_n,j_z_n,in_norm,jn_norm,k_dot_kn;
  if(k_z_n==1.0 && k_x_n==0.0 && k_y_n==0.0){
    i_x_n=1.0;
    i_y_n=0.0;
    i_z_n=0.0;
    j_x_n=0.0;
    j_y_n=1.0;
    j_z_n=0.0;
  }
  else if(k_z_n==-1.0 && k_x_n==0.0 && k_y_n==0.0){
    i_x_n=-1.0;
    i_y_n=0.0;
    i_z_n=0.0;
    j_x_n=0.0;
    j_y_n=-1.0;
    j_z_n=0.0;
  }
  else{
    k_dot_kn = k_z_n;
    i_x_n = -k_dot_kn*k_x_n;
    i_y_n = -k_dot_kn*k_y_n;
    i_z_n = 1 - k_dot_kn*k_z_n;
    in_norm = sqrt(i_x_n*i_x_n + i_y_n*i_y_n + i_z_n*i_z_n);
    i_x_n = i_x_n/in_norm;
    i_y_n = i_y_n/in_norm;
    i_z_n = i_z_n/in_norm;
    j_z_n=0.0;
    if(i_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -i_y_n/i_x_n;
    }
    else if(i_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -i_x_n/i_y_n;
    }
    else if(k_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -k_y_n/k_x_n;
    }
    else if(k_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -k_x_n/k_y_n;
    }
    else{}
    jn_norm = sqrt(j_x_n*j_x_n + j_y_n*j_y_n + j_z_n*j_z_n);
    j_x_n= j_x_n/jn_norm;
    j_y_n = j_y_n/jn_norm;
    j_z_n = j_z_n/jn_norm;
  }
  return i_x_n*x_n + j_x_n*y_n + k_x_n*z_n;
}
double Rotate_yn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n){
  double i_x_n,i_y_n,i_z_n,j_x_n,j_y_n,j_z_n,in_norm,jn_norm,k_dot_kn;
  if(k_z_n==1.0 && k_x_n==0.0 && k_y_n==0.0){
    i_x_n=1.0;
    i_y_n=0.0;
    i_z_n=0.0;
    j_x_n=0.0;
    j_y_n=1.0;
  }
  else{
    k_dot_kn = k_z_n;
    i_x_n = -k_dot_kn*k_x_n;
    i_y_n = -k_dot_kn*k_y_n;
    i_z_n = 1 - k_dot_kn*k_z_n;
    in_norm = sqrt(i_x_n*i_x_n + i_y_n*i_y_n + i_z_n*i_z_n);
    i_x_n = i_x_n/in_norm;
    i_y_n = i_y_n/in_norm;
    i_z_n = i_z_n/in_norm;
    j_z_n = 0.0;
    if(i_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -i_y_n/i_x_n;
    }
    else if(i_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -i_x_n/i_y_n;
    }
    else if(k_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -k_y_n/k_x_n;
    }
    else if(k_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -k_x_n/k_y_n;
    }
    else{}
    jn_norm = sqrt(j_x_n*j_x_n + j_y_n*j_y_n + j_z_n*j_z_n);
    j_x_n= j_x_n/jn_norm;
    j_y_n = j_y_n/jn_norm;
    j_z_n = j_z_n/jn_norm;
  }  
  return i_y_n*x_n + j_y_n*y_n + k_y_n*z_n;
}
double Rotate_zn(double k_x_n, double k_y_n, double k_z_n, double x_n, double y_n, double z_n){
  double i_x_n,i_y_n,i_z_n,j_x_n,j_y_n,j_z_n,in_norm,jn_norm,k_dot_kn;
  if(k_z_n==1.0 && k_x_n==0.0 && k_y_n==0.0){
    i_x_n=1.0;
    i_y_n=0.0;
    i_z_n=0.0;
    j_x_n=0.0;
    j_y_n=1.0;
  }
  else{
    k_dot_kn = k_z_n;
    i_x_n = -k_dot_kn*k_x_n;
    i_y_n = -k_dot_kn*k_y_n;
    i_z_n = 1 - k_dot_kn*k_z_n;
    in_norm = sqrt(i_x_n*i_x_n + i_y_n*i_y_n + i_z_n*i_z_n);
    i_x_n = i_x_n/in_norm;
    i_y_n = i_y_n/in_norm;
    i_z_n = i_z_n/in_norm;
    j_z_n = 0.0;
    if(i_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -i_y_n/i_x_n;
    }
    else if(i_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -i_x_n/i_y_n;
    }
    else if(k_x_n!=0){
      j_y_n = 1.0;
      j_x_n = -k_y_n/k_x_n;
    }
    else if(k_y_n!=0){
      j_x_n = 1.0;
      j_y_n = -k_x_n/k_y_n;
    }
    else{}
    jn_norm = sqrt(j_x_n*j_x_n + j_y_n*j_y_n + j_z_n*j_z_n);
    j_x_n= j_x_n/jn_norm;
    j_y_n = j_y_n/jn_norm;
    j_z_n = j_z_n/jn_norm;
  }
  return i_z_n*x_n + j_z_n*y_n + k_z_n*z_n;
}


