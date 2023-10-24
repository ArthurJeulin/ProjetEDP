#include<time.h>
#include<math.h>
#include<complex.h>
#include<stdio.h>
#include <stdbool.h> 
#include<stdlib.h>
#include<string.h>
#include"subf.h"
#include"VTSW.h"
#include"data.h"

struct state{
  int maillage; // 0 = Uniforme ; 1 = Non Uniforme
  int vitesse ; // 0 = Uniforme ; 1 = Non Uniforme
}defaultState;

int main()
{
struct type_donneesc param;
float **x,**y,**xv,**yv,**vol;
float **T0,**T1;
float **U;     // U definie au centre des facettes est-ouest en x,yv
float **V;     // V definie au centre des facettes nord-sud en xv,y
float **Fadv;
float **Fdiff;
float dt,tf;
int l,N;

float **A;   //matrice A
float *B;    //Vecteur second membre
float **AB;  //matrice bande
float *X;    //vecteur qui contiendra la solution de SOR
int LB;      //Largeur de la bande;

printf("simulation started\n");
param=read_datac();
int tailleMat=(param.nx*param.ny);
printf("-------------------------------------------------------------------\n");
printf(" param.nx : %d  param.ny : %d\n param.Lx : %f m param.Ly : %f m\n param.U0 : %f s^-1 \n param.D : %f m^2/s\n param.Ti : %f ( Initial temperature)\n param.Tg : %f Temperature at the left (Celsius)\n param.Tb : %f Temperature at the bottom (Celsius) \n param.tf : %f final time (s)\n param.Nout : %d number of intermediate unsaved time steps \n param.CFL : %f Courant's number (advection) \n param.R : %f ourier's number (diffusion)\n",
param.nx,param.ny,param.Lx,param.Ly,param.U0,param.D,param.Ti,param.Tg,param.Tb,param.tf,param.Nout,param.CFL,param.R);
printf("-------------------------------------------------------------------\n");
printf("Choix du maillage et de la vitesse : \n");
printf("Definir maillage ===>\n 0 = Uniforme (x,y)\n 1 = Non Uniforme suivant x\n 2 = Non Uniforme suivant y\n 3 = Non Uniforme suivant (x,y)\n Entrer une valeur : ");
scanf("%d", &defaultState.maillage);

if (defaultState.maillage == 0 || defaultState.maillage == 1|| defaultState.maillage == 2|| defaultState.maillage == 3) {
       // printf("Vous avez entré : %d\n", defaultState.maillage);
} 
else {
    printf("Vous avez entré : %d\n", defaultState.maillage);
    printf("La valeur entrée n'est pas bonne.\n");
    exit(1);
}
printf("Definir vitesse ===> \n 0 = Uniforme \n 1 = Non Uniforme \n Entrer une valeur : ");
scanf("%d", &defaultState.vitesse);

if (defaultState.vitesse== 0 || defaultState.vitesse == 1) {
        //printf("Vous avez entré : %d\n", defaultState.vitesse);
} 
else {
    printf("Vous avez entré : %d\n", defaultState.vitesse);
     printf("La valeur entrée n'est ni 1 ni 0.\n");
     exit(1);
}
printf("-------------------------------------------------------------------\n");
//printf("maillage : %d \n vitesse : %d\n",defaultState.maillage,defaultState.vitesse);

// Allocation dynamique
    x=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {x[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    y=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {y[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    xv=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {xv[i]=(float*)malloc((param.ny)*sizeof(float));}
    yv=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {yv[i]=(float*)malloc((param.ny)*sizeof(float));}
    vol=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {vol[i]=(float*)malloc((param.ny)*sizeof(float));}
    T0=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {T0[i]=(float*)malloc((param.ny)*sizeof(float));}
    T1=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {T1[i]=(float*)malloc((param.ny)*sizeof(float));}
    U=(float**)malloc((param.nx+1)*sizeof(float*));
    for (int i=0;i<param.nx+1;i++) {U[i]=(float*)malloc((param.ny)*sizeof(float));}
    V=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {V[i]=(float*)malloc((param.ny+1)*sizeof(float));}
    Fadv=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fdiff=(float**)malloc((param.nx)*sizeof(float*));
    for (int i=0;i<param.nx;i++) {Fdiff[i]=(float*)malloc((param.ny)*sizeof(float));}

meshc(param,x,y,xv,yv,vol,defaultState.maillage,defaultState.vitesse);


initial_conditionc(param,xv,yv,x,y,T0,U,V,defaultState.maillage,defaultState.vitesse);
VTSWriterc(0.0,0,param.nx+1,param.ny+1,x,y,T0,U,V,"ini");

dt= calc_dtc(param.nx,param.ny,x,y,U,param.CFL,V,param.D,param.R);

N=(int)(param.tf/dt);

printf("N=%d\n",N);

printf("Nout=%d\n",param.Nout);

// Computation with the explicit scheme
if(param.i_solver==0)
  {
  printf("Simulation with the explicit solver...\n");
  printf("time advancing...\n");
  for (l=1;l<N;l++)
    {
    printf("Iteration l=%d...\n",l);
    calc_flux_advc(param,x,y,xv,yv,U,V,T0,Fadv);
    calc_flux_diffc(param,x,y,xv,yv,T0,Fdiff);
    advance_timec(param,dt,vol,Fadv,Fdiff,T0,T1);

    if((l%param.Nout)==0)
      {
      VTSWriterc((float)(l)*dt,l,param.nx+1,param.ny+1,x,y,T1,U,V,"int");
      }
    }
  VTSWriterc((float)(N)*dt,N,param.nx+1,param.ny+1,x,y,T1,U,V,"end");
  }

// Computation with the implicit scheme (Gauss method)
if((param.i_solver)==1)
  {
  int NA;
  printf("Simulation with the implicit solver Gauss...\n");

  A=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){A[i]=(float*)malloc((tailleMat)*sizeof(float));}
  if (A==NULL) {printf("Allocation failed");}

  B=(float*)malloc((tailleMat)*sizeof(float*));

  printf("Gauss implemented\n");

  creation_A(param, NA, dt, x,y,xv,yv,vol,A);
  
  //On itère pour chaque pas de temps
  //tf = temps final et dt pas de temps 
  // N=(int)(param.tf/dt);
  for (l=1;l<N;l++)
    {
    printf("Iteration l=%d...\n",l);
    calc_flux_advc(param,x,y,xv,yv,U,V,T0,Fadv);
    //calc_flux_diffc(param,x,y,xv,yv,T0,Fdiff);
    creation_B(param, NA,dt, x, y,xv,yv,vol, Fadv,T0, B);
    
    //A la sortie de GAUSS, la solution se trouve dans B.
    gaussij(param.nx*param.ny, A, B);
    miseajour_T(param,T0,T1,B);
    
    if((l%param.Nout)==0)
      {
      VTSWriterc((float)(l)*dt,l,param.nx+1,param.ny+1,x,y,T1,U,V,"int");
      }
    }
   // (to be completed)

  free(A);
  free(B);
  }


// Computation with the implicit scheme (SOR method)
if((param.i_solver)==2)
  {
  printf("Simulation with the implicit solver SOR...\n");

  A=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){A[i]=(float*)malloc((tailleMat)*sizeof(float));}
  if (A==NULL) {printf("Allocation failed");}

  B=(float*)malloc((tailleMat)*sizeof(float*));
  LB=(2*param.nx)+1;

  AB=(float**)malloc((tailleMat)*sizeof(float *));
  for (int i=0;i<tailleMat;i++){AB[i]=(float*)malloc((LB)*sizeof(float));}
  if (AB==NULL) {printf("Allocation failed for AB");}
  X=(float*)malloc((tailleMat)*sizeof(float*));

  // (to be completed)
  printf("not implemented\n");


  free(A);
  free(B);
  free(AB);
  free(X);
  }

free(x);
free(y);
free(vol);
free(T0);
free(T1);
free(U);
free(V);
free(Fadv);
free(Fdiff);

printf("Simulation done\n");

return 0;

}


