#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
#include"data.h"
#define M_PI       3.14159265358979323846


//*********************************
struct type_donneesc read_datac(){

   struct    type_donneesc param;
   FILE *fp;
   char str[1000], str_copy[1000];
   int iter;
   int read_int;
   float read_float;

   printf("reading data...\n");
   fp = fopen("data.txt", "r");
   // lit les donnees à partir du fichier
      iter = 0;

       while ( fgets (str, 1000, fp) != NULL ) {

           strncpy(str_copy, str, 10);

           if ( iter == 0 || iter == 1 || iter == 10 || iter == 13 || iter == 14 || iter == 15  ){
               read_int = atoi(str_copy);
           }
           else{
               read_float = atof(str_copy);
           }

           switch (iter){

               case 0:
                   param.nx = read_int;
                   break;

               case 1:
                   param.ny = read_int;
                   break;

               case 2:
                   param.Lx = read_float;
                   break;

               case 3:
                   param.Ly = read_float;
                   break;

               case 4:
                   param.U0 = read_float;
                   break;

               case 5:
                   param.D = read_float;
                   break;

               case 6:
                   param.Ti = read_float;
                   break;

               case 7:
                   param.Tg = read_float;
                   break;

               case 8:
                   param.Tb = read_float;
                   break;

               case 9:
                   param.tf = read_float;
                   break;

               case 10:
                   param.Nout = read_int;
                   break;

               case 11:
                   param.CFL = read_float;
                   break;

               case 12:
                   param.R = read_float;
                   break;

               case 13:
                   param.i_mesh = read_int;
                   break;

               case 14:
                   param.i_vit = read_int;
                   break;

               case 15:
                   param.i_solver = read_int;
                   break;

               case 16:
                   param.W = read_float;
                   break;

               case 17:
                   param.R0 = read_float;
                   break;

           }

           iter++;

       }

       fclose(fp);

       return param;

   }


//********************************************************************************************
void meshc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **vol,int maillage, int vitesse) {

int i,j;

printf("creating mesh...\n");
//printf("maillage : %d \n vitesse : %d\n",maillage,vitesse);
// Uniform mesh (streching-compressing flow + rotating flow)
//for (i=0;i<param.nx+1;i++){
//    for (j=0;j<param.ny+1;j++){
//     x[i][j]=-param.Lx+2.*param.Lx*(float)(i)/(float)(param.nx);
//     y[i][j]=-param.Ly+2.*param.Ly*(float)(j)/(float)(param.ny);
//    }
//}

// Uniform mesh (parabolic flow)
/*
x[i][j] correspond au boards de la case
xv[i][j] correspond au maillage horizontal des centres des cases
*/
switch (maillage) {
    case 0:
        for (i=0;i<param.nx+1;i++){
            for (j=0;j<param.ny+1;j++){
            //on calcule x,y regulier
            x[i][j]=param.Lx*(float)(i)/(float)(param.nx);
            y[i][j]=param.Ly*(float)(j)/(float)(param.ny);
            }
        }
        break;
    case 1:
        for (i=0;i<param.nx+1;i++){
            for (j=0;j<param.ny+1;j++){
            //on calcule x,y regulier
            x[i][j]=param.Lx*(float)(i)/(float)(param.nx);
            y[i][j]=param.Ly*(float)(j)/(float)(param.ny);
            //on calcule x irregulier
            x[i][j]= (x[i][j]*x[i][j])/param.Lx;
            }
        }
        break;
    case 2:
        for (i=0;i<param.nx+1;i++){
            for (j=0;j<param.ny+1;j++){
            //on calcule x,y regulier
            x[i][j]=param.Lx*(float)(i)/(float)(param.nx);
            y[i][j]=param.Ly*(float)(j)/(float)(param.ny);
            //on calcule y irregulier
            y[i][j]=(param.Ly)*(1-cos(M_PI*y[i][j]/(2*param.Ly)));
            }
        }
        break;
    case 3: 
        for (i=0;i<param.nx+1;i++){
            for (j=0;j<param.ny+1;j++){
            //on calcule x,y regulier
            x[i][j]=param.Lx*(float)(i)/(float)(param.nx);
            y[i][j]=param.Ly*(float)(j)/(float)(param.ny);
            //on calcule x, y irregulier
            x[i][j]= (x[i][j]*x[i][j])/param.Lx;
            y[i][j]=(param.Ly)*(1-cos(M_PI*y[i][j]/(2*param.Ly)));
            }
        }
        break;
    default:
        printf("Erreur avec maillage.\n");
        exit(1);
        break;
    }


for (i=0;i<param.nx;i++){
    for( j=0;j<param.ny;j++){
          xv[i][j]=0.5*(x[i][j]+x[i+1][j]);
          yv[i][j]=0.5*(y[i][j]+y[i][j+1]);
          vol[i][j]=(x[i+1][j]-x[i][j])*(y[i][j+1]-y[i][j]);
    }
}

   }




//****************************************************************************************************************************
void initial_conditionc(struct type_donneesc param,float **xv,float **yv,float **x,float **y,float **T0,float **U,float **V,int maillage, int vitesse){

    int i,j;

    printf("initial condition...\n");
    //printf("maillage : %d \n vitesse : %d\n",maillage,vitesse);

    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            T0[i][j]=param.Ti; // Condition Initial pour la conduite à t = 0 pour(x,y)
        }
    }

    switch (vitesse) {
        case 0:
            for (i=0;i<param.nx+1;i++){
                for (j=0;j<param.ny;j++){
                    U[i][j] = param.U0 ;
                }
            }
            break;
        case 1:
            for (i=0;i<param.nx+1;i++){
                for (j=0;j<param.ny;j++){
                    //Champ de vitesse selon x
                    U[i][j] =(6*param.U0)*(y[i][j]/(2*param.Ly))*(1-(y[i][j]/(2*param.Ly)));
                }
            }
            break;
        default:
            printf("Erreur avec vitesse.\n");
            exit(1);
            break;
    }

    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny+1;j++){
            V[i][j] = 0.;
        }
    }


      }


//*****************************************************************************************************
float calc_dtc(int nx, int ny,float **x,float **y,float **U,float CFL,float **V,float D,float R) {

 float dx,dy;
 float dt,dt_loc;
 int  i,j;

 printf("computing dt...\n");

 dt = 1.e8;

 for( j=0; j<ny; j++){
   for (i=0; i<nx; i++){
     dx = x[i+1][j]-x[i][j];
     dy = y[i][j+1]-y[i][j];
     dt_loc = 1.0/(fabs(U[i][j])/(CFL*dx)+fabs(V[i][j])/(CFL*dy)+D*(1.0/(dx*dx)+1.0/(dy*dy))/R);
    if (dt_loc < dt) { dt = dt_loc; }
   }
 }

 printf("dt=%f\n",dt);

return dt;

   }



//************************************************************************************************************************************
void calc_flux_advc(struct type_donneesc param,float **x,float **y,float **xv,float **yv,float **U,float **V,float **T0,float **Fadv)
{

    float **Fadv_o, **Fadv_e, **Fadv_n, **Fadv_s;
    float T_amont;
    int i,j,k;
    Fadv_e=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_e[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_o=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_o[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_n=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_n[i]=(float*)malloc((param.ny)*sizeof(float));}
    Fadv_s=(float**)malloc((param.nx)*sizeof(float *));
    for (int i=0;i<param.nx;i++) {Fadv_s[i]=(float*)malloc((param.ny)*sizeof(float));}

    // Advection term

   // East flux
   // In the domain
    for (i=0;i<param.nx-1;i++){
        for (j=0;j<param.ny;j++){
            if (U[i+1][j]>=0.0){
                T_amont = T0[i][j];}
            else{
                T_amont = T0[i+1][j];
            }
            Fadv_e[i][j] = -1.0*U[i+1][j]*T_amont*(y[i+1][j+1]-y[i+1][j]);
        }
    }
   // BC
    i=param.nx-1;
    for (j=0;j<param.ny;j++){
        if (U[i+1][j]>=0.0){
            T_amont = T0[i][j];}
        else{
            printf("BC not implemented 01 \n");
            break;
        }
        Fadv_e[i][j] = -1.0*U[i+1][j]*T_amont*(y[i+1][j+1]-y[i+1][j]);
    }
   //West flux
   //In the domain
    for (i=1;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            Fadv_o[i][j] = -1.0*Fadv_e[i-1][j];
        }
    }
   // BC
    i=0;
    for (j=0;j<param.ny;j++){
        if (U[i][j]<=0.0){
            T_amont = T0[i][j];}
        else{
            T_amont = param.Tg;
        }
        Fadv_o[i][j] = U[i][j]*T_amont*(y[i][j+1]-y[i][j]);
    }
   // North flux
   // In the domain
    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny-1;j++){
            if (V[i][j+1]>=0.0){
                T_amont = T0[i][j];}
            else{
                T_amont = T0[i][j+1];
            }
            Fadv_n[i][j] = -1.0*V[i][j+1]*T_amont*(x[i+1][j+1]-x[i][j+1]);
        }
    }
   // BC
    j=param.ny-1;
    for (i=0;i<param.nx;i++){
        if (V[i][j+1]>=0.0){
            T_amont = T0[i][j];}
       else{
            T_amont = 0.0;
       }
       Fadv_n[i][j] = -1.0*V[i][j+1]*T_amont*(x[i+1][j+1]-x[i][j+1]);
    }
   //South flux
   //In the domain
    for (i=0;i<param.nx;i++){
        for (j=1;j<param.ny;j++){
            Fadv_s[i][j] = -1.0*Fadv_n[i][j-1];
        }
    }
   // BC
    j=0;
    for (i=0;i<param.nx;i++){
        if (V[i][j]<=0.){
            T_amont = T0[i][j];}
        else{
            T_amont = param.Tb;
        }
        Fadv_s[i][j] = V[i][j]*T_amont*(x[i+1][j]-x[i][j]);

    }

   //c) Total summation

    for (i=0;i<param.nx;i++){
        for (j=0;j<param.ny;j++){
            Fadv[i][j] = Fadv_o[i][j] + Fadv_e[i][j] + Fadv_s[i][j] + Fadv_n[i][j];
        }
    }

   }

//*********************************************************************************************************************************
void calc_flux_diffc(struct type_donneesc param,float **x, float **y,float **xv,float **yv,float **T0, float **Fdiff) {

float  **Fdiff_o, **Fdiff_e, **Fdiff_s, **Fdiff_n;
int  i,j,k;
Fdiff_e=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_e[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_o=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_o[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_n=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_n[i]=(float*)malloc((param.ny)*sizeof(float));}
Fdiff_s=(float**)malloc((param.nx)*sizeof(float *));
for (int i=0;i<param.nx;i++) {Fdiff_s[i]=(float*)malloc((param.ny)*sizeof(float));}

// Diffusive term
// West flux
// In the domain
for (i=1;i<param.nx;i++){
    for (j=0;j<param.ny;j++){
        Fdiff_o[i][j] = -1.0*param.D*(T0[i][j]-T0[i-1][j])/(xv[i][j]-xv[i-1][j])*(y[i][j+1]-y[i][j]);
    }
}
//BC
i=0;
for (j=0;j<param.ny;j++){
    Fdiff_o[i][j] = -1.0*param.D*(T0[i][j]-param.Tg)/(2.*(xv[i][j]-x[i][j]))*(y[i][j+1]-y[i][j]);
}
//East flux
// In the domain
for (i=0;i<param.nx-1;i++){
    for (j=0;j<param.ny;j++){
        Fdiff_e[i][j] = -1.0*Fdiff_o[i+1][j];
    }
}
//BC
i=param.nx-1;
for (j=0;j<param.ny;j++){
    Fdiff_e[i][j] = -1.0*Fdiff_o[i][j];
}
// South flux
//In the domain
for (i=0;i<param.nx;i++){
    for (j=1;j<param.ny;j++){
        Fdiff_s[i][j] = -1.0*param.D*(T0[i][j]-T0[i][j-1])/(yv[i][j]-yv[i][j-1])*(x[i+1][j]-x[i][j]);
    }
}
//BC
j=0;
for (i=0;i<param.nx;i++){
    Fdiff_s[i][j] = -1.0*param.D*(T0[i][j]-param.Tb)/(2.*(yv[i][j]-y[i][j]))*(x[i+1][j]-x[i][j]);
}
// North flux
// In the domain
for (i=0;i<param.nx;i++){
    for (j=0;j<param.ny-1;j++){
        Fdiff_n[i][j] = -1.0*Fdiff_s[i][j+1];
    }
}
//BC
j=param.ny-1;
for (i=0;i<param.nx;i++){
    Fdiff_n[i][j] = 0.0;
}

//c) Total summation

for (i=0;i<param.nx;i++){
    for (j=0;j<param.ny;j++){
        Fdiff[i][j] = Fdiff_o[i][j] + Fdiff_e[i][j] + Fdiff_s[i][j] + Fdiff_n[i][j];
    }
}

}

//*********************************************************************************************************************
void  advance_timec(struct type_donneesc param,float dt, float **vol,float **Fadv,float **Fdiff,float **T0,float **T1){

int i,j;

for (i=0;i<param.nx;i++)
{
    for (j=0;j<param.ny;j++)
    {
        T1[i][j] = T0[i][j] + dt/vol[i][j]*(Fadv[i][j] + Fdiff[i][j]);
    }
}
for (i=0;i<param.nx;i++)
{
    for  (j=0;j<param.ny;j++)
    {
        T0[i][j]=T1[i][j];
    }
}
}
// difference horizontal au niveau des bords
// difference horizontal au centre de la cellule de calcul
float deltaX(int i, int j,float **x ){
 return x[i+1][j] - x[i][j];
//deltaX = xv[i+1][j] - xv[i][j];
}
// difference vertical au niveau des bords
// difference vertical au centre de la cellule de calcul
float deltaY(int i, int j,float **y){
    return y[i][j+1] - y[i][j];
//deltaY = yv[i][j+1] - yv[i][j];
}
//*********************************************
void creation_A(struct type_donneesc param,int NA, float dt, float **x,float **y,float **xv,float **yv,float **vol,float **A)
{

    int k, i, j ;
    float a,b,c,d,e;
    printf("param.nx : %d param.ny : %d \n",param.nx,param.ny);
    //Cas dans le domaine
	for (int i=1;i<param.nx-1;++i){
        for (int j=1;j<param.ny-1;++j){
            k = j*param.nx + i;
            printf("k value le domaine: %d\n",k);
            //Calcul des coefficient (a,b,c,d,e,f)
            // a = ((dt/vol[i][j])*(param.D)*( (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/deltaY(i,j-1,yv))+ (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/deltaX(i-1,j,xv))  ));
            // b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
            // c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
            // d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
            // e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
            A[k][k] = 1 + a ;
            A[k][k+1] = 0 ;
            A[k][k-1] = e;
            A[k][k+param.nx] = b;
            A[k][k-param.nx] = c;
            // A[k][k] = 1 ;
            // A[k][k+1] = 2 ;
            // A[k][k-1] = 3;
            // A[k][k+param.nx] = 4;
            // A[k][k-param.nx] = 5;
            

        }
    }
    //Segment [OC]
    /*
    attention a delta(i-1)
    Ti-1 = Te
    e passee dans B
    */
    i = 0; 
	for (int j=1;j<param.ny-1;++j){
        k = j*param.nx + i;
        printf("k value le Segment [OC]: %d\n",k);
        //Calcul des coefficient (a,b,c,d,e,f)
        a = ((dt/vol[i][j])*(param.D)*( (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/deltaY(i,j-1,yv)) + (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/(deltaX(i,j,xv))/2)  ));
        b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
        c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
        d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
        //e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/(deltaX(i,j,xv))/2) ));
        A[k][k] = 1 + a ;
        A[k][k+1] = d ;
        //A[k][k-1] = e;
        A[k][k+param.nx] = b;
        A[k][k-param.nx] = c;
        // A[k][k] = 1  ;
        // A[k][k+1] = 2 ;
        // //A[k][k-1] = e;
        // A[k][k+param.nx] = 4;
        // A[k][k-param.nx] = 5;
    }

    //Segment [OA]
    /*
    attention a deltaY(j-1)
    Tj-1 = Tp
    c passee dans B
    */
    j = 0; 
	for (int i=1;i<param.nx-1;++i){
        k = j*param.nx + i;
        //Calcul des coefficient (a,b,c,d,e,f)
        a = ((dt/vol[i][j])*(param.D)*( (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/(deltaY(i,j,yv)/2)) + (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/deltaX(i-1,j,xv))  ));
        b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
      //  c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
        d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
        e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/(deltaX(i,j,xv))/2) ));
        // A[k][k] = 1 + a ;
        // A[k][k+1] = d ;
        // A[k][k-1] = e;
        // A[k][k+param.nx] = b;
        // //A[k][k-param.nx] = c;
        A[k][k] = 1 +a ;
        A[k][k+1] = d ;
        A[k][k-1] = e;
        A[k][k+param.nx] = b;
        //A[k][k-param.nx] = c;
        
        // A[k][k] = 1  ;
        // A[k][k+1] = 2 ;
        // A[k][k-1] = 3;
        // A[k][k+param.nx] = 4;
        // //A[k][k-param.nx] = c;

        printf("k value le Segment [OA]: %d\n",k);
        
    }
    //Segment [CB]
    /*
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    j = param.ny - 1; 
	for (int i=1;i<param.nx-1;++i){
        k = j*param.nx + i;
        //Calcul des coefficient (a,b,c,d,e,f)
        a = ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv))+ (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/deltaX(i-1,j,xv)) ));
        //b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
        c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
        d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
        e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/(deltaX(i,j,xv))/2) ));
        A[k][k] = 1 + a ;
        A[k][k+1] = d ;
        A[k][k-1] = e;
        //A[k][k+param.nx] = b;
        A[k][k-param.nx] = c;
        // A[k][k] = 1  ;
        // A[k][k+1] = d ;
        // A[k][k-1] = e;
        // //A[k][k+param.nx] = b;
        // A[k][k-param.nx] = c;
       printf("k value le Segment [CB]: %d\n",k);
    }
    //Segment [AB]
    /*
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest  
    */
    i = param.nx-1;
    for (int j=1;j<param.ny-1;++j){
        k = j*param.nx + i;
        printf("k value le Segment [AB]: %d\n",k);
        //Calcul des coefficient (a,b,c,d,e,f)
        a = ((dt/vol[i][j])*(param.D)*( (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
        b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
        c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
        //d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
        //e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
        A[k][k] = 1 + a ;
        //A[k][k+1] = d ;
        //A[k][k-1] = e;
        A[k][k+param.nx] = b;
        A[k][k-param.nx] = c;

    }
    //Coin 0
    /*
    Problème en Ti-1,j et Ti,j-1
    Condition Limite T(x= 0,y = 0) = Te
    */
    i = 0 ;
    j = 0 ;
    k = j*param.nx + i;
   printf("k value le Coin 0: %d\n",k);
    //Calcul des coefficient (a,b,c,d,e,f)
    a = ((dt/vol[i][j])*(param.D)*( (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/(deltaY(i,j,yv)/2))+ (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/(deltaX(i,j,xv)/2))  ));
    b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
    //c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
    d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
    //e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
    A[k][k] = 1 + a ;
    A[k][k+1] = d ;
    //A[k][k-1] = e;
    A[k][k+param.nx] = b;
    //A[k][k-param.nx] = c;
    
    //Coin C
    /*
    Problème en Ti-1,j et Ti,j+1
    Condition Limite T(x= 0,y = Ny-1) = Te
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    i = 0 ;
    j = param.ny - 1;
    k = j*param.nx + i;
    printf("k value le Coin C: %d\n",k);
    //Calcul des coefficient (a,b,c,d,e,f)
    a = ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/(deltaY(i,j-1,yv)))+ (deltaY(i,j,y)/deltaX(i,j,xv)) + (deltaY(i,j,y)/((deltaX(i,j,xv)/2))) ));
    //b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
    c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
    d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
    //e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
    A[k][k] = 1 + a ;
    A[k][k+1] = d ;
    //A[k][k-1] = e;
    //A[k][k+param.nx] = b;
    A[k][k-param.nx] = c;
    
    //Coin A
    /*
    Problème en Ti-1,j et Ti,j+1
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest  
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    i = param.nx - 1 ;
    j = 0;
    k = j*param.nx + i;
    printf("k value le Coin A: %d\n",k);
    //Calcul des coefficient (a,b,c,d,e,f)
    a =   ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) + (deltaX(i,j,x)/(deltaY(i,j,yv)/2) )));
    b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
    // c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
    // d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
    // e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
    A[k][k] = 1 + a ;
    // A[k][k+1] = d ;
    // A[k][k-1] = e;
    A[k][k+param.nx] = b;
    //A[k][k-param.nx] = c;
        
    //Coin B
    /*
    Problème en Ti+1,j et Ti,j+1
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest  
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    printf("k value le Coin B: %d\n",k);
    i = param.nx - 1 ;
    j = param.ny - 1;
    k = j*param.nx + i;
    //Calcul des coefficient (a,b,c,d,e,f)
    a =   ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
    //b = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j,yv)) ));
    c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/deltaY(i,j-1,yv)) ));
    // d = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/deltaX(i,j,xv)) ));
    // e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/deltaX(i-1,j,xv)) ));
    A[k][k] = 1 + a ;
    // A[k][k+1] = d ;
    // A[k][k-1] = e;
    //A[k][k+param.nx] = b;
    A[k][k-param.nx] = c;
}


//**********************************


void creation_B(struct type_donneesc param, int NA, float dt, float **x, float **y,float **xv,float **yv,float **vol, float **Fadv, float **T0, float *B)
    {
    int k;
    float a,b,c,d,e;        
    //Cas dans le domaine
	for (int i=1;i<param.nx-1;++i){
        for (int j=1;j<param.ny-1;++j){
            k = j*param.nx + i;
            //Rappel B juste Vecteur colonne de taille Ny
            B[k] = T0[i][j] +(dt/vol[i][j])*Fadv[i][j] ;
            printf("k value le domaine: => B[%d] : %f\n",k, B[k]);
        }
    }


    //Segment [OC]
    /*
    attention a delta(i-1)
    Ti-1 = Te
    e passe dans B
    */
    int i = 0; 
	for (int j=1;j<param.ny-1;++j){
        k = j*param.nx + i;
        e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/(deltaX(i,j,xv))/2) ));
        B[k] = -e*param.Tg + T0[i][j] +(dt/vol[i][j])*Fadv[i][j];
        printf("k value le Segment [OC]: => B[%d] : %f\n",k,B[k]);
    }
    //Segment [OA]
    /*
    attention a deltaY(j-1)
    Tj-1 = Tp
    c passee dans B
    */
    int j = 0; 
	for (int i=1;i<param.nx-1;++i){
        k = j*param.nx + i;
        c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/(deltaY(i,j,yv)/2)) ));
        B[k] = -c*param.Tb + T0[i][j] +(dt/vol[i][j])*Fadv[i][j];
        printf("k value le Segment [OA]: => B[%d] : %f\n",k,B[k]);
    }
    //Segment [CB]
    /*
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    j = param.ny - 1;
	for (int i=1;i<param.nx-1;++i){
        k = j*param.nx + i;
        B[k] = T0[i][j] +(dt/vol[i][j])*Fadv[i][j] ;
        printf("k value le Segment [CB]: => B[%d] : %f\n",k,B[k]);
        
    }
    //Segment [AB]
    /*
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest 
    */
    i = param.nx-1;
    for (int j=1;j<param.ny-1;++j){
        k = j*param.nx + i;
        B[k] = T0[i][j] +(dt/vol[i][j])*Fadv[i][j] ;
        printf("k value Segment [AB]: => B[%d] : %f\n",k,B[k]);
    }
    //Coin 0
    /*
    Problème en Ti-1,j et Ti,j-1
    Condition Limite T(x= 0,y = 0) = Te
    */
    i = 0 ;
    j = 0 ;
    k = j*param.nx + i;
    //Calcul des coefficient (a,b,c,d,e,f)
    c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/(deltaY(i,j,yv)/2)) ));
    e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/(deltaX(i,j,xv)/2)) ));
    B[k] = -c*param.Tb -e*param.Tg + T0[i][j] +(dt/vol[i][j])*Fadv[i][j];
    printf("k value le Coin 0: => B[%d] : %f\n",k,B[k]);

    //Coin C
    /*
    Problème en Ti-1,j et Ti,j+1
    Condition Limite T(x= 0,y = Ny-1) = Te
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    i = 0 ;
    j = param.ny - 1;
    k = j*param.nx + i;
    //e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,x)/(deltaX(i,j,xv)/2)) ));
    e = - ((dt/vol[i][j])*(param.D)*(  (deltaY(i,j,y)/(deltaX(i,j,xv))/2) ));
    B[k] = -e*param.Tg + T0[i][j] +(dt/vol[i][j])*Fadv[i][j];
    printf("k value le Coin C: => B[%d] : %f\n",k,B[k]);
   
    //Coin A
    /*
    Problème en Ti-1,j et Ti,j+1
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest  
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    i = param.nx - 1 ;
    j = 0;
    k = j*param.nx + i;
    //Calcul des coefficient (a,b,c,d,e,f)
    c = - ((dt/vol[i][j])*(param.D)*(  (deltaX(i,j,x)/(deltaY(i,j,yv)/2) )));
    B[k] = -c*param.Tb + T0[i][j] +(dt/vol[i][j])*Fadv[i][j];
    printf("k value le Coin A: => B[%d] : %f\n",k,B[k]);
    
    //Coin B
    /*
    Problème en Ti+1,j et Ti,j+1
    Condition Limite (dT/dx)(n) = (dT/dx)(n-1)
    On a flux diff Est = - flux diff Ouest  
    Condition Limite dT/dy = 0; Ce qui annule le flux diff Nord
    */
    i = param.nx - 1 ;
    j = param.ny - 1;
    k = j*param.nx + i;
    B[k] =  T0[i][j] +(dt/vol[i][j])*Fadv[i][j];    
    printf("k value le Coin B: => B[%d] : %f\n",k,B[k]);
}


//*****************************************

void miseajour_T(struct type_donneesc param,float **T0,float **T1,float *B)
{
//int  i,j;
int k = 0;
for (int i=0;i<param.nx;i++){
    for (int j=0;j<param.ny;j++){
        k = j*param.nx + i;
        T1[i][j] = B[i*param.ny + j] ;
        T0[i][j] = B[i*param.ny + j] ;
        float result = B[i*param.ny + j]; 
       // printf("B(%d) = %f",k,result);
    }
}

}

//***********************************


