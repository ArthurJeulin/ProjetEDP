#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>
//***********************************
//***********************************
//
//           PASSAGE D'UNE MATRICE PLEINE A UNE MATRICE BANDE
//
// LANGAGE : C
//
//  Donnees : Cette routine cree une matrice bande AB de dimension (NA,LB)
//            a partir d'une matrice pleine A de dimension (NA,NA) telle que:
//           NA  Ordre de la matrice A
//           LB Largeur de la bande = 2 MX + 1;
//            MX est appelee "demi largeur de bande"
//
//             La matrice A est de la forme:
//
//             | A(1,1)       *                  *                    *|
//             |  *                                                    |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |                                                       |
//             |  *                                    *       A(NA,NA)|
//
//
//             Les coefficients de la matrice AB correspondent a:
//
//AB(1,1) ... | AB(1,MX+1)   ...    AB(1,LB)                            |
//             |                                                         |
//             | AB(2,MX)      .....    AB(2,LB)                         |
//             |                                                         |
//             |           *       .....       *                         |
//             |                                                         |
//             |               *       .....       *                     |
//             |                                                         |
//             |                   AB(I,1) ... AB(I,MX+1) ... AB(I,LB)   |
//             |                                                         |
//             |                         *       .....       *           |
//             |                                                         |
//             |                             *       .....       *       |
//             |                                                         |
//             |                               AB(NA,1)  ...  AB(NA,MX+1)| ...  AB(NA,LB)

void create_A_band(int MX,int NA,float **A, float **AB){

    int i,j;
//printf("hala esaie1\n");
    for(i=0;i<NA;i++){
        for(j=0;j<2*MX+1;j++){
            AB[i][j] = 0.0;
        }
    }
//printf("hala esaie2\n");
    for(i=0;i<MX+1;i++){
//printf("hala esaieplz\n");
        for(j=0;j<NA-(i);j++){
//printf("hala esaie10\n");
            AB[j][(MX)+(i)] = A[j][j+(i)];
//printf("hala esaie11\n");
            AB[NA-1-j][(MX)-(i)] = A[NA-1-j][NA-1-j-(i)];
//printf("hala esaie12\n");
        }
   //printf("hala esaie3\n"); 
}
//printf("hala esaie4\n");
       return;

       }
//***************************

//************************************
//************************************
//           RESOLUTION D'UN SYSTEME LINEAIRE
//
//  METHODE : Methode de Sur-Relaxation Successive SOR.
//
//  LANGAGE : C
//
//  Donnees : A  Coefficients de la matrice bande, variable a deux dimensions
//              dont les valeurs numeriques doivent etre fournies conformement
//               au schema ci-dessous.
//              (Les points exterieurs a la bande ne sont pas pris en compte
//              lors du calcul).
//            B  Termes du 2eme membre, variable a un indice
//            X  A la sortie de SOR, la solution se trouve dans X
//           N  Ordre de la matrice A
//            LB Largeur de la bande = 2 MX + 1
//            MX est appelee "demi largeur de bande"
//            R0 precision souhaitee
//           W  parametre de sur-relaxation (0<W<2)
//
//            |                                                       |
// A(1,1) ... | A(1,MX+1)   ...    A(1,LB)                            |
//            |                                                       |
//            | A(2,MX)      .....    A(2,LB)                         |
//            |                                                       |
//            |           *       .....       *                       |
//            |                                                       |
//            |               *       .....       *                   |
//            |                                                       |
//            |                   A(I,1) ... A(I,MX+1) ... A(I,LB)    |
//            |                                                       |
//            |                         *       .....       *         |
//            |                                                       |
//            |                             *       .....       *     |
//            |                                                       |
//            |                                 A(N,1)  ...  A(N,MX+1)| ...  A(N,LB)

void SOR(int MX,int N,float **A,float *B,float R0,float W, float *X){

    int  i,cc,dd,compteur;
    float *Y = malloc(sizeof(int)*(2*MX+N));
    float *Y0= malloc(sizeof(int)*(2*MX+N));
    float ECART,C,D,Z,Z0;
    

    for(i=0;i<2*MX+N;i++){
        Y[i]=0.0;
        Y0[i]=0.0;
    }
    ECART=1.0;
    compteur = 0;
    do  {
        compteur = compteur + 1;
        for(i=MX;i<MX+N;i++){
          // -------> calcul de la somme inférieure appelée C
            C=0.0;
            for(cc=0;cc<MX;cc++){
                C=C+Y[i-1-cc]*A[i-MX][MX-cc-1];
            }
          // -------> calcul de la somme supérieure appelée D
            D=0.0;
            for(dd=0;dd<MX;dd++){
                D=D+Y0[i+dd+1]*A[i-MX][MX+1+dd];
            }
          // -------> calcul de X,k+1
            Y[i]=(1.-W)*Y0[i]+(W/A[i-MX][MX])*(B[i-MX]-C-D);
        }

        // -------> critère d'arrêt
        Z0=0.0;
        Z=0.0;
        for(i=MX;i<MX+N;i++){
            Z=Z+Y[i]*Y[i];
            Z0=Z0+Y0[i]*Y0[i];
        }
      
        ECART=fabs(pow(Z,0.5)-pow(Z0,0.5));
        
        for(i=MX; i<MX+N;i++){
            Y0[i]=Y[i];
        }

    }
    while ((ECART>R0)&&(compteur<=10000));

    if (compteur>=10000) {
        printf("Probleme de convergence du solveur SOR\n");
    }

      // -------> écriture dans X

    for(i=MX;i<MX+N;i++){
        X[i-MX]=Y[i];
    }

    return;

}
//********************

