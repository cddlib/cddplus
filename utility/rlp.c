/* Random LP generator by Komei Fukuda, fukuda@ifor.math.ethz.ch */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rlp.h"

double RowNorm(Amatrix A, long mm, long nn, long i)
{
  long j;
  float norm=0;

  if (i>=1 && i<=mm){
    for (j=2; j<=nn ;j++) {
      norm += A[i-1][j-1]*A[i-1][j-1];
    }
  } else {
    printf("row index i=%ld is outof range\n",i);
  }
  norm= sqrt(norm);
  return norm;
}

double ColNorm(Amatrix A, long mm, long nn, long j)
{
  long i;
  float norm=0;

  if (j>=1 && j<=nn){
    for (i=1; i<=mm ;i++) {
      norm +=  A[i-1][j-1]*A[i-1][j-1];
    }
  } else {
    printf("col index j=%ld is outof range\n",j);
  }
  norm=sqrt(norm);
  return norm;
}


void Normalize(Amatrix A, long m, long n)
{  /* This works only for dual KQ problems  */
  long i, j;
  float shift,a=100;
  Arow c; /* center to be (a,a,a,...,a) */
  
  c=(float *) calloc(n, sizeof shift);
  for (j=2; j<=n; j++){
    c[j-1]=a;  /* center */
  }
   
  for (i=n; i<=m; i++){
    /* Normalize i-th RHS */
    shift=0;
    for (j=2; j<=n; j++){
      shift+=c[j-1]*A[i-1][j-1];
    }
    A[i-1][0]=-(-a*RowNorm(A,m,n,i)+shift);
    /* This is to shift RHS so that it becomes tangent to
    the sphere centered at c with radius a. */
  }
  
  for (j=2; j<=n; j++){
    /* Normalize j-th OBJ value */
    A[m][j-1]=A[m][j-1]*ColNorm(A,m,n,j);
  }
}


void RandomRowVector(Arow row, long nn, int minvalue, int maxvalue)
{
  long j;

  for (j=1; j<=nn ;j++) {
    row[j-1] =  -( rand()%(maxvalue-minvalue+1) + minvalue);
  }
}

void GenerateKuhnQuandt(Amatrix A, long m, long n)
{  /* generate a dual Kuhn-Quandt instance */
  long i, j;
  int maxv=0,minv=-1000;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-10000;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=-10000;
  }
}

void GenerateSignedKuhnQuandt(Amatrix A, long m, long n)
{
  long i, j;
  int maxv=1000,minv=-1000;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-10000;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=10000;
  }
}

void GenerateZeroOne(Amatrix A, long m, long n)
{
  long i, j;
  int maxv=1, minv=-1;

  for (i=1; i<=m; i++){
    if (i<=n-1){
      A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=-1;
    }
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=1;
  }
}

void GenerateDegenerate(Amatrix A, long m, long n)
{
  long i, j;
  int maxv=1000, minv=-1000;

  for (i=1; i<=m; i++){
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=0;
  }
  for (j=1; j<=n; j++){
    if (j==1) A[m][j-1]=0; else A[m][j-1]=10000;
  }
}

void GenerateRandomFeasibility(Amatrix A, long m, long n)
{
  long i, j;
  long maxv=1000, minv=-1000;

  for (i=1; i<=m+1; i++){
    if (i<=n-1){
      if (i==1) A[i-1][0]=-1; else A[i-1][0]=0;
      for (j=2; j<=n; j++){
        if (j==i+1) A[i-1][j-1]=1; else A[i-1][j-1]=0;
      }
    } else {
      RandomRowVector(A[i-1], n, minv, maxv);
      A[i-1][0]=0;
    }
  }
}

void GenerateTotallyRandom(Amatrix A, long m, long n)
{
  long i, j;
  long maxv=1000, minv=-1000;

  for (i=1; i<=m+1; i++){
      RandomRowVector(A[i-1], n, minv, maxv);
      if (i==1) A[i-1][0]=-1; else A[i-1][0]=0;
  }
}

void WriteAmatrix(Amatrix A, long m, long n, unsigned seed, int lptype)
{
  long i, j;

  printf("LP type =  %d   Seed = %u\n", lptype, seed);
  printf("H-representation\nbegin\n  %ld  %ld  real\n", m, n);
  for (i=1; i<=m; i++){
    for (j=1; j<=n; j++){
      printf(" %1.0f", A[i-1][j-1]);
    }
    printf("\n");
  }
  printf("end\n");
  printf("maximize\n");
  for (j=1; j<=n; j++){
    printf(" %1.0f", A[m][j-1]);
  }
  printf("\n");
  if (lptype>=10) printf("feasibility\n");
  printf("!adjacency\n");
}

void WriteAmatrixDual(Amatrix A, long m, long n, unsigned seed, int lptype)
{
  /* This works only when the first (n-1) rows of A is identity matrix
  corresponding to nonnegativity of all variables.
  This applies to dual KQ problem, for example. */
  long i, j;
  long dm, dn;

  dm=m; dn=m-n+2;
  printf("LP type =  %d  Seed = %u\n", lptype, seed);
  printf("H-representation\nbegin\n  %ld  %ld  real\n", dm, dn);
  for (i=1; i<=dm-dn+1; i++){
    for (j=1; j<=dn; j++){
      if (j==1) printf(" %1.0f", -A[m][i]);
      else printf(" %1.0f", -A[n+j-3][i]);
    }
    printf("\n");
  }
  for (i=1; i<=dn-1; i++){
    for (j=1; j<=dn; j++){
      if (j==i+1) printf(" %1d", 1);
      else printf(" %1d", 0);
    }
    printf("\n");
  }
  
  printf("end\nmaximize\n");
  for (j=1; j<=dn; j++){
    printf(" %1.0f", -A[n+j-3][0]);
  }
  printf("\n");
  if (lptype>=10) printf("feasibility\n");
  printf("!adjacency\n");
}

int main(int argc, char *argv[])
{
  long mm,nn,i,j;
  unsigned sd=1;
  int lptype=1;
  Amatrix AA;
  float value;

  if (argc < 3) {
    printf("Too few arguments.\nUsage: rlp m n (seed) (type)\n");
    printf("where type must be either \n 1 (dual Kuhn-Quandt), 2 (Signed dKQ), 3 (-1,0,1),\n 4(normalized dKQ)  0 (sym random)\n 10 (Feasibility) and -i for their dual.\n");
  }
  else{
    mm=atol(argv[1]);
    nn=atol(argv[2]);
    if (nn > NMAX-1 || mm > MMAX-1) 
      printf("Size too large. Recompile with larger MMAX or NMAX\n");
    else {
      if (argc > 3){
        sd=atoi(argv[3]);
        if (argc > 4){
          lptype=atoi(argv[4]);
        }
      }
      srand(sd);
      for (i=0; i<=mm; i++){
        AA[i]=(float *) calloc(nn, sizeof value);
      }
      switch (lptype) {
        case 1:  case -1:
          GenerateKuhnQuandt(AA, mm, nn);
          break;

        case 2:  case -2:
          GenerateSignedKuhnQuandt(AA, mm, nn);
          break;

        case 3:  case -3:
          GenerateZeroOne(AA, mm, nn);
          break;

        case 4:  case -4:
          GenerateKuhnQuandt(AA, mm, nn);
          Normalize(AA, mm, nn);
          break;

        case 0:
          GenerateTotallyRandom(AA, mm, nn);
          break;

        case 10:
          GenerateRandomFeasibility(AA, mm, nn);
          break;

      }
      if (lptype >= 0)
        WriteAmatrix(AA,mm,nn,sd,lptype);
      else
        WriteAmatrixDual(AA,mm,nn,sd,lptype);
      }
  }
  return 1;
}

