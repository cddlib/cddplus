/* rlp.h: Header file for rlp.c 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
*/

#include <stdio.h>
#include <stdlib.h>
#define NMAX 202 
#define MMAX 20002
#ifndef RAND_MAX
#define RAND_MAX 32767
#endif

#define TRUE 1
#define FALSE 0

typedef int boolean;
typedef long rowrange;
typedef long colrange;
typedef float *Arow;
typedef float *Amatrix[MMAX];

void GenerateKuhnQuandt(Amatrix, long, long);
void GenerateSignedKuhnQuandt(Amatrix, long, long);
void GenerateDegenerate(Amatrix, long, long);
void GenerateZeroOne(Amatrix, long, long);
void RandomRowVector(Arow, long, int, int);
void WriteAmatrix(Amatrix,long,long,unsigned, int);
 
