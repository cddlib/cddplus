/* Voronoi preprocessor by Komei Fukuda, fukuda@ifor.math.ethz.ch */
/* This code converts a Polyhedra *.ext file with n d-dimensional points data
   to a Polyhedra ine file representing the Voronoi polyhedron of 
   of the given points.  That is, for each point data (1,p1,p2...,pd)
   \in R^(d+1), output the vector (p1^2+p2^2+...+pd^2, -2 p1, ... , -2pd, 1)
   which means the inequality
  (p1^2+p2^2+...+pd^2) x0 -2 p1 x1  ... -2pd xd + x(d+1) >= 0.
*/

#define FALSE 0
#define TRUE 1
#define NMAX 100

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef double SQmat[NMAX][NMAX];
typedef double Rowvec[NMAX];

void WriteReal(FILE *f, double x)
{
  long ix1,ix2,ix;

  ix1= fabs(x) * 10000. + 0.5;
  ix2= (fabs(x) + 0.5);
  ix2= ix2*10000;
  if ( ix1 == ix2) {
    if (x>0) {
      ix = x + 0.5;
    } else {
      ix = -x + 0.5;
      ix = -ix;
    }
    fprintf(f, " %2ld", ix);
  } else
    fprintf(f, " % .9E", x);
}

main()
{
  char s[30];
  FILE *reading,*writing;
  char inputfile[30],outputfile[30], numbtype[30];
  int found;
  long mm,nn,i,j;
  double value,height;
  SQmat rmat;
  Rowvec aa;

  printf("Voronoi polyhedron preprocessor for cdd\n");
  printf("Input (*.ext) file of given points: ");
  scanf("%s",inputfile);
  printf("Output (*_vo.ine) file for voronoi polyhedron inequalities: ");
  scanf("%s",outputfile);

  reading=fopen(inputfile,"r");
  writing=fopen(outputfile,"w");

  found=FALSE;
  while (!found)
  {
   	  if (fscanf(reading,"%s",s)==EOF) {
   	    printf("improper input format\n");
  	    goto _L99;
  	  }
  	  else if (strncmp(s, "begin", 5)==0) {
  		found=TRUE;
  	  }
  }
  fscanf(reading, "%ld %ld %s", &mm, &nn, numbtype);
  printf("size = %ld x %ld\nNumber Type = %s\n", mm, nn, numbtype);
  fprintf(writing,"H-representation\nbegin\n");
  fprintf(writing," %ld   %ld   real\n",mm,nn+1);
  for (i = 1; i <= mm; i++) {
    for (j = 0; j < nn; j++) {
      fscanf(reading, "%lf", &value);
      aa[j] = value;
    /* printf("a(%3ld,%5ld) = %10.4lf\n",i,j,value); */
    }  /*of j*/
    /* putchar('\n'); */
    height=0;
    for (j = 1; j < nn; j++) {
	  height=height+aa[j]*aa[j];
    }  /*of j*/
    WriteReal(writing,height);
    for (j = 1; j < nn; j++) {
	  WriteReal(writing,-2*aa[j]);
    }  /*of j*/
    fprintf(writing,"  1\n");   
  }  /*of i*/
  fprintf(writing,"end\n");
  fprintf(writing,"incidence\n");
  fprintf(writing,"input_adjacency\n");
  printf("file %s is generated.\n", outputfile);
_L99:;
  fclose(reading);
  fclose(writing);
}


