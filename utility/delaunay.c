/* Delauany preprocessor by Komei Fukuda, fukuda@ifor.math.ethz.ch */
/* This code converts a Polyhedra *.ext file with m d-dimensional points data
   to a Polyhedra ine file containing the vertices and one ray
   of the Delaunay polyhedron.  That is, for each point data (1,p1,p2...,pd)
   \in R^(d+1), output the vector (1, p1, p2,...,pd, p1^2+p2^2+pd^2), and
   ray (0,0,...,1).
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

  printf("Delaunay polyhedron preprocessor for cdd\n");
  printf("Input (*.ext) file of given points: ");
  scanf("%s",inputfile);
  printf("Output (*_de.ext) file for delaunay polyhedron vertices: ");
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
  fprintf(writing,"V-representation\nbegin\n");
  fprintf(writing," %ld   %ld   real\n",mm+1,nn+1);
  for (i = 1; i <= mm; i++) {
    for (j = 0; j < nn; j++) {
      fscanf(reading, "%lf", &value);
      aa[j] = value;
    /* printf("a(%3ld,%5ld) = %10.4lf\n",i,j,value); */
    }  /*of j*/
    /* putchar('\n'); */
    height=0;
    for (j = 0; j < nn; j++) {
      WriteReal(writing,aa[j]);
      if (j) height=height+aa[j]*aa[j];
    }  /*of j*/
    WriteReal(writing,height); 
    fprintf(writing,"\n");   
  }  /*of i*/
  for (j=1; j<= nn; j++) WriteReal(writing,0);
  WriteReal(writing, 1);
  fprintf(writing,"\nend\nhull\n");
  fprintf(writing,"incidence\n");
  printf("file %s is generated.\n", outputfile);
_L99:;
  fclose(reading);
  fclose(writing);
}


