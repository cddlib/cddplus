/* cddrevs.C:  Reverse Search Procedures for cdd.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.72, April 16, 1995
*/

/* cdd.C : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <fstream.h>
#include <strclass.h>
#include "cddtype.h"
#include "cddrevs.h"

extern "C" {
#include "setoper.h"  /* set operation library header (Ver.  April 15,1995 or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
} /* end of extern "C" */

topeOBJECT::topeOBJECT(const topeOBJECT& tope)
{
  long j;
  dim=tope.dim;
  sv = new int[dim]; 
  for (j=1;j<=dim;j++) sv[j-1]=tope.sv[j-1];
}

void topeOBJECT::operator=(const topeOBJECT& tope)
{
  long j;
  delete[] sv;
  dim = tope.dim;
  sv = new int[dim];
  for (j=1;j<=dim;j++) sv[j-1]=tope.sv[j-1];
}

int topeOBJECT::operator[](long j) // return the i-th component of tope
{
  if (j>=1 && j<=dim) return sv[j-1];else return 0; 
}

topeOBJECT operator-(topeOBJECT tope, long j)
// reversing the sign of sv[j-1]
{
  topeOBJECT t=tope;
  if (j>=1 && j<=t.dim) t.sv[j-1]=-t.sv[j-1];
  return t;
}

void topeOBJECT::fwrite(ostream& f)
{
  long j;
  for (j=1; j<=this->dim; j++) 
  { 
    if (this->sv[j-1] > 0) f << " +";
    if (this->sv[j-1] < 0) f << " -";
    if (this->sv[j-1] ==0) f << " 0"; 
  }
  f << "\n"; 
}

int operator==(const topeOBJECT &t1, const topeOBJECT &t2)
{
  long j=1;
  int equal=1;
  if (t1.dim != t2.dim) return 0;
  else
  {
    while (equal && j<=t1.dim)
    {
      if (t1.sv[j-1] != t2.sv[j-1]) equal=0;
      j++;
    }
  }
  return equal; 
}

int operator!=(const topeOBJECT &t1, const topeOBJECT &t2)
{
  if (t1==t2) return 0;
  else return 1; 
}


topeOBJECT f(topeOBJECT v)
{
  long i=1, nexti=0; 
 
  while (nexti==0 && i<=v.dim)
  {
    if (v[i]<0 && Facet_Q(v,i)) nexti=i;
    i++;
  }
  if (nexti==0) return v;
  else return v-nexti;
}

long NeighbourIndex(topeOBJECT u, topeOBJECT v)
{
  long i=1,nindex=0;

  while ( nindex==0 && i <= v.dim)
  {
    if (u==v-i) nindex = i;
    i++;
  }
  return nindex;
}

topeOBJECT Adj(topeOBJECT v, long i)
{
  if (i<=0 || i>v.dim || v[i] <= 0) return v;
  else if (Facet_Q(v, i)) return v-i;
    else return v;
}

void ReverseSearch(ostream &wf, topeOBJECT s, long delta)
{
  topeOBJECT v=s,u=s,w=s;
  long j=0,count=1;
  boolean completed=False;

  cout << "Reverse search starts with #1 object: "; s.fwrite(cout);
  (wf) << "begin\n";
  s.fwrite(wf); 
  while (!completed)
  {
    while (j<delta)
    {
      j=j+1;
      if (debug) cout << "checking " << j << "th neighbour.\n";
      w=Adj(v,j);
      if (debug) {cout << "candidate w is "; w.fwrite(cout);}
      if (w!=v && f(w)==v)
      {
        v=w;   j=0;
        count=count+1; 
        cout << "A new object #" << count << " found: ";
        v.fwrite(cout);
        v.fwrite(wf); 
      }
    }
    if (!(v==s))
    {
      u=v;  v=f(v);
      j=NeighbourIndex(u,v);
      if (j==delta) completed=True;
    }
  }
  (wf) << "end\n";
}

// Facet recognition programs using LPmax (Criss-Cross) code

boolean Facet_Q(topeOBJECT tope, rowrange ii)
{
  colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  myTYPE ov=0, tempRHS=0;  /* LP optimum value */
  long LPiter,testi, i, j;
  boolean answer=True;

  LPsol = new myTYPE[nn];
  LPdsol = new myTYPE[nn]; 
  Conversion=LPmax;
  time(&starttime);
  RHScol=1;
  mm=mm+1;
  OBJrow=mm;
  AA[OBJrow-1]=new myTYPE[nn];
  OBJrow=mm;
  for (i=1; i<=mm; i++) {
    if (tope[i]<0) {
      if (debug) cout << "reversing the signs of " << i << "th inequality\n";
      for (j=1; j<=nn; j++) AA[i-1][j-1]=-AA[i-1][j-1];
    }
  } 
  tempRHS=AA[ii-1][0];
  for (j=1; j<=nn; j++) AA[OBJrow-1][j-1]=-AA[ii-1][j-1]; 
  AA[OBJrow-1][0]=0;
  AA[ii-1][0]=tempRHS+1;   /* relax the ii-th inequality by +1 */
  CrissCrossMaximize(cout, cout, AA, InitialRays, OBJrow, RHScol, 
    &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  if (debug) cout << ii << "-th LP solved with objective value =" << ov << 
    " RHS value = " << tempRHS << "  iter= " << LPiter << "\n";
  if (ov > tempRHS) 
  {
    answer=True;
    if (debug) cout << ii << "-th inequality determines a facet.\n";
  }
  else {
    answer=False;
    if (debug) cout << ii << "-th inequality does not determine a facet.\n";
  }
  AA[ii-1][0]=tempRHS;   /* restore the original RHS */
  delete[] AA[OBJrow-1];
  mm=mm-1;
  for (i=1; i<=mm; i++) {
    if (tope[i]<0) {
      for (j=1; j<=nn; j++)  AA[i-1][j-1]=-AA[i-1][j-1]; /* restore the original data */
    }
  }
  return answer;
}

void FacetListMain(ostream &f, ostream &f_log)
{ 
  rowrange i;
  colrange j;
  rowset subrows; /* rows which define a facet inequality */
  colset allcols;  /* allcols: all column indices */
  rowindex rowequiv;
  rowrange classno;
  topeOBJECT Tope(mm);

  WriteProgramDescription(f);
  (f) << "*Input File:" << inputfile << "(" << minput << "x" << ninput << ")\n";
  WriteRunningMode(f); WriteRunningMode(f_log);
  rowequiv = new long[mm+1];
  FindRowEquivalenceClasses(&classno, rowequiv);
  if (classno<mm) {
    cout << "*There are multiple equivalent rows!!!\n";
    cout << "*You have to remove redundancies before listing facets. \n";
    (f) << "*There are multiple equivalent rows!!!\n";
    (f) << "*You have to remove redundancies before listing facets. \n";
    WriteRowEquivalence(cout, classno, rowequiv);
    WriteRowEquivalence(f, classno, rowequiv);
    goto _L99;
  }
 
  set_initialize(&subrows,mm);
  set_initialize(&allcols,nn); 
  for (j=1;j<=nn;j++) set_addelem(allcols,j);
  if (Inequality==ZeroRHS){
    printf("Sorry, facet listing is not implemented for RHS==0.\n");
    goto _L99;
  }
  for (i=1; i<=mm; i++){
    if (Facet_Q(Tope, i)) {
      if (DynamicWriteOn) cout << "row "<< i << " determines a facet.\n";
      (f) << "row " << i << " determines a facet.\n";
      set_addelem(subrows,i);
    }
    else {
      if (DynamicWriteOn) cout << "row " << i << " does not determine a facet.\n";
      (f) << "row " << i << " does not determine a facet.\n";
    }
  }
  (f) << "* Here is a minimal system representing the same polyhedral set as the input.\n";
  WriteSubMatrixOfAA(f,subrows,allcols,Inequality);
  set_free(&subrows);
  set_free(&allcols); 
_L99:;
  delete[] rowequiv;
}

void TopeListMain(ostream &f, ostream &f_log)
{ 
  rowrange i;
  colrange j;
  rowset subrows; /* rows which define a facet inequality */
  colset allcols;  /* allcols: all column indices */
  rowindex rowequiv;
  rowrange classno;
  topeOBJECT Tope(mm);

  WriteProgramDescription(f);
  (f) << "*Input File:" << inputfile << "(" << minput << "x" << ninput << ")\n";
  WriteRunningMode(f); WriteRunningMode(f_log);
  rowequiv = new long[mm+1];
  FindRowEquivalenceClasses(&classno, rowequiv);
  if (classno<mm) {
    cout << "*There are multiple equivalent rows!!!\n";
    cout << "*You have to remove redundancies before listing topes. \n";
    (f) << "*There are multiple equivalent rows!!!\n";
    (f) << "*You have to remove redundancies before listing topes. \n";
    WriteRowEquivalence(cout, classno, rowequiv);
    WriteRowEquivalence(f, classno, rowequiv);
    goto _L99;
  }

  cout << "the initial tope = "; Tope.fwrite(cout);  

  set_initialize(&subrows,mm);
  set_initialize(&allcols,nn); 
  for (j=1;j<=nn;j++) set_addelem(allcols,j);
  if (Inequality==ZeroRHS){
    printf("Sorry, tope listing is not implemented for RHS==0.\n");
  } else {
    ReverseSearch(f, Tope,mm);
  }
  set_free(&subrows);
  set_free(&allcols); 
_L99:;
  delete[] rowequiv;
}

// end of cddrevs.C

