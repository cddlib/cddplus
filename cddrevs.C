/* cddrevs.C:  Reverse Search Procedures for cdd.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <fstream>
#include <string>
using namespace std;

#include "cddtype.h"
#include "cddrevs.h"

extern "C" {
#include "setoper.h"  /* set operation library header (Ver. May 14,1995 or later) */
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
  long j=0,count=1,fcount=0;
  boolean completed=False;

  cout << "\nReverse search starts with #1 object: ";
  s.fwrite(cout); 
  if (CondensedListOn) {
    wf << "begin\n" << "  *****  " << delta << "  tope_condensed\n";
  }else{
    wf << "begin\n" << "  *****  " << delta << "  tope\n";
  }
  s.fwrite(wf); 
  while (!completed)
  {
    while (j<delta)
    {
      j=j+1;
      w=Adj(v,j);
      if (w!=v && f(w)==v)
      {
        count++;
        if (CondensedListOn){
          cout << "\n #" << count << " r " << j << " f " << fcount; 
          wf << "\n r " << j << " f " << fcount; 
        } else { 
          cout << "\nA new object #" << count << " found: "; w.fwrite(cout);
          wf << "\n"; w.fwrite(wf);
        }
        v=w;   j=0;
        fcount=0;
      }
    }
    if (!(v==s))
    {
      u=v;  v=f(v);
      j=NeighbourIndex(u,v);
      fcount++; 
      if (j==delta) completed=True;
    }
  }
  wf << "\nend";
  cout << "\nNumber ***** of topes = " << count << "\n";
  wf << "\nNumber ***** of topes = " << count << "\n";
}

// Facet recognition programs using Linear Programming

boolean Facet_Q(topeOBJECT tope, rowrange ii)
{
  static colindex nbindex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  static Arow LPdsol;  /*  LP solution and the dual solution (basic var only) */
  static colrange nlast=0;

  if (nlast!=nn){
    if (nlast>0){
      delete[] LPdsol;
    }
    LPdsol = new myTYPE[nn];
    nlast=nn; 
  }
  return Facet_Q2(tope,ii,nbindex,LPdsol);
}

boolean Facet_Q2(topeOBJECT tope, rowrange ii, colindex NBIndex, Arow LPdsol)
{
  /* Before calling this, LPdsol must be initialized with LPdsol = new myTYPE[nn].
     When ii is detected to be non-facet,
     NBIndex returns the nonbasic variables at the evidence solution.
     LPdsol returns the evidence dual solution.
  */
  static Arow LPsol;  /*  LP solution and the dual solution (basic var only) */
  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  colrange s=0;
  myTYPE ov=0, tempRHS=0, purezero=0, pureone=1;  /* ov = LP optimum value */
  long LPiter,testi, i, j;
  boolean answer=True,localdebug=False;
  static colrange nlast=0;
  static Bmatrix BInv;
  static boolean firstcall=True, UsePrevBasis;
  static ConversionType ConversionSave;

  if (nlast!=nn){
    if (nlast>0){
      delete[] LPsol; 
      free_Bmatrix(BInv);
    }
    InitializeBmatrix(BInv);
    LPsol = new myTYPE[nn];
    nlast=nn;
    firstcall=True;
  }
//  if (firstcall || Conversion==TopeListing) 
//    UsePrevBasis=False; else UsePrevBasis=True;
  UsePrevBasis=False;
  ConversionSave=Conversion;
  Conversion=LPmax;
  RHScol=1;
  mm=mm+1;
  OBJrow=mm;
  AA[OBJrow-1]=new myTYPE[nn];
  if (localdebug) cout << "Facet_Q:  create an exra row " << OBJrow << "\n";
  for (i=1; i<=mm; i++) {
    if (tope[i]<0) {
      if (debug) cout << "reversing the signs of " << i << "th inequality\n";
      for (s=0,j=1; j<=nn; j++){
        AA[i-1][j-1]=-AA[i-1][j-1];
      }
    }
  } 
  tempRHS=AA[ii-1][0];
  for (s=0,j=1; j<=nn; j++){
    AA[OBJrow-1][j-1]=-AA[ii-1][j-1];
    if (!firstcall && NBIndex[j]==ii){
      s=j;
      if (localdebug) cout << "Row "<< ii << " is nonbasic" << s << "\n";
    }
  }
  AA[OBJrow-1][0]=purezero;
  AA[ii-1][0]=tempRHS+pureone;   /* relax the ii-th inequality by a large number*/
  if (s>0) GausianColumnPivot2(AA,BInv, ii, s);
  DualSimplexMaximize(cout, cout, AA, BInv, OBJrow, RHScol, UsePrevBasis,
    &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  if (LPStatus!=Optimal){
    if (DynamicWriteOn) cout << "The Dual Simplex failed.  Run the Criss-Cross method.\n";
    CrissCrossMaximize(cout, cout, AA, BInv, OBJrow, RHScol, UsePrevBasis,
    &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  }
  if (localdebug) cout << ii << "-th LP solved with objective value =" << ov << 
    " RHS value = " << tempRHS << "  iter= " << LPiter << "\n";
  if ((ov - tempRHS) > zero) 
  {
    answer=True;
    if (localdebug) cout << ii << "-th inequality determines a facet.\n";
  }
  else {
    answer=False;
    if (localdebug) cout << ii << "-th inequality does not determine a facet.\n";
  }
  AA[ii-1][0]=tempRHS;   /* restore the original RHS */
  for (s=0,j=1; j<=nn; j++){
    if (NBIndex[j]==ii){
      s=j;
    }
  }
  if (s>0){
    if (localdebug) cout << "Row "<< ii << " is nonbasic: basisinv updated " << s << "\n";
    GausianColumnPivot2(AA,BInv, ii, s);
  }
  delete[] AA[OBJrow-1];
  if (localdebug) cout << "Facet_Q:  delete the exra row " << OBJrow << "\n";
  mm=mm-1;
  for (i=1; i<=mm; i++) {
    if (tope[i]<0) {
      for (j=1; j<=nn; j++)  AA[i-1][j-1]=-AA[i-1][j-1]; /* restore the original data */
    }
  }
  firstcall=False;
  Conversion=ConversionSave;
  return answer;
}


void FacetandVertexListMain(ostream &f, ostream &f_log)
{ 
  rowrange i;
  colrange j;
  rowset subrows; /* rows which define a facet inequality */
  colset allcols;  /* allcols: all column indices */
  rowindex rowequiv;
  rowrange classno;
  topeOBJECT Tope(mm);
  Arow LPdsol,center;
  static colindex NBIndex;

  LPdsol = new myTYPE[nn];
  center = new myTYPE[nn];
  WriteProgramDescription(f);
  (f) << "*Input File:" << inputfile << "(" << minput << "x" << ninput << ")\n";
  WriteRunningMode(f); WriteRunningMode(f_log);
  if (Conversion==VertexListing){
    ShiftPointsAroundOrigin(f,f_log, center); 
    /* Shifting the points so that the origin will be in the relative interior of
       their convex hull
    */
  }
  (f) << "* `e` means essential and `r` means redundant.\n";
  if (DynamicWriteOn) cout << "* `e` means essential and `r` means redundant.\n";
  rowequiv = new long[mm+1];
//  FindRowEquivalenceClasses(&classno, rowequiv);
//  if (classno<mm) {
//    cout << "*There are multiple equivalent rows!!!\n";
//    cout << "*You have to remove duplicates before listing facets. \n";
//    (f) << "*There are multiple equivalent rows!!!\n";
//    (f) << "*You have to remove duplicates before listing facets. \n";
//    WriteRowEquivalence(cout, classno, rowequiv);
//    WriteRowEquivalence(f, classno, rowequiv);
//    goto _L99;
//  }

  time(&starttime); 
  set_initialize(&subrows,mm);
  set_initialize(&allcols,nn); 
  for (j=1;j<=nn;j++) set_addelem(allcols,j);
  if (Inequality==ZeroRHS){
    printf("Sorry, facet/vertex listing is not implemented for RHS==0.\n");
    goto _L99;
  }
  (f) << "begin\n";
  if (DynamicWriteOn) (cout) <<"begin\n";
  for (i=1; i<=mm; i++){
    if (Facet_Q2(Tope, i, NBIndex, LPdsol)) {
      if (DynamicWriteOn) cout << i << " e:";
      (f) << i << " e:";
      set_addelem(subrows,i);
    }
    else {
      (f) <<  i << " r:";
      if (DynamicWriteOn) (cout) << i << " r:";
    }
    for (j=1; j<nn; j++){
      (f) << " " << NBIndex[j+1];
      if (LogWriteOn){
        (f) <<"(";  WriteNumber(f,LPdsol[j]); (f) << ")";
      }
    }
    (f) << "\n";
    if (DynamicWriteOn){
      for (j=1; j<nn; j++){
        (cout) << " " << NBIndex[j+1];
        if (LogWriteOn){
          (cout) <<"(";  WriteNumber(cout,LPdsol[j]); (cout) << ")";
        }
      }
      (cout) << "\n";
    }
  }
  (f) << "end\n";
  if (DynamicWriteOn) (cout) <<"end\n";
  time(&endtime);
//  (f) << "* Here is a minimal system representing the same polyhedral set as the input.\n";
//  WriteSubMatrixOfAA(f,subrows,allcols,Inequality);
  WriteTimes(f); WriteTimes(f_log); WriteTimes(cout);
//  set_free(&subrows);
//  set_free(&allcols); 
_L99:;
//  delete[] rowequiv;
//  delete[] LPdsol;
//  delete[] center;
}

void FacetandVertexExternalListMain(ostream &f, ostream &f_log)
{ 
  rowrange i,mmxtn;
  colrange j,nnxtn;
  rowset subrows; /* rows which define a facet inequality */
  colset allcols;  /* allcols: all column indices */
  rowindex rowequiv;
  rowrange classno;
  topeOBJECT Tope(mm);
  Arow LPdsol,center;
  static colindex NBIndex;
  string xtnnumbertype,command;
  myRational rvalue=0;
  myTYPE value=0;
  boolean found,localdebug=False;
  char ch;

  SetReadFileName(xtnfile,'x',"external");
  ifstream reading_xtn(xtnfile);

  if (reading_xtn.is_open()) {
    found=False;
    while (!found)
    {
      if (!reading_xtn.eof()) {
        reading_xtn >> command;
  	if (command=="begin") {
          found=True;
  	}
      }
      else {
  	Error=ImproperInputFormat;
  	goto _L99;
      }
    }
    reading_xtn >> mmxtn;
    reading_xtn >> nnxtn;
    reading_xtn >> xtnnumbertype;

    LPdsol = new myTYPE[nn];
    center = new myTYPE[nn];
    WriteProgramDescription(f);
    (f) << "*Essential File:" << inputfile << "(" << minput << "x" << ninput << ")\n";
    (f) << "*Input File:" << xtnfile << "(" << mmxtn << "x" << nnxtn << ")\n";
    WriteRunningMode(f); WriteRunningMode(f_log);
    if (Conversion==VertexListingExternal){
      ShiftPointsAroundOrigin(f,f_log, center); 
      /* Shifting the points so that the origin will be in the relative interior of
         their convex hull
      */
    }
    (f) << "* `e` means essential and `r` means redundant.\n";

    /* Extrarow to store each line from the external file */
    mm = minput + 1;
    AA[mm-1]= new myTYPE[ninput];

 
    time(&starttime); 
    set_initialize(&subrows,mm);
    set_initialize(&allcols,nn); 
    for (j=1;j<=nn;j++) set_addelem(allcols,j);
    if (Inequality==ZeroRHS){
      printf("Sorry, facet/vertex listing is not implemented for RHS==0.\n");
      goto _L99;
    }
    (f) << "begin\n";
    if (DynamicWriteOn) (cout) <<"begin\n";
    for (i=1; i<=mmxtn; i++){
      for (j = 1; j <= nn; j++) {
        if (xtnnumbertype=="rational" && OutputNumberString=="real"){
          reading_xtn >> rvalue;
          value=myTYPE(rvalue);
        } else {
          reading_xtn >> value;
        }
      	AA[mm-1][j - 1] = value;
	if (localdebug){
           if (xtnnumbertype=="rational" && OutputNumberString=="real")
             cout << "a(" << i << "," << j << ") = " << value << " ("<< rvalue << ")\n";
           else
             cout << "a(" << i << "," << j << ") = " << value  << "\n";
         }
      }  /*of j*/
      while (reading_xtn.get(ch) && ch != '\n') {
        if (localdebug) cout << ch;
      }

      if (Conversion==VertexListingExternal){
        /* Each point must be shifted w.r.t the relative interior of
           their convex hull */
        for (j=2; j<=nn; j++){AA[mm-1][j-1]-=center[j-1];}
        if (localdebug){
          for (j=1; j<=nn; j++){cout << " " << AA[mm-1][j-1];}
          cout << "  " << i << "th point Shifted.\n";
        }
      }
      if (Facet_Q2(Tope, mm, NBIndex, LPdsol)) {
        if (DynamicWriteOn) cout << i << " e:";
        (f) << i << " e:";
        set_addelem(subrows,i);
      }
      else {
        (f) <<  i << " r:";
        if (DynamicWriteOn) (cout) << i << " r:";
      }
      long poscomp_count=0,rindex=0; 
      for (j=1; j<nn; j++){
        (f) << " " << NBIndex[j+1];
        if (LPdsol[j]>zero) {poscomp_count++; rindex=NBIndex[j+1];};
        if (LogWriteOn){
          (f) <<"(";  WriteNumber(f,LPdsol[j]); (f) << ")";
        }
      }
      if (poscomp_count==1 && RowEquivalent_Q(AA[mm-1],AA[rindex-1], nn)) {
        (f) << " =#" << rindex;
      } else {poscomp_count=0;}
      (f) << "\n";
      if (DynamicWriteOn){
        for (j=1; j<nn; j++){
          (cout) << " " << NBIndex[j+1];
          if (LogWriteOn){
            (cout) <<"(";  WriteNumber(cout,LPdsol[j]); (cout) << ")";
          }
        }
      if (poscomp_count==1) cout << " =#" << rindex;
      (cout) << "\n";
      }
    }
    (f) << "end\n";
    if (DynamicWriteOn) (cout) <<"end\n";
    time(&endtime);
    WriteTimes(f); WriteTimes(f_log); WriteTimes(cout);
  _L99:;
    reading_xtn.close();
  } else {
    Error=FileNotFound;
    WriteErrorMessages(cout);
    WriteErrorMessages(f_log);
  }
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
    cout << "*You have to remove duplicates before listing topes. \n";
    (f) << "*There are multiple equivalent rows!!!\n";
    (f) << "*You have to remove duplicates before listing topes. \n";
    WriteRowEquivalence(cout, classno, rowequiv);
    WriteRowEquivalence(f, classno, rowequiv);
    goto _L99;
  }

  time(&starttime);
  set_initialize(&subrows,mm);
  set_initialize(&allcols,nn); 
  for (j=1;j<=nn;j++) set_addelem(allcols,j);
  if (Inequality==ZeroRHS){
    printf("Sorry, tope listing is not implemented for RHS==0.\n");
  } else {
    ReverseSearch(f, Tope,mm);
  }
  time(&endtime);
  WriteTimes(f); WriteTimes(f_log); WriteTimes(cout);
  set_free(&subrows);
  set_free(&allcols); 
_L99:;
  delete[] rowequiv;
}

// end of cddrevs.C

