/* cddpivot.C:  Pivoting Procedures for cdd.C
   written by Komei Fukuda, fukuda@dma.epfl.ch
   Version 0.72a, April 25, 1995 
*/

/* cdd.c : C-Implementation of the double description method for
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
#include "setoper.h"  /* set operation library header (Ver. April 15,1995 or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
} /* end of extern "C" */


void CopyBmatrix(Bmatrix T, Bmatrix TCOPY)
{
  colrange j;

  for (j=0; j < nn; j++) {
    TCOPY[j] = T[j];
  }
}

void free_Bmatrix(Bmatrix T)
{
  colrange j;

  for (j = 0; j < nn; j++) {
    free(T[j]);
  }
}

void SetToIdentity(Bmatrix T)
{
  colrange j1, j2;
  myTYPE t_one=1;
  myTYPE t_zero=0;

  for (j1 = 1; j1 <= nn; j1++) {
    for (j2 = 1; j2 <= nn; j2++) {
      if (j1 == j2)
        T[j1 - 1][j2 - 1] = t_one;
      else
        T[j1 - 1][j2 - 1] = t_zero;
    }
  }
}

void SelectPivot1(Amatrix X, HyperplaneOrderType roworder,
    rowrange rowmax, rowset NopivotRow,
    colset NopivotCol, rowrange *r, colrange *s,
    boolean *selected)
/* Select a position (*r,*s) in the matrix X such that X[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  boolean stop;
  rowrange rtemp;
  rowset rowexcluded;

  stop = False;
  set_initialize(&rowexcluded,mm);
  set_copy(rowexcluded,NopivotRow);
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = False;
  do {
    SelectNextHyperplane(roworder, rowexcluded, &rtemp, &RecomputeRowOrder);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= nn && !*selected) {
        if (!set_member(*s,NopivotCol) && FABS(X[*r - 1][*s - 1]) > zero) {
          *selected = True;
          stop = True;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = True;
    }
  } while (!stop);
  set_free(&rowexcluded);
}

myTYPE TableauEntry(Amatrix X, Bmatrix T,
				rowrange r, colrange s)
/* Compute the (r,s) entry of X.T   */
{
  colrange j;
  myTYPE temp=0;
  boolean localdebug=False;
  
  temp=0;
  if (localdebug) cout << "r = " << r << "   s=" << s << "\n";
  for (j=0; j< nn; j++) {
    if (localdebug) cout << "temp = " << temp << "   j=" << j << "\n";
    if (localdebug) cout << "X[r-1] = " << X[r-1][j] << "   T[j][s-1]=" << T[j][s-1] << "\n";
    temp = temp + X[r-1][j] * T[j][s-1];
  }
  return temp;
}


void SelectPivot2(Amatrix X, Bmatrix T,
    HyperplaneOrderType roworder,
    rowrange rowmax, rowset NopivotRow,
    colset NopivotCol, rowrange *r, colrange *s,
    boolean *selected)
/* Select a position (*r,*s) in the matrix X.T such that (X.T)[*r][*s] is nonzero
   The choice is feasible, i.e., not on NopivotRow and NopivotCol, and
   best with respect to the specified roworder 
 */
{
  boolean stop;
  rowrange i,rtemp;
  rowset rowexcluded;
  myTYPE Xtemp=0;
  boolean localdebug=False;

  if (debug) localdebug=True;
  stop = False;
  set_initialize(&rowexcluded,mm);
  set_copy(rowexcluded,NopivotRow);
  if (localdebug) {
    switch (roworder) {

    case MinIndex:
      cout << "*SelectPivot2: MinIndex\n";
      break;

    case MaxIndex:
      cout << "*SelectPivot2: MaxIndex\n";
      break;

    case MinCutoff:
      cout << "*SelectPivot2: MinCutoff\n";
      break;

    case MaxCutoff:
      cout << "*SelectPivot2: MaxCutoff\n";
    break;

    case MixCutoff:
      cout <<  "*SelectPivot2: MixCutoff\n";
      break;

    case LexMin:
      cout <<  "*SelectPivot2: LexMin\n";
      break;

    case LexMax:
      cout <<  "*SelectPivot2: LexMax\n";
      break;

    case RandomRow:
      cout <<  "*SelectPivot2: Random,  Seed = " << rseed << "\n";
      break;

    case LineShelling:
      cout <<  "*SelectPivot2: LineShelling\n";
      break;
    }
    printf("select pivot2: rowexcluded=");
    set_write(rowexcluded);
  }
  for (rtemp=rowmax+1;rtemp<=mm;rtemp++) {
    set_addelem(rowexcluded,rtemp);   /* cannot pivot on any row > rmax */
  }
  *selected = False;
  do {
    i=1;rtemp=0;
    while (i<=mm && rtemp==0) {  /* EqualitySet vars have highest priorities */
      if (set_member(i,EqualitySet) && !set_member(i,rowexcluded)){
        if (localdebug) printf("marked set %ld chosen as a candidate\n",i);
        rtemp=i;
      }
      i++;
    }
    if (rtemp==0) SelectNextHyperplane(roworder, rowexcluded, &rtemp, &RecomputeRowOrder);
    if (rtemp>=1) {
      *r=rtemp;
      *s=1;
      while (*s <= nn && !*selected) {
        Xtemp=TableauEntry(X,T,*r,*s);
        if (!set_member(*s,NopivotCol) && FABS(Xtemp) > zero) {
          *selected = True;
          stop = True;
        } else {
          (*s)++;
        }
      }
      if (!*selected) {
        set_addelem(rowexcluded, rtemp);
      }
    }
    else {
      *r = 0;
      *s = 0;
      stop = True;
    }
  } while (!stop);
  set_free(&rowexcluded);
}

void GausianColumnPivot1(Amatrix X, rowrange r, colrange s,
				rowrange rowmax)
/* Make a column pivot operation in Amatrix X on position (r,s)  */
{
  long i, j;
  myTYPE Xtemp0=0, Xtemp=0;

  Xtemp0 = X[r - 1][s - 1];
  for (j = 0; j < nn; j++) {
    if (j + 1 != s) {
      Xtemp = X[r - 1][j];
      for (i = 0; i < rowmax; i++) {
        if (i + 1 != r)
        X[i][j] -= X[i][s - 1] * Xtemp / Xtemp0;
      }
      X[r - 1][j] = 0;
    }
  }
  for (i = 0; i < rowmax; i++) {
    if (i + 1 != r)
      X[i][s - 1] /= Xtemp0;
  }
  X[r - 1][s - 1] = 1;
}

void GausianColumnPivot2(Amatrix X, Bmatrix T,
				rowrange r, colrange s)
/* Update the Transformation matrix T with the pivot operation on (r,s) 
   This procedure performs a implicit pivot operation on the matrix X by
   updating the dual basis inverse  T.
 */
{
  long j, j1;
  myTYPE Xtemp0=0, Xtemp=0;
  myTYPE *vec1;

  vec1 = new myTYPE[nn];
  for (j=1; j<=nn; j++) vec1[j-1]=TableauEntry(X, T, r,j);
  Xtemp0 = vec1[s-1];
  for (j = 1; j <= nn; j++) {
    if (j != s) {
      Xtemp = vec1[j-1];
      for (j1 = 1; j1 <= nn; j1++)
        T[j1-1][j-1] -= T[j1-1][s - 1] * Xtemp / Xtemp0;
    }
  }
  for (j = 1; j <= nn; j++)
    T[j-1][s - 1] /= Xtemp0;
  delete[] vec1;
}

void ComputeRank(Amatrix A1, unsigned long *TargetRows, long *rank)
/* Compute the rank of the submatrix of a Amatrix A1 indexed by TargetRows.
   This procedure does not change the matrix A1.
 */
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset NoPivotRow;
  colset ColSelected;
  Bmatrix Btemp;   /* dual basis inverse */
  
  *rank = 0;
  stop = False;
  set_initialize(&NoPivotRow, mm);
  set_compl(NoPivotRow,TargetRows);
  set_initialize(&ColSelected, nn);
  InitializeBmatrix(Btemp);
  SetToIdentity(Btemp);
  if (debug) WriteBmatrix(cout,Btemp);
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, Btemp, MinIndex, mm, NoPivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(NoPivotRow, r);
        set_addelem(ColSelected, s);
        (*rank)++;
        GausianColumnPivot2(A1,Btemp, r, s);
        if (debug) {
          WriteBmatrix(cout,Btemp);
	  printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	}
      } else {
        stop=True;
      }
  } while (!stop);
  set_free(&NoPivotRow);
  set_free(&ColSelected);
  free_Bmatrix(Btemp);
}

void ComputeBInverse(Amatrix A1, long lastrow,
       Bmatrix InvA1, long *rank)
{
  boolean stop, chosen;
  rowrange r;
  colrange s;
  rowset RowSelected;
  colset ColSelected;

  *rank = 0;
  stop = False;
  SetToIdentity(InvA1);
  set_initialize(&RowSelected, mm);
  set_initialize(&ColSelected, nn);
  do {
    SelectPivot2(A1, InvA1, MinIndex, lastrow, RowSelected, ColSelected, &r, &s, &chosen);
    if (chosen) {
      (*rank)++;
      if (debug)
        printf("%3ldth pivot on%3ld, %3ld\n", *rank, r, s);
      GausianColumnPivot2(A1, InvA1, r, s);
      set_addelem(RowSelected, r);
      set_addelem(ColSelected, s);
    } else
      stop = True;
  } while (!stop);
  set_free(&RowSelected);
  set_free(&ColSelected);
}


void FindBasis(Amatrix A1, 
    HyperplaneOrderType roword, 
    rowset RowSelected,colindex ColInd,
    Bmatrix BasisInverse, long *rank)
{
  boolean stop, chosen;
  rowset NopivotRow;
  colset ColSelected;
  rowrange r;
  colrange j,s;

  *rank = 0;
  stop = False;
  for (j=0;j<=nn;j++) ColInd[j]=0;
  set_emptyset(RowSelected);
  set_initialize(&ColSelected, nn);
  set_initialize(&NopivotRow, mm);
  set_copy(NopivotRow,NonequalitySet);
  SetToIdentity(BasisInverse);
  if (debug) WriteBmatrix(cout,BasisInverse);
  if (DynamicWriteOn && !debug){
    cout << "*Initial set of rows:";
  }
  do {   /* Find a set of rows for a basis */
      SelectPivot2(A1, BasisInverse, roword, mm, NopivotRow, ColSelected, &r, &s, &chosen);
      if (debug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(NopivotRow, r);
        set_addelem(ColSelected, s);
        ColInd[s]=r;    /* ColInd[s] stores the corr. row index */
        (*rank)++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
        if (debug) {
          WriteBmatrix(cout,BasisInverse);
          WriteTableau(cout,A1,BasisInverse,NonzeroRHS),
	  printf("%3ldth row added to the initial set (%ldth elem)\n",  r, *rank);
	}
	if (DynamicWriteOn && !debug){
	   cout << " " << r;
	}
      } else {
        stop=True;
      }
      if (*rank==nn) stop = True;
  } while (!stop);
  if (DynamicWriteOn && !debug){
    cout << "\n";
  }
  set_free(&ColSelected);
}


void SelectCrissCrossPivot(Amatrix X, Bmatrix T, rowindex OV,
    long bflag[], rowrange objrow, colrange rhscol,
    rowrange *r, colrange *s,
    boolean *selected, LPStatusType *lps)
{
  boolean colselected=False, rowselected=False;
  rowrange i,k;
  myTYPE val=0;
  
  *selected=False;
  *lps=LPSundecided;
  while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
    for (k=1; k<=mm; k++) {
      i=OV[k];
      if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=TableauEntry(X,T,i,rhscol);
        if (val < -zero) {
          rowselected=True;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=TableauEntry(X,T,objrow,bflag[i]);
        if (val > zero) {
          colselected=True;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=Optimal;
      return;
    }
    else if (rowselected) {
      for (k=1; k<=mm; k++) {
        i=OV[k];
        if (bflag[i] >0) { /* i is nonbasic variable */
          val=TableauEntry(X,T,*r,bflag[i]);
          if (val > zero) {
            colselected=True;
            *s=bflag[i];
            *selected=True;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (k=1; k<=mm; k++) {
        i=OV[k];
        if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=TableauEntry(X,T,i,*s);
          if (val < -zero) {
            rowselected=True;
            *r=i;
            *selected=True;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=DualInconsistent;
    }
    else if (!colselected) {
      *lps=Inconsistent;
    }
  }
}

long InconsistencyDegree(Amatrix X, Bmatrix T, rowrange i,
    long bflag[], rowrange objrow, colrange rhscol)
/* This procedure returns the inconsistency degree
   of the basic row i
*/ 
{ colrange j;
  myTYPE val=0;
  long deg=0;

  if (bflag[i] >0 && TableauEntry(X,T,i,rhscol) < zero) { 
    /* evaluate the following only when i is basic and RHS is negative */
    deg=mm;
    for (j=1; j<=nn; j++) {
      if (j!=rhscol && TableauEntry(X,T,i,j) <= zero) deg++;
    }
  }
  return deg;
}

void OptimizeOrderVector(
    Amatrix X, Bmatrix T, rowindex OV, rowindex L,
    long bflag[], rowrange objrow, colrange rhscol)
/* This procedure modifies the order vector OV according
   to the current L-vector L so that the next Criss-Cross
   pivot improves L as much as possible 
*/
{
  rowrange i,k;
  colrange j;
  myTYPE val=0;
  Amatrix P;  // priority matrix
  
  for (i=1; i<=mm; i++) {
    P[i-1] = new myTYPE[2];
    for (j=1; j<=2; j++) P[i-1][j-1] = 0;
  }
  for (i=1; i<=mm; i++) {
    if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
      P[i-1][0]=-TableauEntry(X,T,i,rhscol);
    }
    else if (bflag[i] >0) { /* i is nonbasic variable */
      P[i-1][0]=TableauEntry(X,T,objrow,bflag[i]);
    }
  }
  QuickSort(OV, 1, mm, P, 2);
  for (i=1; i<=mm; i++) delete[] P[i-1];
}

void SelectCrissCrossDynamicPivot(
    Amatrix X, Bmatrix T, rowindex OV, rowindex L,
    long bflag[], rowrange objrow, colrange rhscol,
    rowrange *r, colrange *s,
    boolean *selected, LPStatusType *lps)
/* This procedure selects a pivot according to order vector OV 
   and modify OV so that L will increase as much as possible
*/ 
{
  boolean colselected=False, rowselected=False;
  rowrange i,k;
  myTYPE val=0;
  
  *selected=False;
  *lps=LPSundecided;
  while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
    for (k=1; k<=mm; k++) {
      i=OV[k];
      if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
        val=TableauEntry(X,T,i,rhscol);
        if (val < -zero) {
          rowselected=True;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=TableauEntry(X,T,objrow,bflag[i]);
        if (val > zero) {
          colselected=True;
          *s=bflag[i];
          break;
        }
      }
    }
    if  ((!rowselected) && (!colselected)) {
      *lps=Optimal;
      return;
    }
    else if (rowselected) {
     for (i=1; i<=mm; i++) {
       if (bflag[i] >0) { /* i is nonbasic variable */
          val=TableauEntry(X,T,*r,bflag[i]);
          if (val > zero) {
            colselected=True;
            *s=bflag[i];
            *selected=True;
            break;
          }
        }
      }
    }
    else if (colselected) {
      for (i=1; i<=mm; i++) {
        if (bflag[i]!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          val=TableauEntry(X,T,i,*s);
          if (val < -zero) {
            rowselected=True;
            *r=i;
            *selected=True;
            break;
          }
        }
      }
    }
    if (!rowselected) {
      *lps=DualInconsistent;
    }
    else if (!colselected) {
      *lps=Inconsistent;
    }
  }
}

void CrissCrossMinimize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
{
   colrange j;
   
   for (j=1; j<=nn; j++)
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   CrissCrossMaximize(f, f_log, A1,BasisInverse, OBJrow, RHScol, 
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter);
   *optvalue=-*optvalue;
   for (j=1; j<=nn; j++){
     dsol[j-1]=-dsol[j-1];
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   }
}

long OrderOf(rowrange p, rowindex OV)
{
  rowrange i=1;

  while (OV[i]!=p && i<=mm) i++;
  if (i>mm) return 0;else return i;
}
 

void UpdateLvector(rowindex L, rowrange p, rowindex OV)
{
  rowrange op, i=1;

  op = OrderOf(p,OV);  /* op is the OV order of a row index p */
  
  for (i=1; i<op; i++) L[OV[i]]=0;
  L[p] = 1;
}

void WriteLvector(ostream &f, rowindex L, rowindex OV)
{
  rowrange i;

/*
  f << " Order    :";
  for (i=1; i<= mm; i++){
    f.width(3); 
    f << i;
  }
*/
  f << " Variable :";  
  for (i=1; i<= mm; i++){
    if (OV[i]!=OBJrow){
      f.width(3); 
      f << OV[i];
    }
  }  
  f << "\n L vector :";
  for (i=1; i<= mm; i++){
    if (OV[i]!=OBJrow){
      f.width(3); 
      f << L[OV[i]];
    }
  }  
  f << "\n";
}

void CrissCrossMaximize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  boolean stop, chosen;
  long rank, ideg;
  rowrange i,r,entering,leaving,keyindex;
  colrange j,s;
  colset ColSelected;
  rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0;
  myTYPE adet=0; /* abs value of the determinant of a basis */
  boolean localdebug=False;
  rowindex L, OrderVec;

  L = new long[mm+1];
  OrderVec = new long[mm+1];
  for (i=0; i<=mm; i++) { L[i]=0; OrderVec[i]=i;} 
  if (debug) localdebug=True;
  if (BasisFlag==NULL || mlast!=mm){
     if (mlast!=mm) free(BasisFlag);   /* called previously with different mm */
     BasisFlag=(long *) calloc(mm+1, sizeof *BasisFlag);  
     /* initialize only for the first time or when a larger space is needed */
     mlast=mm;
  }
  *re=0; *se=0; *iter=0;
  rank = 0;
  stop = False;
  PreOrderedRun=False;
  InitializeBmatrix(InitialRays);
  adet=1.0;
  set_initialize(&Cobasis,mm);
  set_initialize(&Basis,mm);
  set_initialize(&RowSelected, mm);
  set_initialize(&ColSelected, nn);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  for (i=0; i<=mm; i++) BasisFlag[i]=0;
  for (j=0; j<=NMAX; j++) NBIndex[j]=0;
  for (i=1; i<=mm; i++) {
    set_addelem(Basis,i);
    BasisFlag[i]=-1;    /* basic variable has index -1 */
    if (EqualityIndex[i]==-1){
      if (DynamicWriteOn){ 
        cout << "*Warning: strict_inquality option for row " << i << " is ignored for maximization/minimization\n";
      }
    }
  }
  BasisFlag[OBJrow]= 0; /*  BasisFlag of the objective variable is 0, 
    different from other basic variables which have -1 */
  SetToIdentity(BasisInverse);
  if (localdebug) WriteBmatrix(cout,BasisInverse);
  do {   /* Find a LP basis */
      SelectPivot2(A1, BasisInverse, MinIndex
      , mm, RowSelected, ColSelected, &r, &s, &chosen);
      if (localdebug && chosen) printf("Procedure FindBasis: pivot on (r,s) =(%ld, %ld).\n", r, s);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        set_addelem(Cobasis, r);
        set_delelem(Basis,r);
        if (localdebug) {
          cout << "CC: find initial basis: nonbasis = ";   
          set_write(Cobasis); cout << "\n";
          cout << "CC: find initial basis: basis = ";
          set_write(Basis); cout << "\n";
        }
        BasisFlag[r]=s;   /* the nonbasic variable r corresponds to column s */
        NBIndex[s]=r;     /* the nonbasic variable on s column is r */
        if (localdebug) cout << "nonbasic variable " << r << " has index " << BasisFlag[r] << "\n";
        rank++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
        if (localdebug) {
          WriteBmatrix(cout,BasisInverse);
          WriteTableau(cout,A1,BasisInverse,Inequality);
	  cout << r << "th row added to the initial set (" << rank << "elem)\n";
	}
      } else {
        stop=True;
      }
      if (rank==nn-1) stop = True;
  } while (!stop);
  
  stop=False;
  if (OptimizeOrderOn) {
    OptimizeOrderVector(A1, BasisInverse, OrderVec, L, BasisFlag, OBJrow, RHScol);
  } else {
    ComputeRowOrderVector(OrderVec, HyperplaneOrder);
  }
  if (LogWriteOn) {
    WriteLvector(f_log,L,OrderVec);
  }
  if (ShowSignTableauOn){
    WriteSignTableau(cout,A1,BasisInverse,OrderVec, BasisFlag, OBJrow, RHScol);
    WriteLvector(cout,L,OrderVec); cout << "\n";
  }
  do {   /* Criss-Cross Method */
    SelectCrissCrossPivot(A1, BasisInverse, OrderVec, BasisFlag,
       OBJrow, RHScol, &r, &s, &chosen, LPS);
    if (chosen) {
      for (i=1; i<=mm; i++){
        ideg= InconsistencyDegree(A1, BasisInverse, i, BasisFlag, OBJrow, RHScol);
        if (ShowSignTableauOn && ideg > mm) cout << "inconsistency degree of row " << i << " is " << 
          ideg << "\n";
      }
      entering=NBIndex[s];
      leaving=r;
      if (localdebug) {
        cout<< "Procedure Criss-Cross: pivot on (r,s) = (" 
          << r << ", " << s << ")\n";
        cout<< "Procedure Criss-Cross: (leaving, entering) = (" 
          << leaving << ", " << entering << ")\n";
      }
      if (OrderOf(entering,OrderVec) > OrderOf(leaving,OrderVec)) keyindex = entering;
      else keyindex = leaving;
      UpdateLvector(L,keyindex,OrderVec); 
      if (LogWriteOn) WriteLvector(f_log,L,OrderVec);

      set_addelem(Cobasis, leaving);
      set_delelem(Cobasis, entering);
      set_delelem(Basis,leaving);
      set_addelem(Basis,entering);
      BasisFlag[leaving]=s;
      BasisFlag[entering]=-1;
      NBIndex[s]=leaving;
      if (localdebug) {
        cout << "nonbasis = "; set_write(Cobasis);
        cout << "\nbasis = "; set_write(Basis);
        cout << "\nnew nonbasic variable " << leaving << " has index " << BasisFlag[leaving] << "\n";
      }
      GausianColumnPivot2(A1,BasisInverse, r, s);
      (*iter)++;
      if (ShowSignTableauOn){ 
        cout << "Pivot on (" << leaving << ", " << entering << ")\n";
        WriteSignTableau(cout,A1,BasisInverse,OrderVec, BasisFlag, OBJrow, RHScol);
        WriteLvector(cout,L,OrderVec); cout << "\n";
      }
    } else {
      switch (*LPS){
        case Inconsistent: *re=r;
        case DualInconsistent: *se=s;
        default: break;
      }
      stop=True;
    }
  } while(!stop);
  if (debug) cout << "LP solved with " << *iter << " iterations.\n"; 
  switch (*LPS){
  case Optimal:
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,OBJrow,j);
      *optvalue=TableauEntry(A1,BasisInverse,OBJrow,RHScol);
      if (localdebug) cout << "dsol[" << NBIndex[j] << "]= " <<dsol[j-1] << "\n";
    }
    break;
  case Inconsistent:
    if (localdebug) printf("CrissCrossSolve: LP is inconsistent.\n");
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][RHScol-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,*re,j);
      if (localdebug) cout << "dsol " << NBIndex[j] << " " <<dsol[j-1];
    }
    break;
  case DualInconsistent:
    for (j=1;j<=nn; j++) {
      sol[j-1]=BasisInverse[j-1][*se-1];
      dsol[j-1]=-TableauEntry(A1,BasisInverse,OBJrow,j);
      if (localdebug) cout << "dsol " << NBIndex[j] << " " <<dsol[j-1];
    }
    if (localdebug) cout << "CrissCrossSolve: LP is dual inconsistent.\n";
    break;

  default:break;
  }

  delete[] L;
  delete[] OrderVec;
  set_free(&ColSelected);
  set_free(&RowSelected);
  set_free(&Basis);
  set_free(&Cobasis);
}

void InitializeBmatrix(Bmatrix T)
{
  colrange j;

  for (j = 0; j < nn; j++) {
    T[j]= new myTYPE[nn];
  }
}

void ReduceAA(rowset ChosenRow, colset ChosenCol)
/* Set the matrix AA to be the submatrix of AA with chosen 
   rows & columns, and change mm and nn accordingly 
*/
{
  long i,j,inew,jnew,mnew=0,nnew=0;
  Amatrix Acopy;

  mnew=set_card(ChosenRow);
  nnew=set_card(ChosenCol);
  for (i=0; i<mnew; i++){
    Acopy[i]= new myTYPE[nnew];
  }
  inew=0;
  for (i=1; i <= mm; i++) {
    if (set_member(i,ChosenRow)) {
      inew++;
      jnew=0;
      for (j=1; j <= nn; j++) {
        if (set_member(j, ChosenCol)){
          jnew++;
          Acopy[inew-1][jnew-1]=AA[i-1][j-1];
          if (debug) WriteNumber(cout, AA[i-1][j-1]);
        }
      }
      if (debug) cout << "\n";
    }
  }
  for (i=1;i<=mm;i++) {
    if (i<=mnew) 
      set_addelem(ChosenRow,i);
    else
      set_delelem(ChosenRow,i);
  }
  for (j=1;j<=nn;j++) {
    if (j<=nnew) 
      set_addelem(ChosenCol,j);
    else
      set_delelem(ChosenCol,j);
  }
  if (debug) {
    cout << "new row indices:";set_write(ChosenRow);
    cout << "new col indices:";set_write(ChosenCol);
  }
  for (i=0;i<mm;i++){
    delete[] AA[i];
  }
  for (i=0;i<mnew;i++){
    AA[i]=Acopy[i];
  }
  mm=mnew;  nn=nnew;
}

void DualizeAA(Bmatrix T)
/* Set the matrix AA to be the transpose of the matrix [-AA.T | I],
   and change mm and nn accordingly 
*/
{
  long i,j,mnew,nnew;
  Amatrix Acopy;

  mnew=nn+mm;
  nnew=mm;
  for (i=0; i<mm; i++){
    Acopy[i]= new myTYPE[nn];
  }
  if (mnew > MMAX) {
    printf("MMAX  is too small for ray computation. MMAX must be >= %ld.\n",mnew);
    Error = DimensionTooLarge;
    goto _L99;
  }
  if (nnew > NMAX) {
    printf("NMAX  is too small for ray computation. NMAX must be >= %ld.\n",nnew);
    Error = DimensionTooLarge;
    goto _L99;
  }
  for (i=1;i<=mm;i++){
    for (j=1;j<=nn;j++){
      Acopy[i-1][j-1]=TableauEntry(AA,T,i,j);
    }
  }
  for (i=0;i<mm;i++){
    delete[] AA[i];
  }
  for (i=0; i<mnew; i++){
    AA[i]= new myTYPE[nnew];
  }
  for (i=1;i<=nn;i++){
    for (j=1;j<=mm;j++){
      AA[i-1][j-1]=-Acopy[j-1][i-1];
    }
  }
  for (i=1;i<=mm;i++){
    for (j=1;j<=mm;j++){
      if (i==j) AA[nn+i-1][j-1]=1;
      else AA[nn+i-1][j-1]=0;
    }
  }
  for (i=0; i<mm; i++){
    delete[] Acopy[i];
  }
  mm=mnew;  nn=nnew;
  _L99:;
}

void EnlargeAAforInteriorFinding(void)
/* Add an extra column with all minus ones to the matrix AA, 
   add an objective row with (0,...,0,1), and 
   rows & columns, and change mm and nn accordingly 
*/
{
  long i,j,mnew=0,nnew=0;
  Amatrix Acopy;

  mnew=mm+1;
  nnew=nn+1;
  for (i=0; i<mnew; i++){
    Acopy[i]= new myTYPE[nnew];
  }
  for (i=1; i <= mm; i++) {
    for (j=1; j <= nn; j++) {
      Acopy[i-1][j-1]=AA[i-1][j-1];
      if (debug) WriteNumber(cout, AA[i-1][j-1]);
    }
    if (debug) cout << "\n";
  }
  for (i=1;i<=mm;i++) {
    Acopy[i-1][nn]=-1.0;  /* new column with all minus one's */
  }
  for (j=1;j<=nn;j++) {
    Acopy[mm][j-1]=0.0;  /* new row with (0,...,0,1) */
  }
  Acopy[mm][nn]=1.0;  /* new row with (0,...,0,1) */
  for (i=0;i<mm;i++){
    delete[] AA[i];
  }
  for (i=0;i<mnew;i++){
    AA[i]=Acopy[i];
  }
  mm=mnew;  nn=nnew;
}

void RecoverAAafterInteriorFinding(void)
/* Remove the extra column with all minus ones, 
   and remove the objective row with (0,...,0,1), and 
   rows & columns, and change mm and nn accordingly 
*/
{
  long i,j,mnew=0,nnew=0;
  Amatrix Acopy;

  mnew=mm-1;
  nnew=nn-1;
  for (i=0; i<mnew; i++){
    Acopy[i]= new myTYPE[nnew];
  }
  for (i=1; i <= mm; i++) {
    for (j=1; j <= nn; j++) {
      Acopy[i-1][j-1]=AA[i-1][j-1];
      if (debug) WriteNumber(cout, AA[i-1][j-1]);
    }
    if (debug) cout << "\n";
  }
  for (i=0;i<mm;i++){
    delete[] AA[i];
  }
  for (i=0;i<mnew;i++){
    AA[i]=Acopy[i];
  }
  mm=mnew;  nn=nnew;
}

boolean RowEquivalent_Q(Arow a1, Arow a2, colrange n)
/* Check if the row a1 is positive multiple of a2 */
{
  colrange j,j1,j2;
  boolean equivalent=True;
  myTYPE scaler=1, purezero=0;

  j1=1; while (FABS(a1[j1-1]) <=zero && j1<=n) j1++;
  j2=1; while (FABS(a2[j2-1]) <=zero && j2<=n) j2++;
  if (j1!=j2) equivalent=False;
  else {
    if (j1>n) equivalent=True;  /* both are zero vectors */
    else if (a1[j1-1]*a2[j1-1]<purezero) equivalent=False;
      else {
        j=1; scaler=a1[j1-1]/a2[j1-1];
        if (debug) cout << "Checking the equivalence with scaler ="
          << scaler << "\n"; 
        while (equivalent && j<=n) {
          if (FABS(a1[j-1] - a2[j-1]*scaler) > zero ) equivalent = False;
          j++; 
        }
      }
  }
  return equivalent;
}

void FindRowEquivalenceClasses(long *classno, rowindex rindex)
/* Find the classes of equivalent (inequality) rows
   and return the number *classno of classes and
   rowindex rindex[i] which indicates the class number
   it belongs to.  *classno == mm if every row constitutes
   a singleton class (i.e. there is no independence). 
*/
{
  rowrange i,k;
  long rclass=0;
  boolean localdebug=False;

  for (i=1; i<=mm; i++) rindex[i]=0;
  for (i=1; i<=mm; i++)
  {
    if (rindex[i]==0) {
      rclass++; rindex[i]=rclass;
      if (localdebug) cout << "row " << i << 
        " is in the new class " << rclass << ".\n";
      for (k=i+1; k<=mm; k++) {
        if (RowEquivalent_Q(AA[i-1], AA[k-1], nn)) {
          rindex[k]=rclass;
          if (localdebug) cout << "   " << k << " is in the same class.\n";
        }
      } 
    }
  }
  *classno = rclass; 
}

 /* end of cddpivot.C */  
