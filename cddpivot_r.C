/* cddpivot.C:  Pivoting Procedures for cdd.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.c : C-Implementation of the double description method for
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
    delete[] T[j];
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

void SelectPivot1(SignAmatrix X, HyperplaneOrderType roworder,
    rowrange rowmax, rowset NopivotRow,
    colset NopivotCol, rowrange *r, colrange *s,
    boolean *selected)
/* Select a position (*r,*s) in the sign matrix X such that X[*r-1][*s-1] is nonzero
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
        if (!set_member(*s,NopivotCol) && X[*r - 1][*s - 1]!=0 ) {
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

void SetSignAmatrix(SignAmatrix SA, Amatrix X, Bmatrix T)
{
  rowrange i;
  colrange j;

  for (i=1; i<=mm; i++)
    for (j=1; j<=nn; j++)
      SA[i-1][j-1]=myTYPE2sign(TableauEntry(X, T, i, j));
}

myTYPE TableauEntry(Amatrix X, Bmatrix T, rowrange r, colrange s)
/* Compute the (r,s) entry of X.T   */
{
  colrange j;
  myTYPE temp=0;
  
  for (j=0; j< nn; j++) {
    temp += X[r-1][j] * T[j][s-1];
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
  static rowset rowexcluded;
  myTYPE Xtemp=0;
  static rowrange mlast=0;
  boolean localdebug=False;

  if (debug) localdebug=True;
  stop = False;
  if (mlast!=mm){
    if (mlast>0) {
      set_free(&rowexcluded);
    }
    set_initialize(&rowexcluded,mm);
    mlast=mm;
  }
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
}

void GausianColumnPivot1(Amatrix A1, Bmatrix T,
  SignAmatrix X, rowrange r, colrange s)
/* Make a column pivot operation in SignAmatrix X on position (r,s).
   Before calling this procedure, the dual basis inverse T must
   be updated by calling GausianColumnPivot2.
 */
{
  long i, j, count=0;
  int Xtemp0, Xtemp;
  float ratio;

  Xtemp0 = X[r - 1][s - 1];
  for (j = 0; j < nn; j++) {
    if (j + 1 != s) {
      Xtemp = X[r - 1][j];
      if (Xtemp0*Xtemp>0){
        for (i = 0; i < mm; i++) {
          if (i + 1 != r)
            if (X[i][j]*X[i][s - 1]==0) 
              X[i][j]-=X[i][s-1];
            else if (X[i][j]*X[i][s - 1]>0){
              X[i][j]=myTYPE2sign(TableauEntry(A1,T,i+1,j+1));
              count++;
            }
        }
      } else if (Xtemp0*Xtemp<0){
        for (i = 0; i < mm; i++) {
          if (i + 1 != r)
            if (X[i][j]*X[i][s - 1]==0) 
              X[i][j]+=X[i][s-1];
            else if (X[i][j]*X[i][s - 1]<0) {
              X[i][j]=myTYPE2sign(TableauEntry(A1,T,i+1,j+1));
              count++;
            }
        }
      }
      X[r - 1][j] = 0;
    }
  }
  if (Xtemp0<0)
    for (i = 0; i < mm; i++) {
      if (i + 1 != r)
        X[i][s - 1] = - X[i][s - 1];
    }
  X[r - 1][s - 1] = 1;
  ratio = 100*count/((mm-nn+1)*nn);
  // cout << "No of real computations = " << count << " : " << ratio  
  //  << " percent\n"; 
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

void SelectCrissCrossPivot1(SignAmatrix X, rowindex OV,
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
        val=X[i-1][rhscol-1];
        if (val < -zero) {
          rowselected=True;
          *r=i;
          break;
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        val=X[objrow-1][bflag[i]-1];
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
          val=X[*r-1][bflag[i]-1];
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
          val=X[i-1][*s-1];
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

void SelectCrissCrossPivot2(Amatrix X, Bmatrix T, rowindex OV,
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

void SelectDualSimplexPivot(boolean Phase1, Amatrix X, Bmatrix T, rowindex OV, 
    colindex NBIndex, long bflag[], rowrange objrow, colrange rhscol,
    rowrange *r, colrange *s,
    boolean *selected, LPStatusType *lps)
{ /* selects a dual simplex pivot (*r, *s) if the current
     basis is dual feasible and not optimal. If not dual feasible,
     the procedure returns *selected=False and *lps=LPSundecided.
     If Phase1=True, the RHS column will be considered as the negative
     of the column of the largest variable (==mm).  For this case, it is assumed
     that the caller used the auxiliary row (with variable mm) to make the current
     dictionary dual feasible before calling this routine so that the nonbasic
     column for mm corresponds to the auxiliary variable.
  */
  boolean colselected=False, rowselected=False, dualfeasible=True,localdebug=False;
  rowrange i,k;
  colrange j;
  myTYPE val=0, minval=0,rat=0, minrat=0;
  static long lastnn=0;
  static myTYPE purezero=0;
  static Arow rcost;

  if (debug) localdebug=True; 
  if (lastnn != nn){
    if (lastnn>0){
       delete[] rcost;
    }
    rcost=new myTYPE[nn];
    lastnn=nn;
  }
  *r=0; *s=0;
  *selected=False;
  *lps=LPSundecided;
  for (j=1; j<=nn; j++){
    if (j!=rhscol){
      if (localdebug) cout << "checking the column " << j<< "  var" <<NBIndex[j] << "\n"; 
      rcost[j-1]=TableauEntry(X,T,objrow,j);
      if (localdebug) cout << "reduced cost =  " << rcost[j-1] << "\n"; 
      if (rcost[j-1] > zero) dualfeasible=False;
    }
  }
  if (dualfeasible){
    while ((*lps==LPSundecided) && (!rowselected) && (!colselected)) {
      for (i=1; i<=mm; i++) {
        if (localdebug) cout << "checking the row var " << i<< "\n"; 
        if (i!=objrow && bflag[i]==-1) {  /* i is a basic variable */
          if (Phase1){
            val=-TableauEntry(X,T,i,bflag[mm]); // for dual Phase I
          } 
          else {val=TableauEntry(X,T,i,rhscol);}
          if (localdebug) cout << "RHS val =  " << val << "\n"; 
          if (val < minval) {
            *r=i;
            minval=val;
            if (localdebug) cout << "update minval with =" << minval << " *r = " << *r <<"\n";
          }
        }
      }
      if (minval>=-zero) *lps=Optimal;
      else {
        rowselected=True;
        for (j=1; j<=nn; j++){
          val=TableauEntry(X,T,*r,j);
          if (j!=rhscol && val > zero) {
            rat=-rcost[j-1]/val;
            if (localdebug) cout << "new ratio =" << rat << "at col = " << j <<"\n";
            if (*s==0 || rat < minrat){
              minrat=rat;
              *s=j;
              if (localdebug) cout << "update minrat =" << minrat << " *s = " << *s <<"\n";
            }
          }
        }
        if (localdebug) cout << "*s is " << *s << "\n";
        if (*s>0) {colselected=True; *selected=True;}
        else *lps=Inconsistent;
      }
    } /* end of while */
  }
}

void FindDualFeasibleBasis(Amatrix X, Bmatrix T, rowset basis, rowset cobasis, rowindex OV, 
    colindex NBIndex, long bflag[], rowrange objrow, colrange rhscol,
    colrange *s, boolean *found, LPStatusType *lps, long *pivot_no)
{ /* Find a dual feasible basis using Phase I of Dual Simplex method.
     If the problem is dual feasible,
     the procedure returns *found=True, *lps=LPSundecided and a dual feasible
     basis.   If the problem is dual infeasible, this returns
     *found=False, *lps=DualInconsistent and the evidence column *s.
     It still can run infinitely and then stops with maxpivots. 
  */
  boolean phase1, dualfeasible=True,localdebug=False,chosen,stop;
  LPStatusType LPSphase1;
  long pivots_p1=0,maxpivots, maxpivfactor=10;
  rowrange i,k,rtemp,entering,leaving,msvar;
  colrange j,l,mmsave,ms=0,stemp;
  myTYPE val=0,purezero=0,maxcost=-1;
  static long lastnn=0;
  static Arow rcost;

  maxpivots=maxpivfactor*nn;
  *found=True; *lps=LPSundecided; *s=0;
  if (debug) localdebug=True;
  mmsave=mm;
  mm=mm+1;  // increase mm by 1 temporally
  if (lastnn != nn){
    if (lastnn>0){
      delete[] X[mm-1];
      delete[] rcost;
    }
    X[mm-1]= new myTYPE[nn];   // create an auxiliary row
    rcost= new myTYPE[nn];   // objective row work space
    lastnn=nn;
  }

  ms=0;  // ms will be the index of column which has the largest reduced cost
  for (j=1; j<=nn; j++){
    if (j!=rhscol){
      if (localdebug) cout << "checking the column " << j << "  var " <<NBIndex[j] << ".\n"; 
      rcost[j-1]=TableauEntry(X,T,objrow,j);
      if (localdebug) cout << "reduced cost =  " << rcost[j-1] << ".\n"; 
      if (rcost[j-1] > maxcost) {maxcost=rcost[j-1]; ms = j; msvar=NBIndex[ms];}
    }
  }
  if (maxcost > zero) dualfeasible=False;

  if (!dualfeasible){
    for (j=1; j<=nn; j++){
      X[mm-1][j-1]=purezero;
      for (l=1; l<=nn; l++){
        if (NBIndex[l]>0) {
          X[mm-1][j-1]-= X[NBIndex[l]-1][j-1];   // To make the auxiliary row (0,-1,-1,...,-1).
        }
      }
    }
    if (localdebug){
      cout << "Auxiliary row =";
      for (j=1; j<=nn; j++){
        cout << " ( " << j << "):" << TableauEntry(X,T,mm,j); 
      }
      cout << "\n";
//    WriteAmatrix(cout,X,mm,nn,Inequality);
      WriteSignTableau(cout, X, T, OV, bflag, objrow, rhscol);
    }

    if (localdebug){
      cout << "FindDualFeasibleBasis: curruent basis is not dual feasible.\n";
      cout << "because of the column " << ms << " assoc. with var " << NBIndex[ms] 
      << "  dual cost = " << maxcost<< ".\n";
    }

    // Pivot on (mm, ms) so that the dual basic solution becomes feasible
    leaving=mm; entering=msvar;
    set_addelem(cobasis, leaving);
    set_delelem(cobasis, entering);
    set_delelem(basis,leaving);
    set_addelem(basis,entering);
    bflag[leaving]=ms;
    bflag[entering]=-1;
    NBIndex[ms]=leaving;
    if (localdebug) {
      cout << "Phase I ini: nonbasis = "; set_write(cobasis);
      cout << "\nPhase I ini: basis = "; set_write(basis);
      cout << "\nPhase I ini: new nonbasic variable " << leaving << " has index " << bflag[leaving] << "\n";
    }
    GausianColumnPivot2(X,T, mm, ms);
    pivots_p1=pivots_p1+1;
    if (ShowSignTableauOn){
      WriteSignTableau(cout, X, T, OV, bflag, objrow, rhscol);
    }

    phase1=True; stop=False;
    do {   /* Dual Simplex Phase I */
      chosen=False; LPSphase1=LPSundecided;
      SelectDualSimplexPivot(phase1, X, T, OV, NBIndex, bflag,
        objrow, rhscol, &rtemp, &stemp, &chosen, &LPSphase1);
      if (!chosen) {
        /* The current dictionary is terminal.  There are two cases:
           TableauEntry(X,T,objrow,ms) is negative or zero.
           The first case implies dual infeasible,
           and the latter implies dual feasible but mm is still in nonbasis.
           We must pivot in the auxiliary variable mm. */

        myTYPE minval=0;
        rtemp=0;
        for (i=1; i<=mm; i++){
          if (bflag[i]<0) { 
             /* i is basic and not the objective variable */
            val=TableauEntry(X,T,i,ms);  // auxiliary column
            if (val < minval) {
              rtemp=i;
              minval=val;
              if (localdebug) cout << "update minval with =" << minval << " rtemp = " << rtemp <<"\n";
            }
          }
        }

        leaving=rtemp; entering=mm;
        set_addelem(cobasis, leaving);
        set_delelem(cobasis, entering);
        set_delelem(basis,leaving);
        set_addelem(basis,entering);
        bflag[leaving]=ms;
        bflag[entering]=-1;
        NBIndex[ms]=leaving;
        if (localdebug) {
          cout << "Phase I fin: nonbasis = "; set_write(cobasis);
          cout << "\nPhase I fin: basis = "; set_write(basis);
          cout << "\nPhase I fin: new nonbasic variable " << leaving << " has index " << bflag[leaving] << "\n";
        }
        GausianColumnPivot2(X,T, rtemp, ms);
        pivots_p1=pivots_p1+1;
        if (ShowSignTableauOn){
          WriteSignTableau(cout, X, T, OV, bflag, objrow, rhscol);
        }

        if (TableauEntry(X,T,objrow,ms)<-zero){
          if (localdebug){
            cout << "Dual infeasible.\n";
            cout << "obj-ms: " << TableauEntry(X,T,objrow,ms) << " zero = " << zero << "\n";
          }
          *found=False; *lps=DualInconsistent;  *s=ms;
        }
        stop=True;
      } else {
        entering=NBIndex[stemp];
        leaving=rtemp;
        if (localdebug) {
          cout<< "Dual Phase I: pivot on (r,s) = (" 
            << rtemp << ", " << stemp << ")\n";
          cout<< "Dual Phase I: (leaving, entering) = (" 
            << leaving << ", " << entering << ")\n";
        }
        set_addelem(cobasis, leaving);
        set_delelem(cobasis, entering);
        set_delelem(basis,leaving);
        set_addelem(basis,entering);
        bflag[leaving]=stemp;
        bflag[entering]=-1;
        NBIndex[stemp]=leaving;
        if (localdebug) {
          cout << "nonbasis = "; set_write(cobasis);
          cout << "\nbasis = "; set_write(basis);
          cout << "\nnew nonbasic variable " << leaving << " has index " << bflag[leaving] << "\n";
        }
        GausianColumnPivot2(X,T, rtemp, stemp);
        pivots_p1=pivots_p1+1;
        if (pivots_p1 > maxpivots) {
          stop=True; *lps=LPSundecided; *found=False;
          if (localdebug) cout << "Dual simplex phase I failed (cycling?)\n";
        }
        if (ShowSignTableauOn){ 
          WriteSignTableau(cout, X, T, OV, bflag, objrow, rhscol);
          WriteCurrentSolution(cout, X, T, objrow, rhscol, NBIndex);
          cout << "Pivot on (" << leaving << ", " << entering << ")\n";
        }
        if (entering==mm) {
          stop=True; 
          if (localdebug) cout << "Dual Phase I: the auxiliary variable enter the basis, go to phase II\n";
        }
      }
    } while(!stop);
  }
  mm=mmsave;
  *pivot_no=pivots_p1;
}

int myTYPE2sign(myTYPE val)
{ int s;

  if (val > zero) s = 1;
  else if (val < -zero) s = -1;
    else s = 0;
  return s;
}

int long2sign(long val)
{ int s;

  if (val > (long)zero ) s = 1;
  else if (val < (long)zero) s = -1;
    else s = 0;
  return s;
}

void FindRedundancy(SignAmatrix SX,
    rowindex OV, long bflag[], colindex NBIndex,
    rowrange objrow, colrange rhscol, rowset redundancy_found)
/* This procedure finds redundant variables in the
   current sign tableau. 
*/ 
{ rowrange i,h,r; 
  colrange j;
  boolean redundant, localdebug=False;
  long poscount;

  if (localdebug) cout << "FindRedundancy\n";
  for (i=1; i<=mm; i++){
    h=OV[i];
    if (localdebug) cout << "checking the row " << h <<"\n";
    if (bflag[h]<0) { 
      /* h is basic and not the objective variable */
      if (SX[h-1][rhscol-1] >= 0 && !set_member(h, redundancy_found)) {
        j=1; redundant=True; 
        while (redundant && j<=nn){
          if (SX[h-1][j-1] < 0) redundant=False; 
          j++;
        }
        if (redundant) {
          set_addelem(redundancy_found, h);
          if (localdebug) cout << "redundancy of basic " << h << " detected.\n";
        }  
      }
      if ( SX[h-1][rhscol-1] <= 0) {
        poscount=0; 
        for (j=1; j<=nn; j++){
          if (SX[h-1][j-1] > 0) {
            poscount++; 
            r = NBIndex[j];
          } 
        }
        if (poscount==1 && !set_member(r, redundancy_found)) {
          set_addelem(redundancy_found, r);
          if (localdebug) cout << "redundancy of nonbasic " << r << " detected.\n";
        }
      }  
    }
  } /* end of for i */
}

void OptimizeVarOrder(SignAmatrix SX,
    rowindex OV, rowindex L,  long bflag[], 
    colindex NBIndex, rowset Fixed, 
    rowrange objrow, colrange rhscol,
    rowrange LP_size, rowrange B_reduc, colrange N_reduc, int *opt_type)
/* This procedure finds largest subproblems solved, and 
   reorder variables if one can increase the highest nonzero position of
   the current L vector.
   opt_type returns the type (=0,1,2,3) of reordering.  If it is 0,
   no reordering is applied.   
*/ 
{ myTYPE val=0;
  long h,bn1,bn2,bn3;
  long key_digit=0;
  rowrange i,k, b1=0, b2=mm-nn-B_reduc, r2=0, b3=0, b3temp, s3=0; 
  colrange j, n1=0, n2=0, n3=nn-1-N_reduc, n2temp;
  int sval=0;
  Amatrix P;  // priority matrix
 
  *opt_type=0;
  for (i=1; i<=mm; i++) {
    P[i-1] = new myTYPE[2];
    for (j=1; j<=2; j++) P[i-1][j-1] = 0;
  }

  for (i=1; i<=LP_size; i++){
    b3temp=0; n2temp=0;
    h=OV[i];
    if (L[h]==1) key_digit=i;
    switch (long2sign(bflag[h])) { 
      case -1:  /* h is basic and not the objective variable */
       if (SX[h-1][rhscol-1] >=0)  b1++;
        else{ /* SX[h-1][rhscol-1] < 0  */
          for (j=1; j<=nn ; j++){
            if (j!=rhscol && SX[h-1][j-1] <= 0 && 
              !set_member(NBIndex[j], Fixed)) n2temp++;
          }
        }
        break;

      case  1:  /* h is nonbasic */
        j=bflag[h];
        if (SX[objrow-1][j-1] <= 0)  n1++;
        else{ /* SX[objrow-1][j-1] > 0  */
          for (k=1; k<=mm; k++){
            if (bflag[k]==-1 && SX[k-1][j-1] >= 0  &&
              !set_member(k, Fixed)) b3temp++;
          } 
        }
        break;
    } /* end of switch */
    if (n2temp > n2) {n2 =n2temp; r2=h;}
    if (b3temp > b3) {b3 =b3temp; s3=h;}
  } /* end of for i */

/*
  cout << " (b1, n1) = " << b1 << ", " << n1;
  cout << ": (b2, n2) = " << b2 << ", " << n2;
  cout << ": (b3, n3) = " << b3 << ", " << n3 << "\n";
  cout << " r2 = " << r2 << ":   s3 = " << s3 << "\n";
  cout << " key digit = " << key_digit << "\n";
  cout << " fixed variables = "; set_write(Fixed); cout << "\n";
*/

  // bn1=b1*n1; bn2=b2*n2; bn3=b3*n3;
  bn1=b1+n1; bn2=b2+n2; bn3=b3+n3;
  if ((b1+n1 >= key_digit) && bn1>=bn2 && bn1>=bn3){
    if (ShowSignTableauOn) cout << "variable reordering of type I is applied.\n";
    for (i=1; i<=mm; i++) {
      if (bflag[i]==-1) {  /* i is a basic variable */
        P[i-1][0]=-SX[i-1][rhscol-1];
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        P[i-1][0]=SX[objrow-1][bflag[i]-1];
      }
    }
    for (i=LP_size+1; i<=mm; i++) { /* fixed variables should be ordered last */
      P[OV[i]-1][0]=i+2;  /* large enough so that they will be ordered last */
      P[OV[i]-1][1]=i+2;  /* large enough so that they will be ordered last */
    } 
    QuickSort(OV, 1, mm-1, P, 2);
    for (i=0; i<=LP_size; i++) L[OV[i]]=0;
    *opt_type=1;
  }
  else if ((b2+n2 >= key_digit) && bn2>=bn1 && bn2>=bn3){
    if (ShowSignTableauOn) cout << "variable reordering of type II is applied.\n";
    for (i=1; i<=mm; i++) {
      if (bflag[i]==-1) {  /* i is a basic variable */
        if (SX[i-1][rhscol-1]>=0 || i==r2){
          P[i-1][0]=-SX[i-1][rhscol-1];
          P[i-1][1]=0;
        } else{
          P[i-1][0]=-2*SX[r2-1][rhscol-1]; /* large positive number */
          P[i-1][1]=SX[i-1][rhscol-1]; /* negative number */
        }
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        P[i-1][0]=-2*SX[r2-1][rhscol-1]; /* large positive number */
        P[i-1][1]=SX[r2-1][bflag[i]-1]; 
          /* nonbasic var i with nonpositive r2-row comes earliear than pos r2-row */
      }
    }
    for (i=LP_size+1; i<=mm; i++) { /* fixed variables should be ordered last */
      P[OV[i]-1][0]=i+2;  /* large enough so that they will be ordered last */
      P[OV[i]-1][1]=i+2;  /* large enough so that they will be ordered last */
    } 
    QuickSort(OV, 1, mm-1, P, 2);
    for (i=0; i<=LP_size; i++) L[OV[i]]=0;
    *opt_type=2;
  } 
  else if ((b3+n3 >= key_digit) && bn3>=bn1 && bn3>=bn2){ 
    if (ShowSignTableauOn) cout << "variable reordering of type III is applied.\n";
    for (i=1; i<=mm; i++) {
      if (bflag[i]==-1) {  /* i is a basic variable */
        P[i-1][0]=2*SX[objrow-1][bflag[s3]-1]; /* large positive number */
        P[i-1][1]=-SX[i-1][bflag[s3]-1];  /* negative of s3-column */
      }
      else if (bflag[i] >0) { /* i is nonbasic variable */
        if (SX[objrow-1][bflag[i]-1] <= 0 || i==s3){
          P[i-1][0]=SX[objrow-1][bflag[i]-1];  
          P[i-1][1]=0;
        }else{
          P[i-1][0]=2*SX[objrow-1][bflag[s3]-1]; /* large positive number */
          P[i-1][1]=-SX[objrow-1][bflag[i]-1]; 
           /* negative number so that it comes earliear than basic var with negative
              s3-column */
        }
      }
    }
    for (i=LP_size+1; i<=mm; i++) { /* fixed variables should be ordered last */
      P[OV[i]-1][0]=i+2;  /* large enough so that they will be ordered last */
      P[OV[i]-1][1]=i+2;  /* large enough so that they will be ordered last */
    } 
    QuickSort(OV, 1, mm-1, P, 2);
    for (i=0; i<=LP_size; i++) L[OV[i]]=0;
    *opt_type=3;
  }
 
  for (i=1; i<=mm; i++) delete[] P[i-1];
}


void CrissCrossMinimize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, boolean UsePrevBasis, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
{
   colrange j;
   
   for (j=1; j<=nn; j++)
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   CrissCrossMaximize(f, f_log, A1,BasisInverse, OBJrow, RHScol, UsePrevBasis,
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter);
   *optvalue=-*optvalue;
   for (j=1; j<=nn; j++){
     dsol[j-1]=-dsol[j-1];
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   }
}

void DualSimplexMinimize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, boolean UsePrevBasis, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
{
   colrange j;
   
   for (j=1; j<=nn; j++)
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   DualSimplexMaximize(f, f_log, A1,BasisInverse, OBJrow, RHScol, UsePrevBasis, 
     LPS, optvalue, sol, dsol, NBIndex, re,  se, iter);
   *optvalue=-*optvalue;
   for (j=1; j<=nn; j++){
     dsol[j-1]=-dsol[j-1];
     AA[OBJrow-1][j-1]=-AA[OBJrow-1][j-1];
   }
}

long OrderOf(rowrange p, rowindex OV)
{
  rowrange i=0;
  boolean stop=False;

  while (!stop) { 
    i++;
    if (i>mm || OV[i]==p) stop=True;
  }
  if (i>mm) return 0;else return i;
}
 

void UpdateLvector(rowindex L, rowrange p, rowindex OV)
{
  rowrange op, i=1;

  op = OrderOf(p,OV);  /* op is the OV order of a row index p */
  
  for (i=1; i<op; i++) L[OV[i]]=0;
  L[p] = 1;
  if (L[0] < op) L[0]=op; /* L[0] holds the highest 1 position */
}

void WriteLvector(ostream &f, rowindex OV, rowindex L)
{
  rowrange i;

  if (mm <= 30){
    f << "Variable :";
    for (i=1; i<= mm; i++){
      if (OV[i]!=OBJrow){
        f.width(3);
        f << OV[i];
      }
    }
    f << "\nL vector :";
    for (i=1; i<= mm; i++){
      if (OV[i]!=OBJrow){
        f.width(3);
        f << L[OV[i]];
      }
    }
    f << "\n";
  }
  f << "L vector :";
  for (i=mm; i>= 1; i--){
    if (OV[i]!=OBJrow){
      if (L[OV[i]]==1){
        f.width(4); 
        f << i;
      }
    }
  }
  f << "\n";
}

void WriteCurrentSolution(ostream &f, Amatrix A1, Bmatrix Binv, rowrange OBJrow,
   colrange RHScol, colindex NBIndex)
{
  static Arow sol,dsol;
  myTYPE optvalue;
  rowrange i;
  colrange j;
  static long lastnn=0;

  if (lastnn != nn){
    if (lastnn>0){
       delete[] sol; delete[] dsol;
    }
    sol=new myTYPE[nn];
    dsol=new myTYPE[nn];
    lastnn=nn;
  }
  for (j=1;j<=nn; j++) {
    sol[j-1]=Binv[j-1][RHScol-1];
    dsol[j-1]=-TableauEntry(A1,Binv,OBJrow,j);
    optvalue=TableauEntry(A1,Binv,OBJrow,RHScol);
  }
  f << "Primal sol: (";
  for (j=1; j<=nn; j++){
    if (j>1) f << ", " << sol[j-1];
    else f << " " << sol[j-1];
  }
  f << ") ; obj_value=" << optvalue << "\n"; 
}


void CrissCrossMaximize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, boolean UsePrevBasis, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  boolean stop, chosen;
  long rank, ideg;
  long pivots_ds=0, pivots_cc=0;
  rowrange i,r,entering,leaving,keyindex;
  rowrange lp_size, sublp_size, b_reduc, subb_reduc, e_order, l_order;
  colrange n_reduc, subn_reduc, j, s;
  static colset ColSelected;
  static rowset RowSelected,Basis,Cobasis,FixedVariables,subFixedVariables,Redundancy;
  static rowindex BasisFlag;
  static long mlast=0,nlast=0;
  myTYPE adet=0; /* abs value of the determinant of a basis */
  boolean localdebug=False;
  static rowindex L, OrderVec;
  static SignAmatrix SA;
  int opt_type=0; /* reordering type variable for optimize_CC option */

  if (debug) localdebug=True;
  if (mlast!=mm || nlast!=nn){
     if (mlast>0) { /* called previously with different mm */
       if (localdebug) cout << "CCmaximize: deleting the old memory space with mlast = " << mlast << "\n";
       delete[] BasisFlag;
       delete[] L;
       delete[] OrderVec;
       set_free(&ColSelected);
       set_free(&RowSelected);
       set_free(&Basis);
       set_free(&Cobasis);
       set_free(&FixedVariables);
       set_free(&subFixedVariables);
       set_free(&Redundancy);
       for (i=0; i<mlast; i++) delete[] SA[i];
     }
     if (localdebug) cout << "CCmaximize: allocating a new memory space with mm = " << mm<< "\n";
     BasisFlag=new long[mm+1];  
     L = new long[mm+1];
     OrderVec = new long[mm+1];
     set_initialize(&Cobasis,mm);
     set_initialize(&Basis,mm);
     set_initialize(&RowSelected, mm);
     set_initialize(&ColSelected, nn);
     set_initialize(&FixedVariables, mm);
     set_initialize(&subFixedVariables, mm);
     set_initialize(&Redundancy, mm);
     for (i=0; i<mm; i++) SA[i]=new int[nn];
     /* initialize only for the first time or when a larger space is needed */
     mlast=mm;nlast=nn;
  }
  *re=0; *se=0; *iter=0;
  rank = 0;
  lp_size=mm-1; b_reduc=0; n_reduc=0; 
    /* these sizes will change when some basic or nonbasic variables are fixed */ 
  stop = False;
  PreOrderedRun=False;
  adet=1.0;
  set_emptyset(Cobasis);
  set_emptyset(Basis);
  set_emptyset(RowSelected);
  set_emptyset(ColSelected);
  set_emptyset(FixedVariables);
  set_emptyset(subFixedVariables);
  set_emptyset(Redundancy);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  for (j=0; j<=nn; j++) NBIndex[j]=0;
  for (i=0; i<=mm; i++) {
    L[i]=0; OrderVec[i]=i; 
  }
  for (i=1; i<=mm; i++) {
    BasisFlag[i]=0;
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
  if (ManualPivotOn){ 
    ManualPivot(f,f_log,A1,BasisInverse,OBJrow,RHScol,LPS,optvalue,sol,dsol,NBIndex, re,se,iter);
    for (i=1; i<=mm; i++){ set_addelem(Basis,i);}
    for (j=1; j<=nn; j++){
      i=NBIndex[j];
      set_addelem(Cobasis,i);
      set_delelem(Basis,i);
      BasisFlag[i]=j;
    }
  } 
  else do {   /* Find a LP basis */
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
  ComputeRowOrderVector(OrderVec, HyperplaneOrder);
  if (localdebug) {
    WriteLvector(f_log,OrderVec,L);
  }
  SetSignAmatrix(SA,A1,BasisInverse);
  if (ShowSignTableauOn){
    if (nn <= 11) WriteSignAmatrix(cout,SA,OrderVec, BasisFlag, OBJrow, RHScol);
    WriteLvector(cout,OrderVec,L); 
  }
  do {   /* Criss-Cross Method */
    if (LPsolver==CombMaxImprove){
      opt_type=0; 
      sublp_size=lp_size;
      subb_reduc=b_reduc;
      subn_reduc=n_reduc;
      set_copy(subFixedVariables,FixedVariables);
      while (opt_type==0 && sublp_size>0){
        // cout << "Order optimization applied to subproblem size = " << sublp_size << "\n"; 
        OptimizeVarOrder(SA, OrderVec, L, BasisFlag, NBIndex, 
          subFixedVariables, OBJrow, RHScol, sublp_size, subb_reduc, subn_reduc, &opt_type);
        if (opt_type==0){
          i=sublp_size; chosen=False;
          while (i>0 && !chosen){
            if (BasisFlag[OrderVec[i]]<0) subb_reduc++;
            else subn_reduc++;
            set_addelem(subFixedVariables,OrderVec[i]);
            if (L[OrderVec[i]]==1) chosen=True;
            i--;
          }
          sublp_size = i;
        }else{
          if (ShowSignTableauOn){
            WriteSignAmatrix(cout,SA,OrderVec, BasisFlag, OBJrow, RHScol);
            WriteDictionary(cout, A1, BasisInverse, OrderVec, BasisFlag, NBIndex, OBJrow, RHScol);
            WriteLvector(cout,OrderVec,L);
          }
        } 
      }
    }
    chosen=False; *LPS=LPSundecided;
    if (!chosen && *LPS==LPSundecided) {
      if (LPsolver==CombMaxImprove || SignPivotOn) {
        SelectCrissCrossPivot1(SA, OrderVec, BasisFlag, 
          OBJrow, RHScol, &r, &s, &chosen, LPS);
        FindRedundancy(SA, OrderVec, BasisFlag, NBIndex,
          OBJrow, RHScol, Redundancy);
        if (chosen) pivots_cc=pivots_cc+1;
      } else {
        SelectCrissCrossPivot2(A1, BasisInverse, OrderVec, BasisFlag, 
          OBJrow, RHScol, &r, &s, &chosen, LPS);
        if (chosen) pivots_cc=pivots_cc+1;
      }
    }
    if (chosen) {
      entering=NBIndex[s];
      leaving=r;
      if (localdebug) {
        cout<< "Procedure Criss-Cross: pivot on (r,s) = (" 
          << r << ", " << s << ")\n";
        cout<< "Procedure Criss-Cross: (leaving, entering) = (" 
          << leaving << ", " << entering << ")\n";
      }
      e_order=OrderOf(entering,OrderVec);
      l_order=OrderOf(leaving,OrderVec);
      if (e_order > l_order) {
        keyindex = entering;
        UpdateLvector(L,keyindex,OrderVec);
        if (e_order==lp_size) {
          L[0]=0;
          lp_size--;  /* lp size is reduced */
          b_reduc++;  /* the last active variable is fixed to basic */
          set_addelem(FixedVariables, entering);
        }
      }else {
        keyindex = leaving;
        UpdateLvector(L,keyindex,OrderVec);
        if (l_order==lp_size){
          L[0]=0;
          lp_size--;  /* lp size is reduced */
          n_reduc++;  /* the last active variable is fixed to basic */
          set_addelem(FixedVariables, leaving);
        }
      }
      if (localdebug){
        WriteLvector(f_log,OrderVec,L);
        // f_log << "cobasis: "; set_fwrite(f_log,Cobasis); f_log << "\n";
      }
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
      if (SignPivotOn) GausianColumnPivot1(A1,BasisInverse, SA, r, s);
      (*iter)++;
      if (localdebug) WriteCurrentSolution(f_log, A1, BasisInverse, OBJrow, RHScol, NBIndex);
      if (ShowSignTableauOn){ 
        WriteCurrentSolution(cout, A1, BasisInverse, OBJrow, RHScol, NBIndex);
        cout << "Pivot on (" << leaving << ", " << entering << ")\n";
        if (SignPivotOn){
          WriteSignAmatrix(cout,SA,OrderVec, BasisFlag, OBJrow, RHScol);
        }else{
          WriteSignTableau(cout,A1,BasisInverse,OrderVec, BasisFlag, OBJrow, RHScol);
          WriteDictionary(cout, A1, BasisInverse, OrderVec, BasisFlag, NBIndex, OBJrow, RHScol);
        }
        WriteLvector(cout, OrderVec, L);
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
  if (localdebug){
    WriteLvector(f_log,OrderVec,L);
    // f_log << "cobasis: "; set_fwrite(f_log,Cobasis); f_log << "\n";
  }
  if (DynamicWriteOn && LogWriteOn){
     cout << "LP solved with " << *iter << " pivots. (ds piv#=" << 
       pivots_ds << ", cc piv#=" << pivots_cc << ")\n"; 
  }
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
}


void DualSimplexMaximize(ostream &f, ostream &f_log,
   Amatrix A1,Bmatrix BasisInverse, 
   rowrange OBJrow, colrange RHScol, boolean UsePrevBasis, LPStatusType *LPS,
   myTYPE *optvalue, Arow sol, Arow dsol, colindex NBIndex,
   rowrange *re, colrange *se, long *iter)
/* 
When LP is inconsistent then *re returns the evidence row.
When LP is dual-inconsistent then *se returns the evidence column.
*/
{
  boolean stop, chosen, phase1, found;
  long rank;
  long pivots_ds=0, pivots_p1=0, pivots_pc=0, maxpivots, maxpivfactor=70;
  rowrange i,r,entering,leaving;
  colrange j,s;
  static colset ColSelected;
  static rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0,nlast=0;
  myTYPE adet=1; /* abs value of the determinant of a basis */
  boolean localdebug=False;

  *LPS=LPSundecided;
  maxpivots=maxpivfactor*nn;  /* maximum pivots to be performed before cc pivot is applied. */
  if (debug) localdebug=True;
  if (mlast!=mm || nlast!=nn){
     if (mlast>0) { /* called previously with different mm */
       if (localdebug) cout << "DualSimplex: deleting the old memory space with mlast = " << mlast << "\n";
       delete[] BasisFlag;
       set_free(&ColSelected);
       set_free(&RowSelected);
       set_free(&Basis);
       set_free(&Cobasis);
     }
     if (localdebug) cout << "DualSimplex: allocating a new memory space with mm = " << mm<< "\n";
     BasisFlag=new long[mm+2];  // one more element for an auxiliary variable
     set_initialize(&Cobasis,mm+1);
     set_initialize(&Basis,mm+1);
     set_initialize(&RowSelected, mm+1);
     set_initialize(&ColSelected, nn);
     /* initialize only for the first time or when a larger space is needed */
     mlast=mm;nlast=nn;
  }
  /* Initializing control variables. */
  BasisFlag[OBJrow]= 0; /*  BasisFlag of the objective variable is 0, 
    different from other basic variables which have -1 */
  set_emptyset(RowSelected);
  set_emptyset(ColSelected);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  *re=0; *se=0; *iter=0;

  PreOrderedRun=True;
  if (UsePrevBasis){
    if (s=BasisFlag[OBJrow]>0){ /* OBJrow is in cobasis, and we must pivot it in. */
      for (j=0;j<=mm+1;j++){if (j!=s) set_addelem(ColSelected,j);}
      rank=nn-2;
      if (localdebug) cout << "UsePrevBasis but the current basis does not contain OBJrow.\n";
      set_copy(RowSelected,Cobasis);
    } else {
      if (localdebug){
         cout << "UsePrevBasis. The current basis contains OBJrow and thus an LP basis.\n";
        cout << " Current nonbasis = ";   
        set_write(Cobasis); cout << "\n";
      }
      goto _L88;
    }
  } else {
    if (localdebug) cout << "Do not UsePrevBasis.\n";
    rank = 0;
    set_emptyset(Cobasis);
    set_emptyset(Basis);
    for (j=0; j<=nn; j++) NBIndex[j]=0;
    for (i=1; i<=mm+1; i++) {
      BasisFlag[i]=0;
      set_addelem(Basis,i);
      BasisFlag[i]=-1;    /* basic variable has index -1 */
      if (EqualityIndex[i]==-1){
        if (DynamicWriteOn){ 
          cout << "*Warning: strict_inquality option for row " << i << " is ignored for maximization/minimization\n";
        }
      }
    }
    SetToIdentity(BasisInverse);
  }

  stop = False;
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
        cout << "DSimplex: find initial basis: nonbasis = ";   
        set_write(Cobasis); cout << "\n";
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

_L88:
  stop=False;
  if (ShowSignTableauOn){
    WriteSignTableau(cout, A1, BasisInverse, OrderVector, BasisFlag, OBJrow, RHScol);
    WriteDictionary(cout, A1, BasisInverse, OrderVector, BasisFlag, NBIndex, OBJrow, RHScol);
  }

  FindDualFeasibleBasis(A1, BasisInverse, Basis, Cobasis, OrderVector, NBIndex,BasisFlag,
      OBJrow, RHScol, &s, &found, LPS, &pivots_p1);
  *iter+=pivots_p1;
  if (*LPS!=LPSundecided && !found){
     *se=s;   
     // No dual feasible basis is found, and thus DualInconsistent.  Output the evidence column.
  }
  else {
    do {   /* Dual Simplex Method */
      chosen=False; *LPS=LPSundecided; phase1=False;
      if (pivots_ds<maxpivots) {
        SelectDualSimplexPivot(phase1, A1, BasisInverse, OrderVector, NBIndex, BasisFlag,
        OBJrow, RHScol, &r, &s, &chosen, LPS);
      }
      if (chosen) pivots_ds=pivots_ds+1;
      if (!chosen && *LPS==LPSundecided) {  
        /* In principle this should not be executed because we already have dual feasibility
           attained and dual simplex pivot should have been chosen.  This might occur
           under floating point computation, or the case of cycling.
        */
        SelectCrissCrossPivot2(A1, BasisInverse, OrderVector, BasisFlag, 
          OBJrow, RHScol, &r, &s, &chosen, LPS);
        if (chosen) pivots_pc=pivots_pc+1;
      }
      if (chosen) {
        entering=NBIndex[s];
        leaving=r;
        if (localdebug) {
          cout<< "Phase 1: pivot on (r,s) = (" 
            << r << ", " << s << ")\n";
          cout<< "Phase 1: (leaving, entering) = (" 
            << leaving << ", " << entering << ")\n";
        }
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
          WriteCurrentSolution(cout, A1, BasisInverse, OBJrow, RHScol, NBIndex);
          WriteDictionary(cout, A1, BasisInverse, OrderVector, BasisFlag, NBIndex, OBJrow, RHScol);
          cout << "Pivot on (" << leaving << ", " << entering << ")\n";
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
  }
  if (DynamicWriteOn && LogWriteOn){
     cout << "LP solved with " << *iter << " pivots. (ds piv#=" << 
       pivots_ds << ", p1 piv#=" << pivots_p1;
     if (pivots_pc > 0) cout << ", cc piv#=" << pivots_pc;
     cout << ")\n";
  }
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
    if (localdebug) printf("DualSimplexSolve: LP is inconsistent.\n");
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
    if (localdebug) cout << "DualSimplexSolve: LP is dual inconsistent.\n";
    break;

  default:break;
  }
}

void ManualPivot(ostream &f, ostream &f_log,
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
  rowrange i,r,entering,leaving;
  rowrange lp_size;
  colrange j, s;
  static colset ColSelected;
  static rowset RowSelected,Basis,Cobasis;
  static rowindex BasisFlag;
  static long mlast=0,nlast=0;
  myTYPE adet=0; /* abs value of the determinant of a basis */
  static rowindex L, OrderVec;
  static SignAmatrix SA;
  int opt_type=0; /* reordering type variable for optimize_CC option */

  if (mlast!=mm || nlast!=nn){
     if (mlast>0) { /* called previously with different mm */
       delete[] BasisFlag;
       delete[] OrderVec;
       set_free(&ColSelected);
       set_free(&RowSelected);
       set_free(&Basis);
       set_free(&Cobasis);
       for (i=0; i<mlast; i++) delete[] SA[i];
     }
     BasisFlag=new long[mm+1];  
     OrderVec = new long[mm+1];
     set_initialize(&Cobasis,mm);
     set_initialize(&Basis,mm);
     set_initialize(&RowSelected, mm);
     set_initialize(&ColSelected, nn);
     for (i=0; i<mm; i++) SA[i]=new int[nn];
     /* initialize only for the first time or when a larger space is needed */
     mlast=mm;
  }
  *re=0; *se=0; *iter=0;
  rank = 0;
  lp_size=mm-1; 
  stop = False;
  PreOrderedRun=False;
  adet=1.0;
  set_emptyset(Cobasis);
  set_emptyset(Basis);
  set_emptyset(RowSelected);
  set_emptyset(ColSelected);
  set_addelem(RowSelected, OBJrow);
  set_addelem(ColSelected, RHScol);
  for (j=0; j<=nn; j++) NBIndex[j]=0;
  for (i=0; i<=mm; i++) {
    OrderVec[i]=i; 
  }
  for (i=1; i<=mm; i++) {
    BasisFlag[i]=0;
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
  do {   /* Find a LP basis */
      SelectPivot2(A1, BasisInverse, MinIndex
      , mm, RowSelected, ColSelected, &r, &s, &chosen);
      if (chosen) {
        set_addelem(RowSelected, r);
        set_addelem(ColSelected, s);
        set_addelem(Cobasis, r);
        set_delelem(Basis,r);
        BasisFlag[r]=s;   /* the nonbasic variable r corresponds to column s */
        NBIndex[s]=r;     /* the nonbasic variable on s column is r */
        rank++;
        GausianColumnPivot2(A1,BasisInverse, r, s);
      } else {
        stop=True;
      }
      if (rank==nn-1) stop = True;
  } while (!stop);
  
  stop=False;
  ShowSignTableauOn=True;
  SignPivotOn=True;
  ComputeRowOrderVector(OrderVec, HyperplaneOrder);
  SetSignAmatrix(SA,A1,BasisInverse);
  if (ShowSignTableauOn){
    WriteSignAmatrix(cout,SA,OrderVec, BasisFlag, OBJrow, RHScol);
    WriteDictionary(cout, A1, BasisInverse, OrderVec, BasisFlag, NBIndex, OBJrow, RHScol);
  }
  do {   /* Manual Povot */
    chosen=False;
    while (!stop && !chosen){ 
      cout << "Select pivot (r,s): r  s = ";
      cin >> r;  cin >> entering;
      if (r==0) stop=True;
      else if (!set_member(r,Basis)){cout << r << " is not in the basis\n";}
      else if (!set_member(entering,Cobasis)) {cout << entering << " is not in the cobasis\n";}
      else{
        chosen=True; s=BasisFlag[entering];
      }  
    } 
    if (chosen) {
      // entering=NBIndex[s];
      leaving=r;
      set_addelem(Cobasis, leaving);
      set_delelem(Cobasis, entering);
      set_delelem(Basis,leaving);
      set_addelem(Basis,entering);
      BasisFlag[leaving]=s;
      BasisFlag[entering]=-1;
      NBIndex[s]=leaving;
      GausianColumnPivot2(A1,BasisInverse, r, s);
      GausianColumnPivot1(A1,BasisInverse, SA, r, s);
      (*iter)++;
      WriteCurrentSolution(cout, A1, BasisInverse, OBJrow, RHScol, NBIndex);
      cout << "Pivot on (" << leaving << ", " << entering << ")\n";
      // f_log << "cobasis: "; set_fwrite(f_log,Cobasis); f_log << "\n";
      WriteSignAmatrix(cout,SA,OrderVec, BasisFlag, OBJrow, RHScol);
      WriteDictionary(cout, A1, BasisInverse, OrderVec, BasisFlag, NBIndex, OBJrow, RHScol);
    }
  } while(!stop);
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

void ShiftPointsAroundOrigin(ostream &f, ostream &f_log, Arow center)
/* Shifting all points in AA so that the origin
   lies in the (relative) interior of its convex hull
*/
{
  long i,j,maxpoint,minpoint;
  myTYPE maxvalue=0, minvalue=0, two=2;
  boolean localdebug=False;

  /* compute a center */
  f << "* Computing 2d coordinate-wise extremal points:\n* min(1) max(1) ... min(d) max(d) =";
  if (DynamicWriteOn) cout << "* Computing 2d coordinate-wise extremal points:\n* min(1) max(1) ... min(d) max(d) =";
  for (j=1; j <= nn; j++) {
    center[j-1]=0;
    i=1;
    maxpoint=i;maxvalue=AA[maxpoint-1][j-1];
    minpoint=i;minvalue=AA[minpoint-1][j-1];
    for (i=2; i <= mm; i++) {
      if (maxvalue < AA[i-1][j-1]) {
        maxpoint=i;maxvalue=AA[maxpoint-1][j-1];
      }
      if (minvalue > AA[i-1][j-1]) {
        minpoint=i;minvalue=AA[minpoint-1][j-1];
      }
    }
    center[j-1]=(maxvalue+minvalue)/two;
    f << " " << minpoint << " " << maxpoint;
    if (DynamicWriteOn) cout << " " << minpoint << " " << maxpoint;
  }
  f << "\n* Shifting all points w.r.t. the center (of 2d points) =\n* Point0 =";
  for (j=2; j <= nn; j++) {
    WriteNumber(f, center[j-1]);
  }
  f << "\n";
  if (DynamicWriteOn) {
    cout << "\n* Shifting all points w.r.t. the center (of 2d points) =\n* Point0 =";
    for (j=2; j <= nn; j++) {
      WriteNumber(cout, center[j-1]);
    }
    cout << "\n";
  }

  /* do shifting */
  for (i=1; i <= mm; i++) {
    for (j=2; j <= nn; j++) {
      AA[i-1][j-1]-=center[j-1];
      if (localdebug) WriteNumber(cout, AA[i-1][j-1]);
    }
    if (localdebug) cout << "\n";
  }
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
  for (i=0; i < mnew; i++) {
    for (j=0; j < nnew; j++) {
      Acopy[i][j]=AA[i][j];
      if (debug) WriteNumber(cout, AA[i][j]);
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

void EnlargeAAforZeroRHSLP(void)
/* Add an extra column with all zeros to the matrix AA, 
   and change nn accordingly 
*/
{
  long i,j,nnew=0;
  Amatrix Acopy;

  nnew=nn+1;
  // Objective row mm  AA[mm-1] is correctly stored already. No change.
  for (i=1; i<=mm-1; i++){
    Acopy[i-1]= new myTYPE[nnew];
  }
  for (i=1; i <= mm-1; i++) {
    for (j=2; j <= nnew; j++) {
      Acopy[i-1][j-1]=AA[i-1][j-2];
      if (debug) WriteNumber(cout, AA[i-1][j-2]);
    }
    if (debug) cout << "\n";
  }
  for (i=1;i<=mm-1;i++) {
    Acopy[i-1][0]=0;  /* new column with all minus one's */
  }
  for (i=1;i<=mm-1;i++){
    delete[] AA[i-1];
  }
  for (i=1;i<=mm-1;i++){
    AA[i-1]=Acopy[i-1];
  }
  nn=nnew;
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
