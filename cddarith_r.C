/* cddarith.C:  Arithmetic Procedures for cdd.C
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


myTYPE FABS(myTYPE x)
{
  myTYPE pure_zero=0;

  if (x < pure_zero) return -x; else return x;
}

void SetInequalitySets(rowindex id)
{
  long i;
  static long mprev=0;
  boolean localdebug=False;
  
  if (mprev!=mm){
    if (localdebug) printf("SetInequalitySets: initializing inequality sets.\n");
    if (GroundSet!=NULL) set_free(&GroundSet);
    if (EqualitySet!=NULL) set_free(&EqualitySet);
    if (NonequalitySet!=NULL) set_free(&NonequalitySet);
    set_initialize(&EqualitySet, mm);
    set_initialize(&NonequalitySet, mm);
    set_initialize(&GroundSet, mm);
  }
  if (localdebug && mprev==mm) printf("SetInequalitySets: Resetting inequality sets.\n");
  set_emptyset(GroundSet);
  set_emptyset(EqualitySet);
  set_emptyset(NonequalitySet);  
  for (i = 1; i <= mm; i++){
    set_addelem(GroundSet, i);
    if (id[i]==1) set_addelem(EqualitySet,i);
    if (id[i]==-1) set_addelem(NonequalitySet,i);
  }
  mprev=mm;
}


void CheckAdjacency1(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  long rank;

  *adjacent = True;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (debug)
    printf("Check adjacency\n");
  if (set_card(Face)< nn - 2) {
    *adjacent = False;
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = True;
  	return;
  }
  ComputeRank(AA,Face,&rank);
  if (rank < nn - 2){
    *adjacent = False;
  }
}

void CheckAdjacency2(RayRecord **RP1, RayRecord **RP2,
			    boolean *adjacent)
{
  RayRecord *TempRay;
  boolean localdebug=False;

  if (debug) localdebug=True;
  *adjacent = True;
  set_int(Face1, (*RP1)->ZeroSet, (*RP2)->ZeroSet);
  set_int(Face, Face1, AddedHyperplanes);
  if (localdebug){
    printf("Check adjacency of\n");
    WriteRayRecord(cout, *RP1);
    WriteRayRecord(cout, *RP2);    
  }
  if (set_card(Face)< nn - 2) {
    *adjacent = False;
    if (localdebug) {
      printf("non adjacent: set_card(face) %ld < %ld = nn.\n",
        set_card(Face),nn);
    }
    return;
  }
  else if (NondegAssumed) {
  	*adjacent = True;
  	return;
  }
  TempRay = FirstRay;
  while (TempRay != NULL && *adjacent) {
    if (TempRay != *RP1 && TempRay != *RP2) {
    	set_int(Face1, TempRay->ZeroSet, AddedHyperplanes);
      	if (set_subset(Face, Face1)) *adjacent = False;
    }
    TempRay = TempRay->Next;
  }
}

void Eliminate(RayRecord **Ptr)
{
  /*eliminate the record pointed by Ptr^.Next*/
  RayRecord *TempPtr;

  if (debug) {
    printf("            Delete:");
    WriteRayRecord(cout, (*Ptr)->Next);
  }
  TempPtr = (*Ptr)->Next;
  (*Ptr)->Next = (*Ptr)->Next->Next;
  if (TempPtr == FirstRay)   /*Update the first pointer*/
    FirstRay = (*Ptr)->Next;
  if (TempPtr == LastRay)   /*Update the last pointer*/
    LastRay = *Ptr;
  delete[] TempPtr->Ray;          /* free the ray vector memory */
  delete TempPtr->ARay;         /* free the Axray value memory */
  set_free(&(TempPtr->ZeroSet));  /* free the ZeroSet memory */
  delete TempPtr;   /* free the RayRecord structure memory */
  RayCount--; 
}


void SelectNextHyperplane0(unsigned long *excluded, rowrange *hnext)
{
  /*A natural way to choose the next hyperplane.  Simply the largest index*/
  long i;
  boolean determined;

  i = mm;
  determined = False;
  do {
    if (set_member(i, excluded))
      i--;
    else
      determined = True;
  } while (!determined && i>=1);
  if (determined) 
    *hnext = i;
  else
    *hnext = 0;
}

void SelectNextHyperplane1(unsigned long *excluded, rowrange *hnext)
{
  /*Natural way to choose the next hyperplane.  Simply the least index*/
  long i;
  boolean determined;

  i = 1;
  determined = False;
  do {
    if (set_member(i, excluded))
      i++;
    else
      determined = True;
  } while (!determined && i<=mm);
  if (determined) 
    *hnext = i;
  else 
    *hnext=0;
}

void SelectNextHyperplane2(unsigned long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmin, fi=0;   /*feasibility and infeasibility numbers*/

  infmin = RayCount + 1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (inf < infmin) {
	infmin = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (DynamicWriteOn) {
    cout << "*infeasible rays (min) = " << infmin << ",  #feas rays =" << fi << "\n";
  }
}

void SelectNextHyperplane3(unsigned long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with maximum infeasibility*/
  long i, fea, inf, infmax, fi=0;   /*feasibility and infeasibility numbers*/

  infmax = -1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (inf > infmax) {
	infmax = inf;
	fi = fea;
	*hnext = i;
      }
    }
  }
  if (DynamicWriteOn) {
    cout << "*infeasible rays (max) = " << infmax << ",  #feas rays =" << fi << "\n";
  }
}

void SelectNextHyperplane4(unsigned long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane with the most unbalanced cut*/
  long i, fea, inf, max, tmax, fi=0, infi=0;
      /*feasibility and infeasibility numbers*/

  max = -1;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      FeasibilityIndices(&fea, &inf, i);
      if (fea <= inf)
        tmax = inf;
      else
        tmax = fea;
      if (tmax > max) {
        max = tmax;
        fi = fea;
        infi = inf;
        *hnext = i;
      }
    }
  }
  if (!DynamicWriteOn)
    return;
  if (max == fi) {
    cout << "*infeasible rays (min) = " << infi << ",  #feas rays =" << fi << "\n";
  } else {
    cout << "*infeasible rays (max) = " << infi << ",  #feas rays =" << fi << "\n";
  }
}

void SelectNextHyperplane5(unsigned long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-min*/
  long i, minindex;
  myTYPE *v1, *v2;

  minindex = 0;
  v1 = NULL;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
	  v2 = AA[i - 1];
      if (minindex == 0) {
	    minindex = i;
	    v1=v2;
      } else if (LexSmaller(v2,v1,nn)) {
        minindex = i;
	    v1=v2;
      }
    }
  }
  *hnext = minindex;
}


void SelectNextHyperplane6(unsigned long *excluded, rowrange *hnext)
{
  /*Choose the next hyperplane which is lexico-max*/
  long i, maxindex;
  myTYPE *v1, *v2;

  maxindex = 0;
  v1 = NULL;
  for (i = 1; i <= mm; i++) {
    if (!set_member(i, excluded)) {
      v2= AA[i - 1];
      if (maxindex == 0) {
        maxindex = i;
        v1=v2;
      } else if (LexLarger(v2, v1, nn)) {
        maxindex = i;
        v1=v2;
     }
    }
  }
  *hnext = maxindex;
}

long Partition(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  myTYPE *x;
  long i,j,ovi;
  
  x=A[OV[p]-1];
  i=p-1;
  j=r+1;
  while (True){
    do{
      j--;
    } while (LexLarger(A[OV[j]-1],x,nmax));
    do{
      i++;
    } while (LexSmaller(A[OV[i]-1],x,nmax));
    if (i<j){
      ovi=OV[i];
      OV[i]=OV[j];
      OV[j]=ovi;
    }
    else{
      return j;
    }
  }
}

void QuickSort(rowindex OV, long p, long r, Amatrix A, long nmax)
{
  long q;
  
  if (p < r){
    q = Partition(OV, p, r, A, nmax);
    QuickSort(OV, p, q, A, nmax);
    QuickSort(OV, q+1, r, A, nmax);
  }
}


#ifndef RAND_MAX 
#define RAND_MAX 32767 
#endif

void RandomPermutation(rowindex OV, long t, unsigned int seed)
{
  long k,j,ovj;
  double u,xk,r,rand_max=(double) RAND_MAX;
  boolean localdebug=False;

  srand(seed);
  for (j=t; j>1 ; j--) {
    r=rand();
    u=r/rand_max;
    xk=j*u +1;
    k=(long)xk;
    if (localdebug) printf("u=%lg, k=%ld, r=%lg, randmax= %lg\n",u,k,r,rand_max);
    ovj=OV[j];
    OV[j]=OV[k];
    OV[k]=ovj;
    if (localdebug) printf("row %ld is exchanged with %ld\n",j,k); 
  }
}

void UpdateRowOrderVector(unsigned long *PriorityRows)
/* Update the RowOrder vector to shift selected rows
in highest order.
*/
{
  rowrange i,j,k,j1=0,oj=0;
  long rr;
  boolean found, localdebug=False;
  
  if (debug) localdebug=True;
  found=True;
  rr=set_card(PriorityRows);
  if (localdebug) set_write(PriorityRows);
  for (i=1; i<=rr; i++){
    found=False;
    for (j=i; j<=mm && !found; j++){
      oj=OrderVector[j];
      if (set_member(oj, PriorityRows)){
        found=True;
        if (localdebug) printf("%ldth in sorted list (row %ld) is in PriorityRows\n", j, oj);
        j1=j;
      }
    }
    if (found){
      if (j1>i) {
        /* shift everything lower: OV[i]->OV[i+1]..OV[j1-1]->OV[j1] */
        for (k=j1; k>=i; k--) OrderVector[k]=OrderVector[k-1];
        OrderVector[i]=oj;
        if (localdebug){
          printf("OrderVector updated to:\n");
          for (j = 1; j <= mm; j++) printf(" %2ld", OrderVector[j]);
          printf("\n");
        }
      }
    } else {
      printf("UpdateRowOrder: Error.\n");
      goto _L99;
    }
  }
_L99:;
}

void SelectPreorderedNext(unsigned long *excluded, rowindex OV, rowrange *hnext)
{
  rowrange i,k;
  
  *hnext=0;
  for (i=1; i<=mm && *hnext==0; i++){
    k=OV[i];
    if (!set_member(k, excluded)) *hnext=k ;
  }
}

void SelectNextHyperplane(HyperplaneOrderType ho, 
         unsigned long *excluded, rowrange *hh, boolean *RefreshOrderVector)
{
  if (PreOrderedRun){
    if (debug) {
      printf("debug SelectNextHyperplane: Use PreorderNext\n");
    }
    SelectPreorderedNext(excluded, OrderVector, hh);
  }
  else {
    if (debug) {
      printf("debug SelectNextHyperplane: Use DynamicOrderedNext\n");
    }

    switch (ho) {

    case MaxIndex:
      SelectNextHyperplane0(excluded,hh);
      break;

    case MinIndex:
      SelectNextHyperplane1(excluded,hh);
      break;

    case MinCutoff:
      SelectNextHyperplane2(excluded,hh);
      break;

    case MaxCutoff:
      SelectNextHyperplane3(excluded, hh);
      break;

    case MixCutoff:
      SelectNextHyperplane4(excluded, hh);
      break;

    default:
      SelectNextHyperplane0(excluded,hh);
      break;
    }
  }
}


myTYPE AValue(myTYPE *p, rowrange i)
{
  /*return the ith component of the vector  A x p */
  colrange j;
  myTYPE temp=0;

  temp = 0;
  for (j = 0; j < nn; j++)
    temp = temp + AA[i - 1][j] * p[j];
  return temp;
}

void StoreRay1(myTYPE *p, RayRecord *RR, boolean *feasible)
{  /* Original ray storing routine when RelaxedEnumeration is False */
  rowrange i,k,fii=mm+1;
  colrange j;
  myTYPE temp=0;
  boolean localdebug=False;

  if (debug) localdebug=True;
  *feasible = True;
  set_initialize(&(RR->ZeroSet),mm);
  if (localdebug) cout << "StoreRay1: RR->ARay will be set.\n"; 
  *(RR->ARay) = 0;
  if (localdebug) cout << "StoreRay1: RR->ARay is set to " << RR->ARay <<"\n"; 
  for (j = 0; j < nn; j++)
    RR->Ray[j] = p[j];
  if (localdebug) cout << "StoreRay1: RR->Ray[j] are set.\n"; 
  for (i = 1; i <= mm; i++) {
    k=OrderVector[i];
    temp = AValue(p, k);
    if (FABS(temp) <= zero)
      set_addelem(RR->ZeroSet, k);
    if (temp < -zero){
      *feasible = False;
      if (fii>mm) fii=i;  /* the first violating inequality index */
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
  if (localdebug) 
    if (fii<= mm)
      printf("StoreRay1:store ray with fii= %ld (row =%ld)\n", fii,OrderVector[fii]);
    else
      printf("StoreRay1:store ray with fii= %ld (feasible=%d)\n", fii,*feasible);    
}

void StoreRay2(myTYPE *p, RayRecord *RR, 
  boolean *feasible, boolean *weaklyfeasible)
   /* Ray storing routine when RelaxedEnumeration is True.
       weaklyfeasible is true iff it is feasible with
       the strict_inequality conditions deleted. */
{
  rowrange i,k,fii=mm+1;
  colrange j;
  myTYPE temp=0;
  boolean localdebug=False;

  if (debug) localdebug=True;
  *feasible = True;
  *weaklyfeasible = True;
  set_initialize(&(RR->ZeroSet),mm);
  *(RR->ARay) = 0;
  for (j = 0; j < nn; j++)
    RR->Ray[j] = p[j];
  for (i = 1; i <= mm; i++) {
    k=OrderVector[i];
    temp = AValue(p, k);
    if (FABS(temp) <= zero){
      set_addelem(RR->ZeroSet, k);
      if (EqualityIndex[k]==-1) 
        *feasible=False;  /* strict inequality required */
    }
    if (temp < -zero){
      *feasible = False;
      if (fii>mm && EqualityIndex[k]>=0) {
        fii=i;  /* the first violating inequality index */
        *weaklyfeasible=False;
      }
    }
  }
  RR->FirstInfeasIndex=fii;
  RR->feasible = *feasible;
  if (localdebug) {
    if (fii<= mm)
      printf("StoreRay2:store ray with fii= %ld (row =%ld)\n", fii,OrderVector[fii]);
    else
      printf("StoreRay2:store ray with fii= %ld (weaklyfeasible=%d)\n", fii,*weaklyfeasible);    
    if (*weaklyfeasible) WriteRayRecord(cout,RR);
  }
}



void ConditionalAddEdge(RayRecord *Ray1, RayRecord *Ray2, RayRecord *ValidFirstRay)
{
  long it,it_row,fii1,fii2,fmin,fmax;
  boolean adjacent,lastchance;
  RayRecord *TempRay,*Rmin,*Rmax;
  AdjacencyRecord *NewEdge;
  boolean localdebug=False;
  rowset ZSmin, ZSmax;
  
  fii1=Ray1->FirstInfeasIndex;
  fii2=Ray2->FirstInfeasIndex;
  if (fii1<fii2){
    fmin=fii1; fmax=fii2;
    Rmin=Ray1;
    Rmax=Ray2;
  }
  else{
    fmin=fii2; fmax=fii1;
    Rmin=Ray2;
    Rmax=Ray1;
  }
  ZSmin = Rmin->ZeroSet;
  ZSmax = Rmax->ZeroSet;
  if (localdebug) {
    printf("ConditionalAddEdge: FMIN = %ld (row%ld)   FMAX=%ld\n",
      fmin,OrderVector[fmin], fmax);
  }
  if (fmin==fmax){
    if (localdebug) printf("ConditionalAddEdge: equal FII value-> No edge added\n");
  }
  else if (set_member(OrderVector[fmin],ZSmax)){
    if (localdebug) printf("ConditionalAddEdge: No strong separation -> No edge added\n");
  }
  else {  /* the pair will be separated at the iteration fmin */
    lastchance=True;
    /* flag to check it will be the last chance to store the edge candidate */
    set_int(Face1, ZSmax, ZSmin);
    count_int++;
    if (localdebug){
      printf("Face: ");
      for (it=1; it<=mm; it++) {
        it_row=OrderVector[it];
        if (set_member(it_row, Face1)) printf("%ld ",it_row);
      }
      printf("\n");
    }
    for (it=Iteration+1; it < fmin && lastchance; it++){
      it_row=OrderVector[it];
      if (EqualityIndex[it_row]>=0 && set_member(it_row, Face1)){
        lastchance=False;
        count_int_bad++;
        if (localdebug){
          printf("There will be another chance iteration %ld (row %ld) to store the pair\n", it, it_row);
        }
      }
    }
    if (lastchance){
      adjacent = True;
      count_int_good++;
      /* adjacent checking */
      set_int(Face, Face1, AddedHyperplanes);
      if (localdebug){
        printf("Check adjacency\n");
        printf("AddedHyperplanes: "); set_write(AddedHyperplanes);
        printf("Face: ");
        for (it=1; it<=mm; it++) {
          it_row=OrderVector[it];
          if (set_member(it_row, Face)) printf("%ld ",it_row);
        }
        printf("\n");
      }
      if (set_card(Face)< nn - 2) {
        adjacent = False;
      }
      else if (NondegAssumed) {
    	adjacent = True;
      }
      else{
        TempRay = ValidFirstRay;  /* the first ray for adjacency checking */
        while (TempRay != NULL && adjacent) {
          if (TempRay != Ray1 && TempRay != Ray2) {
            set_int(Face1, TempRay->ZeroSet, AddedHyperplanes);
            if (set_subset(Face, Face1)) {
              if (localdebug) set_write(Face1);
              adjacent = False;
            }
          }
          TempRay = TempRay->Next;
        }
      }
      if (adjacent){
        if (localdebug) printf("The pair is adjacent and the pair must be stored for iteration %ld (row%ld)\n",
          fmin, OrderVector[fmin]);
        NewEdge=(struct AdjacencyRecord *) malloc(sizeof *NewEdge);
        NewEdge->Ray1=Rmax;  /* save the one remains in iteration fmin in the first */
        NewEdge->Ray2=Rmin;  /* save the one deleted in iteration fmin in the second */
        NewEdge->Next=NULL;
        EdgeCount++; 
        TotalEdgeCount++;
        if (Edges[fmin]==NULL){
          Edges[fmin]=NewEdge;
          if (localdebug) printf("Create a new edge list of %ld\n", fmin);
        }else{
          NewEdge->Next=Edges[fmin];
          Edges[fmin]=NewEdge;
        }
      }
    }
  }
}

void CreateInitialEdges(void)
{
  RayRecord *Ptr1, *Ptr2;
  rowrange fii1,fii2;
  long count=0;
  boolean localdebug=False;
  
  if (FirstRay ==NULL || LastRay==NULL){
    printf("Error found: CreateInitialEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=FirstRay;
  while(Ptr1!=LastRay && Ptr1!=NULL){
    fii1=Ptr1->FirstInfeasIndex;
    Ptr2=Ptr1->Next;
    while(Ptr2!=NULL){
      fii2=Ptr2->FirstInfeasIndex;
      count++;
      if (localdebug) printf("CreateInitialEdges: edge %ld \n",count);
      if (fii1!=fii2) ConditionalAddEdge(Ptr1,Ptr2,FirstRay);
      Ptr2=Ptr2->Next;
    }
    Ptr1=Ptr1->Next;
  }
_L99:;  
}


void UpdateEdges(RayRecord *RRbegin, RayRecord *RRend)
/* This procedure must be called after the ray list is sorted
   by EvaluateARay2 so that FirstInfeasIndex's are monotonically
   increasing.
*/
{
  RayRecord *Ptr1, *Ptr2begin, *Ptr2;
  rowrange fii1;
  boolean ptr2found,quit,localdebug=False;
  long count=0,pos1, pos2;
  float workleft,prevworkleft=110.0,totalpairs;

  totalpairs=(ZeroRayCount-1.0)*(ZeroRayCount-2.0)+1.0;
  Ptr2begin = NULL; 
  if (RRbegin ==NULL || RRend==NULL){
    if (1) printf("Warning: UpdateEdges called with NULL pointer(s)\n");
    goto _L99;
  }
  Ptr1=RRbegin;
  pos1=1;
  do{
    ptr2found=False;
    quit=False;
    fii1=Ptr1->FirstInfeasIndex;
    pos2=2;
    for (Ptr2=Ptr1->Next; !ptr2found && !quit; Ptr2=Ptr2->Next,pos2++){
      if  (Ptr2->FirstInfeasIndex > fii1){
        Ptr2begin=Ptr2;
        ptr2found=True;
      }
      else if (Ptr2==RRend) quit=True;
    }
    if (ptr2found){
      quit=False;
      for (Ptr2=Ptr2begin; !quit ; Ptr2=Ptr2->Next){
        count++;
        if (localdebug) printf("UpdateEdges: edge %ld \n",count);
        ConditionalAddEdge(Ptr1,Ptr2,RRbegin);
        if (Ptr2==RRend || Ptr2->Next==NULL) quit=True;
      }
    }
    Ptr1=Ptr1->Next;
    pos1++;
    workleft = 100.0 * (ZeroRayCount-pos1) * (ZeroRayCount - pos1-1.0) / totalpairs;
    if (ZeroRayCount>=500 && DynamicWriteOn && pos1%10 ==0 && prevworkleft-workleft>=10 ) {
      cout << "*Work of iteration " << Iteration << "(/" << mm << "):" 
         << pos1 << "/" << ZeroRayCount << " => ";
      cout.width(4); cout << workleft; cout << "% left\n";
      prevworkleft=workleft;
    }    
  }while(Ptr1!=RRend && Ptr1!=NULL);
_L99:;  
}

void FreeDDMemory(void)
{
  RayRecord *Ptr, *PrevPtr;
  long count;
  boolean localdebug=False;
  
  PrevPtr=ArtificialRay;
  count=0;
  for (Ptr=ArtificialRay->Next; Ptr!=NULL; Ptr=Ptr->Next){
    if (!PostAnalysisOn) delete[] PrevPtr->Ray;
    delete PrevPtr->ARay;
    free(PrevPtr->ZeroSet);
    delete PrevPtr;
    count++;
    PrevPtr=Ptr;
  };
  LastRay=NULL;
  FirstRay=NULL;
  ArtificialRay=NULL;
  if (localdebug) printf("%ld ray storage spaces freed\n",count);
  
  set_free(&InitialHyperplanes);
  set_free(&AddedHyperplanes);
  set_free(&GroundSet);
  set_free(&Face);
  set_free(&Face1);
  free(OrderVector);
  
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  WeaklyFeasibleRayCount = 0;
  VertexCount = 0;
}

void Normalize(myTYPE *V)
{
  long j;
  myTYPE min=0, temp=0;

  min = 1.0e+20;
  for (j = 0; j < nn; j++) {
    temp = FABS(V[j]);
    if (temp > zero && temp < min)
      min = temp;
  }
  for (j = 0; j < nn; j++)
    V[j] /= min;
}


void ZeroIndexSet(myTYPE *x, rowset ZS)
{
  rowrange i;
  myTYPE temp=0;

  set_emptyset(ZS);
  for (i = 1; i <= mm; i++) {
    temp = AValue(x, i);
    if (FABS(temp) <= zero)
      set_addelem(ZS, i);
  }
}


void FindInitialRays(rowset InitHyperplanes,
			    Bmatrix InitRays, colindex PivRow, boolean *found)
{
  Bmatrix BInverse;
  rowset CandidateRows;
  long i, j, rank;
  HyperplaneOrderType roworder=LexMin;

  *found = False;
  set_initialize(&CandidateRows, mm);
  if (InitBasisAtBottom==True) {
    roworder=MaxIndex;
    PreOrderedRun=False;
  }
  else PreOrderedRun=True;
  for (i = 1; i <= mm; i++)
    if (!set_member(i,NonequalitySet)) set_addelem(CandidateRows, i);
    /*all rows not in NonequalitySet are candidates for initial cone*/
  if (DynamicWriteOn)
    printf("*Computing an initial set of rays\n");
  InitializeBmatrix(BInverse);
  FindBasis(AA, roworder, InitHyperplanes, PivRow, BInverse, &rank);
  if (debug) {
    printf("FindInitialBasis: InitHyperplanes=");
    set_write(InitHyperplanes);
    for (j = 1; j <= nn; j++)
      printf("Pivotrow[%ld] = %ld \n", j, PivRow[j]);
    printf("nn = %ld, rank = %ld\n",nn,rank);
  }

/* The case of rank==nn-1 is not treated properly and deleted
  on Jun 8, 1996.  Check cddbag960531
*/
/*if (rank < nn-1) {  */
  if (rank <= nn-1) {
    if (debug) WriteBmatrix(cout,BInverse);
    Error = LowColumnRank;
    return;
  }
  if (!set_subset(EqualitySet,InitHyperplanes)) {
    Error = DependentMarkedSet;
    return;
  }
  *found = True;
  if (debug) {
    WriteBmatrix(cout,BInverse);
  }
  /* free_Bmatrix(InitRays);  */
  CopyBmatrix(BInverse,InitRays);
  if (debug) WriteBmatrix(cout,InitRays);
  set_free(&CandidateRows);
  if (HyperplaneOrder==MaxCutoff||HyperplaneOrder==MinCutoff||HyperplaneOrder==MixCutoff){
    PreOrderedRun=False;
  } else PreOrderedRun=True;
}

void CheckEquality(RayRecord **RP1, RayRecord **RP2, boolean *equal)
{
  long j;
  myTYPE two=2;

  if (debug)
    printf("Check equality of two rays\n");
  *equal = True;
  j = 1;
  while (j <= nn && *equal) {
    if (FABS((*RP1)->Ray[j - 1] - (*RP2)->Ray[j - 1]) > two * zero)
      *equal = False;
    j++;
  }
  if (*equal)
    printf("Equal records found !!!!\n");
}

void CreateNewRay(RayRecord *Ptr1, RayRecord *Ptr2, rowrange ii)
{
  /*Create a new ray by taking a linear combination of two rays*/
  colrange j;
  myTYPE v1=0, v2=0;
  myTYPE *vec1;

  vec1 = new myTYPE[nn];
  v1 = FABS(AValue(Ptr1->Ray, ii));
  v2 = FABS(AValue(Ptr2->Ray, ii));
  for (j = 0; j < nn; j++)
    vec1[j] = Ptr1->Ray[j] * v2 + Ptr2->Ray[j] * v1;
  Normalize(vec1);
  AddRay(vec1);
  if (debug){
    WriteRayRecord(cout,Ptr1);
    WriteRayRecord(cout,Ptr2);  
    printf("create a new ray by eliminating %ld:\n",ii);
    WriteRayRecord(cout,LastRay);
  }
  delete[] vec1;
}

void EvaluateARay1(rowrange i)
/* Evaluate the ith component of the vector  A x RD.Ray 
    and rearrange the linked list so that
    the infeasible rays with respect to  i  will be
    placed consecutively from First 
 */
{
  colrange j;
  myTYPE temp=0;
  RayRecord *Ptr, *PrevPtr, *TempPtr;

  Ptr = FirstRay;
  PrevPtr = ArtificialRay;
  if (PrevPtr->Next != Ptr) {
    printf("Error.  Artificial Ray does not point to FirstRay!!!\n");
  }
  while (Ptr != NULL) {
    temp = 0;
    for (j = 0; j < nn; j++)
      temp = temp + AA[i - 1][j] * Ptr->Ray[j];
    *(Ptr->ARay) = temp;
    // if ( temp <= -zero && Ptr != FirstRay) {
    if ( temp < -zero && Ptr != FirstRay) {
      /* printf("Moving an infeasible record w.r.t. %ld to FirstRay\n",i); */
      if (Ptr==LastRay) LastRay=PrevPtr;
      TempPtr=Ptr;
      Ptr = Ptr->Next;
      PrevPtr->Next = Ptr;
      ArtificialRay->Next = TempPtr;
      TempPtr->Next = FirstRay;
      FirstRay = TempPtr;
    }
    else {
      PrevPtr = Ptr;
      Ptr = Ptr->Next;
    }
  }
}

void EvaluateARay2(rowrange i)
/* Evaluate the ith component of the vector  A x RD.Ray 
   and rearrange the linked list so that
   the infeasible rays with respect to  i  will be
   placed consecutively from First. Also for all feasible rays,
   "positive" rays and "zero" rays will be placed consecutively.
 */
{
  colrange j;
  myTYPE temp=0;
  RayRecord *Ptr, *NextPtr;
  boolean zerofound=False,negfound=False,posfound=False;

  PosHead=NULL;ZeroHead=NULL;NegHead=NULL;
  PosLast=NULL;ZeroLast=NULL;NegLast=NULL;
  Ptr = FirstRay;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    Ptr->Next=NULL;     /* then clear the Next pointer */
    temp = 0;
    for (j = 0; j < nn; j++)
      temp = temp +  AA[i - 1][j] * Ptr->Ray[j];
    *(Ptr->ARay) = temp;
    if ( temp < -zero) {
      if (!negfound){
        negfound=True;
        NegHead=Ptr;
        NegLast=Ptr;
      }
      else{
        Ptr->Next=NegHead;
        NegHead=Ptr;
      }
    }
    else if (temp > zero){
      if (!posfound){
        posfound=True;
        PosHead=Ptr;
        PosLast=Ptr;
      }
      else{  
        Ptr->Next=PosHead;
        PosHead=Ptr;
       }
    }
    else {
      if (!zerofound){
        zerofound=True;
        ZeroHead=Ptr;
        ZeroLast=Ptr;
      }
      else{
        Ptr->Next=ZeroHead;
        ZeroHead=Ptr;
      }
    }
    Ptr=NextPtr;
  }
  /* joining three neg, pos and zero lists */
  if (negfound){                 /* -list nonempty */
    FirstRay=NegHead;
    if (posfound){               /* -list & +list nonempty */
      NegLast->Next=PosHead;
      if (zerofound){            /* -list, +list, 0list all nonempty */
        PosLast->Next=ZeroHead;
        LastRay=ZeroLast;
      } 
      else{                      /* -list, +list nonempty but  0list empty */
        LastRay=PosLast;      
      }
    }
    else{                        /* -list nonempty & +list empty */
      if (zerofound){            /* -list,0list nonempty & +list empty */
        NegLast->Next=ZeroHead;
        LastRay=ZeroLast;
      } 
      else {                      /* -list nonempty & +list,0list empty */
        LastRay=NegLast;
      }
    }
  }
  else if (posfound){            /* -list empty & +list nonempty */
    FirstRay=PosHead;
    if (zerofound){              /* -list empty & +list,0list nonempty */
      PosLast->Next=ZeroHead;
      LastRay=ZeroLast;
    } 
    else{                        /* -list,0list empty & +list nonempty */
      LastRay=PosLast;
    }
  }
  else{                          /* -list,+list empty & 0list nonempty */
    FirstRay=ZeroHead;
    LastRay=ZeroLast;
  }
  ArtificialRay->Next=FirstRay;
  LastRay->Next=NULL;
}

void DeleteNegativeRays(void)
/* Eliminate the infeasible rays with respect to  i  which
   are supposed to be consecutive from the head of the RayRecord list,
   and sort the zero list assumed to be consecutive at the
   end of the list.
 */
{
  rowrange fii,fiitest;
  myTYPE temp=0;
  RayRecord *Ptr, *PrevPtr, *NextPtr, *ZeroPtr1, *ZeroPtr0;
  boolean found, completed, zerofound=False,negfound=False,posfound=False;
  boolean localdebug=False;
  
  PosHead=NULL;ZeroHead=NULL;NegHead=NULL;
  PosLast=NULL;ZeroLast=NULL;NegLast=NULL;

  /* Delete the infeasible rays  */
  PrevPtr= ArtificialRay;
  Ptr = FirstRay;
  if (PrevPtr->Next != Ptr) 
    printf("Error at DeleteNegativeRays: ArtificialRay does not point the FirstRay.\n");
  completed=False;
  while (Ptr != NULL && !completed){
    if ( *(Ptr->ARay) < -zero ){
      Eliminate(&PrevPtr);
      Ptr=PrevPtr->Next;
    }
    else{
      completed=True;
    }
  }
  
  /* Sort the zero rays */
  Ptr = FirstRay;
  ZeroRayCount=0;
  while (Ptr != NULL) {
    NextPtr=Ptr->Next;  /* remember the Next record */
    temp = *(Ptr->ARay);
    if (localdebug) cout << "Ptr->ARay :" << temp << "\n";
    if ( temp < -zero) {
      if (!negfound){
        printf("Error: An infeasible ray found after their removal\n");
        negfound=True;
      }
    }
    else if (temp > zero){
      if (!posfound){
        posfound=True;
        PosHead=Ptr;
        PosLast=Ptr;
      }
      else{  
        PosLast=Ptr;
       }
    }
    else {
      ZeroRayCount++;
      if (!zerofound){
        zerofound=True;
        ZeroHead=Ptr;
        ZeroLast=Ptr;
        ZeroLast->Next=NULL;
      }
      else{/* Find a right position to store the record sorted w.r.t. FirstInfeasIndex */
        fii=Ptr->FirstInfeasIndex; 
        found=False;
        ZeroPtr1=NULL;
        for (ZeroPtr0=ZeroHead; !found && ZeroPtr0!=NULL ; ZeroPtr0=ZeroPtr0->Next){
          fiitest=ZeroPtr0->FirstInfeasIndex;
          if (fiitest >= fii){
            found=True;
          }
          else ZeroPtr1=ZeroPtr0;
        }
        /* printf("insert position found \n %d  index %ld\n",found, fiitest); */
        if (!found){           /* the new record must be stored at the end of list */
          ZeroLast->Next=Ptr;
          ZeroLast=Ptr;
          ZeroLast->Next=NULL;
        }
        else{
          if (ZeroPtr1==NULL){ /* store the new one at the head, and update the head ptr */
            /* printf("Insert at the head\n"); */
            Ptr->Next=ZeroHead;
            ZeroHead=Ptr;
          }
          else{                /* store the new one inbetween ZeroPtr1 and 0 */
            /* printf("Insert inbetween\n");  */
            Ptr->Next=ZeroPtr1->Next;
            ZeroPtr1->Next=Ptr;
          }
        }
        /*
        Ptr->Next=ZeroHead;
        ZeroHead=Ptr;
        */
      }
    }
    Ptr=NextPtr;
  }
  /* joining the pos and zero lists */
  if (posfound){            /* -list empty & +list nonempty */
    FirstRay=PosHead;
    if (zerofound){              /* +list,0list nonempty */
      PosLast->Next=ZeroHead;
      LastRay=ZeroLast;
    } 
    else{                        /* 0list empty & +list nonempty */
      LastRay=PosLast;
    }
  }
  else{                          /* +list empty & 0list nonempty */
    FirstRay=ZeroHead;
    LastRay=ZeroLast;
  }
  ArtificialRay->Next=FirstRay;
  LastRay->Next=NULL;
}

void FeasibilityIndices(long *fnum, long *infnum, rowrange i)
{
  /*Evaluate the number of feasible rays and infeasible rays*/
  /*  w.r.t the hyperplane  i*/
  colrange j;
  myTYPE temp=0;
  RayRecord *Ptr;

  *fnum = 0;
  *infnum = 0;
  Ptr = FirstRay;
  while (Ptr != NULL) {
    temp = 0;
    for (j = 0; j < nn; j++)
      temp = temp + AA[i - 1][j] * Ptr->Ray[j];
    if (temp >= -zero)
      (*fnum)++;
    else
      (*infnum)++;
    Ptr = Ptr->Next;
  }
}

boolean LexSmaller(myTYPE *v1, myTYPE *v2, long nmax)
{ /* nmax is the size of vectors v1,v2 */
  boolean determined, smaller;
  colrange j;

  smaller = False;
  determined = False;
  j = 1;
  do {
    if (FABS(v1[j - 1]-v2[j - 1])>zero) {
      if (v1[j - 1] < v2[j - 1]) {
	    smaller = True;
	  }
      determined = True;
    } else
      j++;
  } while (!(determined) && (j <= nmax));
  return smaller;
}

boolean LexLarger(myTYPE *v1, myTYPE *v2, long nmax)

{ /* nmax is the size of vectors v1,v2 */
  boolean determined, larger;
  colrange j;

  larger = False;
  determined = False;
  j = 1;
  do {
    if (FABS(v1[j - 1]-v2[j - 1])>zero) {
      if (v1[j - 1] > v2[j - 1]) {
	    larger = True;
	  }
      determined = True;
    } else
      j++;
  } while (!(determined) && (j <= nmax));
  return larger;
}

void CopyArow(myTYPE *vcopy, myTYPE *v, long nmax)
{
 colrange j;

  for (j = 1; j <= nmax; j++) {
    vcopy[j-1] = v[j-1];
  }
}

void AddNewHyperplane1(rowrange hnew)
/* This procedure 1 must be used with PreorderedRun=False 
   This procedure is the most elementary implementation of
   DD and can be used with any type of ordering, including
   dynamic ordering of rows, e.g. MaxCutoff, MinCutoff.
   The memory requirement is minimum because it does not
   store any adjacency among the rays.
*/
{
  RayRecord *RayPtr0, *RayPtr1, *RayPtr2, *RayPtr2s, *RayPtr3;
  long pos1, pos2;
  double prevprogress, progress;
  myTYPE value1=0, value2=0;
  boolean adj, equal, completed;

  EvaluateARay1(hnew);  
   /*Check feasibility of rays w.r.t. hnew 
     and put all infeasible ones consecutively */

  RayPtr0 = ArtificialRay;   /*Pointer pointing RayPrt1*/
  RayPtr1 = FirstRay;        /*1st hnew-infeasible ray to scan and compare with feasible rays*/
  value1 = *(FirstRay->ARay);
  if (value1 > -zero) {
    if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
    goto _L99;        /* Sicne there is no hnew-infeasible ray and nothing to do */
  }
  else {
    RayPtr2s = RayPtr1->Next;/* RayPtr2s must point the first feasible ray */
    pos2=1;
    while (RayPtr2s!=NULL && *(RayPtr2s->ARay) <= -zero) {
      RayPtr2s = RayPtr2s->Next;
      pos2++;
    }
  }
  if (RayPtr2s==NULL) {
    FirstRay=NULL;
    ArtificialRay->Next=FirstRay;
    RayCount=0;
    VertexCount=0;
    if (Inequality==ZeroRHS && Conversion==IneToExt){
      CompStatus=AllFound;
    } else {
      CompStatus=RegionEmpty;
    }
     goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  RayPtr2 = RayPtr2s;   /*2nd feasible ray to scan and compare with 1st*/
  RayPtr3 = LastRay;    /*Last feasible for scanning*/
  prevprogress=-10.0;
  pos1 = 1;
  completed=False;
  while ((RayPtr1 != RayPtr2s) && !completed) {
    value1 = *(RayPtr1->ARay);
    value2 = *(RayPtr2->ARay);
    if (debug) {
      WriteRayRecord2(cout, RayPtr1);
      WriteRayRecord2(cout, RayPtr2);
      cout << "check feasibility " <<  value1 << "  " << value2 << "\n";
    }
    CheckEquality(&RayPtr1, &RayPtr2, &equal);
    if ((value1 > zero && value2 < -zero) || (value2 > zero && value1 < -zero)) {
      switch (AdjacencyTest) {

      case Algebraic:
		CheckAdjacency1(&RayPtr1, &RayPtr2, &adj);
		break;

      case Combinatorial:
		CheckAdjacency2(&RayPtr1, &RayPtr2, &adj);
		break;
      }
      if (adj)
		CreateNewRay(RayPtr1, RayPtr2, hnew);
    }
    if (RayPtr2 != RayPtr3) {
      RayPtr2 = RayPtr2->Next;
      continue;
    }
    if (value1 < -zero || equal) {
      Eliminate(&RayPtr0);
      RayPtr1 = RayPtr0->Next;
      RayPtr2 = RayPtr2s;
    } else {
      completed=True;
    }
    pos1++;
    progress = 100.0 * ((double)pos1 / pos2) * (2.0 * pos2 - pos1) / pos2;
    if (progress-prevprogress>=10 && pos1%10==0 && DynamicWriteOn) {
      cout << "*Progress of iteration " << Iteration << "(/" << mm << "):" <<
        pos1 << "/" << pos2 << " => ";
      cout.width(4); cout << progress; cout << "% done\n";
      prevprogress=progress;
    }
  }
  if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
  _L99:;
}

void AddNewHyperplane2(rowrange hnew)
/* This procedure must be used under PreOrderedRun mode */
{
  RayRecord *RayPtr0,*RayPtr1, *RayPtr2;
  AdjacencyRecord *EdgePtr, *EdgePtr0;
  long pos1;
  rowrange fii1, fii2;
  boolean localdebug=False;

  EvaluateARay2(hnew);
   /* Check feasibility of rays w.r.t. hnew 
      and sort them. ( -rays, +rays, 0rays)*/

  if (PosHead==NULL && ZeroHead==NULL) {
    FirstRay=NULL;
    ArtificialRay->Next=FirstRay;
    RayCount=0;
    VertexCount=0;
    if (Inequality==ZeroRHS && Conversion==IneToExt){
      CompStatus=AllFound;
    } else {
      CompStatus=RegionEmpty;
    }
    goto _L99;   /* All rays are infeasible, and the computation must stop */
  }
  
  if (localdebug){
    pos1=0;
    printf("(pos, FirstInfeasIndex, A Ray)=\n");
    for (RayPtr0=FirstRay; RayPtr0!=NULL; RayPtr0=RayPtr0->Next){
      pos1++;
      printf("(%ld,%ld,",pos1,RayPtr0->FirstInfeasIndex);
      WriteNumber(cout,*(RayPtr0->ARay)); 
      printf(") ");
   }
    printf("\n");
  }
  
  if (ZeroHead==NULL) ZeroHead=LastRay;

  EdgePtr=Edges[Iteration];
  while (EdgePtr!=NULL){
    RayPtr1=EdgePtr->Ray1;
    RayPtr2=EdgePtr->Ray2;
    fii1=RayPtr1->FirstInfeasIndex;   
    CreateNewRay(RayPtr1, RayPtr2, hnew);
    fii2=LastRay->FirstInfeasIndex;
    if (fii1 != fii2) ConditionalAddEdge(RayPtr1,LastRay,PosHead);
    EdgePtr0=EdgePtr;
    EdgePtr=EdgePtr->Next;
    free(EdgePtr0);
    EdgeCount--;
  }
  Edges[Iteration]=NULL;
  
  DeleteNegativeRays();
    
  set_addelem(AddedHyperplanes, hnew);

  if (Iteration<mm){
    if (ZeroHead!=NULL && ZeroHead!=LastRay){
      if (ZeroRayCount>200 && DynamicWriteOn) printf("*New edges being scanned...\n");
      UpdateEdges(ZeroHead,LastRay);
    }
    if (localdebug) printf("*Edges currently stored = %ld\n", EdgeCount);
  }

  if (RayCount==WeaklyFeasibleRayCount) CompStatus=AllFound;
_L99:;
}


void AddRay(myTYPE *p)
{  
  boolean feasible, weaklyfeasible;

  if (FirstRay == NULL) {
    FirstRay = new RayRecord;
    FirstRay->Ray = new myTYPE[nn];
    FirstRay->ARay = new myTYPE;
    if (debug)
      printf("Create the first ray pointer\n");
    LastRay = FirstRay;
    ArtificialRay->Next = FirstRay;
  } else {
    LastRay->Next = new RayRecord;
    LastRay->Next->Ray = new myTYPE[nn];
    LastRay->Next->ARay = new myTYPE;
    if (debug)
      printf("Create a new ray pointer\n");
    LastRay = LastRay->Next;
  }
  LastRay->Next = NULL;
  RayCount++;
  TotalRayCount++;
  if (RelaxedEnumeration){
    StoreRay2(p, LastRay, &feasible, &weaklyfeasible);
    if (weaklyfeasible) WeaklyFeasibleRayCount++;
  } else {
    StoreRay1(p, LastRay, &feasible);
    if (feasible) WeaklyFeasibleRayCount++;
    /* weaklyfeasible is equiv. to feasible in this case. */
  }
  if (!feasible) return;
  else {
    FeasibleRayCount++;
    if (FABS(LastRay->Ray[0]) > zero && Inequality == NonzeroRHS)
      VertexCount++;
    if (DynamicRayWriteOn) {
      WriteRayRecord(cout,LastRay);
    }
  }
  if (DynamicWriteOn) {
    if (TotalRayCount % 100 == 0) {
      cout << "*Rays (Total, Currently Active, Feasible) = " 
        << TotalRayCount << " " << RayCount << " " << FeasibleRayCount << "\n";
    }
  }
  if (PostAnalysisOn){
       delete[] LastRay->Ray;
  }
}

void AddArtificialRay(void)
{  
  long j;
  boolean feasible;

  if (ArtificialRay != NULL) {
    printf("Warning !!!  FirstRay in not nil.  Illegal Call\n");
    return;
  }
  ArtificialRay = new RayRecord;
  ArtificialRay->Ray = new myTYPE[nn];
  ArtificialRay->ARay = new myTYPE;
  if (debug)
    printf("Create the artificial ray pointer\n");
  myTYPE *zerovec;
  zerovec = new myTYPE[nn];
  for (j = 0; j < nn; j++)
    zerovec[j] = 0;
  StoreRay1(zerovec, ArtificialRay, &feasible);
  ArtificialRay->Next = NULL;
}

void ComputeRowOrderVector(rowindex OV, HyperplaneOrderType ho)
{
  long i,itemp,j;
 
  OV[0]=0;
  switch (ho){
  case MaxIndex:
    for(i=1; i<=mm; i++) OV[i]=mm-i+1;
    break;

  case MinIndex: 
    for(i=1; i<=mm; i++) OV[i]=i;
    break;

  case LexMin: case MinCutoff: case MixCutoff: case MaxCutoff:
    for(i=1; i<=mm; i++) OV[i]=i;
    QuickSort(OV, 1, mm, AA, nn);
    break;

  case LexMax:
    for(i=1; i<=mm; i++) OV[i]=i;
    QuickSort(OV, 1, mm, AA, nn);
    for(i=1; i<=mm/2;i++){   /* just reverse the order */
      itemp=OV[i];
      OV[i]=OV[mm-i+1];
      OV[mm-i+1]=itemp;
    }
    break;

  case RandomRow:
    for(i=1; i<=mm; i++) OV[i]=i;
    if (rseed<=0) rseed=1;
    RandomPermutation(OV, mm, rseed);
    break;

  case LineShelling:
    myTYPE *vec1, *vec2;
    vec1 = new myTYPE[nn];
    vec2 = new myTYPE[nn]; 
    for(i=1; i<=mm; i++) OV[i]=i;
    vec1[0]=1;
    vec2[0]=0;
    if (rseed<=0) rseed=1;
    srand(rseed);
    for(j=2; j<=nn; j++){
      vec1[j-1]=0;
      vec2[j-1]=nn-j+1;
      /* vec2[j-1]=rand(); */
    }
    LineShellingOrder(OV, vec1, vec2);
    delete[] vec1;
    delete[] vec2;
    break;
  }
}

void LineShellingOrder(rowindex OV, myTYPE *z, myTYPE *d)
/* find the shelling ordering induced by a point 
   z (interior point, i.e. A z > 0) and a direction vector  d */
{
  long i,j;
  myTYPE temp1=0,temp2=0,infinity=(myTYPE)1e+20;
  boolean localdebug=False;
  myTYPE* beta;
  
  beta=new myTYPE[mm];
  for (i=1; i<= mm; i++) beta[i-1]=AA[i-1][0]; /* store the first column in beta */
  for (i=1; i<= mm; i++){
    temp1 = 0;
    temp2 = 0;
    for (j = 1; j <= nn; j++){
      temp1 = temp1 + AA[i - 1][j-1] * z[j-1];
      temp2 = temp2 + AA[i - 1][j-1] * d[j-1];
    }
    if (FABS(temp1)>zero) AA[i-1][0]=temp2/temp1;  
    else if (temp1*temp2 > zero) AA[i-1][0]= infinity;
    else AA[i-1][0]= -infinity;
     /* use the first column of AA tentatively */
  }
  if (localdebug) 
    for (i=1; i<= mm; i++){
      cout << "set AA[" << i << "] =" << AA[i-1][0] << "\n";
    }
  QuickSort(OV, 1, mm, AA, 1);
  for (i=1; i<= mm; i++) {
    AA[i-1][0]=beta[i-1]; 
     /* restore the first column of AA */ 
    if (localdebug) cout << "restore AA[" << i << "] =" << AA[i-1][0] << "\n";;
  }
}

/* end of cddarith.C */
