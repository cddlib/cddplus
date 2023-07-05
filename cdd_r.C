/* cdd.C: Main program of the sofware cdd+
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd+ : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

/*  This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
*/

/* The first version C0.21 was created on November 10,1993 
   with Dave Gillespie's p2c translator 
   from the Pascal program pdd.p written by Komei Fukuda. 
*/


#include <fstream>
#include <string>
using namespace std;

#include "cddtype.h"
#include "cddrevs.h"

extern "C" {
#include "setoper.h" 
  /* set operation library header (May 14, 1995 version or later) */
#include "cdddef.h"
#include "cdd.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

/* #include <profile.h>    THINK C PROFILER */
/* #include <console.h>    THINK C PROFILER */
} /* end of extern "C"  */


long minput, ninput;   /*size of input data [b -A] */
long mm, nn;   /*size of the homogenous system to be solved by dd*/
long projdim;  /*dimension of orthogonal preprojection */
colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1, SubproblemRowSet, RedundantRowSet;
rowrange Iteration, hh;
rowindex OrderVector;  /* the permutation vector to store a preordered row indeces */
rowindex EqualityIndex;  
  /* ith component is 1 if it is equality, -1 if it is strict inequality, 0 otherwise. */
rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
  TotalRayCount, VertexCount, ZeroRayCount;
long EdgeCount, TotalEdgeCount;
long count_int=0,count_int_good=0,count_int_bad=0;
boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, 
  ShowSignTableauOn=False, 
  PostAnalysisOn=False, CondensedListOn=False, SignPivotOn=False, 
  ManualPivotOn=False, debug;
Amatrix AA;
Bmatrix InitialRays;
colindex InitialRayIndex; /* 0 if the corr. ray is for generator of an extreme line */ 
colrange RHScol;   /* LP RHS column */
rowrange OBJrow;   /* LP OBJ row */
LPStatusType LPStatus;
Arow LPcost;  /* LP cost vector to be maximized  */
RayRecord *ArtificialRay, *FirstRay, *LastRay;
RayRecord *PosHead, *ZeroHead, *NegHead, *PosLast, *ZeroLast, *NegLast;
AdjacencyRecord *Edges[MMAX];  /* adjacency relation storage for iteration k */
boolean RecomputeRowOrder, found, inputsuccessful;
HyperplaneOrderType HyperplaneOrder;
AdjacencyTestType AdjacencyTest;
NumberType Number;
string InputNumberString, OutputNumberString;
InequalityType Inequality;
boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
boolean RestrictedEnumeration; /* Restricted enumeration switch (True if it is restricted on the intersection of EqualitySet hyperplanes) */
boolean RelaxedEnumeration; /* Relaxed enumeration switch (True if NonequalitySet inequalities must be satisfied with strict inequality) */
boolean VerifyInput; /* Verification switch for the input data */
boolean PreOrderedRun; 
  /* True if the rows are ordered before execution & all necessary adjacencies are stored */
CompStatusType CompStatus;  /* Computation Status */
ConversionType Conversion;
SpecialConversionType SpecialConversion;  /* rowdecomposition, rowsubproblem */
LPsolverType LPsolver;
IncidenceOutputType IncidenceOutput;
AdjacencyOutputType AdjacencyOutput;
ErrorType Error;
FileInputModeType FileInputMode;
DataFileType inputfile,ifilehead,ifiletail,
  outputfile,projfile,incfile,adjfile,logfile,dexfile,verfile,xtnfile;
time_t starttime, endtime;
unsigned int rseed=1;  /* random seed for random row permutation */

myTYPE zero=ZERO;      /*Rational or floating zero*/
boolean Round_Output;    /* rounding option for floating-point output. */
int output_digits=OUTPUTDIGITS;  /* Float digits for output.  Does not affect the computation. */

void DefaultOptionSetup(void)
{
  debug = False;
  DynamicWriteOn = True;
  DynamicRayWriteOn = True;
  LogWriteOn = False;
  HyperplaneOrder = LexMin;
  AdjacencyTest = Combinatorial;
  NondegAssumed = False;
  RecomputeRowOrder=True;
  PreOrderedRun=True;
  VerifyInput=False;
  Conversion = IneToExt;
  SpecialConversion = Nothing;
  LPsolver = DualSimplex;   
 
  IncidenceOutput = IncOff;
  AdjacencyOutput = AdjOff;
  InitBasisAtBottom = False;
  Round_Output=True;
}

void DDMain(ostream &f,ostream &f_log)
{
  Iteration = nn + 1;
  while (Iteration <= mm) {
    SelectNextHyperplane(HyperplaneOrder, WeaklyAddedHyperplanes, &hh, &RecomputeRowOrder);
    if (DynamicWriteOn) {
      cout << "*----------  Iteration = " << Iteration << " :   add  row # " << hh << "  ----------\n";
    }
    if (set_member(hh,NonequalitySet)){  /* Skip the row hh */
      if (DynamicWriteOn) {
        cout << "*The row # " << hh << " should be inactive and thus skipped.\n";
      }
      set_addelem(WeaklyAddedHyperplanes, hh);
    } else {
      if (PreOrderedRun)
        AddNewHyperplane2(hh);
      else
        AddNewHyperplane1(hh);
      set_addelem(AddedHyperplanes, hh);
      set_addelem(WeaklyAddedHyperplanes, hh);
    }
    (f_log) << Iteration << " " <<  hh << " " << TotalRayCount << " " <<
      RayCount << " " << FeasibleRayCount << "\n";
    if (CompStatus==AllFound||CompStatus==RegionEmpty) {
      set_addelem(AddedHyperplanes, hh);
      goto _L99;
    }
    Iteration++;
  }
  _L99:;
}

void Initialization(int ARGC, char *ARGV[])
/* Initialization of global variables */
{
  long i;

  Error=None;

  CompStatus=InProgress;
  if (ARGC>1){
    FileInputMode=Auto;
    strcpy(inputfile,ARGV[1]);
  }
  else{
    FileInputMode=Manual;
  }
}

void DDEnumerate(ostream &f, ostream &f_log)
{
  DDInit();
  time(&starttime);
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting(f, f_log);
    DDMain(f, f_log);
    WriteExtFile(f, f_log);
    if (IncidenceOutput == OutputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'i', "incidence");
      ofstream writing_ocd(incfile);
      WriteIncidenceFile(writing_ocd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_ocd.close();
    }
    if (IncidenceOutput == InputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'n', "input_incidence");
      ofstream writing_icd(incfile);
      WriteInputIncidenceFile(writing_icd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_icd.close();
    }
    if (AdjacencyOutput == OutputAdjacency || AdjacencyOutput == IOAdjacency){
      SetWriteFileName(adjfile, 'a', "adjacency");
      ofstream writing_adj(adjfile);
      if (DynamicWriteOn) printf("Writing the adjacency file %s...\n",adjfile);
      WriteAdjacencyFile(writing_adj);
      writing_adj.close();
      if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
    }
    if (AdjacencyOutput == InputAdjacency || AdjacencyOutput == IOAdjacency){
      SetWriteFileName(adjfile, 'j', "input_adjacency");
      ofstream writing_iad(adjfile);
      if (DynamicWriteOn) printf("Writing the input_adjacency file %s...\n",adjfile);
      WriteInputAdjacencyFile(writing_iad);
      writing_iad.close();
      if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
    }
    FreeDDMemory();
  } else {
    WriteExtFile(f, f_log);
    WriteErrorMessages(f);
  }
}

void DecompositionCore(ostream &f, ostream &f_log)
{
  DDInit();
  time(&starttime);
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting(f, f_log);
    DDMain(f, f_log);
    time(&endtime);
    WriteDecompResult(f, f_log);
    FreeDDMemory();
  } else {
    time(&endtime);
    WriteDecompResult(f, f_log);
    WriteErrorMessages(cout);
    WriteErrorMessages(f);
  }
}

void DDRowDecomposition(ostream &f, ostream &f_log)
{
  rowrange i,k;
  long FeasibleRaySum=0;
  time_t starttime_save;
  
  time(&starttime_save);
  SetWriteFileName(dexfile, 'd', "row-decomposition");
  ofstream writing_dex(dexfile);

  RestrictedEnumeration=True;
  RelaxedEnumeration=True;
  for (i = 0; i <= mm; i++) EqualityIndex[i]=0;
  for (k = 1; k <= mm-nn+2; k++){
    EqualityIndex[k]=1;   /* Equality for k-th inequality */
    if (k>=2) EqualityIndex[k-1]=-1;  /* Strict inequality for 1,2,...,(k-1)st inequalities */
    if (DynamicWriteOn) {
      writing_dex << "* Decomposition problem number = " << k << "(/" << (mm-nn+2) << ")\n";
      cout << "* Decomposition problem number = " << k << "(/" << (mm-nn+2) << ")\n";
    }
    DecompositionCore(writing_dex, f_log);
    FeasibleRaySum=FeasibleRaySum+FeasibleRayCount;
  }
  switch (Inequality) {
  case ZeroRHS:
    writing_dex << "*Total outputs = " << FeasibleRaySum << " " << (nn+1) 
       << "  " << OutputNumberString << "\n";
    cout << "*Total outputs = " << FeasibleRaySum << " " << (nn+1) 
       << "  " << OutputNumberString << "\n";
    break;
  case NonzeroRHS:
    writing_dex << "*Total outputs = " << FeasibleRaySum << " " << nn 
       << "  " << OutputNumberString << "\n";
    cout << "*Total outputs = " << FeasibleRaySum << " " << nn 
       << "  " << InputNumberString << "\n";
    break;
  }
  DDInit();
  for (i = 0; i <= mm; i++) {
    EqualityIndex[i]=0;
    set_addelem(AddedHyperplanes,i);
  }
  CompileDecompResult(writing_dex);
  starttime=starttime_save;
  RestrictedEnumeration=False;
  RelaxedEnumeration=False;
  WriteExtFile(f, f_log);
  if (IncidenceOutput == OutputIncidence || IncidenceOutput == IOIncidence){
    SetWriteFileName(incfile, 'i', "incidence");
    ofstream writing_ocd(incfile);
    WriteIncidenceFile(writing_ocd);
    if (DynamicWriteOn) printf("closing the file %s\n",incfile);
    writing_ocd.close();
  }
  if (IncidenceOutput == InputIncidence || IncidenceOutput == IOIncidence){
    SetWriteFileName(incfile, 'n', "input_incidence");
    ofstream writing_icd(incfile);
    WriteInputIncidenceFile(writing_icd);
    if (DynamicWriteOn) printf("closing the file %s\n",incfile);
    writing_icd.close();
  }
  if (AdjacencyOutput == OutputAdjacency || AdjacencyOutput == IOAdjacency){
    SetWriteFileName(adjfile, 'a', "adjacency");
    ofstream writing_adj(adjfile);
    if (DynamicWriteOn) printf("Writing the adjacency file %s...\n",adjfile);
    WriteAdjacencyFile(writing_adj);
    writing_adj.close();
    if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
  }
  writing_dex.close();
  FreeDDMemory();
}

void SolveRowSubproblem(ostream &f, ostream &f_log)
{
/* 
  This is not complete.  This part should be very similar to
  DDenumerate.  The only difference is that this first read the 
  rows specified by SubProblemRowSet, set mm, nn, etc properly
  and run DDenumerate.  For this, AmatrixInput must be modified
  so that it does not read the matrix data for this mode: 
    SpecialConversion=RowSubproblemSolve.
  Then it is necessary to make a new AmatrixInput program which
  reads a specified set of rows for the matrix AA.  
*/
}

void PreProjection(ostream &f, ostream &f_log)
{
  rowset subrows1,subrows2,DBrows;
  colset subcols1,subcols2;  /* subcols1:projvars,  subcols2:rest */
  rowrange i;
  colrange j,k;
  colindex pivrow;
  Bmatrix DBinv;  /* dual basis matrix inverse */
  long DBrank;
 
  time(&starttime);
  set_initialize(&subrows1,mm);
  set_initialize(&subrows2,mm);
  set_initialize(&DBrows,mm);
  set_initialize(&subcols1,nn);  /* subcol1 : projvar & RHS columns */
  set_initialize(&subcols2,nn);  /* subcol2 : remaining columns */
  SetWriteFileName(projfile, 'p', "preprojection variable subsystem");
  ofstream writing_proj(projfile);
  for (j=1;j<=nn;j++){
    if (set_member(j,projvars) || (j==1 && Inequality==NonzeroRHS))
      set_addelem(subcols1,j);
    else
      set_addelem(subcols2,j);
  }
  for (i=1; i<=mm; i++) set_addelem(subrows1,i);
  if (DynamicWriteOn){
    WriteSubMatrixOfAA(cout,subrows1,subcols1,Inequality);
  }
  WriteSubMatrixOfAA(writing_proj,subrows1,subcols1,Inequality);
  Inequality=ZeroRHS;
  ReduceAA(subrows1,subcols2);
    /* Extract the submatrix of AA index by subcols2. 
       subcols2 is changed to a consecutive sequence starting from 1 */
  if (debug) {
    WriteAmatrix(cout,AA,mm,nn,NonzeroRHS);
    WriteAmatrix(f,AA,mm,nn,NonzeroRHS);
  }
  PreOrderedRun=False;
  InitializeBmatrix(DBinv);
  FindBasis(AA,MinIndex,DBrows,pivrow,DBinv,&DBrank);
    /* DBrows stores the rows associated with a dual basis */
  if (debug){
    printf("rank of the new (deletion col) matrix is %ld\n", DBrank);
    printf("dual basis rows ="); set_write(DBrows);
    for (j=1;j<=nn;j++) (f) << "pivot row at col " << j << " = " << pivrow[j] << "\n";
  }
  set_diff(subrows2,subrows1,DBrows); 
    /* subrows2 stores the rows not in DBrows */
  for (j=1; j<=nn;j++){
    if (pivrow[j]==0) {
      set_delelem(subcols2,j);
      (f) << "Warning: col " << j << " is a linear combination of the other colums. The column linear dependency must be deleted for ray computation\n";
      for (k=j; k<=nn-1; k++){ /* shifting all pivrow information */
        pivrow[j]=pivrow[j+1];
      }
      pivrow[nn]=0;
      nn--;
    }
  }
  if (debug)  {
    printf("rows for ray enumeration:");set_write(subrows2);
    printf("cols for ray enumeration:");set_write(subcols2);
  }
  ReduceAA(subrows2,subcols2); 
    /* subrows2 is changed to a consecutive sequence starting from 1 */
  DualizeAA(DBinv);
  if (Error==DimensionTooLarge) goto _L99;
  if (debug) {
    WriteAmatrix(cout,AA,mm,nn,ZeroRHS);
    WriteAmatrix(f,AA,mm,nn,ZeroRHS);
  }
  if (DynamicWriteOn) {
    WriteRunningMode(cout);
  }
  DDInit();
  FindInitialRays(InitialHyperplanes, InitialRays, InitialRayIndex, &found);
  if (found) {
    InitialDataSetup();
    InitialWriting(f, f_log);
    DDMain(f,f_log);
    WriteProjResult(f, f_log, pivrow);
    if (IncidenceOutput == OutputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'i', "incidence");
      ofstream writing_ocd(incfile);
      WriteIncidenceFile(writing_ocd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_ocd.close();
    }
    if (IncidenceOutput == InputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'n', "input_incidence");
      ofstream writing_icd(incfile);
      WriteInputIncidenceFile(writing_icd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_icd.close();
    }
    if (AdjacencyOutput == OutputAdjacency || AdjacencyOutput == IOAdjacency){
      SetWriteFileName(adjfile, 'a', "adjacency");
      ofstream writing_adj(adjfile);
      if (DynamicWriteOn) printf("Writing the adjacency file %s...\n",adjfile);
      WriteAdjacencyFile(writing_adj);
      writing_adj.close();
      if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
    }
    FreeDDMemory();
  } else {
    _L99:;
    WriteErrorMessages(cout);
    WriteErrorMessages(f);
  }
  set_free(&subrows1);
  set_free(&subrows2);
  set_free(&DBrows);
  set_free(&subcols1);
  set_free(&subcols2);
}



void CompileDecompResult(ofstream &f)
{
  long i,j,k;
  myTYPE value=0;
  long mray,nray;
  char numbtype[wordlenmax],command[wordlenmax];
  boolean localdebug=False;
  myTYPE* vec;
  
  vec = new myTYPE[mm];
  AddArtificialRay();
  if ((f).is_open()){
    (f).close();
    if (DynamicWriteOn) printf("closing the file %s\n",dexfile);
  }
  ifstream reading_dex(dexfile);
  for (i=1; i<=mm-nn+2;i++){
    found=False;
    while (!found)
    {
      if (reading_dex.eof()) {
       Error=ImproperInputFormat;
       goto _L99;
      }
      else {
        reading_dex >> command;
        if (strncmp(command, "begin", 5)==0) {
          found=True;
        }
      }
    }
    reading_dex >> mray;
    reading_dex >> nray;
    reading_dex >> numbtype;
    if (localdebug) printf("decomp size = %ld x %ld\nNumber Type = %s\n", mray, nray, numbtype);
    for (k=1; k<=mray;k++){
      for (j=1; j<=nray; j++){
        reading_dex >> value;
        if (Inequality==NonzeroRHS) {
          vec[j - 1] = value;
        } else if (j>=2) {
          vec[j - 2] = value;
        }
        if (localdebug) WriteNumber(cout, value);
      }
      if (localdebug) printf("\n");
      AddRay(vec);
    }
  }
_L99:;
}

void  PostAnalysisMain(ifstream &f, ostream &f_log)
{
  rowrange i,k;
  
  time(&starttime);
  DDInit();
  if (f.is_open()) {
    ReadExtFile(f);
    RestrictedEnumeration=False;
    RelaxedEnumeration=False;
    for (i=1; i<=mm; i++) set_addelem(AddedHyperplanes, i);
    if (IncidenceOutput == OutputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'i', "incidence");
      ofstream writing_ocd(incfile);
      WriteIncidenceFile(writing_ocd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_ocd.close();
    }
    if (IncidenceOutput == InputIncidence || IncidenceOutput == IOIncidence){
      SetWriteFileName(incfile, 'n', "input_incidence");
      ofstream writing_icd(incfile);
      WriteInputIncidenceFile(writing_icd);
      if (DynamicWriteOn) printf("closing the file %s\n",incfile);
      writing_icd.close();
    }
    if (AdjacencyOutput == OutputAdjacency || AdjacencyOutput == IOAdjacency){
      SetWriteFileName(adjfile, 'a', "adjacency");
      ofstream writing_adj(adjfile);
      if (DynamicWriteOn) printf("Writing the adjacency file %s...\n",adjfile);
      WriteAdjacencyFile(writing_adj);
      writing_adj.close();
      if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
    }
    if (AdjacencyOutput == InputAdjacency || AdjacencyOutput == IOAdjacency){
      SetWriteFileName(adjfile, 'j', "input_adjacency");
      ofstream writing_iad(adjfile);
      if (DynamicWriteOn) printf("Writing the input_adjacency file %s...\n",adjfile);
      WriteInputAdjacencyFile(writing_iad);
      writing_iad.close();
      if (DynamicWriteOn) printf("closing the file %s\n",adjfile);
    }
    FreeDDMemory();
  } else {
    Error=FileNotFound;
    WriteErrorMessages(cout);
    WriteErrorMessages(f_log);
  }
}

void LPInit(void)
{
  rowrange i;

  Error=None;
  CompStatus=InProgress;
  if (debug) WriteAmatrix(cout,AA,mm,nn, Inequality);
  OrderVector = new long[mm+3];/* two more element for auxiliary variables */
  for (i=1; i<=mm+3; i++) OrderVector[i-1]=i-1; 
}


void DDInit(void)
{
  colrange j;

  Error=None;
  CompStatus=InProgress;
  SetInequalitySets(EqualityIndex);
  set_initialize(&InitialHyperplanes,mm);
  set_initialize(&AddedHyperplanes,mm);
  set_initialize(&WeaklyAddedHyperplanes,mm);
  set_initialize(&Face, mm);   /* used in CheckAdjacency  */
  set_initialize(&Face1, mm);  /* used in CheckAdjacency  */
  OrderVector=(long *)calloc(mm+1, sizeof *OrderVector);
  if (debug) WriteAmatrix(cout,AA,mm,nn, Inequality);
  ComputeRowOrderVector(OrderVector, HyperplaneOrder);
  RecomputeRowOrder=False;
  InitializeBmatrix(InitialRays);
  if (debug) WriteBmatrix(cout,InitialRays);
  RayCount = 0;
  TotalRayCount = 0;
  FeasibleRayCount = 0;
  WeaklyFeasibleRayCount = 0;
  VertexCount = 0;
  EdgeCount=0; /* active edge count */
  TotalEdgeCount=0; /* active edge count */
}

void InitialDataSetup(void)
{
  long j, r;
  rowset ZSet;
  myTYPE *vec1, *vec2;

  vec1 = new myTYPE[nn];
  vec2 = new myTYPE[nn];
  RecomputeRowOrder=False;
  ArtificialRay = NULL;
  FirstRay = NULL;
  LastRay = NULL;
  set_initialize(&ZSet,mm);
  AddArtificialRay();
  Iteration = nn;   /*Initially,we have already  nn  hyperplanes */
  set_copy(AddedHyperplanes, InitialHyperplanes);
  set_copy(WeaklyAddedHyperplanes, InitialHyperplanes);
  UpdateRowOrderVector(InitialHyperplanes);
  for (r = 1; r <= nn; r++) {
    for (j = 0; j < nn; j++){
      vec1[j] = InitialRays[j][r-1];
      vec2[j] = -InitialRays[j][r-1];
    }
    Normalize(vec1);
    Normalize(vec2);
    ZeroIndexSet(vec1, ZSet);
    if (set_subset(EqualitySet, ZSet)){
      if (debug) {
        printf("add an initial ray with zero set:");
        set_write(ZSet);
      }
      AddRay(vec1);
      if (InitialRayIndex[r]==0) {
        AddRay(vec2);
        if (debug) {
          printf("and add its negative also.\n");
        }
      }
    }
  }
  CreateInitialEdges();
  set_free(&ZSet);
}


void LPMain(ostream &f, ostream &f_log)
{
  static colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  static Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  static Bmatrix BasisInv;

  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  myTYPE ov=0;  /* LP optimum value */
  long LPiter;

  if (Inequality==ZeroRHS){
    // printf("Sorry, LP optimization is not implemented for RHS==0.\n");
    // goto _L99;
    EnlargeAAforZeroRHSLP();
  }
  time(&starttime);
  LPsol = new myTYPE[nn];
  LPdsol = new myTYPE[nn];
  InitializeBmatrix(BasisInv);
  boolean UsePrevBasis=False;
  if (Conversion==LPmax){
    if (LPsolver==DualSimplex){
      DualSimplexMaximize(f, f_log, AA, BasisInv, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
    } else {
      CrissCrossMaximize(f, f_log, AA, BasisInv, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
    }
  }
  else if (Conversion==LPmin){
    if (LPsolver==DualSimplex){
      DualSimplexMinimize(f, f_log, AA, BasisInv, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
    } else {
      CrissCrossMinimize(f, f_log, AA, BasisInv, OBJrow, RHScol, UsePrevBasis,
        &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
    }
  }
  WriteLPResult(f, LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
  if (DynamicWriteOn)
    WriteLPResult(cout,LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
_L99:;
}

void InteriorFindMain(ostream &f, ostream &f_log, boolean *found)
{
  static colindex NBIndex;  /* NBIndex[s] stores the nonbasic variable in column s */ 
  static Arow LPsol, LPdsol;  /*  LP solution and the dual solution (basic var only) */
  static Bmatrix BasisInv;

  rowrange re;  /* evidence row when LP is inconsistent */
  colrange se;  /* evidence col when LP is dual-inconsistent */
  myTYPE ov=0;  /* LP optimum value */
  long LPiter;
  boolean UsePrevBasis=False;

  *found = False;
  if (Inequality==ZeroRHS){
    printf("Sorry, find_interior is not implemented for RHS==0.\n");
    goto _L99;
  }
  EnlargeAAforInteriorFinding();
  InitializeBmatrix(BasisInv);
  LPsol = new myTYPE[nn];
  LPdsol = new myTYPE[nn]; 
  time(&starttime);
  OBJrow=mm; RHScol=1;
  DualSimplexMaximize(f, f_log, AA, BasisInv, OBJrow, RHScol, UsePrevBasis,
    &LPStatus, &ov, LPsol, LPdsol,NBIndex, &re, &se, &LPiter);
  WriteLPResult(f, LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
  if (LPStatus==Optimal || LPStatus==DualInconsistent) *found=True;
  if (DynamicWriteOn)
    WriteLPResult(cout,LPStatus, ov, LPsol, LPdsol, NBIndex, re, se, LPiter);
  RecoverAAafterInteriorFinding();
_L99:;
}

void CheckConversionConsistency(ostream &f, ostream &f_log)
{
  switch (Conversion) {
  case ExtToIne: /* facet enumeration is chosen */
    if (strcmp(ifiletail,"ext")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ext'' for hull computation.\n";
      f << "*Warning: Extension of inputfile name should ``.ext'' for hull computation.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ext'' for hull computation.\n";
    }
    break;

  case IneToExt: /* vertex enumeration is chosen */
    if (strcmp(ifiletail,"ine")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ine'' for vertex/ray enumeration.\n";
      f << "*Warning: Extension of inputfile name should be ``.ine'' for vertex/ray enumeration.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ine'' for vertex/ray enumeration.\n";
    }
    break;
    
  case LPmax:  case LPmin:  case InteriorFind:    /* LP is chosen */
    if (strcmp(ifiletail,"ine")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ine'' for LP optimization/interior_find.\n";
      f << "*Warning: Extension of inputfile name should be ``.ine'' for LP optimization/interior_find.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ine'' for LP optimization/interior_find.\n";
    }
    break;

  case FacetListing: case FacetListingExternal: 
           /* Facet Listing is chosen */
    if (strcmp(ifiletail,"ine")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ine'' for facet_listing.\n";
      f << "*Warning: Extension of inputfile name should be ``.ine'' for facet_listing.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ine'' for facet_listing.\n";
    }
    break;

   case VertexListing:  case VertexListingExternal:
        /* Vertex Listing with ExternalFile is chosen */
    if (strcmp(ifiletail,"ext")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ext'' for vertex_listing.\n";
      f << "*Warning: Extension of inputfile name should be ``.ext'' for vertex_listing.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ext'' for vertex_listing.\n";
    }
    break;
  
  case TopeListing:           /* TopeListing is chosen */
    if (strcmp(ifiletail,"ine")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ine'' for tope_listing.\n";
      f << "*Warning: Extension of inputfile name should be ``.ine'' for tope_listing.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ine'' for tope_listing.\n";
    }
    break;

  case Projection:            /* preprojection is chosen */
    if (strcmp(ifiletail,"ine")!=0){
      cout << "*Warning: Extension of inputfile name should be ``.ine'' for preprojection.\n";
      f << "*Warning: Extension of inputfile name should be ``.ine'' for preprojection.\n";
      f_log << "*Warning: Extension of inputfile name should be ``.ine'' for preprojection.\n";
    }
    break;

  default: break;
  }
}

int main(int argc, char *argv[])
{
  OutputHeading();
  DefaultOptionSetup();
  Initialization(argc, argv);
  AmatrixInput(&inputsuccessful);

  if (inputsuccessful) {
    SetWriteFileName(logfile,'l',"log");
    ofstream writing_log(logfile);
    if (VerifyInput){
      SetWriteFileName(verfile,'v',"input verification");
      ofstream writing_ver(verfile);
      WriteSolvedProblem(writing_ver);
      writing_ver.close();
      if (DynamicWriteOn) printf("closing the file %s\n",verfile);
     }
    if (DynamicWriteOn) {
      WriteRunningMode(cout);
    }
    SetWriteFileName(outputfile,'o',"output");
    if (PostAnalysisOn){
      /* Post analysis is chosen */
      ifstream reading_ext(outputfile);
      PostAnalysisMain(reading_ext, writing_log);
    }
    else {
      ofstream writing(outputfile);
      CheckConversionConsistency(writing,writing_log);
      switch (Conversion) {
      case ExtToIne: case IneToExt: /* vertex/facets enumeration is chosen */
        switch (SpecialConversion){
        case RowDecomposition:
          DDRowDecomposition(writing,writing_log);
          break;
        case RowSubproblemSolve:
          SolveRowSubproblem(writing,writing_log);
          break; 
        case Nothing:
          DDEnumerate(writing, writing_log);
          break;
        }
        break;
    
      case LPmax:  case LPmin:      /* LP is chosen */
        LPInit();
        LPMain(writing, writing_log);
        break;

      case FacetListing: case VertexListing:
                          /* Facet or Vertex Listing is chosen */
        LPInit();
        FacetandVertexListMain(writing, writing_log);
        break;

      case FacetListingExternal: case VertexListingExternal:
        /* Facet or Vertex Listing with ExternalFile is chosen */
        LPInit();
        FacetandVertexExternalListMain(writing, writing_log);
        break;
  
      case TopeListing:           /* TopeListing is chosen */
        LPInit();
        TopeListMain(writing, writing_log);
        break;

      case Projection:            /* preprojection is chosen */
        PreProjection(writing, writing_log);
        break;

      case InteriorFind:      /* Interior point search is chosen */
        boolean found;
        LPInit();
        InteriorFindMain(writing, writing_log, &found);
        break;

      default: break;
      }
      if (writing.is_open()) {
        writing.close();
        if (DynamicWriteOn) printf("closing the file %s\n",outputfile);
      }   
    }
    if (writing_log.is_open()) {
      writing_log.close();
      if (DynamicWriteOn) printf("closing the file %s\n",logfile);
    }
  } else {
    ofstream writing("cdd.error");
    WriteErrorMessages(writing);
    writing.close();
  }

}


/* end of cdd.C */
