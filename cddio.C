/* cddio.C:  Basic Input and Output Procedures for cdd.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.73, September 6, 1995 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <fstream.h>
#include <libc.h>
#include <strclass.h>
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

void SetInputFile(boolean *success)
{
  boolean opened=False,stop, quit=False;
  int i,dotpos=0, semipos=0;
  char ch;
  char *tempname;
  
  *success=False;
  while (!opened && !quit) {
    if (FileInputMode!=Auto){
      printf("\n>> Input file (*.ine) : ");
      scanf("%s",inputfile);
      ch=getchar();
    }
    stop=False;
    for (i=0; i<filenamelen && !stop; i++){
      ch=inputfile[i];
      switch (ch) {
        case '.': 
          dotpos=i+1;
          break;
        case ';':  case ' ':  case '\0':  case '\n':  case '\t':     
          if (ch==';'){
            semipos=i+1;
            FileInputMode=SemiAuto;   
            /* semicolon at the end of the filename
            -> output file names will be creeated with defaults. */
          }
          stop=True;
          tempname=(char*)calloc(filenamelen,sizeof ch);
          strncpy(tempname, inputfile, i);
          strcpy(inputfile,tempname);
          break;
      }
    }
    if (dotpos>0){
      strncpy(ifilehead, inputfile, dotpos-1);
    }else{
      strcpy(ifilehead, inputfile);
    }
    if (debug){
      printf("inputfile name: %s\n", inputfile);  
      printf("inputfile name head: %s\n", ifilehead);  
      printf("semicolon pos: %d\n", semipos);
    }  
    ifstream ftemp(inputfile);
    if (ftemp) {
      if (DynamicWriteOn) printf("input file %s is open\n", inputfile);
      opened=True;
      *success=True;
      ftemp.close();
    }
    else{
      printf("The file %s not found\n",inputfile);
      if (FileInputMode==Auto) {
        quit=True;
      }
    }
  }
}

void SetWriteFileName(DataFileType fname, char cflag, char *fscript)
{
  boolean quit=False;
  char *extension;
  
  switch (cflag) {
    case 'o':
      extension=".ext";break;   /* vertex and ray output file; general output file */
    case 'a':
      extension=".adj";break;   /* adjacency file */
    case 'i':
      extension=".icd";break;   /* incidence file */
    case 'j':
      extension=".iad";break;   /* input adjacency file */
    case 'l':
      extension=".ddl";break;   /* log file */
    case 'd':
      extension=".dex";break;   /* decomposition output */
    case 'p':
      extension="sub.ine";break;  /* preprojection sub inequality file */
    case 'v':
      extension=".solved";break;  /* verify_input file */
    default:
      extension=".xxx";break;
  }
  if (FileInputMode==Manual){
    while (!quit) {
      printf("\n>> %s file name (*%s)   : ",fscript, extension);
      scanf("%s",fname);
      if (fname[0]==';'|| fname[0]<'0'){
        quit=True;  /* default file name */
      } 
      else if (strcmp(inputfile, fname)!=0){
        goto _L99;
      }
      else {
        printf("%s file %s must have a name different from inputfile.\n",fscript,fname);
      }
    }
  }
  /* Auto or SemiAuto FileInput */
  strcpy(fname,ifilehead); 
  strcat(fname,extension); 
  if (strcmp(inputfile, fname)==0) {
    strcpy(fname,inputfile); 
    strcat(fname,extension); 
  }
_L99:;
  if (DynamicWriteOn) printf("Open %s file %s.\n",fscript,fname);
}

void SetNumberType(string line)
{
  if (line== "integer") {
    Number = Integer;
    InputNumberString="integer";
#ifdef RATIONAL
    OutputNumberString="rational";
#else
    OutputNumberString="real";
#endif
    return;
  }
  else if (line== "rational") {
    Number = Rational;
    InputNumberString="rational";
    OutputNumberString="rational";
#ifndef RATIONAL
    Error=ImproperExecutable; /* the executable cannot handle rational data */ 
#endif
    return;
  }
  else if (line== "real") {
    Number = Real;
    InputNumberString="real";
    OutputNumberString="real";
#ifdef RATIONAL
    Error=ImproperExecutable; /* the executable cannot handle real data */ 
#endif
    return;
  }
  else { 
    Number=Unknown;
    Error=ImproperInputFormat;
  }
}

void ProcessCommandLine(ifstream &f, string line)
{
  colrange j;
  long var,msize;
  myTYPE cost=0;

  if (debug) cout << line << "\n";
  if (line=="dynout_off") {
    DynamicRayWriteOn = False;
    return;
  }
  if (line== "stdout_off") {
    DynamicRayWriteOn = False;
    DynamicWriteOn = False;
    return;
  }
  if (line== "logfile_on") {
    LogWriteOn = True;
    return;
  }
  if (line== "logfile_off") {
    LogWriteOn = False;
    return;
  }
  if (line== "hull") {
    Conversion = ExtToIne;
    return;
  }
  if (line== "incidence") {
    IncidenceOutput = IncSet;
    return;
  }
  if (line== "#incidence") {
    IncidenceOutput = IncCardinality;
    return;
  }
  if (line== "adjacency") {
    if (AdjacencyOutput==InputAdjacency || AdjacencyOutput==IOAdjacency)
      AdjacencyOutput=IOAdjacency;
    else
      AdjacencyOutput = OutputAdjacency;
    return;
  }
  if (line== "input_adjacency") {
    if (AdjacencyOutput==OutputAdjacency || AdjacencyOutput==IOAdjacency)
      AdjacencyOutput = IOAdjacency;
    else 
      AdjacencyOutput = InputAdjacency;
    return;
  }
  if (line== "postanalysis") {
    PostAnalysisOn=True;
    return;
  }
/*  algebraic option is not efficient in most cases and deleted from Version 051 */
/*
  if (line== "algebraic") {
    AdjacencyTest = Algebraic;
    return;
  }
*/
  if (line== "nondegenerate") {
    NondegAssumed = True;
    return;
  }
  if (line== "minindex") {
    HyperplaneOrder = MinIndex;
    return;
  }
  if (line== "maxindex") {
    HyperplaneOrder = MaxIndex;
    return;
  }
  if (line== "mincutoff") {
    HyperplaneOrder = MinCutoff;
    return;
  }
  if (line== "maxcutoff") {
    HyperplaneOrder = MaxCutoff;
    return;
  }
  if (line== "mixcutoff") {
    HyperplaneOrder = MixCutoff;
    return;
  }
  if (line== "lexmin") {
    HyperplaneOrder = LexMin;
    return;
  }
  if (line== "lexmax") {
    HyperplaneOrder = LexMax;
    return;
  }
  if (line== "random") {
    HyperplaneOrder = RandomRow;
    (f) >> rseed;
    if (rseed <= 0) rseed = 1;
    return;
  }
  if (line== "lineshell") {
    HyperplaneOrder = LineShelling;
    return;
  }
  if (line== "initbasis_at_bottom") {
    InitBasisAtBottom = True;
    return;
  }
  if (line== "verify_input") {
    VerifyInput = True;
    return;
  }
  if (line== "debug") {
    debug = True;
    return;
  }
  if ((line== "partial_enum" || line=="equality") 
    && RestrictedEnumeration==False) {
    (f) >> msize;
    for (j=1;j<=msize;j++) {
      (f) >> var;
      EqualityIndex[var]=1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Partial Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    RestrictedEnumeration=True;
    return;
  }
  if (line== "strict_ineq" && RelaxedEnumeration==False) {
    (f) >> msize;
    for (j=1;j<=msize;j++) {
      (f) >> var;
      EqualityIndex[var]=-1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Partial Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    RelaxedEnumeration=True;
    return;
  }
  if (line== "preprojection" && Conversion != Projection) {
    set_initialize(&projvars,nn);
    (f) >> projdim;
    if (debug) printf("dimension of projection = %ld  in variables:\n",projdim);
    for (j=1;j<=projdim;j++) {
      (f) >> var;
      if (debug) printf(" %ld",var);
      if (Inequality==NonzeroRHS)
        set_addelem(projvars,var+1);
      else
        set_addelem(projvars,var);
    }
    Conversion=Projection;
    return;
  }
  if (line== "maximize" && Conversion != LPmax) {
    if (debug) printf("linear maximization is chosen.\n");
    LPcost=new myTYPE[ninput];
    for (j=0;j<ninput;j++) {
      (f) >> cost;
      LPcost[j]=cost;
      if (debug) cout << " cost[" << j << "] = "<< LPcost[j] <<"\n";
    }
    Conversion=LPmax;
    if (debug) {
      printf("\n");
    }
    return;
  }
  if (line== "minimize" && Conversion != LPmin) {
    if (debug) printf("linear minimization is chosen.\n");
    LPcost=new myTYPE[ninput];
    for (j=0;j<ninput;j++) {
      (f) >> cost;
      LPcost[j]=cost;
      if (debug) cout << " cost[" << j << "] = "<< LPcost[j] <<"\n";
    }
    Conversion=LPmin;
    if (debug) {
      printf("\n");
    }
    return;
  }
  if (line== "find_interior" && Conversion != InteriorFind) {
    if (debug) printf("Interior finding option is chosen.\n");
    Conversion=InteriorFind;
    return;
  }
  if (line== "facet_listing" && Conversion != LPmax && Conversion != LPmin) {
    if (debug) printf("facet_listing is chosen.\n");
    Conversion=FacetListing;
    return;
  }
  if (line== "tope_listing" && Conversion != LPmax && Conversion != LPmin) {
    if (debug) printf("tope_listing is chosen.\n");
    Conversion=TopeListing;
    return;
  }
  if (line== "condensed_listing") {
    CondensedListOn=True;
    return;
  }
  if (line== "signpivot") {
    SignPivotOn=True;
    return;
  }
  if (line== "CMIalgorithm") {
    if (DynamicWriteOn) printf("Use the combinatorial maximum improvement method.\n");
    LPsolver=CombMaxImprove; SignPivotOn=True;
    return;
  }
  if (line== "dual_simplex") {
    if (DynamicWriteOn) printf("Use the dual simplex method.\n");
    LPsolver=DualSimplex;
    return;
  }
  if (line== "criss-cross") {
    if (DynamicWriteOn) printf("Use the crss-cross method.\n");
    LPsolver=CrissCross;
    return;
  }
  if (line== "manual_pivot") {
    ManualPivotOn=True;
    return;
  }
  if (line== "show_tableau" && !ShowSignTableauOn) {
    if (DynamicWriteOn) printf("Sign tableau will be output to CRT.\n");
    ShowSignTableauOn=True;
    return;
  }
  if (line== "row_decomp" && !RowDecomposition) {
    if (DynamicWriteOn) printf("Row decomposition is chosen.\n");
    RowDecomposition=True;
    return;
  }
}


void WriteRayRecord(ostream &f, RayRecord *RR)
{
  long j;
  myTYPE scaler=0;

  if (Inequality == ZeroRHS) {
    (f) << " " << 0;
    for (j = 0; j < nn; j++)
      WriteNumber(f, RR->Ray[j]);
  } 
  else {
    scaler = FABS(RR->Ray[0]);
    if (scaler > zero) {
      if (RR->Ray[0] > 0) {
        if (Conversion == IneToExt) {
	  (f) << " " << 1;
	  for (j = 1; j < nn; j++)
	    WriteNumber(f, RR->Ray[j] / scaler);
        } 
        else {
	      /* hull computation is chosen */
          for (j = 0; j < nn; j++)
	        WriteNumber(f, RR->Ray[j]);
	    }
      }
      else {
        /* hull computation must have been chosen, since RR->Ray[0] is negative */
	    for (j = 0; j < nn; j++)
	      WriteNumber(f, RR->Ray[j]);
      }
    } 
    else {
      (f) << " " << 0;
      for (j = 1; j < nn; j++)
        WriteNumber(f, RR->Ray[j]);
    }
  }
  if (IncidenceOutput==IncCardinality) {
    (f) << " : " << set_card(RR->ZeroSet);
  }
  (f).put('\n');
}


void WriteRayRecord2(ostream &f, RayRecord *RR)
{
  long j;

  (f) << " Ray = ";
  for (j = 0; j < nn; j++)
    (f) << RR->Ray[j];
  putchar('\n');
  (f) << " ZeroSet =";
  set_fwrite(f, RR->ZeroSet);
  (f).put('\n');
}

void WriteSubMatrixOfAA(ostream &f, rowset ChosenRow, colset ChosenCol,
      InequalityType ineq)
{
  long i,j;

  (f) << "begin\n";
  if (ineq==ZeroRHS)
    (f) << "  " <<  set_card(ChosenRow)  << "  " << (set_card(ChosenCol)+1) << "  " << InputNumberString << "\n";
  else
    (f) << "  " <<  set_card(ChosenRow) << "  " << set_card(ChosenCol) << "  " << InputNumberString << "\n";
  for (i=1; i <= mm; i++) {
    if (set_member(i,ChosenRow)) {
      if (ineq==ZeroRHS){  /* If RHS==0, output zero first */
        WriteNumber(f, 0);
        for (j=1; j <= nn; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteNumber(f, AA[i-1][j-1]);
          }
        }
      }
      else {
        for (j=1; j <= nn; j++) {
          if (set_member(j, ChosenCol)){ 
            WriteNumber(f, AA[i-1][j-1]);
          }
        }
      }
      (f) << "\n";
    }
  }
  (f) << "end\n"; 
}

void WriteAmatrix(ostream &f, Amatrix A, long rowmax, long colmax,
      InequalityType ineq)
{
  long i,j;

  (f) << "begin\n";
  if (ineq==ZeroRHS)
    (f) << "  " << rowmax << "  " << colmax+1 <<  "  " << InputNumberString << "\n";
  else
    (f) << "  " << rowmax << "  " << colmax  << "  " << InputNumberString << "\n";
  for (i=1; i <= rowmax; i++) {
    if (ineq==ZeroRHS)
      WriteNumber(f, 0);  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j <= colmax; j++) {
      WriteNumber(f, A[i-1][j-1]);
    }
    (f) << "\n";
  }
  (f) << "end\n";
}


void WriteNumber(ostream &f, myTYPE x)
{
#ifndef RATIONAL
  long ix1,ix2,ix;

  ix1= (long) (FABS(x) * 10000. + 0.5);
  ix2= (long) (FABS(x) + 0.5);
  ix2= (long) (ix2*10000);
  if ( ix1 == ix2) {
    if (x>0) {
      ix = (long) (x + 0.5);
    } else {
      ix = (long) (-x + 0.5);
      ix = -ix;
    }
    (f) << " " << ix;
  } else {
    (f).setf(ios::scientific,ios::floatfield);
    (f) << " " << x;
  }
#else
  (f) << " " << x;
#endif
}

void WriteTableau(ostream &f, Amatrix X, Bmatrix T, InequalityType ineq)
/* Write the tableau  X.T   */
{
  colrange j;
  rowrange i;
  
  (f) << "begin\n";
  if (ineq==ZeroRHS)
    (f) << "  " << mm << "  " << (nn+1) << "  real\n";
  else
    (f) << "  " << mm << "  " << nn << "  real\n";
  for (i=1; i<= mm; i++) {
    if (ineq==ZeroRHS)
      (f) << " " <<  0;  /* if RHS==0, the column is not explicitely stored */
    for (j=1; j<= nn; j++) {
      (f) << " " << TableauEntry(X,T,i,j);
    }
    (f) << "\n";
  }
  (f) << "end\n";
}

char Sign(myTYPE val)
{
  if (val > zero) return '+';
  else if ( val < -zero) return '-';
    else return '0'; 
}

char intSign(int val)
{
  if (val > 0) return '+';
  else if ( val < 0) return '-';
    else return '0'; 
}

void WriteSignAmatrix(ostream &f, SignAmatrix X,
 rowindex OV, long bflag[], rowrange objrow, colrange rhscol)
/* Write the sign matrix X. 
   This works properly only for Nonhomogeneous inequality  */
{
  colrange j;
  rowrange i,k,l;
  
  (f) << "      g|";
  for (l=1; l<= mm; l++) {
    j=OV[l];
    if  (bflag[j] >0) { /* j is nonbasic variable */
      (f).width(3);
      (f) << j;
    }
  }
  (f) << "\n";
  for (k=1; k<= mm; k++) {
    i=OV[k];
    if (bflag[i]==0) {  /* i the objective variable */
      (f) << "   f";
      (f) << "  " << intSign(X[i-1][0]) << "|";   
      for (l=1; l<= mm; l++) {
        j=OV[l];
        if  (bflag[j] >0) { /* j is nonbasic variable */
          (f) << "  " << intSign(X[i-1][bflag[j]-1]);
        }
      }
      (f) << "\n";
    }
  }
  (f) << "  =====+";
  for (j=1; j<=nn; j++) { 
    (f) << "===";
  }
  (f) << "\n";
  for (k=1; k<= mm; k++) {
    i=OV[k];
    if (bflag[i]!=0 && bflag[i]==-1) {  /* i is a basic variable */
      (f).width(4); 
      (f) << i;   
      (f) << "  " << intSign(X[i-1][0]) << "|";   
      for (l=1; l<= mm; l++) {
        j=OV[l];
        if  (bflag[j] >0) { /* j is nonbasic variable */
          (f) << "  " << intSign(X[i-1][bflag[j]-1]);
        }
      }
      (f) << "\n";
    }
  }
  (f) << "\n";
}

void WriteSignTableau(ostream &f, Amatrix X, Bmatrix T,
 rowindex OV, long bflag[], rowrange objrow, colrange rhscol)
/* Write the sign matrix of tableau  X.T.
   This works properly only for Nonhomogeneous inequality  */
{
  colrange j;
  rowrange i,k,l;
  
  (f) << "      g|";
  for (l=1; l<= mm; l++) {
    j=OV[l];
    if  (bflag[j] >0) { /* j is nonbasic variable */
      (f).width(3);
      (f) << j;
    }
  }
  (f) << "\n";
  for (k=1; k<= mm; k++) {
    i=OV[k];
    if (bflag[i]==0) {  /* i the objective variable */
      (f) << "   f";
      (f) << "  " << Sign(TableauEntry(X,T,i,1)) << "|";   
      for (l=1; l<= mm; l++) {
        j=OV[l];
        if  (bflag[j] >0) { /* j is nonbasic variable */
          (f) << "  " << Sign(TableauEntry(X,T,i,bflag[j]));
        }
      }
      (f) << "\n";
    }
  }
  (f) << "  =====+";
  for (j=1; j<=nn; j++) { 
    (f) << "===";
  }
  (f) << "\n";
  for (k=1; k<= mm; k++) {
    i=OV[k];
    if (bflag[i]!=0 && bflag[i]==-1) {  /* i is a basic variable */
      (f).width(4); 
      (f) << i;   
      (f) << "  " << Sign(TableauEntry(X,T,i,1)) << "|";   
      for (l=1; l<= mm; l++) {
        j=OV[l];
        if  (bflag[j] >0) { /* j is nonbasic variable */
          (f) << "  " << Sign(TableauEntry(X,T,i,bflag[j]));
        }
      }
      (f) << "\n";
    }
  }
  (f) << "\n";
}

void WriteBmatrix(ostream &f, Bmatrix T)
{
  colrange j1, j2;

  for (j1 = 0; j1 < nn; j1++) {
    for (j2 = 0; j2 < nn; j2++) {
      f << " " << T[j1][j2];
    }  /*of j2*/
    f.put('\n');
  }  /*of j1*/
  f.put('\n');
}

void WriteIncidence(ostream &f, RayRecord *RR)
{
  rowset cset;
  long zcar;

  set_initialize(&cset,mm);
  zcar = set_card(RR->ZeroSet);
  switch (IncidenceOutput) {

  case IncCardinality:
    (f) << " : " << zcar;
    break;

  case IncSet:
    if (mm - zcar >= zcar) {
      (f) << " " << zcar << " : ";
      set_fwrite(f, RR->ZeroSet);
    } else {
      set_diff(cset, GroundSet, RR->ZeroSet);
      (f) << " " << (zcar - mm) << " : ";
      set_fwrite(f, cset);
    }
    break;

  case IncOff:
    break;
  }
  (f).put('\n');
  set_free(&cset);
}

void WriteInputIncidence(ostream &f, Aincidence Aicd, rowrange i)
{
  set_type cset;
  long zcar;

  set_initialize(&cset,RayCount);
  zcar = set_card(Aicd[i-1]);
  if (RayCount - zcar >= zcar) {
    (f) << " " << zcar << " :";
    set_fwrite(f, Aicd[i-1]);
  } else {
    set_compl(cset, Aicd[i-1]);
    (f) << " " << (zcar - RayCount) << " :";
    set_fwrite(f, cset);
  }
  (f).put('\n');
  set_free(&cset);
}

void OutputHeading(void)
{
  cout << "* cdd+: Double Description Method:" << DDVERSION << "\n";
  cout << "* "<< COPYRIGHT << "\n";
  cout << "* "<< ARITHMETIC << "\n";
  cout << "----------------------------------------------\n";
  cout << " Enumeration of all vertices and extreme rays\n";
  cout << " of a convex polyhedron P={ x : b - A x >= 0}\n";
  cout << " Use hull option for convex hull computation!\n";
  cout << "----------------------------------------------\n";
}

void WriteProgramDescription(ostream &f)
{
  (f) << "* cdd+: Double Description Method in C++:" << DDVERSION << "\n";
  (f) << "* "<< COPYRIGHT << "\n";
  (f) << "* "<< ARITHMETIC << "\n";
}

void WriteRunningMode(ostream &f)
{
  if (Conversion==IneToExt || Conversion==ExtToIne){ 
    switch (HyperplaneOrder) {

    case MinIndex:
      (f) << "*HyperplaneOrder: MinIndex\n";
      break;

    case MaxIndex:
      (f) << "*HyperplaneOrder: MaxIndex\n";
      break;

    case MinCutoff:
      (f) << "*HyperplaneOrder: MinCutoff\n";
      break;

    case MaxCutoff:
      (f) << "*HyperplaneOrder: MaxCutoff\n";
    break;

    case MixCutoff:
      (f) << "*HyperplaneOrder: MixCutoff\n";
      break;

    case LexMin:
      (f) << "*HyperplaneOrder: LexMin\n";
      break;

    case LexMax:
      (f) << "*HyperplaneOrder: LexMax\n";
      break;

    case RandomRow:
      (f) << "*HyperplaneOrder: Random,  Seed = " << rseed << "\n";
      break;

    case LineShelling:
      (f) << "*HyperplaneOrder: LineShelling\n";
      break;
    }
    if (NondegAssumed) {
      (f) << "*Degeneracy preknowledge for computation: NondegenerateAssumed\n";
    }
    else {
      (f) << "*Degeneracy preknowledge for computation: None (possible degeneracy)\n";
    }
  }
  switch (Conversion) {
    case ExtToIne:
      (f) << "*Hull computation is chosen.\n";
      break;
    
    case IneToExt:
      (f) << "*Vertex/Ray enumeration is chosen.\n";
      break;
    
    case LPmax:
      (f) << "*LP (maximization) is chosen.\n";
      break;

    case LPmin:
      (f) << "*LP (minimization) is chosen.\n";
      break;

    case FacetListing:
      (f) << "*Facet listing is chosen.\n";
      break;

    case TopeListing:
      (f) << "*Tope listing is chosen.\n";
      break;

    case Projection:
      (f) << "*Preprojection is chosen.\n";
      break;
  
    default: break;
  }
  if (RestrictedEnumeration) {
    (f) << "*The equality option is chosen.\n* => Permanently active rows are:";
    set_fwrite(f,EqualitySet);
    (f) << "\n";
  }
  if (RelaxedEnumeration) {
    (f) << "*The strict_inequality option is chosen.\n* => Permanently nonactive rows are:";
    set_fwrite(f,NonequalitySet);
    (f) << "\n";
  }
  if (PostAnalysisOn) {
    (f) << "*Post analysis is chosen.\n";
  }
}

void WriteRunningMode2(ostream &f)
{
  long j;

  if (Conversion==IneToExt || Conversion==ExtToIne){ 
    switch (HyperplaneOrder) {

    case MinIndex:
      (f) << "minindex\n";
      break;

    case MaxIndex:
      (f) <<  "maxindex\n";
      break;

    case MinCutoff:
      (f) <<  "mincutoff\n";
      break;

    case MaxCutoff:
      (f) <<  "maxcutoff\n";
    break;

    case MixCutoff:
      (f) <<  "mixcutoff\n";
      break;

    case LexMin:
      (f) <<  "lexmin\n";
      break;

    case LexMax:
      (f) <<  "lexmax\n";
      break;

    case RandomRow:
      (f) <<  "random  " << rseed << "\n";
      break;

    case LineShelling:
      (f) <<  "lineshelling\n";
      break;
    }
    if (NondegAssumed) {
      (f) <<  "nondegenerate\n";
    }
  }
  switch (Conversion) {
    case ExtToIne:
      (f) <<  "hull\n";
      break;
    
    case LPmax:
      (f) <<  "maximize\n";
      for (j=0; j<ninput; j++) WriteNumber(f,LPcost[j]);
      (f) << "\n";
      break;

    case LPmin:
      (f) <<  "minimize\n";
      for (j=0; j<ninput; j++) WriteNumber(f,LPcost[j]);
      (f) << "\n";
      break;

    case FacetListing:
      (f) << "facet_listing\n";
      break;

    case TopeListing:
      (f) << "tope_listing\n";
      break;

    case Projection:
      (f) <<  "*Preprojection is chosen.\n";
      break;
  
    default: break;
  }
  if (RestrictedEnumeration) {
    (f) <<  "equality ";
    set_fwrite(f,EqualitySet);
    (f) << "\n";
  }
  if (RelaxedEnumeration) {
    (f) <<  "strict_inequality ";
    set_fwrite(f,NonequalitySet);
    (f) << "\n";
  }
  if (PostAnalysisOn) {
    (f) << "postanalysis\n";
  }
}

void WriteCompletionStatus(ostream &f)
{
  if (Iteration<mm && CompStatus==AllFound) {
    (f) << "*Computation completed at Iteration " << Iteration << ".\n";
  } 
  if (CompStatus == RegionEmpty) {
    (f) << "*Computation completed at Iteration " << Iteration << " because the region found empty.\n";
  }
}


void WriteTimes(ostream &f)
{ 
  long ptime,ptime_sec,ptime_minu, ptime_hour;
  
  /* ptime=difftime(endtime,starttime); */   /* This function is ANSI standard, but not available sometime */
  ptime=endtime-starttime;      /* This is to replace the line above, but it may not give correct time in seconds */ 
  ptime_hour=ptime/3600;
  ptime_minu=(ptime-ptime_hour*3600)/60;
  ptime_sec=ptime%60;
  (f) << "*Computation starts     at " << asctime(localtime(&starttime));
  (f) << "*            terminates at " << asctime(localtime(&endtime));
  (f) << "*Total processor time = " << ptime << " seconds\n";
  (f) << "*                     = " <<  ptime_hour << "h " << ptime_minu << "m " << ptime_sec << "s\n";
}

boolean InputAdjacentQ(Aincidence Aicd, 
  set_type RedundantSet, set_type DominantSet, 
  rowrange i1, rowrange i2)
/* Before calling this function, RedundantSet must be 
   a set of row indices whose removal results in a minimal
   nonredundant system to represent the input polyhedron,
   DominantSet must be the set of row indices which are
   active at every extreme points/rays.
*/
{
  boolean adj=True;
  rowrange i;
  static set_type common;
  static lastRayCount=0;

  if (lastRayCount!=RayCount){
    if (lastRayCount >0) set_free(&common);
    set_initialize(&common, RayCount);
    lastRayCount=RayCount;
  }
  if (set_member(i1, RedundantSet) || set_member(i2, RedundantSet)){
    adj=False;
    goto _L99;
  }
  if (set_member(i1, DominantSet) || set_member(i2, DominantSet)){
  // dominant inequality is considered adjacencent to all others.
    adj=True;
    goto _L99;
  }
  set_int(common, Aicd[i1-1], Aicd[i2-1]);
  i=0;
  while (i<mm && adj==True){ 
    i++; 
    if (i!=i1 && i!=i2 && !set_member(i, RedundantSet) &&
        !set_member(i, DominantSet) && set_subset(common,Aicd[i-1])){
      adj=False;
    }
  }
_L99:;
  return adj;
} 

void WriteInputAdjacencyFile(ostream &f)
{
  RayRecord *RayPtr1, *RayPtr2;
  rowrange i,k, r1, r2, elem;
  colrange j;
  long pos1, pos2, degree, scard;
  boolean adj, redundant;
  node *headnode, *tailnode, *newnode, *prevnode;
  Aincidence Aicd;
  rowset Ared;  /* redundant inequality set */
  rowset Adom;  /* dominant inequality set */

  WriteProgramDescription(f);
  for(i=1; i<=mm; i++) set_initialize(&(Aicd[i-1]),RayCount); 
  set_initialize(&Ared, mm); 
  set_initialize(&Adom, mm); 
  headnode=NULL; tailnode=NULL;
  switch (Conversion) {
  case IneToExt:
    (f) << "*Adjacency List of input (=inequalities/facets)\n";
    break;

  case ExtToIne:
    (f) << "*Adjacency List of input (=vertices/rays)\n";
    break;

  default:
    break;
  }
  (f) <<  "*cdd input file : " << inputfile << " (" << minput << " x " << ninput << ")\n";
  (f) <<  "*cdd output file: " << outputfile << "\n";
  if (RayCount==0){
    goto _L99;
  }
  LastRay->Next=NULL;
  for (RayPtr1=FirstRay, pos1=1;RayPtr1 != NULL; RayPtr1 = RayPtr1->Next, pos1++){
    scard=set_card(RayPtr1->ZeroSet);
    elem = 0; i = 0;
    while (elem < scard && i < mm){
      i++;
      if (set_member(i,RayPtr1->ZeroSet)) {
        set_addelem(Aicd[i-1],pos1);
        elem++;
      }
    }
  }
  for (i=1; i<=mm; i++){
    if (set_card(Aicd[i-1])==RayCount){
      f << "*row " << i << " is dominating.\n";
      if (DynamicWriteOn) cout << "*row " << i << " is dominating.\n";
      set_addelem(Adom, i);
    }  
  }
  for (i=mm; i>=1; i--){
    if (set_card(Aicd[i-1])==0){
      redundant=True;
      f << "*row " << i << " is redundant;dominated by all others.\n";
      if (DynamicWriteOn) cout << "*row " << i << " is redundant;dominated by all others.\n";
      set_addelem(Ared, i);
    }else {
      redundant=False;
      for (k=1; k<=mm; k++) {
        if (k!=i && !set_member(k, Ared)  && !set_member(k, Adom) && 
            set_subset(Aicd[i-1], Aicd[k-1])){
          if (!redundant){
            f << "*row " << i << " is redundant;dominated by:"; 
            if (DynamicWriteOn) cout << "*row " << i << " is redundant;dominated by:"; 
            redundant=True;
          }
          f << " " << k;
          if (DynamicWriteOn) cout << " " << k;
          set_addelem(Ared, i);
        }
      }
      if (redundant){
        f << "\n";
        if (DynamicWriteOn) cout << "\n";
      }
    }
  }

  (f) << "begin\n";
  (f) << "  " << mm << "\n";
  for (i=1; i<=mm; i++){
    degree=0;
    for (j=1; j<=mm; j++){
      if (i!=j && InputAdjacentQ(Aicd, Ared, Adom, i, j)) {
        degree++;
        if (degree==1){
          newnode=new node;
          newnode->key=j;
          newnode->next=NULL;
          headnode=newnode;
          tailnode=newnode;
        }
        else{
          newnode=new node;
          newnode->key=j;
          newnode->next=NULL;
          tailnode->next=newnode;
          tailnode=newnode;
        }
      }
    } /* end of j */
    (f) << " " << i << " " << degree << " :";
    for (newnode=headnode; newnode!=NULL; newnode=newnode->next, delete prevnode){
      prevnode=newnode;
      (f) << " " << newnode->key;
    }
    headnode=NULL;
    (f) << "\n";
  } /* end of i */
  (f) << "end\n";

// Input Incidence information won't be output if the following lines
// are commented out 
//  (f) << "\n*Incidences of input and output (dual of incidence information)\n";
//  (f) << "*After <begin> three numbers are m1, m and output_size,\n";
//  (f) << "*where m1 is m+1 (for vertex/ray enumeration) or m (for convex hull).\n";
//  (f) << "begin\n";
//  (f) << "  " << mm << "  " << minput << "  " << RayCount << "\n";
//  for (i=1; i<= mm ; i++) {
//    WriteInputIncidence(f, Aicd, i);
//  }
//  (f) << "end\n";
// up to here.

_L99:;
  for(i=1; i<=mm; i++) set_free(&(Aicd[i-1]));  
}

void WriteAdjacencyFile(ostream &f)
{
  RayRecord *RayPtr1, *RayPtr2;
  long pos1, pos2, degree;
  boolean adj;
  node *headnode, *tailnode, *newnode, *prevnode;

  WriteProgramDescription(f);
  headnode=NULL; tailnode=NULL;
  switch (Conversion) {
  case IneToExt:
    (f) << "*Adjacency List of output (=vertices/rays)\n";
    break;

  case ExtToIne:
    (f) << "*Adjacency List of output (=inequalities=facets)\n";
    break;

  default:
    break;
  }
  (f) <<  "*cdd input file : " << inputfile << " (" << minput << " x " << ninput << ")\n";
  (f) <<  "*cdd output file: " << outputfile << "\n";
  (f) << "begin\n";
  (f) << "  " << RayCount << "\n";
  if (RayCount==0){
    goto _L99;
  }
  LastRay->Next=NULL;
  for (RayPtr1=FirstRay, pos1=1;RayPtr1 != NULL; RayPtr1 = RayPtr1->Next, pos1++){
    for (RayPtr2=FirstRay, pos2=1,degree=0; RayPtr2 != NULL; RayPtr2 = RayPtr2->Next, pos2++){
      if (RayPtr1!=RayPtr2){
        CheckAdjacency2(&RayPtr1, &RayPtr2, &adj);
        if (adj) {
          degree++;
          if (degree==1){
            newnode=(node *)malloc(sizeof *newnode);
            newnode->key=pos2;
            newnode->next=NULL;
            headnode=newnode;
            tailnode=newnode;
          }
          else{
            newnode=(node *)malloc(sizeof *newnode);
            newnode->key=pos2;
            newnode->next=NULL;
            tailnode->next=newnode;
            tailnode=newnode;
          }
        }
      }
    }
    (f) << " " << pos1 << " " << degree << " :";
    for (newnode=headnode; newnode!=NULL; newnode=newnode->next, free(prevnode)){
      prevnode=newnode;
      (f) << " " << newnode->key;
    }
    (f) << "\n";
  }
_L99:;
  (f) << "end\n";
}

void WriteIncidenceFile(ostream &f)
{
  RayRecord *TempPtr;

  WriteProgramDescription(f);
  switch (Conversion) {
    case IneToExt:
       (f) << "*Incidences of output(=vertices/rays) and input (=hyperplanes)\n";
       (f) << "*   for each output, #incidence and the set of hyperplanes containing it\n";
       (f) << "*   or its complement with its cardinality with minus sign\n";
      break;
    case ExtToIne:
       (f) << "*Incidences of output(=facets) and input (=points)\n";
       (f) << "*   for each output, #incidence and the set of points lying on it\n";
       (f) << "*   or its complement with its cardinality with minus sign\n";
      break;

    default:
      break;
  }
  (f) << "*cdd input file : " << inputfile << "  (" << minput << " x " << ninput << ")\n";
  (f) << "*cdd output file: " << outputfile << "\n";
  (f) << "*After <begin>, three numbers are output_size, m and m1,\n";
  (f) << "*where m1 is m+1 (for vertex/ray enumeration) or m (for convex hull).\n";
  (f) << "begin\n";
  (f) << "  " << FeasibleRayCount << "  " << minput << "  " << mm << "\n";
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteIncidence(f, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  (f) << "end\n";
}

void WriteExtFile(ostream &f, ostream &f_log)
{
  RayRecord *TempPtr;

  time(&endtime);
  WriteProgramDescription(f);
  (f) << "*Input File:" << inputfile << "(" << minput << "x" << ninput << ")\n";
  WriteRunningMode(f);
  WriteCompletionStatus(f);
  WriteTimes(f);
  if (DynamicWriteOn){
    WriteCompletionStatus(cout);
    WriteTimes(cout);
  }
  if (Conversion == ExtToIne)
    (f) << "*Since hull computation is chosen, the output is a minimal inequality system\n";
  (f) << "*FINAL RESULT:\n";
  if (DynamicWriteOn)
    printf("*Computation complete.\n");
  if (Conversion == IneToExt) {
    if (Inequality==ZeroRHS && Conversion == IneToExt){
      (f)<< "*Number of Rays =" << FeasibleRayCount << "\n";
      (f)<< "*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n";
      if (DynamicWriteOn){
        cout << "*Number of Rays =" << FeasibleRayCount << "\n";
        cout << "*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n";
      }
    }else{
      if (DynamicWriteOn)
        cout << "*Number of Vertices =" << VertexCount << ", Rays =" 
          << (FeasibleRayCount - VertexCount) << "\n";
      (f)<< "*Number of Vertices =" << VertexCount << ", Rays =" 
        << (FeasibleRayCount - VertexCount) << "\n";
    }
  } else {
    if (DynamicWriteOn)
      cout << "*Number of Facets =" << FeasibleRayCount << "\n";
    (f)<< "*Number of Facets =" << FeasibleRayCount << "\n";
  }
  (f)<< "begin\n";
  switch (Inequality) {
  case ZeroRHS:
    (f)<<  FeasibleRayCount << "  " << (nn + 1) << "  " << OutputNumberString << "\n";
    break;
  case NonzeroRHS:
    (f)<<  FeasibleRayCount << "  " << nn <<  "  " << OutputNumberString << "\n";
    break;
  }

 if (DynamicWriteOn)
   cout << "*Writing the output file " << outputfile << "...\n";
  
TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(f, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  (f)<< "end\n";
  if (Conversion == IneToExt) (f)<< "hull\n";

// write the final comments on the log file.
  (f_log) << "end\n";
  WriteRunningMode(f_log);
  WriteCompletionStatus(f_log);
  if (PreOrderedRun) {
    (f_log) << "*set_intersection total#, effective#, loss# = "
      <<  count_int << "  " << count_int_good << "  " << count_int_bad << "\n";
  }
  WriteTimes(f_log);
}

void WriteDecompResult(ostream &f, ostream &f_log)
{
  RayRecord *TempPtr;
  static int callnumber=0;

  callnumber++;
  if (callnumber==1) {
    WriteProgramDescription(f);
    (f) << "*Input File:" << inputfile << "(" <<minput << "x" << ninput <<")\n";
    if (Conversion == ExtToIne)
      (f) << "*Since hull computation is chosen, the output is a minimal inequality system\n";
  }
  (f) << "begin\n";
  switch (Inequality) {
  case ZeroRHS:
    (f) << " " << FeasibleRayCount << "  " << (nn + 1)  << "  " << OutputNumberString << "\n";
    break;
  case NonzeroRHS:
    (f) << " " << FeasibleRayCount << "  " << nn  << "  " << OutputNumberString << "\n";
    break;
  }
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    if (TempPtr->feasible) {
      WriteRayRecord(f, TempPtr);
    }
    TempPtr = TempPtr->Next;
  }
  (f) << "end\n";
  if (LogWriteOn) {
    (f_log) << "end\n";
    WriteCompletionStatus(f_log);
    WriteTimes(f_log);
  }
  if (DynamicWriteOn) {
    WriteCompletionStatus(cout);
    WriteTimes(cout);
  }
}


void WriteProjResult(ostream &f, ostream &f_log, long *dbrow)
{
  RayRecord *TempPtr;

  time(&endtime);
  WriteProgramDescription(f);
  f << "*Input File: " << inputfile << "  (" << minput << " x " << ninput << ")\n";
  WriteRunningMode(f);
  WriteCompletionStatus(f);
  WriteTimes(f);
  f << "*FINAL RESULT:\n";
  if (DynamicWriteOn)
    printf("*Computation complete.\n");
  if (DynamicWriteOn)
     cout << "*Number of Vertices =" << VertexCount << ",   Rays =" 
       <<  (FeasibleRayCount-VertexCount) << "\n";
  f << "*Number of Vertices =" << VertexCount << ",   Rays =" 
     <<  (FeasibleRayCount-VertexCount) << "\n";
  f << "begin\n";
  f << " " << FeasibleRayCount << "  " << (mm+1) << " " << OutputNumberString <<"\n";
  TempPtr = FirstRay;
  while (TempPtr != NULL) {
    WriteProjRayRecord(f, TempPtr, dbrow);
    TempPtr = TempPtr->Next;
  }
  f << "end\n";
  if (DynamicWriteOn) {
    WriteCompletionStatus(cout);
    WriteTimes(cout);
  }
  f_log << "end\n";
  WriteRunningMode(f_log);
  WriteCompletionStatus(f_log);
  WriteTimes(f_log);
}

void InitialWriting(ostream &f, ostream &f_log)
{
  WriteProgramDescription(f_log);
  (f_log) << "*Input File: " << inputfile << "(" <<minput << " x " << ninput << ")\n";
  (f_log) << "*Initial set of hyperplanes: ";
  set_fwrite(f_log, AddedHyperplanes);
  (f_log) << "\n";
  
  (f_log) << "begin\n";
  (f_log) << " " << (mm-nn) << "  " <<  5 << "\n";
}

void WriteErrorMessages(ostream &f)
{
  switch (Error) {

  case LowColumnRank:
    if (Conversion==IneToExt) {
      (f)<< "*Input Error: Input matrix (b, -A) is not column full rank => no vertices.\n";
      cout<< "*Input Error: Input matrix (b, -A) is not column full rank => no vertices.\n";
      break;
    } else {
      (f)<< "*Input Error: The polytope (convex hull) is not full dimensional.\n";
      cout<< "*Input Error: The polytope (convex hull) is not full dimensional.\n";
      break;
    }
 
  case DimensionTooLarge:
    (f)<< "*Input Error: Input matrix is too large:\n";
    (f)<< "*Please increase MMAX and/or NMAX in the source code and recompile.\n";
    cout<< "*Input Error: Input matrix is too large:\n";
    cout<< "*Please increase MMAX and/or NMAX in the source code and recompile.\n";
    break;

  case DependentMarkedSet:
    (f)<< "*Input Error: Active (equality) rows are linearly dependent.\n";
    (f)<< "*Please select independent rows for equality restriction.\n";
    cout<< "*Input Error: Active (equality) rows are linearly dependent.\n";
    cout<< "*Please select independent rows for equality restriction.\n";
    break;

  case FileNotFound:
    (f)<< "*Input Error: Specified input file does not exist.\n";
    cout<< "*Input Error: Specified input file does not exist.\n";
    break;

  case ImproperExecutable:
    if (Number == Rational) {
      (f) << "*Improper Executable Error: the current executable cannot handle rational data.\n";
      (f) << "*Please use cdd+ (cddr+) compiled with -DRATIONAL compiler option.\n";
      cout << "*Improper Executable Error: the current executable cannot handle rational data.\n";
      cout << "*Please use cdd+ (cddr+) compiled with -DRATIONAL compiler option.\n";
    } else {  /* Number == Real */
      (f)<< "*Improper Executable Error: the current execubtale cannot handle floating-point data.\n";
      (f)<< "*Please use cdd+ (cddf+) compiled without -DRATIONAL compiler option.\n";
      cout << "*Improper Executable Error: the current execubtale cannot handle floating-point data.\n";
      cout << "*Please use cdd+ (cddf+) compiled without -DRATIONAL compiler option.\n";
    }
    break;
 
  case ImproperInputFormat:
    (f)<< "*Input Error: Input format is not correct.\n";
    (f)<< "*Format:\n";
    (f)<< " begin\n";
    (f)<< "   m   n  NumberType(real, rational or integer)\n";
    (f)<< "   b  -A\n";
    (f)<< " end\n";
    cout<< "*Input Error: Input format is not correct.\n";
    cout<< "*Format:\n";
    cout<< " begin\n";
    cout<< "   m   n  NumberType(real, rational or integer)\n";
    cout<< "   b  -A\n";
    cout<< " end\n";
    break;

  case None:
    (f)<< "*No Error found.\n";
    cout<< "*No Error found.\n";
    break;
  }
}

void WriteLPResult(ostream &f, LPStatusType LPS, myTYPE optval,
   Arow sol, Arow dsol, colindex NBIndex, rowrange re, colrange se,
   long iter)
{
  long j;

  time(&endtime);
  WriteProgramDescription(f);
  (f) << "*cdd LP Result\n";  
  (f) << "*cdd input file : " << inputfile << "  (" << minput << " x " << ninput << ")\n";
  switch (LPsolver){
    case DualSimplex:
      (f) << "*LP solver: Dual Simplex\n"; break;

    case CrissCross:
      (f) << "*LP solver: Criss-Cross Method\n"; break;

    case CombMaxImprove:
      (f) << "*LP solver: Combinatorial Maximum Improvement Algorithm\n"; break;

  }
    (f) << "*LP status: a dual pair (x, y) of optimal solutions found.\n";
  if (Conversion==LPmax)
    (f) << "*maximization is chosen.\n";  
  else if (Conversion==LPmin)
    (f) << "*minimization is chosen.\n";
  else if (Conversion==InteriorFind){
    (f) << "*inerior point computation is chosen.\n";
    (f) << "*the following is the result of solving the LP:\n";
    (f) << "*   maximize      x_{d+1}\n";
    (f) << "*   s.t.    A x + x_{d+1}  <=  b.\n";
    (f) << "*Thus, the optimum value is zero     if the polyhedron has no interior point.\n";
    (f) << "*      the optimum value is negative if the polyhedron is empty.\n";
    (f) << "*      the LP is dual inconsistent   if the polyhedron admits unbounded inscribing balls.\n";
  }
  if (Conversion==LPmax||Conversion==LPmin){
    (f) << "*Objective function is\n";  
    for (j=0; j<nn; j++){
      if (j>0 && AA[OBJrow-1][j]>=0 ) (f) << " +";
      if (j>0 && (j % 5) == 0) (f) << "\n";
      WriteNumber(f, AA[OBJrow-1][j]);
      if (j>0) (f) << " X[" << j << "]";
    }
    (f) << "\n";
  }

  switch (LPS){
  case Optimal:
    (f) << "*LP status: a dual pair (x, y) of optimal solutions found.\n";
    (f) << "begin\n";
    (f) << "  primal_solution\n";
    for (j=1; j<nn; j++) {
      (f) << "  " << j << " : ";
      WriteNumber(f,sol[j]);
      (f) << "\n";
    }
    (f) << "  dual_solution\n";
    for (j=1; j<nn; j++){
      (f) << "  " << NBIndex[j+1] << " : ";
      WriteNumber(f,dsol[j]);
      (f) << "\n";
    }
    (f) << "  optimal_value : " << optval << "\n";
    (f) << "end\n";
    break;

  case Inconsistent:
    (f) << "*LP status: LP is inconsistent.\n";
    (f) << "*The nonnegative combination of original inequalities with\n";
    (f) << "*the following coefficients will prove the inconsistency.\n";
    (f) << "begin\n";
    (f) << "  dual_direction\n";
    (f) << "  " << re << " : ";
    WriteNumber(f,1); 
    (f) << "\n";
    for (j=1; j<nn; j++){
      (f) << "  " << NBIndex[j+1] << " : ";
      WriteNumber(f,dsol[j]);
      (f) << "\n";
    }
    (f) << "end\n";
    break;

  case DualInconsistent:
    (f) << "LP status: LP is dual inconsistent.\n";
    (f) << "*The linear combination of columns with\n";
    (f) << "*the following coefficients will prove the dual inconsistency.\n";
    (f) << "*(It is also an unbounded direction for the primal LP.)\n";
    (f) << "begin\n";
    (f) << "  primal_direction\n";
    for (j=1; j<nn; j++) {
      (f) << "  "<< j << " : ";
      WriteNumber(f,sol[j]);
      (f) << "\n";
    }
    (f) << "end\n";
    break;

  default:
    break;
  }
  (f) << "*number of pivot operations = " << iter << "\n";
  WriteTimes(f);
}


void WriteSolvedProblem(ostream &f)
{
  (f) << "*cdd input file : " << inputfile << "  ( " << minput << " x " << ninput << ")\n";
  (f) << "*The input data has been interpreted as the following.\n";
  WriteAmatrix(f,AA,minput,nn,Inequality);
  WriteRunningMode2(f);
}

void AmatrixInput(boolean *successful)
{
  long i,j;
  myTYPE value=0;
  long value1,value2;
  boolean found=False,decided=False,fileopened, localdebug=False;
  string numbtype, stemp; 
  char ch;
  string command="", inputst="";

  *successful = False;
  SetInputFile(&fileopened);

  if (!fileopened){
    Error=FileNotFound;
  } 
  else {
    ifstream reading(inputfile);
    while (!reading.eof() && inputst!="begin")
    {
      reading >> inputst;
    }
    reading >> minput;
    reading >> ninput;
    reading >> numbtype;
    cout << "size = " << minput << " x " << ninput << "\nnumber type = " << numbtype << "\n";
    SetNumberType(numbtype);
    if (Error==ImproperExecutable) {
      goto _L99;
    } 
    Inequality=ZeroRHS;
    for (i=1; i<=minput && !decided; i++) {
      reading >> value;
      if (FABS(value) > zero) {
        Inequality = NonzeroRHS;
        decided=True;
      }
      for (j=2; j<=ninput; j++) {
        reading >> value;
      }
      if (localdebug) printf("remaining data to be skipped:");
      while (reading.get(ch) && ch != '\n') {
        if (localdebug) cout << ch;
      }
      if (localdebug) putchar('\n');
    }
    if (Inequality==NonzeroRHS) {
      printf("Nonhomogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
      nn = ninput;
    }
    else {
      printf("Homogeneous system with  m = %5ld  n = %5ld\n", minput, ninput);
      nn = ninput-1;
    }
    if (nn > NMAX || minput > MMAX) {
      Error = DimensionTooLarge;
      goto _L99;
    }
  
    RestrictedEnumeration = False;
    RelaxedEnumeration = False;
    EqualityIndex=(long *)calloc(minput+2, sizeof *EqualityIndex);
    for (i = 0; i <= minput+1; i++) EqualityIndex[i]=0;

    while (!reading.eof()) {
      reading >> command;
      if (debug) cout << command << "\n";
      ProcessCommandLine(reading,command);
    } 
    reading.close();

    reading.open(inputfile); 
    found=False;
    while (!found)
    {
      if (!reading.eof()) {
        reading >> command;
  	if (command=="begin") {
          found=True;
  	}
      }
      else {
  	Error=ImproperInputFormat;
  	goto _L99;
      }
    }
    reading >> value1;
    reading >> value2;
    reading >> stemp;
    for (i = 1; i <= minput; i++) {
      AA[i-1]=new myTYPE[ninput];
      for (j = 1; j <= ninput; j++) {
        reading >> value;
        if (Inequality==NonzeroRHS) 
      	  AA[i-1][j - 1] = value;
        else if (j>=2) {
          AA[i-1][j - 2] = value;
	}
	if (debug) cout << "a(" << i << "," << j << ") = " << value << "\n";
      }  /*of j*/
      if (debug) printf("comments to be skipped:");
      while (reading.get(ch) && ch != '\n') {
        if (debug) cout << ch;
      }
      if (debug) putchar('\n');
    }  /*of i*/
    if (reading.eof()) {
      Error=ImproperInputFormat;
      goto _L99;
    }
    else{
      reading >> command;
      if (command!="end") {
        if (debug) cout << "'end' missing or illegal extra data:" << command << "\n";
   	  Error=ImproperInputFormat;
  	   goto _L99;
      }
    }
  
    switch (Conversion) {
    case IneToExt:
      if (Inequality==NonzeroRHS){
        mm = minput + 1;
        AA[mm-1]= new myTYPE[ninput];
        for (j = 1; j <= ninput; j++) {   /*artificial row for x_1 >= 0*/
        if (j == 1)
          AA[mm - 1][j - 1] = 1;
        else
          AA[mm - 1][j - 1] = 0;
        }
      } else {
        mm = minput;
      }
    break;

    case ExtToIne:
      mm = minput;
      break;

    case LPmax:  case LPmin:
      mm = minput + 1;
      OBJrow=mm;
      RHScol=1L;
      AA[mm-1]= new myTYPE[ninput];
      for (j = 1; j <= ninput; j++) {   /*objective row */
 	 AA[mm - 1][j - 1] = LPcost[j-1];
      }
      break;

    default:
      mm = minput;
    }
    SetInequalitySets(EqualityIndex);
    *successful = True;
    _L99:
    reading.close();
  }
}

void WriteProjRayRecord(ostream &f, RayRecord *RR, long *dbrow)
{
  long i,j,k;
  static rowset dbset;
  myTYPE* vec;

  set_initialize(&dbset,mm);  
  vec= new myTYPE[mm];
  for (j = 1; j <= mm-nn; j++){
    i=dbrow[j];
    set_addelem(dbset,i);
    if (debug)  cout << "index " << i << " is added to dbset\n";
    vec[i-1]=0;
    for (k=1; k<=nn; k++) {
      vec[i-1]+= (RR->Ray[k-1])*AA[j-1][k-1];
    }
    if (debug) cout << "vec[" << i-1 << "]= " << vec[i-1] << "\n";
  }
  i=1;
  for (j = 1; j <= mm; j++){
    if (!set_member(j,dbset)){
      vec[j-1]=RR->Ray[i-1];
      i++;
    }
  }
  (f)<< 0;
  for (j = 0; j < mm; j++)
    WriteNumber(f, vec[j]);
  (f).put('\n');
}


void WriteRowEquivalence(ostream &f, long classno, rowindex rowequiv)
{
  rowrange i,k;
  long cn=0;
 
  f << "*Classification of rows (rows in the same class are equivalent)\n"; 
  for (i=1; i<=mm; i++)
  {
    if (cn<rowequiv[i]) {
      cn=rowequiv[i];
      f << "class " << cn << ": " << i ;
      for (k=i+1; k<=mm; k++){
        if (rowequiv[k]==cn) f << ", " << k;
      }
      f << "\n";
    } 
  }
}

/* end of cddio.c */

