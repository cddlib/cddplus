/* cddio.C:  Basic Input and Output Procedures for cdd.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C++-Implementation of the double description method for
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
      strcpy(ifiletail, &(inputfile[dotpos]));
    }else{
      strcpy(ifilehead, inputfile);
    }
    if (debug){
      printf("inputfile name: %s\n", inputfile);  
      printf("inputfile name head: %s\n", ifilehead);  
      printf("inputfile name tail: %s\n", ifiletail);  
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
  DataFileType newname;
  
  switch (cflag) {
    case 'o':
      switch (Conversion) {
        case ExtToIne:
          extension=".ine"; break;     /* output file for ine data */
        case IneToExt:   case Projection:
          extension=".ext"; break;     /* output file for ext data */
        case FacetListing:  case FacetListingExternal:
          extension=".fis"; break;     /* output file for facet_listing */
        case VertexListing: case VertexListingExternal:
          extension=".vis"; break;     /* output file for vertex_listing */
        case TopeListing:
          extension=".tis"; break;     /* output file for tope_listing */
        case LPmax:  case LPmin:  case InteriorFind:
          extension=".lps"; break;     /* output file for LPmax, LPmin, interior_find */
        default:
        extension=".xxx";break;
      }
      break;

    case 'a':         /* decide for output adjacence */
      if (Conversion==IneToExt)
        extension=".ead";       /* adjacency file for ext data */
      else
        extension=".iad";       /* adjacency file for ine data */
      break;
    case 'i':         /* decide for output incidence */
      if (Conversion==IneToExt)
        extension=".ecd";       /* ext incidence file */
      else
        extension=".icd";       /* ine incidence file */
      break;
    case 'n':         /* decide for input incidence */
      if (Conversion==IneToExt)
        extension=".icd";       /* ine incidence file */
      else
        extension=".ecd";       /* ext incidence file */
      break;
    case 'j':        /* decide for input adjacence */
      if (Conversion==IneToExt)
        extension=".iad";       /* ine adjacency file */
      else
        extension=".ead";       /* ext adjacency file */
      break;
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
  ifstream reading_tmp(fname);
  if (reading_tmp.is_open()) {
    if (!PostAnalysisOn || cflag !='o') {  
    /* rename the file only when it is not for postanalysis inputfiles nor external file */
      strcpy(newname,fname);  
      strcat(newname,".old");
      rename(fname,newname);
      if (DynamicWriteOn){
        printf("Default output file %s exists.\n",fname);
        printf("Caution: The old file %s is renamed to %s!\n",fname,newname);
        printf("Create a new file %s.\n",fname);
      }
    }
  }
  reading_tmp.close();
  if (DynamicWriteOn) printf("Open %s file %s.\n",fscript,fname);
}

void SetReadFileName(DataFileType fname, char cflag, char *fscript)
{
  boolean quit=False;
  char *extension;
  DataFileType newname;
  
  switch (cflag) {
    case 'o':
      if (Conversion==ExtToIne)
        extension=".ine";       /* output file for ine data */
      else
        extension=".ext";       /* output file for ext data, and for any other conversions */
      break;
    case 'x':
      if (Conversion==FacetListingExternal)
        extension=".ine.external";       /* external ine file */
      else
        extension=".ext.external";       /* external ext file */
      break;
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
    // the executable cannot handle rational data so it converts to real 
    // modified on 1996-02-18 
    Number = Real;
    InputNumberString="rational";
    OutputNumberString="real";
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
  static long var,msize, outdig=-1;
  myTYPE cost=0, zero_input=-1, purezero=0;

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
    if (IncidenceOutput==InputIncidence || IncidenceOutput==IOIncidence)
      IncidenceOutput=IOIncidence;
    else
      IncidenceOutput = OutputIncidence;
    return;
  }
  if (line== "#incidence") {
    IncidenceOutput = IncCardinality;
    return;
  }
  if (line== "input_incidence") {
    if (IncidenceOutput==OutputIncidence || IncidenceOutput==IOIncidence)
      IncidenceOutput=IOIncidence;
    else
      IncidenceOutput = InputIncidence;
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
  if ((line== "partial_enumeration" || line=="equality" || line=="linearity") 
    && RestrictedEnumeration==False) {
    (f) >> msize;
    for (j=1;j<=msize;j++) {
      (f) >> var;
      EqualityIndex[var]=1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Projection is cancelled because it cannot be compatible with Partial Enumeration.\n");
      Conversion=IneToExt;
    }
    RestrictedEnumeration=True;
    return;
  }
  if (line== "strict_inequality" && RelaxedEnumeration==False) {
    (f) >> msize;
    for (j=1;j<=msize;j++) {
      (f) >> var;
      EqualityIndex[var]=-1;
    }
    printf("\n");
    if (Conversion==Projection) {
      printf("Warning: Projection is cancelled because it is not compatible with Partial Enumeration.\n");
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
  if (line== "solve_rowsubproblem" && SpecialConversion!=RowSubproblemSolve) {
    set_initialize(&SubproblemRowSet,minput+1);
    set_initialize(&RedundantRowSet,minput+1);
    if (debug) printf("solve_rowsubproblem is chosen.\n");

    (f) >> msize;
    for (j=1;j<=msize;j++) {
      (f) >> var;
      set_addelem(SubproblemRowSet,var);
    }
    SpecialConversion=RowSubproblemSolve;
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
  if (line== "facet_listing" && Conversion != LPmax && Conversion != FacetListing) {
    if (debug) printf("facet_listing is chosen.\n");
    Conversion=FacetListing;
    return;
  }
  if (line== "facet_listing_external" && Conversion != LPmax && Conversion != FacetListing) {
    if (debug) printf("facet_listing_external is chosen.\n");
    Conversion=FacetListingExternal;
    return;
  }
  if (line== "vertex_listing" && Conversion != LPmax && Conversion != VertexListing) {
    if (debug) printf("vertex_listing is chosen.\n");
    Conversion=VertexListing;
    return;
  }
  if (line== "vertex_listing_external" && Conversion != LPmax && Conversion != VertexListingExternal) {
    if (debug) printf("vertex_listing_external is chosen.\n");
    Conversion=VertexListingExternal;
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
  if (line== "row_decomposition" && SpecialConversion!=RowDecomposition) {
    if (DynamicWriteOn) printf("Row decomposition is chosen.\n");
    SpecialConversion=RowDecomposition;
    return;
  }
  if (line== "round_output_off" && Round_Output) {
    if (DynamicWriteOn) printf("No rounding of output numbers.\n");
    Round_Output=False;
    return;
  }
  if (line== "zero_tolerance" && zero_input<purezero) {
    (f) >> zero_input;
    if (zero_input > purezero){
      zero=zero_input;
      if (DynamicWriteOn)  cout << "zero_tolerance is reset to " << zero_input << ".\n";
    } 
    else cout << "Warning: float_zero must be positive. Use the default =" << zero
      << ".\n";
    return;
  }
  if (line== "output_digits" && outdig==-1) {
    (f) >> outdig;
    if (outdig > 0){
      output_digits=outdig;
      if (DynamicWriteOn) cout << "output_digits is reset to " << output_digits << ".\n";
    } 
    else {
      if (DynamicWriteOn) cout << "Warning: output_digits must be >= 1. Use the default =" 
      << output_digits << ".\n";
    }
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
      if (RR->Ray[0] > zero) {
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

  (f) << "H-representation\n" << "begin\n";
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
    (f) << "  " << rowmax << "  " << colmax+1 <<  "  " << OutputNumberString << "\n";
  else
    (f) << "  " << rowmax << "  " << colmax  << "  " << OutputNumberString << "\n";
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
  long ix1,ix2,ix, i;
  static myTYPE sig=100000.;

  if (Round_Output){
    ix1= (long) (FABS(x) * sig + 0.5);
    ix2= (long) (FABS(x) + 0.5);
    ix2= (long) (ix2* (long)sig);
    if (ix1 == ix2) {
      if (x>0) {
        ix = (long) (x + 0.5);
      } else {
        ix = (long) (-x + 0.5);
        ix = -ix;
      }
      (f) << " " << ix;
    } else {
      (f).setf(ios::scientific,ios::floatfield);
      (f).precision(output_digits);
      (f) << " " << x;
    }
  } else {
    (f).setf(ios::scientific,ios::floatfield);
    (f).precision(output_digits);
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

void WriteDictionary(ostream &f, Amatrix X, Bmatrix T,
 rowindex OV, long bflag[], colindex nbindex, rowrange objrow, colrange rhscol)
/* Write the sign matrix of tableau  X.T.
   This works properly only for Nonhomogeneous inequality  */
{
  colrange j;
  rowrange i,k,l;
  
  (f) << "\\multicolumn{2}{r}{g}";
  for (j=2; j<= nn; j++) {
    (f).width(3);
    (f) << "&\\multicolumn{1}{r}{" << nbindex[j]<< "}\n";
  }
  (f) << "\\\\ \\cline{2-"<< (nn+1) << "}\n";
  (f) << "     g|";
  for (j=2; j<= nn; j++) {
    (f).width(3);
    (f) << " " << nbindex[j];
  }
  (f) << "\\\\\n";
  for (i=1; i<= mm; i++) {
    if (bflag[i]==0) {  /* i the objective variable */
      (f) << "   f";
      (f) << " &" << TableauEntry(X,T,i,1) << " ";   
      for (j=2; j<= nn; j++) {
        (f).width(3);
        (f) << "&" << TableauEntry(X,T,i,j);
      }
      (f) << "\\\\ \\cline{2-" << (nn+1) << "}\n";
    }
  }
  (f) << "  -----+";
  for (j=2; j<=nn; j++) { 
    (f) << "-----";
  }
  (f) << "\\\\\n";
  for (i=1; i<= mm; i++) {
    if (bflag[i]!=0 && bflag[i]==-1) {  /* i is a basic variable */
      (f).width(2); 
      (f) << " " << i;   
      (f).width(3);
      (f) << "&" << TableauEntry(X,T,i,1) << " ";   
      for (j=2; j<= nn; j++) {
        (f).width(3);
        (f) << "&" << TableauEntry(X,T,i,j);
      }
      (f) << "\\\\\n";
    }
  }
  (f) << "\\cline{2-" << (nn+1) << "}\n";
  (f) << "\\multicolumn{" << nn+1 << "}{r}{}\\\\\n";
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

  case OutputIncidence: case IOIncidence: 
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
  }
  if (NondegAssumed) {
    (f) << "*Degeneracy preknowledge for computation: NondegenerateAssumed\n";
  }
  else {
    (f) << "*Degeneracy preknowledge for computation: None (possible degeneracy)\n";
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

    case VertexListing:
      (f) << "*Vertex listing is chosen.\n";
      break;

    case FacetListingExternal:
      (f) << "*Facet listing external is chosen.\n";
      break;

    case VertexListingExternal:
      (f) << "*Vertex listing external is chosen.\n";
      break;

    case TopeListing:
      (f) << "*Tope listing is chosen.\n";
      break;

    case Projection:
      (f) << "*Preprojection is chosen.\n";
      break;
  
    default: break;
  }
  switch (AdjacencyOutput) {
    case IOAdjacency:
      (f) <<  "*Output adjacency file is requested.\n";
      (f) <<  "*Input adjacency file is requested.\n";
      break;
    
    case InputAdjacency:
      (f) <<  "*Input adjacency file is requested.\n";
      break;

    case OutputAdjacency:
      (f) <<  "*Output adjacency file is requested.\n";
      break;

    default: break;
  }
  switch (IncidenceOutput) {
    case IOIncidence:
      (f) <<  "*Output incidence file is requested\n";
      (f) <<  "*Input incidence file is requested.\n";
      break;
    
    case InputIncidence:
      (f) <<  "*Input incidence file is requested.\n";
      break;

    case OutputIncidence:
      (f) <<  "*Output incidence file is requested.\n";
      break;

    case IncCardinality:
      (f) <<  "*Incidence cardinality output is requested.\n";
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
#ifndef RATIONAL
  (f) << "*Zero tolerance = " << zero << "\n";
#endif
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

    case VertexListing:
      (f) << "vertex_listing\n";
      break;

    case FacetListingExternal:
      (f) << "facet_listing_external\n";
      break;

    case VertexListingExternal:
      (f) << "vertex_listing_external\n";
      break;

    case TopeListing:
      (f) << "tope_listing\n";
      break;

    case Projection:
      (f) <<  "*Preprojection is chosen.\n";
      break;
  
    default: break;
  }
  switch (AdjacencyOutput) {
    case IOAdjacency:
      (f) <<  "adjacency\n";
      (f) <<  "input_adjacency\n";
      break;
    
    case InputAdjacency:
      (f) <<  "input_adjacency\n";
      break;

    case OutputAdjacency:
      (f) <<  "adjacency\n";
      break;

    default: break;
  }
  switch (IncidenceOutput) {
    case IOIncidence:
      (f) <<  "incidence\n";
      (f) <<  "input_incidence\n";
      break;
    
    case InputIncidence:
      (f) <<  "input_incidence\n";
      break;

    case OutputIncidence:
      (f) <<  "incidence\n";
      break;

    case IncCardinality:
      (f) <<  "#incidence\n";
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

void WriteRunningMode0(ostream &f)
{
  switch (Conversion) {
    case ExtToIne: case VertexListing: case VertexListingExternal:
      (f) << "V-representation\n";
      break;
    
    case IneToExt: case LPmax: case LPmin: case FacetListing: 
    case FacetListingExternal: case TopeListing: case Projection:
      (f) << "H-representation\n";
      break;

    default: break;
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
  static long lastRayCount=0;

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

_L99:;
  for(i=1; i<=mm; i++) set_free(&(Aicd[i-1]));  
}

void WriteInputIncidenceFile(ostream &f)
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
  (f) <<  "*cdd input file : " << inputfile << " (" << minput << " x " << ninput << ")\n";
  (f) <<  "*cdd output file: " << outputfile << "\n";

  switch (Conversion) {
  case IneToExt:
    (f) << "*Incidence of input (=inequalities/facets) w.r.t. output (=vertices/rays).\n";
    break;

  case ExtToIne:
    (f) << "**Incidence of input (=vertices/rays) w.r.t. output (=inequalities/facets).\n";
    break;

  default:
    break;
  }

  for(i=1; i<=mm; i++) set_initialize(&(Aicd[i-1]),RayCount); 
  set_initialize(&Ared, mm); 
  set_initialize(&Adom, mm); 
  headnode=NULL; tailnode=NULL;

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
      set_addelem(Adom, i);
    }  
  }
  for (i=mm; i>=1; i--){
    if (set_card(Aicd[i-1])==0){
      redundant=True;
      f << "*row " << i << " is redundant;dominated by all others.\n";
      set_addelem(Ared, i);
    }else {
      redundant=False;
      for (k=1; k<=mm; k++) {
        if (k!=i && !set_member(k, Ared)  && !set_member(k, Adom) && 
            set_subset(Aicd[i-1], Aicd[k-1])){
          if (!redundant){
            f << "*row " << i << " is redundant;dominated by:"; 
            redundant=True;
          }
          f << " " << k;
          set_addelem(Ared, i);
        }
      }
      if (redundant){
        f << "\n";
      }
    }
  }
  (f) << "begin\n";
  (f) << "  " << minput << "  " << RayCount;
  switch (Conversion) {
  case IneToExt:
    (f) << "  " << RayCount << "\n";
    break;

  case ExtToIne:
    (f) << "  " << RayCount+1 << "\n";
    break;

  default:
    break;
  }
  for (i=1; i<= mm ; i++) {
    WriteInputIncidence(f, Aicd, i);
  }
  (f) << "end\n";

_L99:;
  for(i=1; i<=mm; i++) set_free(&(Aicd[i-1]));  
}

void WriteAdjacencyFile(ostream &f)
{
  RayRecord *RayPtr1, *RayPtr2;
  long pos1, pos2, degree;
  boolean adj,localdebug=False;
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
  if (localdebug) cout << "  " << RayCount << "\n";
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
    if (localdebug) cout << " " << pos1 << " " << degree << " :";
    for (newnode=headnode; newnode!=NULL; newnode=newnode->next, free(prevnode)){
      prevnode=newnode;
      (f) << " " << newnode->key;
      if (localdebug) cout << " " << newnode->key;
    }
    (f) << "\n";
    if (localdebug) cout << "\n";
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
      (f)<< "*Number of Rays = " << FeasibleRayCount << "\n";
      (f)<< "*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n" <<"V-representation\n";
      if (DynamicWriteOn){
        cout << "*Number of Rays = " << FeasibleRayCount << "\n";
        cout << "*Caution!: the origin is a vertex, but cdd does not output this trivial vertex\n"<<"V-representation\n";
      }
    }else{
      if (DynamicWriteOn)
        cout << "*Number of Vertices = " << VertexCount << ", Rays =" 
          << (FeasibleRayCount - VertexCount) << "\n" << "V-representation\n";
      (f)<< "*Number of Vertices =" << VertexCount << ", Rays =" 
        << (FeasibleRayCount - VertexCount) << "\n" << "V-representation\n";
    }
  } else {
    if (DynamicWriteOn)
      cout << "*Number of Facets = " << FeasibleRayCount << "\n" << "H-representation\n";
    (f)<< "*Number of Facets = " << FeasibleRayCount << "\n" << "H-representation\n";
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
  f << "V-representation\n" << "begin\n";
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
      if (j>0 && AA[OBJrow-1][j]>=zero ) (f) << " +";
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
  WriteRunningMode0(f);
  WriteAmatrix(f,AA,minput,nn,Inequality);
  WriteRunningMode2(f);
}


void AmatrixInput(boolean *successful)
{
  long i,j;
  myRational rvalue=0;
  myTYPE value=0;    // modified on 1996-02-18 
  long value1,value2;
  boolean found=False,decided=False,fileopened, newformat=False, localdebug=False;
  string numbtype, stemp; 
  char ch;
  string command="", inputst="";

  *successful = False;
  SetInputFile(&fileopened);

  if (!fileopened){
    Error=FileNotFound;
  } 
  else {
     {
	ifstream reading(inputfile);
	while (!reading.eof() && inputst!="begin")
	   {
	      reading >> inputst;
	      if (inputst=="V-representation"){
		 Conversion = ExtToIne; newformat=True;
	      } else if (inputst=="H-representation"){
		 Conversion =IneToExt; newformat=True;
	      }
	   }
	reading >> minput;
	reading >> ninput;
	reading >> numbtype;
	cout << "size = " << minput << " x " << ninput << "\nnumber type = " << numbtype << "\n";
	SetNumberType(numbtype);
	if (Error==ImproperExecutable) {
	   reading.close();
	   return;
	} 
	Inequality=ZeroRHS;
	for (i=1; i<=minput && !decided; i++) {
	   if (InputNumberString=="rational" && OutputNumberString=="real"){
	      reading >> rvalue;
	      value=myTYPE(rvalue);
	   } else {
	      reading >> value;
	   }
	   if (FABS(value) > zero) {
	      Inequality = NonzeroRHS;
	      decided=True;
	   }
	   for (j=2; j<=ninput; j++) {
	      if (InputNumberString=="rational" && OutputNumberString=="real"){
		 reading >> rvalue;
	      } else {
		 reading >> value;
	      }
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
	   reading.close();
	   return;
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
     }

    ifstream reading(inputfile);
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
	reading.close();
	return;
      }
    }
    reading >> value1;
    reading >> value2;
    reading >> stemp;
    for (i = 1; i <= minput; i++) {
      AA[i-1]=new myTYPE[ninput];
      for (j = 1; j <= ninput; j++) {
        if (InputNumberString=="rational" && OutputNumberString=="real"){
          reading >> rvalue;
          value=myTYPE(rvalue);
        } else {
          reading >> value;
        }
        if (Inequality==NonzeroRHS) 
      	  AA[i-1][j - 1] = value;
        else if (j>=2) {
          AA[i-1][j - 2] = value;
	}
	if (debug){
           if (InputNumberString=="rational" && OutputNumberString=="real")
             cout << "a(" << i << "," << j << ") = " << value << " ("<< rvalue << ")\n";
           else
             cout << "a(" << i << "," << j << ") = " << value  << "\n";
         }
      }  /*of j*/
      if (debug) printf("comments to be skipped:");
      while (reading.get(ch) && ch != '\n') {
        if (debug) cout << ch;
      }
      if (debug) putchar('\n');
    }  /*of i*/
    if (reading.eof()) {
      Error=ImproperInputFormat;
      reading.close();
      return;
    }
    else{
      reading >> command;
      if (command!="end") {
        if (debug) cout << "'end' missing or illegal extra data:" << command << "\n";
   	  Error=ImproperInputFormat;
	  reading.close();
	  return;
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

void ReadExtFile(ifstream &f)
{
  long i,j,k;
  myRational rvalue=0;
  myTYPE value=0;
  myRational mvalue=0;
  long mray,nray;
  string numbtype, command;
  boolean localdebug=False,found;
  myTYPE* vec;
  
  vec = new myTYPE[mm];
  AddArtificialRay();
  if (!f.is_open()) {
    Error=ImproperInputFormat;
    goto _L99;
  };
  found=False;
  while (!found)
  {
    if (f.eof()) {
     Error=ImproperInputFormat;
     goto _L99;
    }
    else {
      f >> command;
      if (command=="begin") {
        found=True;
      }
    }
  }
  f >> mray;
  f >> nray;
  f >> numbtype;
  if (DynamicWriteOn){ 
    cout << "ext object size =" << mray << "x" <<nray << "\n";
    cout << "Number Type = " << numbtype << "\n";
  }
  for (k=1; k<=mray;k++){
    for (j=1; j<=nray; j++){
        if (OutputNumberString=="real" && numbtype=="rational"){
          f >> rvalue;
          value=myTYPE(rvalue);
        } else {
          f >> value;
        }
      if (Inequality==NonzeroRHS) {
        vec[j - 1] = value;
      } else if (j>=2) {
        vec[j - 2] = value;
      }
      if (localdebug) WriteNumber(cout, value);
    }
    if (localdebug) printf("\n");
    if (DynamicRayWriteOn) cout << "#" << k << ":";
    AddRay(vec);
  }
_L99:;
}

/* end of cddio.c */

