/* cdd.h: Header file for cdd.C 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define COPYRIGHT   "Copyright (C) 1999, Komei Fukuda, fukuda@ifor.math.ethz.ch"
#define DDVERSION   "Version 0.77(August 19, 2003)"

#ifdef RATIONAL
 #ifdef GMP
 #define ARITHMETIC  "Compiled for Rational Exact Arithmetic with GMP"
 #else
 #define ARITHMETIC  "Compiled for Rational Exact Arithmetic with G++"
 #endif
#else
#define ARITHMETIC  "Compiled for Floating-Point Arithmetic"
#endif

#include <ctime>
#include <iostream>
#include <fstream>

using std::string;
using std::ifstream;
using std::ofstream;


#ifdef  RATIONAL
typedef Rational myTYPE;
#else
typedef double myTYPE;
#endif  // RATIONAL
typedef Rational myRational;
typedef int boolean;
typedef long rowrange;
typedef long colrange;
typedef set_type rowset;  /* set_type defined in setoper.h */
typedef set_type colset;
typedef long *rowindex;   
    /* rowindex should be intialized to be an array of [mm+1] components */
typedef long colindex[NMAX+1];
typedef myTYPE *Amatrix[MMAX];
typedef myTYPE *Arow;
typedef myTYPE *Bmatrix[NMAX];
typedef set_type Aincidence[MMAX];
typedef int *SignAmatrix[MMAX]; /* Sign (+1, 0 ,-1)-matrix of Amatrix type */
typedef char DataFileType[filenamelen];
typedef char LineType[linelenmax];
typedef char WordType[wordlenmax];

typedef struct RayRecord {
  myTYPE *Ray;
  rowset ZeroSet;
  rowrange FirstInfeasIndex;  /* the first inequality the ray violates */
  boolean feasible;  /* flag to store the feasibility */
  myTYPE *ARay;   /* temporary area to store some row of A*Ray */
  struct RayRecord *Next;
} RayRecord;

typedef struct AdjacencyRecord {
  RayRecord *Ray1, *Ray2;
  struct AdjacencyRecord *Next;
} AdjacencyRecord;

typedef struct node {long key; struct node *next;} node;

typedef enum {
  Combinatorial, Algebraic
} AdjacencyTestType;

typedef enum {
  MaxIndex, MinIndex, MinCutoff, MaxCutoff, MixCutoff,
   LexMin, LexMax, RandomRow, LineShelling
} HyperplaneOrderType;

typedef enum {
  Real, Rational, Integer, Unknown
} NumberType;

typedef enum {
  ZeroRHS, NonzeroRHS
} InequalityType;

typedef enum {
  IneToExt, ExtToIne, Projection, 
  LPmax, LPmin, FacetListing, VertexListing, TopeListing, InteriorFind,
  FacetListingExternal,VertexListingExternal
} ConversionType;

typedef enum {
  RowDecomposition, RowSubproblemSolve, Nothing
} SpecialConversionType;

typedef enum {
  CrissCross,DualSimplex,CombMaxImprove
} LPsolverType;

typedef enum {
  IncOff=0, IncCardinality, InputIncidence, OutputIncidence, IOIncidence
} IncidenceOutputType;

typedef enum {
  AdjOff=0, OutputAdjacency, InputAdjacency, IOAdjacency
} AdjacencyOutputType;

typedef enum {
  Auto, SemiAuto, Manual
} FileInputModeType;   
   /* Auto if a input filename is specified by command arguments */

typedef enum {
  DimensionTooLarge, LowColumnRank, ImproperInputFormat, DependentMarkedSet, 
  ImproperExecutable, FileNotFound, None
} ErrorType;

typedef enum {
  InProgress, AllFound, RegionEmpty
} CompStatusType;

typedef enum {
  LPSundecided, Optimal, Inconsistent, DualInconsistent, Unbounded, DualUnbounded
} LPStatusType;

extern long minput, ninput;   /*size of input data [b -A] */
extern long mm, nn;   /*size of the homogenous system to be solved by dd*/
extern long projdim;  /*dimension of orthogonal preprojection */
extern colset projvars;   /*set of variables spanning the space of preprojection, 
     i.e. the remaining variables are to be removed*/
extern rowset EqualitySet, NonequalitySet, GroundSet, Face, Face1,SubproblemRowSet, RedundantRowSet;
extern rowrange Iteration, hh;
extern rowindex OrderVector;
extern rowindex EqualityIndex;  
extern rowset AddedHyperplanes, WeaklyAddedHyperplanes, InitialHyperplanes;
extern long RayCount, FeasibleRayCount, WeaklyFeasibleRayCount,
 TotalRayCount, VertexCount,ZeroRayCount;
extern long EdgeCount,TotalEdgeCount;
extern long count_int,count_int_good,count_int_bad;
extern boolean DynamicWriteOn, DynamicRayWriteOn, LogWriteOn, 
  ShowSignTableauOn, PostAnalysisOn, 
  CondensedListOn, SignPivotOn, ManualPivotOn, debug;
extern Amatrix AA;
extern Bmatrix InitialRays;
extern colindex InitialRayIndex;
extern Arow LPcost;
extern LPStatusType LPStatus;
extern colrange RHScol;   /* LP RHS column */
extern rowrange OBJrow;   /* LP OBJ row */
extern RayRecord *ArtificialRay, *FirstRay, *LastRay;
extern RayRecord *PosHead, *ZeroHead, *NegHead, *PosLast, *ZeroLast, *NegLast;
extern AdjacencyRecord *Edges[MMAX];  /* adjacency relation storage for iteration k */
extern boolean RecomputeRowOrder, inputsuccessful;
extern HyperplaneOrderType HyperplaneOrder;
extern AdjacencyTestType AdjacencyTest;
extern NumberType Number;
extern string InputNumberString, OutputNumberString;
extern InequalityType Inequality;
extern boolean NondegAssumed;   /* Nondegeneacy preknowledge flag */
extern boolean InitBasisAtBottom;  /* if it is on, the initial Basis will be selected at bottom */
extern boolean RestrictedEnumeration; /* Restricted enumeration Switch (True if it is restricted on the intersection of EqualitySet hyperplanes) */
extern boolean RelaxedEnumeration; /* Relaxed enumeration Switch (True if NonequalitySet inequalities must be satisfied with strict inequality) */
extern boolean VerifyInput; /* Verification switch for the input data */
extern boolean PreOrderedRun; 
extern CompStatusType CompStatus;     /* Computation Status */
extern ConversionType Conversion;
extern SpecialConversionType SpecialConversion;  /* rowdecomposition, rowsubproblem */
extern LPsolverType LPsolver;
extern IncidenceOutputType IncidenceOutput;
extern AdjacencyOutputType AdjacencyOutput;
extern ErrorType Error;
extern FileInputModeType FileInputMode;
extern DataFileType inputfile,ifilehead,ifiletail,
     outputfile,projfile,incfile,adjfile,logfile,dexfile,verfile,xtnfile;
extern time_t starttime, endtime;
extern unsigned int rseed;

extern myTYPE zero;    /*Rational or floating zero*/
extern boolean Round_Output;    /* rounding option for floating-point output. */
extern int output_digits;  /* Float digits for output.  Does not affect the computation. */

void SetInputFile(boolean *);
void SetWriteFileName(DataFileType, char, char *);
void SetReadFileName(DataFileType, char, char *);

myTYPE FABS(myTYPE);
void SetNumberType(string);
void ProcessCommandLine(ifstream &, string);
void AmatrixInput(boolean *);
void SetInequalitySets(rowindex);
void RandomPermutation(rowindex, long, unsigned int);
void QuickSort(rowindex, long, long, Amatrix, long);
myTYPE AValue(myTYPE *, rowrange );
void WriteIncidence(ostream &, RayRecord *);
void StoreRay1(myTYPE *, RayRecord *, boolean *);
void StoreRay2(myTYPE *, RayRecord *, boolean *, boolean *);
void AddRay(myTYPE *);
void AddArtificialRay(void);
void ConditionalAddEdge(RayRecord *Ray1, RayRecord *Ray2, RayRecord *ValidFirstRay);
void CreateInitialEdges(void);
void UpdateEdges(RayRecord *RRbegin, RayRecord *RRend);
void FreeDDMemory(void);
void Normalize(myTYPE *);
void ZeroIndexSet(myTYPE *, rowset);
void CopyBmatrix(Bmatrix T, Bmatrix TCOPY);
void SelectPivot1(Amatrix, HyperplaneOrderType,
   rowrange, rowset, colset, rowrange *, colrange *,boolean *);
myTYPE TableauEntry(Amatrix, Bmatrix T, rowrange, colrange);
void WriteTableau(ostream &,Amatrix, Bmatrix T, InequalityType);
char Sign(myTYPE);
int myTYPE2sign(myTYPE);
int long2sign(long);
void WriteSignAmatrix(ostream &f, SignAmatrix X,
  rowindex OV, long bflag[], rowrange objrow, colrange rhscol);
void WriteSignTableau(ostream &f, Amatrix X, Bmatrix T,
  rowindex OV, long bflag[], rowrange objrow, colrange rhscol);
void WriteDictionary(ostream &f, Amatrix X, Bmatrix T,
  rowindex OV, long bflag[], colindex NBIndex, rowrange objrow, colrange rhscol);
void OutputTableau(Amatrix, Bmatrix T,InequalityType);
void SelectPivot2(Amatrix, Bmatrix T,
   HyperplaneOrderType,rowrange, rowset, colset,
   rowrange *, colrange *,boolean *);
void GausianColumnPivot1(Amatrix, Bmatrix, SignAmatrix, rowrange, colrange);
void GausianColumnPivot2(Amatrix, Bmatrix,  rowrange, colrange);
void InitializeBmatrix(Bmatrix T);
void SetToIdentity(Bmatrix T);
void free_Bmatrix(Bmatrix T);
void WriteBmatrix(ostream &, Bmatrix T);
void ReduceAA(rowset, colset);
void DualizeAA(Bmatrix T);
void EnlargeAAforInteriorFinding(void);
void EnlargeAAforZeroRHSLP(void);
void RecoverAAafterInteriorFinding(void);
void ShiftPointsAroundOrigin(ostream &, ostream &, Arow);
boolean RowEquivalent_Q(Arow a1, Arow a2, colrange n);
void FindRowEquivalenceClasses(long *classno, rowindex);
void WriteSubMatrixOfAA(ostream &, rowset, colset, InequalityType);
void WriteAmatrix(ostream &, Amatrix, long, long, InequalityType);
void WriteErrorMessages(ostream &);
void ComputeRank(Amatrix, unsigned long *, long *);
void ComputeBInverse(Amatrix, long, Bmatrix InvA1, long *);
void FindBasis(Amatrix, HyperplaneOrderType, rowset, long *,
   Bmatrix BasisInverse, long *);
void SelectCrissCrossPivot(Amatrix, Bmatrix, rowindex,
  long, rowrange,colrange,rowrange *,colrange *,
  boolean *, LPStatusType *);
void SelectDualSimplexPivot(boolean, Amatrix, Bmatrix, rowindex,
    colindex, long, rowrange, colrange,
    rowrange *, colrange *, boolean *, LPStatusType *);
void CrissCrossMinimize(ostream &, ostream &, Amatrix,Bmatrix,
  rowrange, colrange, boolean,
  LPStatusType *, myTYPE *, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void CrissCrossMaximize(ostream &, ostream &, Amatrix,Bmatrix,
  rowrange, colrange, boolean, 
  LPStatusType *, myTYPE *, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void DualSimplexMinimize(ostream &, ostream &, Amatrix,Bmatrix,
  rowrange, colrange, boolean,
  LPStatusType *, myTYPE *, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void DualSimplexMaximize(ostream &, ostream &, Amatrix,Bmatrix,
  rowrange, colrange, boolean,
  LPStatusType *, myTYPE *, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void ManualPivot(ostream &, ostream &, Amatrix,Bmatrix BasisInverse,
  rowrange, colrange, 
  LPStatusType *, myTYPE *optvalue, Arow, Arow, colindex,
  rowrange *, colrange *, long *);
void WriteLPResult(ostream &, LPStatusType, myTYPE, 
  Arow, Arow, colindex, rowrange, colrange, long);
void FindInitialRays(rowset, Bmatrix, colindex, boolean *);
void CheckAdjacency1(RayRecord **, RayRecord **,boolean *);
void CheckAdjacency2(RayRecord **, RayRecord **,boolean *);
void CheckEquality(RayRecord **, RayRecord **, boolean *);
void Eliminate(RayRecord **);
void CreateNewRay(RayRecord *, RayRecord *, rowrange);
void EvaluateARay1(rowrange);
void EvaluateARay2(rowrange);
void FeasibilityIndices(long *, long *, rowrange);
boolean LexSmaller(myTYPE *, myTYPE *, long);
boolean LexLarger(myTYPE *, myTYPE *, long);
void LineShellingOrder(rowindex, myTYPE *, myTYPE *);
void CompileDecompResult(ofstream &);
void ReadExtFile(ifstream &);
void CopyArow(myTYPE *, myTYPE *, long);
void ComputeRowOrderVector(rowindex OV, HyperplaneOrderType ho);
void UpdateRowOrderVector(rowset PriorityRows);
void SelectNextHyperplane(HyperplaneOrderType, unsigned long *, rowrange *, boolean *);
void SelectPreorderedNext(long *excluded, rowindex, rowrange *hnext);
void AddNewHyperplane1(rowrange);
void AddNewHyperplane2(rowrange);
void WriteRunningMode(ostream &);
void WriteRunningMode0(ostream &);
void WriteRunningMode2(ostream &);
void WriteCompletionStatus(ostream &);
void WriteTimes(ostream &);
void WriteProgramDescription(ostream &);
void WriteSolvedProblem(ostream &);
void WriteNumber(ostream &, myTYPE);
void WriteRayRecord(ostream &, RayRecord *);
void WriteRayRecord2(ostream &, RayRecord *);
void WriteExtFile(ostream &, ostream &);
void WriteIncidenceFile(ostream &);
void WriteInputIncidenceFile(ostream &);
void WriteInputAdjacencyFile(ostream &);
void WriteAdjacencyFile(ostream &);
void WriteProjResult(ostream &, ostream &, long *);
void WriteProjRayRecord(ostream &, RayRecord *, long *);
void InitialWriting(ostream &, ostream &);
void WriteDecompResult(ostream &f, ostream &f_log);
void WriteRowEquivalence(ostream &f, long classno, rowindex rowequiv);
void OutputHeading(void);
void InitialDataSetup(void);
void LPInit(void);
void DDInit(void);
void LPMain(ostream &, ostream &);
void PostAnalysisMain(ifstream &, ostream &);
void InteriorFindMain(ostream &, ostream &, boolean *);
void WriteCurrentSolution(ostream &, Amatrix, Bmatrix, rowrange, colrange, colindex);
void CheckConversionConsistency(ostream &, ostream &);

// procedures in cddrevs.C
boolean Facet_Q(topeOBJECT, rowrange);
boolean Facet_Q2(topeOBJECT, rowrange,colindex,Arow);
void ReverseSearch(ostream &, topeOBJECT, long delta);
void FacetandVertexListMain(ostream &, ostream &);
void FacetandVertexExternalListMain(ostream &, ostream &);
void TopeListMain(ostream &, ostream &);

/* end of cdd.h */


