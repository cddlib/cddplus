/* header file for cdd.C and cddrevs.C
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron
   P= {x :  b - A x >= 0}.
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

class topeOBJECT {
  long dim;
  int * sv;  /* {+1,-1}^dim sign vector */

public:
  topeOBJECT(long d) {long j; dim=d; sv=new int[d];for (j=1; j<=d; j++) sv[j-1]=1;}
  ~topeOBJECT() {delete[] sv;}
  topeOBJECT(const topeOBJECT& tope);
  void operator=(const topeOBJECT&);
  int operator[](long i); // return the i-th component (saved as i-1 comp) of tope 
  friend topeOBJECT operator-(topeOBJECT, long); // reversing the sign of sv[j-1]
  void fwrite(ostream&);
  friend int operator==(const topeOBJECT &t1, const topeOBJECT &t2);
  friend int operator!=(const topeOBJECT &t1, const topeOBJECT &t2);
  friend topeOBJECT f(topeOBJECT);
  friend long NeighbourIndex(topeOBJECT, topeOBJECT);
  friend topeOBJECT Adj(topeOBJECT, long);
};
/* {+1,-1}^dim vector representing a full-dimentional region of
the arrangement of hyperplanes associated with Ax <= b */

// end of cddrevs.h

