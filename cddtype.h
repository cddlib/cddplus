/* cddtype.h: Arithmetic type header file for cdd.C 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.73, September 6, 1995 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/


#ifdef	RATIONAL
#include <Rational.h>
#define ZERO 0
#else
#define	ZERO 1.0E-5
#endif	RATIONAL

/* end of cddtype.h */
