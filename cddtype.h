/* cddtype.h: Arithmetic type header file for cdd.C 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/


#ifdef	GMP
// using GMP wrapper of Polymake (polymake@math.tu-berlin.de)
#include <Rational.h>
#else
// using GNU g++ lib Rational library
#include <Rational.h>
#endif	// GMP

#ifdef	RATIONAL
#define ZERO 0
#define OUTPUTDIGITS 0
#else
#define ZERO 1.E-6
#define OUTPUTDIGITS 8
#endif	// RATIONAL

/* ZERO is the default value for the sign recognition.
   This should not be modified.  It can be controlled by
   float_zero option with caution.
   OUTPUTDIGITS is the number of decimal digits for each floating
   point output. It can be controlled by output_digits option with caution.
*/

/* end of cddtype.h */
