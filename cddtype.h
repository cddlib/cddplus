/* cddtype.h: Arithmetic type header file for cdd.C 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.75, November 30, 1997 
*/

/* cdd.C : C++-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#include <Rational.h>

#ifdef	RATIONAL
#define ZERO 0
#define OUTPUTDIGITS 0
#else
#define ZERO 1.E-6
#define OUTPUTDIGITS 8
#endif	RATIONAL

/* ZERO is the default value for the sign recognition.
   This should not be modified.  It can be controlled by
   float_zero option with caution.
   OUTPUTDIGITS is the number of decimal digits for each floating
   point output. It can be controlled by output_digits option with caution.
*/

/* end of cddtype.h */
