/* cdddef.h:  Definition file for cdd.C 
   written by Komei Fukuda, fukuda@ifor.math.ethz.ch
   Version 0.77, August 19, 2003 
*/

/* cdd.C : C-Implementation of the double description method for
   computing all vertices and extreme rays of the polyhedron 
   P= {x :  b - A x >= 0}.  
   Please read COPYING (GNU General Public Licence) and
   the manual cddman.tex for detail.
*/

#define MMAX      50002  /* USER'S CHOICE: max row size of A plus two */
#define NMAX      101   /* USER'S CHOICE: max column size of A plus one */

#define rowsetsize MMAX   /* The size of the column index set */
#define colsetsize NMAX   /* The size of the row index set */

#define datawidth       10
#define filenamelen     256 
#define wordlenmax      128 
#define linelenmax      256

#define False 0
#define True 1


/* end of cdddef.h */

