/* Header file for setoper.C  */

/* setoper.h: 
 * A set operation library 
 * created by Komei Fukuda, April 3, 1995
 * Last modified, June 19, 1996
 * Modified for newer C++ standard, April 20, 2003  
 */

#include <iostream>

using std::ostream;

typedef unsigned long *set_type;  /* set type definition */

long set_blocks(unsigned long len);
void set_initialize(set_type *setp,long len);
void set_free(set_type *set);
void set_emptyset(set_type set);
void set_copy(set_type setcopy,set_type set);
void set_addelem(set_type set, long elem);
void set_delelem(set_type set, long elem);
void set_int(set_type set,set_type set1,set_type set2);
void set_uni(set_type set,set_type set1,set_type set2);
void set_diff(set_type set,set_type set1,set_type set2);
void set_compl(set_type set,set_type set1);
int set_subset(set_type set1,set_type set2);
int set_member(long elem, set_type set);
long set_card(set_type set);
void set_write(set_type set);
void set_fwrite(ostream &, set_type set);
void set_binwrite(set_type set);
void set_fbinwrite(ostream &,set_type set);

/* End of File: setoper.h */
/* End of File: setoper.h */

