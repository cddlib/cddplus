/* Header file for setoper.C  */

/* setoper.C: 
 * A set operation library 
 * created by Komei Fukuda, April 3, 1995
 * last modified, April 15, 1995
 */

#include <fstream.h>

#define SETBITS 32       /* Important Constant: Number of bits in a long integer */
#include <stdio.h>
#include <stdlib.h>

typedef unsigned long *set_type;  /* set type definition */

/* Definitions for optimized set_card function 
   by David Bremner bremner@cs.mcgill.ca  
 */

typedef unsigned char set_card_lut_t;

#define LUTBLOCKS(set) (((set[0]-1)/SETBITS+1)*(sizeof(long)/sizeof(set_card_lut_t)))

static unsigned char set_card_lut[]={
0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};
/* End of Definitions for optimized set_card */

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
void set_fwrite(ofstream, set_type set);
void set_binwrite(set_type set);
void set_fbinwrite(ofstream,set_type set);

/* End of File: setoper.h */
/* End of File: setoper.h */

