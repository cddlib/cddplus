/* setoper.C:
 A set operation library 
 created by Komei Fukuda, April 3, 1995
 Last modified, June 19, 1996
 Modified for newer C++ standard, October 21, 2007 
*/

#include <iostream>
#include <climits>

using std::cout;

extern "C" {
#include "setoper.h"
}

#define SETBITS (sizeof(long) * CHAR_BIT)
/* (Number of chars in a long) * (number of bits in a char) */

/* Definitions for optimized set_card function 
   by David Bremner bremner@cs.mcgill.ca  
*/
/* Caution!!!
   Bremner's technique depends on the assumption that CHAR_BIT == 8.
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

long set_blocks(unsigned long len)
{
	long blocks=1L;
	
	if (len>0) blocks=((long)len-1)/SETBITS+2;
	// cout << "length =" << len << "   blocks =" << blocks <<"\n";
	return blocks;
}

void set_initialize(set_type *setp, long len)
/* Make a set with a given bit lengths  */
{
	long i,forlim1;
	
	forlim1=set_blocks(len);
	(*setp)=new unsigned long[forlim1];
	(*setp)[0]=(unsigned long) len;  /* size of the ground set */
	// cout << " size of the set is set to " << (*setp)[0] << "\n";
	for (i=1; i<forlim1; i++)
		(*setp)[i]=0U;
}

void set_free(set_type *set)
/* Free the space created with the set pointer set*/
{
    delete[] (*set);
}

void set_emptyset(set_type set)
/* Set set to be the emptyset  */
{
	long i,forlim;
	
	forlim=set_blocks(set[0])-1;
	for (i=1; i<=forlim; i++)
		set[i]=0;
}


void set_copy(set_type setcopy,set_type set)
/* Copy the set set[] to setcopy[] with setcopy[] length */
{
	long i,forlim;

	forlim=set_blocks(setcopy[0])-1;
	for (i=1; i<=forlim; i++)
		setcopy[i]=set[i];
}

void set_addelem(set_type set, long elem)
/* add elem only if it is within the set[] range */
{
	long i,j;
	unsigned long change;
	unsigned long one=1U;
	
	if (elem<=set[0])    
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change= one << j;  /* put 1 in jth position */
		set[i]=set[i] | change;
	}
}

void set_delelem(set_type set, long elem)
/* delete elem only if it is within the set[] range */
{
	long  i,j;
	unsigned long change;
	unsigned long one=1U;	 
	
	if (elem<=set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		change=one << j;  /* put 1 in jth position */
		set[i]=(set[i] | change) ^ change;
	}
}

void set_int(set_type set,set_type set1,set_type set2)
/* Set intersection, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;
	
	forlim=set_blocks(set[0])-1;
	for (i=1;i<=forlim;i++)
		set[i]=(set1[i] & set2[i]);
}

void set_uni(set_type set,set_type set1,set_type set2)
/* Set union,assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] | set2[i];
}

void set_diff(set_type set,set_type set1,set_type set2)
/* Set difference se1/set2, assuming set1 and set2 have the same length as set */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]=set1[i] & (~set2[i]);
}

void set_compl(set_type set,set_type set1)
/* set[] will be set to the complement of set1[] */
{
	long  i,forlim;

	forlim=set_blocks(set[0])-1;	
	for (i=1;i<=forlim;i++)
		set[i]= ~set1[i];
}

int set_subset(set_type set1,set_type set2)
/* Set containment check, set1 <= set2 */
{
	int  yes=1;
	long i,forlim;
	
	forlim=set_blocks(set2[0])-1;
	for (i=1;i<=forlim && yes;i++)
		if ((set1[i] | set2[i])!=set2[i])
			yes=0;
	return yes;
}

int set_member(long elem, set_type set)
/* Set membership check, elem in set */
{
	int  yes=0;
	long  i,j;
	unsigned long testset;
	unsigned long one=1U;	 
	
	if (elem<=set[0])
	{
		i=(elem-1)/SETBITS+1;
		j=(elem-1)%SETBITS;
		testset=set[i] | (one<<j);   /* add elem to set[i] */
		if (testset==set[i])
			yes=1;
	}
	return yes;
}

long set_card(set_type set)
/* set cardinality, modified by David Bremner bremner@cs.mcgill.ca
   to optimize for speed. */
{
  unsigned long block;
  long car=0;
  set_card_lut_t *p;
  
  p=(set_card_lut_t *)&set[1];
  for (block=0; block< LUTBLOCKS(set);block++) {
    car+=set_card_lut[p[block]];
  }
  return car;
}

/* old cardinality code 
long set_card(set_type set)
{
	long elem,car=0;
	
	for (elem=1; elem<=set[0]; elem++) {
		if (set_member(elem,set)) car++;
    }
	return car;
}
*/ 

void set_write(set_type set)
{
	long elem;
	
	for (elem=1;elem<=set[0];elem++)
	{
		if (set_member(elem,set))
			cout << " " << elem;
	}
}

void set_fwrite(ostream &f,set_type set)
{
	long elem;
	
	for (elem=1;elem<=set[0];elem++)
	{
		if (set_member(elem,set))
			f << " " << elem;
	}
}

void set_binwrite(set_type set)
{
	long i,j;
	long forlim;
	unsigned long e1,e2;
	
	cout << "ground set size = " << set[0] << "\n";
	forlim=set_blocks(set[0])-1;
	for (i=forlim;i>=1;i--)
	{
		e1=e2=set[i];
		for (j=SETBITS-1;j>=0;j--)
		{
			e1=(e1>>j);
                        cout.width(1);
			cout << e1;
			e1=e2-(e1<<j);
			e2=e1;
		}
		cout << " ";
	}
}


void set_fbinwrite(ostream &f,set_type set)
{
	long i,j;
	long forlim;
	unsigned long e1,e2;
	
	f << "ground set size = " << set[0] << "\n";
	forlim=set_blocks(set[0])-1;
	for (i=forlim;i>=1;i--)
	{
		e1=e2=set[i];
		for (j=SETBITS-1;j>=0;j--)
		{
			e1=(e1>>j);
                        f.width(1);
			f << e1;
			e1=e2-(e1<<j);
			e2=e1;
		}
		f << " ";
	}
}

/* End of the library:  setoper.C  */

