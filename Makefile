# Makefile for cdd+ compiler.

# You must use GNU g++ compiler of version 2.6.0 or higher
CC = gcc

# Location of g++-library archive file libg++.a
LIBDIR = /usr/local/lib

# Location of g++-library include files
INCLUDEDIR = /usr/local/lib/g++-include
#INCLUDEDIR = /usr/local/include/g++

# Compiler options
#CFLAGS = -g -O -I$(INCLUDEDIR)
CFLAGS = -O -I$(INCLUDEDIR)
#CFLAGS = -g -I$(INCLUDEDIR)

########## You shouldn't have to change anything after this point ##########

LDFLAGS = -L$(LIBDIR)

cddio.o: cddio.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -c cddio.C

cddarith.o: cddarith.C cdd.h cdddef.h cddtype.h 
	$(CC) $(CFLAGS) -c cddarith.C

cddpivot.o: cddpivot.C cdd.h cdddef.h cddtype.h 
	$(CC) $(CFLAGS) -c cddpivot.C

cddrevs.o: cddrevs.C cddrevs.h cdd.h cdddef.h cddtype.h 
	$(CC) $(CFLAGS) -c cddrevs.C

cdd.o: cdd.C cdd.h cdddef.h cddtype.h cddrevs.h
	$(CC) $(CFLAGS) -c cdd.C

setoper.o: setoper.C
	$(CC) $(CFLAGS) -c setoper.C

cddio_r.o: cddio.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -c -o cddio_r.o -DRATIONAL cddio.C

cddarith_r.o: cddarith.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -c -o cddarith_r.o -DRATIONAL cddarith.C

cddpivot_r.o: cddpivot.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -c -o cddpivot_r.o -DRATIONAL cddpivot.C

cddrevs_r.o: cddrevs.C cdd.h cddrevs.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -c -o cddrevs_r.o -DRATIONAL cddrevs.C

cdd_r.o: cdd.C cdd.h cdddef.h cddtype.h cddrevs.h
	$(CC) $(CFLAGS) -c -o cdd_r.o -DRATIONAL cdd.C

cddr+: cdd_r.o cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o
	$(CC) $(CFLAGS) $(LDFLAGS) cdd_r.o cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o -o cddr+ -lg++

cddf+: cdd.o cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o
	$(CC) $(CFLAGS) $(LDFLAGS) cdd.o cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o -o cddf+ -lg++

clean:
	rm -rf core a.out cddf+ cddr+ *.o *~

all: cddr+ cddf+

