# Makefile for cdd+ compilation with gcc-3.*.

# You must have GMP installed and use GNU g++ compiler of version 3.* .
CC = gcc
#CC = /usr/local/bin/gcc
#CC = /bin/cc

# Location of gnu c++ library.  
LIBDIR = /usr/lib 
#LIBDIR = /usr/local/lib 

# Location of gnu gmp library libgmp.a  
GMPLIBDIR = /usr/lib 
#GMPLIBDIR = /usr/local/lib 

# Location of gnu gmp-library include file gmp.h 
GMPINCLUDEDIR = /usr/include
#GMPINCLUDEDIR = /usr/local/include

# Compiler optimization/debug options
#OPTFLAGS = -g -static -O
#OPTFLAGS = -g -static -pg -O
OPTFLAGS = -static -O3

########## You shouldn't have to change anything after this point ##########

RATLIB = gmp
GMPFLAG = -DGMP
RATOBJ = gmp_init.o Integer.o Rational.o
RATEXE = cddr+

CFLAGS = $(OPTFLAGS) -I$(INCLUDEDIR) -I$(GMPINCLUDEDIR) -I. $(GMPFLAG)

LDFLAGS = -L$(LIBDIR) -L$(GMPLIBDIR)

LIBS = -lstdc++ -l$(RATLIB)

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
	$(CC) $(CFLAGS) -DRATIONAL -c -o cddio_r.o cddio.C

cddarith_r.o: cddarith.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -DRATIONAL -c -o cddarith_r.o cddarith.C

cddpivot_r.o: cddpivot.C cdd.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -DRATIONAL -c -o cddpivot_r.o cddpivot.C

cddrevs_r.o: cddrevs.C cdd.h cddrevs.h cdddef.h cddtype.h
	$(CC) $(CFLAGS) -DRATIONAL -c -o cddrevs_r.o cddrevs.C

cdd_r.o: cdd.C cdd.h cdddef.h cddtype.h cddrevs.h
	$(CC) $(CFLAGS) -DRATIONAL -c -o cdd_r.o cdd.C

gmp_init.o: gmp_init.cc
	$(CC) $(CFLAGS) -c -o gmp_init.o gmp_init.cc 

Integer.o: Integer.cc
	$(CC) $(CFLAGS) -c -o Integer.o Integer.cc 

Rational.o: Rational.cc 
	$(CC) $(CFLAGS) -c -o Rational.o Rational.cc 

$(RATEXE): cdd_r.o $(RATOBJ) cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o 
	$(CC) $(CFLAGS) $(LDFLAGS) cdd_r.o cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o $(RATOBJ) -o $(RATEXE) $(LIBS)

cddf+: cdd.o $(RATOBJ) cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o
	$(CC) $(CFLAGS) $(LDFLAGS) cdd.o cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o $(RATOBJ) -o cddf+ $(LIBS)

clean:
	rm -rf core a.out cddf+ cddr+* *.o *~

all: $(RATEXE) cddf+

