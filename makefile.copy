# Makefile for cdd+ compilation with gcc-2.8.*.

# You must use GNU g++ compiler of version 2.8.0 or higher
# Use Makefile.2.7 when your compiler is older than gcc-2.8.0.
CC = gcc
#CC = /usr/local/bin/gcc
#CC = /bin/cc

# Location of gnu c++ library libstdc++.a.  
LIBDIR = /usr/local/lib 

# Location of gnu g++-library include files.  
INCLUDEDIR = /usr/local/include/g++

# Select the second line to use GMP instead of gnu g++ Rational
#GMPUSED = TRUE
GMPUSED = FALSE

# Location of gnu gmp library libgmp.a  
GMPLIBDIR = /usr/local/lib 

# Location of gnu gmp-library include file gmp.h 
GMPINCLUDEDIR = .
#GMPINCLUDEDIR = /usr/local/include

# Compiler optimization/debug options
#OPTFLAGS = -g -O
#OPTFLAGS = -g -pg -O
OPTFLAGS = -O3

########## You shouldn't have to change anything after this point ##########

ifeq ($(GMPUSED),TRUE)
RATLIB = gmp
GMPFLAG = -DGMP
RATOBJ = gmp_integer.o gmp_rational.o
RATEXE = cddr+_gmp
else
RATLIB = g++
GMPFLAG = 
RATOBJ = 
RATEXE = cddr+_g++
endif

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

gmp_integer.o: gmp_integer.cc gmp_integer.h
	$(CC) $(CFLAGS) -c -o gmp_integer.o gmp_integer.cc 

gmp_rational.o: gmp_rational.cc gmp_rational.h 
	$(CC) $(CFLAGS) -c -o gmp_rational.o gmp_rational.cc 

$(RATEXE): cdd_r.o $(RATOBJ) cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o gmp_rational.o
	$(CC) $(CFLAGS) $(LDFLAGS) cdd_r.o cddrevs_r.o cddio_r.o cddarith_r.o cddpivot_r.o setoper.o $(RATOBJ) -o $(RATEXE) $(LIBS)

cddf+: cdd.o $(RATOBJ) cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o
	$(CC) $(CFLAGS) $(LDFLAGS) cdd.o cddrevs.o cddio.o cddarith.o cddpivot.o setoper.o $(RATOBJ) -o cddf+ $(LIBS)

clean:
	rm -rf core a.out cddf+ cddr+* *.o *~

all: $(RATEXE) cddf+

