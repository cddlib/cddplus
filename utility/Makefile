# Makefile for utitlity programs of cdd-061.

# Select ANSI C compiler
#CC = /usr/local/bin/gcc
CC = /usr/bin/gcc

# Location of include files
INCLUDEDIR = .

# Compiler options
#CFLAGS = -g -O -I$(INCLUDEDIR)
#CFLAGS = -g -pg -O -I$(INCLUDEDIR)
CFLAGS = -O3 -I$(INCLUDEDIR)

########## You shouldn't have to change anything after this point ##########

rlp: rlp.c rlp.h
	$(CC) $(CFLAGS) -o rlp rlp.c

voronoi: voronoi.c
	$(CC) $(CFLAGS) -o voronoi voronoi.c

delaunay: delaunay.c
	$(CC) $(CFLAGS) -o delaunay delaunay.c


clean:
	rm -rf core a.out voronoi delaunay rlp *.o *~

all: voronoi delaunay rlp

