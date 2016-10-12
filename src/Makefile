CLAPACKINC=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1/INCLUDE
CLAPACKLIB=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1
F2CLIBS=/nfs/users/nfs_n/nk5/Applications/CLAPACK-3.1.1.1/F2CLIBS
GSLLIBS=/nfs/users/nfs_n/nk5/Applications/gsl-1.16/lib
INTELLIB=/software/intel-tools-2015/lib/intel64

CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -I/usr/include -I/usr/include/gsl -fpic -g -O2
LDFLAGS := $(LDFLAGS) -lf2c -lblas -lgslcblas -llapack -lgsl -lm -lz -lpthread

PROGRAM = rasqual
OBJS = usage.o parseVCF.o sort.o nbglm.o nbem.o util.o

.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): main.c $(OBJS)
	$(CC) $(CFLAGS) -o $(PROGRAM) $^ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

usage.c:
	sh makeUsage.sh

install:
	mkdir -p ../bin
	mv rasqual ../bin

clean:
	rm *.o usage.c
	rm -rf ../bin
