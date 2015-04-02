LAPACK=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE
CLAPACKLIB=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1
F2C=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS

CC=gcc
CFLAGS = -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include -I/usr/include/gsl -fpic  -g -O2

PROGRAM = rasqual
OBJS = src/usage.o src/parseVCF.o src/sort.o src/nbglm.o src/nbem.o src/util.o



.SUFFIXES: .c .o

.PHONY: all
all: $(PROGRAM)

$(PROGRAM): src/main.c $(OBJS)
	#sh src/makeUsage.sh
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include -I/usr/include/gsl -fpic  -g -O2 -c src/usage.c -o src/usage.o
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include -I/usr/include/gsl -fpic  -g -O2 -c src/parseVCF.c -o src/parseVCF.o
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include -I/usr/include/gsl -L$(F2C)  -fpic  -g -O2 -c src/sort.c -o src/sort.o
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include/gsl  -L$(F2C)       -fpic  -g -O2 -c src/nbglm.c -o src/nbglm.o
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include/gsl  -L$(F2C)     -fpic  -g -O2 -c src/nbem.c -o src/nbem.o
	#$(CC) -std=gnu99 -I$(F2C) -I$(LAPACK) -I/usr/include/gsl  -L$(F2C)      -fpic  -g -O2 -c src/util.c -o src/util.o
	$(CC) -L$(F2C) -I$(F2C) -I$(LAPACK) -o $(PROGRAM) $^ -I/usr/include/gsl -L$(CLAPACKLIB)  -lf2c -lblas -lgslcblas  -llapack -lgsl -lm -lz -ltmglib -lblas -lf2c

.c.o:
	$(CC) $(CFLAGS) -c $<

usage.c:
	sh src/makeUsage.sh

convert:
	/software/R-3.0.0/bin/R --vanilla --quiet --args test/Y.txt test/K.txt < R/txt2bin.R

clean:
	rm src/*.o src/usage.c
