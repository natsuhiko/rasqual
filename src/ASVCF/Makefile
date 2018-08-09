CLAPACKINC=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE
CLAPACKLIB=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1
F2CLIBS=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS

CC=gcc
CFLAGS := $(CFLAGS) -std=gnu99 -I/usr/include -fpic -g -O2
#LDFLAGS := $(LDFLAGS)  -lf2c -lblas -lgslcblas  -llapack -lgsl -lm -lz -ltmglib -lblas -lf2c
LDFLAGS := $(LDFLAGS)  -lm -lz

#PROGRAM = rasqual.test2
#OBJS = usage.o parseVCF.o sort.o nbglm.o nbem.o util.o



.SUFFIXES: .c .o

.PHONY: all
all: qcFilterBam.o pasteFiles.o zpaste.o countAS.o
	$(CC) -o qcFilterBam qcFilterBam.o $(LDFLAGS)
	$(CC) -o pasteFiles pasteFiles.o $(LDFLAGS)
	$(CC) -o zpaste zpaste.o $(LDFLAGS)
	$(CC) -o countAS countAS.o $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) -c $<

clean:
	rm qcFilterBam pasteFiles zpaste countAS *.o
