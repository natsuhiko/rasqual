
LAPACK=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE
F2C=/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS

grep -v "^#" usage.txt | sed 's/\\$/\\\\/g' | awk -F '\n' 'BEGIN{print "#include<stdio.h>\nvoid usage(){"}{print "\tfprintf(stderr, \""$1"\\n\");"};END{print "}"}' > usage.c

gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include -I/usr/include/gsl -fpic  -g -O2 -c usage.c -o usage.o
gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include -I/usr/include/gsl -fpic  -g -O2 -c parseVCF.c -o parseVCF.o
gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include -I/usr/include/gsl -I$F2C -L$F2C  -fpic  -g -O2 -c sort.c -o sort.o 
gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include/gsl  -I$F2C -L$F2C       -fpic  -g -O2 -c nbglm.c -o nbglm.o 
gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include/gsl    -I$F2C -L$F2C     -fpic  -g -O2 -c nbem.c -o nbem.o 
gcc -std=gnu99 -I$F2C -I$LAPACK -I/usr/include/gsl  -I$F2C -L$F2C      -fpic  -g -O2 -c util.c -o util.o 

gcc -L$F2C -I$F2C -I$LAPACK -o rasqual usage.o nbem.o nbglm.o sort.o parseVCF.o util.o main.c -I$F2C -I/usr/include/gsl -L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1  -lf2c -lblas -lgslcblas  -llapack -lgsl -lm -lz -ltmglib -lblas -lf2c

rm nbem.o nbglm.o sort.o parseVCF.o util.o usage.c usage.o

