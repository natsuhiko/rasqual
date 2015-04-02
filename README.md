# rasqual development
Robust Allele Specific Quantification and quality controL

## How to build & install

	RASQUALDIR=/path/to/rasqualdir/
	cd $RASQUALDIR/src
	export CFLAGS="-I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE -I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	export LDFLAGS="-L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1 -L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	make
	make install

## Getting started

	cd $RASQUALDIR
	
