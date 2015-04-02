#/bin/bash

# converting read count and offset files into binary data
${RHOME}/R --vanilla --quiet --args data/Y.txt data/K.txt < R/txt2bin.R > log

# mapping eQTLs for C11orf21 and TSPAN11 genes

tabix test/chr11.gz 11 | ./rasqual -y test/Y.bin -k test/K.bin -n 24 -l 1000 -m 1000 -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 -z -j 1 -t
