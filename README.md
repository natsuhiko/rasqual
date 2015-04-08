# RASQUAL
RASQUAL (Robust Allele Specific Quantification and quality controL) maps QTLs using population and allele-specific signals for any sequenced based cellular tratis.

## How to build & install

Please make sure CLAPACK is installed in your environment.  Otherwise you can get the library at http://www.netlib.org/clapack/.  To build and install RASQUAL, firstly go to the _source_ directory (src), then set environment variables appropriately to point the CLAPACK library as follows.  Finally use "make" to build and install RASQUAL which will be installed in "$RASQUALDIR/bin".

	RASQUALDIR=/path/to/rasqualdir/
	cd $RASQUALDIR/src
	export CFLAGS="-I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/INCLUDE -I/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	export LDFLAGS="-L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1 -L/nfs/users/nfs_n/nk5/team170/Natsuhiko/CLAPACK-3.1.1.1/F2CLIBS"
	make
	make install

## First QTL mapping with RASQUAL

Mapping QTLs with RASQUAL is straightforward.  You just need to prepare two data files:

1. fragment (read) count table with sample specific offset (usually library size)
2. VCF file with imputed SNP genotypes and allele-specific counts

For example, the following commands give you the first expression QTL mapping for two genes (C11orf21 and TSPAN32) with RASQUAL:

    # make sure tabix is installed in your environment
    tabix data/chr11.gz 11 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 409 -m 63 \
        -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 
        -t -f C11orf21 -z
    tabix data/chr11.gz 11 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 2 -l 409 -m 61 \
        -s 2323227,2323938,2324640,2325337,2328175,2329966,2330551,2331219,2334884,2335715,2338574,2339093 \
        -e 2323452,2324188,2324711,2325434,2328220,2330040,2330740,2331248,2334985,2337897,2338755,2339430 \
        -t -f TSPAN32 -z

Here RASQUAL requies SNP information in the VCF file as standard input and the count table and sample specific offset as binary files (Y.bin and K.bin, respectively).  Sample size (*N*=24) is given by **-n** option and the feature ID is given by **-j** option (only two genes exist in this example, thereby j=1,2).  To save memory usage, you are also required to count the numbers of testing SNPs and feature SNPs a priori (**-l** and **-m**, respectively).  RASQUAL also requires feature start and end positions (as comma separated values for more than one positions, i.e., union of exons in this example) for inputs (**-s** and **-e**, respectively).  As a default, RASQUAL outputs QTL mappint results for all tested SNPs, but you can specify to get only the lead QTL SNP (**-t** option).  In the output, you can also specify the feature name by **-f** option.  To take account of genotype uncertainty, imputation quality score (R square value) is used in this example (**-z** option; see the section below).  

## Data preparation

You can find an example expression data for C11orf21 and TSPAN32 genes in the _data_ directory.  There are two files: read count table (Y.txt) and sample specific offset data (K.txt), both has to be organized feature by sample format (i.e., row: gene; col: sample).  The following R script allows you to convert the count and offset files (text) into binary format used in RASQUAL software.  The script essentially converts a table data into a vector of values (double precision).

	cd $RASQUALDIR
	export RHOME=/software/R-3.0.0/bin/
	$RHOME/R --vanilla --quiet --args data/Y.txt data/K.txt < R/txt2bin.R > log

You also need to prepare custom VCF files containing the allele specific counts at all SNPs.  The files have to contain an additional subfield, "AS", located within the genotype field consisting of two integers, the reference and alternative allele counts, separated by a comma.  For example, sample **_i_** is heterozygous at a SNP and has 1 and 10 reads overlapping with the reference and alternative alleles respectively, the genotype field for the sample becomes

	... FORMAT ... Sample_i ...
	... GL:AS  ... 0|1:1,10 ...

You can find the example VCF (chr11.gz) in the _data_ directory.  Note that, genotypes have to be phased for QTL mapping with RASQUAL.  Currently, SNP genotypes are only used to map QTLs, but short INDELs and some form of structural variations will be able to use shortly.

## Genotype uncertainty

To maximize the ability of RASQUAL, we recommend to incorporate uncertainty in imputed genotypes.  There are 4 options: 

1. **Allelic probability** (AP) 

    Allelic probability can be obtained from the standard 2-step imputation scheme, where genotypes are first phased then imputed on haplotype-by-haplotype basis.  The custom AP field consists of the allelic probabilities (in Log10 scale) of two haplotypes from one individual, separated by a comma:

        ... FORMAT ... Sample_i      ...
        ... AP:AS  ... 0.0,-5.0:1,10 ...

2. **Genotype likelihood** (GL)

    You may also use genotype likelihood from conventional genotype imputation in conjunction with phased genotype data:

        ... FORMAT    ... Sample_i                ...
        ... GT:GL:AS  ... 0|1:-4.0,-0.1,-0.6:1,10 ...

3. **Dosage** (DS)

    Genotype dosage can also be utilized to take account of genotype uncertainty:

        ... FORMAT    ... Sample_i     ...
        ... GT:DS:AS  ... 0|1:1.1:1,10 ...

4. **Imputation quality scoare** (R squre value; RSQ)

    Imputation methods often provide a quality score for each SNP locus that approximates squared correlation coefficient between true and observed genotypes (*e.g.*, *R*^2 from MaCH or Beagle; *I*^2 from IMPUTE2 ).  RASQUAL can convert the score into genotyping error rate to handle uncertainly:

        ... INFO            FORMAT ... Sample_i ...
        ... ...;RSQ=0.9;... GL:AS  ... 0|1:1,10 ...

    If you want to prioritize RSQ, you need to specify **-z** option for RASQUAL (see above example).

5. **Population allele frequency**

    Even if genotype information is not available, RASQUAL can run in the SNP free mode (**--population-allele-frequency** option).  In this case, you just need to provide the population allele frequency (not minor allele frequency!) for the SNP in VCF files.  Note that, this mode is only feasible when the causal SNP is a feature SNP (*e.g.*, ChIP-seq QTL).

        ... INFO           FORMAT ... Sample_i ...
        ... ...;AF=0.4;... AS     ... 1,10     ...

We strongly recommend to use AP for QTL mapping.  If there are multiple subfields of AP, GL and DS, then AP is prioritized than GL and DS and GL is prioritized than DS.  If neither of GL, DS nor AP is provided, GL is used as AP.

## Offset calculation

Sample specific offset terms can be calculated from the count table.  See the script makeOffset.R in the *R* directory.  The usage is:

	# Does not work!
	$RHOME/R --vanilla --quiet --args data/your.Y.txt < R/makeOffset.R > log

## Covariates

There are usually several confounding factors in the real data, which affects count data and reduces power to detect QTLs (such as sequencing batch, sample preparation date etc.).  RASQUAL can handle covariates as an input (**-x** option).  The following is the same eQTL mapping example above, but with covariates:

    tabix data/chr11.gz 11 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 409 -m 63 \
        -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 -z -t -f C11orf21 \
        -x data/X.bin

The covariate file is based on a sample-by-variable table (see X.txt in the *data* directory).  You can convert the file into binary format by using the same R script:

	$RHOME/R --vanilla --quiet --args data/Y.txt data/K.txt data/X.txt < R/txt2bin.R > log

Those confunding factors are not often observed but can be captured by principal component analysis (PCA).  We applied PCA onto log FPKMs with and without permutation and picked up the first several components whose contribution rates are greater than those from permutation result as covariates for subsequent analyses.  A sample code is also available in the *R* directory:

	# Does not work!
	$RHOME/R --vanilla --quiet --args data/your.Y.txt < R/makeCovariates.R > log


