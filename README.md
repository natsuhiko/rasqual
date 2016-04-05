# RASQUAL
RASQUAL (Robust Allele Specific QUAntification and quality controL) maps QTLs for sequenced based cellular traits by combining population and allele-specific signals.

## How to build & install

**Please make sure CLAPACK and GSL (GNU Scientific Library) are installed in your environment** (if you don't have them, then see below for installation tips).  GSL is usually installed in /usr directory.  Please check /usr/include/gsl and /usr/lib/libgsl.a are existing.

To build and install RASQUAL, firstly go to the _source_ directory (*src*), then set environment variables appropriately to point to the CLAPACK library and GSL.  Finally use "make" to build and install RASQUAL which will be installed in "$RASQUALDIR/bin".

	RASQUALDIR=/path/to/rasqualdir/
	cd $RASQUALDIR/src
	# Not run!  Please export your environment.
	export CFLAGS="-I/path/to/your/CLAPACK-*.*.*.*/INCLUDE -I/path/to/your/CLAPACK-*.*.*.*/F2CLIBS -I/path/to/your/gsl-*.*/gsl"
	export LDFLAGS="-L/path/to/your/CLAPACK-*.*.*.* -L/path/to/your/CLAPACK-*.*.*.*/F2CLIBS -I/path/to/your/gsl-*.*/lib"
	make
	make install

## QTL mapping with RASQUAL

RASQUAL needs two input data files:

1. A fragment (read) count table, with sample specific offsets (such as library size)
2. A VCF file with \*phased\* SNP genotypes and allele-specific counts. 

An example of each of these files can be found in the data directory. In the usage example below, RASQUAL takes SNP information from a tabix-indexed VCF file as standard input, while the count table and sample specific offsets are binary files (Y.bin and K.bin, respectively). Tabix-indexing is not strictly necessary but allows for genotype and allelic count information to be accessed quickly from the command line.
Using the example data files, you can use the following commands to map expression QTLs for two genes (C11orf21 and TSPAN32) using RASQUAL:

    # make sure tabix is installed in your environment
    cd $RASQUALDIR
    tabix data/chr11.gz 11:2315000-2340000 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 378 -m 62 \
        -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 \
        -t -f C11orf21 -z
    tabix data/chr11.gz 11:2315000-2340000 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 2 -l 378 -m 60 \
        -s 2323227,2323938,2324640,2325337,2328175,2329966,2330551,2331219,2334884,2335715,2338574,2339093 \
        -e 2323452,2324188,2324711,2325434,2328220,2330040,2330740,2331248,2334985,2337897,2338755,2339430 \
        -t -f TSPAN32 -z

Sample size (in this example, *N*=24) is given by **-n** option, the feature ID is given by **-j** option (only two genes exist in this example, thereby j=1,2).  You need to provide the number of testing SNPs and feature SNPs in the *cis*-window *a priori* (**-l** and **-m**, respectively).  RASQUAL also requires the feature start and end positions (as comma separated values for more than one positions, e.g. such as for a union of exons in this example) as inputs (**-s** and **-e**, respectively).  By default, RASQUAL outputs QTL mapping results for all tested SNPs, but you can also specify only the lead QTL SNP (**-t** option).  In the output, you can also specify the feature name by **-f** option.  To take account of genotype uncertainty, imputation quality score (R square value) are used in this example (**-z** option; see the section below).  

## Output

On output, RASQUAL provides the following values for each tested SNP:

1. Feature ID*
2. rs ID*
3. Chromosome*
4. SNP position*
5. Ref allele
6. Alt allele
7. Allele frequency (not MAF!)*
8. HWE Chi-square statistic
9. Imputation quality score (IA)
10. Benjamini-Hochberg Q-value
11. Chi square statistic (2 x log Likelihood ratio)*
12. Effect size (Pi)
13. Sequencing/mapping error rate (Delta)
14. Reference allele mapping bias (Phi)
15. Overdispersion
16. SNP ID within the region
17. No. of feature SNPs
18. No. of tested SNPs*
19. No. of iterations for null hypothesis
20. No. of iterations for alternative hypothesis
21. Random location of ties (tie lead SNP; only useful with **-t** option)
22. Log likelihood of the null hypothesis
23. Convergence status (0=success)
24. Squared correlation between prior and posterior genotypes (fSNPs)
25. Squared correlation between prior and posterior genotypes (rSNP)

You may need columns with (*) for the downstream analysis.

## Data preparation

You can find an example expression data for C11orf21 and TSPAN32 genes in the _data_ directory.  There are two files: a read/fragment count table (Y.txt) and sample specific offset data (K.txt), both have to be organised in feature by sample format (*i.e.*, row: gene; col: sample).  We have included an R script to allow you to convert the count and offset files (text) into binary format for us by RASQUAL.  The script converts a table data into a vector of double precision values.

	cd $RASQUALDIR
	RHOME=/software/R-3.0.0/bin/
	$RHOME/R --vanilla --quiet --args data/Y.txt data/K.txt < R/txt2bin.R > log

You will also need to prepare custom VCF files containing the allele specific counts of your target cellular trait at all SNPs.  The files have to contain an additional subfield, "AS", located within the genotype field consisting of two integers, the reference and alternative allele counts, separated by a comma.  For example, sample **_i_** is heterozygous at a SNP and has 1 and 10 reads overlapping the reference and alternative alleles respectively, the genotype field for the sample becomes

	... FORMAT ... Sample_i ...
	... GL:AS  ... 0|1:1,10 ...

An example VCF file (chr11.gz) can be found in the _data_ directory.  Note that, phased genotypes are required for QTL mapping with RASQUAL.  Currently, SNP genotypes are only used to map QTLs, but short INDELs and some form of structural variations will be able to use shortly.

## Genotype uncertainty

To maximise the ability of RASQUAL, we recommend to incorporate uncertainty in imputed genotypes.  There are 4 options: 

1. **Allelic probability** (AP) 

    Allelic probability can be obtained from the standard 2-step imputation scheme, where genotypes are first phased then imputed on haplotype-by-haplotype basis.  The custom AP field consists of the allelic probabilities (in Log10 scale) of two haplotypes from one individual, separated by a comma:

        ... FORMAT ... Sample_i      ...
        ... AP:AS  ... 0.0,-5.0:1,10 ...

2. **Genotype likelihood** (GL)

    You may also use genotype likelihood (in Log10 scale) from conventional genotype imputation in conjunction with phased genotype data:

        ... FORMAT    ... Sample_i                ...
        ... GT:GL:AS  ... 0|1:-4.0,-0.1,-0.6:1,10 ...

    Likewise, you can also provide the genotype likelihood in nominal scale [0-1] with GP FORMAT:

        ... FORMAT    ... Sample_i                ...
        ... GT:GP:AS  ... 0|1:0.01,0.99,0.0:1,10 ...

3. **Dosage** (DS)

    Genotype dosage can also be utilised to take account of genotype uncertainty:

        ... FORMAT    ... Sample_i     ...
        ... GT:DS:AS  ... 0|1:1.1:1,10 ...

4. **Imputation quality score** (R square value; RSQ)

    Imputation methods often provide a quality score for each SNP locus that approximates squared correlation coefficient between true and observed genotypes (*e.g.*, *R*^2 from MaCH or Beagle; *I*^2 from IMPUTE2 ).  RASQUAL can convert the score into genotyping error rate to handle uncertainly:

        ... INFO            FORMAT ... Sample_i ...
        ... ...;RSQ=0.9;... GL:AS  ... 0|1:1,10 ...

    If you want to prioritise RSQ, you need to specify **-z** option for RASQUAL (see above example).

5. **Population allele frequency**

    If genotype information is not available, RASQUAL can run in SNP free mode (**--population-allele-frequency** option).  In this case, you just need to provide the population allele frequency (not minor allele frequency!) for the SNP in VCF files.  Note that, this mode is only feasible when the causal SNP is a feature SNP (*e.g.*, ChIP-seq QTL).

        ... INFO           FORMAT ... Sample_i ...
        ... ...;AF=0.4;... AS     ... 1,10     ...

We strongly recommend using the AP field for QTL mapping.  If there are multiple subfields of AP, GL and DS, AP is prioritised over both GL and DS, and GL is prioritised over DS.  If neither of GL, DS nor AP is provided, GL is used as AP.

## Offset calculation

Sample specific offset terms (K.txt) can be calculated from the count table.  See the script makeOffset.R in the *R* directory.  The usage is:

	# Not run!
	$RHOME/R --vanilla --quiet --args data/your.Y.txt < R/makeOffset.R > log
	# With GC content; not run!
	$RHOME/R --vanilla --quiet --args data/your.Y.txt data/gcc.txt < R/makeOffset.R > log

Note that you need to prepare a GC content file (gcc.txt in this example) to apply GC correction for the read count at each feature.  The file is a vector of GC% values for all features as a text file (separated by either, a comma, a tab or a line break).  In order to obtain the GC% for each feature, we normally extract the reference sequence overlapping with the feature annotation, count G/C bases and then divide the count by the total feature length.

## Covariates

Real data is usually affected by hidden confounding factors, such as sequencing batch, sample preparation date etc, that can reduce power to detect QTLs. RASQUAL handles covariates as an input (**-x** option).  The following is the same eQTL mapping example above, but with covariates:

    tabix data/chr11.gz 11:2315000-2340000 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 378 -m 62 \
        -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 \
        -z -t -f C11orf21 \
        -x data/X.bin

The covariate file is based on a sample-by-variable table (see X.txt in the *data* directory).  You can convert the file into binary format for RASQUAL by using the same R script:

	$RHOME/R --vanilla --quiet --args data/Y.txt data/K.txt data/X.txt < R/txt2bin.R > log

Those confounding factors are not often observed but can be captured by principal component analysis (PCA).  In our example, we applied PCA onto log FPKMs with and without permutation and picked the first N components whose explained variances are greater than those from permutation result as covariates for subsequent analyses.  A sample code is also available in the *R* directory:

	# Not run!
	$RHOME/R --vanilla --quiet --args data/your.Y.txt data/your.K.txt < R/makeCovariates.R > log

Note that, the result of PCA is always sensitive to few outliers (just one or two), which explain almost all variance in the data.  Using those PCs as covariates hurts your QTL mapping result.  We strongly recommend to spend some time to explore and clean up your data first.

## Installation tips for CLAPACK and GSL

You first need to get the latest library from http://www.netlib.org/clapack/.  Then, compile it like

	tar zxvf clapack.tgz
	cd CLAPACK-3.2.1
	mv make.inc.example make.inc
	make

When it has been done, you will find three archive files in the directory which need to be either linked or renamed, such as

	ln -s lapack_LINUX.a liblapack.a
	ln -s tmglib_LINUX.a libtmglib.a
	ln -s blas_LINUX.a libblas.a

before compiling RASQUAL.

You may also need to get GSL (GNU Scientific Library) from http://www.gnu.org/software/gsl/ (if it is not installed).  Then, compile it like

	tar zxvf gsl-*.*.tar.gz
	cd gsl-*.*
	./configure --prefix=$PWD
	make
	make install

## Creating a VCF file with AS counts

We provide a useful script to create a VCF file with AS counts from a master VCF file and a set of BAM files.  Before using the script you need to compile some C codes called from the script:

	cd $RASQUALDIR/src/ASVCF
	make

The command to produce the VCF with AS counts is:

	sh $RASQUALDIR/src/ASVCF/createASVCF.sh bam.list.txt master.vcf.gz

which creates *master.vcf.new.gz* in the same directory.  The *master.vcf.gz* must be tabix indexed and the *bam.list.txt* is a text file which contains absolute path to your set of BAMs from which AS counts are produced:

	# bam.list.txt
	/path/to/your/bam/sample1.bam 
	/path/to/your/bam/sample2.bam 
	/path/to/your/bam/sample3.bam
	...

The order of the samples **MUST** be the same as that in the master VCF.  Before using the script, please make sure the latest tabix (http://www.htslib.org/doc/tabix.html) is installed in your environment. 

Note that, our script doesn't filter out any AS read by means of QC criteria (depth of coverage, mapping quality, etc.).  You may also want to filter out some reads a priori, using GATK ASEReadCounter (https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php).

## Permutation test

The permutation test is also implemented in RASQUAL.  The **-r/--random-permutation** option generates a random permutation for each feature to break the correlation between genotype and total feature count as well as AS counts.

## Tips to speed up QTL mapping

### Multithreading

RASQUAL is now multithreaded in order to speed up execution times, which requires the **pthread** library.  You just need to specify the additional option **--n-thereads** to use this function.  Using the above example, you can use the following commands to map the eQTL with 10 threads:

    tabix data/chr11.gz 11:2315000-2340000 | bin/rasqual -y data/Y.bin -k data/K.bin -n 24 -j 1 -l 409 -m 63 \
        -s 2316875,2320655,2321750,2321914,2324112 -e 2319151,2320937,2321843,2323290,2324279 \
        -t -f C11orf21 -z --n-threads 10

### Filters for fSNPs

To maximize power to detect QTLs, RASQUAL uses all fSNPs with MAF>0.0, pHWE>0.0 and imputation quality score RSQ>0.0.  However, RASQUAL takes ages to map QTLs with a number of fSNPs in a feature (e.g., long genes).  Therefore you may want to reduce the number of fSNPs with additional filters.  We introduced the following new options **--minor-allele-frequency-fsnp**, **--imputation-quality-fsnp** and **--hardy-weinberg-pvalue-fsnp** to eliminate some of fSNPs which are possibly not so informative.
        
## Conditional analysis

To map subsidiary QTLs conditional on the lead QTL variant(s) can be performed with the following option:

    bin/rasqual ... -k2 rs0001:0.1,rs0002:0.2 ...
    
You may introduce any number of variants with thier effect sizes (Pi values) as comma separated values where each variant ID and its Pi value have to be connected by colon (:).

## Warnings

To save the memory, each variant ID in the VCF file must be shorter than 100 characters; otherwise a buffer overflow happens.
