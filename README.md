EPIQ
====

This is the beta version of EPIQ - An efficient tool for detecting SNP-SNP epistatic interactions for quantitative traits.
See Arkin2014.pdf for a description of the algorithm.

How to Instal:
--------------
Requirments: Linux, Unix or Mac OS X; g++ version 4.6.3 or higher.
Download files to a local directory. From the terminal, cd to the "AlgorithmRelease" directory and type "make all". Thish should create an executable named "projections".

How to run:
-----------

1. Input: EPIQ requires 3 files: A genotypes file, a phenotypes file, and a map file.
     * Genotypes file: a file containing only 0's or 1's. The file contains one row per individual, one column per SNP. No whitespaces or missing data are allowd.
     * Phenotypes file: contains one row per individual. Each row contains one number.
     * A map file: a 4 column file, in the format of PLINK's map files.
2. Convert files to epi format: 
3. Usage: 
