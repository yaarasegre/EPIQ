EPIQ
====

This is the beta version of EPIQ - An efficient tool for detecting SNP-SNP epistatic interactions for quantitative traits.
See Arkin2014.pdf for a description of the algorithm.

How to Install:
--------------
Requirements: Linux, Unix or Mac OS X; g++ version 4.6.3 or higher.
Download files to a local directory. From the terminal, cd to the "AlgorithmRelease" directory and type "make all". This should create an executable named "projections".

How to run:
-----------

1. Input: EPIQ requires 5 files: 
     * Genotypes file: a file containing only 0's or 1's. The file contains one row per individual, one column per SNP. No whitespaces or missing data are allowed.
     * Phenotypes file: contains one row per individual. Each row contains one number.
     * A map file: a 4 column file, in the format of PLINK's map files.
     * An edges file: defines how SNPs are distributed to bins, according to the minor allele frequencies. Contains B+1 values, separated by a comma, where B is the number of bins.
     * A thresholds file: A CSV, containing B rows and B columns, where B is the number of bins. Each value is the threshold used for the corresponding pair of bins. 
     
2. Converting genotype files to EPIQ format: 
	* Use PLINK to convert data to a "raw" format: 
	
        	plink --bfile <fname> --recodeA --noweb --out temp

	* Edit the raw file output from the last command to get an EPIQ formatted file: 

		A. Remove first 6 cols:
		
        	cut -f 7- -d ' ' temp.raw  > data.epi

		B. Remove blanks and convert to binary. For example, replace 2's with 1's:
		
	    	find . -name data.epi | xargs perl -pi -e 's/ //g'
	    	find . -name data.epi | xargs perl -pi -e 's/2/1/g'

		C. Resolve missing data. For example:
		
	    	find . -name data.epi | xargs perl -pi -e 's/NA/0/g'

3. Usage: 

            ./projections -p \<phenotype file\> -f \<genotypes file\>  -F  \<thresholds file\> -e \<edges file\> 
            -n \<# of individuals\>  -L \<# of iterations\> -C \<chi-squared threshold\> -M \<map file\>

