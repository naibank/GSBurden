# GSBurden

GSBurden is an R package for genotype-phenotype gene-set enrichment analysis of several genomic variants (CNVs, SNVs and Indels). The package offers tools to do global burden analysis (testing of the number of genes impacted by CNVs or the number of CNVs), gene-set burden analysis (testing of the number of genes within a gene-set) and loci testing (testing of a genomic locus).

In current version, only CNVs analysis can be done using the package (SNVs analysis can also be done but the pipeline is being tested). A choice of three models (linear regression, logistic regression and ordinal logistic regression) is automatically selected based on the given outcome variable in the data. Permutation test is done as default but can be disabled to make it faster. BiasedUrn label permutation could also be utilized to incorporate covariates in the label permutation. 

The analysis can be done as demonstrated in example.R

# Installation
1. GSBurden requires R version 3.3 or later
2. GSBurden requires GenomicRanges (>= 1.30.3), BiasedUrn (>= 1.07), MASS (>= 7.3-49), ordinal (>= 2018.8-25), fastmatch (1.1-8), speedglm (0.3-5)
3. Install the package through R script as below 

devtools::install_github("naibank/GSBurden")

**Any data preprocessing steps of CNVs/SNVs data are not included in the library. 
**For CNVs, it is advisible that large CNVs (i.e. size > 3MB) and CNVs significantly overlapped segmental duplications (i.e. > 70% overlap) be excluded from the analysis. 

# Notes on CNV loci test
The test is done in a gene-based fashion. However, as CNVs oftenly overlap multiple genes, the loci are merged if 80% (found to generated FWER adjusted p-value that is approximate to that of permutation FWER) of the overlapped CNVs are common between loci. See image below for an example.

![Loci test](https://github.com/naibank/GSBurden/blob/master/material/loci%20test.jpg)
