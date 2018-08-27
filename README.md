# GSBurden

GSBurden is an R package for gene-set enrichment analysis of human variants (CNVs, SNVs and Indels). The package offers tools to do global burden analysis (testing of number of genes impacted by CNVs), gene-set burden analysis (testing of number of genes in a gene-set) and loci testing (testing of association loci).

In current version, only CNVs analysis can be done using the package. A choice of two models (linear regression and logistic regression) is automatically selected based on outcome variable in the data. Permutation test is done as default but can be disabled to make it faster. BiasedUrn label permutation is utilized to incorporate covariates in label permutation. 

The analysis can be done as demonstrated in example.R
