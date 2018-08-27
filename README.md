# GSBurden

GSBurden is an R package for analysis of human variants (CNVs, SNVs and Indels). The package offers tools to do global burden analysis (testing of number of genes impacted by CNVs), gene-set burden analysis (testing of number of genes in a gene-set) and loci testing (testing of association loci).

In current version, only CNVs analysis is offered. A choice of two models (linear regression and logistic regression) is automatically selected based on outcome variable in the data. Permutation test is done as default but can be disable. BiasedUrn label permutation is utilized to incorporate covariates in label permutation. 
