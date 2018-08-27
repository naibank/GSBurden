#########################################################
#### Install GSBurden ####
library(devtools)
install_github("naibank/GSBurden")
library(GSBurden)
#### read gene set ####
#### text file format (table with at least 2 columns; enzid and gene set name)
#### default separator is "\t" with header. If other separator is used, please define using parameter "sep".
gs <- readGeneset("gs.table.tsv")

#### read cnvs ####
#### any file format but need to have following information; 
#### sample ID, location(chr, start, end) and type of CNVs (optional)
#### each row represent one CNV
cnv.in <- read.delim("CNV.tsv", stringsAsFactors = F)
cnv.in <- cnv.in[which(cnv.in$CNV %in% c("Loss", "Gain") & !is.na(cnv.in$status)), ]
cnvs <- CNVSet(cnv.in$SID, cnv.in$chr, cnv.in$start, cnv.in$end, cnv.in$CNV)

#### read gene annotation ####
#### any file format but need to have following information;
#### loction(chr, start, end) genesymbol and Entrez gene ID)
#### each row represent a gene or an exon
gene.in <- read.delim("gene.tsv", stringsAsFactors = F, header = F)
genes <- GeneAnnotation(gene.in$V6, gene.in$V1, gene.in$V2, gene.in$V3, gene.in$V5)

#### get gene sets matrix ####
#### transform CNV data into gene-set matrix
cnv.matrix <- getCNVGSMatrix(cnvs, genes, gs)

#### add outcome variable, covariates to the matrix ####
sample.info <- unique(cnv.in[, c("SID", "Cohort", "status", "gender", "C1", "C2", "C3", "C4")])
cnv.matrix <- merge(cnv.matrix, sample.info, by.x = "sample", by.y = "SID",
                    all.x = T)
cnv.matrix$status <- ifelse(cnv.matrix$status == "Affected", 1, 0)
covariates <- c("gender", "C1", "C2", "C3")

#### perform global burden and gene set burden test ####
global.test.out <- CNVGlobalTest(cnv.matrix, "status", covariates)
burden.test.out <- CNVBurdenTest(cnv.matrix, gs, "status", covariates)

#### select significant geneset for loci test (optional) ####
geneset <- gs[c("hi015", "TADA1000", "PhHs_MindFun_ADX", "hi035")]

#### perform loci testing. In case you want to test all genes, remove geneset from the parameters ####
loci.test.out <- CNVLociTest(cnvs, cnv.matrix, genes, "status", covariates, geneset)

#### write output to files ####
write.table(global.test.out, "global.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(burden.test.out, "gs.burden.tsv", sep="\t", row.names=F, quote=F, col.names=T)
write.table(loci.test.out, "loci.test.tsv", sep="\t", row.names=F, quote=F, col.names=T)
