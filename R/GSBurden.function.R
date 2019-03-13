#' readGeneset Function
#'
#' This function to read gene set information to perform gene set burden test for CNVs and SNVs.
#' @param pathfile path to gene set file. Either 1) A text file of two-column table; 1st column is Entrez gene ids and 2nd column is Gene set name. or 2) R object file, contain list variable. Each element is a vector contains Entrez gene ids. Gene set name is required as names of list elements.
#' @keywords gene set, entrez gene ids, gene set analysis
#' @export 
#' @examples
#' readGeneset('test.geneset.RData')
readGeneset <- function(pathfile, sep="\t", header=T){
  out <- tryCatch(
    {
      return(get(load(pathfile)))
    },
    error=function(cond) {
    },
    warning=function(cond) {
      dt <- read.delim(pathfile, sep = sep, header = header, stringsAsFactors = F)
      names(dt) <- c("enzid", "gsname")
      
      list.out <- list()
      for(gs in unique(dt$gsname)){
        list.out[[gs]] <- na.omit(dt$enzid[which(dt$gsname == gs)])
      }
      return(list.out)
    },
    finally={
      message("Read gene set file successfully")
    }
  )
  
  return(out)
}

#'
#' This function to get CNV in the format used by GSBurden.
#' @param sample a vector of sample IDs
#' @param chr a vector of chromosome name of CNVs
#' @param start a vector of start locations of CNVs
#' @param end a vector of end locations of CNVs
#' @param type a vector of type of CNVs i.e. c("Loss", "Gain", "Loss", "Loss")
#' @keywords Entrez id
#' @export
CNVSet <- function(sample, chr, start, end, type = "-"){
  if(length(type) == 1){
    message("found only one type of CNVs")
    if(sd(c(length(sample), length(chr), length(start), length(end))) != 0){
      stop("Supplied parameters do not have same length")
    }else{
      return(data.frame(sample, chr, start, end, type, stringsAsFactors = F))
    }
  }else{
    if(sd(c(length(sample), length(chr), length(start), length(end), length(type))) != 0){
      stop("Supplied parameters do not have same length")
    }else{
      return(data.frame(sample, chr, start, end, type, stringsAsFactors = F))
    }
  }
  
}

#'
#' This function to get coordinate of gene for annotating CNVs/SNVs in the format used by GSBurden.
#' @param enzid a vector of Entrez gene IDs
#' @param chr a vector of chromosome name of CNVs
#' @param start a vector of start locations of CNVs
#' @param end a vector of end locations of CNVs
#' @keywords Entrez id
#' @export
GeneAnnotation <- function(enzid, chr, start, end, gsymbol = "-"){
  if(sd(c(length(enzid), length(chr), length(start), length(end))) != 0){
    stop("Supplied parameters do not have same length")
  }else{
    return(data.frame(enzid, chr, start, end, gsymbol, stringsAsFactors = F))
  }
}

#'
#' This function is for internal use by other functions in the package.
#' @param enzid Entrez id
#' @param geneset geneset
#' @keywords GSBurden
#' @export
getGenesetCount <- function(enzid, geneset){
  return(sapply(sapply(geneset, is.element, enzid), sum))
}

#'
#' This function is for internal use by other functions in the package.
#' @param enzid Entrez id
#' @param geneset geneset
#' @keywords GSBurden
#' @export
getGenesetList <- function(enzid, geneset){
  gs.sum <- sapply(sapply(geneset, is.element, enzid), sum)
  
  return(paste(names(gs.sum)[gs.sum > 0], collapse = ";"))
}

#'
#' This function is to get a matrix of gene set count by sample.
#' @param cnv.table CNV table
#' @param annotation.table gene annotation table
#' @param geneset gene set object
#' @keywords GSBurden
#' @export
getCNVGSMatrix <- function(cnv.table, annotation.table, geneset){
  all.out <- data.frame()
  cnvtypes <- unique(cnv.table$type)
  for(cnvtype in cnvtypes){
    this.cnv.table <- cnv.table[cnv.table$type == cnvtype, ]
    cnv.g <- GenomicRanges::GRanges(this.cnv.table$chr, IRanges::IRanges(this.cnv.table$start, this.cnv.table$end), "*")
    annotation.g <- GenomicRanges::GRanges(annotation.table$chr, IRanges::IRanges(annotation.table$start, annotation.table$end), "*")
    
    olap <- data.frame(IRanges::findOverlaps(cnv.g, annotation.g))
    olap$sample <- this.cnv.table$sample[olap$queryHits]
    olap$enzid <- annotation.table$enzid[olap$subjectHits]
    
    cnvCount <- table(this.cnv.table$sample[unique(olap$queryHits)])
    
    olap <- unique(olap[, c("sample", "enzid")])
    geneCount <- table(olap$sample)
    temp <- aggregate(enzid ~ sample, olap, c)
    gsMatrix <- data.frame(t(sapply(temp$enzid, getGenesetCount, geneset)))
    gsMatrix$sample <- temp$sample
    
    cnvCount <- data.frame(sample = names(cnvCount), cnv_count = as.numeric(cnvCount))
    geneCount <- data.frame(sample = names(geneCount), gene_count = as.numeric(geneCount))
    
    dt.out <- merge(cnvCount, geneCount, by = "sample", all.x = T)
    dt.out <- merge(dt.out, gsMatrix, by = "sample", all.x = T)
    dt.out[is.na(dt.out)] <- 0
    
    if(length(cnvtypes) > 1){
      names(dt.out)[-1] <- paste(names(dt.out)[-1], cnvtype, sep="_")
    }
    
    if(nrow(all.out) == 0){
      all.out <- dt.out
    }else{
      all.out <- merge(all.out, dt.out, by = "sample", all = T)
    }

    all.out[is.na(all.out)] <- 0
  }
  
  message("Transform gene set count matrix successfully")
  return(all.out)
}

#'
#' This function is to get a matrix of gene set count by sample but separate duplication by disruption.
#' @param cnv.table CNV table
#' @param annotation.table gene annotation table
#' @param appris.cds appris principal isoform cds. Each row represent an exon but only coding sequence.
#' @param appris.nodisrupt appris principal isoform promoter and stop codon.
#' @param geneset gene set object
#' @keywords GSBurden Separate duplication by disruption
#' @export
getCNVDupGSMatrix <- function(cnv.table, annotation.table, appris.cds, appris.nodisrupt, geneset){
  this.cnv.table <- cnv.table
  all.out <- data.frame()
  cnv.g <- GenomicRanges::GRanges(this.cnv.table$chr, IRanges::IRanges(this.cnv.table$start, this.cnv.table$end), "*")
  annotation.g <- GenomicRanges::GRanges(annotation.table$chr, IRanges::IRanges(annotation.table$start, annotation.table$end), "*")
  appris.cds.g <- GenomicRanges::GRanges(appris.cds$chr, IRanges::IRanges(appris.cds$start, appris.cds$end),"*")
  appris.nodisrupt.g <- GenomicRanges::GRanges(appris.nodisrupt$chr, IRanges::IRanges(appris.nodisrupt$start, appris.nodisrupt$end),"*")
  
  olap.nodisrupt <- data.frame(IRanges::findOverlaps(cnv.g, appris.nodisrupt.g))
  olap.nodisrupt$elementsize <- appris.nodisrupt$end[olap.nodisrupt$subjectHits]-appris.nodisrupt$start[olap.nodisrupt$subjectHits] + 1
  olap.nodisrupt$olapsize <- IRanges::width(IRanges::pintersect(cnv.g[olap.nodisrupt$queryHits], appris.nodisrupt.g[olap.nodisrupt$subjectHits]))
  olap.nodisrupt <- olap.nodisrupt[olap.nodisrupt$elementsize == olap.nodisrupt$olapsize, ]
  olap.nodisrupt$gsymbol <- appris.nodisrupt$gsymbol[olap.nodisrupt$subjectHits]
  olap.nodisrupt$key <- paste(olap.nodisrupt$queryHits, olap.nodisrupt$gsymbol, sep=":")
  
  olap.appris <- data.frame(IRanges::findOverlaps(cnv.g, appris.cds.g))
  olap.appris$gsymbol <- appris.cds$gsymbol[olap.appris$subjectHits]
  olap.appris$key <- paste(olap.appris$queryHits, olap.appris$gsymbol, sep=":")
  olap.appris <- olap.appris[!olap.appris$key %in% olap.nodisrupt$key, ]
  
  # olap.appris$size <- GenomicRanges::width(GenomicRanges::pintersect(cnv.g[olap.appris$queryHits], 
  #                                                                    appris.g[olap.appris$subjectHits]))
  # olap.appris$possible.size <- GenomicRanges::width(appris.g[olap.appris$subjectHits])
  # olap.appris.nodis <- olap.appris[olap.appris$size >= olap.appris$possible.size, ]
  # olap.appris <- olap.appris[olap.appris$size < olap.appris$possible.size, ]
  
  olap <- data.frame(IRanges::findOverlaps(cnv.g, annotation.g))
  
  for(i in c("NoDisruptDup", "DisruptDup")){
    if(i == "NoDisruptDup"){
      th.olap <- olap #olap[!olap$queryHits %in% olap.appris$queryHits, ]
      info.table <- annotation.table
      th.olap$enzid <- info.table$enzid[th.olap$subjectHits]
      th.olap$gsymbol <- info.table$gsymbol[th.olap$subjectHits]
      th.olap <- th.olap[!paste(th.olap$queryHits, th.olap$gsymbol, sep=":") %in% 
                           paste(olap.appris$queryHits, olap.appris$gsymbol, sep=":"), ]
    }else{
      th.olap <- olap.appris
      th.olap <- merge(th.olap, annotation.table[, c("gsymbol", "enzid")], by = "gsymbol", all.x = T)
    }
    
    th.olap$sample <- this.cnv.table$sample[th.olap$queryHits]
    cnvCount <- table(th.olap$sample[!duplicated(th.olap$queryHits)])
    
    th.olap <- unique(th.olap[, c("sample", "enzid")])
    geneCount <- table(th.olap$sample)
    temp <- aggregate(enzid ~ sample, th.olap, c)
    gsMatrix <- data.frame(t(sapply(temp$enzid, getGenesetCount, geneset)))
    gsMatrix$sample <- temp$sample
    
    th.cnvCount <- data.frame(sample = names(cnvCount), cnv_count = as.numeric(cnvCount))
    th.geneCount <- data.frame(sample = names(geneCount), gene_count = as.numeric(geneCount))
    
    dt.out <- merge(th.cnvCount, th.geneCount, by = "sample", all.x = T)
    dt.out <- merge(dt.out, gsMatrix, by = "sample", all.x = T)
    dt.out[is.na(dt.out)] <- 0
    
    names(dt.out)[-1] <- paste(names(dt.out)[-1], i, sep="_")
    
    if(nrow(all.out) == 0){
      all.out <- dt.out
    }else{
      all.out <- merge(all.out, dt.out, by = "sample", all = T)
    }
    
    all.out[is.na(all.out)] <- 0
  }
  
  message("Transform gene set count matrix successfully")
  return(all.out)
}

#'
#' This function is to test for global burden of gene count.
#' @param cnv.matrix CNV table containing samples' outcome, gene count, covariates (optional) and cnv count (optional)
#' @param label variable name that is used as outcome
#' @param covariates list of covariates to be used in the model
#' @param correctCNVCount logical value to indicate whether the burden will be corrected for CNV count or not
#' @param standardizeCoefficient logical value to indicate whether coefficient will be standardized or not
#' @keywords GSBurden
#' @export
#' 
CNVGlobalTest <- function(cnv.matrix, label, covariates, correctCNVCount = T, standardizeCoefficient = T){
  
  distinct.prefixes <- names(cnv.matrix)[grep("gene_count", names(cnv.matrix))]
  distinct.prefixes <- gsub(sprintf("%s_", "gene_count"), "", distinct.prefixes)
  
  model = "lm"
  if(length(unique(cnv.matrix[, label])) == 2){
    message("dichotomous outcome variable detected. Logistic regression is being done ...")
    model = "glm"
  }else if(is.numeric(cnv.matrix[, label])){
    message("continuous outcome variable detected. Linear regression is being done ...")
    model = "lm"
  }else if(is.factor(cnv.matrix[, label])){
    message("ordinal outcome variable detected. Ordinal regression is being done ...")
    model = "clm"
  }else{
    stop("Non dichotomous or continuous or ordinal outcome variable detected. The burden test cannot be run")
  }
  
  test.out <- data.frame()
  for(cnvtype in distinct.prefixes){
    feature <- sprintf("gene_count_%s", cnvtype)
    cnvcount <- sprintf("cnv_count_%s", cnvtype)
    
    this.covariates <- covariates
    if(correctCNVCount){
      this.covariates <- c(this.covariates, cnvcount)
    }
    
    ref.term <- sprintf("%s ~ %s", label, paste(this.covariates, collapse = " + "))
    add.term <- sprintf("%s + %s", ref.term, feature)
    
    if(standardizeCoefficient){
      cnv.matrix[, feature] <- scale(cnv.matrix[, feature])
    }
    
    if(model == "lm"){
      ref.model <- lm(ref.term, cnv.matrix)
      add.model <- lm(add.term, cnv.matrix)
      ano <- anova(ref.model, add.model, test = "Chisq"))
    }else if(model == "glm"){
      ref.model <- glm(ref.term, cnv.matrix, family = binomial(link = "logit"))
      add.model <- glm(add.term, cnv.matrix, family = binomial(link = "logit"))
      ano <- anova(ref.model, add.model, test = "Chisq"))
    }else{
      ref.model <- ordinal::clm(ref.term, data = cnv.matrix)
      add.model <- ordinal::clm(add.term, data = cnv.matrix)
      ano <- anova(ref.model, add.model))
    }
    
    names(ano)[length(names(ano))] <- "pvalue"
    pvalue <- ano$pvalue[2]
    coefficient <- add.model$coefficients[feature]
    intervals <- confint(add.model)
    upperbound <- intervals[feature, "97.5 %"]
    lowerbound <- intervals[feature, "2.5 %"]
    
    temp.out <- data.frame("global" = "Gene", "type" = cnvtype, "coefficient" = coefficient, lowerbound, upperbound, "pvalue" = pvalue)
    test.out <- rbind(test.out, temp.out)
  }
  
  return(test.out)
}


#'
#' This function is to test for geneset burden of gene count.
#' @param cnv.matrix CNV table containing samples' outcome, gene count, covariates (optional) and cnv count (optional)
#' @param geneset a list of gene set, which each element contains a list of Entrez gene IDs
#' @param label variable name that is used as outcome
#' @param covariates list of covariates to be used in the model
#' @param correctCNVCount logical value to indicate whether the burden will be corrected for CNV count or not
#' @param standardizeCoefficient logical value to indicate whether coefficient will be standardized or not
#' @param BiasedUrn logical value to indicate whether BaisedUrn will be used to permute label or not
#' @keywords GSBurden
#' @export
#' 
CNVBurdenTest <- function(cnv.matrix, geneset, label, covariates, correctGlobalBurden = T, standardizeCoefficient = T, 
                          permutation = T, nperm = 100, BiasedUrn = F){
  
  distinct.prefixes <- names(cnv.matrix)[grep(names(geneset)[1], names(cnv.matrix))]
  distinct.prefixes <- gsub(sprintf("%s_", names(geneset)[1]), "", distinct.prefixes)
  
  model = "lm"
  if(length(unique(cnv.matrix[, label])) == 2){
    message("dichotomous outcome variable detected. Logistic regression is being done ...")
    model = "glm"
  }else if(is.numeric(cnv.matrix[, label])){
    message("continuous outcome variable detected. Linear regression is being done ...")
    model = "lm"
  }else if(is.factor(cnv.matrix[, label])){
    message("ordinal outcome variable detected. Ordinal regression is being done ...")
    model = "clm"
  }else{
    stop("Non dichotomous or continuous outcome variable detected. The burden test cannot be run")
  }
  
  if(permutation){
    ref.term <- sprintf("%s ~ %s", label, paste(covariates, collapse = " + "))
    message(sprintf("Permuting sample labels for %s times", nperm))
    if(BiasedUrn){
      lm.odds <- glm(ref.term, cnv.matrix, family = binomial (link = "logit"))
      d.odds <- exp(lm.odds$linear.predictors)

      n.case <- sum(cnv.matrix[, label])
      n.all <- length(cnv.matrix[, label])
      
      perm.hg <- BiasedUrn::rMFNCHypergeo(nran = nperm, m = rep(1, n.all), n = n.case, odds = d.odds)
    }else{
      perm.hg <- data.frame(1:nrow(cnv.matrix), nrow(cnv.matrix):1)
      for(i in 1:nperm){
        permuted <- sample(cnv.matrix[, label])
        perm.hg <- cbind(perm.hg, permuted)
      }
      
      perm.hg <- perm.hg[, -c(1:2)]
    }
    
  }
  
  perm.test.pvalues <- data.frame()
  test.out <- data.frame()
  for(cnvtype in distinct.prefixes){
    for(this.gs in names(geneset)){
      out.message <- sprintf("Testing %s", this.gs)
      if(length(distinct.prefixes) > 1){
        out.message <- sprintf("%s for %s CNVs", out.message, cnvtype)
      }
      message(out.message)
      
      feature <- sprintf("%s_%s", this.gs, cnvtype)
      global <- sprintf("gene_count_%s", cnvtype)
      
      this.covariates <- covariates
      if(correctGlobalBurden){
        this.covariates <- c(this.covariates, global)
      }
      
      ref.term <- sprintf("%s ~ %s", label, paste(this.covariates, collapse = " + "))
      add.term <- sprintf("%s + %s", ref.term, feature)
      
      if(standardizeCoefficient & mean(cnv.matrix[, feature]) != 0 & sd(cnv.matrix[, feature]) != 0){
        cnv.matrix[, feature] <- scale(cnv.matrix[, feature])
      }
      
      if(model == "lm"){
        ref.model <- lm(ref.term, cnv.matrix)
        add.model <- lm(add.term, cnv.matrix)
        ano <- anova(ref.model, add.model, test = "Chisq")
      }else if(model == "glm"){
        ref.model <- glm(ref.term, cnv.matrix, family = binomial(link = "logit"))
        add.model <- glm(add.term, cnv.matrix, family = binomial(link = "logit"))
        ano <- anova(ref.model, add.model, test = "Chisq")
      }else{
        ref.model <- ordinal::clm(ref.term, data = cnv.matrix)
        add.model <- ordinal::clm(add.term, data = cnv.matrix)
        ano <- anova(ref.model, add.model)
      }
      
      names(ano)[length(names(ano))] <- "pvalue"
      pvalue <- ano$pvalue[2]
      coefficient <- add.model$coefficients[feature]
      conf <- confint(add.model)
      coeff.l <- conf[feature, 1]
      coeff.u <- conf[feature, 2]
      
      temp.out <- data.frame("geneset" = this.gs, "type" = cnvtype, "coefficient" = coefficient, 
                             "coeff.upper" = coeff.u, "coeff.lower" = coeff.l, "pvalue" = pvalue)
      test.out <- rbind(test.out, temp.out)
      
      if(permutation){
        for(iperm in 1:nperm){
          cnv.matrix$outcome.perm <- perm.hg[, iperm]
          ref.perm.term <- sprintf("outcome.perm ~ %s", paste(this.covariates, collapse = " + "))
          add.perm.term <- sprintf("%s + %s", ref.perm.term, feature)
          
          if(model == "lm"){
            ref.perm.model <- lm(ref.perm.term, cnv.matrix)
            add.perm.model <- lm(add.perm.term, cnv.matrix)
            ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
          }else if(model == "glm"){
            ref.perm.model <- glm(ref.perm.term, cnv.matrix, family = binomial(link = "logit"))
            add.perm.model <- glm(add.perm.term, cnv.matrix, family = binomial(link = "logit"))
            ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
          }else{
            ref.perm.model <- ordinal::clm(ref.perm.term, data = cnv.matrix)
            add.perm.model <- ordinal::clm(add.perm.term, data = cnv.matrix)
            ano.perm <- anova(ref.perm.model, add.perm.model)
          }
          
          names(ano.perm)[length(names(ano.perm))] <- "pvalue"
          coeff <- add.perm.model$coefficients[feature]
          perm.test.pvalues <- rbind(perm.test.pvalues, data.frame("cnvtype" = cnvtype, "pvalue" = ano.perm$pvalue[2], coeff))
        }
      }
    }
  }
  
  if(permutation){
    message("Calculating permutation-based FDR")
    test.out$permFDR <- 1
    for(i in 1:nrow(test.out)){
      rec <- test.out[i, ]
      this.perm <- perm.test.pvalues[perm.test.pvalues$cnvtype == rec$type, ]
      
      actual <- sum(test.out$pvalue[test.out$type == rec$type] <= rec$pvalue)/nrow(test.out[test.out$type == rec$type,])
      perm <- sum(this.perm$pvalue <= rec$pvalue)/nrow(this.perm)
      
      fdr <- ifelse(perm/actual > 1, 1, perm/actual)

      test.out$permFDR[i] <- fdr
    }
  }
  
  test.out <- test.out[order(test.out$pvalue), ]
  for(i in 1:nrow(test.out)){
    test.out$permFDR[i] <- min(test.out$permFDR[i:nrow(test.out)])
  }
  
  list.out <- list(test.out, perm.test.pvalues)
  names(list.out) <- c("Test", "Permutation.Test")
  return(list.out)
}


#'
#' This function is for loci test.
#' @param cnv.table CNV table
#' @param cnv.matrix CNV matrix containing samples' outcome, gene count, covariates (optional) and cnv count (optional)
#' @param annotation.table gene annotation table
#' @param label variable name that is used as outcome
#' @param covariates list of covariates to be used in the model
#' @param geneset a list of gene set, which each element contains a list of Entrez gene IDs. If geneset parameter is specified, only genes in geneset will be tested
#' @param permutation logical, doing permutation FDR or not
#' @param nperm number of permutation to be done
#' @param nsubject minimum number of subjects with CNVs impacting the loci, required for a locus to be tested
#' @param BiasedUrn logical value to indicate whether BaisedUrn will be used to permute label or not
#' @keywords GSBurden
#' @export
#' 
CNVLociTest <- function(cnv.table, cnv.matrix, annotation.table, label, covariates, geneset = list(),
                        permutation = T, nperm = 100, nsubject = 3, BiasedUrn = F){
  cnvtypes <- unique(cnv.table$type)
  all.cnv.table <- cnv.table
  model = "lm"
  if(length(unique(cnv.matrix[, label])) == 2){
    message("dichotomous outcome variable detected. Logistic regression is being done ...")
    model = "glm"
  }else if(is.numeric(cnv.matrix[, label])){
    message("continuous outcome variable detected. Linear regression is being done ...")
    model = "lm"
  }else if(is.factor(cnv.matrix[, label])){
    message("ordinal outcome variable detected. Ordinal regression is being done ...")
    model = "clm"
  }else{
    stop("Non dichotomous or continuous outcome variable detected. The burden test cannot be run")
  }
  
  
  if(length(geneset) != 0){
    annotation.table <- annotation.table[annotation.table$enzid %in% unlist(geneset), ]
    # map.geneset <- sapply(enzid, GSBurden::getGenesetList, geneset)
    # map.geneset <- data.frame(enzid, "geneset" = map.geneset) 
  }
  
  if(permutation){
    ref.term <- sprintf("%s ~ %s", label, paste(covariates, collapse = " + "))
    message(sprintf("Permuting sample labels for %s times", nperm))
    if(BiasedUrn){
      lm.odds <- glm(ref.term, cnv.matrix, family = binomial (link = "logit"))
      d.odds <- exp(lm.odds$linear.predictors)
      
      n.case <- sum(cnv.matrix[, label])
      n.all <- length(cnv.matrix[, label])
      
      perm.hg <- BiasedUrn::rMFNCHypergeo(nran = nperm, m = rep(1, n.all), n = n.case, odds = d.odds)
    }else{

      perm.hg <- data.frame(1:nrow(cnv.matrix), nrow(cnv.matrix):1)
      for(i in 1:nperm){
        permuted <- sample(cnv.matrix[, label])
        perm.hg <- cbind(perm.hg, permuted)
      }
      
      perm.hg <- perm.hg[, -c(1:2)]
    }
    
  }
  
  final.out <- data.frame()
  for(cnvtype in cnvtypes){
    message(sprintf("Testing %s CNVs...", cnvtype))
    
    cnv.table <- all.cnv.table[all.cnv.table$type == cnvtype, ]
  
    cnv.g <- GenomicRanges::GRanges(cnv.table$chr, IRanges::IRanges(cnv.table$start, cnv.table$end), "*")
    annotation.g <- GenomicRanges::GRanges(annotation.table$chr, IRanges::IRanges(annotation.table$start, annotation.table$end), "*")
    olap <- data.frame(IRanges::findOverlaps(cnv.g, annotation.g))
    table <- aggregate(queryHits ~ subjectHits, olap, paste, collapse = ",")
    
    ########### generate loci ############
    new.loci <- data.frame()
    for(i in unique(table$queryHits)){
      temp.loci <- table[table$queryHits == i, ]
      this.samples <- unique(cnv.table$sample[as.numeric(as.character(strsplit(i, ",")[[1]]))])
      if(length(this.samples) >= nsubject){
        sample <- paste(this.samples, collapse = ",")
        temp.loci <- annotation.table[temp.loci$subjectHits, ]
        temp.loci$enzid <- paste(unique(temp.loci$enzid), collapse = ",")
        temp.loci$gsymbol <- paste(unique(temp.loci$gsymbol), collapse = ",")
        temp.loci$start <- min(temp.loci$start)
        temp.loci$end <- max(temp.loci$end)
        temp.loci$sample <- sample
        temp.loci <- unique(temp.loci)
        
        new.loci <- rbind(new.loci, temp.loci)
      }
    }
    
    current.count <- 1
    all.count <- nrow(new.loci)
    dt.out <- data.frame()
    for(iloci in 1:nrow(new.loci)){
      temp.out <- data.frame()
      this.loci <- new.loci[iloci, ]
      
      dt.temp <- cnv.matrix[, c("sample", label, covariates)]
      dt.temp$gene_count <- 0
      sample.with.loci <- strsplit(this.loci$sample,",")[[1]]
      dt.temp$gene_count[which(dt.temp$sample %in% sample.with.loci)] <- 1
            
      ref.term <- sprintf("%s ~ %s", label, paste(covariates, collapse = " + "))
      add.term <- sprintf("%s + %s", ref.term, "gene_count")
            
      if(model == "lm"){
        ref.model <- lm(ref.term, dt.temp)
        add.model <- lm(add.term, dt.temp)
        ano <- anova(ref.model, add.model, test = "Chisq")
      }else if(model == "glm"){
        ref.model <- glm(ref.term, dt.temp, family = binomial(link = "logit"))
        add.model <- glm(add.term, dt.temp, family = binomial(link = "logit"))
        ano <- anova(ref.model, add.model, test = "Chisq")
      }else{
        ref.model <- ordinal::clm(ref.term, data = dt.temp)
        add.model <- ordinal::clm(add.term, data = dt.temp)
        ano <- anova(ref.model, add.model)
      }
            
      names(ano)[length(names(ano))] <- "pvalue"
      pvalue <- ano$pvalue[2]
      coefficient <- add.model$coefficients["gene_count"]
      
      temp.out <- data.frame("enzid" = this.loci$enzid, "chr" = this.loci$chr, "start" = this.loci$start,
                             "end" = this.loci$end, "gsymbol" = this.loci$gsymbol, 
                             "type" = cnvtype, "coefficient" = coefficient, "pvalue" = pvalue, 
                             "sampleid" = this.loci$sample, 
                             "sampleclass" = paste(dt.temp$status[dt.temp$gene_count == 1], collapse = ","),
                             stringsAsFactors = F)
      
      if(permutation){
        temp.out[, sprintf("perm.pvalue.n%s", 1:nperm)] <- 1
        for(iperm in 1:nperm){
          dt.temp$outcome.perm <- perm.hg[, iperm]
          ref.perm.term <- sprintf("outcome.perm ~ %s", paste(covariates, collapse = " + "))
          add.perm.term <- sprintf("%s + %s", ref.perm.term, "gene_count")
          
          if(model == "lm"){
            ref.perm.model <- lm(ref.perm.term, dt.temp)
            add.perm.model <- lm(add.perm.term, dt.temp)
            ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
          }else if(model == "glm"){
            ref.perm.model <- glm(ref.perm.term, dt.temp, family = binomial(link = "logit"))
            add.perm.model <- glm(add.perm.term, dt.temp, family = binomial(link = "logit"))
            ano.perm <- anova(ref.perm.model, add.perm.model, test = "Chisq")
          }else{
            ref.perm.model <- ordinal::clm(ref.perm.term, data = dt.temp)
            add.perm.model <- ordinal::clm(add.perm.term, data = dt.temp)
            ano.perm <- anova(ref.perm.model, add.perm.model)
          }
          
          names(ano.perm)[length(names(ano.perm))] <- "pvalue"
          temp.out[, sprintf("perm.pvalue.n%s", iperm)] <- ano.perm$pvalue[2]
        }
      }
            
      dt.out <- rbind(dt.out, temp.out)
      
      current.count <- current.count + 1
      
      if((current.count%%10) == 0){
        print(sprintf("%s Loci testing - %s out of %s tests were done", cnvtype, current.count, all.count))
        flush.console()   
      }
    }
    message(sprintf("Loci testing done", all.count, all.count))
    
    message("Calculate FDR ...")
    dt.out.merge <- mergeLoci(dt.out, "pvalue")
    
    if(permutation){
      perm.pvalue <- c()
      for(i in 1:nperm){
        message(sprintf("Merging permution #%s", i))
        perm.merge <- mergeLoci(dt.out, sprintf("perm.pvalue.n%s", i))
        perm.pvalue <- c(perm.pvalue, perm.merge$pvalue)
      }
      
      dt.out.merge$permFDR <- 1
      for(i in 1:nrow(dt.out.merge)){
        actual <- sum(dt.out.merge$pvalue <= dt.out.merge$pvalue[i])/nrow(dt.out.merge)
        perm <- sum(perm.pvalue <= dt.out.merge$pvalue[i])/length(perm.pvalue)
        
        dt.out.merge$permFDR[i] <- perm/actual
        dt.out.merge$permFDR[i] <- ifelse(dt.out.merge$permFDR[i] > 1, 1, dt.out.merge$permFDR[i])
      }
      
      rownames(dt.out.merge) <- NULL
      dt.out.merge <- dt.out.merge[, c(2:4, 6, 7, 9:10, 1, 5, 8)]
      
      dt.out.merge <- dt.out.merge[order(dt.out.merge$pvalue), ]
      for(i in 1:nrow(dt.out.merge)){
        dt.out.merge$permFDR[i] <- min(dt.out.merge$permFDR[i:nrow(dt.out.merge)])
      }
    }
    
    final.out <- rbind(final.out, dt.out.merge)
  }
  return(final.out)
}

#'
#' This function is for internal use by other functions in the package.
#' @param test.table test result table
#' @param pvalue.column name of column contain p-values to be used
#' @keywords GSBurden 
#' @export
mergeLoci <- function(test.table, pvalue.column){
  test.table <- test.table[order(test.table[, pvalue.column]), ]
  test.out <- test.table[, c(names(test.table)[c(1:7,9)], pvalue.column)]
  names(test.out)[which(names(test.out) == pvalue.column)] <- "pvalue"
  
  mergable <- T
  
  while(mergable){
    test.out <- test.out[order(test.out$pvalue),]
    
    t.g <- GenomicRanges::GRanges(test.out$chr, IRanges::IRanges(test.out$start, test.out$end), "*")
    olap.t <- data.frame(IRanges::findOverlaps(t.g, t.g)) 
    olap.t <- olap.t[olap.t$queryHits != olap.t$subjectHits, ]
    olap.t$pair <- paste(pmin(olap.t$queryHits, olap.t$subjectHits), pmax(olap.t$queryHits, olap.t$subjectHits), sep=",")
    olap.t <- olap.t[!duplicated(olap.t$pair), ]
    
    mergable <- F
    if(nrow(olap.t) > 0){
      if(pvalue.column == "pvalue")
        message(sprintf("%s loci pair/s can potentially be merged", nrow(olap.t)))
      mergable <- T
      
      olap.t$pvalue <- pmin(test.out$pvalue[olap.t$queryHits], test.out$pvalue[olap.t$subjectHits])
      
      remove.test <- c()
      for(i in 1:nrow(olap.t)){
        
        olap.rec <- olap.t[i, ]
        subject1 <- unique(strsplit(as.character(test.out$sampleid[olap.rec$queryHits]), ",")[[1]])
        subject2 <- unique(strsplit(as.character(test.out$sampleid[olap.rec$subjectHits]), ",")[[1]])
        
        if(length(intersect(subject1, subject2))/length(union(subject1, subject2)) >= 0.75){
          enzid <- paste(unique(strsplit(paste(test.out$enzid[olap.rec$queryHits], test.out$enzid[olap.rec$subjectHits],sep=","),",")[[1]]), collapse = ",")
          chr <- test.out$chr[olap.rec$queryHits]
          start <- min(test.out$start[olap.rec$queryHits], test.out$start[olap.rec$subjectHits])
          end <- max(test.out$end[olap.rec$queryHits], test.out$end[olap.rec$subjectHits])
          gsymbol<- paste(unique(strsplit(paste(test.out$gsymbol[olap.rec$queryHits], test.out$gsymbol[olap.rec$subjectHits],sep=","),",")[[1]]), collapse = ",")
          type <- test.out$type[olap.rec$queryHits]
          sampleid <- paste(unique(strsplit(paste(test.out$sampleid[olap.rec$queryHits], test.out$sampleid[olap.rec$subjectHits],sep=","),",")[[1]]), collapse = ",")
          pvalue <- which.min(test.out$pvalue[c(olap.rec$queryHits, olap.rec$subjectHits)])
          coefficient <- test.out$coefficient[c(olap.rec$queryHits, olap.rec$subjectHits)][which.min(test.out$pvalue[c(olap.rec$queryHits, olap.rec$subjectHits)])]
          pvalue <- min(test.out$pvalue[c(olap.rec$queryHits, olap.rec$subjectHits)])
          
          merge.test <- data.frame(enzid, chr, start, end, gsymbol, type, coefficient, sampleid, pvalue)
          test.out <- rbind(test.out, merge.test)
          
          remove.test <- c(remove.test, olap.rec$queryHits, olap.rec$subjectHits)
        }
      }
      
      remove.test <- unique(remove.test)
      
      if(length(remove.test) == 0){
        mergable = F
      }else{
        test.out <- unique(test.out[-c(remove.test),])
      }
      
    }
    
  }
  
  message("Loci testing done!")
  return(test.out)
}
