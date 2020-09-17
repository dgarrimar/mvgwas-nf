#!/usr/bin/env Rscript

#### Prepare transcript expression file 

## 1. Load libraries and arguments

library(optparse)
library(data.table)
library(seqminer)
library(mlm)
  
option_list <- list(
    make_option(c("-p", "--phenotypes"), type = "character", 
                help = "Phenotype file (indId, [phenotypes])", metavar = "FILE"),
    make_option(c("-c", "--covariates"), type = "character",
                help = "Covariate file (indId, [covariates])", 
                metavar = "FILE"),
    make_option(c("-g", "--genotypes"), type = "character", 
                help = "Genotype file (indexed VCF)", metavar = "FILE"),
    make_option(c("-r", "--region"), type = "character",
                help = "Genomic region", metavar = "CHARACTER"),
    make_option(c("-m", "--min_value"), type = "numeric", default = 0,
                help = "Minimum value of the phenotype [default %default]", metavar = "NUMERIC"),
    make_option(c("-n",	"--min_nb_ind_geno"), type = "numeric", default = 10,
       	       	help = "Minimum number of individuals per genotype group"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output summary stats", metavar = "FILE"),
    make_option(c("-s", "--seed"), type = "numeric", help = "Set seed for random processes",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "[default %default]", 
                default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

tabix.read.table.nochecknames <- function (tabixFile, tabixRange, col.names = TRUE, stringsAsFactors = FALSE) 
{
  stopifnot(seqminer:::local.file.exists(tabixFile))
  stopifnot(all(isTabixRange(tabixRange)))
  tabixFile <- path.expand(tabixFile)
  storage.mode(tabixFile) <- "character"
  storage.mode(tabixRange) <- "character"
  header <- .Call("readTabixHeader", tabixFile, PACKAGE = "seqminer")
  body <- .Call("readTabixByRange", tabixFile, tabixRange, 
                PACKAGE = "seqminer")
  body <- do.call(rbind, strsplit(body, "\t"))
  body <- as.data.frame(body, stringsAsFactors = FALSE)
  if (ncol(body) > 0) {
    for (i in 1:ncol(body)) {
      body[, i] <- utils::type.convert(body[, i], as.is = !stringsAsFactors)
    }
    num.col <- ncol(body)
    header <- header[nchar(header) > 0]
    if (length(header) == 0 || !col.names) {
      colNames <- paste0("V", 1L:num.col)
    }
    else {
      hdrLine <- header[length(header)]
      hdrLine <- sub("^#", "", hdrLine)
      # colNames <- make.names(strsplit(hdrLine, "\t")[[1]])
      colNames <- strsplit(hdrLine, "\t")[[1]]
      if (length(colNames) > ncol(body)) {
        colNames <- colNames[1:ncol(body)]
      }
      else if (length(colNames) < ncol(body)) {
        tmpNames <- paste0("V", 1L:num.col)
        tmpNames[1:length(colNames)] <- colNames
        colNames <- tmpNames
      }
    }
    colnames(body) <- colNames
  }
  body
}

## 2. Input files: transcript expression, sample groups, gene location
  
pheno.f <- opt$phenotypes   
cov.f <- opt$covariates
geno.f <- opt$genotypes
out.f <- opt$output
region <- opt$region
min.val <- opt$min_value

if ( is.null(pheno.f) || is.null(geno.f) || is.null(cov.f) || is.null(region) || is.null(out.f) ){
    print_help(opt_parser)
    stop("Missing/not found input files", call.= FALSE)
}

pheno.df <- read.table(pheno.f, header = TRUE, as.is = TRUE, sep = "\t")
cov.df <- read.table(cov.f, header = TRUE, as.is = TRUE, sep = "\t")
colnames(pheno.df)[1] <- colnames(cov.df)[1] <- "ID" 
rownames(pheno.df) <- pheno.df$ID 
rownames(cov.df) <- cov.df$ID
pheno.df$ID <- cov.df$ID <- NULL

geno.df <- tabix.read.table.nochecknames(geno.f, region)

subset.ids <- Reduce(intersect, list(colnames(geno.df),
                                     rownames(pheno.df), rownames(cov.df))) 
if (length(subset.ids) == 0) {
  stop("No common samples between genotype, covariate and transcript files.")
}

pheno.df <- pheno.df[subset.ids, , drop = FALSE]
cov.df <- cov.df[subset.ids, , drop = FALSE] 

colnames(geno.df)[1:5] <- c("chr", "pos", "variant", "REF", "ALT")
geno.df <- geno.df[,colnames(geno.df) %in% c("chr", "pos", "variant", "REF", "ALT", subset.ids)]
geno.df[,-c(1:5)] <- apply(geno.df[,-c(1:5)], 2, function(x){y <- gsub("([012.]/[012.]):.*","\\1", x)
sapply(y, function(z){switch(z,
                             # Unphased
                             "./." = NA,
                             "0/0" = "0",
                             "1/0" = "1",
                             "0/1" = "1",
                             "1/1" = "2",
                             # Phased
                             ".|." = NA,
                             "0|0" = "0",
                             "1|0" = "1",
                             "0|1" = "1",
                             "1|1" = "2",
                             # Haploid (X,Y,MT)
                             "." = NA,
                             "0" = "0",
                             "1" = "1"
                             )})
})

## 3. Pre-process phenotypes

all.na <- unlist(lapply(pheno.df, function(x){all(is.na(x))}))
if (opt$verbose && sum(all.na) > 0){
  warning(sprintf("Phenotypes with NA values for all individuals were removed: %s",
                  paste(names(which(all.na)), collapse = ", ")))
}

pheno.na <- apply(pheno.df, 1, function(x){any(is.na(x))})
if(sum(pheno.na) / nrow(pheno.df) > 0.05){
  stop("More than 5% of the individuals contain NA values for at least one phenotype")  
}

allzero <- apply(pheno.df, 2, function(x){ all(x < min.val) })
if (any(allzero)){
  out <- names(which(allzero))
  pheno.df <- pheno.df[, !colnames(pheno.df)%in%out]
  warning(sprintf("Phenotype(s) %s is(are) removed due to value < %s in all common individuals.", 
                  paste(out, collapse = ", "), min.val) )
}

## 4. Pre-process covariates

all.na <- unlist(lapply(cov.df, function(x){all(is.na(x))}))
if (opt$verbose && sum(all.na) > 0){
  warning(sprintf("Covariates with NA values for all individuals were removed: %s",
                  paste(names(which(all.na)), collapse = ", ")))
}

cov.df <- cov.df[, !all.na, drop = FALSE]
for(i in 1:ncol(cov.df)){
  typ <- class(cov.df[, i])
  if(typ == "character"){
    cov.df[, i] <- as.factor(cov.df[, i])
  } else if (typ == "numeric" || typ == "integer"){
    next
  } else {
    stop ("Covariates should be either 'numeric' or 'character'")
  }
}

types <- unlist(lapply(cov.df, class))
if (opt$verbose) {
  message("Covariate types:\n", 
          paste(names(types), types, sep=": ", collapse = ", "))
}

cov.na <- apply(cov.df, 1, function(x){any(is.na(x))})
if(sum(cov.na) / nrow(cov.df) > 0.05){
  stop("More than 5% of the individuals contain NA values for at least one covariate")  
}

multiclass <- apply(cov.df, 2, function(x){length(table(x)) > 1})
cov.df <- cov.df[, multiclass, drop = FALSE]

if (opt$verbose && any(!multiclass)){
  message("\t", "Covariates removed due to only one value: ",
          paste(names(multiclass)[!multiclass], collapse = ", "))
}

if (ncol(cov.df) > 1){
  vifs <- car::vif(lm(pheno.df[, 1] ~ ., data = cov.df))
  if (opt$verbose){
    message("\t", "Covariates VIF - ", 
            paste(names(vifs), round(vifs, 2), sep = ": ", collapse = ", "))
  }
  if (any(vifs > 5)){
    warning("Check multicollinearity. VIF > 5 for some covariates:", "\n",
            paste(names(vifs), round(vifs, 2), sep = ": ", collapse = ", "))
  }  
}  

## 5. Filter SNPs

check.genotype <- function(geno.df, min.nb.ind.geno = 10)
{
  apply(geno.df, 1, function(geno.snp){
    if(sum(is.na(geno.snp)) > 0.05*length(as.numeric(geno.snp))){
      return("Missing genotype in > 5% individuals")
    }
    geno.snp.t <- table(geno.snp[!is.na(geno.snp)])
    if (length(geno.snp.t) < 3) {
      return("One or two genotype groups")                    
    }
    if (sum(geno.snp.t >= min.nb.ind.geno) < length(geno.snp.t)) {
      return(sprintf("Not all the groups with >%s samples", min.nb.ind.geno))        
    }
    return("PASS")
  })
}

snps.to.keep <- check.genotype(geno.df[, subset.ids], min.nb.ind.geno = opt$min_nb_ind_geno)

if(opt$verbose){
  snps.to.keep.t <- table(snps.to.keep)
  message("\t", paste(names(snps.to.keep.t), snps.to.keep.t, sep = ": ", collapse=", "))
}

if(any(snps.to.keep == "PASS")){
  geno.df <- geno.df[snps.to.keep == "PASS", ]
}

out.df <- c()
for (p in geno.df$pos){
  snp <- subset(geno.df, pos == p)
  rec <- snp[, !colnames(snp)%in%subset.ids]
  snp <- as.numeric(snp[, subset.ids])
  
  mvfit <- tryCatch(mlm(as.matrix(pheno.df) ~ ., data = data.frame(cov.df, "GT" = snp)),
                    error = function(e) NULL)
  if(is.null(mvfit)){
    warning(sprintf("SNP %s skipped",  subset(geno.df, pos == p)$variant))
    next
  }
  out.df <- rbind(out.df, c(t(rec), mvfit$aov.tab[nrow(mvfit$aov.tab)-1,c(5,6)]))
}

## 6. Write output

write.table(out.df, file = out.f, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

#### END

