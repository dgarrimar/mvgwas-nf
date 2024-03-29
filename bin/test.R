#!/usr/bin/env Rscript

## 1. Load libraries and arguments, define functions

library(optparse)
library(data.table)
library(seqminer)
library(manta)
  
option_list <- list(
    make_option(c("-p", "--phenotypes"), type = "character", 
                help = "Pre-processed phenotype file (as generated by 'preprocess.R')", 
                metavar = "FILE"),
    make_option(c("-c", "--covariates"), type = "character",
                help = "Pre-processed covariate file (as generated by 'preprocess.R')", 
                metavar = "FILE"),
    make_option(c("-g", "--genotypes"), type = "character", 
                help = "Genotype file (indexed VCF)", metavar = "FILE"),
    make_option(c("-r", "--region"), type = "character",
                help = "Genomic region", metavar = "CHARACTER"),
    make_option(c("-t", "--transform"), type = "character", default = 'none',
                help = "Phenotype transformation: 'none', 'log', 'sqrt' [default %default]",
                metavar = "CHARACTER"),
    make_option(c("-i", "--interaction"), type = "character", default = 'none',
                help = "Test for interaction with a covariate [default %default]",
                metavar = "CHARACTER"),
    make_option(c("-n",	"--min_nb_ind_geno"), type = "numeric", default = 10,
       	       	help = "Minimum number of individuals per genotype group [default %default]",
                metavar = "NUMERIC"),
    make_option(c("-o", "--output"), type = "character",
                help = "Output summary stats", metavar = "FILE"),
    make_option(c("-s", "--seed"), type = "numeric", help = "Set seed for random processes [default %default]",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "[default %default]", 
                default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

tabix.read.table.nochecknames <- function (tabixFile, tabixRange, 
                                           col.names = TRUE, 
                                           stringsAsFactors = FALSE) {
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
        } else {
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

check.genotype <- function(geno.df, min.nb.ind.geno = 10) {
    apply(geno.df, 1, function(geno.snp) {
        if (sum(is.na(geno.snp)) > 0.05*length(as.numeric(geno.snp))) {
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

## 2. Input files
  
pheno.f <- opt$phenotypes   
cov.f <- opt$covariates
geno.f <- opt$genotypes
out.f <- opt$output
region <- opt$region

set.seed(opt$seed)

if (is.null(pheno.f) || is.null(geno.f) || is.null(cov.f) || is.null(region) || 
    is.null(out.f)) {
    print_help(opt_parser)
    stop("Missing/not found I/O files", call.= FALSE)
}

pheno.df <- data.frame(fread(pheno.f, header = TRUE, sep = "\t"), row.names = 1)
cov.df <- data.frame(fread(cov.f, header = TRUE, sep = "\t"), row.names = 1)
geno.df <- tabix.read.table.nochecknames(geno.f, region)

subset.ids <- rownames(pheno.df)

cn <- c("chr", "pos", "variant", "REF", "ALT") 
colnames(geno.df)[1:5] <- cn
geno.df <- geno.df[, colnames(geno.df) %in% c(cn, subset.ids)]
geno.df[, -c(1:5)] <- apply(geno.df[, -c(1:5)], 2, function(x){
    y <- gsub("([012.]/[012.]):.*","\\1", x)
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

## 3. Filter SNPs
snps.to.keep <- check.genotype(geno.df[, subset.ids], min.nb.ind.geno = opt$min_nb_ind_geno)

if (opt$verbose) {
    snps.to.keep.t <- table(snps.to.keep)
    message("\t", paste(names(snps.to.keep.t), snps.to.keep.t, sep = ": ", collapse=", "))
}

## 4. Test & write output
if (any(snps.to.keep == "PASS")) {
    geno.df <- geno.df[snps.to.keep == "PASS", ]
    out.df <- c()
    Y <- as.matrix(pheno.df)
    if (opt$interaction == "none") {
        for (p in geno.df$pos) {
            snp <- subset(geno.df, pos == p)
            rec <- snp[, !colnames(snp)%in%subset.ids]
            snp <- as.numeric(snp[, subset.ids])
              
            mvfit <- tryCatch(manta(Y ~ ., data = data.frame(cov.df, "GT" = snp), type = "I", subset = "GT", transform = opt$transform),
                                error = function(e) NULL)
            if (is.null(mvfit)) {
                warning(sprintf("SNP %s skipped",  subset(geno.df, pos == p)$variant))
                next
            }
            out.df <- rbind(out.df, c(t(rec), mvfit$aov.tab[1, 4:6]))
        }
    } else {
        INT <- paste0(opt$interaction, ":GT")
        for (p in geno.df$pos) {
            snp <- subset(geno.df, pos == p)
            rec <- snp[, !colnames(snp)%in%subset.ids]
            snp <- as.numeric(snp[, subset.ids])
            Data <- data.frame(cov.df, "GT" = snp)
            fm <- as.formula(paste("Y ~", paste0(c(colnames(Data), INT), collapse = "+")))
            mvfit <- tryCatch(manta(fm,  data = data.frame(cov.df, "GT" = snp), type = "II", transform = opt$transform, 
                                    subset = c(opt$interaction, "GT", INT)),
                              error = function(e) NULL)
            if (is.null(mvfit)) {
                warning(sprintf("SNP %s skipped",  subset(geno.df, pos == p)$variant))
                next
            }
            out.df <- rbind(out.df, c(t(rec), mvfit$aov.tab[1:3, 4:6]))
        }
    }
    fwrite(out.df, file = out.f, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
} 

#### END

