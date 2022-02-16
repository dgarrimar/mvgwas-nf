#!/usr/bin/env Rscript

## 1. Load libraries and arguments

library(optparse)
library(data.table)
library(car)

option_list <- list(
    make_option(c("-p", "--phenotypes"), type = "character", 
                help = "Phenotype file (indId, [phenotypes])", metavar = "FILE"),
    make_option(c("-c", "--covariates"), type = "character",
                help = "Covariate file (indId, [covariates])", 
                metavar = "FILE"),
    make_option(c("-g", "--genotypes"), type = "character", 
                help = "Genotype file (indexed VCF)", metavar = "FILE"),
    make_option(c("-m", "--min_value"), type = "numeric", default = 0,
                help = "Minimum value of the phenotype [default %default]", metavar = "NUMERIC"),
    make_option("--out_pheno", type = "character",
                help = "Output pre-processed phenotypes", metavar = "FILE"),
    make_option("--out_cov", type = "character",
                help = "Output pre-processed covariates", metavar = "FILE"),
    make_option(c("-s", "--seed"), type = "numeric", help = "Set seed for random processes",
                metavar = "NUMERIC", default = 123),
    make_option(c("-v", "--verbose"), action = "store_true", 
                help = "[default %default]", 
                default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

## 2. Input files

pheno.f <- opt$phenotypes   
cov.f <- opt$covariates
geno.f <- opt$genotypes
min.val <- opt$min_value
out_pheno.f <- opt$out_pheno
out_cov.f <- opt$out_cov

set.seed(opt$seed)

if ( is.null(pheno.f) || is.null(geno.f) || is.null(cov.f) || 
     is.null(out_pheno.f) || is.null(out_cov.f) ){
    print_help(opt_parser)
    stop("Missing/not found I/O files")
}

pheno.df <- fread(pheno.f, header = TRUE, data.table = FALSE, sep = "\t")
cov.df <- fread(cov.f, header = TRUE, data.table = FALSE, sep = "\t")
colnames(pheno.df)[1] <- colnames(cov.df)[1] <- "ID" 
rownames(pheno.df) <- as.character(pheno.df$ID)
rownames(cov.df) <- as.character(cov.df$ID)
pheno.df$ID <- cov.df$ID <- NULL

geno.df <- fread(cmd = sprintf("zcat %s | head -n 10000", geno.f), skip = "#CHROM", nrows = 1,
                 data.table = FALSE, header = TRUE, check.names = FALSE)

subset.ids <- Reduce(intersect, list(colnames(geno.df),
                                     rownames(pheno.df), 
                                     rownames(cov.df))) 
if (length(subset.ids) == 0) {
    stop("No common samples between genotype, covariate and transcript files")
}

pheno.df <- pheno.df[subset.ids, , drop = FALSE]
cov.df <- cov.df[subset.ids, , drop = FALSE] 

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
    warning(sprintf("Phenotype(s) %s is(are) removed due to value < %s in all common individuals", 
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

## 5. Write pre-processed phenotypes and covariates

fwrite(pheno.df, file = out_pheno.f, row.names = TRUE, quote = FALSE, sep = "\t")
fwrite(cov.df, file = out_cov.f, row.names = TRUE, quote = FALSE, sep = "\t")

#### END

