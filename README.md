# mvgwas-nf

[![nextflow](https://img.shields.io/badge/nextflow-%E2%89%A520.04.1-blue.svg)](http://nextflow.io)
[![CI-checks](https://github.com/guigolab/sqtlseeker2-nf/actions/workflows/ci.yaml/badge.svg)](https://github.com/guigolab/sqtlseeker2-nf/actions/workflows/ci.yaml)

A pipeline for multi-trait genome-wide association studies (GWAS) using [MANTA](https://github.com/dgarrimar/manta).

The pipeline performs the following analysis steps:

* Split genotype file 
* Preprocess phenotype and covariate data
* Test for association between phenotypes and genetic variants
* Collect summary statistics

The pipeline uses [Nextflow](http://www.nextflow.io) as the execution backend. Please check [Nextflow documentation](http://www.nextflow.io/docs/latest/index.html) for more information.

## Requirements

- Unix-like operating system (Linux, MacOS, etc.)
- Java 8 or later 
- [Docker](https://www.docker.com/) (v1.10.0 or later) or [Singularity](http://singularity.lbl.gov) (v2.5.0 or later)

## Quickstart (~2 min)

1. Install Nextflow:
    ```
    curl -fsSL get.nextflow.io | bash
    ```

2. Make a test run:
    ```
    nextflow run dgarrimar/mvgwas-nf -with-docker
    ```

**Notes**: move the `nextflow` executable to a directory in your `$PATH`. Set `-with-singularity` to use Singularity instead of Docker. 

(*) Alternatively you can clone this repository:
```
git clone https://github.com/dgarrimar/mvgwas-nf
cd mvgwas-nf
nextflow run mvgwas.nf -with-docker
```

## Pipeline usage

Launching the pipeline with the `--help` parameter shows the help message:

```
nextflow run mvgwas.nf --help
```

```
N E X T F L O W  ~  version 20.04.1
Launching `mvgwas.nf` [amazing_roentgen] - revision: 56125073b7

mvgwas-nf: A pipeline for multivariate Genome-Wide Association Studies
==============================================================================================
Performs multi-trait GWAS using using MANTA (https://github.com/dgarrimar/manta)

Usage:
nextflow run mvgwas.nf [options]

Parameters:
--pheno PHENOTYPES          phenotype file
--geno GENOTYPES            indexed genotype VCF file
--cov COVARIATES            covariate file
--l VARIANTS/CHUNK          variants tested per chunk (default: 10000)
--t TRANSFOMATION           phenotype transformation: none, sqrt, log (default: none)
--i INTERACTION             test for interaction with a covariate: none, <covariate> (default: none)
--ng INDIVIDUALS/GENOTYPE   minimum number of individuals per genotype group (default: 10)
--dir DIRECTORY             output directory (default: result)
--out OUTPUT                output file (default: mvgwas.tsv)
```

## Input files and format

`mvgwas-nf` requires the following input files:

* **Genotypes.** 
[bgzip](http://www.htslib.org/doc/bgzip.html)-compressed and indexed [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) genotype file.

* **Phenotypes.**
Tab-separated file with phenotype measurements (quantitative) for each sample (i.e. *n* samples x *q* phenotypes).
The first column should contain sample IDs. Columns should be named.

* **Covariates.**
Tab-separated file with covariate measurements (quantitative or categorical) for each sample (i.e. *n* samples x *k* covariates). 
The first column should contain sample IDs. Columns should be named. 

Example [data](data) is available for the test run.

## Pipeline results

An output text file containing the multi-trait GWAS summary statistics (default: `./result/mvgwas.tsv`), with the following information:

* `chr`: chromosome
* `pos`: position
* `snp`: variant ID
* `REF`: reference allele
* `ALT`: alternative allele
* `F`: pseudo-F statistic
* `r2`: fraction of variance explained by the variant
* `pv`: P-value

The output folder and file names can be modified with the `--dir` and `--out` parameters, respectively.
