/*
 * Multivariate Genome Wide Association Studies (MVGWAS) 
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.pheno = 'data/phenotypes.tsv'
params.geno = 'data/genotypes.vcf.gz'
params.cov = 'data/covariates.tsv'
params.l = 10000
params.dir = 'result'
params.out = 'mvgwas.tsv'
params.help = false

/*
 *  Print usage and help
 */

if (params.help) {
  log.info ''
  log.info 'Multivariate Genome-Wide Association Studies (MVGWAS)'
  log.info '======================================================================='
  log.info 'Performs multivariate GWAS given a set of phenotypes and genotypes'
  log.info ''
  log.info 'Usage: '
  log.info '    nextflow run mvgwas.nf [options]'
  log.info ''
  log.info 'Parameters:'
  log.info ' --pheno PHENOTYPES          phenotype file (default: phenotypes.tsv)'
  log.info ' --geno GENOTYPES            indexed genotype VCF file (default: genotypes.vcf.gz)'
  log.info ' --cov COVARIATES            covariate file (default: covariates.tsv)'
  log.info ' --l VARIANTS/CHUNK          variants tested per chunk (default: 10000)'
  log.info ' --dir DIRECTORY             output directory (default: result)'
  log.info ' --out OUTPUT                output file (default: mvgwas.tsv)'
  log.info ''
  exit(1)
}

/*
 *  Print parameter selection
 */

log.info ''
log.info 'Parameters'
log.info '------------------'
log.info "Phenotype data               : ${params.pheno}"
log.info "Genotype data                : ${params.geno}"
log.info "Covariates                   : ${params.cov}"
log.info "Variants/chunk               : ${params.l}"
log.info "Output directory             : ${params.dir}"
log.info "Output file                  : ${params.out}"
log.info ''


/*
 *  Preprocess VCF
 */

process split {
 
    input:
   
    file(vcf) from file(params.geno)
    file(index) from file("${params.geno}.tbi")    

    output:
    
    file("chunk*") into chunks_ch    

    script:
    """
    bcftools query -f '%CHROM\t%POS\n' $vcf > positions
    split -d -a 10 -l ${params.l} positions chunk
    """
}

/*
 *  GWAS: mlm testing
 */

process mvgwas {

    input:

    file pheno from file(params.pheno)
    file cov from file(params.cov)
    file(vcf) from file(params.geno)
    file(index) from file("${params.geno}.tbi")
    each file (chunk) from chunks_ch

    output:

    file('sstats.txt') into sstats_ch

    script:
    """
    if [[ \$(cut -f1 $chunk | sort | uniq -c | wc -l) == 2 ]]; then
        k=1
        cut -f1 $chunk | sort | uniq | while read chr; do
        region=\$(paste <(grep "^\$chr" $chunk | head -1) <(grep "^\$chr" $chunk | tail -1 | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\$k.txt --verbose
        ((k++))
    done
    cat sstats.*.txt > sstats.txt
    else
        region=\$(paste <(head -1 $chunk) <(tail -1 $chunk | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.txt --verbose
    fi
    """
}

sstats_ch.collectFile(name: "${params.out}").set{pub_ch}

/*
 * Summary stats
 */

process end {

   publishDir "${params.dir}"     

   input:
   file(out) from pub_ch

   output:
   file(out) into end_ch

   script:
   """
   sed -i "1 s/^/chr\tpos\tsnp\tREF\tALT\tr2\tpv\\n/" ${out}
   """
}

