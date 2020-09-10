/*
 * Multivariate Genome Wide Association Studies (MVGWAS) 
 * Diego Garrido Martín 
 */

/*
 *  Define parameters
 */

// General params
params.pheno = null
params.geno = null
params.cov = null
params.l = 500
params.ng = 10
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
  log.info ' --ng INDIVIDUALS/GENOTYPE   minimum number of individuals per genotype group (default: 10)'        
  log.info ' --dir DIRECTORY             output directory (default: result)'
  log.info ' --out OUTPUT                output file (default: mvgwas.tsv)'
  log.info ''
  exit(1)
}

/*
 * Check mandatory parameters
 */

if (!params.pheno) {
    exit 1, "Phenotype file not specified."
} else if (!params.geno){
    params.help
    exit 1, "Genotype not specified."
} else if (!params.cov){
    exit 1, "Covariate file not specified."
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
log.info "Individuals/genotype         : ${params.ng}" 
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

    file('sstats.*.txt') into sstats_ch

    script:
    """
    chunknb=\$(basename $chunk | sed 's/chunk//')
    if [[ \$(cut -f1 $chunk | sort | uniq -c | wc -l) == 2 ]]; then
        k=1
        cut -f1 $chunk | sort | uniq | while read chr; do
        region=\$(paste <(grep "^\$chr" $chunk | head -1) <(grep "^\$chr" $chunk | tail -1 | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\$k.tmp --min_nb_ind_geno ${params.ng} --verbose
        ((k++))
    done
    cat sstats.*.tmp > sstats.\${chunknb}.txt
    else
        region=\$(paste <(head -1 $chunk) <(tail -1 $chunk | cut -f2) | sed 's/\t/:/' | sed 's/\t/-/')
        test.R --phenotypes $pheno --covariates $cov --genotypes $vcf --region "\$region" --output sstats.\${chunknb}.txt --min_nb_ind_geno ${params.ng} --verbose
    fi
    """
}

sstats_ch.collectFile(name: "${params.out}", sort: { it.name }).set{pub_ch}

/*
 * Summary stats
 */

process end {

   publishDir "${params.dir}", mode: 'copy'     

   input:
   file(out) from pub_ch

   output:
   file(out) into end_ch

   script:
   """
   sed -i "1 s/^/chr\tpos\tsnp\tREF\tALT\tr2\tpv\\n/" ${out}
   """
}

