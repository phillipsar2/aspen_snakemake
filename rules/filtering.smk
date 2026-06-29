# (17) Index gvcfs
rule index_vcfs:
    input:
        "data/vcf/gatk/genotyped/{geno}_p{geno_ploidy}.g.vcf.gz"
    output:
        "data/vcf/gatk/genotyped/{geno}_p{geno_ploidy}.g.vcf.gz.tbi"
    shell:
        "bcftools index -t {input}"

# (18) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
rule diagnostics:
    input:
        vcf = "data/vcf/gatk/genotyped/{geno}_p{geno_ploidy}.g.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "reports/filtering/{geno}_p{geno_ploidy}.table"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F QD -F DP -F MQ -F ReadPosRankSum \
        -O {output}
        """

##########

# (20) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "data/raw/vcf_bpres/{cov}/all.AG.{cov}.{chr}.raw.snps.vcf.gz"
    output:
        "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.snps.vcf"
    run:
        shell("gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 30.0\" --filter-name \"MQ30\" \
        -O {output}")


# (21) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.snps.vcf"
    output:
        "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf"
    run: 
        shell("gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}")
        

# (22) Exclude clones
# Clones were identified by looking at the kinship matrix using all genotypes.
rule exclude_clones:
    input:
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"],
        clones = "clones.list"
    output:
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.vcf",
    shell:
        """
        vcftools --vcf {input.vcf} --remove {input.clones} --stdout > {output.vcf}
        """ 

# (23) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
#        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.{cov}.{chr}.filtered.nocall.vcf",
        vcf = "data/processed/filtered_snps_bpres/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    output:
#        dp = "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.table"
        dp = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.table"
    run:
        shell("gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output.dp}")

# (24) Filter by genotype depth and missingness
# p = probability of the given read depth
# miss = percent missing data threshold at a site
# min = minimum depth for a genotype at a site
rule filter_depth:
    input:
#        vcf = "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.table"
        vcf = "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.table"
    output:
        ## Output used to select sites in ANGSD
#        "reports/filtering/depth/{cov}/all.AG.{cov}.{chr}.filtered.nocall.0.99_0.2.txt"
        "reports/filtering/depth/{cov}/all.AG.noclones.{cov}.{chr}.filtered.nocall.0.99_0.2.txt"
    params:
#        p = "{p}",
#        miss = "{miss}"
        p = "0.99",
        miss = "0.2",
        min = "1"
    shell:
        "Rscript scripts/genoDPfilter.R {input.vcf} -q {params.p} -m {params.miss} --min {params.min}"
