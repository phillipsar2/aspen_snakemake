# make bam list
#rule bamlist:
#    input:
#       bams = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{bam}.dedup.bam", bam = BAM)
#    output:
#       bamlist = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{date}.bamlist.txt"
#    params:
#       path = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/"
#    shell:
#       "ls {params.path}*bam > {output}"

# (7) Call snps initially with bcftools to identify variable sites
# Call on each bam independently, as bcftools doesn't do joint calling anyway
# default only sites with max 250 reads considered at each positin, this is way above the max coverage
# -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
# -r output for only the given region
# --annotate FORMAT/AD,FORMAT/DP give allele and genotype depths
# # -b {input.bamlist} -r {params.chr} \
rule mpileup:
    input:
        ref = config["data"]["reference"]["genome"],
#        bamlist = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/aspen1222.bamlist.txt"
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{bam}.dedup.bam"
    output:
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{bam}.raw.vcf.gz"
    shell:
        """
        /global/scratch/users/arphillips/toolz/bcftools/bcftools mpileup -Ou \
        -f {input.ref} {input.bam} \
        -A --annotate FORMAT/AD,FORMAT/DP --threads 2 | \
        /global/scratch/users/arphillips/toolz/bcftools/bcftools call -mv -Oz -o {output}
        /global/scratch/users/arphillips/toolz/bcftools/bcftools index -t {output}
        """

# (8) Merge sample vcfs into larger vcf by region
rule vcflist:
    input:
        vcfs = expand("/global/scratch/users/arphillips/data/vcf/wgs_aspen.{bam}.raw.vcf.gz", bam = BAM)
    output:
        vcflist = "/global/scratch/users/arphillips/data/vcf/vcflist.raw.txt"
    params:
        path = "/global/scratch/users/arphillips/data/vcf/"
    shell:
        "ls {params.path}*.raw.vcf.gz > {output}"

rule bcftools_merge:
    input:
        list = "/global/scratch/users/arphillips/data/vcf/vcflist.raw.txt"
    output:
        "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.raw.merged.vcf.gz"
    params:
        region = "{region}"
    shell:
        "/global/scratch/users/arphillips/toolz/bcftools/bcftools merge -l {input} -r {params} -Oz -o {output}"

# (9) Extract SNPs from each vcf by region
rule get_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.raw.merged.vcf.gz"
    output:
         "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.snps.vcf.gz"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
        gatk SelectVariants \
        -R {input.ref} \
        -V {input.vcf} \
        -select-type SNP \
        -O {output}
        """


# (10) Filtering diagnostics - extract variant quality scores
# Roughly following suggestions in https://evodify.com/gatk-in-non-model-organism/
# Extract alternate base quality (QUAL), mapping quality (MQ), and total depth at the site (DP)
rule diagnostics:
    input:
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.snps.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{region}.table"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS -F QUAL -F DP -F MQ \
        -O {output}
        """

# (11) Hard filter SNPs
# https://gatk.broadinstitute.org/hc/en-us/articles/360035531112?id=23216#2
# https://gatk.broadinstitute.org/hc/en-us/articles/360037499012?id=3225

# Hard filter for mapping quality and base quality
rule filter_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{region}.snps.vcf.gz"
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.filtered.snps.vcf.gz"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
        gatk VariantFiltration \
        -V {input.vcf} \
        -filter \"QUAL < 30.0\" --filter-name \"QUAL30\" \
        -filter \"MQ < 30.0\" --filter-name \"MQ30\" \
        -O {output}
        """

# (12) Filter SNPs to only biallelic sites and exclude the sites that failed the hard filter
rule filter_nocall:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.filtered.snps.vcf.gz"
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.filtered.nocall.vcf.gz"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
        gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}
        """

# (13) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.filtered.nocall.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.{region}.filtered.nocall.table"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
        gatk VariantsToTable \
        -R {input.ref} \
        -V {input.vcf} \
        -F CHROM -F POS \
        -GF DP \
        -O {output}
        """

# (14) Fitlter by depth of each genotype at each site
rule filter_depth:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.filtered.nocall.vcf.gz",
        ref = config["data"]["reference"]["genome"]
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    output:
        dp = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.depth.{min_dp}dp{max_dp}.vcf.gz"
    params:
        min = "{min_dp}",
        max = "{max_dp}"
    shell:
        """
        gatk VariantFiltration \
        -R {input.ref} \
        -V {input.vcf} \
        -G-filter \"DP < {params.min} || DP > {params.max} \" \
        -G-filter-name \"DP_{params.min}-{params.max}\" \
        --set-filtered-genotype-to-no-call true -O {output.dp}
        """

# (15) Filter snps for genotype missingness (10%)
rule depth_nocall:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.depth.{min_dp}dp{max_dp}.vcf.gz",
    output:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf.gz",
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered true --max-nocall-fraction 0.1 -O {output}"

# (16) Exclude genotypes that don't meet quality thresholds
rule geno_filt:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf.gz",
        exc = "metadata/samples_to_drop_from_vcf.txt"
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.goodg.{min_dp}dp{max_dp}.vcf.gz"
    shell:
        "bcftools view -Oz -S ^{input.exc} {input.vcf} > {output}"
