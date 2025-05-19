# (7) make bam list
rule bamlist:
    input:
       bams = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{bam}.dedup.bam", bam = BAM)
    output:
       bamlist = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{date}.bamlist.txt"
    params:
       path = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/"
    shell:
       "ls {params.path}*bam > {output}"

# (8) Call snps initially with bcftools to identify variable sites
# default only sites with max 250 reads considered at each positin, this is way above the max coverage
# -v option asks to output variant sites only (this is sufficient for the analyses we want to run)
# -r output for only the given region
# --annotate FORMAT/AD,FORMAT/DP give allele and genotype depths
rule mpileup:
    input:
        ref = config["data"]["reference"]["genome"],
#        bamlist = expand("/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{date}.bamlist.txt", date = DATE) 
        bamlist = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/bamlist.05062025.txt"
    output:
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{chr}.raw.vcf.gz"
    params:
        chr = "{chr}"
#    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/bcftools.yaml"
#    benchmark:
#         "/global/scratch/users/arphillips/benchmarks/{chr}.mpileup.benchmark.txt"
    shell:
        """
        /global/scratch/users/arphillips/toolz/bcftools/bcftools mpileup -Ou -f {input.ref} -b {input.bamlist} -r {params.chr} -A --annotate FORMAT/AD,FORMAT/DP --threads 10 | /global/scratch/users/arphillips/toolz/bcftools/bcftools call -mv -Oz -o {output.vcf}
        /global/scratch/users/arphillips/toolz/bcftools/bcftools index -t {output}
        """

# (9) Extract SNPs from each vcf
rule get_snps:
    input:
        ref = config["data"]["reference"]["genome"],
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{chr}.raw.vcf.gz"
    output:
         "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{chr}.snps.vcf.gz"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
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
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{chr}.snps.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/reports/filtering/wgs_aspen.{chr}.table"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
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
        vcf = "/global/scratch/users/arphillips/data/vcf/wgs_aspen.{chr}.snps.vcf.gz"
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.filtered.snps.vcf"
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
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.filtered.snps.vcf"
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.filtered.nocall.vcf"
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        """
        gatk SelectVariants -V {input.vcf} --exclude-filtered true  --restrict-alleles-to BIALLELIC -O {output}
        """

# (13) Extract genotype depth across samples to determine DP cutoff
rule depth:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
    output:
        "/global/scratch/users/arphillips/reports/filtering/depth/wgs_aspen.{chr}.filtered.nocall.table"
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
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.filtered.nocall.vcf",
        ref = config["data"]["reference"]["genome"]
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    output:
        dp = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.depth.{min_dp}dp{max_dp}.vcf"
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
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.depth.{min_dp}dp{max_dp}.vcf",
    output:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.nocall.{min_dp}dp{max_dp}.vcf",
#    conda: "/global/home/users/arphillips/.conda/envs/gatk"
    shell:
        "gatk SelectVariants -V {input} --exclude-filtered true --max-nocall-fraction 0.1 -O {output}"

