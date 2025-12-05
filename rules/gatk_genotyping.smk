# (1) Merge filtered vcf files
#rule merge_filt:
#    input:
#        vcfs = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf.gz", region = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
#    output:
#        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz"
#    params:
#        list = "/global/scratch/users/arphillips/data/processed/filtered_snps/vcflist.{min_dp}dp{max_dp}.txt",
#        path = "/global/scratch/users/arphillips/data/processed/filtered_snps/",
#        suffix = "wgs_aspen.Chr*.nocall.{min_dp}dp{max_dp}.vcf.gz" 
#    shell:
#        """
#        ls {params.path}*{params.suffix} > {params.list} 
#        /global/scratch/users/arphillips/toolz/bcftools/bcftools concat -f {params.list} --threads 10 -Oz -o {output}
#        """

# (2) Subset vcf to a single sample
## had to sort and index the vcf
#rule select_sample:
#    input:
#       vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.sorted.vcf.gz",
#        ref = config["data"]["reference"]["genome"]
#    output: 
#        temp("/global/scratch/users/arphillips/data/processed/filtered_snps/{geno}.10dp90.vcf.gz")
#    params:
#        geno = "{geno}"
#    shell:
#        """
#        gatk SelectVariants \
#        -R {input.ref} \
#        -V {input.vcf} \
#        --sample-name {params.geno} \
#        -O {output}
#        """

# (1) Haplotype caller
rule haplotype_caller:
    input:
        ref = config["data"]["reference"]["genome"], 
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{geno}.dedup.bam" 
    output:
        temp("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.{region}.g.vcf.gz")
    params:
        region = "{region}"
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        --input {input.bam} \
        --output {output} \
        --reference {input.ref} \
        -ERC BP_RESOLUTION \
        -G StandardAnnotation \
        -L {params.region} 
        """

# (2) Individually genotype with correct ploidy level specified
## Merge vcfs for each genotype then do genotyping
rule indv_geno:
    input:
        vcf = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.{region}.g.vcf.gz",
#        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/{geno}.10dp90.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
       temp("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.{region}.ploidy{geno_ploidy}.g.vcf.gz")
    params:
        ploidy = "{geno_ploidy}",
        tmpdir =  "/global/scratch/users/arphillips/tmp/joint_geno/{geno}.{region}",
        pre = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}",
        list = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.list"
    shell:
        """
        mkdir -p {params.tmpdir}
        ls {params.pre}*.g.vcf.gz > {params.list}
        gatk MergeVcfs \
        -I {params.list} \
        -O {params.pre}.g.vcf.gz

        gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {input.pre}.g.vcf.gz \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --tmp-dir {params.tmpdir} \
        --sample-ploidy {params.ploidy} \
        -O {output}
        rm -rf {params.tmpdir}
        """

#(3) Merge samples into large VCF
rule mergevcfs:
    input:
        vcfs = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.{region}.ploidy{geno_ploidy}.vcf.gz", zip,  geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY, region = CHR)
    output:
        "/global/scratch/users/arphillips/data/vcf/gatk/called/wgs_aspen.all.genos.{region}.vcf.gz"
    params:
        dir = "/global/scratch/users/arphillips/data/vcf/gatk/called/",
        list = "/global/scratch/users/arphillips/data/vcf/gatk/called/genotyped.vcf.list" 
    shell:
        """
        ls {params.dir}.*dp*ploidy*.vcf.gz > {params.list}
        gatk MergeVcfs \
        -I {params.list} \
        -O {output}
        """


# Genotyping quickly with GATK
# Gentoyping assumbing all are diploid or triploid to see how results change
rule joint_geno:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.nocall.{min_dp}dp{max_dp}.vcf",
#        dir = directory("data/interm/combined_database_bpres/{REF}/{interval}"),
        ref = config["data"]["reference"]["genome"]
    output:
        "data/processed/genotyped/wgs_aspen.{chr}.genos.{min_dp}dp{max_dp}.ploidy{genotype_ploidy}.vcf.gz"
    params:
#        db = "gendb://data/interm/combined_database_bpres/{REF}/{interval}",
#        region = "data/processed/scattered_intervals/{REF}/{interval}-scattered.interval_list",
        ploidy = "{genotype_ploidy}",
        tmp =  "/global/scratch/users/arphillips/tmp/joint_geno/{chr}"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/gatk.yaml"
    shell:
        """
        mkdir -p {params.tmp}
        gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {input.vcf} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --tmp-dir {params.tmp} \
        --sample-ploidy {params.ploidy} \
        -O {output}
        rm -rf {params.tmp}
        """
