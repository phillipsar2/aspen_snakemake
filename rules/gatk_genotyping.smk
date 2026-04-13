#ruleorder: haplotype_caller > merge_gvcfs > indv_geno > exclude_MNPs > genomicsdb > big_gvcf

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
# Need to specify ploidy here, can't be undone in later steps.
# --max-mnp-distance 0 is important for excluding MNPs for GenomicsDBImport
rule haplotype_caller:
    input:
        ref = config["data"]["reference"]["genome"], 
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{geno}.dedup.bam" 
    output:
        "data/vcf/gatk/called/{geno}_p{geno_ploidy}.{region}.haplo.g.vcf.gz"
    params:
        region = "{region}",
        ploidy = "{geno_ploidy}"
    shell:
        """
        gatk --java-options "-Xmx4g" HaplotypeCaller \
        --input {input.bam} \
        --output {output} \
        --reference {input.ref} \
        -ERC BP_RESOLUTION \
        -G StandardAnnotation \
        -L {params.region} \
        -ploidy {params.ploidy}
        """

# (2) Merge VCFs
## Merge vcfs for each genotype then do genotyping
## database approach doesn't work for the scale of data
rule merge_gvcfs:
    input: 
        vcf = expand("data/vcf/gatk/called/{{geno}}_p{{geno_ploidy}}.{region}.haplo.g.vcf.gz", region = CHR),
    output:
        "data/vcf/gatk/merged/{geno}_p{geno_ploidy}.g.vcf.gz"
    params:
        pre = "data/vcf/gatk/called/{geno}",
        list = "data/vcf/gatk/called/{geno}.list"
    shell:
        """
        ls {params.pre}*.g.vcf.gz > {params.list}
        gatk MergeVcfs \
        -I {params.list} \
        -O {output}
        rm {params.list}
        """

# (3) Individually genotype with correct ploidy level specified
## Merge vcfs for each genotype then do genotyping
### double bracket masks geno wildcard in input
rule indv_geno:
    input:
        vcf = "data/vcf/gatk/merged/{geno}_p{geno_ploidy}.g.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
       "data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz"
    params:
        ploidy = "{geno_ploidy}",
        tmpdir =  "global/scratch/users/arphillips/tmp/joint_geno/{geno}",
    shell:
        """
        mkdir -p {params.tmpdir}
        gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {input.vcf} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --include-non-variant-sites \
        --tmp-dir {params.tmpdir} \
        --sample-ploidy {params.ploidy} \
        -O {output}
        rm -rf {params.tmpdir}
        """

# (4) Merge vcfs
#> cat list.txt | while read file
#> do
#> bcftools index -t "${file}"
#> done

rule bcftools_merge:
    input:
        vcfs = expand("data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz", zip,  geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY)
    output:
        "data/vcf/gatk/called/wgs_aspen.all.genos.{region}.g.vcf.gz"
    params:
        chr = "{region}",
        vcfs = lambda wildcards, input: input.vcfs 
    shell:
        "bcftools merge {params.vcfs} -m all -r {params.chr} --threads 5 -Oz -o {output}"


