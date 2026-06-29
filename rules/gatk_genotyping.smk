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
        -ploidy {params.ploidy} \
        --max-mnp-distance 0
        """

# (2) Merge VCFs
## Merge vcfs for each genotype then do genotyping
## database approach doesn't work for the scale of data
rule merge_gvcfs:
    input: 
        vcf = expand("data/vcf/gatk/called/{{geno}}_p{{geno_ploidy}}.{region}.haplo.g.vcf.gz", region = CHR),
    output:
        "data/vcf/gatk/merged/{geno}_p{geno_ploidy}.merged.g.vcf.gz"
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
rule indv_geno:
    input:
        vcf = "data/vcf/gatk/merged/{geno}_p{geno_ploidy}.merged.g.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
       "data/vcf/gatk/genotyped/{geno}_p{geno_ploidy}.g.vcf.gz"
    params:
        ploidy = "{geno_ploidy}",
        tmpdir =  "global/scratch/users/arphillips/tmp/joint_geno/{geno}",
    shell:
        """
        gatk IndexFeatureFile -I {input.vcf}
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

