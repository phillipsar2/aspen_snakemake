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
        "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.{region}.haplo.g.vcf.gz"
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
rule merge_gvcfs:
    input: 
        vcf = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{{geno}}_p{{geno_ploidy}}.{region}.g.vcf.gz",region = CHR),
    output:
        "/global/scratch/users/arphillips/data/vcf/gatk/merged/{geno}_p{geno_ploidy}.haplo.g.vcf.gz"
    params:
        pre = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}",
        list = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}.list"
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
        vcf = "/global/scratch/users/arphillips/data/vcf/gatk/merged/{geno}_p{geno_ploidy}.g.vcf.gz",
        ref = config["data"]["reference"]["genome"]
    output:
       "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz"
    params:
        ploidy = "{geno_ploidy}",
        tmpdir =  "/global/scratch/users/arphillips/tmp/joint_geno/{geno}",
        pre = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}"
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
        vcfs = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz", zip,  geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY)
    output:
        "/global/scratch/users/arphillips/data/vcf/gatk/called/wgs_aspen.all.genos.{region}.g.vcf.gz"
    params:
        chr = "{region}",
        vcfs = lambda wildcards, input: input.vcfs 
    shell:
        "bcftools merge {params.vcfs} -m all -r {params.chr} --threads 5 -Oz -o {output}"


#rule exclude_MNPs:
#    input:
#        vcf = "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.g.vcf.gz",
#        ref = config["data"]["reference"]["genome"]
#    output:
#       "/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.nomnps.g.vcf.gz"
#    shell:
#        "gatk SelectVariants -R {input.ref} -V {input.vcf} --select-type-to-exclude MNP -O {output}"

# (3) Create a database with GenomicsDB
#rule split_intervals:
#    input:
#        ref = config["data"]["reference"]["genome"]
#    output:
#        int = expand("/global/scratch/users/arphillips/data/processed/scattered_intervals/{intv}-scattered.interval_list", intv = INTERVALS)
#    params:
#        regions = config["data"]["reference"]["contigs"],
#        dir = "/global/scratch/users/arphillips/data/processed/scattered_intervals/"
#    shell:
#        "gatk SplitIntervals -R {input.ref} -L {params.regions} --scatter-count 200 -O {params.dir}"

#rule genomicsdb:
#    input:
#        gvcfs = expand("/global/scratch/users/arphillips/data/vcf/gatk/called/{geno}_p{geno_ploidy}.nomnps.g.vcf.gz", zip,  geno = GENOTYPE, geno_ploidy = GENOTYPE_PLOIDY),
#        region = "/global/scratch/users/arphillips/data/processed/scattered_intervals/{intv}-scattered.interval_list",
#        map = "/global/scratch/users/arphillips/data/vcf/gatk/called/aspen.sample_map"
#    output:
#        directory("/global/scratch/users/arphillips/data/interm/combined_database_bpres/{intv}")
#    params:
#        tmp = "/global/scratch/users/arphillips/tmp/genomicsdbimport/{intv}"
#    wildcard_constraints:
#        intv = r"\d{4}"
#    shell:
#        """
#        mkdir -p {params.tmp}
#        gatk --java-options \"-Xmx90g -Xms90g\" \
#        GenomicsDBImport \
#        --genomicsdb-workspace-path {output} \
#        --batch-size 50 \
#        --reader-threads 8 \
#        --sample-name-map {input.map} \
#        --intervals {input.region} --tmp-dir {params.tmp}
#        rm -rf {params.tmp}
#        """

# (4) Conslidate into large gvcf
#rule big_gvcf:
#    input:
#        ref = config["data"]["reference"]["genome"],
##        db = "gendb://global/scratch/users/arphillips/data/interm/combined_database_bpres/{interval}",
#        db2 = "/global/scratch/users/arphillips/data/interm/combined_database_bpres/{intv}"
#    output:
#        "/global/scratch/users/arphillips/data/vcf/gatk/genotyped/wgs_aspen.all.genos.{intv}.g.vcf.gz"
#    params:
#        db = "gendb://global/scratch/users/arphillips/data/interm/combined_database_bpres/{intv}"
#    wildcard_constraints:
#        intv = r"\d{4}"
#    shell:
#        "gatk SelectVariants -R {input.ref} -V {params.db} -O {output}"

