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
    run:
        shell("gatk GenotypeGVCFs \
        -R {input.ref} \
        -V {params.vcf} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        --tmp-dir {params.tmp} \
        --sample-ploidy {params.ploidy} \ 
        -O {output}")
