# Call SNPs with freebayes
# -g 1000 don't consider SNPs with coverage more than 1000
# --min-alternate-fraction 0 don't require a minimum alt allele read fraction
rule freebayes:
    input:
        ref = config["data"]["reference"]["genome"], 
        bam = "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/bams/{geno}.dedup.bam" 
    output:
        temp("/global/scratch/users/arphillips/data/vcf/freebayes/called/{geno}.{ploidy}.vcf.gz")
    params:
#        region = "{region}",
        ploidy = "{ploidy}"
    shell:
        """
        freebayes -f {input.ref} \
#        -r {params.region} \
        -p {input.ploidy} \
        -g 1000 \
        --min-alternate-fraction 0 \
        {input.bam} \
        > {output.vcf}
        """
