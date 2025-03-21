# (19) Genotype with updog
# Additional filters applied that exclude poor quality genotypes.
# --ploidy is diploid or triploid
rule updog:
    input:
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{chr}.nocall.10dp75.vcf",
#        meta = "/global/scratch/projects/fc_moilab/aphillips/spectral_aspen/aspendatasite-levelprocessed30Mar2020.csv" 
    output:
        "/global/scratch/users/arphillips/data/updog/updog.genomat.{ploidy}.{chr}.{date}.txt"
    params:
        outdir = "/global/scratch/users/arphillips/data/updog",
        ploidy = "{ploidy}",
        chr = "{chr}"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/updog.yaml"
    shell:
        """
        Rscript scripts/updog.R {input.vcf} \
        --ploidy {params.ploidy} --cores 10 --outdir {params.outdir} --chr {params.chr}
        """
