# () Merge vcf files
rule merge_filt:
    input:
        vcfs = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf", region = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz"
    params:
        list = "/global/scratch/users/arphillips/data/processed/filtered_snps/vcflist.{min_dp}dp{max_dp}.txt",
        path = "/global/scratch/users/arphillips/data/processed/filtered_snps/",
        suffix = ".nocall.{min_dp}dp{max_dp}.vcf" 
    shell:
        """
        ls {params.path}*{params.suffix} > {params.list} 
        /global/scratch/users/arphillips/toolz/bcftools/bcftools concat -f {params.list} --threads 10 -Oz -o {output}
        """

# (19) Determine ploidy with gbs2ploidy
rule gbs2ploidy:
    input:
        expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz", min_dp = MIN_DP, max_dp = MAX_DP)
    output:
        "/global/scratch/users/arphillips/data/gbs2ploidy/{bam}.propOut.csv"
    params:
        geno = "{bam}",
        tmp_dir = "/global/scratch/users/arphillips/tmp/gbs2ploidy/{bam}",
        temp = "/global/scratch/users/arphillips/tmp/gbs2ploidy/{bam}/{bam}.vcf.gz"
    conda: "/global/home/users/arphillips/.conda/envs/stuff_in_r"
    shell:
        """
        mkdir -p {params.tmp_dir}
        /global/scratch/users/arphillips/toolz/bcftools/bcftools view -Oz -s {params.geno} {input} >> {params.temp}
        Rscript scripts/gbs2ploidy.R {params.temp} --out {output}
        rm -rf {params.tmp_dir}
        """

# (19) Genotype with updog
# Additional filters applied that exclude poor quality genotypes.
# --ploidy is diploid or triploid
rule updog:
    input:
#        vcf = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz", min_dp = MIN_DP, max_dp = MAX_DP)
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.10dp90.vcf.gz",
        meta = "/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions.csv"
    output:
        "/global/scratch/users/arphillips/data/updog/updog.genomat.{ploidy}.{chr}.1206.txt"
    params:
        tmp_dir = "/global/scratch/users/arphillips/tmp/updog/{ploidy}",
        tmp = "/global/scratch/users/arphillips/tmp/updog/{ploidy}/wgs_aspen.{chr}.nocall.10dp90.vcf.gz",
        outdir = "/global/scratch/users/arphillips/data/updog",
        ploidy = "{ploidy}",
        chr = "{chr}"
    shell:
        """
        mkdir -p {params.tmp_dir}
        bcftools view -r {params.chr} {input.vcf} -Oz -o {params.tmp}
        Rscript scripts/updog.R {params.tmp} \
        --ploidy {params.ploidy} --cores 10 --outdir {params.outdir} --chr {params.chr} --meta {input.meta}
        rm -f {params.tmp}
        """
