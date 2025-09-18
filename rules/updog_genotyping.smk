# () Merge vcf files
rule merge_filt:
    input:
        vcfs = expand("/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.nocall.{min_dp}dp{max_dp}.vcf.gz", region = REGION, min_dp = MIN_DP, max_dp = MAX_DP)
    output:
        "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.all.nocall.{min_dp}dp{max_dp}.vcf.gz"
    params:
        list = "/global/scratch/users/arphillips/data/processed/filtered_snps/vcflist.{min_dp}dp{max_dp}.txt",
        path = "/global/scratch/users/arphillips/data/processed/filtered_snps/",
        suffix = "wgs_aspen.Chr*.nocall.{min_dp}dp{max_dp}.vcf.gz" 
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
        vcf = "/global/scratch/users/arphillips/data/processed/filtered_snps/wgs_aspen.{region}.goodg.10dp90.vcf.gz",
        meta = "/global/scratch/users/arphillips/data/gbs2ploidy/flow_cyt_predictions.csv"
    output:
        "/global/scratch/users/arphillips/data/updog/updog.genomat.{ploidy}.{region}.txt"
    params:
        ploidy = "{ploidy}",
        prefix = "/global/scratch/users/arphillips/data/updog/updog.genomat.{ploidy}.{region}"
    conda: "stuff_in_r"
    shell:
        """
        n=$(zcat {input.vcf} | grep -vc "#")
        if [ "$n" -eq 0 ]; then
            touch {output}
        else
            Rscript scripts/updog.R {input.vcf} --ploidy {params.ploidy} --cores 4 --prefix {params.prefix} --meta {input.meta}
        fi
        """
