# (1) Call snps initially with bcftools to identify variable sites
# Call on each bam independently, as bcftools doesn't do joint calling anyway
# default only sites with max 250 reads considered at each positin, this is way above the max coverage
# -v option asks to output variant sites only (this is sufficient as we only need heterozygous sites)
# -r output for only the given region
# --annotate FORMAT/AD,FORMAT/DP give allele and genotype depths
# # -b {input.bamlist} -r {params.chr} \
rule mpileup:
    input:
        ref = config["data"]["reference"]["genome"],
        bam = "data/bams/{bam}.dedup.bam"
    output:
        vcf = temp("data/vcf/bcftools/raw/wgs_aspen.{bam}.raw.vcf.gz")
    shell:
        """
        /global/scratch/users/arphillips/toolz/bcftools/bcftools mpileup -Ou \
        -f {input.ref} {input.bam} \
        -A --annotate FORMAT/AD,FORMAT/DP --threads 2 | \
        /global/scratch/users/arphillips/toolz/bcftools/bcftools call -mv -Oz -o {output}
        /global/scratch/users/arphillips/toolz/bcftools/bcftools index -t {output}
        """

# (2) Merge vcf files
rule merge_filt:
    input:
        vcfs = expand("data/vcf/bcftools/raw/wgs_aspen.{bam}.raw.vcf.gz", bam = SAMPLE)
    output:
        temp("data/bcftools/raw/wgs_aspen.all.raw.gz")
    params:
        list = "data/vcf/bcftools/raw/vcflist.txt",
        pattern = "data/vcf/bcftools/raw/wgs_aspen.*.raw.vcf.gz"
    shell:
        """
        ls {params.path}*{params.suffix} > {params.list}
        bcftools merge -l {params.pattern} --threads 10 -Oz -o {output}
        """

# (3) Filter vcf all together
# -v snps keep only SNPs
# -m2 -M2 maximum and minimum allele count is 2 (biallelic)
# QUAL and MQ score greater than 40
# Genotype depth greater than 10 and less than 90
rule bcftools_filt:
    input:
        vcf = "data/bcftools/raw/wgs_aspen.all.raw.gz"
    output:
        temp("data/bcftools/filtered/wgs_aspen.all.snps.gz")
    shell:
        "bcftools view -v snps -m2 -M2 -i 'QUAL>=40 && MQ >=40 && FMT/DP>=10 && FMT/DP<=90' {input.vcf} -Oz -o {output}" 

# (19) Determine ploidy with gbs2ploidy
rule gbs2ploidy:
    input:
        "data/bcftools/filtered/wgs_aspen.all.snps.gz"
    output:
        "/data/gbs2ploidy/{bam}.propOut.csv"
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
