configfile: "profiles/config.yaml"

from random import randint

SAMPLE = glob_wildcards("/global/scratch/users/arphillips/raw/jgi_wgs/{sample}.fastq.gz").sample
#SAMPLE = ["53044.3.564626.CTAGGGCCGC-TGGCTCTGTT", "52882.3.482409.CGCATGAT-CGCATGAT", "53044.4.564661.GATTTGGACT-GTGATGGCTC"]
#SAMPLE = ["52896.1.487686.GATGCACTAT-GATGCACTGT", "52865.3.474771.TGCTTGGT-TGCTTGGT", "52864.3.474633.AGCAAGCA-AGCAAGCA"]
#print(SAMPLE)

# =================================================================================================
#     Target Rules
# =================================================================================================
rule all:
    input:
        #fastqc = expand("/global/scratch/users/arphillips/qc/fastqc/{sample}_fastqc.zip", sample = SAMPLE),
        fastp = expand("/global/scratch/users/arphillips/data/trimmed/{sample}.trim.fastq.gz" , sample = SAMPLE),
        #bwa_prep = "/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.gz.0123"
        bam = expand("/global/scratch/users/arphillips/data/interm/mark_dups/{sample}.dedup.bam", sample = SAMPLE),
        bamqc = expand("/global/scratch/users/arphillips/reports/bamqc/{sample}_stats/qualimapReport.html", sample = SAMPLE)

# =================================================================================================
#     Rule Modules
# =================================================================================================
include: "rules/mapping.smk"
