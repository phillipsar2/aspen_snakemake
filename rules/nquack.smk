# Determining ploidy with nQuack


# (1) Data pre-processing
# Running on each bam in paralell,rather than using a for loop
#
# -U output regions that don't meet filtering criteria
# -L select regions in repeat bed file
# Have to skip the repetitive region filtering for now as we don't have an annotation:
# #  samtools view {input.bam} -h -U -L {input.bed} | \
# #      samtools view {input.bam} -b -q 10 | \
#        samtools mpileup --no-BAQ --ff UNMAP,DUP -A -Q 0 -q 0 /dev/stdin > {output.txt}

# (2) Prepare - convert data to txt file
rule prepare:
    input:
        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam"
    output:
        txt = "/global/scratch/users/arphillips/data/nquack/prepared/{sample}.dedup.txt"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/nquack.yaml"
    params:
        bam = "{sample}.dedup",
        tmp = "/global/scratch/users/arphillips/tmp/nquack_prep/{sample}/"
    shell:
        """
        Rscript scripts/nquack_02prep.R {params.bam} --tmp {params.tmp}
        ""


# (3) Process
rule process:
    input:
        txt = "/global/scratch/users/arphillips/data/nquack/prepared/{sample}.dedup.txt"
    output:
        csv = "/global/scratch/users/arphillips/data/nquack/processed/{sample}.dedup.csv"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/nquack.yaml"
    params:
        samp = "{sample}.dedup"
    shell:
        "Rscript scripts/nquack_03proc.R {params.samp}"

# (4) Model inference
rule infer:
    input:
        csv = "/global/scratch/users/arphillips/data/nquack/processed/{sample}.dedup.csv"
    output:
        csv = "/global/scratch/users/arphillips/data/nquack/model_inference/{sample}.dedup.csv"
    params:
        samp = "{sample}.dedup"
    shell:
        "Rscript scripts/nquack_04infer.R {params.samp}"
