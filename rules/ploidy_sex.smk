# Extract avg coverage in TOZ19 region
rule toz19:
    input:
        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam"
    output:
        "/global/scratch/users/arphillips/data/toz19/{bam}.chr13.cov.txt"
    conda: "/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        samtools depth {input.bam} -r ""Chr13"" | \
        cut -f2,3 > {output}
        """
