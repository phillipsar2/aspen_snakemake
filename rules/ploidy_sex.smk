# Extract avg coverage in TOZ19 region
rule toz19:
    input:
        bam = "/global/scratch/users/arphillips/data/interm/mark_dups/{bam}.dedup.bam"
    output:
        "/global/scratch/users/arphillips/data/toz19/{bam}.toz19.cov.txt"
    conda: "/global/home/users/arphillips/aspen/aspen_snakemake/envs/samtools.yaml"
    shell:
        """
        samtools index {input.bam}
        samtools depth {input.bam} -r ""Potre.1MX.sc0049:320490-324352"" | \
        cut -f3 | \
        awk ''{{sum += $1}} END {{print sum/NR}}'' > {output}
        """
