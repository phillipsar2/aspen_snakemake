#!/bin/sh
#SBATCH --job-name=minimap
#SBATCH --account=ac_moilab
#SBATCH --partition=savio3_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --output /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/minimap_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/minimap_%j.err
#SBATCH --chdir /global/home/users/arphillips/aspen/aspen_snakemake


# Align TOZ19 ref to reference genome
# -a output to SAM

ref="/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.gz"
toz="ref/TOZ19.fasta"
out="ref/TOZ19_minimap.paf"

~/toolz/minimap2/minimap2 $ref $toz > $out
