#!/bin/sh
#SBATCH --job-name=blasnt
#SBATCH --account=ac_moilab
#SBATCH --partition=savio3_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --output /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/blastn_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/blastn_%j.err
#SBATCH --chdir /global/home/users/arphillips/aspen/aspen_snakemake

module load bio/blast-plus/2.13.0-gcc-11.4.0

# Align TOZ19 ref to reference genome

gz="/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta.gz"
ref="/global/scratch/projects/fc_moilab/projects/aspen/genome/mex_genome/genome.1MX.fasta"
toz="ref/TOZ19.fasta"
db="ref/mex_genome_db"
out="ref/TOZ19_blastn.txt"

gunzip -c $gz > $ref

# Create a database of your reference contigs to blast against
makeblastdb -in $ref -out $db -dbtype nucl -title "TOZ19" -parse_seqids

# Run the blast
blastn -query $toz -db $db \
	-outfmt '7 qseqid sseqid length qlen slen qstart qend sstart send evalue' \
	-out $out -max_target_seqs 10
