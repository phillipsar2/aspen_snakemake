#!/bin/sh
#SBATCH --job-name=blasnt
#SBATCH --account=fc_mel
#SBATCH --partition=savio3_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --output /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/blastn_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/blastn_%j.err
#SBATCH --chdir /global/home/users/arphillips/aspen/aspen_snakemake

module load bio/blast-plus/2.13.0-gcc-11.4.0

# Align TOZ19 ref to reference genome

ref="/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_release/Populus_tremuloi
des_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta"
toz="ref/toz19_ref/Potri.019G047300.fasta"
db="ref/CAM1604/CAM1604_db"
out="ref/CAM1604/TOZ19_blastn.txt"

# Create a database of your reference contigs to blast against
makeblastdb -in $ref -out $db -dbtype nucl -title "TOZ19" -parse_seqids

# Run the blast
blastn -query $toz -db $db \
	-outfmt '7 qseqid sseqid length qlen slen qstart qend sstart send evalue' \
	-out $out -max_target_seqs 10
