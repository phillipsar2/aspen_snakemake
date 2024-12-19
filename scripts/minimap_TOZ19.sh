#!/bin/sh
#SBATCH --job-name=minimap
#SBATCH --account=fc_mel
#SBATCH --partition=savio3_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --output /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/minimap_%j.out
#SBATCH --error /global/home/users/arphillips/aspen/aspen_snakemake/slurm_log/minimap_%j.err
#SBATCH --chdir /global/home/users/arphillips/aspen/aspen_snakemake


# Align TOZ19 ref to reference genome
# -a output to SAM

ref="/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta"
#ref="/global/scratch/projects/fc_moilab/PROJECTS/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP2.mainGenome.fasta"
#toz="ref/toz19_ref/Potri.019G047300.fasta"
toz="ref/toz19_ref/TOZ19.fasta"
out="ref/CAM1604/TOZ19_minimap.HAP1.paf"

~/toolz/minimap2/minimap2 $ref $toz > $out
