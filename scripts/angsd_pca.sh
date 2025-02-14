#!/bin/sh
#SBATCH --job-name=pca
#SBATCH --account=ac_moilab
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00
#SBATCH --output /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/pca_%j.out
#SBATCH --error /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/pca_%j.err
#SBATCH --chdir /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake

module load angsd

bamlist="/global/scratch/users/arphillips/data/interm/mark_dups/bamlist.txt"
ref="/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta"
pre="/global/scratch/users/arphillips/data/angsd/pca/singlepca.3dp70"
pre_chr="/global/scratch/users/arphillips/data/angsd/pca/singlepca.3dp70.chr2"

ls /global/scratch/users/arphillips/data/interm/mark_dups/*bam > $bamlist

#angsd \
#-bam $bamlist -P 56 \
#-GL 2 \
#-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 \
#-doMajorMinor 4 \
#-doCounts 1 \
#-setMinDepthInd 6 -setMaxDepthInd 75 \
#-minInd 150 \
#-ref $ref \
#-doCov 1 \
#-doIBS 1 \
#-makeMatrix 1 \
#-out $pre


## non-optional arguments: doIBS, majorminor, docount, output, docov, intToMajorMinorAA
## majorminor can be determined from bams or genotype likelihoods (-doMajorMinor 4, use ref allele as major allele, requires -ref)

angsd \
-bam $bamlist -P 56 \
-GL 2 \
-uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 2 \
-doMajorMinor 4 \
-doCounts 1 \
-setMinDepthInd 6 -setMaxDepthInd 75 \
-r Chr02 \
-minInd 150 \
-ref $ref \
-doCov 1 \
-doIBS 1 \
-makeMatrix 1 \
-out $pre_chr
