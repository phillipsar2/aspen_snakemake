#!/bin/bash
#SBATCH --job-name=mappingCAM
#SBATCH --account=co_moilab
#SBATCH --partition=savio4_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --time=3-00:00:00
#SBATCH --output  /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/mappingCAM_%j.out
#SBATCH --error  /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/mappingCAM_%j.err
#SBATCH --chdir  /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake

module load bio/samtools bio/gatk java bio/bcftools

### Mapping CAM 1604 short reads to reference ###

ref=/global/scratch/projects/fc_moilab/projects/aspen/genome/CAM1604/Populus_tremuloides_var_CAM1604-4_HAP1_V2_release/Populus_tremuloides_var_CAM1604-4/sequences/Populus_tremuloides_var_CAM1604-4_HAP1.mainGenome.fasta
sample=KWQH_PCRfree_1_1_TGGATCGA_Populus_tremuloides_CAM_1604-4_I1373
r1=/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/genome/KWQH_PCRfree_1_1_TGGATCGA_Populus_tremuloides_CAM_1604-4_I1373_L1_R1.fastq.gz
r2=/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/genome/KWQH_PCRfree_1_1_TGGATCGA_Populus_tremuloides_CAM_1604-4_I1373_L1_R2.fastq.gz
tmp=/global/scratch/users/arphillips/temp/$sample

# Uncompress then map
#/global/scratch/users/arphillips/toolz/bwa-mem2/bwa-mem2 index $ref

#/global/scratch/users/arphillips/toolz/bwa-mem2/bwa-mem2 mem -t 20 $ref $r1 $r2 | \
#samtools view -Sb > /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.mapped.bam

#echo "successfully mapped reads"

# Sort
#mkdir -p $tmp
#samtools sort -T $tmp -@ 1 /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.mapped.bam > /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.sorted.bam
#rm -rf $tmp

#echo "sorted bam"

# add read groups
#mkdir -p $tmp
#gatk --java-options ""-Xmx4G"" AddOrReplaceReadGroups \
#-I /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.sorted.bam \
#-O /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.rg.bam \
#-RGID $sample \
#-RGLB lib1 \
#-RGPL illumina \
#-RGPU unit1 \
#-RGSM $sample \
#--TMP_DIR $tmp \
#--CREATE_INDEX true
#rm -rf $tmp

#echo "mapped bam"

# call variants
bcftools mpileup -Ou \
-f $ref /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.rg.bam \
-A --annotate FORMAT/AD,FORMAT/DP --threads 20 | \
bcftools call -mv -Oz -o /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.vcf.gz
echo "called variants"

bcftools index -t /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/data/interm/mapped_bam/$sample.vcf.gz
echo "indexed vcf"

