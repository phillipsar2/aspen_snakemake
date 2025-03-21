#!/bin/sh
#SBATCH --job-name=getpop
#SBATCH --account=ac_moilab
#SBATCH --partition=savio3_htc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:00
#SBATCH --output /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/otherpop_%j.out
#SBATCH --error /global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/slurm_log/otherpop_%j.err

## SRA

#pop_list="/global/scratch/projects/fc_moilab/aphillips/aspen_snakemake/metadata/other_poplar_sra.txt"

cd /global/scratch/users/arphillips/raw/other_poplars/

#while read line; do
#/global/scratch/users/arphillips/toolz/sratoolkit.3.2.0-centos_linux64/bin/prefetch $line
#/global/scratch/users/arphillips/toolz/sratoolkit.3.2.0-centos_linux64/bin/fastq-dump --gzip --skip-technical --readids --split-3 $line
#done < $pop_list

## GSA

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050703/CRR050703_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050703/CRR050703_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050665/CRR050665_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050665/CRR050665_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050664/CRR050664_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050664/CRR050664_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050663/CRR050663_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050663/CRR050663_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050662/CRR050662_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050662/CRR050662_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050661/CRR050661_f1.fastq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050661/CRR050661_r2.fastq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050702/CRR050702_f1.fq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050702/CRR050702_r2.fq.gz

wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050700/CRR050700_f1.fq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050700/CRR050700_r2.fq.gz

#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050701/CRR050701_f1.fq.gz
#wget https://download.cncb.ac.cn/gsa3/CRA001510/CRR050701/CRR050701_r2.fq.gz
