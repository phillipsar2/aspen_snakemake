# Variant calling in quaking aspen WGS

Author: Alyssa Phillips

Snakemake (v8.27.1) pipeline for variant calling in a large quaking aspen WGS dataset.
The directory roughly follows a CookieCutter directory structure.

## Running the pipeline

`module load anaconda3 gatk samtools`
`conda activate grenepipe`
`rm -r .snakemake/metadata .snakemake/log .snakemake/slurm_logs`
`snakemake --executor slurm --profile profiles/ --use-conda --rerun-triggers input`

## Project organization
<pre>
├── README.md  
|  
├── rules  
|    ├── calling.smk  
|    ├── ploidy_sex.smk  <- rules to determine ploidy and sex
|    ├── filtering.smk  
|    └── mapping.smk	 <- rules for aligning raw data to the reference  	
|  
├─  environment.yml  
├─  scripts  
│    ├── README.md  
│    └── filtering       <- Custom scripts for variant filtering  
|  
├── data  
│    ├── raw 		 <- Original data dump  
│    ├── trimmed         <- Trimmed fastqs
│    ├── genome 	 <- Reference genome  
│    ├── interim  	 <- Intermediate files in read mapping and SNP calling
|    ├── ref/toz19		 <- Files with for determining sex of samples  
│    └── processed	 <- Final vcfs for analysis  
|  
├── reports 		 <- Generated analyses as HTML, PDF, or .txt.  
├── qc 			 <- Quality check output for raw data  
├── Snakefile  
└── profiles/config.yaml	<- Specifies job resources, snakemake setting  
</pre>

## Workflow overview

1. Pre-processing of reads
* Assess read quality with fastqc
* Trim reads with fastp and re-evaluate quality. Reads are trimmed via sliding windows (4 bp windows, min quality of 15) and automated detection of adapters.

2. Mapping
* Map reads to the reference with bwa-mem2
* Sort, add read groups, and deduplicate BAM files with samtools and GATK.
* Assess mapping quality with qualimap's bamqc
* Assess the DNA damage AdDeam
* As some genotypes were sequenced mutliple times, BAM files from multiple high-quality runs were merged with `samtools merge` and processed as described above

3. Extracting plastid genome reads
* Reads were independently mapped to a P. tremuloides chloroplast genome (MW376839.1) using bwa-mem2
* `samtools sort` was used to sort BAM files
* `samtools fastq` was used to convert BAMs back to fastq files, excluding unmapped reads and singletons. Pair-end reads were seperated into two files

4. Variant calling and filtering
* Variants are initially called with bcftools mpileup. Quality of SNPs is assed between each filtering step.
	Raw SNPs: 90,102,90  
* Variants are hard filtered for biallelic sites, MQ > 40, and QUAL > 40.
	SNPs after hard filtering: 61,364,088
* Variants are then filtered for: 
        10 < DP < 90 & less than 10% missing data: 544,839 
* After filtering, 87 samples (`metadata/samples_to_drop_from_vcf.txt`) were excluded on the critera of: genotype was duplicately sequenced (kept the better of the two), mean coverage < 10, percent of missing SNPs > 20%, ploidy could not be determined (see below), percent of reads mapped < 85% 
* Multiple LD filters were assessed:
	+ one SNP per 500 bp
	+ LD thinning with PLINK with r^2 threshold of 0.1
	+ Estimation of LD with `ldsep`

5. Ploidy determination
* Ploidy was determined using `gbs2ploidy`
* Heterozygous sites from the filtered vcf were used (without LD filter applied, before dropping genotypes)
* The allele ratio with the highest posterior probability was used to assign ploidy (`data/gbs2ploidy/flow_cyt_predictions.csv`)

6. Genotyping
* `updog` was used for genotyping with dosage specified by the called ploidy 

X. Population structure
* 

X. Determining sex
* The TOZ19 sex locus region was identified by mapping the *P. trichocarpa* genomic sequence to the reference with minimap2 (the reference is Male)
* The sex-linked region was identified by mapping the *P. tremuloides* concensus sequence from Pakull et al. (2014) to the reference with minimap2 and blastn using the scripts `scripts/blastn_TOZ19.sh` and  `scripts/minimap_TOZ19.sh`
* Samtools depth is used to extract the per-bp coverage of the TOZ19 region and then mean coverage of the region is calculated.
