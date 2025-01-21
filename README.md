# Variant calling in quaking aspen WGS

Author: Alyssa Phillips

Snakemake pipeline for variant calling in a large quaking aspen WGS dataset.
The directory roughly follows a CookieCutter directory structure.

## Running the pipeline

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
* Likely need to merge bam files from multiple sequencing runs of the same genotype

3. Variant calling and filtering
* Variants are initially called with bcftools mpileup. Quality of SNPs is assed between each filtering step.
	Raw SNPs: 
* Variants are hard filtered for biallelic sites, MQ > 40, and QUAL > 30

4. Determining sex
* The TOZ19 sex locus region was identified by mapping the *P. trichocarpa* genomic sequence to the reference with minimap2 (the reference is Male)
* The sex-linked region was identified by mapping the *P. tremuloides* concensus sequence from Pakull et al. (2014) to the reference with minimap2 and blastn using the scripts `scripts/blastn_TOZ19.sh` and  `scripts/minimap_TOZ19.sh`
* Samtools depth is used to extract the per-bp coverage of the TOZ19 region and then mean coverage of the region is calculated.
