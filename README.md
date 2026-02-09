# Variant calling in quaking aspen WGS

Author: Alyssa Phillips

Snakemake (v9.14.5) pipeline for variant calling in a large quaking aspen WGS dataset.
The directory roughly follows a CookieCutter directory structure.

## Running the pipeline

`module load anaconda3 bio/gatk java bio/samtools bio/bcftools`
`conda activate grenepipe`
`module load java bio/gatk bio/bcftools bio/samtools`
`rm -r .snakemake/metadata .snakemake/log .snakemake/slurm_logs`
`snakemake --executor slurm --profile profiles/ --use-conda --rerun-triggers input`

## Project organization
<pre>
├── README.md  
|  
├── rules 
|    ├── mapping.smk	       		<- processing and aligning raw data to the reference
|    ├── calling.smk           		<- variant calling and filtering with bcftools 
|    ├── updog_genotyping.smk  		<- determine ploidy and test rules for updog genotyping
|    ├── gatk_genotyping.smk 		<- variant calling and filtering with GATK
|    ├── mapping_otherpoplar.smk	<- draft rules for mapping other poplar species
|    ├── ploidy_sex.smk				<- finding TOZ19	
|    └── nquack.smk					<- running nQuack for ploidy determination  	
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

1. Pre-processing of reads (`rules/mapping.smk`)
* Assess read quality with fastqc
* Trim reads with fastp and re-evaluate quality. Reads are trimmed via sliding windows (4 bp windows, min quality of 15) and automated detection of adapters.

2. Mapping (`rules/mapping.smk`)
* Map reads to the reference with bwa-mem2
* Sort, add read groups, and deduplicate BAM files with samtools and GATK.
* Assess mapping quality with qualimap's bamqc
* Assess the DNA damage AdDeam
* As some genotypes were sequenced mutliple times, BAM files from multiple high-quality runs were merged with `samtools merge` and processed as described above. This step will removed going forward and replace with only keeping the highest quality BAM. This prevents accidentally merging samples not from the same genotype/clone. 

3. Extracting plastid genome reads (`rules/plastid.smk`)
* Reads were independently mapped to a P. tremuloides chloroplast genome (MW376839.1) using bwa-mem2
* `samtools sort` was used to sort BAM files
* `samtools fastq` was used to convert BAMs back to fastq files, excluding unmapped reads and singletons. Pair-end reads were seperated into two files

4. Ploidy determination
* Variants were intitally called with samtools mpileup, as it is efficient and fast. (`rules/calling.smk`)
* Ploidy was determined using `gbs2ploidy` following recommended protocol. (`rules/updog_genotyping.smk`)
* Heterozygous sites were selected and filtered (MQ > 40, QUAL > 40, DP > 10, DP < 90, bialleleic) 
* The allele ratio with the highest posterior probability was used to assign ploidy (`data/gbs2ploidy/flow_cyt_predictions.csv`)

5. Genotyping and variant filtering (`rules/gatk_genotyping.smk`)
* Variants were called with GATK haplotype caller, per chromosome per genotype to reduce runtime.
* GVCFs were merged for each genotype and GATK GenotypeGVCFs was used to genotype each individual seperately, specifying their ploidy determined with gbs2ploidy.  
* BCFtools merge was used to combine GVCFs into 1 Mb regions. GATK SelectVariants was used to extract SNPs:
	Raw SNPs: 49,979,789 (n = 1,100)
* Variants were hard filtered for biallelic sites, MQ > 40, and QUAL > 40.
	SNPs after hard filtering: 35,819,008
* Variants were then filtered for min genotype depth, genotype depth max three times the mean, and missing data at each variant: 
        10 < DP < 90 & less than 10% missing data: 26,589,983  
* After filtering, 87 samples (`metadata/samples_to_drop_from_vcf.txt`) were excluded on the critera of: genotype was duplicately sequenced (kept the better of the two), mean coverage < 10, percent of missing SNPs > 20%, ploidy could not be determined (see below), percent of reads mapped < 85% 
* Multiple LD filters were assessed:
	+ one SNP per 500 bp
	+ LD thinning with PLINK with r^2 threshold of 0.1
	+ Estimation of LD with `ldsep`

X. Population structure
* 

X. Determining sex
* The TOZ19 sex locus region was identified by mapping the *P. trichocarpa* genomic sequence to the reference with minimap2 (the reference is Male) (`scripts/minimap_TOZ19.sh`)
* The sex-linked region was identified by mapping the *P. tremuloides* concensus sequence from Pakull et al. (2014) to the reference with minimap2 and blastn using the scripts `scripts/blastn_TOZ19.sh` and  `scripts/minimap_TOZ19.sh`
* Samtools depth is used to extract the per-bp coverage of the TOZ19 region and then mean coverage of the region is calculated (`rules/ploidy_sex.smk`). This didn't work very well.
