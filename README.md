# Variant calling in quaking aspen WGS

Author: Alyssa Phillips

Snakemake pipeline for variant calling in a large quaking aspen WGS dataset.
Directory roughly follows a CookieCutter directory structure.

├── README.md <- The top-level README.md for Andropogon or polyploidy enthusiasts that want to replicate these analyses.   
|
├── rules
|   ├── calling.smk
|   ├── determine_ploidy.smk
|   ├── filtering.smk
|   └── mapping.smk
|
├─  environment.yml 
├─  scripts
│   ├── README.md
│   ├── filtering <- custom scripts for variant filtering
│   ├── common_garden <- scripts related to linear mixed models & phenotype data
│   └── popgen <- analysis and plotting scripts for pop gen analyses
|
├── data
│   ├── raw 		<- Original data dump
│   ├── genome 		<- reference genome
│   ├── interim  	<- intermediate files in read mapping and SNP calling
│   └── processed	<- final vcfs for analysis
|
├── reports 		<- Generated analyses as HTML, PDF, or .txt.
├── qc 			<- quality check output for raw data
├── Snakefile
└── config.yaml 
