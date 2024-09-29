# Fastq Processing and Genome Alignment Script

## Description
This Bash script is designed to automate the preprocessing and alignment of FASTQ files using the STAR aligner. It handles genome file downloading,
adapter sequence fetching, and genome indexing for human reference genomes (hg19 and hg38). Additionally, it uses Trimmomatic for trimming low-quality
reads and adapter sequences and performs genome alignment using STAR, outputting both BAM files and gene count tables. The script supports paired-end
and single-end FASTQ files.

## Requirements
- STAR: A fast RNA-seq aligner.
- Trimmomatic: A tool for trimming Illumina sequencing data.
- FastQC: A quality control tool for high throughput sequence data.
- Picard: A tool for manipulating BAM files.
- gtfToGenePred: A UCSC tool for converting GTF files into gene prediction formats.
- wget: For downloading files from the internet.

## Dependencies Installation
You can install required tools using conda or apt:

- Using conda:  
conda install -c bioconda star trimmomatic fastqc picard wget

- Or using apt:  
sudo apt-get install star trimmomatic fastqc picard wget

## Usage
The script takes in the directory containing FASTQ files and the reference genome version (either hg19 or hg38) as input arguments.

Basic Usage:  
./script.sh <FASTQ_DIR> <GENOME> <CROP>  

where,
- <FASTQ_DIR>: Directory containing the FASTQ files.
- <GENOME>: Genome version (hg19 or hg38).
- <CROP>: Crop length for trimming reads (default: 50).

Example:  
./script.sh /path/to/fastq hg19 50

__Default Parameters:__
If no arguments are provided, the script uses default values:
- FASTQ Directory: ./
- Genome: hg19
- Crop: 50

Outputs:
- Trimmed Reads: FASTQ files with trimmed reads.
- Quality Reports: FastQC reports for raw and trimmed FASTQ files.
- Aligned Reads: BAM files containing aligned reads.
- Gene Count Tables: Count tables with gene expression levels.

## Features
- Automatic Genome Download: Downloads reference genome FASTA and GTF files for hg19 or hg38 from UCSC.
- Genome Indexing: If genome indexes are not available, they are generated using STAR.
- Adapter Sequence Download: Automatically downloads Illumina TruSeq adapter sequences.
- FastQC Analysis: Runs FastQC on both raw and trimmed reads.
- Read Trimming: Uses Trimmomatic to remove low-quality bases and adapter contamination.
- Genome Alignment: Uses STAR to align reads to the specified reference genome.

## Workflow
- Directory Setup: The script sets up directories for genomes, indexes, adapters, and other necessary files in the $HOME/GENOMEDIR folder.
- Genome Download: Downloads the necessary reference genome (FASTA) and annotation (GTF) files from UCSC.
- Adapter Sequences: Downloads adapter sequences required for trimming.
- Genome Indexing: If not available, indexes the genome using STAR.
- Trimming: Performs read trimming using Trimmomatic to remove adapters and low-quality bases.
- Alignment: Aligns the trimmed reads to the reference genome using STAR.
- Gene Count Table: Outputs gene expression counts from the alignment.

Example Output Files
- *_trimmed_R1.fastq: Trimmed read 1 FASTQ file.
- *_trimmed_R2.fastq: Trimmed read 2 FASTQ file (paired-end data).
- Aligned.sortedByCoord.out.bam: Sorted BAM file of aligned reads.
- ReadsPerGene.out.tab: Table of read counts per gene.
