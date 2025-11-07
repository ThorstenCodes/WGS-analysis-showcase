#!/bin/bash

# Script demonstrating to run WGS analysis and call germline variants in a human WGS paired end reads 2 X 100bp
# Using GATK4

# Before starting generate a folder structure in your project folder:ter
# mkidr reads aligned_reads scripts data results


######################################################################################################################################################
########################################################## Downloading Files for Demonstration (TO BE DOWNLOADED ONLY ONCE) ##########################
######################################################################################################################################################

# Change all paths in this and the next section to your corresponding project diretory!!!!
# Through the if clause will run only once!!!!! Delete supporting_files folder completely if you need to run it again!
if [ ! -d "/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files" ]; then

    echo "supporting_files directory not found — creating and downloading files..."
    mkdir -p Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files


    # download data
    wget -P /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
    wget -P /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


    echo "Run Prep files..."

######################################################################################################################################################
################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################
######################################################################################################################################################


    # download reference files
    wget -P /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/hg38.fa.gz

    # index ref - .fai file before running haplotype caller
    samtools faidx /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/hg38.fa


    # ref dict - .dict file before running haplotype caller
    gatk CreateSequenceDictionary R=/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/hg38.fa O=/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/hg38.dict


    # download known sites files for BQSR from GATK resource bundle
    wget -P /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    wget -P /Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

else
    echo "supporting_files directory already exists — skipping downloads."
fi

######################################################################################################################################################
################################################### Variant Calling Steps ############################################################################
######################################################################################################################################################

# directories
ref="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/hg38.fa"
known_sites="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/aligned_reads"
reads="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/reads"
results="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/results"
data="/Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/data"


# -----------------------------
# STEP 1: QC - Run fastqc
# -----------------------------

echo "STEP 1: QC - Run fastqc"

#fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
#fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, as quality of reads look okay.

# ----------------------------------------
# STEP 2: Map to reference using BAM-MEM
# ---------------------------------------

echo "STEP 2: Map to refenece $(basename "$ref" .fa) genome using BWA-MEM

# BWA index reference
bwa index ${ref}

# BWA alignment
bwa mem -t 4 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.filt.fastq.gz ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634_paired.sam
