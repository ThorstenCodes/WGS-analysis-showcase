#!/bin/bash

# Script demonstrating to run WGS analysis and call germline variants in a human WGS paired end reads 2 X 100bp
# Using GATK4

# Structure is established in Dockerfile

# directories
BASE="/GENOMICS_WGS"
ref="$BASE/supporting_files/hg38/hg38.fa"
known_sites="$BASE/supporting_files/hg38/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="$BASE/aligned_reads"
reads="$BASE/reads"
results="$BASE/results"
data="$BASE/data"

######################################################################################################################################################
########################################################## Downloading Files for Demonstration (TO BE DOWNLOADED ONLY ONCE) ##########################
######################################################################################################################################################

# Change all paths in this and the next section to your corresponding project diretory!!!!
# Through the if clause will run only once!!!!! Delete supporting_files folder completely if you need to run it again!
if [ ! -f "$BASE/reads/ftp:/SRR062634_1.filt.fastq.gz" ]; then

    echo "FASTQ files not found — downloading files..."
  #  mkdir -p Users/thorsten/code/ThorstenCodes/Bioinformatics_TK/Projects/Genomics_WGS/supporting_files


    # download data
    wget -P $BASE/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_1.filt.fastq.gz
    wget -P $BASE/reads ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/phase3/data/HG00096/sequence_read/SRR062634_2.filt.fastq.gz


    echo "Run Prep files..."

######################################################################################################################################################
################################################### Prep files (TO BE GENERATED ONLY ONCE) ##########################################################
######################################################################################################################################################


    # download reference files
    wget -P $BASE/supporting_files/hg38/ https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip $BASE/supporting_files/hg38/hg38.fa.gz

    # index ref - .fai file before running haplotype caller
    samtools faidx $BASE/supporting_files/hg38/hg38.fa


    # ref dict - .dict file before running haplotype caller
    gatk CreateSequenceDictionary R=$BASE/supporting_files/hg38/hg38.fa O=$BASE/supporting_files/hg38/hg38.dict


    # download known sites files for BQSR from GATK resource bundle
    wget -P $BASE/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
    wget -P $BASE/supporting_files/hg38/ https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

else
    echo "supporting_files directory already exists — skipping downloads."
fi

######################################################################################################################################################
################################################### Variant Calling Steps ############################################################################
######################################################################################################################################################


# -----------------------------
# STEP 1: QC - Run fastqc
# -----------------------------

echo "STEP 1: QC - Run fastqc"

fastqc ${reads}/SRR062634_1.filt.fastq.gz -o ${reads}/
fastqc ${reads}/SRR062634_2.filt.fastq.gz -o ${reads}/

# No trimming required, as quality of reads look okay.

# ----------------------------------------
# STEP 2: Map to reference using BAM-MEM
# ---------------------------------------

echo "STEP 2: Map to refenece $(basename "$ref" .fa) genome using BWA-MEM"

# BWA index reference
bwa index ${ref}

# BWA alignment
bwa mem \
  -t 4 \
  -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" \
  ${ref} \
  ${reads}/SRR062634_1.filt.fastq.gz \
  ${reads}/SRR062634_2.filt.fastq.gz > ${aligned_reads}/SRR062634_paired.sam

# ------------------------------------------
# STEP 3: Mark Duplicates and Sort - GATK4
# ------------------------------------------

echo "STEP 3: Mark Duplicates and Sort - GATK4"

gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam

# ------------------------------------------
# STEP 4: Base Quality Requalibration
# ------------------------------------------

echo "STEP 4: Base quality recalibration"

# 1. Build the model for requalibration using the variants in the reference
gatk BaseRecalibrator \
    -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam \
    -R ${ref} --known-sites ${known_sites} \
    -O ${data}/recal_data.table

# 2. Adjust Base qualitiy Score
gatk ApplyBQSR \
    -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam \
    -R ${ref} --bqsr-recal-file ${data}/recal_data.table \
    -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam

# ----------------------------------------------------
# STEP 5: Collect Alignment and Insert Size Metrices
# ----------------------------------------------------

echo "STEP 5: Collect Alignment and Insert Size Metrices"

gatk CollectAlignmentSummaryMetrices \
      R=${ref} \
      I=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam \
      O=${aligned_reads}/alignment_metrics.txt


gatk CollectInsertSizeMetrics \
    INPUT=${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam
    OUTPUT=${aligned_reads}/Insert_size_metrics.txt \
    HISTOGRAM_FILE=${aligned_reads}/insert_size_histogram.pdf

# ----------------------------------------------------
# STEP 6: Call Variants - gatk haplotype caller
# ----------------------------------------------------

echo "STEP 6: Call Variants - gatk haplotype caller"

gatk HaplotypeCaller \
    -R ${ref} \
    -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam \
    -O ${results}/raw_variants.vcf

# extract SNP and Indels

gatk SelectVariants \
    -R ${ref} \
    -V ${results}/raw_variants.vcf \
    --select_type SNP \
    -O ${results}/raw_snps.vcf

gatk SelectVariants \
    -R ${ref} \
    -V ${results}/raw_variants.vcf \
    --select_type INDEL \
    -O ${results}/raw_indels.vcf
