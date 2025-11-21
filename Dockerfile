FROM broadinstitute/gatk:4.6.2.0


COPY scripts /scripts
COPY

WORKDIR /GENOMICS_WGS

RUN mkdir -p supporting_files/hg38 \
    reads \
    aligned_reads \
    results \
    data
RUN apt-get update && apt-get install bwa gatk4 samtools fastqc
