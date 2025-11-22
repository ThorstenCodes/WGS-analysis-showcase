FROM broadinstitute/gatk:4.6.2.0


COPY scripts /scripts

WORKDIR /GENOMICS_WGS

RUN mkdir -p supporting_files/hg38 \
    reads \
    aligned_reads \
    results \
    data
RUN apt-get update && apt-get install bwa samtools fastqc

CMD bash scripts/WGS_Analysis.sh
