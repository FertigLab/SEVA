#!/bin/bash

# Using python version 2.7.6

ID=${1}

# run the mapsplice command

python MapSplice_multi_threads_2.0.1.9_unix/bin/mapsplice_multi_thread.py  --all-chromosomes-files bowtie_ref/indv_chr/HG19.fasta  --pairend -X 2 -Q fq   --chromosome-dir bowtie_ref/indv_chr/ -x bowtie_ref/HG19.index -1 ${ID}_R1_001.fastq -2 ${ID}_R2_001.fastq -o ${ID}

# change directory to mapsplice directory

cd ${ID}

# Add read groups and sort by coordinate - takes about an hour per file, appeared to have stopped after ~20 minutes
java -Xmx4G -jar picard-tools-1.93/AddOrReplaceReadGroups.jar INPUT=alignments.sam OUTPUT=sorted_genome_alignments.bam RGSM=${ID} RGID=${ID} RGLB=${ID}  RGPL=illumina RGPU=illumina  VALIDATION_STRINGENCY=SILENT SO=coordinate
# replaced ^ with the mapsplice --bam option

# Flagstat - statistics such as total number of reads, reads with QC failure flag set, number of duplicates, percentage mapped, etc.
samtools flagstat sorted_genome_alignments.bam > sorted_genome_alignments.bam.flagstat
# ^ altered to accomodate with the mapsplice --bam option and not using AddOrReplaceReadGroups
#samtools flagstat alignments.bam > sorted_genome_alignments.bam.flagstat


# Index - create index to allow for rapid access in future queries
samtools index sorted_genome_alignments.bam

## Sort by chromosome to allow followup steps to function
perl sort_bam_by_reference_and_name.pl --input sorted_genome_alignments.bam --output sorted_by_chr_read.bam --temp-dir . --samtools samtools

# Translate to transcriptome coords
java -Xmx4G -jar ubu-1.2-jar-with-dependencies.jar sam-xlate --single --bed unc_hg19.bed  --in sorted_by_chr_read.bam  --out transcriptome_alignments.bam  --order rsem_ref/hg19_M_rCRS_ref.transcripts.fa --xgtags --reverse

# Filter indels, large inserts, zero mapping quality from transcriptome bam

java -Xmx4G -jar ubu-1.2-jar-with-dependencies.jar sam-filter --in transcriptome_alignments.bam  --out transcriptome_alignments_filtered.bam --strip-indels --max-insert 10000 --mapq 1

# RSEM summary

rsem-calculate-expression  --paired-end --bam --estimate-rspd -p 8 transcriptome_alignments_filtered.bam rsem_ref/hg19_M_rCRS_ref  rsem


# Strip trailing tabs from rsem.isoforms.results
perl strip_trailing_tabs.pl --input rsem.isoforms.results  --temp orig.isoforms.results

# Prune isoforms from gene quant file keep only genes*
mv rsem.genes.results orig.genes.results
sed /^uc0/d orig.genes.results > rsem.genes.results

# Normalize gene quant changed -c 2 to 5
perl quartile_norm.pl -c 5 -q 75 -t 1000 -o rsem.genes.normalized_results rsem.genes.results

# Normalize isoform quant changed -c 2 to 5
perl quartile_norm.pl -c 5 -q 75 -t 300 -o rsem.isoforms.normalized_results rsem.isoforms.results



