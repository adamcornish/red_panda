#!/bin/sh

sample=test
input_bam_file=data/aln/test.bam
data_dir=data
genome=hg38
genome_fasta=/inbre-storage-b/adamc/bcbio/genomes/Hsapiens/hg38/seq/hg38.fa
annotation=/inbre-storage-b/adamc/bcbio/genomes/Hsapiens/hg38/rnaseq/ref-transcripts.gtf

mkdir -p $data_dir/{aln,vcf/unfiltered}

# Run samtools mpileup (v. 0.1.19 - 1.8)
samtools mpileup -B -t AD -v -u $input_bam_file -f $genome_fasta | grep -vP "\\s<\\*>\\s" > $data_dir/vcf/unfiltered/$sample.vcf

# Run bcftools mpileup (v. 1.9)
bcftools mpileup -B -a AD,DP $input_bam_file -f $genome_fasta | grep -vP "\\s<\\*>\\s" > $data_dir/vcf/unfiltered/$sample.vcf

# Run Red Panda
perl ../src/red_panda.pl -q data/quant/$sample.sf -v $data_dir/vcf/unfiltered/$sample.vcf -s $sample -d $data_dir -g $genome -A $annotation -V
