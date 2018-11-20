sample=example
genome=/path/to/genome/genome.fa
bam_file=/path/to/bam/file/example.bam
samtools mpileup -B -t AD -v -u $bam_file -f $genome | grep -vP "\\s<\\*>\\s" > $sample.vcf
