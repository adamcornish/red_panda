dir=/path/to/output/dir
sample=example
genome=mm10
annotation=/path/to/genome/annotation.gtf
perl red_panda.pl -v $dir/$sample.vcf -s $sample -d data -g $genome -A $annotation -V
