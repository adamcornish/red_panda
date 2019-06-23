# Overview
This will run tests to make sure that everything is set up correctly on your environment.

Before the tests can be run, please make sure that you have most recent version of the fasta sequence and annotation for either the human or mouse genome from Ensembl. 

## Genome and Annotation

For Red Panda to run correctly, the fasta sequence and annotation of your genome, must be downloaded. 

Fasta sequence for the human genome (GRCh38) can be found here: ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/

Fasta sequence for the mouse genome (GRCm38) can be found here: ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/

Annotation for the human genome (release v. 96) can be found here: ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens

Annotation for the mouse genome (release v. 96) can be found here: ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus

## Edit the script to suit your environment

### Variables
At the top of the run_test.sh shell script, a number of variables have been included, two of which need to be changed to fit your environment:

    genome_fasta=/path/to/your/genome/fasta/file.fa
    annotation=/path/to/your/annotation/file.gtf
    
### Samools/Bcftools version
As of right now, Red Panda works with samtools versions 0.1.19–1.8. In version 1.9, mpileup was moved to the bcftools package, so depending on which version you have installed in your environment, you will need to uncomment the appropriate line of the shell script. If you have 1.8 or older installed, uncomment line 13, but if you have version 1.9 installed, uncomment line 15.

## Run
Once you have ensured that you have the genome and annotation, run the script "run_test.sh" to ensure everything is working.
