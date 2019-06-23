# Overview
This will run tests to make sure that everything is set up correctly on your environment.

Before the tests can be run, please make sure that you have most recent version of the fasta sequence and annotation for either the human or mouse genome from Ensembl. 

## Genome and Annotation

For Red Panda to run correctly, the fasta sequence and annotation of your genome, must be downloaded. 

Fasta sequence for the human genome (GRCh38) can be found here: ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/

Fasta sequence for the mouse genome (GRCm38) can be found here: ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/

Annotation for the human genome (release v. 96) can be found here: ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens

Annotation for the mouse genome (release v. 96) can be found here: ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus

## Edit the variables to suit your environment
At the top of the run_test.sh shell script, a number of variables have been included, two of which need to be changed to fit your environment:

    genome_fasta=/path/to/your/genome/fasta/file.fa
    annotation=/path/to/your/annotation/file.gtf
    
## Run

Once you have ensured that you have the genome and annotation, run the script "run_test.sh" to ensure everything is working.
