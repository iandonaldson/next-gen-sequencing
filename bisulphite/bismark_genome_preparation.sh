#!/bin/bash


# bismark_genome_preparation

# usage
# modify the parameters section below and then ./bismark_genome_preparation.sh


BISMARK_PATH=/data/home/wgw057/tools/bismark/bismark_v0.16.3
GENOME_DIR=/data/WHRI-GenomeCentre/shares/Projects/NGS_Projects/DNA_Sequencing/Elliot_Evan/GC-EE-5295/Analysis/Ian_Analysis/Elliot_Genome

module load bowtie2


${BISMARK_PATH}/bismark_genome_preparation --verbose --bowtie2 ${GENOME_DIR}

exit

Notes

Bismark User Guide
http://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide.pdf
${BISMARK_PATH}/bismark_genome_preparation --help
${BISMARK_PATH}/bismark --help

Input to the application is a directory path 'GENOME_DIR' which is expected to contain one or more fasta sequence files
with the .fa or .fasta extension. 

the default index will be made using (and compatible with) bowtie2 - if you want to use bowtie for the alignment, a separate index must be made 
specifying --bowtie 

A new directory is made 'Bisulfite_Genome' in the same directory as GENOME_DIR.

Progress is sent to stdout and can be captured to a log.
