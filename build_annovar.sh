#!/bin/bash

#build_annovar.sh

HUMANDB="/data/home/wgw057/scratch/annovardb/human"

mkdir -p ${HUMANDB}

module load annovar

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene ${HUMANDB}/

annotate_variation.pl -buildver hg19 -downdb cytoBand ${HUMANDB}/

annotate_variation.pl -buildver hg19 -downdb genomicSuperDups ${HUMANDB}/ 

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all ${HUMANDB}/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2014oct ${HUMANDB}/

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 ${HUMANDB}/ 

annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ljb26_all ${HUMANDB}/
