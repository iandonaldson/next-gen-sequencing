#!/usr/bin/bash

#prepare a reference genome for bwa alignment

#learning resources
#Prepare a reference for use with BWA and GATK
#https://www.broadinstitute.org/gatk/guide/topic?name=tutorials

#set up directory
mkdir reference
cd reference

###
#retrieve a reference genome
#
#begin browsing at http://hgdownload.cse.ucsc.edu/downloads.html
#drill down to human http://hgdownload.cse.ucsc.edu/downloads.html#human
#drill down to hg38/GRCh38 
#(most recent though many are still using hg19/GRCh37)
#and choose the chromosomes view http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/
#also see README for this directory and README at http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/

mkdir hg38
cd hg38

for i in $(seq 1 22) X Y M; do 
	echo $i; 
	wget --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chr${i}.fa.gz"; 
done
#3 min

du -sh #905 MB
ls -lSh #sort files by size

#check the checksums
wget --timestamping "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/md5sum.txt"
grep -P "chr.{1,2}.fa.gz" md5sum.txt | sort > md5sum.chr.txt
md5sum *.fa.gz | sort > download.md5
diff -s download.md5 md5sum.chr.txt #should be identical

gunzip *.gz
du -sh #3.0 GB

#make a single reference file
for i in $(seq 1 22) X Y M; do 
	cat chr${i}.fa >> hg38.fa; 
done

#store individual chromosome refs for possible later use
gzip *.fa


#test creation of indices on a small chromosome
#generate the bwa index
man bwa #see the index tool
bwa index -a bwtsw chr22.fa 
#about 45 sec
#five index files are created that will be used by bwa
#chr22.fa.amb
#chr22.fa.ann
#chr22.fa.bwt
#chr22.fa.pac
#chr22.fa.sa

#generate the fasta file index
samtools faidx chr22.fa 
#less than 1 second
#one file is created
#chr22.fa.fai - 1 line per fasta sequence in the indexed file where each line details the
#name of the sequence given in the fasta file, its size, location, basesPerLine and bytesPerLine.
#for example
#name,   size,     location, basesPerLine and bytesPerLine.
#chr22   50818468  7         50               51

#generate the sequence dictionary
java -jar $PICARD CreateSequenceDictionary REFERENCE=chr22.fa OUTPUT=chr22.dict 
#less than 1 second
#this creates a file called reference.dict formatted like a SAM header, describing the contents of the reference fasta file. it looks like this:
#@HD     VN:1.4  SO:unsorted
#@SQ     SN:chr22        LN:50818468     M5:221733a2a15e2de66d33e73d126c5109     UR:file:/home/ian/ngs/reference/hg38/chr22.fa


#run the indexing steps for real and track the time
#commented in indivdual steps to get times 
START=$(date +%s.%N)
#bwa index -a bwtsw hg38.fa
#2822 seconds = 47 minutes
#samtools faidx hg38.fa
#20 seconds
java -jar $PICARD CreateSequenceDictionary REFERENCE=hg38.fa OUTPUT=hg38.dict 
#13 seconds
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF






###
#alternatives
#
#it is also possible to retrieve a resource bundle from broad here
#http://gatkforums.broadinstitute.org/discussion/1215/how-can-i-access-the-gsa-public-ftp-server
#location: ftp.broadinstitute.org
#username: gsapubftp-anonymous
#password: <blank>

#sftp gsapubftp-anonymous@ftp.broadinstitute.org
#or point your browser to 
#ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle


#Sources
#https://arvados.org/projects/arvados
#
