#!/usr/bin/bash


#learning resources
#Install all software packages required to follow the GATK Best Practices.
#https://www.broadinstitute.org/gatk/guide/topic?name=tutorials

cd ~
mkdir -parents Downloads
cd Downloads

###
# bwa - Burrows-Wheeler transformation aligner
#
#main site: http://bio-bwa.sourceforge.net/
#code: https://github.com/lh3/bwa
#download: http://sourceforge.net/projects/bio-bwa/files/ latest tarball or bwakit (precompiled)
#documentation: man bwa or http://bio-bwa.sourceforge.net/bwa.shtml
#questions: bio-bwa-help@sourceforge.net and http://seqanswers.com/ and https://www.biostars.org/
#mailing list signup: https://lists.sourceforge.net/lists/listinfo/bio-bwa-help
wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2
tar --bzip2 -xvf bwa-0.7.12.tar.bz2
cd bwa-0.7.12 
make
#there is no make target for install so
sudo mkdir --parents --mode 755 /usr/local/bin /usr/local/share/man/man1
sudo install --preserve-timestamps bwa  /usr/local/bin
sudo install --preserve-timestamps --mode 644 bwa.1 /usr/local/share/man/man1
bwa
man bwa
cd ..

###
# samtools
#
#main site: http://www.htslib.org/
#old site: http://samtools.sourceforge.net/  <== documentation is here
#code: https://github.com/samtools/samtools/archive/develop.zip
#download site: http://sourceforge.net/projects/samtools/files/samtools/1.2/
cd ~
mkdir -parents Downloads
cd Downloads

wget http://sourceforge.net/projects/samtools/files/samtools/1.2/samtools-1.2.tar.bz2
wget http://sourceforge.net/projects/samtools/files/samtools/1.2/bcftools-1.2.tar.bz2
wget http://sourceforge.net/projects/samtools/files/samtools/1.2/htslib-1.2.1.tar.bz2
tar --bzip2 -xvf samtools-1.2.tar.bz2
cd samtools-1.2 
less README
less INSTALL
#make has a dependency on zlib development libs
sudo yum install zlib-devel  
make
#by default, install is in /usr/local
sudo make install
#or specify your own directory
#make prefix=~/ngs/usr/bin install
samtools
man samtools
cd ..

###
# tabix - destributed on samtools site
# http://www.htslib.org/download/
# now distributed as part of htslib - includes tzbix, bgzip and htsfile

# historical only - see http://www.htslib.org/download/ above instead
# http://sourceforge.net/projects/samtools/files/tabix/
# cd ~
# mkdir -parents Downloads
# cd Downloads
# wget http://sourceforge.net/projects/samtools/files/tabix/tabix-0.2.6.tar.bz2
# tar --bzip2 -xvf tabix-0.2.6.tar.bz2
# cd tabix-0.2.6
# make
# no target for install so...
# sudo install --preserve-timestamps tabix  /usr/local/bin
# sudo install --preserve-timestamps bgzip  /usr/local/bin
# sudo install --preserve-timestamps --mode 644 tabix.1 /usr/local/share/man/man1
# tabix
# man tabix #includes bgzip man


###
# picard
#
#main site: http://broadinstitute.github.io/picard/
#mailing list: use the samtools list
#code: http://github.com/broadinstitute/picard
#usage info (after installing as described below): java -jar $PICARD
#also online here: http://broadinstitute.github.io/picard/command-line-overview.html#Overview
#picard can generate a large number of metrics described here:
#http://broadinstitute.github.io/picard/picard-metric-definitions.html
cd ~
mkdir -parents Downloads
cd Downloads

#install pre-compiled picard tools
wget https://github.com/broadinstitute/picard/releases/download/1.131/picard-tabixtools-1.131.zip
unzip picard-tools-1.131.zip
cd picard-tools-1.131
java -jar picard.jar #will display usage info
#put it in the same place as all other execs
sudo cp *.jar /usr/local/bin/.
java -jar /usr/local/bin/picard.jar #will display usage info
cd ..
#make it a little easier to access
echo "export PICARD='/usr/local/bin/picard.jar'" >> ~/.bashrc
#check it
bash
echo $PICARD
java -jar $PICARD
#in case you don't understand this, see this tutorial 
https://www.digitalocean.com/community/tutorials/how-to-read-and-set-environmental-and-shell-variables-on-a-linux-vps

###
# note on java options
man java
# picard tools have been built expecting (optimized for 2g of RAM) hence...
# java -Xmx2g ... sets the maximum memory available to the JVM

###
#build picard tools from source
#
#this is optional and easy but confusing due to the way the source is packaged

#install or update ant
sudo yum install ant

#build the samtools htsjdk.jar (high-throughput sequencing java development kit)
cd ~
mkdir -parents Downloads
cd Downloads

wget https://github.com/samtools/htsjdk/tarball/master
tar -xvf master
cd samtools-htsjdk-6d2c9dd
ant htsjdk-jar
ls -lt dist #note htsjdk-1.132.jar
cd ..

#build picard from source
#get most recent release archive here: https://github.com/broadinstitute/picard/releases/latest
cd ~
mkdir -parents Downloads
cd Downloads

wget https://github.com/broadinstitute/picard/archive/1.131.zip
#mv the htsjdk built above to this directory (required for build of picard tools)
mv samtools-htsjdk-6d2c9dd/ picard-1.131/htsjdk
cd picard-1.131/
ant -lib lib/ant package-commands
ls -lt dist #picard.jar is here
cd ..

###
# GATK
#
cd ~
mkdir -parents Downloads
cd Downloads
# signup then download from here
# https://www.broadinstitute.org/gatk/download/


cd ~
mkdir -p gatk
cd gatk
mv ../Downloads/GenomeAnalysisTK-3.4-0.tar.bz2 .
tar --bzip2 -xvf GenomeAnalysisTK-3.4-0.tar.bz2
sudo cp GenomeAnalysisTK.jar /usr/local/bin/.
#make it a little easier to access
echo "export GATK='/usr/local/bin/GenomeAnalysisTK.jar'" >> ~/.bashrc
#check it
bash
echo $GATK
java -jar $GATK -h


###
# IGV and igvtools
#
cd ~
mkdir -parents Downloads
cd Downloads
#download from https://www.broadinstitute.org/software/igv/download

cd ~
mkdir -p igv

mv Downloads/IGV_2.3.52.zip igv/.
mv Downloads/igvtools_2.3.52.zip igv/.

unzip igv/IGV_2.3.52.zip
unzip igv/igvtools_2.3.52.zip
rm -rf igv/*.zip

mv igv/IGV_2.3.52 ~/igv
mv igv/IGVTools ~/.

#to start IGV
cd ~/igv
. igv.sh


###
# Download and install AsperaConnect web plugin
# 
# http://downloads.asperasoft.com/connect2/
# documentation link on the same page or
# http://d3gcli72yxqn2z.cloudfront.net/connect/docs/user/linux/en/pdf/Connect_User_3.6.0_Linux.pdf
sudo yum install GLIB
cd ~/Downloads
tar -zxvf apera-connect*
sh aspera-connect*
#will appear as a button in firefox on aspera relevant pages
#to execute manually use
#~/.aspera/connect/bin/asperaconnect
#test - see directions in documentation
#http://demo.asperasoft.com/aspera/user/?B=%2Faspera-test-dir-tiny
#also try here
#http://www.ncbi.nlm.nih.gov/public/

##############
#command line client
#http://downloads.asperasoft.com/en/downloads/50
#sudo sh ascp-install-3.5.4.102989-linux-64.sh
#man ascp
#lots of discussion on biostars
##############

#aspera calculator
#http://asperasoft.com/performance-calculator/



###
# Download data resources
# http://downloads.asperasoft.com/connect2/

cd ~
mkdir -p resources

#dbSNP141
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/snp141Common.txt.gz

#dbSNP in SNV format
#see: http://www.ncbi.nlm.nih.gov/variation/docs/human_variation_vcf/
#ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_GRCh38/VCF/ #note the build subdirectory
wget ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_GRCh38/VCF/All.vcf.gz
gunzip All.vcf.gz
mv All.vcf All.known.human.snps.hg38.vcf

#GATK resource bundle
## - these resources are for use with hg19
## - there is no GATK resource bundle for hg38
## 
GATK_FTP=ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/2.8/hg19

#download library file for GATK's VariantRecalibrator
for f in  hapmap_3.3.hg19.sites.vcf.gz \
          1000G_omni2.5.hg19.sites.vcf.gz \
          1000G_phase1.indels.hg19.sites.vcf.gz \
          Mills_and_1000G_gold_standard.indels.hg19.vcf.gz ; do \
                echo $f; \
                curl --remote-name ${GATK_FTP}/${f}; \
done


###
# Download sample data for analysis from 1000G
#
#
#start here
# http://www.ebi.ac.uk/ena
# http://www.ebi.ac.uk/ena/data/warehouse/search
# Taxon: 9606
# Library layout: PAIRED
# Instrument mode:l Illumina HiSeq 2000
# Instrument platform: ILLUMINA
# Library selection: RNADOM
# Library strategy: WGS
# Library source: Genomic
# Free text search: 1000G to find individuals from 1000 genomes project
#
#chose this one - only 8 GB each - 19B bases - individual NA20362
#http://www.ebi.ac.uk/ena/data/view/ERX232571
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR257/ERR257982/ERR257982_2.fastq.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR257/ERR257982/ERR257982_1.fastq.gz
#there is also this one used in the cbw course - 50 GB each file - 130B bases
#http://www.ebi.ac.uk/ena/data/view/ERX237515
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997/ERR262997_1.fastq.gz
#ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR262/ERR262997/ERR262997_2.fastq.gz

#try retrieving small files 
cd /storage/1000G
START=$(date +%s.%N)
#wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR257/ERR257982/ERR257982_2.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR257/ERR257982/ERR257982_1.fastq.gz
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF
#56 min for 16 GB - rate averages 5 Mb/s

#now try with aspera
#http://www.ebi.ac.uk/ena/browse/read-download#downloading_files_aspera
#http://download.asperasoft.com/download/docs/ascp/3.5.2/html/index.html
#man ascp
#ascp -QT -l 300m -i <aspera connect installation directory>/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:<file or files to download> <download location>
#
#-Q enable fair transfer rate to play nicely with other traffic
#-T disable encryption for max throughput
#-l max rate
#-k2 recover transfer after fatal interupt
#-L directory for a log file (default is /var/log/messages)
#--overwrite=always overwrite pre-existing files of the same name
#-i private key file - typically in $HOME/.ssh/id_[algorithm]
#  locate asperaweb_id_dsa.openssh
#  /home/ian/.aspera/connect/etc/asperaweb_id_dsa.openssh

cd /storage/1000G
START=$(date +%s.%N)
#ascp -T -l 300m -k2 -L /storage/1000G --overwrite=always -i /home/ian/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR257/ERR257982/ERR257982_2.fastq.gz .
ascp -QT -l 50m -k2 -L /storage/1000G --overwrite=always -i /home/ian/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/ERR257/ERR257982/ERR257982_1.fastq.gz .
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF
#rate averages 55 Mb/s - for each file transfer was 100% after 5 min but took until 21 minutes to complete?
#recommended settings - https://support.asperasoft.com/entries/20151798-Aspera-and-Concurrency


###
# Download VCF data for 1000G project

#VCF data is distributed per chromosome for all subjects or for whole genome sequence (wgs) as a summary
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/
#	ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz      #209 MB
#	ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi
#       ...
#	ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz 		#2GB
#	ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz.tbi
#
#lots of additional data (including functional annotation) can be found here
#ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/

#example use of 1000G aspera site - see http://www.1000genomes.org/data#DataAccess
#ascp -i bin/aspera/etc/asperaweb_id_dsa.putty -Tr -Q -l 100M -L- fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz ./

cd /storage/1000G

#ascp does not appear to support globbing so multiple files are downloaded using a for-in loop 
START=$(date +%s.%N)
for FILENAME in \
ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz \
ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz.tbi \
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz.tbi \
ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz; do
ascp -T -l 300m -k2 \
-L /storage/1000G \
--overwrite=always \
-i /home/ian/.aspera/connect/etc/asperaweb_id_dsa.openssh \
fasp-g1k@fasp.sra.ebi.ac.uk:/vol1/ftp/release/20130502/$FILENAME .;
done
END=$(date +%s.%N)
DIFF=$(echo "$END - $START" | bc)
echo $DIFF
#347 s = 6 min

###
#a short introduction to working with large files
FILENAME=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz
#take a look
zless $FILENAME
#how many lines
zcat $FILENAME | wc -l # 84 M
#uncompressed size ; i.e. how many characters
zcat $FILENAME | wc -c # 13 GB
#extract a line
# do this first to avoid surprises - returns just 1 record - cf time for tabix below
zgrep rs367896724 $FILENAME | wc -l 
time zgrep rs367896724 $FILENAME
#real	1m1.647s  <=="wall time"
#user	1m5.360s  <==amount of CPU time spent in user-mode code (outside the kernel)
#sys	0m3.415s  <==amount of CPU time spent in kernel/supervisor-mode code (inside the kernel)
#u+s=CPU time from all threads and can exceed real time
#r time includes time waiting for other processes to complete such as read/wrote that will not be included in u+s

#tabix - depends on the presence of a precomputed index file called FILENAME.tbi 
tabix $FILENAME 1:10,000,000-20,000,000 > tmp
tabix $FILENAME 1:10000000-20000000 > tmp
time tabix $FILENAME 1:10000000-20000000 > tmp
#real	0m0.289s
#user	0m0.264s
#sys	0m0.025s
wc -l tmp 
#294767 - enough said


#looking at individual data
FILENAME=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
zless $FILENAME

tabix $FILENAME 22:16050075-16050075
#the index file either does not exist or is older than the vcf file. Please reindex.
#this is because the tbi file was downloaded first - 
#you could do
stat $FILENAME.tbi
touch -m $FILENAME.tbi
#then the above tabix command will work
#but lets re-index anyway
mv $FILENAME.tbi $FILENAME.tbi.bk
time tabix -p vcf -s 1 -b 2 -e 2 -S 0 -c '#' $FILENAME
#30 s

tabix $FILENAME 22:16050075-16050075 > tmp; wc -l
#1 line
#how many columns?
tabix $FILENAME 22:16050075-16050075 | t2r > tmp2; wc -l tmp2
#2516
#how many individuals have at least 1 allele?
grep -c '0|1' tmp2
grep -n '0|1' tmp2
#only one individual in column 1447
RESULT=$(grep -n '0|1' tmp2)

#tabix works on remote files !!
tabix -h ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz 2:39967768-39967768
#in this case, the corresponding .tbi index file will be downloaded to the local directory and tabix will only return the portion of the file corresponding to query

#but tabix is limited to extracting parts of a vcf file (or other chromosome:location based file)
#by location - consider if you have to extract by some other feature
#
#example - try to pick out indels from a vcf file
FILENAME=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
FILENAME=ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.vcf.gz
zgrep -v "^#" $FILENAME | awk '{if(length($4) != length($5)) {print $0}}' | less
#immediate results from the top of the file...but how long does it take to complete entire file?
time zgrep -v "^#" $FILENAME | awk '{if(length($4) != length($5)) {print $0}}' > tmp
#real	1m1.712s - not bad
#subsections of vcf files based on individual id, population etc. can be retrieved using vcftools (just like how samtools can help retrieve subsections of bam or sam files).
#1000G provides a web interface to these tools called  Data Slicer Tool on their tool page (http://browser.1000genomes.org/tools.html) but its probably better to install and learn these tools yourself

FILENAME=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

#generate 
vcftools --gzvcf $FILENAME
#reports number of individuals and sites

#subset just the indels like above
time vcftools --gzvcf $FILENAME --keep-only-indels --recode --recode-INFO-all --out INDELs_only
#1m24.626s - almost as good as 'grep and awk' method above

#snp filter for rs587593704 (can also exclude or specify multiple in file)
FILENAME=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
time vcftools --gzvcf $FILENAME --snp rs587593704 --recode --recode-INFO-all --out rs557400692_only
#28 s - seems a little slow

#individual filter for HG00118 (can also exclude or specify multiple in file)
FILENAME=ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
time vcftools --gzvcf $FILENAME --indv HG00118 --recode --recode-INFO-all --out rs557400692_only
#8 min

#select a random number of individuals up to some maximum 

#output as a matrix using --012
mkdir tmp_dir
time vcftools --gzvcf $FILENAME --012 --temp tmp_dir --out matrix
###does not work --- 
# Writing 012 matrix files ... Error:  Could not open temporary file.
ulimit -n 3000 #to allow many multiple files
vcftools --gzvcf $FILENAME --012 --temp tmp_dir --out matrix
#now it works
#568 sec

###
# vcftools
# http://vcftools.sourceforge.net/
# http://www.1000genomes.org/faq/how-do-i-get-sub-section-vcf-file

cd ~
mkdir -parents Downloads
cd Downloads
wget http://sourceforge.net/projects/vcftools/files/vcftools_0.1.12b.tar.gz
gunzip vcftools_0.1.12b.tar.gz
tar -xvf vcftools_0.1.12b.tar
cd tar vcftools_0.1.12b
make
#there is no make target for install so
sudo mkdir --parents --mode 755 /usr/local/bin /usr/local/share/man/man1
sudo install --preserve-timestamps bin/vcf*  /usr/local/bin
sudo install --preserve-timestamps bin/fill*  /usr/local/bin
sudo install --preserve-timestamps perl/Vcf.pm  /usr/lib64/perl5
sudo install --preserve-timestamps perl/VcfStats.pm  /usr/lib64/perl5
sudo install --preserve-timestamps perl/FaSlice.pm  /usr/lib64/perl5
sudo install --preserve-timestamps --mode 644 bin/man1/vcftools.1 /usr/local/share/man/man1
man vcftools
#or http://vcftools.sourceforge.net/docs.html



###
# pIRS (profile based Illumina pair-end Reads Simulator
#
#developed for de novo data simulation. It uses empirical distribution to reproduce Illumina pair-end reads with #real distribution of substitution sequencing errors, quality values and GC%-depth bias.
#For future information, please refer to https://github.com/galaxy001/pirs.
#http://bioinformatics.oxfordjournals.org/content/28/11/1533
#https://github.com/galaxy001/pirs/archive/master.zip

cd ~/Downloads
wget https://github.com/galaxy001/pirs/archive/master.zip
unzip master.zip
cd pirs-master
#may have to install several libs (BOTH static and devel versions) - see the INSTALL readme.
make
q
#put it in the same place as all other execs
sudo cp src/pirs/pirs /usr/local/bin/.
#check it - will show usage info
pirs

#store the profiles data under resources
mv pirs ~/ngs/resources/.

echo "export PIRS=~/ngs/resources/pirs" >> ~/.bashrc
exec bash
echo $PIRS


###
# Annotation tools
###

###
# ClinVar
# where the rubber hits the road
# spend some quality time here: ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/README.txt

cd ~/Downloads
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
wget ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz.tbi
# zless clinvar.vcf.gz
mkdir -p ~/ngs/resources/clinvar
cp clinvar.vcf.gz ~/ngs/resources/clinvar/.
gunzip ~/ngs/resources/clinvar/clinvar.vcf.gz
mv clinvar.vcf.gz ~/ngs/resources/clinvar/.
mv clinvar.vcf.gz.tbi ~/ngs/resources/clinvar/.
echo "export CLINVARDIR=~/ngs/resources/clinvar" >> ~/.bashrc



###
# SnpEff
#
# http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip

cd ~/Downloads
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
cd snpEff
#get the protocols example data - see http://snpeff.sourceforge.net/protocol.html
wget http://sourceforge.net/projects/snpeff/files/protocols.zip
unzip protocols.zip
rm protocols.zip
cd ..

#i am going to install snpEff in the ngs/resources directory
#snpEff automatically downloads required databases into a subdirectory of snpEff called data

mv snpEff ~/ngs/resources/.
echo "export SNPEFF=~/ngs/resources/snpEff/snpEff.jar" >> ~/.bashrc
echo "export SNPSIFT=~/ngs/resources/snpEff/SnpSift.jar" >> ~/.bashrc
echo "export SNPEFFDIR=~/ngs/resources/snpEff" >> ~/.bashrc
exec bash
echo $SNPEFF
echo $SNPSIFT
echo $SNPEFFDIR

# see example usage in annotation.sh


#Sources
# Google na12878_r1.fq.gz
# https://arvados.org/attachments Single_Sample_SNV_Pipeline.txt

