###
# SnpEff
#
# http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip

# download and install as described in setup.sh


#usage - getting list of command line params
java -jar $SNPEFF -h
java -jar $SNPSIFT -h

#what databases are available - #20,000!
java -jar $SNPEFF databases > database_list
grep "Homo_sapiens" database_list
grep "Homo_sapiens" database_list | wc -l 
#only 16 databases related to human - mostly different genome versions
#like hg19 and hg19kg (known genes), hg38 etc. and different versions of
#GRCh37 or GRCh38
#note that databases are specific to different versions of snpEff - yikes

###
#expect to spend some quality time with the snpEff and SnpSift manuals
###
#spend some quality time here: http://snpeff.sourceforge.net/SnpEff_manual.html
#spend some quality time here: http://snpeff.sourceforge.net/SnpSift.html

#sample usage
mkdir -p ~/ngs/annotation_demo
cd ~/ngs/annotation_demo
time java -Xmx4g -jar $SNPEFF -v GRCh37.75 $SNPEFFDIR/examples/test.chr22.vcf > test.chr22.ann.vcf
#real	2m44.529s
#user	2m23.154s
#sys	0m16.015s
#the majority of real time was spent downloading data

#open snpEff_summary.html with a browser
#then look at gene specific details
head snpEff_genes.txt
head -n 2 snpEff_genes.txt | t2r
wc -l snpEff_genes.txt 
#3455 snpEff_genes.txt

#and the modified VCF file now contains an ANN field
#the contents of the ANN field (16 pipe-delimited fields) is described here
#http://snpeff.sourceforge.net/SnpEff_manual.html#input
#see sections on ANN field and EFF field
#N.B. - multiple effects can be listed for each ANN and each EFF (separated by ,)
#because effects will be listed for upstream and downstream genes and transcripts in the case that
#a SNP occurs between two genes - for a total of four effects
#this will be listed in two ANN entries (separated by a comma) and each ANN entry will have two effects - these two effect types will be listed in the second ANN field separated by a '&' character
#i know...a mess ... i only want to note the arity of the data here
#thankfully, there is SnpSift to help troll through all this
#and if you want to extract a list of genes that have certain effect types (e.g. codon_change) in a given region type (e.g. exon) you are better off using the snpEff_genes.txt output file

#there are also additional sources of annotation not covered here (see http://snpeff.sourceforge.net/SnpEff_manual.html#effNc)



#example protocols
#http://snpeff.sourceforge.net/protocol.html


less test.chr22.ann.vcf
grep -v "^##" test.chr22.ann.vcf | wc -l #15415

#how many variants with high impact
grep -E "^[^##].*HIGH" test.chr22.ann.vcf | wc -l #251

###
# 
#annotate snps having clinical significance

#clinvar publishes results in vcf format that are used by the snpeff demos
zless $SNPEFFDIR/protocols/db/clinvar.vcf.gz
#most recent clinvar is here as installed in setup.sh
ll $CLINVARDIR
java -jar $SNPSIFT -h
java -jar $SNPSIFT annotate -h

cd ~/ngs/annotation_demo
java -Xmx1g -jar $SNPSIFT \
    annotate \
    -v \
    $CLINVARDIR/clinvar.vcf.gz \
    test.chr22.ann.vcf \
    > test.chr22.ann.clinvar.vcf

#there are no clinical variants in this example
grep CLNDBN test.chr22.ann.clinvar.vcf

#instead, try this one
cd ~/ngs/annotation_demo
java -Xmx1g -jar $SNPSIFT \
    annotate \
    -v \
    $CLINVARDIR/clinvar.vcf.gz \
    $SNPEFFDIR/protocols/ex1.ann.cc.clinvar.filtered.vcf \
    > ex1.cc.ann.clinvar.filtered.vcf
#and one entry (cariant) is annotated
grep CLNDBN ex1.cc.ann.clinvar.filtered.vcf

#but it is really hard to read...so you can use this utility
grep CLNDBN ex1.cc.ann.clinvar.filtered.vcf | $SNPEFFDIR/scripts/vcfInfoOnePerLine.pl 



####continue here
