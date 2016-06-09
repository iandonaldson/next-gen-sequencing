#use VCFTools perl module to combine multiple VCF files (one per sample) into one VCF file
#vcf-merge will intelligently merge SNVs - mutiple insertions/deletions at the same site
#will be merged into a single variant line with the correct default allele and list of alternative alleles
#some INFO field results should be viewed with caution - they may only be copies of data in the original files and
#may not reflect newly created and re-ordered alternate alleles; espcially AD, PL, VF
#the function of the -c and -d options are unclear in the documentation and did nothing in my limited tests

#vcf-merge
#https://vcftools.github.io/perl_module.html#vcf-merge
#vcf-merge requires compressed and indexed files

#collect vcfs in one directory
mdir combined_vcf
cp *_S*.vcf combined_vcf/.

#vcf-merge requires compressed and indexed files
module load htslib
cd combined_vcf
for THIS in *.vcf; do bgzip ${THIS}; tabix -p vcf ${THIS}.gz; done

#set up to use vcftools perl module installed in my home directory
export PERL5LIB=~/tools/vcftools/vcftools/src/perl/
~/tools/vcftools/vcftools/src/perl/vcf-merge *.vcf.gz > combined.vcf

#compress and make index
bgzip combined.vcf
tabix -p vcf combined.vcf.gz
