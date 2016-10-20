#!/bin/bash


# collate.sh
#
# collect columns from multiple output files into a key by sample matrix
# based on join
# author: Ian Donaldson - i.donaldson@qmul.ac.uk
#

# usage
# review the params below then
# ./collate.sh
#
# set these parameters
PROJECT_DIR=/data/WHRI-GenomeCentre/idonaldson/bisulfite_dev
INPUT_DIR=${PROJECT_DIR}/results_bismark_consol
OUTPUT_DIR=${PROJECT_DIR}/results_collate_unmeth
#multiple columns in the input file may be specified to create a compound key column
KEY_COL="1,2"
#only one column in the input file may be specified as the input column
DATA_COL=6
#samples without a data entry for some key will be assigned an NA value (character or string)
NA="NA"
#an extension string can be removed from the name of uncompressed input files to get the sample name
#so if the input files are named like "sample_name.bismark.cov.gz", you could specify ".bismark.cov"
EXTENSION=".bismark.cov"
#these working directories can be left as is
TMP_DIR=${OUTPUT_DIR}/tmp
KEY_FILE=${TMP_DIR}/keys


###
# no changes required beyond this point
###

#all output files will be written to a single sub-directory
if [ ! -e ${OUTPUT_DIR} ]; then mkdir -p ${OUTPUT_DIR}; fi
if [ ! -e ${TMP_DIR} ]; then mkdir -p ${TMP_DIR}; fi

#cp data files to be consolidated to a tmp dir
find ${INPUT_DIR} -name *.bismark.cov.gz | xargs cp -t ${TMP_DIR}/. 
gunzip ${TMP_DIR}/*


# step 1
#rewrite all input data files in two columns (key and data) format
#this allows for creation of compound keys by pasting together multiple columns specified by KEY_COL
for THIS_FILE in $(ls ${TMP_DIR}); do
     cut --output-delimiter '|' -f ${KEY_COL}  ${TMP_DIR}/${THIS_FILE} > ${TMP_DIR}/key.${THIS_FILE} 
     cut                        -f ${DATA_COL} ${TMP_DIR}/${THIS_FILE} > ${TMP_DIR}/data.${THIS_FILE}
     paste ${TMP_DIR}/key.${THIS_FILE} ${TMP_DIR}/data.${THIS_FILE} > ${TMP_DIR}/key_data.${THIS_FILE}
     rm ${TMP_DIR}/data.${THIS_FILE}
     rm ${TMP_DIR}/key.${THIS_FILE}
     rm ${TMP_DIR}/${THIS_FILE}
done

# step 2
#collect a non-redundant list of keys from all input files
#keep the file sorted and unique in place to keep mem requirements low
echo "!key.sample" > ${KEY_FILE}
for THIS_FILE in $(ls ${TMP_DIR}/key_data.*); do
     cut -f 1 ${THIS_FILE} >> ${KEY_FILE}
     sort -u -o ${KEY_FILE} ${KEY_FILE}
done

# step 3
#a. add an entry corresponding to the header (see key.sample in KEY_FILE)
#b. rewrite the files again so that all keys are present in each with NA's for missing data
for THIS_FILE in $(ls  ${TMP_DIR}/key_data.*); do
        # a
        FILE_NAME=$(basename ${THIS_FILE} ${EXTENSION})
        COL_NAME=${FILE_NAME:9}
        echo -e "!key.sample\t${COL_NAME}" >> $THIS_FILE
        # b
	join -1 1 -2 1 -t $'\t' -a 1 -a 2 -o 0,2.2 -e "${NA}" ${KEY_FILE} <(sort ${THIS_FILE}) > "${TMP_DIR}/complete_${FILE_NAME}"
        #wc -l "${TMP_DIR}/complete_${FILE_NAME}"
done


# step 4
COLLATED=${TMP_DIR}/collated
cp ${KEY_FILE} ${COLLATED}
TMP_COLLATED=${TMP_DIR}/tmp_collated
for THIS_FILE in $(ls ${TMP_DIR}/complete_*); do
	join -1 1 -2 1 -t $'\t' ${COLLATED} <(sort ${THIS_FILE}) > ${TMP_COLLATED}
      mv ${TMP_COLLATED}  ${COLLATED}
done


exit


### implementing a full-outer join using join in bash
the script above does a full outer join on multiple output files that have a common set of keys but where many of the result files will be missing results for some keys

the result files are from bismark methylation extractor and the files are named sample_name.bismark.cov.gz

the script follows 4 main steps
first - make a list of unique keys found in any of the files
second - rewrite each file as just two columns (a key and a data result)
third - rewrite files with all keys (and add an NA to those samples that are missing data for a given key)
fourth - full-outer join to create a matrix on all keys on all samples



###
