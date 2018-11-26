#!/bin/bash

# run mixcr
# usage: bash mixcr.sh (fastq_dir) (fastq_basename) (loci) (run_dir) (log_file) (SLURM_ARRAY_TASK_ID)


### get arguments

# input
fastq_dir=$1						# fastq path
fastq_basename=$2       			# fastq basename
run_dir=$3           				# output path
loci=$4								# chain
log_file=$5         				# log file for all sample processing progress
SLURM_ARRAY_TASK_ID=$6				# task id
MIXCR_JAR=$7						# path to mixcr .jar file


# output

# executables
# mixcr=/cga/wu/jingsun/software/mixcr-2.1.5/mixcr.jar


### set environments

# source /broad/software/scripts/useuse
# #reuse -q .java-jdk-1.8.0_92-x86-64
# reuse -q .java-jdk-1.8.0_121-x86-64

### start data processing

echo
echo "[ALIGN TCR READS WITH MIXCR]"
echo "-------------------------------------------------"
echo
echo

start_time=$(date +%s)


### make directories

echo ">> MAKE DIRECTORIES"
echo

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ mkdir ${run_dir}mixcr/"
if [[ ! -d ${run_dir}mixcr/ ]]; then
	mkdir ${run_dir}mixcr/
fi
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ mkdir ${run_dir}mixcr/${fastq_basename}_${loci}/"
if [[ -d ${run_dir}mixcr/${fastq_basename}_${loci}/ ]]; then
	rm -r ${run_dir}mixcr/${fastq_basename}_${loci}/
fi
mkdir ${run_dir}mixcr/${fastq_basename}_${loci}/

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ cd ${run_dir}mixcr/${fastq_basename}_${loci}/"
cd ${run_dir}mixcr/${fastq_basename}_${loci}/

echo
echo


### run mixcr

echo ">> RUN MIXCR"
echo
# java -Xmx5g
mixcr_align="java -Xmx5g -jar $MIXCR_JAR align --chains $loci --species hsa --save-description ${fastq_dir}${fastq_basename}_${loci}_L001_R1_001.fastq alignments.vdjca"
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $mixcr_align"
$mixcr_align
echo

mixcr_assemble="java -Xmx5g -jar $MIXCR_JAR assemble --index index_file alignments.vdjca clones.clns"
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $mixcr_assemble"
$mixcr_assemble
echo

mixcr_alignment_pretty="java -Xmx5g -jar $MIXCR_JAR exportAlignmentsPretty alignments.vdjca alignments.pretty.txt"
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $mixcr_alignment_pretty"
$mixcr_alignment_pretty
echo

mixcr_alignment="java -Xmx5g -jar $MIXCR_JAR exportAlignments --preset full -readId -descrR1 -cloneIdWithMappingType index_file alignments.vdjca alignments.txt"
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $mixcr_alignment"
$mixcr_alignment
echo

mixcr_clone="java -Xmx5g -jar $MIXCR_JAR exportClones --chains $loci --preset full -readIds index_file clones.clns clones.txt"
echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $mixcr_clone"
$mixcr_clone
echo



### end data processing

echo "-------------------------------------------------"
echo

wc -l clones.txt

echo "-------------------------------------------------"
echo

end_time=$(date +%s)
elapsed_time="$(expr ${end_time} - ${start_time})"

start_time_format=`date -d @$start_time`
end_time_format=`date -d @$end_time`
elapsed_time_format=`echo | awk -v D=$elapsed_time '{printf "%02d:%02d:%02d\n", D/(60*60), D%(60*60)/60, D%60}'`

echo "ELAPSED TIME: $elapsed_time_format"
echo

echo -e "$SLURM_ARRAY_TASK_ID\t${fastq_basename}_${loci}\t$start_time_format\t$end_time_format\t$elapsed_time_format\tsuccessfully completed" >> ${log_file}

