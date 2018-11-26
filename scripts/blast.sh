#!/bin/bash

##updated for not broad server

# run blast
# usage: blast.sh (fastq_dir) (fastq_basename) (run_dir) (target_gene) (log_file) (sge_task_id)

#set the path to the nucleotide databases made from the primer fasta files
# blast_database_dir=/cga/wu/juliet/blast_database

### get arguments

# input
fastq_dir=$1						# fastq path
fastq_basename=$2       			# fastq basename
run_dir=$3           				# output path
target_gene=$4						# if there are target genes
log_file=$5         				# log file for all sample processing progress
sge_task_id=$6						# task id
blast_database_dir=$7				# path to directory where the blast databases are found


### start data processing

echo
echo "[RUN BLAST]"
echo "-------------------------------------------------"
echo
echo

start_time=$(date +%s)


### make directories

echo ">> MAKE DIRECTORIES"
echo

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ mkdir ${run_dir}fasta/"
if [[ ! -d ${run_dir}fasta/ ]]; then
	mkdir ${run_dir}fasta/
fi

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ mkdir ${run_dir}blast/"
if [[ ! -d ${run_dir}blast/ ]]; then
	mkdir ${run_dir}blast/
fi

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ mkdir ${run_dir}blast/${fastq_basename}/"
if [[ -d ${run_dir}blast/${fastq_basename}/ ]]; then
	rm -rf ${run_dir}blast/${fastq_basename}/
fi
mkdir ${run_dir}blast/${fastq_basename}/


echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ cd ${run_dir}blast/${fastq_basename}/"
cd ${run_dir}blast/${fastq_basename}/

echo
echo


### convert fastq to fasta

echo ">> CONVERT FASTQ TO FASTA"
echo

reads=(R1 R2)
for which_read in ${reads[@]}
do
	echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $convert_fastq2fasta"

	gunzip -c ${fastq_dir}${fastq_basename}_L001_${which_read}_001.fastq.gz | paste - - - -  |awk '{print ">"$1 ; print $3;}' > ${run_dir}fasta/${fastq_basename}_L001_${which_read}_001.fasta
	echo
done

echo
echo

### blast primers

echo ">> BLAST PRIMERS"
echo

# blast std format
## qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore

# blast R1 to TCR variable primer
run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R1_001.fasta -query_loc 1-25 -strand plus -task blastn-short -db ${blast_database_dir}TRV_primer.fasta -out R1_vs_TRV.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

eval $run_blast
echo

# blast R2 to TCR constant primer
run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R2_001.fasta -query_loc 1-35 -strand plus -task blastn-short -db ${blast_database_dir}TRC_primer.fasta -out R2_vs_TRC.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

eval $run_blast
echo

if [[ $target_gene == "Y" ]]; then


	# blast R1 to target gene forward primer
	run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R1_001.fasta -query_loc 1-25 -strand plus -task blastn-short -db ${blast_database_dir}target_gene_primer_forward.fasta -out R1_vs_target_gene_forward_primer.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

	echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

	eval $run_blast
	echo

	# blast R2 to target gene reverse primer
	run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R2_001.fasta -query_loc 1-35 -strand plus -task blastn-short -db ${blast_database_dir}target_gene_primer_reverse.fasta -out R2_vs_target_gene_reverse_primer.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

	echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

	eval $run_blast
	echo

	#blast R2 to target gene full length reverse amplicon
	run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R2_001.fasta -query_loc 1-220 -strand plus -task blastn -db ${blast_database_dir}target_gene_reverse.fasta -out R2_vs_target_gene_reverse.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

	echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

	eval $run_blast
	echo

	# blast R1 to target gene full length forward amplicon
	run_blast="blastn -query ${run_dir}fasta/${fastq_basename}_L001_R1_001.fasta -query_loc 1-220 -strand plus -task blastn -db ${blast_database_dir}target_gene_forward.fasta -out R1_vs_target_gene_forward.blast -outfmt \"6 std qlen slen sstrand\" -max_target_seqs 1"

	echo "`date '+[%a %b %d %H:%M:%S %Z %Y]'` $ $run_blast"

	eval $run_blast
	echo

fi

echo
echo


### end data processing

echo "-------------------------------------------------"
echo

end_time=$(date +%s)
elapsed_time="$(expr ${end_time} - ${start_time})"

start_time_format=`date -d @$start_time`
end_time_format=`date -d @$end_time`
elapsed_time_format=`echo | awk -v D=$elapsed_time '{printf "%02d:%02d:%02d\n", D/(60*60), D%(60*60)/60, D%60}'`

echo "ELAPSED TIME: $elapsed_time_format"
echo

echo -e "$SLURM_ARRAY_TASK_ID\t$fastq_basename\t$start_time_format\t$end_time_format\t$elapsed_time_format\tsuccessfully completed" >> ${log_file}

