#!/bin/bash
#SBATCH -c 20 #Number of cores
#SBATCH -N 1
#SBATCH -t 15:00:00 # Runtime in minutes
#SBATCH --mem=20000 #Memory per node in MB (see also --mem-per-cpu)
#SBATCH -p medium
#SBATCH -o /Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/master/master.out
#SBATCH -J master_sample_data
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jforman@broadinstitute.org
### define parameters

# batch job parameters
source `which env_parallel.bash`
env_parallel --citation 'will cite'
echo ____________________
echo
echo start setup
BLAST_DATABASE_DIR=/Users/jforman/Documents/WuLab/TCR_protocol_files/blast_database/
MIXCR_JAR=/Users/jforman/bin/mixcr-2.1.5/mixcr.jar
fastq_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/data/sample_data/
mixcr_fastq_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/fastq/tcr/
run_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/
target_gene_exist=N
umi_exist=Y
sc_data=N
run_info=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/run_info.csv
>${run_dir}logs/master/master.out

# sample arguments
blast_fastq_list_path=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/results/fastq_basename4blast.list
mixcr_fastq_list_path=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/results/fastq_basename4mixcr.list
src_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/scripts/
echo setup completed
echo
echo ____________________
echo
echo start map
python3 ${src_dir}analyze_tcr.py --info $run_info --step map
sample_range=`python ${src_dir}get_parallel_range.py --fastq_list_path $blast_fastq_list_path`
blast_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/blast/

blast_log_file=${blast_log_dir}array_jobs.log
blast_parse_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/blast_parse/

blast_parse_log_file=${blast_parse_log_dir}array_jobs.log
separate_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/separate/

separate_log_file=${separate_log_dir}array_jobs.log
mixcr_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/mixcr/

mixcr_log_file=${mixcr_log_dir}array_jobs.log
merge_TRBV_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/merge_TRBV/

merge_TRBV_log_file=${mixcr_log_dir}array_jobs.log
mixcr_umi_log_dir=/Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data/logs/mixcr_umi/

mixcr_umi_log_file=${mixcr_umi_log_dir}array_jobs.log
if [[ ! -e $blast_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $blast_log_file
fi
if [[ ! -e $blast_parse_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $blast_parse_log_file
fi
if [[ ! -e $separate_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $separate_log_file
fi
if [[ ! -e $mixcr_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $mixcr_log_file
fi
if [[ ! -e $merge_TRBV_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $merge_TRBV_log_file
fi
if [[ ! -e $mixcr_umi_log_file ]]; then
    echo -e "slurm_no	sample_id	start_time	end_time	elapsed_time	status" > $mixcr_umi_log_file
fi
echo map completed
echo
echo ____________________
echo
echo start blast
python3 ${src_dir}print_description.py --info $run_info --step blast
env_parallel 'blast_fastq_basename=`awk "NR=={}" $blast_fastq_list_path | awk '\''{print $1}'\''` && bash ${src_dir}blast.sh "$fastq_dir" "$blast_fastq_basename" "$run_dir" "$target_gene_exist" "$blast_log_file" {} "$BLAST_DATABASE_DIR" &> ${blast_log_dir}${blast_fastq_basename}.out' ::: $sample_range
wait
echo blast completed
echo
echo ____________________
echo
echo start blast parse
python3 ${src_dir}print_description.py --info $run_info --step blast_parse
env_parallel 'blast_parse_fastq_basename=`awk "NR=={}" $blast_fastq_list_path | awk '\''{print $1}'\''` && python3 ${src_dir}parse_blast_results.py --run_dir $run_dir --fastq_basename $blast_parse_fastq_basename --target_gene_exist $target_gene_exist --umi_exist $umi_exist --log_path $blast_parse_log_file --sge_task_id {} &> ${blast_parse_log_dir}${blast_parse_fastq_basename}.out' ::: $sample_range
wait
echo blast parse completed
echo
echo ____________________
echo
echo start blast count
python3 ${src_dir}analyze_tcr.py --info $run_info --step blast_count
wait
echo blast count completed
echo
echo ____________________
echo
echo start separate
python3 ${src_dir}print_description.py --info $run_info --step separate
env_parallel 'separate_fastq_basename=`awk "NR=={}" $blast_fastq_list_path | awk '\''{print $1}'\''` && python3 ${src_dir}separate_fastq.py --run_dir $run_dir --fastq_dir $fastq_dir --fastq_basename $separate_fastq_basename --umi_exist $umi_exist --log_path $separate_log_file --sge_task_id {} &> ${separate_log_dir}${separate_fastq_basename}.out' ::: $sample_range
wait
echo separate completed
echo
echo ____________________
echo
echo start mixcr
python3 ${src_dir}print_description.py --info $run_info --step mixcr
mixcr_sample_range=`python ${src_dir}get_parallel_range.py --fastq_list_path $mixcr_fastq_list_path`
env_parallel 'mixcr_fastq_basename=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $1}'\''` && loci=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $2}'\''` && echo $mixcr_fastq_basename && bash ${src_dir}mixcr.sh "$mixcr_fastq_dir" "$mixcr_fastq_basename" "$run_dir" "$loci" "$mixcr_log_file" {} "$MIXCR_JAR" &> ${mixcr_log_dir}${mixcr_fastq_basename}.out' ::: $mixcr_sample_range
wait
echo mixcr completed
echo
echo ____________________
echo
echo start merge_TRBV
python3 ${src_dir}print_description.py --info $run_info --step merge_TRBV
env_parallel 'merge_TRBV_fastq_basename=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $1}'\''` && loci=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $2}'\''` && echo $merge_TRBV_fastq_basename && python3 ${src_dir}merge_TRBV.py --run_dir $run_dir --fastq_basename $merge_TRBV_fastq_basename --log_path $merge_TRBV_log_file --sge_task_id {} --loci $loci &> ${mixcr_log_dir}${mixcr_fastq_basename}.out' ::: $mixcr_sample_range
wait
echo merge_TRBV completed
echo
echo ____________________
echo
if [[ $umi_exist == "Y" ]]; then
echo start mixcr umi
python3 ${src_dir}print_description.py --info $run_info --step mixcr_umi
env_parallel 'mixcr_umi_fastq_basename=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $1}'\''` && loci=`awk "NR=={}" $mixcr_fastq_list_path | awk '\''{print $2}'\''` && python3 ${src_dir}count_umi.py --run_dir $run_dir --fastq_basename $mixcr_umi_fastq_basename --loci $loci --log_path $mixcr_umi_log_file --sge_task_id {} &> ${mixcr_umi_log_dir}${mixcr_umi_fastq_basename}.out' ::: $mixcr_sample_range
wait
echo mixcr umi completed
echo
echo ____________________
echo
fi
if [[ $sc_data == "N" ]]; then
echo start compare clonotype
python3 ${src_dir}compare_clonotype.py --run_dir /Users/jforman/Documents/WuLab/TCR_protocol_files/out/sample_data --collapse_rules_path /Users/jforman/Documents/WuLab/TCR_protocol_files/scripts/collapse_rules.txt
wait
echo compare clonotype completed
echo
echo ____________________
echo
fi
if [[ $sc_data == "Y" ]]; then
echo start mixcr parse
python3 ${src_dir}print_description.py --info $run_info --step mixcr_parse
python3 ${src_dir}analyze_tcr.py --info $run_info --step mixcr_parse
wait
echo mixcr parse completed
echo
echo ____________________
echo
fi
echo 'PIPELINE FINISHED :) :) :)'
echo
