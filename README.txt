TCR Sequencing Analysis Pipeline

System Requirements:

This pipeline has been tested on MacOS High Sierra Version 10.13.6. It requires the following software dependencies, and has been tested with the versions indicated:

BLAST version NCBI-BLAST-2.2.30 (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.30/)
MiXCR version 2.1.5 (https://github.com/milaboratory/mixcr/releases/tag/v2.1.5)
Java JDK version 8 (https://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html)
GNU parallel version 20180722 (http://git.savannah.gnu.org/cgit/parallel.git/)
Scripts and reference files from the Github repository at https://github.com/julietforman/rhTCRseq

It also requires Python 3 with the following modules:

Biopython (https://biopython.org/wiki/SourceCode)
os
sys
argparse
re
csv
pandas
time
collections
gzip
matplotlib

Installation guide:

BLAST, MiXCR, Java, and GNU Parallel can be installed from the provided URLs. The scripts and reference files specific to this pipeline are available from https://github.com/julietforman/rhTCRseq.git. They can be downloaded via a web browser by navigating to this URL and selecting "Clone or download" > "Download ZIP". Alternatively, the Github repository can be cloned onto the user's machine using the command "git clone https://github.com/julietforman/rhTCRseq.git". Cloning this repository should take no more than a few minutes.

Demo:

A demo dataset is available in the Github repository in the data/sample_run folder. The pipeline can be run on this data according to the instructions provided below. The pipeline should take approximately 10-15 minutes to run on this demo dataset.

The expected output is available in the Github repository in the out/sample_run folder.

Pipeline steps:

This pipeline is designed to analyze data from TCR or target gene sequencing runs. The pipeline comprises several steps:

1: Setup

This step is completed by a call to the script make_index_list.py, which generates two files of information about the run that will be referred to in later steps:
	<ROOT_DIR>out/<RUN_NAME>/index_list.txt
		This file stores information about each well and its corresponding index 1 and index 2.
	<ROOT_DIR>out/<RUN_NAME>/run_info.csv
		This file lists metadata about the run, including the name of the run, paths to the run directory, data directory, and index_list.txt file, and whether the run is single cell, has UMIs, and/or has target genes.

2: Map

This step is completed by a call to the function mapFastq2Well in analyze_tcr.py. It maps each fastq_basename (the beginning of the file name for the fastq files) to a well of the plate, and generates two files:
	<ROOT_DIR>out/<RUN_NAME>/results/fastq_basename4blast.list
		This file contains a list of fastq_basenames to be used as input in future steps.
	<ROOT_DIR>out/<RUN_NAME>/results/fastq_basename2well.list
		This file records the mapping between fastq_basenames and wells, as well as the read count within each well.

3: Blast

This step runs once for each fastq_basename, via calls to blast.sh. This script blasts each read against the relevant primers/amplicons and returns the best hit in order to determine which gene/locus each read came from. It generates one file per blast call with the results of that call. 
	For TCR runs, the outputs for each fastq_basename will be:
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R1_vs_TRV.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R2_vs_TRC.blast

	For target_gene runs, the outputs for each fastq_basename will be:
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R1_vs_TRV.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R1_vs_target_gene_forward_primer.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R2_vs_target_gene_reverse.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R1_vs_target_gene_forward.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R2_vs_TRC.blast
		<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/R2_vs_target_gene_reverse_primer.blast

4: Blast parse

This step runs once for each fastq_basename by calling parse_blast_results.py. It goes through the blast outputs, filters them according to some quality thresholds, and identifies which reads had hits to the TRA locus, TRB locus, or particular target genes. This step generates two output files for each fastq_basename:
	<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/fastq2tcr.txt
		This file contains a mapping between read sequence ids and the locus/target gene to which that read aligns
	<ROOT_DIR>out/<RUN_NAME>/blast/<fastq_basename>/fastq2tcr_count.txt
		This file contains a read count for each locus/target gene.

5: Blast count

This step is run by calling the countBlastResults function in analyze_tcr.py. This step counts the number of reads aligned to each locus/target gene for each well, and generates tables with the total number and the proportion of reads for each locus/target gene for each well:
	<ROOT_DIR>out/<RUN_NAME>/results/read_count.list
		This file lists the total number of reads and number of reads from each locus/target gene for each well.
	<ROOT_DIR>out/<RUN_NAME>/results/read_count.table
		This file contains one table for the total across all genes, one table for each locus/target gene, and one table for the all of the target genes cumulatively. The table with the total lists the total number of reads for each well. Each of the other tables list the proportion of reads in each well that come from that gene.

6: Separate

For each fastq_basename, this step separates the reads aligned to TRA and the reads aligned to TRB into different files. The outputs for each fastq_basename are:
	<ROOT_DIR>out/<RUN_NAME>/fastq/tcr/<fastq_basename>_TRA_L001_R1_001.fastq
		fastq file containing all of the sequences from the fastq_basename R1 file that align to TRA
	<ROOT_DIR>out/<RUN_NAME>/fastq/tcr/<fastq_basename>_TRB_L001_R1_001.fastq
		fastq file containing all of the sequences from the fastq_basename R1 file that align to TRB
Additionally, for runs with UMIs, there will be two additional files:
	<ROOT_DIR>out/<RUN_NAME>/fastq/umi/<fastq_basename>_TRA_umi.txt
		file containing read ids matched to their UMI sequences for TRA
	<ROOT_DIR>out/<RUN_NAME>/fastq/umi/<fastq_basename>_TRB_umi.txt
		file containing read ids matched to their UMI sequences for TRB

7: Mixcr

This step runs the Mixcr program (via the script mixcr.sh) for each locus (TRA and TRB) in each fastq_basename. Mixcr aligns each read to reference V, D, J, and C genes of the TCR receptor, and then assembles the aligned reads into clonotypes. This step generates four output files per locus per fastq_basename:
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/alignments.vdjca
		This is a binary file that contains the alignments of the reads to reference genes
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/alignments.txt
		This is a human readable file that contains the alignments of the reads to reference genes
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/alignments.pretty.txt
		This is another human readable file containing the alignments of the reads, which indicates the positions within the reads of the different genes and regions in the TCR, and shows the read and the reference sequences side by side.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/clones.clns
		This is a binary file that contains the clonotypes assembled from the reads.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/clones.txt
		This is a human readable file that lists the information about each assembled clonotype.

8: Merge TRBV

Some of the V genes at the TRB locus to which Mixcr aligns reads and assigns clonotypes are indistinguishable based on the locations of the primers for those genes. This step edits the clonotypes returned by Mixcr to replace indistinguishable V gene assignments when they appear as the first V gene for that clonotype. Each indistinguishable V gene is replaced by a V gene identifier for the particular group of indistinguishable V genes. There are two groups of indistinguishable V genes: TRBV6-2, TRBV6-3, TRBV6-5, and TRBV6-6 are replaced with TRBV6-2/3/5/6, and TRBV12-3 and TRBV12-4 are replaced with TRBV12-3/4. This step runs once on the TRB clonotypes for each fastq_basename via the merge_TRBV.py script. This script generates two output files for each fastq_basename:
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/clones_pre_combine_TRBV6_TRBV12.txt
		This file is a renamed copy of the original clones.txt file returned by Mixcr, with the original V gene assignments intact.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/clones.txt
		This file is an updated version of clones.txt that has the indistinguishable V genes replaced by the identifier for their group whenever they appear as the first V gene for a particular clonotype.

### FOR BULK RNA RUNS ONLY ###
9: Mixcr UMI

This step runs once for each locus for each fastq_basename via a call to the script count_umi.py. This step identifies the UMIs that go with each read in each clonotype. Within each clonotype, any UMIs that have only one read supporting them are removed, and any clonotypes with no remaining UMIs are eliminated. Then all remaining clonotypes are sorted into groups with the same V and J gene assignment pair. Within each group with matching V and J gene pairs, clonotypes with at least 95% nucleotide sequence identity are clustered together. Within each of these clusters, the clonotypes are divided into groups that have the same length cdr3 amino acid sequence, and for each group, the majority clonotype with the most unique UMIs is identified and the rest of the clonotypes in the group are collapsed into that majority clonotype, which in turn gains all of the reads and UMIs from the clonotypes that are subsumed. Finally, clonotypes that are pseudogenes or ORFs or have only one read in support are deleted. This step generates several output files for each locus for each fastq_basename:
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/collapse_stat.txt
		This file records the total number of reads and UMIs across all clonotypes, and the number of clonotypes present before collapse.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/clone_stat.csv
		This file records the fastq_basename, locius, total reads, total clonotypes, total UMIs, average reads/umi, average UMIs/clone, and percentage of clonotypes collapsed.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>/removed_clones.txt")
		This file records the clonotypes that were deleted after collapse due to being pseudogenes or ORFs or having no supporting UMIs remaining after filtering.
	<ROOT_DIR>out/<RUN_NAME>/results/mixcr/<fastq_basename>_<locus>_clones.txt
		This file records all of the information about the clonotypes after all collapsing and filtering steps.
	<ROOT_DIR>out/<RUN_NAME>/mixcr/<fastq_basename>_<locus>_clones_umi.txt
		This file records the UMI count, clone_count, and clone_fraction of each clonotype.

### FOR BULK RNA RUNS ONLY ###
10: Compare Clonotype

This step compares clonotypes between different replicate wells of the same sample. For each sample, this step collects all the clonotypes from all replicates. Then these clonotypes are sorted into groups with the same J gene and CDR3 sequence pair. Then all the clonotypes in each group are collapsed to the clonotype within that group that has the the highest UMI count. If there is a tie for the highest UMI count, the clonotype with the longest V gene sequence is chosen, as defined by the position of the primer for the V gene of that clonotype. Once these clonotypes have been collapsed, any non-productive clonotypes and clonotypes with a read count below a threshold value are removed. This step generates many output files:
These files contain information about all of the samples:
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/all_clonotype_count_TRA.csv
    	This file is a table that records the number of clonotypes with each TRA V gene in each sample. For each sample, one row is added to the table with the number of clonotypes per V gene in that sample.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/all_clonotype_count_TRB.csv
    	This file is a table that records the number of clonotypes with each TRB V gene in each sample. For each sample, one row is added to the table with the number of clonotypes per V gene in that sample.
These files contain information about a particular sample, and this step generates one of each of the following files for each sample in the run:
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/<sample>_<locus>_clones_count.csv
    	For each clonotype, this file records the number of reads assigned to that clonotype for each replicate well in the sample, as well as the number of replicates in which the clonotype appears, and the total number of reads across all replicates.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/<sample>_<locus>_clones_umi_count.csv
    	For each clonotype, this file records the number of UMIs supporting that clonotype for each replicate well in the sample, as well as the number of replicates in which the clonotype appears, and the total number of UMIs across all replicates.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/<sample>_<locus>_observation_across_replicates.csv
    	This file records the distribution of the number of replicates in which each clonotype appears. For each number of replicates, it records the number of clonotypes that appear in that number of replicates within the sample.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/<sample>_<locus>_clonotype_count.csv
    	This file records the number of clonotypes within the sample that have each of the potential V genes.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr/<sample>_<locus>_clonotype_count.png
    	This figure is a barplot of the number of clonotypes per V gene in the sample.

### FOR SINGLE CELL RNA RUNS ONLY ###
11: Mixcr Parse

This step examines clonotypes for each locus in each well. Within each well, it collapses clonotypes with the same V and J genes based on sequence similarity.


    <ROOT_DIR>out/<RUN_NAME>/results/collapse_stat_all.txt
    	This file contains records of the pairs of clonotypes that were collapsed for each locus and well, as well as the initial number of clonotypes in that well and the number that were collapsed.
    <ROOT_DIR>out/<RUN_NAME>/results/clone_stat_all.txt
    	This file records the number of clonotypes left after collapse and the percentage of clonotypes collapsed for each locus in each well.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr_clonotype_per_well.list
    	This file contains information about the read count, top three clonotypes, and more detailed information about the top clonotype for each locus in each well.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr_clonotype_<locus>_across_well.table")
    	Two versions of this file are generated: one for TRA and once for TRB. There is one row in each file for each well. Each row records the clone fraction of the top clonotype in that well in each other well in the run.
    <ROOT_DIR>out/<RUN_NAME>/results/mixcr_clonotype_recurrence.list
    	This file records each cdr3 sequence (both TRA and TRB) and lists the wells in which that cdr3 is found, as well as the total number of wells in which it appears. It also records each TRA/TRB clonotype pair found and lists the TRA cdr3 sequence and the TRB cdr3 sequence, as well as the number of wells and list of wells in which that pair appears.

Instructions for use:

The required input for the rhTCRseq pipeline includes the sequencing reads in fastq.gz files separated by well, with four files per well. Each filename for each well should begin with the fastq_basename prefix for that well, and the four filenames should end in _L001_R1_001.fastq.gz, _L001_R2_001.fastq.gz, _L001_I1_001.fastq.gz, and _L001_I2_001.fastq.gz. Also required is the SampleSheet.csv file for the run.

1: Install MiXCR-2.1.5, NCBI Blast, Java, and GNU Parallel versions specified above.

2: Set up directory structure:
	choose a ROOT_DIR in which all requirements for the pipeline will go
	create a directory <ROOT_DIR>/scripts (add to this folder all of the files from the Github repository folder rhTCRseq/scripts) 
	create a directory <ROOT_DIR>/data
	create a directory <ROOT_DIR>/out
	create a directory <ROOT_DIR>/blast_database (add to this folder all of the blast database files from the Github repository folder rhTCRseq/blast_database)

3: Prepare to analyze a particular run:
	open the config file <ROOT_DIR>/scripts/config.py and edit the variables in the file to match your run
	create a directory in <ROOT_DIR>/data/ with the same name as RUN_NAME in config.py (place the fastq.gz files from the sequencing run in this folder)
	create a directory in <ROOT_DIR>/out/ with the same name as RUN_NAME in config.py (place SampleSheet.csv for the run into this folder)

After these steps are complete, the directory structure should look like this:

ROOT_DIR
	blast_database
		TRV_primer.fasta
		TRV_primer.fasta.nhr
		TRV_primer.fasta.nin
		TRV_primer.fasta.nog
		TRV_primer.fasta.nsd
		TRV_primer.fasta.nsi
		TRV_primer.fasta.nsq
		TRC_primer.fasta
		TRC_primer.fasta.nhr
		TRC_primer.fasta.nin
		TRC_primer.fasta.nog
		TRC_primer.fasta.nsd
		TRC_primer.fasta.nsi
		TRC_primer.fasta.nsq
		target_gene_primer_forward.fasta
		target_gene_primer_forward.fasta.nhr
		target_gene_primer_forward.fasta.nin
		target_gene_primer_forward.fasta.nog
		target_gene_primer_forward.fasta.nsd
		target_gene_primer_forward.fasta.nsi
		target_gene_primer_forward.fasta.nsq
		target_gene_primer_reverse.fasta
		target_gene_primer_reverse.fasta.nhr
		target_gene_primer_reverse.fasta.nin
		target_gene_primer_reverse.fasta.nog
		target_gene_primer_reverse.fasta.nsd
		target_gene_primer_reverse.fasta.nsi
		target_gene_primer_reverse.fasta.nsq
		target_gene_reverse.fasta
		target_gene_reverse.fasta.nhr
		target_gene_reverse.fasta.nin
		target_gene_reverse.fasta.nog
		target_gene_reverse.fasta.nsd
		target_gene_reverse.fasta.nsi
		target_gene_reverse.fasta.nsq
		target_gene_forward.fasta
		target_gene_forward.fasta.nhr
		target_gene_forward.fasta.nin
		target_gene_forward.fasta.nog
		target_gene_forward.fasta.nsd
		target_gene_forward.fasta.nsi
		target_gene_forward.fasta.nsq
	data
		<RUN_NAME>
			<fastq.gz files>
	out
		<RUN_NAME>
			SampleSheet.csv
	scripts
		blast.sh
		count_umi.py
		separate_fastq.py
		collapse_rules.txt
		parse_blast_results.py
		analyze_tcr.py
		compare_clonotype.py
		merge_TRBV.py
		mixcr.sh
		get_parallel_range.py
		print_description.py
		config.py
		run_pipeline.sh
		make_index_list.py

4: Run the pipeline:
	run <ROOT_DIR>/scripts/run_pipeline.sh







