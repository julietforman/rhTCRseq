#!/usr/bin/env python

"""
Description: parse BLAST results

Usage: python parse_blast_result.py

Jing Sun
Sep, 2017
"""


import os
import sys
import argparse
import re
import csv
import pandas
import time
from collections import defaultdict

nested_dict = lambda: defaultdict(nested_dict)


# default parameters
TRV_PIDENT_CUTOFF = 95
TRV_LENGTH_CUTOFF = 15
TRV_EVALUE_CUTOFF = 0.001

TRC_PIDENT_CUTOFF = 95
TRC_LENGTH_CUTOFF = 15
TRC_EVALUE_CUTOFF = 0.001

TG_PIDENT_CUTOFF = 95
TG_LENGTH_CUTOFF = 150
TG_EVALUE_CUTOFF = 0.001



def main():

    parser = argparse.ArgumentParser(description="Parse BLAST results")

    parser.add_argument("--run_dir",
                        dest="run_dir",
                        default=None,
                        type=str,
                        help="path to the run")
    parser.add_argument("--fastq_basename",
                        dest="fastq_basename",
                        default=None,
                        type=str,
                        help="fastq file basename")
    parser.add_argument("--target_gene_exist",
                        dest="target_gene_exist",
                        default=None,
                        type=str,
                        help="if there are target genes in library")
    parser.add_argument("--umi_exist",
                        dest="umi_exist",
                        default=None,
                        type=str,
                        help="if there are umi in library")
    parser.add_argument("--log_path",
                        dest="log_path",
                        default=None,
                        type=str,
                        help="log file")
    parser.add_argument("--sge_task_id",
                        dest="sge_task_id",
                        default=None,
                        type=str,
                        help="array task id")

    args = parser.parse_args()

    run_dir = args.run_dir
    fastq_basename = args.fastq_basename
    target_gene_exist = args.target_gene_exist
    umi_exist = args.umi_exist
    log_path = args.log_path
    sge_task_id = args.sge_task_id


    # start data processing
    start_time = time.time()
    print("\n" +
          "[PARSE BLAST RESULTS]" + "\n" +
          "-" * 50 + "\n\n")


    # read blast result
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] read tcr blast results")

    blast_result_path = os.path.join(run_dir, "blast", fastq_basename, "R1_vs_TRV.blast")
    id2trv, trv2id = parseBlastResult(blast_result_path, "N", TRV_PIDENT_CUTOFF, TRV_LENGTH_CUTOFF, TRV_EVALUE_CUTOFF)

    blast_result_path = os.path.join(run_dir, "blast", fastq_basename, "R2_vs_TRC.blast")
    id2trc, trc2id = parseBlastResult(blast_result_path, umi_exist, TRC_PIDENT_CUTOFF, TRC_LENGTH_CUTOFF, TRC_EVALUE_CUTOFF)


    # find reads with consistent trv and trc
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] find reads with consistent trv and trc blast results")

    trv_id = set(id2trv.keys())
    trc_id = set(id2trc.keys())

    id2tcr = {}
    tcr2id = nested_dict()
    tcr2id["TRA"] = {}
    tcr2id["TRB"] = {}
    #for id in trv_id & trc_id:
    for id in trv_id:
        v = id2trv[id][:1]

        loci = "TR" + v
        id2tcr[id] = loci
        tcr2id[loci][id] = ""

    tcr_id = set(id2tcr.keys())

    # print fastq id to tcr mapping
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] print read - tcr blast results mapping")

    map_path = os.path.join(run_dir, "blast", fastq_basename, "fastq2tcr.txt")
    with open(map_path, "w") as map_file:
        map_file.write("\t".join(["seq_id", "loci"]) + "\n")
        for id in tcr_id:
            map_file.write(id + "\t" + id2tcr[id] + "\n")
    map_file.close()

    # print read count aligned to each tcr allele
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] print tcr read counts" + "\n")

    count_path = os.path.join(run_dir, "blast", fastq_basename, "fastq2tcr_count.txt")
    with open(count_path, "w") as count_file:
        count_file.write("\t".join(["loci", "read_count"]) + "\n")
        for loci in tcr2id:
            count_file.write(loci + "\t" + str(len(tcr2id[loci])) + "\n")
    count_file.close()

    log_file = open(log_path, "a")

    # if there are target gene
    if target_gene_exist == "Y":

        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] read target gene blast results")

        blast_result_forward_path = os.path.join(run_dir, "blast", fastq_basename, "R1_vs_target_gene_forward.blast")
        id2forward, forward2id = parseBlastResult(blast_result_forward_path, "N", TG_PIDENT_CUTOFF, TG_LENGTH_CUTOFF, TG_EVALUE_CUTOFF)

        target_gene2reads = {}

        for read in id2forward.keys():
            gene = id2forward[read]
            if gene in target_gene2reads:
                    target_gene2reads[gene].append(read)
            else:
                target_gene2reads[gene] = [read]

        count_path = os.path.join(run_dir, "blast", fastq_basename, "fastq2tg_count.txt")
        with open(count_path, 'w') as count_file:
            count_file.write("\t".join(["target_gene", "read_count"]) + "\n")
            for gene in target_gene2reads:
                count_file.write(gene + "\t" + str(len(target_gene2reads[gene])) + "\n")



    # end data processing
    end_time = time.time()
    elapsed_time = end_time - start_time

    start_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(start_time))
    end_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(end_time))
    elapsed_time_format = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print("\n" +
          "-" * 50 + "\n" +
          "ELAPSED TIME: " + elapsed_time_format + "\n")

    log_file.write(sge_task_id + "\t" + fastq_basename + "\t" +
                       start_time_format + "\t" + end_time_format + "\t" + elapsed_time_format + "\t" +
                       "successfully completed" + "\n")
    log_file.close()

def name2primer(map_path):
    name2forward_primer = {}
    name2reverse_primer = {}
    with open(map_path) as map_file:
        map_data = pandas.read_table(map_file, sep="\t", header=0)
        for index, row in map_data.iterrows():
            name = row["name"]
            forward_primer_name = row["forward primer name"]
            reverse_primer_name = row["reverse primer name"]
            name2forward_primer[name] = forward_primer_name
            name2reverse_primer[name] = reverse_primer_name

    return name2forward_primer, name2reverse_primer


def parseBlastResult(blast_result_path, umi_exist, pident_cutoff, length_cutoff, evalue_cutoff):

    id2gene = {}
    gene2id = nested_dict()

    blast_header = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                    "qlen", "slen", "sstrand"]

    with open(blast_result_path) as blast_result_file:
        blast_result = pandas.read_table(blast_result_file, sep = "\t", header = None, names = blast_header)
        for index, row in blast_result.iterrows():
            id = row["qseqid"]
            gene = row["sseqid"]

            pident = row["pident"]
            length = row["length"]

            qstart = row["qstart"]
            sstart = row["sstart"]

            evalue = row["evalue"]

            if pident < pident_cutoff:
                continue
            if evalue > evalue_cutoff:
                continue
            if length < length_cutoff:
                continue

            id2gene[id] = gene
            gene2id[gene][id] = ""

    return id2gene, gene2id

if __name__ == "__main__": main()