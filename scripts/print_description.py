#!/usr/bin/env python

"""
Description: Get range (ex: 4 samples in file returns 1 2 3 4) for GNU Parallel

Usage: python get_parallel_range.py

Roger Jin
July, 2018
"""


import argparse
from collections import defaultdict
import csv

nested_dict = lambda: defaultdict(nested_dict)


def main():

    parser = argparse.ArgumentParser(description="Print Step Descriptions")
    parser.add_argument("--step",
                        dest="step",
                        default=None,
                        type=str,
                        help="step to run [map, blast, blast_parse, separate, mixcr, merge_TRBV, mixcr_umi, mixcr_parse]")
    parser.add_argument("--info",
                        dest="run_info_path",
                        default=None,
                        type=str,
                        help="file of run information")

    args = parser.parse_args()
    step = args.step
    info = args.run_info_path
    step_description = {"map": "Step 1: Map FASTQ file to well",
                        "blast": "Step 2.1: Generate BLAST wrapper",
                        "blast_parse": "Step 2.2: Generate BLAST parse wrapper",
                        "blast_count": "Step 2.3: Count BLAST results aligned to TRA, TRB and target genes",
                        "separate": "Step 3: Separate FASTQ files to TRA and TRB",
                        "mixcr": "Step 4.1: Generate MiXCR wrapper",
                        "merge_TRBV": "Step 4.1.5: Merge indistinguishable TRBV regions",
                        "mixcr_umi": "Step 4.2: Generate UMI count wrapper",
                        "mixcr_parse": "Step 4.3: Get clonotype per cell",
                        "master": "One script for all."
                        }
    print("\n" +
          "#" + "-" * 50 + "\n\n" +
          step_description[step] + "\n\n")
    # read run info
    run_list = []
    with open(info) as run_info_file:
        for row in csv.DictReader(run_info_file, skipinitialspace=True):
            run_list = [{k: v for k, v in row.items()}]
    run_info_file.close()

    # process each run
    for run in run_list:
        # print(run)
        run_name = run["run_name"]

        index_path = run["index_path"]
        # print('index_path', index_path)
        target_gene_exist = run["target_gene"]
        umi_exist = run["umi"]
        sc_data = run["single"]
        if "other" in run:
            other_method = run["other"]
        else:
            other_method = "N"

        fastq_dir = run["fastq_dir"]
        run_dir = run["run_dir"]

        print("          [run name]  " + run_name + "\n" +
              "        [index path]  " + index_path + "\n" +
              " [target gene exist]  " + target_gene_exist + "\n" +
              "         [umi exist]  " + umi_exist + "\n" +
              "       [single-cell]  " + sc_data + "\n" +
              "        [input path]  " + fastq_dir + "\n" +
              "       [output path]  " + run_dir + "\n")

    print("#" + "-" * 50 + "\n")


if __name__ == "__main__": main()

