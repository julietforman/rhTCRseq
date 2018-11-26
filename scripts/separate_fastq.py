#!/usr/bin/env python

"""
Description: separate FASTQ files

Usage: python separate_fastq.py

Jing Sun
Nov, 2017
"""


import os
import sys
import argparse
import re
import csv
import pandas
import time
from collections import defaultdict
import gzip

nested_dict = lambda: defaultdict(nested_dict)



def main():

    parser = argparse.ArgumentParser(description="Separate FASTQ files")

    parser.add_argument("--run_dir",
                        dest="run_dir",
                        default=None,
                        type=str,
                        help="path to the run")
    parser.add_argument("--fastq_dir",
                        dest="fastq_dir",
                        default=None,
                        type=str,
                        help="path to fastq files")
    parser.add_argument("--fastq_basename",
                        dest="fastq_basename",
                        default=None,
                        type=str,
                        help="fastq file basename")
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
    fastq_dir = args.fastq_dir
    fastq_basename = args.fastq_basename
    umi_exist = args.umi_exist
    log_path = args.log_path
    sge_task_id = args.sge_task_id


    # start data processing
    start_time = time.time()
    print("\n" +
          "[SEPARATE FASTQ FILES]" + "\n" +
          "-" * 50 + "\n\n")


    # read fastq id - tcr mapping
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] READ SEQ ID - TCR MAPPING")

    id2tcr = {}
    map_path = os.path.join(run_dir, "blast", fastq_basename, "fastq2tcr.txt")
    with open(map_path) as map_file:
        map_data = pandas.read_table(map_file, sep = "\t", header = 0)
        for index, row in map_data.iterrows():
            id2tcr[row["seq_id"]] = row["loci"]

    map_file.close()
    # read raw fastq, read1
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] READ RAW FASTQ READ1")

    fastq_tcr = defaultdict(list)
    fastq_path = os.path.join(fastq_dir, fastq_basename + "_L001_R1_001.fastq.gz")
    with gzip.open(fastq_path, "tr") as fastq_file:
        fastq_block = ""
        loci = ""
        for count, row in enumerate(fastq_file, start=1):
            if count % 4 == 1:
                fastq_block = ""
                loci = ""

                id = re.search("(@.+?)\s", row).group(1)
                if id in id2tcr:
                    loci = id2tcr[id]

            fastq_block += row

            if count % 4 == 0:
                if loci == "":
                    continue
                fastq_tcr[loci].append(fastq_block)
    fastq_file.close()

    # separate fastq by tra and trb
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] SEPARATE FASTQ BY TRA AND TRB")

    fastq_tcr_dir = os.path.join(run_dir, "fastq", "tcr")
    if not os.path.exists(fastq_tcr_dir):
        os.makedirs(fastq_tcr_dir)

    fastq_tcr_path = {"TRA": os.path.join(fastq_tcr_dir, fastq_basename + "_TRA_L001_R1_001.fastq"),
                      "TRB": os.path.join(fastq_tcr_dir, fastq_basename + "_TRB_L001_R1_001.fastq")}

    for loci in fastq_tcr_path:
        if len(fastq_tcr[loci]) == 0:
            continue

        fastq_path = fastq_tcr_path[loci]
        with open(fastq_path, "w") as fastq_file:
            fastq_file.write("".join(fastq_block for fastq_block in fastq_tcr[loci]))
        fastq_file.close()


    if umi_exist == "Y":

        # read umi, read2
        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] READ UMI READ2")

        fastq_umi_dir = os.path.join(run_dir, "fastq", "umi")
        if not os.path.exists(fastq_umi_dir):
            os.makedirs(fastq_umi_dir)

        fastq_umi = defaultdict(list)
        fastq_path = os.path.join(fastq_dir, fastq_basename + "_L001_R2_001.fastq.gz")
        with gzip.open(fastq_path, "rt") as fastq_file:
            id = ""
            umi = ""
            loci = ""
            for count, row in enumerate(fastq_file, start = 1):
                if count % 4 == 1:
                    loci = ""

                    id = re.search("(@.+?)\s", row).group(1)
                    if id in id2tcr:
                        loci = id2tcr[id]

                if count % 4 == 2:
                    umi = row[0:7]

                if count % 4 == 0:
                    if loci == "":
                        continue
                    fastq_umi[loci].append(id + "\t" + umi)

        fastq_file.close()

        # print umi by tra and trb
        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] PRINT UMI BY TRA AND TRB")

        fastq_umi_path = {"TRA": os.path.join(fastq_umi_dir, fastq_basename + "_TRA_umi.txt"),
                          "TRB": os.path.join(fastq_umi_dir, fastq_basename + "_TRB_umi.txt")}

        for loci in fastq_umi_path:
            if len(fastq_umi[loci]) == 0:
                continue

            umi_path = fastq_umi_path[loci]
            with open(umi_path, "w") as umi_file:
                umi_file.write("\n".join(fastq_umi for fastq_umi in fastq_umi[loci]))
            umi_file.close()


    # end data processing
    end_time = time.time()
    elapsed_time = end_time - start_time

    start_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(start_time))
    end_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(end_time))
    elapsed_time_format = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print("\n" +
          "-" * 50 + "\n" +
          "ELAPSED TIME: " + elapsed_time_format + "\n")

    with open(log_path, "a") as log_file:
        log_file.write(sge_task_id + "\t" + fastq_basename + "\t" +
                       start_time_format + "\t" + end_time_format + "\t" + elapsed_time_format + "\t" +
                       "successfully completed" + "\n")
    log_file.close()


if __name__ == "__main__": main()