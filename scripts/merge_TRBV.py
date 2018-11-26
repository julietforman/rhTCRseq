#!/usr/bin/env python


import os
import argparse
import time

import sys
import re
import csv
import pandas


def main():

    parser = argparse.ArgumentParser(description="merge indistinguishable V segments")

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
    parser.add_argument("--loci",
                        dest="loci",
                        default=None,
                        type=str,
                        help="loci")

    args = parser.parse_args()

    run_dir = args.run_dir
    fastq_basename = args.fastq_basename
    sge_task_id = args.sge_task_id
    log_path = args.log_path
    loci = args.loci

    if loci == "TRA":
        print("skipping TRA")
        return

    # start data processing
    start_time = time.time()
    print("\n" +
          "[PARSE BLAST RESULTS]" + "\n" +
          "-" * 50 + "\n\n")


    # read clones file
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] read clones.txt file")

    clone_path = os.path.join(run_dir, "mixcr", fastq_basename + "_TRB", "clones.txt")

    clone_file = open(clone_path, 'r')
    clone_lines = clone_file.readlines()
    clone_file.close()


    # rename original clones file
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] rename original clones file to clones_pre_combine_TRBV6_TRBV12.txt")

    old_clone_path = os.path.join(run_dir, "mixcr", fastq_basename + "_TRB", "clones_pre_combine_TRBV6_TRBV12.txt")
    os.rename(clone_path, old_clone_path)


    # merge indistinguishable V hits
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] write clones.txt with indistinguishable V hits merged")

    clone_file = open(clone_path, 'w')
    clone_file.write(clone_lines[0].rstrip("\n"))
    for line in clone_lines[1:]:
        fields = line.rstrip("\n").split("\t")
        clone_file.write("\n" + "\t".join(fields[:5]) + "\t")

        v_hits = fields[5].split(",")
        found6 = False
        found12 = False
        v_hit_out = []
        for v_hit in v_hits:
            v_label, v_tail = v_hit.split("*")[:2]
            if v_label in ["TRBV6-2", "TRBV6-3", "TRBV6-5", "TRBV6-6"]:
                if found6:
                    continue
                else:
                    v_label = "TRBV6-2/3/5/6"
                    v_hit_out.append("*".join([v_label, v_tail]))
                    found6 = True
            elif v_label in ["TRBV12-3", "TRBV12-4"]:
                if found12:
                    continue
                else:
                    v_label = "TRBV12-3/4"
                    v_hit_out.append("*".join([v_label, v_tail]))
                    found12 = True
            else:
                v_hit_out.append(v_hit)

        clone_file.write(",".join(v_hit_out) + "\t")
        clone_file.write("\t".join(fields[6:]))
    clone_file.close()

    # end data processing
    end_time = time.time()
    elapsed_time = end_time - start_time

    start_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(start_time))
    end_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(end_time))
    elapsed_time_format = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print("\n" +
          "-" * 50 + "\n" +
          "ELAPSED TIME: " + elapsed_time_format + "\n")

    log_file = open(log_path, "a")
    log_file.write(sge_task_id + "\t" + fastq_basename + "\t" +
                       start_time_format + "\t" + end_time_format + "\t" + elapsed_time_format + "\t" +
                       "successfully completed" + "\n")
    log_file.close()


if __name__ == "__main__": main()






























