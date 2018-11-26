#!/usr/bin/env python

from collections import defaultdict
import sys
import os
from config import *


path = ROOT_DIR + "out"
file = "%s/%s/SampleSheet.csv"%(path, RUN_NAME)
out = "%s/%s/index_list.txt"%(path, RUN_NAME)
bulk = 1
no_name = 1
samples = defaultdict(list)

o = open(out, "w")

info = 0
with open(file, "r") as f:
   for line in f:
      if line.startswith("Sample_ID"):
         info = 1
         o.write("well\tIndex1\tIndex2\n")
         continue
      
      if info == 1:
         if ',' in line:
            fields = line.split(",")
         elif ';' in line:
            fields = line.split(";")
         well = fields[0]
         idx1 = fields[5]
         idx2 = fields[7]
         
         o.write("%s\t%s\t%s\n"%(well, idx1, idx2))
         
         if no_name == 0:
            sample = well[1:]
         else:
            sample = fields[1]
         samples[sample].append(well)

o.close()

# Make wells to compare if bulk
compare_file = "%s/%s/wells_to_compare.txt"%(path, RUN_NAME)
with open(compare_file, "w") as cf:
   cf.write("group\twells\tlocus\n")
   
   for sample in sorted(samples.keys()):
      cf.write("%s\t%s\tTRA,TRB\n"%(sample,",".join(samples[sample])))


#make run_info.csv file
run_info_path = "%s/%s/run_info.csv"%(path, RUN_NAME)
run_info_file = open(run_info_path, "w")
run_info_file.write("run_name,run_dir,fastq_dir,index_path,target_gene,umi,single\n")
run_info_file.write(RUN_NAME + ",")
run_info_file.write(path + "/" + RUN_NAME + ",")
run_info_file.write(ROOT_DIR + "data/" + RUN_NAME + ",")
run_info_file.write(path + "/" + RUN_NAME + "/index_list.txt,")
run_info_file.write(TARGET_GENE + ",")
run_info_file.write(UMI + ",")
run_info_file.write(SINGLE)

if __name__ == "__main__":
   print(RUN_NAME)
