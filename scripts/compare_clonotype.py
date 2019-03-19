#!/usr/bin/env python

"""
Description: compare clonotypes across cells

Usage: python compare_clonotype.py

Jing Sun
Nov, 2017

edited: Juliet Forman Aug. 2018
"""


import os
import sys
import argparse
import re
import csv
import pandas
import time
from collections import defaultdict
import matplotlib
matplotlib.use('agg') #added
import matplotlib.pyplot as plt

# PSEUDOGENES_ORFS = ["TRBV3-2", "TRBV5-7", "TRBV6-7", "TRBV7-1", "TRBV7-5", "TRBV12-1", "TRBV12-2"]

nested_dict = lambda: defaultdict(nested_dict)


## mainly for sctcr data
CLONE_COUNT_CUTOFF = 2


def main():
    
    parser = argparse.ArgumentParser(description="Compare clonotypes")

    parser.add_argument("--run_dir",
                        dest="run_dir",
                        default=None,
                        type=str,
                        help="path to the run")
    parser.add_argument("--clone_dir",
                        dest="clone_dir",
                        default=None,
                        type=str,
                        help="path to with clones")
    parser.add_argument("--compare_list",
                        dest="compare_list",
                        default="wells_to_compare.txt",
                        type=str,
                        help="name of file for comparing")
    parser.add_argument("--collapse_rules_path",
                        dest="collapse_rules_path",
                        default=None,
                        type=str,
                        help="path to collapse rules file")

    args = parser.parse_args()

    run_dir = args.run_dir
    collapse_rules_path = args.collapse_rules_path
    if args.clone_dir == None:
        clone_dir = os.path.join(run_dir, "results", "mixcr")
    else:
        clone_dir = args.clone_dir
    
    compare_list_path = os.path.join(run_dir, args.compare_list)
    all_v_hits = defaultdict(dict)
    
    # read fastq - well mapping
    well2fastq_basename = {}
    map_path = os.path.join(run_dir, "results", "fastq_basename2well.list")
    with open(map_path) as map_file:

        map_data = pandas.read_table(map_file, sep = "\t", header = 0)
        for index, row in map_data.iterrows():
            well2fastq_basename[row["well"]] = row["fastq_basename"]

    map_file.close()


    # read list of wells to compare
    with open(compare_list_path) as compare_list_file:

        compare_list_data = pandas.read_table(compare_list_file, sep = "\t", header = 0)
        for index, row in compare_list_data.iterrows():
            group = row["group"]
            wells = row["wells"].split(",")
            locus = row["locus"].split(",")
            
            # Remove white space in group
            tmp = group.split(" ")
            group = "_".join(tmp)
            
            fastq_basenames = []
            for well in wells:
                if well in well2fastq_basename:
                    fastq_basenames.append(well2fastq_basename[well])
                else:
                    fastq_basenames.append("empty_file")

            for loci in locus:
                all_v_hits[group][loci] = compareClonotype(clone_dir, group, fastq_basenames, loci, collapse_rules_path)
        

        # Combine clone counts
        av_regions = ["TRAV1-1","TRAV1-2","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7","TRAV8-1","TRAV8-2","TRAV8-3","TRAV8-4","TRAV8-6","TRAV9-1","TRAV9-2","TRAV10","TRAV12-1","TRAV12-2","TRAV12-3","TRAV13-1","TRAV13-2","TRAV14DV4","TRAV16","TRAV17","TRAV18","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRAV24","TRAV25","TRAV26-1","TRAV26-2","TRAV27","TRAV29DV5","TRAV30","TRAV34","TRAV35","TRAV36DV7","TRAV38-1","TRAV38-2DV8","TRAV39","TRAV40","TRAV41"]
        bv_regions = ["TRBV2","TRBV3-1","TRBV3-2","TRBV4-1","TRBV4-2","TRBV4-3","TRBV5-1","TRBV5-4","TRBV5-5","TRBV5-6","TRBV5-7","TRBV5-8","TRBV6-1","TRBV6-2","TRBV6-3","TRBV6-4","TRBV6-5","TRBV6-6","TRBV6-7","TRBV6-8","TRBV6-9","TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-5","TRBV7-6","TRBV7-7","TRBV7-8","TRBV7-9","TRBV9","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1","TRBV11-2","TRBV11-3","TRBV12-2","TRBV12-3","TRBV12-4","TRBV12-5","TRBV13","TRBV14","TRBV15","TRBV16","TRBV18","TRBV19","TRBV20-1","TRBV24-1","TRBV25-1","TRBV27","TRBV28","TRBV29-1","TRBV30"]
    
        combine_file_a = open(os.path.join(clone_dir, "all_clonotype_count_TRA.csv"), "w")
        combine_file_a.write("Sample\t%s\n"%("\t".join(av_regions)))
        combine_file_b = open(os.path.join(clone_dir, "all_clonotype_count_TRB.csv"), "w")
        combine_file_b.write("Sample\t%s\n"%("\t".join(bv_regions)))
        
        for sample in all_v_hits:
            for chain in all_v_hits[sample]:
                if chain == "TRA":
                    combine_file_a.write("%s\t%s\n"%(sample, "\t".join(map(str,all_v_hits[sample][chain]))))
                if chain == "TRB":
                    combine_file_b.write("%s\t%s\n"%(sample, "\t".join(map(str,all_v_hits[sample][chain]))))
                        
    compare_list_file.close()


def getLenDict(collapse_rules_path):
    collapse_rules_file = open(collapse_rules_path, 'r')
    out_dict = {}
    for line in collapse_rules_file:
        fields = line.rstrip().split("\t")
        pos = fields[1]
        for v_hit in fields[2:]:
            if "DV" in v_hit:
                parts = v_hit.split("/")
                v_hit = parts[0] + parts[1]
            out_dict[v_hit] = pos

    return out_dict

def getPseudogenesORFs():

    collapse_rules_path = "/cga/wu/juliet/scripts/collapse_rules.txt"
    collapse_rules_file = open(collapse_rules_path, 'r')
    out_list = []
    for line in collapse_rules_file:
        fields = line.rstrip().split("\t")
        pos = fields[1]
        for v_hit in fields[2:]:
            if "DV" in v_hit:
                parts = v_hit.split("/")
                v_hit = parts[0] + parts[1]
            out_list.append(v_hit)

    return out_list
   

def compareClonotype(clone_dir, group, fastq_basenames, loci, collapse_rules_path):

    j_cdr3_collapse_dict = {} #dict of collapsed clonotypes and what to replace them with
    same_j_cdr3_hits = defaultdict(list)
    clonotype_umi_counts = defaultdict(int)
    fastq_basename_umi_counts = {}
    fastq_basename_read_counts = {}
    for fastq_basename in fastq_basenames:
        fastq_basename_umi_counts[fastq_basename] = {}
        fastq_basename_read_counts[fastq_basename] = {}
        clone_path = os.path.join(clone_dir, fastq_basename + "_" + loci + "_clones.txt")
        if not os.path.exists(clone_path):
                continue
        with open(clone_path) as clone_file:
            clone_data = pandas.read_table(clone_file, sep = "\t", header = 0)
            for index, row in clone_data.iterrows():

                clonotype = row["clonotype"]

                # get basic counts for this clone
                clone_count = row["clone_count"]
                clone_fraction = row["clone_fraction"]

                # get allele for v/j/c regions
                v_hit = row["v_hit"]
                j_hit = row["j_hit"]
                c_hit = row["c_hit"]

                # get sequences
                clone_sequence = row["clone_sequence"]
                cdr3 = row["cdr3"]
                umi_count = row["umi_count"]
                umis = row["umis"]

                j_cdr3_hit = "%s_%s"%(j_hit,cdr3)


                same_j_cdr3_hits[j_cdr3_hit].append({"clonotype": clonotype,
                                "clone_count": clone_count,
                                "clone_fraction": clone_fraction,
                                "clone_sequence": clone_sequence,
                                "v_hit": v_hit, "j_hit": j_hit, "c_hit": c_hit,
                                "cdr3": cdr3,
                                "umi_count": umi_count,
                                "umis": umis})

                clonotype_umi_counts[clonotype] += umi_count

                if j_cdr3_hit not in fastq_basename_umi_counts[fastq_basename]:
                    fastq_basename_umi_counts[fastq_basename][j_cdr3_hit] = []

                fastq_basename_umi_counts[fastq_basename][j_cdr3_hit].extend(eval(umis))

                if j_cdr3_hit not in fastq_basename_read_counts[fastq_basename]:
                    fastq_basename_read_counts[fastq_basename][j_cdr3_hit] = 0

                fastq_basename_read_counts[fastq_basename][j_cdr3_hit] += clone_count

                

    len_dict = getLenDict(collapse_rules_path) #get lengths of different v regions

    functional_gene_list = getPseudogenesORFs()

    # identify majority clonotype(s) by umi count
    

    for j_cdr3_hit in same_j_cdr3_hits:
        if len(same_j_cdr3_hits[j_cdr3_hit]) == 1: #continue if only one clonotype has this j_cdr3_hit
            continue

        majorityClonotypes = []
        high_umi_count = 0
        backup = ""
        for clone in same_j_cdr3_hits[j_cdr3_hit]:
            clonotype = clone["clonotype"]
            umi_count = clonotype_umi_counts[clonotype]
            v_hit = clone["v_hit"]
            # if v_hit in PSEUDOGENES_ORFS:
            if v_hit not in functional_gene_list
                print("found pseudogene/orf: ", v_hit)
                backup = clonotype
                continue
            if umi_count > high_umi_count:
                high_umi_count = umi_count
                majorityClonotypes = [clonotype]
            elif umi_count == high_umi_count:
                majorityClonotypes.append(clonotype)

        if majorityClonotypes == []:
            majorityClonotypes.append(backup)

        if len(majorityClonotypes) == 1:
            collapseClonotype = majorityClonotypes[0]

        else:
            collapseClonotype = ""
            low_pos = float('inf')
            for clonotype in majorityClonotypes:
                v_hit = clonotype.split(",")[0]
                pos = int(len_dict[v_hit])
                if pos < low_pos:
                    low_pos = pos
                    collapseClonotype = clonotype

        new_dict = {}
        new_clone_count = 0
        new_clone_fraction = 0
        new_umis = []
        new_umi_count = 0

        #print("collapseClonotype: ", collapseClonotype)
        for clone in same_j_cdr3_hits[j_cdr3_hit]:
            clonotype = clone["clonotype"]
            # print("umis: ", clone["umis"])

            new_clone_count += clone["clone_count"]
            new_clone_fraction += clone["clone_fraction"]
            umis = eval(clone["umis"])
            # print("eval(umis): ", umis)
            # print("first umi: ", umis[0])
            new_umis.extend(umis)
            if clonotype == collapseClonotype:
                new_dict["clone_sequence"] = clone["clone_sequence"]
                new_dict["v_hit"] = clone["v_hit"]
                new_dict["j_hit"] = clone["j_hit"]
                new_dict["c_hit"] = clone["c_hit"]
                new_dict["cdr3"] = clone["cdr3"]
                new_dict["clonotype"] = clonotype
        # print("new_umis: ", new_umis)
        new_umi_count = len(set(new_umis))
        # print("new_umi_count: ", new_umi_count)

        new_dict["clone_count"] = new_clone_count
        new_dict["clone_fraction"] = new_clone_fraction
        new_dict["umi_count"] = new_umi_count


        for clone in same_j_cdr3_hits[j_cdr3_hit]:
            clonotype = clone["clonotype"]
            j_cdr3_collapse_dict[clonotype] = new_dict







    # read clonotype list
    clone_dict = nested_dict()
    for fastq_basename in fastq_basenames:

        clone_path = os.path.join(clone_dir, fastq_basename + "_" + loci + "_clones.txt")
        if not os.path.exists(clone_path):
                continue
        with open(clone_path) as clone_file:

            clone_data = pandas.read_table(clone_file, sep = "\t", header = 0)

            for index, row in clone_data.iterrows():

                clonotype = row["clonotype"]

                if clonotype in j_cdr3_collapse_dict:

                    clone_count = j_cdr3_collapse_dict[clonotype]["clone_count"]
                    clone_fraction = j_cdr3_collapse_dict[clonotype]["clone_fraction"]
                    v_hit = j_cdr3_collapse_dict[clonotype]["v_hit"]
                    j_hit = j_cdr3_collapse_dict[clonotype]["j_hit"]
                    cdr3 = j_cdr3_collapse_dict[clonotype]["cdr3"]
                    umi_count = j_cdr3_collapse_dict[clonotype]["umi_count"]
                    new_clonotype = j_cdr3_collapse_dict[clonotype]["clonotype"]

                    # replace collapsed clonotypes with the correct clonotype info in dataframe
                    clone_data.at[index, "clone_count"] = clone_count
                    clone_data.at[index, "clone_fraction"] = clone_fraction
                    clone_data.at[index, "v_hit"] = v_hit
                    clone_data.at[index, "j_hit"] = j_hit
                    clone_data.at[index, "c_hit"] = c_hit
                    clone_data.at[index, "umi_count"] = umi_count
                    clone_data.at[index, "clonotype"] = new_clonotype

            # re-sort dataframe to reflect new clone_fractions
            clone_data.sort_values(["clone_fraction"], ascending=False)


            for index, row in clone_data.iterrows():

                clonotype = row["clonotype"]

                #skip any repeated clonotypes within the fastq_basename (due to collapsing)
                if clonotype in clone_dict:
                    if fastq_basename in clone_dict[clonotype]:
                        continue

                # only productive clonotypes
                if re.search("\_|\*", row["cdr3"]):
                    continue

                # stop when the clonotype with read count is lower than cutoff
                if row["clone_count"] < CLONE_COUNT_CUTOFF:
                    continue


                clone_dict[row["clonotype"]][fastq_basename] = {"clone_count": row["clone_count"],
                                                                "clone_fraction": row["clone_fraction"],
                                                                "v_hit": row["v_hit"],
                                                                "j_hit": row["j_hit"],
                                                                "cdr3": row["cdr3"],
                                                                "umi_count": row["umi_count"]}


        clone_file.close()



    # print clone count per group
    clone_path = os.path.join(clone_dir, str(group) + "_" + loci + "_clones_count.csv")
    with open(clone_path, "w+") as clone_file:

        clone_file.write(",".join(["v_hit", "j_hit", "cdr3"] + fastq_basenames + ["freq", "count_sum"]) + "\n")
        for clonotype in clone_dict:
            clone_file.write(clonotype)

            v_hit, j_hit, cdr3 = clonotype.split(",")
            j_cdr3_hit = "%s_%s"%(j_hit,cdr3)

            clone_n = 0
            clone_sum = 0
            for fastq_basename in fastq_basenames:
                if fastq_basename not in clone_dict[clonotype]:
                    clone_count = 0
                else:
                    # clone_count = clone_dict[clonotype][fastq_basename]["clone_count"]
                    clone_count = fastq_basename_read_counts[fastq_basename][j_cdr3_hit]

                    clone_n = clone_n + 1
                    clone_sum = clone_sum + clone_count

                clone_file.write("," + str(clone_count))

            clone_file.write("," + str(clone_n) + "," + str(clone_sum) + "\n")

    clone_file.close()


    # print clone umi count per group
    clone_path = os.path.join(clone_dir, str(group) + "_" + loci + "_clones_umi_count.csv")
    v_hits = defaultdict(int)
    clone_n_frequencies = defaultdict(int)
    with open(clone_path, "w+") as clone_file:

        clone_file.write(",".join(["v_hit", "j_hit", "cdr3"] + fastq_basenames + ["freq", "count_sum"]) + "\n")
        for clonotype in clone_dict:
            clone_file.write(clonotype)
            
            v_hit, j_hit, cdr3 = clonotype.split(",")
            j_cdr3_hit = "%s_%s"%(j_hit,cdr3)
            v_hits[v_hit] += 1
            
            clone_n = 0
            clone_sum = 0
            for fastq_basename in fastq_basenames:
                if fastq_basename not in clone_dict[clonotype]:
                    clone_count = 0
                else:
                    #clone_count = clone_dict[clonotype][fastq_basename]["umi_count"]
                    # print("j_cdr3 umis in fastq_basename: ", fastq_basename_umi_counts[fastq_basename][j_cdr3_hit])
                    clone_count = len(set(fastq_basename_umi_counts[fastq_basename][j_cdr3_hit]))

                    clone_n = clone_n + 1
                    clone_sum = clone_sum + clone_count
                    
                clone_file.write("," + str(clone_count))

            clone_file.write("," + str(clone_n) + "," + str(clone_sum) + "\n")
            clone_n_frequencies[clone_n] += 1
    clone_file.close()

    
    clone_n_freq = os.path.join(clone_dir, str(group) + "_" + loci + "_observation_across_replicates.csv")
    with open(clone_n_freq, "w") as n_clone_file:
        n_clone_file.write("How many times clone is observed,Count\n")
        for count in clone_n_frequencies:
            n_clone_file.write("%d,%d\n"%(count, clone_n_frequencies[count]))

    # Make plot
    png = os.path.join(clone_dir, str(group) + "_" + loci + "_clonotype_count.png")
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.25)
    fig.subplots_adjust(top=0.90)
    ax = fig.add_subplot(111)
    # labels for bars
    if loci == "TRA":
        v_regions = ["TRAV1-1","TRAV1-2","TRAV2","TRAV3","TRAV4","TRAV5","TRAV6","TRAV7","TRAV8-1","TRAV8-2","TRAV8-3","TRAV8-4","TRAV8-6","TRAV9-1","TRAV9-2","TRAV10","TRAV12-1","TRAV12-2","TRAV12-3","TRAV13-1","TRAV13-2","TRAV14DV4","TRAV16","TRAV17","TRAV18","TRAV19","TRAV20","TRAV21","TRAV22","TRAV23DV6","TRAV24","TRAV25","TRAV26-1","TRAV26-2","TRAV27","TRAV29DV5","TRAV30","TRAV34","TRAV35","TRAV36DV7","TRAV38-1","TRAV38-2DV8","TRAV39","TRAV40","TRAV41"]
    if loci == "TRB":
        v_regions = ["TRBV2","TRBV3-1","TRBV3-2","TRBV4-1","TRBV4-2","TRBV4-3","TRBV5-1","TRBV5-4","TRBV5-5","TRBV5-6","TRBV5-7","TRBV5-8","TRBV6-1","TRBV6-2/3/5/6","TRBV6-4","TRBV6-7","TRBV6-8","TRBV6-9","TRBV7-2","TRBV7-3","TRBV7-4","TRBV7-5","TRBV7-6","TRBV7-7","TRBV7-8","TRBV7-9","TRBV9","TRBV10-1","TRBV10-2","TRBV10-3","TRBV11-1","TRBV11-2","TRBV11-3","TRBV12-2","TRBV12-3/4","TRBV12-5","TRBV13","TRBV14","TRBV15","TRBV16","TRBV18","TRBV19","TRBV20-1","TRBV24-1","TRBV25-1","TRBV27","TRBV28","TRBV29-1","TRBV30"]
    
    tick_label = v_regions
    x_values = range(1, len(tick_label) + 1)
    y_values = []
    v_hits_for_file = dict()

    for v_region in v_regions:
        if v_region in v_hits:
            v_hits_for_file[v_region] = v_hits[v_region]
            y_values.append(v_hits[v_region])
        else:
            v_hits_for_file[v_region] = 0
            y_values.append(0)

    
    # Make file for V region frequencies
    v_freq_path = os.path.join(clone_dir, str(group) + "_" + loci + "_clonotype_count.csv")
    with open(v_freq_path, "w") as v_freq_file:
        v_freq_file.write("v_hit,Count\n")
        for v_hit in sorted(v_hits_for_file):
            v_freq_file.write("%s,%s\n"%(v_hit, v_hits_for_file[v_hit]))
        v_freq_file.write("Total,%d\n"%(sum(v_hits_for_file.values())))


    ax.text(0.98, 1, '%d Unique Clonotypes'%(sum(v_hits.values())),
        verticalalignment='top', horizontalalignment='right',
        transform=ax.transAxes, fontsize=11)

    
    plt.bar(x_values, y_values, color = ['blue'], align='center', alpha = 0.7)
    plt.xticks(x_values, tick_label, fontsize=9, rotation='vertical')
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    # naming the x-axis
    plt.xlabel('All Productive %s Genes'%(loci))
    # naming the y-axis
    plt.ylabel('Number of Unique Clonotypes')


    plt.savefig(png)

    plt.close(fig)




    return y_values
if __name__ == "__main__": main()