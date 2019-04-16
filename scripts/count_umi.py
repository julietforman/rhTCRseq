#!/usr/bin/env python

"""
Description: get unique UMI count per clonetype

Usage: python count_umi.py

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
# import os.path
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

# PSEUDOGENES_ORFS = ["TRBV3-2", "TRBV5-7", "TRBV6-7", "TRBV7-1", "TRBV7-5", "TRBV12-1", "TRBV12-2"]

# functional_gene_list = getPseudogenesORFs(collapse_rules_path)

nested_dict = lambda: defaultdict(nested_dict)

# def write_clone_collapse(old_hit, hit, sf1):
#     sf1.write(",".join(["clone_count", "clone_fraction", "v_hit", "j_hit", "c_hit",
#                     "clone_sequence", "cdr3", "umi_count", "clonotype"]) + "\n")
#     sf1.write(",".join([str(old_hit["clone_count"]),
#                     str(old_hit["clone_fraction"]),
#                     old_hit["v_hit"],
#                     old_hit["j_hit"],
#                     old_hit["c_hit"],
#                     old_hit["clone_sequence"],
#                     old_hit["cdr3"],
#                     str(old_hit["umi_count"]),
#                     old_hit["clonotype"]]) + "\n")
#     sf1.write(",".join([str(hit["clone_count"]),
#                     str(hit["clone_fraction"]),
#                     hit["v_hit"],
#                     hit["j_hit"],
#                     hit["c_hit"],
#                     hit["clone_sequence"],
#                     hit["cdr3"],
#                     str(hit["umi_count"]),
#                     hit["clonotype"]]) + "\n")

# def remove_umis(umi_list):
#     umi_sorted = sorted(umi_list,key=umi_list.count,reverse=True)
#     used = set()
#     umi_unique = [x for x in umi_sorted if x not in used and (used.add(x) or True)]
    
#     umi_count = 0
#     umis_to_keep = list()
#     for umi in umi_unique:
#         if umi_list.count(umi) > 1:
#            umis_to_keep.append(umi)
#            umi_count += 1
    
#     umis_included = [x for x in umi_list if x in umis_to_keep]
#     return umi_count, umis_included

# def count_umi(umi_list, sf):
#     max = 1
#     umi_sorted = sorted(umi_list,key=umi_list.count,reverse=True)
#     used = set()
#     umi_unique = [x for x in umi_sorted if x not in used and (used.add(x) or True)]
    
#     umi_count = 0
#     total_freq = 0
#     low_freq_umis = list()
#     high_freq_umis = list()
#     true_umi = list()
#     umi_read_counts = dict()
    
#     for umi in umi_unique:
#         freq = umi_list.count(umi) / float(len(umi_list))
#         total_freq += freq
#         umi_read_counts[umi] = umi_list.count(umi)
        
#         if (umi_count == 0 and total_freq > 0.10):
#             umi_count +=1
#             high_freq_umis.append(umi)
#             true_umi.append(umi)
#         else:
#             if (total_freq > 0.10):
#                 low_freq_umis.append(umi)
#             else:
#                 umi_count +=1
#                 high_freq_umis.append(umi)
#                 true_umi.append(umi)
        
#     for low_umi in low_freq_umis:
#         collaps = True
#         for high_umi in high_freq_umis:
#             collaps = True
#             n = 0
#             for c1, c2 in zip(low_umi, high_umi):
#                 if c1 != c2:
#                     n += 1
#                     if n > max:
#                         collaps = False
#                         break

#             if collaps == True:
#                 break
            
#         if collaps == False:
#             umi_count +=1
#             high_freq_umis.append(low_umi)
#             true_umi.append(low_umi)
#         else:
#             sf.write("Umi %s collapsed with %s\n"%(low_umi, high_umi))
#     umis_included = [x for x in umi_list if x in true_umi]
        
#     return umi_count,umis_included

def initClusterList(clonotype_list):
    '''input is list of dictionaries which each hold info about one clonotype. Output is a list of clusters. Each cluster
    is a tuple with two elements. The first is a list of clone_ids, and the second is a set of umis represented in all of
    the clone_ids in the cluster. At initialization, each clonotype is in its own cluster.'''
    out_list = []
    for clonotype in clonotype_list:
        clone_id = clonotype["clone_id"]
        clone_id_list = [clone_id]
        umi_set = set(clonotype["umis"])
        out_list.append((clone_id_list, umi_set))
    return out_list

def mergeClusters(cluster_list, index_list):
    '''takes as input a cluster_list and a list of indices of clusters to be merged together. Returns an updated cluster list
    where the clusters at the positions in index_list have been merged into one.'''
    out_list = []
    newClusterCloneIDs = []
    newClusterUmiSet = set()
    for i in range(len(cluster_list)):
        if i in index_list:
            clone_id_list = cluster_list[i][0]
            umi_set = cluster_list[i][1]
            newClusterCloneIDs.extend(clone_id_list)
            newClusterUmiSet.update(umi_set)
        else:
            out_list.append(cluster_list[i])
    if newClusterCloneIDs != []:
        out_list.append((newClusterCloneIDs, newClusterUmiSet))
    return out_list

def collapseClusters(cluster_list, clone_dict, functional_gene_list):

    for cluster in cluster_list:
        clone_id_list = cluster[0]
        # clone_umi_set = cluster[1]
        if len(clone_id_list) == 1: #if there is only one clone in the cluster, add it to output dict and continue
            continue

        cdr3_len_dict = {}

        for clone_id in clone_id_list: # group each clone by the aa length of its cdr3
            len_cdr3 = len(clone_dict[clone_id]["cdr3"])
            if len_cdr3 not in cdr3_len_dict:
                cdr3_len_dict[len_cdr3] = []
            cdr3_len_dict[len_cdr3].append(clone_id)

        for len_cdr3 in cdr3_len_dict: # collapse clones with same length cdr3
            # find the clonotype with the most unique UMIs (ties broken arbitrarily)
            majorityClonotype = ""
            backup = ""
            high_umi_count = 0
            for clone_id in cdr3_len_dict[len_cdr3]:
                umi_count = clone_dict[clone_id]["umi_count"]
                v_hit = clone_dict[clone_id]["v_hit"]
                # if v_hit in PSEUDOGENES_ORFS:
                if v_hit not in functional_gene_list:
                    backup = clone_id
                    continue
                if umi_count > high_umi_count:
                    high_umi_count = umi_count
                    majorityClonotype = clone_id


            if majorityClonotype == "":
                majorityClonotype = backup

            # merge all clonotypes in cluster to the majority clonotype and add merged clone to output dict
            clone_dict = mergeClonotypes(majorityClonotype, cdr3_len_dict[len_cdr3], clone_dict)

    return clone_dict



def mergeClonotypes(majority_clone_id, clone_id_list, clone_dict):
    new_dict = clone_dict[majority_clone_id]

    new_clone_count = 0
    new_clone_fraction = 0
    new_umis = []
    new_umi_count = 0

    for clone_id in clone_id_list:
        new_clone_count += clone_dict[clone_id]["clone_count"]
        new_clone_fraction += clone_dict[clone_id]["clone_fraction"]
        new_umis.extend(clone_dict[clone_id]["umis"])
        if clone_id != majority_clone_id:
            del clone_dict[clone_id]

    new_umi_count = len(set(new_umis))

    new_dict["clone_count"] = new_clone_count
    new_dict["clone_fraction"] = new_clone_fraction
    new_dict["umis"] = new_umis
    new_dict["umi_count"] = new_umi_count

    clone_dict[majority_clone_id] = new_dict

    return clone_dict


# def getLenDict():
#     collapse_rules_path = "/cga/wu/juliet/scripts/collapse_rules.txt"
#     collapse_rules_file = open(collapse_rules_path, 'r')
#     out_dict = {}
#     for line in collapse_rules_file:
#         fields = line.rstrip().split("\t")
#         pos = fields[1]
#         for v_hit in fields[2:]:
#             if "DV" in v_hit:
#                 parts = v_hit.split("/")
#                 v_hit = parts[0] + parts[1]
#             out_dict[v_hit] = pos

#     return out_dict


def getPseudogenesORFs(collapse_rules_path):

    # collapse_rules_path = "/cga/wu/juliet/scripts/collapse_rules.txt"
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


def getHammingDist(umi1, umi2):
    if umi1 == umi2:
        return 0

    count = 0
    for i in range(len(umi1)):
        if umi1[i] != umi2[i]:
            count += 1

    return count

def mergeUmiClusters(clusters, index_list):
    outClusters = []
    newCluster = []
    for i in range(len(clusters)):
        if i not in index_list:
            outClusters.append(clusters[i])
        else:
            newCluster.extend(clusters[i])
    outClusters.append(newCluster)
    return outClusters


# def collapseUmiClusters(clusters, umi_list):
#     out_list = []
#     for cluster in clusters:

#         #get umi with highest count, and total count of all umis in cluster
#         majorityUmi, totalCount = getMajorityUmi(cluster, umi_list)

#         #add majority umi to list once for each time a umi in the cluster appeared in the original umi list
#         for i in range(totalCount):
#             out_list.append(majorityUmi)

#     # filter out umis with only one read before returning
#     return filterUmis(out_list)

# def filterUmis(umi_list):
#     out_list = []
#     for umi in umi_list:
#         if umi_list.count(umi) > 1:
#             out_list.append(umi)
#     return out_list

# def getMajorityUmi(cluster, umi_list):
#     majorityUmi = ""
#     highCount = 0
#     totalCount = 0
#     for umi in cluster:
#         count = umi_list.count(umi)
#         totalCount += count
#         if count > highCount:
#             majorityUmi = umi
#             highCount = count

#     return majorityUmi, totalCount

# def getCloneList(cloneList):
#     out_clone_list = []
#     for clonotype in cloneList:
#         clone_id = clonotype["clone_id"]
#         out_clone_list.append(clone_id)
#     return out_clone_list


def main():

    parser = argparse.ArgumentParser(description="Count unique UMI per clonotype")

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
    parser.add_argument("--loci",
                        dest="loci",
                        default=None,
                        type=str,
                        help="loci")
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
    parser.add_argument("--collapse_identity",
                        dest="collapse_identity",
                        default=0.95,
                        type=float,
                        help="CDR3 nuc identity threshold for collaple of clones")
    parser.add_argument("--collapse_rules_path",
                        dest="collapse_rules_path",
                        default=None,
                        type=str,
                        help="path to collapse rules file")

    args = parser.parse_args()
    
    collapse_identity = args.collapse_identity
    run_dir = args.run_dir
    fastq_basename = args.fastq_basename
    loci = args.loci
    log_path = args.log_path
    sge_task_id = args.sge_task_id
    collapse_rules_path = args.collapse_rules_path

    functional_gene_list = getPseudogenesORFs(collapse_rules_path)


    # start data processing
    start_time = time.time()
    print("\n" +
          "[GET UNIQUE UMI COUNT]" + "\n" +
          "-" * 50 + "\n\n")
    
    # read umi per read
    print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] read umi per read")
    
    umi_path = os.path.join(run_dir, "fastq", "umi", fastq_basename + "_" + loci + "_umi.txt")
    if os.path.isfile(umi_path):
        with open(umi_path) as umi_file:
            umi_data = pandas.read_table(umi_file, sep="\t", header=None)
        umi_file.close()
        # read clonotype
        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] get umi count per clonotype")
    
        clone_path = os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci, "clones.txt")
        
        os.makedirs(os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci), exist_ok=True)
        stat_file1 = os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci, "collapse_stat.txt")
        stat_file2 = os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci, "clone_stat.csv")
        # same_vj_hits = defaultdict(list)
        clone_dict = defaultdict(list)
        read_counts = 0
        # total_umis = 0
        sf1 = open(stat_file1, "w")
        sf2 = open(stat_file2, "w")
        sf2.write("Fastq_basename,loci,Total reads,Total clones,Total UMIs,Reads/umi,UMIs/clone,% clones collapsed\n")
        with open(clone_path) as clone_file:
            clone_data = pandas.read_table(clone_file, sep = "\t", header = 0, keep_default_na = False)
            for index, row in clone_data.iterrows():
                clone_id = row["cloneId"]

                # get basic counts for this clone
                clone_count = row["cloneCount"]
                clone_fraction = row["cloneFraction"]
    
                # get one allele per v/j/c region
                v_hit = re.search("(%s.+?)\*" % loci, row["allVHitsWithScore"]).group(1)
                j_hit = re.search("(%s.+?)\*" % loci, row["allJHitsWithScore"]).group(1)
                vj_hit = "%s_%s"%(v_hit,j_hit)
                
                c_hit = ""
                if row["allCHitsWithScore"]:
                    c_hit = re.search("(%s.+?)\*" % loci, row["allCHitsWithScore"]).group(1)
    
                # get sequences
                clone_sequence = row["clonalSequence"]
                cdr3 = row["aaSeqCDR3"]

                j_cdr3_hit = "%s_%s"%(j_hit,cdr3)
    
                # process umi
                reads = row["reads"]
                if re.search(",", str(reads)):
                    read_list = reads.split(",")
                    umi_all = umi_data[1][list(map(int, read_list))]
                    umi_list = umi_all.tolist()
                    umis = umi_list
                    umi_count = len(set(umi_list))
                else:
                    umi_count = 1
                    read_list = [reads]
                    umi_all = umi_data[1][list(map(int, read_list))]
                    umis = umi_all.tolist()

                
                # total_umis += umi_count
                
                
                # same_vj_hits[vj_hit].append({"clone_count": clone_count,
                #                 "clone_fraction": clone_fraction,
                #                 "clone_sequence": clone_sequence,
                #                 "v_hit": v_hit, "j_hit": j_hit, "c_hit": c_hit,
                #                 "cdr3": cdr3,
                #                 "umi_count": umi_count,
                #                 "clonotype": "\"" + v_hit + "," + j_hit + "," + cdr3 + "\"",
                #                 "umis": umis, "clone_id": clone_id})

                clone_dict[clone_id] = {"clone_count": clone_count,
                                "clone_fraction": clone_fraction,
                                "clone_sequence": clone_sequence,
                                "v_hit": v_hit, "j_hit": j_hit, "c_hit": c_hit,
                                "cdr3": cdr3,
                                "umi_count": umi_count,
                                "clonotype": "\"" + v_hit + "," + j_hit + "," + cdr3 + "\"",
                                "umis": umis, "clone_id": clone_id}
    
        clone_file.close()
        clones = len(clone_dict.keys())
        
        # Compare and collaps clonal types based on similiarty in CDR3
        if clones != 0:
            # sf1.write("\nClones collapsed based on CDR3 amino acids sequence:\n\n")

            #remove umis with only one read
            remove_clone_list = []
            for clone_id in clone_dict:
                umi_list = clone_dict[clone_id]["umis"]
                remove_umi_list = []
                for umi in umi_list:
                    if umi_list.count(umi) == 1:
                        remove_umi_list.append(umi)
                for umi in remove_umi_list:
                    umi_list.remove(umi)
                if len(umi_list) == 0:
                    remove_clone_list.append(clone_id)
                clone_dict[clone_id]["umis"] = umi_list
                clone_dict[clone_id]["umi_count"] = len(set(umi_list))
            for clone_id in remove_clone_list:
                del clone_dict[clone_id]


            same_vj_hits = defaultdict(list)
            for clone_id in clone_dict:
                v_hit = clone_dict[clone_id]["v_hit"]
                j_hit = clone_dict[clone_id]["j_hit"]
                vj_hit = "%s_%s"%(v_hit,j_hit)

                same_vj_hits[vj_hit].append(clone_dict[clone_id])


            #for clonotypes with same V and J: cluster clonotypes with shared UMIs
            for vj_combo in same_vj_hits:
                if len(same_vj_hits[vj_combo]) == 1: #continue if only one clonotype has this vj_combo
                    continue

                cluster_list = initClusterList(same_vj_hits[vj_combo])

                for clonotype in same_vj_hits[vj_combo]:
                    cluster_index_list_to_merge = [] # list to store indices of clusters to be merged
                    clone_id = clonotype["clone_id"]
                    # umi_set = set(clonotype["umis"])
                    clone_sequence_1 = clonotype["clone_sequence"]

                    for i in range(len(cluster_list)):
                        # cluster_umi_set = cluster_list[i][1]
                        # if not umi_set.isdisjoint(cluster_umi_set):
                        #     cluster_index_list_to_merge.append(i)
                        clone_list = cluster_list[i][0]
                        for clone in clone_list:
                            clone_sequence_2 = clone_dict[clone]["clone_sequence"]
                            alignment_score = pairwise2.align.globalxx(clone_sequence_1, clone_sequence_2, score_only = True)
                            
                            if len(clone_sequence_1) >= len(clone_sequence_2):
                                identity = (alignment_score/len(clone_sequence_1))
                            else:
                                identity = (alignment_score/len(clone_sequence_2))

                            if identity >= collapse_identity:
                                cluster_index_list_to_merge.append(i)
                                break

                    cluster_list = mergeClusters(cluster_list, cluster_index_list_to_merge)

                clone_dict = collapseClusters(cluster_list, clone_dict, functional_gene_list)
            

            # for each clonotype, collapse umis within one hamming distance, remove umis with only one read
            # for clone_id in clone_dict:

            #     umi_list = clone_dict[clone_id]["umis"]
            #     umi_set = set(umi_list)
            #     clusters = [[umi] for umi in set(umi_list)]

            #     for umi1 in umi_set:
            #         cluster_index_list_to_merge = [] # list of indices of umi clusters to merge
            #         for i in range(len(clusters)):
            #             cluster = clusters[i]
            #             for umi2 in cluster:
            #                 dist = getHammingDist(umi1, umi2)
            #                 if dist <= 1:
            #                     cluster_index_list_to_merge.append(i)
            #                     continue
            #         clusters = mergeUmiClusters(clusters, cluster_index_list_to_merge)

            #     umi_list = collapseUmiClusters(clusters, umi_list)

            #     clone_dict[clone_id]["umis"] = umi_list
            #     clone_dict[clone_id]["umi_count"] = len(set(umi_list))

            # remove clones with no umis after filtering, remove clones from pseudogenes or ORFs
            remove_clone_list = []
            removedClonesPath = os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci, "removed_clones.txt")
            removedClonesFile = open(removedClonesPath, 'w')
            removedClonesFile.write("v_hit\tj_hit\tcdr3\treason removed\n")
            for clone_id in clone_dict:
                v_hit = clone_dict[clone_id]["v_hit"]
                j_hit = clone_dict[clone_id]["j_hit"]
                cdr3 = clone_dict[clone_id]["cdr3"]
                umi_list = clone_dict[clone_id]["umis"]

                #if only one read or if V hit is pseudogene/orf, put clone in remove list
                if len(umi_list) == 1: #if only one left for this clone after filtering (note: each umi in umi_list represents one read)
                    removedClonesFile.write(v_hit + "\t" + j_hit + "\t" + cdr3 + "\tonly 1 read\n")
                    remove_clone_list.append(clone_id)
                elif v_hit in PSEUDOGENES_ORFS:
                    removedClonesFile.write(v_hit + "\t" + j_hit + "\t" + cdr3 + "\tpseudogene/ORF\n")
                    remove_clone_list.append(clone_id)

            removedClonesFile.close()

            # delete removed clones from clone_dict
            for clone_id in remove_clone_list:
                del clone_dict[clone_id]

    
        # print clonetype with unique umi count
        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] print umi count per clonotype")
    
        os.makedirs(os.path.join(run_dir, "results", "mixcr"), exist_ok=True)
        clone_path = os.path.join(run_dir, "results", "mixcr", fastq_basename + "_" + loci + "_clones.txt")
        with open(clone_path, "w") as clone_file:
            clone_file.write("\t".join(["clone_count", "clone_fraction", "v_hit", "j_hit", "c_hit",
                                        "clone_sequence", "cdr3", "umi_count", "clonotype", "umis"]) + "\n")
            umi_dict = {}
            if clones != 0:
                newlist = sorted(clone_dict.values(), key=lambda k: k['clone_fraction'], reverse=True)
                collapsed_clone_list = newlist
                for clone in collapsed_clone_list:
                    clone_file.write("\t".join([str(clone["clone_count"]),
                                                str(clone["clone_fraction"]),
                                                clone["v_hit"],
                                                clone["j_hit"],
                                                clone["c_hit"],
                                                clone["clone_sequence"],
                                                clone["cdr3"],
                                                str(clone["umi_count"]),
                                                clone["clonotype"],
                                                str(clone["umis"])]) + "\n")
                    
                    # Save umi, clone count and clone fraction
                    umi_dict[clone["cdr3"]] = "%d\t%d\t%f"%(clone["umi_count"], clone["clone_count"], clone["clone_fraction"])
                    
                # Write stat files
                total_umi_count = sum(k['umi_count'] for k in clone_dict.values())
                read_counts = sum(k['clone_count'] for k in clone_dict.values())
                no_clones = len(clone_dict.keys())
                if total_umi_count > 0 and no_clones > 0 and clones > 0:
                    sf2.write("%s,%s,%d,%d,%d,%0.2f,%0.2f,%0.2f\n"%(fastq_basename, loci, read_counts, no_clones, total_umi_count, int(read_counts)/float(total_umi_count), int(total_umi_count)/float(no_clones), (int(clones - len(clone_dict.keys()))/float(clones))*100))
                
                sf1.write("Total read count:%d\n"%(read_counts))
                sf1.write("Total number of UMIs across all clones:%d\n"%(total_umi_count))
                sf1.write("Total number of clones before collapse:%d\n"%(clones))

        clone_file.close()
        
        # print unique umi count per clonotype
        print("[" + time.strftime("%a %b %d %H:%M:%S %Z %Y") + "] print umi count per clonotype")
    
        umi_count_path = os.path.join(run_dir, "mixcr", fastq_basename + "_" + loci, "clones_umi.txt")
        with open(umi_count_path, "w") as umi_count_file:
            umi_count_file.write("cloneId\tumi_count\tclone_count\tclone_fraction\n")
            umi_count_file.write("\n".join(str(id) + "\t" + str(umi_dict[id]) for id in umi_dict))
    
        umi_count_file.close()
        
        sf1.close()
        sf2.close()
    
    # end data processing
    end_time = time.time()
    elapsed_time = end_time - start_time

    start_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(start_time))
    end_time_format = time.strftime("%a %b %d %H:%M:%S %Z %Y", time.localtime(end_time))
    elapsed_time_format = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))

    print("\n" +
          "-" * 50 + "\n" +
          "ELAPSED TIME: " + elapsed_time_format + "\n")
    
    if os.path.isfile(umi_path):
        with open(log_path, "a") as log_file:
            log_file.write(sge_task_id + "\t" + fastq_basename + "\t" +
                       start_time_format + "\t" + end_time_format + "\t" + elapsed_time_format + "\t" +
                       "successfully completed" + "\n")
    else:
        with open(log_path, "a") as log_file:
            log_file.write(sge_task_id + "\t" + fastq_basename + "\t" +
                       start_time_format + "\t" + end_time_format + "\t" + elapsed_time_format + "\t" +
                       "successfully completed - no UMIs" + "\n")
    log_file.close()
if __name__ == "__main__": main()