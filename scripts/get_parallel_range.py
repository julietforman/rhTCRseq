#!/usr/bin/env python

"""
Description: Get range (ex: 4 samples in file returns 1 2 3 4) for GNU Parallel

Usage: python get_parallel_range.py

Roger Jin
July, 2018
"""


import argparse
from collections import defaultdict

nested_dict = lambda: defaultdict(nested_dict)


def count_samples(fastq_list_path):
    count = 0
    with open(fastq_list_path) as f:
        for line in f:
            if line!='':
                count+=1

    return count

def get_parallel_indices(n):
    #if n = 4 returns 1 2 3 4 - purpose is to feed into gnu parallel
    s = ''
    for i in range(n):
        s+=str(i+1) + ' '
    return s

def get_parallel_range(fastq_list_path):
    count = count_samples(fastq_list_path)
    return get_parallel_indices(count)

def main():

    parser = argparse.ArgumentParser(description="Get String for Parallel")
    parser.add_argument("--fastq_list_path",
                        dest="list_path",
                        default=None,
                        type=str,
                        help="path to the list")

    args = parser.parse_args()
    list_path = args.list_path
    print(get_parallel_range(list_path))

if __name__ == "__main__": main()

