#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO
import collections
import time
import argparse
from multiprocessing import Pool

def parse_option(parser):
    parser.add_argument('-i', '--input', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, required=True)
    parser.add_argument('-t', '--threads', type=str, default=1,)
    args = parser.parse_args()
    return args

def calc(record):
    length = len(str(record.seq))
    frequencies = collections.Counter(str(record.seq))
    GC = frequencies['C'] + frequencies['G']
    return (record.id, length, GC)

def __resolve_and_merge_results(results_list):
    resolved = []
    try:
        for element in results_list:
            if isinstance(element, mp.pool.ApplyResult):
                element = element.get()
                print (element)
            if isinstance(element, list):
                for sub_element in element:
                    resolved.append(sub_element)
            else:
                resolved.append(element)
    except TypeError:
        resolved = results_list
    return resolved

parser = argparse.ArgumentParser(description='')
options = parse_option(parser)

contigs = []
report = ''
now = time.process_time()

with open(options.input) as fd, Pool(int(options.threads)) as pool:
    gc_for_sequence = pool.map(
        calc,
        (record for record in SeqIO.parse(fd, 'fasta')),
        chunksize=1000,
        )

time_stamp = time.process_time() - now
msg = 'Script took ' + str(round(time_stamp, 4)) + ' seconds'
print (msg)
print (gc_for_sequence)
