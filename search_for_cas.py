#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

from tqdm import tqdm

import os
import argparse
from Bio import SeqIO

#import subprocess as sp
#from termcolor import colored
#from glob import glob
#from Bio.Seq import Seq
#from Bio.Alphabet import generic_dna
#from glob import glob
#from tqdm import tqdm
#import collections
#import re

## Parse options
def parse_option(parser):
    parser.add_argument("-i", '--input', type=str, help='Archivo FASTA a analizar')
    parser.add_argument("-o", '--output', type=str, help='Directorio de salida')
    
    parser.add_argument("-v", '--verbose', action='store_true')
    
    args = parser.parse_args()
    return args

## Look for spacers
def look_for_spacers (seq, options):
    run_vmatch2 (seq, options)
#	run_fuzznuc
#	extract_spacers

def run_vmatch2 (seq, options):
    os.makedirs (f'/tmp/casper_search/current_contig/')
    os.chdir(f'/tmp/casper_search/current_contig/')
    SeqIO.write(seq, f'contig.fasta', "fasta")
    seq_id = seq.id
    cmd = f' mkvtree2 -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1 -db contig.fasta'
    exec(cmd, verbose=verbose)
    cmd = f'vmatch2 -l 23 25 60 -s leftseq -evalue 0.001 -absolute -noevalue -noscore \
          -noidentity -sort sd -best 10000 -selfun ~/Programas/CRISPRCasFinder/sel392v2.so \
          contig.fasta > {options.out}/{seq_id}_vmatch2.txt'
    exec(cmd, verbose=verbose)

#def run_fuzznuc ():

#def extract_spacers ():

## Look for CAS
#def look_for_cas ():
#	run_hmmsearch
#	find_best_cas

#def run_hmmsearch ():

#def find_best_cas ():

## Wrap everything

## Main
def main():
    parser = argparse.ArgumentParser(description='')
    options = parse_option(parser)
    os.makedirs(options.output)
    os.makedirs (f'/tmp/casper_search')
    with open(options.input, "rU") as in_fasta:
        for seq in SeqIO.parse(in_fasta, "fasta"):
            look_for_spacers(seq, options)
