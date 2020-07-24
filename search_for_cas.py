#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

from tqdm import tqdm

import os
import subprocess
import argparse
import shutil
from Bio import SeqIO

# import subprocess as sp
# from termcolor import colored
# from glob import glob
# from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
# from glob import glob
# from tqdm import tqdm
# import collections
# import re

## Parse options
def parse_option(parser):
    parser.add_argument('-i', '--input', type=str, help='Archivo FASTA a analizar')
    parser.add_argument('-o', '--output', type=str, help='Directorio de salida')
    parser.add_argument('-v', '--verbose', action='store_true')
    
    args = parser.parse_args()
    return args

## Look for spacers
def look_for_spacers (seq, options):
    run_vmatch2 (seq, options)
#	run_fuzznuc
#	extract_spacers

def run_vmatch2 (seq, options):
    work_dir = '/tmp/casper_search/current_contig/'
    verbose = options.verbose
    if not os.path.exists(work_dir):
        os.makedirs (work_dir)
    os.chdir(work_dir)
        
    SeqIO.write(seq, 'contig.fasta', 'fasta')
    seq_id = seq.id
    
    cmd = 'mkvtree -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1 -db contig.fasta'
    subprocess.run(cmd.split())
    
    cmd = 'vmatch2 -l 23 25 60 -s leftseq -evalue 0.001 -absolute -noevalue -noscore \
          -noidentity -sort sd -best 10000 -selfun /media/2tb_Viejo/CRISPRCasFinder/sel392v2.so \
          contig.fasta'
    salida = subprocess.getoutput(cmd)
    print (seq_id, len(salida))
    print (salida)

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
    current_dir = os.getcwd()
    if not os.path.exists(options.output):
        os.makedirs(options.output)
    if not os.path.exists("/tmp/casper_search"):
        os.makedirs ("/tmp/casper_search")
    with open(options.input, "rU") as in_fasta:
        for seq in SeqIO.parse(in_fasta, "fasta"):
            look_for_spacers(seq, options)
    shutil.rmtree("/tmp/casper_search") 

if __name__ == "__main__":
    main()
