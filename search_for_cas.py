#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

from tqdm import tqdm

import os
import argparse
import Bio.SeqIO

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
	
	parser.add_argument("-p", '--primers', type=str, help='fasta con primers')
	parser.add_argument("-v", '--verbose', action='store_true')
	parser.add_argument("-r", '--reference', help="fasta reference", default="/ref/MN996528.fna")
	parser.add_argument("-a", '--annotation', help="fasta reference", default="/app/snpEff/data/covid19/genes.gbk")
	
	args = parser.parse_args()
	return args

## Look for spacers
def look_for_spacers (seq, options):
	run_vmatch2
	run_fuzznuc
	extract_spacers

def run_vmatch2 ():
	cmd = f'vmatch2 ... {sample}.fasta"'
	exec(cmd, verbose=verbose)

def run_fuzznuc ():

def extract_spacers ():

## Look for CAS
def look_for_cas ():
	run_hmmsearch
	find_best_cas

def run_hmmsearch ():

def find_best_cas ():

## Wrap everything


## Main
def main():
	parser = argparse.ArgumentParser(description='')
	options = parse_option(parser)
	in_fasta = bpio.read(option.input, "fasta")
	for seq in in_fasta:
		look_for_spacers(seq, options)
