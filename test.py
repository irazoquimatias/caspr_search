#!/usr/bin/env python3
from tqdm import tqdm
import argparse
import os
import subprocess as sp
import Bio.SeqIO as bpio
from termcolor import colored
from glob import glob
from Bio.Seq import Seq
from glob import glob
from tqdm import tqdm
import collections
import re

seq1 = "AAAMMMNNNNNNNAAA"
length = len(seq1)
new_length = length
Ns = 0
n_iter = re.finditer('n+', seq1.lower())
for match in n_iter:
    if match.span()[0] == 0 or match.span()[1]==length:
        new_length = new_length - match.span()[1] + match.span()[0]
        continue 
    Ns += match.span()[1] - match.span()[0]
print1 = "Seq1\tNs: "+str(Ns)+"\tLen: "+str(length)+"\tNew_Len: "+str(new_length)
#print (print1)

seq2 = "NNNMMMNNNNOOONNN"
length = len(seq2)
new_length = length
Ns = 0
n_iter = re.finditer('n+', seq2.lower())
for match in n_iter:
    if match.span()[0] == 0 or match.span()[1]==length:
        print (new_length, match.span()[1], match.span()[0])
        new_length = new_length - match.span()[1] + match.span()[0]
        continue 
    Ns += match.span()[1] - match.span()[0]
print2 = "Seq2\tNs: "+str(Ns)+"\tLen: "+str(length)+"\tNew_Len: "+str(new_length)
#print (print2)

seq3 = "A-ATT-CCC-GG"
seq4 = "AAATTTCCCGGG"
i = 0
while i < len(seq3):
	report_codon = ""
	for a in range(i, i+3):
		if seq3[a] != seq4[a]:
			report_codon += seq4[a].upper()
		else:
			report_codon += seq4[a].lower()
	print (i, i+2, report_codon)
	i += 3
