#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import os
import subprocess
import argparse
import shutil
from Bio import SeqIO
import re

## Parse options
def parse_option(parser):
    parser.add_argument('-i', '--input', type=str, help='Archivo FASTA a analizar')
    parser.add_argument('-o', '--output', type=str, help='Directorio de salida')
    parser.add_argument('-s', '--so', type=str, help='Archivo so requerido por vmatch')
    parser.add_argument('-m', '--mismatch', type=int, default=0, help='Numero de mismatches permitidos en direct repeat')

    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()
    return args

## Look for spacers
def look_for_spacers (record, options):
    (positions, dr) = run_vmatch2 (record, options)
    extract_spacers (record, positions, options)

def run_vmatch2 (record, options):
    work_dir = '/tmp/casper_search/current_contig/'
    verbose = options.verbose
    if not os.path.exists(work_dir):
        os.makedirs (work_dir)
    os.chdir(work_dir)

    SeqIO.write(record, 'contig.fasta', 'fasta')

    cmd = 'mkvtree2 -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1 -db contig.fasta'
    subprocess.run(cmd.split(' '))

    cmd = 'vmatch2 -l 23 25 60 -s leftseq -evalue 0.001 -absolute -noevalue -noscore \
          -noidentity -sort sd -best 10000 -selfun %s \
          contig.fasta' % (options.so)
    salida = subprocess.getoutput(cmd)

    best_score = 0
    best_positions = []
    best_dr = ''
    for l in salida.split('\n'):
        l = l.rstrip()
        if re.match('^>', l) or not l:
            continue
        (score, positions) = run_fuzznuc(l, options)
        if score > best_score:
            best_score = score
            best_positions = positions
            best_dr = l
    return (best_positions, best_dr)

def run_fuzznuc (pattern, options):
    score = 0
    positions = []
    cmd = 'fuzznuc -sequence contig.fasta -pattern %s -pmismatch %s -outfile fuzznuc.txt' \
    % (pattern, options.mismatch)
    subprocess.run(cmd.split(' '), stderr = subprocess.PIPE)
    with open('fuzznuc.txt', 'rU') as fuzz_report:
        for line in fuzz_report:
            if re.match('^\s+\d', line):
                columns = line.split()
                columns[4] = 0 if columns[4] == '.' else columns[4]
                score = score + len(columns[5]) - int(columns[4])
                positions.append((columns[0],columns[1]))
    return (score, positions)

def extract_spacers (record, positions, options):
    outfasta = options.output + '/CRISPR.fasta'
    spacers = ''
    if os.path.isfile(outfasta):
        fasta = open(outfasta, 'a')
    else:
        fasta = open(outfasta, 'w')
    for i in range (0, len(positions)-2):
        start = int(positions[i][1])+2 # +1 por el indice, +1 porque el spacer empieza despues del repeat
        end = int(positions[i+1][0])
        spacer = str(record.seq)[start:end]
        spacers = spacers + '>' + record.id + '_' + str(start) + '_' + str(end)
        spacers = spacers + '\n' + spacer + '\n'
    fasta.write(spacers)

## Look for CAS
def look_for_cas (record, options):
    run_prodigal (record, options)
#	run_hmmsearch
#	find_best_cas

def run_prodigal (record, options):
    if os.path.isfile('contig.fasta'):
        print("esta")
    else:
        print("no esta")

#def run_hmmsearch ():

#def find_best_cas ():

## Wrap everything

## Main
def main():
    #bool_dependencies = check_dependencies()
    #if not bool_dependencies:
        #install_dependencies()
    parser = argparse.ArgumentParser(description='')
    options = parse_option(parser)
    current_dir = os.getcwd()
    if not re.match('^/', options.output):
        options.output = current_dir + '/' + options.output
    if os.path.exists(options.output):
        print("ERROR: la carpeta de salida existe")
        exit
    os.makedirs(options.output)
    if not os.path.exists("/tmp/casper_search"):
        os.makedirs ("/tmp/casper_search")
    with open(options.input, "rU") as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            look_for_spacers(record, options)
            look_for_cas(record, options)
    shutil.rmtree("/tmp/casper_search")

if __name__ == "__main__":
    main()
