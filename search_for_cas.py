#!/usr/bin/env python3

__author__ = "Matias Irazoqui"

import os
import subprocess
import argparse
import shutil
import time
from Bio import SeqIO
import re

## Global variables ##
script_dir = os.path.dirname(os.path.abspath(__file__))
hmmer_dir = script_dir + '/data/profiles/'
temp_dir = '/tmp/casper_search'

## Check dependencies ##
#def check_dependencies()
#def install_dependencies()

## Parse options ##
def parse_option(parser):
    parser.add_argument('-i', '--input', type=str, required=True, help='Archivo FASTA a analizar')
    parser.add_argument('-o', '--output', type=str, required=True, help='Directorio de salida')
    parser.add_argument('-s', '--so', type=str, required=True, help='Archivo so requerido por vmatch')
    parser.add_argument('--minspacer', type=int, default=17, help='Tamano minimo del spacer')
    parser.add_argument('--maxspacer', type=int, default=50, help='Tamano maximo del spacer')
    parser.add_argument('-m', '--mismatch', type=int, default=0, help='Numero de mismatches permitidos en direct repeat')
    parser.add_argument('--meta', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    return args

## Look for spacers ##
def look_for_spacers (record, options):
    (positions, dr) = run_vmatch2 (record, options)
    extract_spacers (record, positions, options)

def run_vmatch2 (record, options):
    work_dir = temp_dir #+ record.id
    verbose = options.verbose
    if not os.path.exists(work_dir):
        os.makedirs (work_dir)
    os.chdir(work_dir)

    SeqIO.write(record, 'contig.fasta', 'fasta')

    cmd = 'mkvtree2 -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1 -db contig.fasta'
    subprocess.run(cmd.split())

    cmd = 'vmatch2 -l %s %s -s leftseq -evalue 0.001 -absolute -noevalue -noscore \
          -noidentity -sort sd -best 10000 -selfun %s \
          contig.fasta' % (options.minspacer, options.maxspacer, options.so)
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
    subprocess.run(cmd.split(), stderr = subprocess.PIPE)
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

## Look for CAS ##
def look_for_cas (record, options):
    run_prodigal (options)
    run_hmmsearch (options)
#	find_best_cas

def run_prodigal (options):
    cmd = "prodigal -a translation.faa -i contig.fasta -o translation.gbk"
    if options.meta: cmd += " -p meta"
    subprocess.run(cmd.split(), stderr=subprocess.PIPE)

def run_hmmsearch (options):
    cmd = "hmmsearch --domtblout search.txt $i cas14a_proteins.faa
#def find_best_cas ():

## Wrap everything ##

## Main ##
def main():
    #bool_dependencies = check_dependencies()
    #if not bool_dependencies:
        #install_dependencies()
    parser = argparse.ArgumentParser(description='')
    options = parse_option(parser)
    current_dir = os.getcwd()
    if not re.match('^/', options.output):
        options.output = current_dir + '/' + options.output
    if os.path.exists(options.output) and os.listdir(options.output):
        print("ERROR: la carpeta de salida no esta vacia")
        exit()
    if not os.path.exists(options.output):
        os.makedirs(options.output)
    if not os.path.exists(temp_dir):
        os.makedirs (temp_dir)
    with open(options.input, "rU") as in_fasta:
        for record in SeqIO.parse(in_fasta, "fasta"):
            now = time.process_time()
            msg = "Sequence " + record.id + " started at " + str(time.asctime())
            print (msg)
            look_for_spacers(record, options)
            time_stamp = time.process_time() - now
            msg = "look_for_spacers took " + str(time_stamp)
            print (msg)
            look_for_cas(record, options)
            time_stamp = time.process_time() - now
            msg = "look_for_cas took " + str(time_stamp)
            print (msg)
     #shutil.rmtree(temp_dir)

if __name__ == "__main__":
    main()
