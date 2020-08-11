#!/usr/bin/env python3

__author__ = 'Matias Irazoqui'

import os
import subprocess
import argparse
import re
import shutil
import time
import pprint
from Bio import SeqIO

## Global variables ##
current_dir = os.getcwd()
script_dir = os.path.dirname(os.path.abspath(__file__))
so_path = script_dir + '/data/sel392v2.so'
profile_path = script_dir + '/data/originales.hmm'
temp_dir = '/tmp/caspr_search'

## Check dependencies ##
#def check_dependencies()
#def install_dependencies()

## Parse options ##
def parse_option(parser):
    parser.add_argument('-i', '--input', type=str, required=True, help='Archivo FASTA a analizar')
    parser.add_argument('-o', '--output', type=str, required=True, help='Directorio de salida')
    parser.add_argument('-s', '--so', type=str, default=so_path, help='Archivo so requerido por vmatch')
    parser.add_argument('--minspacer', type=int, default=17, help='Tamano minimo del spacer')
    parser.add_argument('--maxspacer', type=int, default=50, help='Tamano maximo del spacer')
    parser.add_argument('-m', '--mismatch', type=int, default=0, help='Numero de mismatches permitidos en direct repeat')
    parser.add_argument('--meta', action='store_true')
    parser.add_argument('--profiles', default=profile_path, help='Archivo con los perfiles a buscar')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    args.profiles = os.path.abspath(args.profiles)
    return args

## Look for spacers ##
def look_for_spacers (record, options):
    dr = run_vmatch2 (record, options)
    extract_spacers (record, dr['positions'], options)
    return dr

def run_vmatch2 (record, options):
    work_dir = temp_dir + '/' + record.id
    verbose = options.verbose
    if not os.path.exists(work_dir):
        os.makedirs (work_dir)
    os.chdir(work_dir)

    SeqIO.write(record, 'contig.fasta', 'fasta')

    cmd = 'mkvtree -dna -pl -lcp -suf -tis -ois -bwt -bck -sti1 -db contig.fasta'
    subprocess.run(cmd.split())

    cmd = 'vmatch -l 25 %s %s -s leftseq -evalue 0.001 -absolute -noevalue -noscore \
          -noidentity -sort sd -best 10000 -selfun %s \
          contig.fasta' % (options.minspacer, options.maxspacer, options.so)
    salida = subprocess.getoutput(cmd)

    best = { 'id': record.id,
             'score': 0,
             'mismatchs': 0,
             'positions': [],
             'dr': '' }
    cont = 0
    for line in salida.split('\n'):
        line = line.rstrip()
        if re.match('^>', line) or not line:
            continue
        (score, positions, mismatchs) = run_fuzznuc(line, options, cont)
        cont += 1
        if score > best['score']:
            best['score'] = score
            best['positions'] = positions
            best['dr'] = line
            best['mismatchs'] = mismatchs
    return best

def run_fuzznuc (pattern, options, cont):
    score = 0
    mismatchs = 0
    positions = []
    cmd = 'fuzznuc -sequence contig.fasta -pattern %s -pmismatch %s -outfile %s' \
    % (pattern, options.mismatch, str('fuzznuc'+str(cont)+'.txt'))
    subprocess.run(cmd.split(), stderr = subprocess.PIPE)
    with open(str('fuzznuc'+str(cont)+'.txt'), 'r') as fuzz_report:
        for line in fuzz_report:
            if re.match('^\s+\d', line):
                columns = line.split()
                columns[4] = 0 if columns[4] == '.' else columns[4]
                score = score + len(columns[5]) - int(columns[4])
                if int(columns[4]) > 0: mismatchs += 1
                positions.append((int(columns[0]),int(columns[1])))
    return (score, positions, mismatchs)

def extract_spacers (record, positions, options):
    outfasta = options.output + '/CRISPR.fasta'
    spacers = ''
    if os.path.isfile(outfasta):
        fasta = open(outfasta, 'a')
    else:
        fasta = open(outfasta, 'w')
    for i in range (0, len(positions)-2):
        start = positions[i][1]+2 # +1 por el indice, +1 porque el spacer empieza despues del repeat
        end = positions[i+1][0]
        spacer = str(record.seq)[start:end]
        spacers = spacers + '>' + record.id + '_' + str(start) + '_' + str(end)
        spacers = spacers + '\n' + spacer + '\n'
    fasta.write(spacers)

## Look for CAS ##
def look_for_cas (record, options):
    run_prodigal (options)
    run_hmmsearch (options)
    genes = find_best_cas ()
    return genes

def run_prodigal (options):
    cmd = 'prodigal -a translation.faa -i contig.fasta -o translation.gbk'
    if options.meta: cmd += ' -p meta'
    subprocess.run(cmd.split(), stderr=subprocess.PIPE)

def run_hmmsearch (options):
    cmd = 'hmmsearch --domtblout search.txt %s translation.faa' % (options.profiles)
    subprocess.run(cmd.split(), stderr=subprocess.PIPE, stdout=subprocess.PIPE)

def find_best_cas ():
    genes = {}
    strand = 0
    if os.path.exists('search.txt'):
        with open('search.txt', 'r') as search:
            for line in search:
                if not re.match('^#', line):
                    fields = line.split()
                    if not fields[0] in genes:
                        partial = re.match('.*partial\=(\d+)', fields[29])
                        strand=int(fields[27])
                        genes[fields[0]] = {'len': int(fields[2]),
                                            'start': int(fields[23]),
                                            'end': int(fields[25]),
                                            'partial': partial.group(1),
                                            'annotation': '',
                                            'evalue': 9999,
                                            'aligned': 0,
                                            'score':999,
                                            'annotations': {} }
                    if not fields[3] in genes[fields[0]]['annotations']:
                        genes[fields[0]]['annotations'][fields[3]] = {'evalue': float(fields[6]),
                                                                 'aligned': ((int(fields[18]) - int(fields[17]) + 1)/int(fields[2])) }
                    else:
                        genes[fields[0]]['annotations'][fields[3]]['evalue'] += float(fields[6])
                        genes[fields[0]]['annotations'][fields[3]]['aligned'] += float((int(fields[18]) - int(fields[17]) + 1)/int(fields[2]))
    else:
        genes['no_genes'] = {'len': 0,
                             'start': 0,
                             'end': 1,
                             'partial': 00,
                             'annotation': 'No genes predicted',
                             'evalue': 9999,
                             'aligned': 0,
                             'score':999,
                             'annotations': {} }
        return (genes, strand)
    for g in genes:
        for a in genes[g]['annotations']:
            score = genes[g]['annotations'][a]['evalue']/genes[g]['annotations'][a]['aligned']
            if score < genes[fields[0]]['score'] and genes[g]['annotations'][a]['aligned'] >= 0.7:
                genes[g]['annotation'] = a
                genes[g]['evalue'] = genes[g]['annotations'][a]['evalue']
                genes[g]['aligned'] = genes[g]['annotations'][a]['aligned']
                genes[g]['score'] = score
    return (genes, strand)

## Write output ##
def write_output(results, options):
    outfile = options.output + '/Results.tsv'
    report = 'Contig\t# repeats\tLength repeats\tStrand\tStart repeats\tEnd repeats\t# repeats w/mismatchs\
              \t#CRISPRs\tAvg length CRISPR\tStart Cas\tEnd Cas\tCas cassette\n'
    for c in results:
        start_dr = 0
        end_dr =0
        if results[c]['strand'] == 1:
            start_dr = results[c]['crisprs']['positions'][0][0]
            end_dr = results[c]['crisprs']['positions'][-1][1]
        elif results[c]['strand'] == -1:
            start_dr = results[c]['crisprs']['positions'][-1][1]
            end_dr = results[c]['crisprs']['positions'][0][0]

        number_repeats = len(results[c]['crisprs']['positions'])
        length = 0
        avg_length = 0
        if number_repeats > 0:
            for i in range(number_repeats-1, 1, -1):
                length += results[c]['crisprs']['positions'][i][0] - results[c]['crisprs']['positions'][i-1][1] + 1
            avg_length = length / (number_repeats - 1)

        cas_cassette = []
        strand = 0
        for cas in results[c]['genes']:
            cas_cassette.append((results[c]['genes'][cas]['annotation'], results[c]['genes'][cas]['start'], results[c]['genes'][cas]['end'], results[c]['genes'][cas]['partial']))
        if results[c]['strand']==1:
            cas_cassette.sort(key=lambda x:x[1])
        else:
            cas_cassette.sort(key=lambda x:x[1], reverse=True)

        report += c + '\t' + str(number_repeats) + '\t' + str(len(results[c]['crisprs']['dr'])) + '\t' + str(results[c]['strand']) + '\t'
        report += str(start_dr) + '\t' + str(end_dr) + '\t' + str(results[c]['crisprs']['mismatchs']) + '\t'
        report += str(number_repeats - 1) + '\t' + str(round(avg_length,2))
        if cas_cassette:
            report += '\t' + str(cas_cassette[0][1]) + '\t' + str(cas_cassette[-1][2]) + '\t'
            for cas in cas_cassette:
               report += cas[0]
               if cas[3] != '00': report += '*' # If partial, add '*' to the name of the gene
               report += ';'
            report.strip(';')
        report += '\n'

    with open (outfile, 'w') as o:
        o.write(report)

## Main ##
def main():
    #bool_dependencies = check_dependencies()
    #if not bool_dependencies:
        #install_dependencies()
    parser = argparse.ArgumentParser(description='')
    options = parse_option(parser)

    if not re.match('^/', options.output):
        options.output = current_dir + '/' + options.output
    if os.path.exists(options.output) and os.listdir(options.output):
        print('ERROR: la carpeta de salida no esta vacia')
        exit()
    if not os.path.exists(options.output):
        os.makedirs(options.output)
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    os.makedirs (temp_dir)

    results = {}
    with open(options.input, 'r') as in_fasta:
        for record in SeqIO.parse(in_fasta, 'fasta'):
            now = time.process_time()
            msg = 'Sequence ' + record.id + ' started at ' + str(time.asctime())
            if options.verbose: print (msg)
            
            crisprs = look_for_spacers(record, options)
            time_stamp = time.process_time() - now
            msg = 'look_for_spacers took ' + str(round(time_stamp,4)) + ' seconds'
            if options.verbose: print (msg)
            
            now = time.process_time()
            (genes,strand) = look_for_cas(record, options)
            time_stamp = time.process_time() - now
            msg = 'look_for_cas took ' + str(round(time_stamp, 4)) + ' seconds'
            if options.verbose: print (msg)
            
            if not 'no_genes' in genes and len(crisprs['positions']) != 0:
                results[record.id] = {'crisprs': crisprs,
                                      'genes': genes,
                                      'strand': strand }
    write_output(results, options)
    #write_cas_fasta()
    #shutil.rmtree(temp_dir)

if __name__ == '__main__':
    main()
