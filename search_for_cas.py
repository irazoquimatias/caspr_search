#!/usr/bin/env python3

__author__ = 'Matias Irazoqui'

import os
import subprocess
import argparse
import re
import shutil
import time
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool

## Global variables ##
current_dir = os.getcwd()
script_dir = os.path.dirname(os.path.abspath(__file__))
so_path = script_dir + '/data/sel392v2.so'
profile_path = script_dir + '/data/profiles.hmm'
temp_dir = '/tmp/caspr_search'
list_effectors = ['cas12d', 'cas12e', 'cas12i', 'cas13d', 'cas14a',
                 'cas14b', 'cas14c', 'cas14d', 'cas14e', 'cas14f',
                 'cas14g', 'cas14h', 'cas14u', 'cas13b_0', 'cas12c',
                 'cas12a_1']

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
    parser.add_argument('-t', '--threads', type=int, default=1, help='Numero de procesadores')
    parser.add_argument('-k', '--keeptemp', action='store_true')
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()
    args.profiles = os.path.abspath(args.profiles)
    return args

## Search everything
def search_everything (record):
    this_temp_dir = temp_dir + '/' + record.id
    dr = look_for_spacers(record)
    genes = look_for_cas(record)
    strand = 0
    for g in genes:
        if genes[g]['effector']:
            strand = genes[g]['strand']
            genes[g]['sequence'] = extract_protein (genes[g]['id'], this_temp_dir)
    if not options.keeptemp:
        shutil.rmtree(this_temp_dir)
    return (record.id, dr, genes, strand)

def extract_protein (gene_id, temp_dir):
    translation = temp_dir + '/translation.faa'
    with open(translation, "rU") as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if record.id == gene_id:
                return str(record.seq)

## Look for spacers ##
def look_for_spacers (record):
    dr = run_vmatch2 (record)
    extract_spacers (record, dr['positions'])
    return dr

def run_vmatch2 (record):
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
        (score, positions, mismatchs) = run_fuzznuc(line, cont)
        cont += 1
        if score > best['score']:
            best['score'] = score
            best['positions'] = positions
            best['dr'] = line
            best['mismatchs'] = mismatchs
    return best

def run_fuzznuc (pattern, cont):
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

def extract_spacers (record, positions):
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
def look_for_cas (record):
    run_prodigal ()
    run_hmmsearch ()
    genes = find_best_cas ()
    return genes

def run_prodigal ():
    cmd = 'prodigal -a translation.faa -i contig.fasta -o translation.gbk'
    if options.meta: cmd += ' -p meta'
    subprocess.run(cmd.split(), stderr=subprocess.PIPE)

def run_hmmsearch ():
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
                        genes[fields[0]] = {'id': '',
                                            'length': int(fields[2]),
                                            'start': int(fields[23]),
                                            'end': int(fields[25]),
                                            'strand': int(fields[27]),
                                            'partial': partial.group(1),
                                            'annotation': '',
                                            'evalue': 9999,
                                            'aligned': 0,
                                            'score':999,
                                            'effector':'FALSE',
                                            'annotations': {} }
                    if not fields[3] in genes[fields[0]]['annotations']:
                        genes[fields[0]]['annotations'][fields[3]] = { 'id': fields[0],
                                                                 'length': int(fields[2]),
                                                                 'evalue': float(fields[6]),
                                                                 'aligned': ((int(fields[18]) - int(fields[17]) + 1)/int(fields[2])),
                                                                 'hmm_length': int(fields[5]),
                                                                 'hmm_start': int(fields[15]),
                                                                 'hmm_end': int(fields[16])}
                    else:
                        genes[fields[0]]['annotations'][fields[3]]['evalue'] += float(fields[6])
                        genes[fields[0]]['annotations'][fields[3]]['aligned'] += float((int(fields[18]) - int(fields[17]) + 1)/int(fields[2]))
                        if int(fields[15]) < genes[fields[0]]['annotations'][fields[3]]['hmm_start']:
                            genes[fields[0]]['annotations'][fields[3]]['hmm_start'] = int(fields[15])
                        if int(fields[16]) > genes[fields[0]]['annotations'][fields[3]]['hmm_end']:
                            genes[fields[0]]['annotations'][fields[3]]['hmm_end'] = int(fields[16])
    else:
        genes['no_genes'] = {'id': 'None_found',
                             'length': 0,
                             'strand': 0,
                             'start': 0,
                             'end': 1,
                             'partial': 00,
                             'annotation': 'No genes predicted',
                             'evalue': 9999,
                             'aligned': 0,
                             'score':999,
                             'effector':'FALSE',
                             'annotations': {} }
        return (genes)
    for g in genes:
        for a in genes[g]['annotations']:
            score = genes[g]['annotations'][a]['evalue']/genes[g]['annotations'][a]['aligned']
            if score < genes[fields[0]]['score']:
                if a in list_effectors:
                    genes[g]['effector'] = 'TRUE'
                genes[g]['annotation'] = a
                genes[g]['id'] = genes[g]['annotations'][a]['id']
                genes[g]['length'] = genes[g]['annotations'][a]['length']
                genes[g]['evalue'] = genes[g]['annotations'][a]['evalue']
                genes[g]['aligned'] = genes[g]['annotations'][a]['aligned']
                genes[g]['hmm_length'] = genes[g]['annotations'][a]['hmm_length']
                genes[g]['hmm_start'] = genes[g]['annotations'][a]['hmm_start']
                genes[g]['hmm_end'] = genes[g]['annotations'][a]['hmm_end']
                genes[g]['score'] = score
    return (genes)

## Write output ##
def write_output(results):
    outfile = options.output + '/Results.tsv'
    outfasta = options.output + '/Effector.fasta'
    fasta_effectors = ''
    report = 'Contig\t# repeats\tLength repeats\tStrand\tStart repeats\tEnd repeats\t'
    report += 'Sequence repeat\t# repeats w/mismatchs\t#CRISPRs\tAvg length CRISPR\t'
    report += 'Start Cas\tEnd Cas\tCas cassette\tEfector\tLargo efector\t'
    report += 'Length HMM efector\tStart HMM efector\tEnd HMM efector\t'
    report += 'Sequence efector\n'
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
        ids = []
        effector = []
        strand = 0
        for cas in results[c]['genes']:
            if results[c]['genes'][cas]['id'] not in ids:
                cas_cassette .append((results[c]['genes'][cas]['annotation'], results[c]['genes'][cas]['start'],
                                      results[c]['genes'][cas]['end'], results[c]['genes'][cas]['partial']))
                ids.append(results[c]['genes'][cas]['id'])
            if results[c]['genes'][cas]['effector'] == 'TRUE':
                effector = [results[c]['genes'][cas]['annotation'], results[c]['genes'][cas]['length'],
                            results[c]['genes'][cas]['hmm_length'], results[c]['genes'][cas]['hmm_start'],
                            results[c]['genes'][cas]['hmm_end'], results[c]['genes'][cas]['sequence']]
        if results[c]['strand']==1:
            cas_cassette.sort(key=lambda x:x[1])
        else:
            cas_cassette.sort(key=lambda x:x[1], reverse=True)

        report += c + '\t' + str(number_repeats) + '\t' + str(len(results[c]['crisprs']['dr'])) + '\t' + str(results[c]['strand']) + '\t'
        report += str(start_dr) + '\t' + str(end_dr) + '\t' + results[c]['crisprs']['dr'] + '\t'
        report += str(results[c]['crisprs']['mismatchs']) + '\t'
        report += str(number_repeats - 1) + '\t' + str(round(avg_length,2)) + '\t'

        if cas_cassette:
            report += str(cas_cassette[0][1]) + '\t' + str(cas_cassette[-1][2]) + '\t'
            for cas in cas_cassette:
               report += cas[0]
               if cas[3] != '00': report += '*' # If partial, add '*' to the name of the gene
               report += ';'
            report.strip(';')

        if effector:
            report += '\t' + effector[0] + '\t' + str(effector[1]) + '\t' + str(effector[2]) + '\t'
            report += str(effector[3]) + '\t' + str(effector[4]) + '\t' + str(effector[5]) + '\t'
            fasta_effectors += '>' + c + '_' + effector[0] + '\n' + str(effector[5]) + '\n'
        else:
            report += 'No effector found\t0\t\t\t\t\t\t'
        report += '\n'

    with open (outfile, 'w') as o:
        o.write(report)
    with open (outfasta, 'w') as of:
        of.write(fasta_effectors)

## Main ##
def main():
    #bool_dependencies = check_dependencies()
    #if not bool_dependencies:
        #install_dependencies()
    parser = argparse.ArgumentParser(description='')
    global options
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

    now = time.process_time()
    msg = 'Process started at ' + str(time.asctime())
    if options.verbose: print (msg)

    with open(options.input) as fin, Pool(int(options.threads)) as pool:
        # [0]: record.id; [1] repeats; [2]: genes, [3]: strand
        contigs = pool.map(
            search_everything,
            (record for record in SeqIO.parse(fin, 'fasta')),
            chunksize=50,
            )

    time_stamp = time.process_time() - now
    msg = 'Process ended at '  + str(time.asctime()) + '\n'
    msg += 'Process took ' + str(time_stamp) + ' seconds'
    if options.verbose: print (msg)

    results = {}
    for c in contigs:
        if not 'no_genes' in c[2] and len(c[1]['positions']) != 0:
            results[c[0]] = {'crisprs': c[1], 'genes': c[2], 'strand': c[3] }
    write_output(results)
    if not options.keeptemp:
        shutil.rmtree(temp_dir)

if __name__ == '__main__':
    main()
