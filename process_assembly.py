#!/usr/bin/env python3

__author__ = "Ezequiel Sosa - Matias Irazoqui"

from tqdm import tqdm
import argparse
import os
import subprocess as sp
import Bio.SeqIO as bpio
from termcolor import colored
from glob import glob
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from glob import glob
from tqdm import tqdm
import collections
import re


## Functions
def exec(cmd, verbose=False):
    if verbose:
        print(cmd)
    sp.call(cmd, shell=True)


def sample_file(item):
    fasta_path = f'/tmp/{item.id.replace("/", "_").replace(" ", "_").replace("|", "_")}.fasta'
    bpio.write(item, fasta_path, "fasta")
    return fasta_path


def get_codon(pos, genes_dict, gene):
    pos_in_genome = pos + int(genes_dict[gene]["gapped_coding_location"][0][0])
    cant_aa = 0
    for cr in genes_dict[gene]["gapped_coding_location"]:
        cant_aa += int((int(cr[0]) - int(genes_dict[gene]["gapped_coding_location"][0][0]) + 1) / 3)
        if int(cr[0]) <= pos_in_genome and pos_in_genome < int(cr[1]):
            pos_aa = cant_aa + (int((pos_in_genome - int(cr[0])) / 3)) + 1
            rest = (pos_in_genome - int(cr[0]) + 1) / 3 - int((pos_in_genome - int(cr[0]) + 1) / 3)
            if rest == 0:
                return (pos - 2, pos + 1, pos_aa)
            elif rest > 0.5:
                return (pos - 1, pos + 2, pos_aa)
            else:
                return (pos, pos + 3, pos_aa)

def fix_msa (genes, gene):
    i = 0
    for loc in genes[gene]["gapped_coding_location"]:
        while i < (int(loc[1]) - int(loc[0])):
            codon = genes[gene]["dna_seq"][i:i+3]
            if '-' in codon:
                extra_seq = ''
                start = i
                i += 3
                next_codon = genes[gene]["dna_seq"][i:i+3]
                while next_codon == '---':
                    extra_seq += next_codon
                    i += 3
                    next_codon = genes[gene]["dna_seq"][i:i+3]
                if re.search( '\w\w-', codon) and re.search( '--\w', next_codon):
                    new_codon = codon[0:2] + next_codon[-1]
                    genes[gene]["dna_seq"] = genes[gene]["dna_seq"][:start] + new_codon + extra_seq + "---" + genes[gene]["dna_seq"][i+3:]
                elif re.search( '\w--', codon) and re.search( '-\w\w', next_codon):
                    new_codon = codon[0] + next_codon[1:3]
                    genes[gene]["dna_seq"] = genes[gene]["dna_seq"][:start] + "---" + extra_seq + new_codon + genes[gene]["dna_seq"][i+3:]
            i += 3


###main

parser = argparse.ArgumentParser(description='')
parser.add_argument("-i", '--fasta', type=str, help='directorio de fastas no alineados')
parser.add_argument("-o", '--dir_salida', type=str, help='Directorio donde se ponen los reportes de salida')
parser.add_argument("-p", '--primers', type=str, help='fasta con primers')
parser.add_argument("-v", '--verbose', action='store_true')
parser.add_argument("-r", '--reference', help="fasta reference", default="/ref/MN996528.fna")
parser.add_argument("-a", '--annotation', help="fasta reference", default="/app/snpEff/data/covid19/genes.gbk")
args = parser.parse_args()

"""
1) % de cada ORF con datos (entiéndase: bases nucleotídicas)
2) % NNNs en cada ORF
3) gaps encontrados en los ORF
4) Codones START atrasados
5) Codones de STOP prematuros
6) Detección de mutaciones en los primers que se usan para hacer las qPCR de diagnóstico
"""
din = args.fasta
out = args.dir_salida
primers = args.primers
verbose = args.verbose

assert os.path.exists(din), f"no existe el directorio de entrada {din}"
if not os.path.exists(out):
    os.makedirs(out)

assert os.path.exists(args.reference), f"no existe la referencia {args.reference}"
assert os.path.exists(args.annotation), f"no existe la anotacion {args.annotation}"

if os.path.isdir(din):
    pbar = tqdm(glob(f'{din}/*.fasta'))
    sample_fn = lambda item: os.path.basename(item).split(".fasta")[0]
    sample_file = lambda item: item
else:
    pbar = tqdm(list(bpio.parse(din, "fasta")))
    sample_fn = lambda item: item.id.replace("/", "_").replace(" ", "_").replace("|", "_")

general_report = "Secuencia\tReferencia usada\tLargo (nt)\t# Ns\t%Ns\tGenes completos\tGenes incompletos\tGenes no encontrados\t# Mut. (nt)\t% Mut.\t# Mut. Sin.\t# Mut. No Sin.\n"
genes_report = "Gen\tMuestra\tLargo\tInicio\tFin\tBases_identicas\tMismatches\tInserciones\tDeleciones\tNs\t%Ns\tMutaciones_sinonimas\tMutaciones_no_sinonimas\n"
record = bpio.read(args.annotation, "genbank")
ref_id = record.id
sample_to_warn = {}

with open(f"{out}/primers.txt", "w") as hw:
    hw.write("#sample\tprimer\tpos\tpos_primer\tvariant\n")

for item in pbar:
    fasta_file = sample_file(item)
    sample = sample_fn(item)
    pbar.set_description(fasta_file)
    msa = "/tmp/aln.fasta"
    genomes = {}
    ref_genes = {}
    sample_genes = {}

    if not os.path.exists(f"{out}/samples/{sample}"):
        os.makedirs(f"{out}/samples/{sample}")

    cmd = f'cat  {args.reference} {fasta_file} > /tmp/seqs.fasta'
    exec(cmd, verbose=verbose)
    cmd = f'mafft --auto /tmp/seqs.fasta > {msa} 2>>"{out}/samples/{sample}/log.txt"'
    exec(cmd, verbose=verbose)

    ## Read reference
    for feature in record.features:
        if feature.type == "CDS" or feature.type == "mat_peptide":
            locations = re.findall('\[(\d+)\:(\d+)\]', str(feature.location))
            ref_genes[feature.qualifiers['product'][0]] = {
                "length": (int(locations[-1][1]) - int(locations[0][0])),
                "coding_location": locations,  # positions of coding regions (in case there's splicing)
                "coding_regions": len(locations)  # number of coding regions
            }

    ##  Read MSA ##
    msa_dict = bpio.to_dict(bpio.parse(msa, "fasta"))
    assert ref_id in msa_dict, f"ref {ref_id} not found in msa ({','.join(msa_dict)})"

    for seq_record in msa_dict.values():
        list_gaps = []
        list_ns = []
        length = len(str(seq_record.seq))
        Ns = 0
        if seq_record.id != ref_id:
            sample_id = seq_record.id

        gap_iter = re.finditer('-+', str(seq_record.seq))
        for match in gap_iter:
            list_gaps.append((match.span()[0], match.span()[1]))
            length = length - match.span()[1] + match.span()[0]
        unamb_length = length
        n_iter = re.finditer('n+', str(seq_record.seq).lower())
        for match in n_iter:
            if match.span()[0] == 0 or match.span()[1] >= length:
                unamb_length = unamb_length - match.span()[1] + match.span()[0]
                continue 
            Ns += match.span()[1] - match.span()[0]
            if (match.span()[1] - match.span()[0] >= 20):
                list_ns.append((match.span()[0], match.span()[1]))

        genomes[seq_record.id] = {"Ns": Ns,
                                  "gaps": list_gaps,
                                  "n_islands": list_ns,
                                  "length": unamb_length,
                                  "seq": str(seq_record.seq).lower(),
                                  "mut": 0,
                                  "syn_mut": 0,
                                  "non_syn_mut": 0,
                                  "complete_genes": [],
                                  "incomplete_genes": [],
                                  "absent_genes": []
                                  }
        del seq_record

    ## Find coding regions in MSA ##
    for gene in ref_genes:
        if genomes[ref_id]["gaps"]:
            ref_genes[gene]["gapped_coding_location"] = []
            for m in ref_genes[gene]["coding_location"]:
                new_start = int(m[0])
                new_end = int(m[1])
                for gap_start, gap_end in (reversed(genomes[ref_id]['gaps'])):
                    if new_start >= gap_end:
                        new_start += gap_end - gap_start
                        new_end += gap_end - gap_start
                    elif new_start < gap_start and new_end >= gap_end:
                        new_end += gap_end - gap_start
                ref_genes[gene]["gapped_coding_location"].append((new_start, new_end))
        else:
            ref_genes[gene]["gapped_coding_location"] = ref_genes[gene]["coding_location"]

    ## Extract genes ##
    for gene in ref_genes:
        sample_genes[gene] = {}
        start = int(ref_genes[gene]["gapped_coding_location"][0][0])
        end = int(ref_genes[gene]["gapped_coding_location"][-1][1])
        for genome in genomes:
            if genome == ref_id:
                ref_genes[gene]["dna_seq"] = genomes[genome]["seq"][start:end]
                frequencies = collections.Counter(ref_genes[gene]["dna_seq"])
                sample_genes[gene]["insertions"] = frequencies['-']
            else:
                sample_genes[gene]["dna_seq"] = genomes[genome]["seq"][start:end]
                sample_genes[gene]["gapped_coding_location"] = ref_genes[gene]["gapped_coding_location"]
                frequencies = collections.Counter(sample_genes[gene]["dna_seq"])
                sample_genes[gene]["Ns"] = frequencies['n']
                sample_genes[gene]["deletions"] = frequencies['-']
                if sample_genes[gene]["Ns"] < 20:
                    genomes[sample_id]["complete_genes"].append(gene)
                elif sample_genes[gene]["Ns"] < ref_genes[gene]["length"]:
                    genomes[sample_id]["incomplete_genes"].append(gene)
                else:
                    genomes[sample_id]["absent_genes"].append(gene)
                sample_genes[gene]["ident"] = 0
                sample_genes[gene]["mut"] = 0
                sample_genes[gene]["non_syn_mut"] = {}
                sample_genes[gene]["syn_mut"] = {}
                fix_msa(ref_genes, gene)
                fix_msa(sample_genes, gene)

    ## Compare genes ##
    for gene in ref_genes:
        ref_genes[gene]["gaps"] = []
        for i in range(0, len(ref_genes[gene]["dna_seq"])):
            bool_island = 0
            if ref_genes[gene]["dna_seq"][i] != sample_genes[gene]["dna_seq"][i]:
                for island in genomes[sample_id]["n_islands"]:
                    if island[0] <= (i + int(ref_genes[gene]["coding_location"][0][0])) and (
                        i + int(ref_genes[gene]["coding_location"][0][0])) < island[1]:
                        bool_island = 1
                if bool_island == 1:
                    continue
                genomes[sample_id]["mut"] += 1
                sample_genes[gene]["mut"] += 1
                (start, end, aa_pos) = get_codon(i, ref_genes, gene)

                if "-" not in str(ref_genes[gene]["dna_seq"][start:end]):
                    ref_aa = str(Seq(str(ref_genes[gene]["dna_seq"][start:end])).translate(gap="-"))
                elif str(ref_genes[gene]["dna_seq"][start:end]) == "---":
                    while re.search('^-+$', str(ref_genes[gene]["dna_seq"][start:end])):
                        start -= 3
                        ref_aa = str(Seq(str(ref_genes[gene]["dna_seq"][start:end])).translate(gap="-"))
                        if not (start, start + 3) in ref_genes[gene]["gaps"]:
                            ref_genes[gene]["gaps"].append((start, start + 3))
                    ref_aa = ref_aa.replace("-", "")
                else:
                    ref_aa = "*"
                    ins_codon = ""
                    for a in range(start, end):
                        if str(ref_genes[gene]["dna_seq"][a]) != str(sample_genes[gene]["dna_seq"][a]):
                            ins_codon += str(sample_genes[gene]["dna_seq"][a]).upper()
                        else:
                            ins_codon += str(sample_genes[gene]["dna_seq"][a])
                    warn = sample_id+"\t"+gene+"\t"+str(start)+"-"+str(end-1)+"\t"+ins_codon+"\tinsercion"
                    sample_to_warn[warn] = {"bool" : 1}
                    
                if "-" not in str(sample_genes[gene]["dna_seq"][start:end]):
                    sample_aa = Seq(str(sample_genes[gene]["dna_seq"][start:end])).translate()
                    if sample_aa == "*":
                        sample_aa = "STOP"
                elif str(sample_genes[gene]["dna_seq"][start:end]) == "---" :
                    sample_aa = "-"
                else:
                    sample_aa = "*"
                    warn = sample_id+"\t"+gene+"\t"+str(start)+"-"+str(end-1)+"\t"+str(sample_genes[gene]["dna_seq"][start:end])+"\tdelecion"
                    sample_to_warn[warn] = {"bool" : 1}
                if ref_aa != sample_aa:
                    mut = ref_aa[0] + str(aa_pos) + sample_aa[0]
                    genomes[sample_id]["non_syn_mut"] += 1
                    sample_genes[gene]["non_syn_mut"][aa_pos - len(ref_genes[gene]["gaps"])] = {"ref": str(ref_aa),
                                                                                                    "alt": str(sample_aa)}
                else:
                    sample_genes[gene]["syn_mut"][int(i + 1 - (len(ref_genes[gene]["gaps"]) * 3))] = {
                        "ref": ref_genes[gene]["dna_seq"][i], "alt": sample_genes[gene]["dna_seq"][i]}
                    mut = (ref_genes[gene]["dna_seq"][i] + str(i + 1) + sample_genes[gene]["dna_seq"][i]).upper()
                    genomes[sample_id]["syn_mut"] += 1
            else:
                sample_genes[gene]["ident"] += 1

    ## Printing outputs ##
    general_report += sample_id + "\t" + ref_id + "\t"
    general_report += str(genomes[sample_id]["length"]) + "\t"
    general_report += str(genomes[sample_id]["Ns"]) + "\t"
    general_report += str(round(genomes[sample_id]["Ns"] / int(genomes[sample_id]["length"])*100, 2)) + "\t"
    general_report += str(len(genomes[sample_id]["complete_genes"])) + "\t"
    general_report += str(len(genomes[sample_id]["incomplete_genes"])) + "\t"
    general_report += str(len(genomes[sample_id]["absent_genes"])) + "\t"
    general_report += str(genomes[sample_id]["mut"]) + "\t"
    general_report += str(round(genomes[sample_id]["mut"] / genomes[sample_id]["length"]*100, 2)) + "\t"
    general_report += str(genomes[sample_id]["syn_mut"]) + "\t"
    general_report += str(genomes[sample_id]["non_syn_mut"]) + "\n"

    for gene in ref_genes:
        genes_report += gene + "\t" + sample_id + "\t" + str(ref_genes[gene]["length"]) + "\t"
        genes_report += ref_genes[gene]["coding_location"][0][0] + "\t" + ref_genes[gene]["coding_location"][-1][
            1] + "\t"
        genes_report += str(sample_genes[gene]["ident"]) + "\t"
        genes_report += str(sample_genes[gene]["mut"]) + "\t"
        genes_report += str(sample_genes[gene]["insertions"]) + "\t"
        genes_report += str(sample_genes[gene]["deletions"]) + "\t"
        genes_report += str(sample_genes[gene]["Ns"]) + "\t"
        genes_report += str(round(sample_genes[gene]["Ns"]/ref_genes[gene]["length"]*100, 2)) + "\t"
        for m in sample_genes[gene]["syn_mut"]:
            genes_report += str(sample_genes[gene]["syn_mut"][m]["ref"]) + str(m) + str(
                sample_genes[gene]["syn_mut"][m]["alt"]) + ";"
        genes_report = genes_report.rstrip(';')
        genes_report += "\t"
        for m in sample_genes[gene]["non_syn_mut"]:
            genes_report += str(sample_genes[gene]["non_syn_mut"][m]["ref"]) + str(m) + str(
                sample_genes[gene]["non_syn_mut"][m]["alt"]) + ";"
        genes_report = genes_report.rstrip(';')
        genes_report += "\n"

    ## Primers
    if primers:
        assert os.path.exists(primers), f"no existe el archivo de primers{primers}"
        cmd = f'python3 /app/script/MSAMap.py -r MN996528.1  -i {msa} >> "{out}/samples/{sample}/{sample}_variants.vcf" 2>>"{out}/samples/{sample}/log.txt"'
        exec(cmd, verbose=verbose)

        cmd = f'java -jar /app/snpEff/snpEff.jar -stats {out}/samples/{sample}/snpEff.html covid19  "{out}/samples/{sample}/{sample}_variants.vcf" > "{out}/samples/{sample}/{sample}_ann.vcf" 2>>"{out}/samples/{sample}/log.txt"'
        exec(cmd, verbose=verbose)

        with open(f'{out}/samples/{sample}/{sample}_ann.vcf') as h, open(f'{out}/samples/{sample}/{sample}_vars.txt',
                                                                         "w") as hw:
            hw.write("pos ref alt annotation gene hgvs_c hgvs_p\n")
            for line in h:
                if not line.startswith("#"):
                    contig, pos, _, ref, alt, _, _, info = line.split("\t")[:8]
                    ann_str = info.split("ANN=")[1].split(",")[0]

                    (alt, annotation, impact, gene, geneid, feature_type, feature_id, transcript_biotype,
                     rank_div_total, hgvs_c, hgvs_p, c_dna_pos, cds_pos, (aa_pos_aa_len), dist_to_feature,
                     errors) = ann_str.split("|")
                    hw.write(" ".join([pos, ref, alt, annotation, gene, hgvs_c, hgvs_p]) + "\n")

        cmd = f'blastn -task "blastn-short"  -query {primers} -db {args.reference} -qcov_hsp_perc 100 -outfmt "6 sseqid sstart send qseqid" > /tmp/{sample}primers_raw.bed 2>/dev/null'
        exec(cmd, verbose=verbose)
        with open(f"{out}/samples/{sample}/primers.bed", "w") as h:
            for line in open(f'/tmp/{sample}primers_raw.bed'):
                vec = line.split("\t")
                if int(vec[1]) > int(vec[2]):
                    vec[1], vec[2] = vec[2], vec[1]
                h.write("\t".join(vec))
        cmd = f'bedtools  intersect -a "{out}/samples/{sample}/{sample}_ann.vcf" -b "{out}/samples/{sample}/primers.bed" -wb > "/tmp/{sample}primers.intersect" 2>>"{out}/samples/{sample}/log.txt"'
        exec(cmd, verbose=verbose)
        with open(f'/tmp/{sample}primers.intersect') as h, open(f'{out}/samples/{sample}/{sample}_primers.txt',
                                                                "w") as hw:
            hw.write("#sample\tprimer\tpos_ref\tpos_primer\tvariant\n")
            for line in h:
                vec = line.strip().split()
                # MN996528.1      3037    .       C       T       225     .       DP=26;VDB=0.343005;SGB=-0.690438;MQSB=1;MQ0F=0;AC=1;AN=1;DP4=0,0,5,12;MQ=60     GT:PL   1:255,0 MN996528.1      3028    3047    p2
                pos, ref, alt, primer, pstart, pend = int(vec[1]), vec[3], vec[4], vec[-1], int(vec[-3]), int(vec[-2])
                vcf_end = pos + len(ref) - 1
                start, end = max(pos, pstart), min(vcf_end, pend)
                pos_primer = start - pstart + 1
                if len(ref) > len(alt):
                    alt = f'NotFound -> Deletion at {pos} of {len(ref) - len(alt)}bps'
                    pos_primer = "-"
                elif len(ref) < len(alt):
                    alt = f'Insertion at {pos} of {len(alt) - len(ref)}bps'
                else:
                    if len(ref) == 1:
                        alt = f'{ref}->{alt}'
                    else:
                        if "N" in alt:
                            alt = "Presence of Ns in the region"
                        else:
                            alt = f'{ref[max(pos, start) - pos:end - pos + 1 ]}->{alt[max(pos, pstart) - pstart:end - pstart + 1 ]}'
                hw.write("\t".join([sample, primer, str(start), str(pos_primer), alt]) + "\n")

        cmd = f'grep -v "^#" "{out}/samples/{sample}/{sample}_primers.txt" >> {out}/primers.txt'
        exec(cmd, verbose=verbose)

with open(f"{out}/output_general.tsv", 'w') as general_output, open(f"{out}/output_genes.tsv", 'w') as genes_output:
    print(general_report, file=general_output)
    print(genes_report, file=genes_output)

if sample_to_warn:
    import sys
    sys.stderr.write ("WARNING: Hubo muestras con inserciones/deleciones menores a un triplete. Para más información ver el archivo log_warnings.txt\n")
    warnings_report = "Muestra\tGen\tPosicion\tCodon\tTipo\n"
    for i in sample_to_warn:
        warnings_report += str(i) + "\n"
    with open(f"{out}/log_warnings.txt", 'w') as warning_log:
        print(warnings_report, file=warning_log)
