#!/usr/bin/env python3

"""
> split_bismark_cov_by_context.py <

As Bismark covs produced by bismark_methylation_extractor with the --CX flag
has methylated Cs in all three contexts (CpG, CHG, CHH) mushed together
in the output file, this script will split them up by context.

Produces three output files by replacing the original .cov into .CpG.cov,
.CHG.cov and .CHH.cov.
"""
import argparse
import csv
import re

import numpy as np

import parse_fasta

parser = argparse.ArgumentParser(description="""
As Bismark covs produced by bismark_methylation_extractor with the --CX flag
has methylated Cs in all three contexts (CpG, CHG, CHH) mushed together
in the output file, this script will split them up by context.""")

parser.add_argument('genome_fasta', metavar='fasta_file',
                    type=argparse.FileType('r'),
                    help='genome FASTA file.')
parser.add_argument('bismark_cov', metavar="cov_filename",
                    type=argparse.FileType('r'),
                    help="Bismark .cov filename.")
args = parser.parse_args()

def reverse_complement(seq):
    seq = seq.replace('U', 'T')

    translation_from = 'AaTtGgCcYyRrSsWwKkMmBbDdHhVvNn'
    translation_to   = 'TtAaCcGgRrYySsWwMmKkVvHhDdBbNn'
    translation_table = str.maketrans(translation_from, translation_to)

    seq = seq[::-1].translate(translation_table)

    return seq

# read sequences
genome_sequences = parse_fasta.get_all_sequences(args.genome_fasta, 'fasta')

# genomic context is denoted in a NumPy array as follows:
#   0: not a cytosine on Watson/Crick strand
#   1: cytosine in CpG context on Watson/Crick strand
#   2: cytosine in CHG context on Watson/Crick strand
#   3: cytosine in CHH context on Watson/Crick strand
genomic_context = {}
for s in genome_sequences:
    genomic_context[s] = np.zeros(len(genome_sequences[s]), np.int8)
    
    scaf_watson = genome_sequences[s]
    scaf_crick = reverse_complement(genome_sequences[s])
    
    # get CpG
    for x in re.finditer('CG', scaf_watson):
        genomic_context[s][x.span()[0]] = 1
    
    for x in re.finditer('CG', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 1
    
    # get CHG
    for x in re.finditer('C[A|C|T]G', scaf_watson):
        genomic_context[s][x.span()[0]] = 2
    
    for x in re.finditer('C[A|C|T]G', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 2
    
    # get CHH--note the lookahead assertion: by default, Python only does non-
    # overlapping matches, which causes stuff like CCCC to report only 1 CHH.
    # lookahead assertion would correctly report 2 CHHs in CCCC.
    for x in re.finditer('(?=C[A|C|T][A|C|T])', scaf_watson):
        genomic_context[s][x.span()[0]] = 3
    
    for x in re.finditer('(?=C[A|C|T][A|C|T])', scaf_crick):
        genomic_context[s][-x.span()[0] - 1] = 3

# parse cov file and split them into component files
tsv_reader = csv.reader(args.bismark_cov, delimiter='\t')

output_file = {}
cpg_output_file = args.bismark_cov.name[:-4] + '.cpg' + args.bismark_cov.name[-4:]
output_file[1] = open(cpg_output_file, 'w')

chg_output_file = args.bismark_cov.name[:-4] + '.chg' + args.bismark_cov.name[-4:]
output_file[2] = open(chg_output_file, 'w')

chh_output_file = args.bismark_cov.name[:-4] + '.chh' + args.bismark_cov.name[-4:]
output_file[3] = open(chh_output_file, 'w')

for row in tsv_reader:
    if not row: continue
    
    scaf = row[0]
    pos = int(row[1]) - 1   # convert to 0-based numbering
    
    if genomic_context[scaf][pos]:
        print (*row, sep='\t', file=output_file[genomic_context[scaf][pos]])

for n in range(1, 4):
    output_file[n].close()
