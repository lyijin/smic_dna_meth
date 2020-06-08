#!/usr/bin/env python3

"""
> parse_mean_coverages.py <

Calculate per-gene, per-exon AND per-intron mean coverages in the form
    gene \t exon/intron \t exon/intron number \t exon/intron reversed \t RPKM

All numbers sum up to 1 million, so that genes can be compared between
individual files.
"""
import csv
import glob
import gzip

import numpy as np

import natural_sort
import parse_gff3

GENOME_GFF3_FILE = '../../raw_data/Smic.genome.annotation.gff3'
GENOME_SCAF_LEN_FILE = './smic_genome_scaffold_lengths.tsv'
HISAT2_SAMTOOLS_DEPTH_FILES = glob.glob('../../rnaseq/hisat2/*.depth.tsv.gz')

# read coordinates of genes and exons from .gff3 file.
scaffold_gff3 = parse_gff3.parse_gff3(open(GENOME_GFF3_FILE), 'exon')

# as genes might contain overlapping isoforms, the longest isoform is chosen,
# if multiples exist.
scaffold_gff3 = parse_gff3.pick_longest_mRNA(scaffold_gff3)

# make sure features in all mRNAs are sorted properly (for exon numbering).
scaffold_gff3 = parse_gff3.sort_features(scaffold_gff3)

# create dictionary which has all exons/introns as the key, and
# (start_coord, end_coord) as the value
all_ei = {}
for scaf in scaffold_gff3:
    for gene_id in scaffold_gff3[scaf]:
        sole_mRNA = list(scaffold_gff3[scaf][gene_id].mRNAs.values())[0]
        
        exon_coords = sole_mRNA.details['exon']
        for n, e in enumerate(exon_coords):
            k = f'{gene_id}_Exon_{n+1}_{n-len(exon_coords)}'
            all_ei[k] = (scaf, min(e), max(e))
        
        intron_coords = [(x[1], y[0]) for x, y in
                         zip(exon_coords[:-1], exon_coords[1:])]
        for n, i in enumerate(intron_coords):
            k = f'{gene_id}_Intron_{n+1}_{n-len(intron_coords)}'
            all_ei[k] = (scaf, min(i), max(i))

# get scaffold lengths for smic
scaf_lengths = {}
tsv_reader = csv.reader(open(GENOME_SCAF_LEN_FILE), delimiter='\t')
for row in tsv_reader:
    scaf_lengths[row[0]] = int(row[1])

# iterate through depth files, and calculate normalised exon/intron coverages
for hsd in HISAT2_SAMTOOLS_DEPTH_FILES:
    # create empty coverage array to store samtools depth output
    samtools_depth = {}
    for scaf in scaf_lengths:
        samtools_depth[scaf] = np.zeros(scaf_lengths[scaf], np.int32)
    
    # parse samtools depth coverage values into numpy array
    tsv_reader = csv.reader(gzip.open(hsd, 'rt'), delimiter='\t')
    for row in tsv_reader:
        cov = int(row[2])
        if not cov: continue
        
        scaf = row[0]
        pos = int(row[1]) - 1
        samtools_depth[scaf][pos] = cov
    
    # calculate mean coverage in every exon/intron
    prenorm_cov = {}
    for ei in all_ei:
        scaf, startpos, endpos = all_ei[ei]
        
        cov_slice = samtools_depth[scaf][startpos:endpos]
        nonzero_covs = cov_slice[cov_slice > 0]
        prenorm_cov[ei] = np.mean(nonzero_covs) if np.any(nonzero_covs) else 0
    
    total_coverage = sum(prenorm_cov.values())
    scaling_factor = 1000000 / total_coverage
    
    # print normalised exon/intron expression
    with open(hsd + '.ei_norm.tsv', 'w') as o:
        for ei in natural_sort.natural_sort(all_ei):
            print (*ei.split('_'), round(prenorm_cov[ei] * scaling_factor, 3),
                   sep='\t', file=o)
