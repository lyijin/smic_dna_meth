#!/usr/bin/env python3

"""
> filter_pos.four_criteria.py <

Script checks the "compiled_coverage_meth_unmeth.pre-filtered.tsv" file, and 
filters for positions that have median coverage >= 10 (criteria I) and present
in ALL samples (i.e. unmethylated + methylated reads > 5) (criteria II). To
make sure that lowly methylated loci are truly methylated, at least one of the 
treatments needs to have >= 1 methylated read in all replicates (criteria III).

The script reads in positions considered bona fide ("filter_miscalled_Cs.py" 
ran on "all.pre-filtered.merged.cov"), then filters out positions that are not
in this file. (criteria IV).

Afterwards, the script reads in merged_unfilt_covs/*.cov files, filters out 
unwanted positions, and saves the filtered files as ${a/cov/filt.cov} files.
"""
import csv
import glob
import gzip
import statistics
import sys

import numpy as np

def lazy_int(string):
    string = string.strip()
    return int(string) if string else 0

def calc_cols(input_list):
    return sum([1 for x in input_list if x])

# bit topsy-turvy, but start with criteria IV first!
print ('Checking positions that are considered bona fide (criteria IV)...', 
       file=sys.stderr)

# bona_fide_pos[scaf][pos]: 1 == bona fide; 0 == not bona fide
bona_fide_pos = {}
tsv_reader = csv.reader(open('../all.bona_fide_meth_pos.annot.c20.cov'), delimiter='\t')
for row in tsv_reader:
    scaf = row[0]
    pos = int(row[1])
    
    if scaf not in bona_fide_pos:
        scaf_size = int(scaf.split('size')[-1])
        # scaf_size, due to gap filling etc, tend to be MERELY SUGGESTIVE of
        # the actual length. the scaffolds are usually longer than labelled.
        # to counter this: create an array twice the expected length
        bona_fide_pos[scaf] = np.zeros(scaf_size * 2, dtype=np.int8)
    
    bona_fide_pos[scaf][pos] = 1

# check compiled, pre-filtered, meth & unmeth reads
print ('Checking full compiled table...', file=sys.stderr)
tsv_reader = csv.reader(gzip.open('compiled_coverage.pre-filt.meth_unmeth.tsv.gz', 'rt'), 
                        delimiter='\t')
next(tsv_reader)
for row in tsv_reader:
    if not row: continue
    
    per_sample_cov = [lazy_int(i) + lazy_int(j)
                      for i, j in zip(row[2::2], row[3::2])]
    
    # criteria I: coverage >= 5 for all treatments
    if not all([1 if x >= 5 else 0 for x in per_sample_cov]): continue
    
    # criteria II: median coverage >= 10
    if not statistics.median(per_sample_cov) >= 10: continue
    
    # criteria III: all replicates in at least one treatment have methylation
    per_sample_meth = [lazy_int(i) for i in row[2::2]]
    
    # for this dataset, there are 5 treatments, each with 3 replicates:
    # - 26A, B, C
    # - 26D, E, F
    # - 29D, E, F
    # - 32A, B, C
    # - 32D, E, F
    # in the per_sample_meth array, the treatments have indices of
    # 0,1,2; 3,4,5; 6,7,8; ...
    # following lines are hacky, but whatever works!
    one_treat_all_reps_meth = False
    for n in range(0, 15, 3):
        if all(per_sample_meth[n:n+3]): one_treat_all_reps_meth = True
    
    if not one_treat_all_reps_meth: continue
    
    scaf = row[0]
    pos = int(row[1])
    
    if scaf in bona_fide_pos:
        # positions passing n criteria would have a value of n
        bona_fide_pos[scaf][pos] += 3

# start filtering for correct positions
unfilt_files = sorted(glob.glob('subsampled_covs/*.subsamp.cov.gz'))

for u in unfilt_files:
    print ('Processing {}...'.format(u), file=sys.stderr)
    
    tsv_reader = csv.reader(gzip.open(u, 'rt'), delimiter='\t')
    output_file = u.split('/')[-1].replace('cov.gz', 'filt.cov')
    with open(output_file, 'w') as o:
        for row in tsv_reader:
            if not row: continue
            
            scaf = row[0]
            pos = int(row[1])
            
            if scaf in bona_fide_pos:
                if bona_fide_pos[scaf][pos] == 4:
                    print ('\t'.join(row), file=o)
