#!/usr/bin/env python3

"""
> tally_overlaps.py <

Tally numbers for genes that are over/under/not-differentially expressed
vs. meth/unmethylated.
"""
import gzip

import pandas as pd

# define dictionary that maps shortened smic gene names to full gene names
full_genename = {}
for line in gzip.open('smic_full_gene_names.txt.gz', 'rt'):
    fg = line.strip()
    full_genename[fg.split(' ')[0]] = fg

# get list of methylated genes
meth_genes = []
for line in open('smic_1pos.txt'):
    meth_genes.append(line.strip().split(' ')[0])

# read in expression TPMs: [0, inf)
data = pd.read_table('../../rnaseq/sleuth_results.26Cv32C.tsv',
                     index_col=0)

# drop NAs (i.e., restrict tally to expressed genes)
data = data.dropna(how='any')

# split dataset into methylated vs. unmethylated genes
meth_data = data[data.index.isin(meth_genes)]
sig_meth_data = meth_data[meth_data['qval'] < 0.05]
sig_up_meth_data = sig_meth_data[sig_meth_data['b'] > 0]
sig_down_meth_data = sig_meth_data[sig_meth_data['b'] < 0]

unmeth_data = data[~data.index.isin(meth_genes)]
sig_unmeth_data = unmeth_data[unmeth_data['qval'] < 0.05]
sig_up_unmeth_data = sig_unmeth_data[sig_unmeth_data['b'] > 0]
sig_down_unmeth_data = sig_unmeth_data[sig_unmeth_data['b'] < 0]

# yep, sufficient info to start printing tallies out
print('# expressed, methylated genes:', len(meth_data), sep='\t')
print ('# sig overexpressed at 32C, methylated genes:',
       len(sig_up_meth_data), sep='\t')
print ('# not sig diff expressed at 32C, methylated genes:',
       len(meth_data) - len(sig_meth_data), sep='\t')
print ('# sig underexpressed at 32C, methylated genes:',
       len(sig_down_meth_data), sep='\t')
print ()
print('# expressed, unmethylated genes:', len(unmeth_data), sep='\t')
print ('# sig overexpressed at 32C, unmethylated genes:',
       len(sig_up_unmeth_data), sep='\t')
print ('# not sig diff expressed at 32C, unmethylated genes:',
       len(unmeth_data) - len(sig_unmeth_data), sep='\t')
print ('# sig underexpressed at 32C, unmethylated genes:',
       len(sig_down_unmeth_data), sep='\t')

# dump the gene names into plain text files for topGO analysis
print (*[full_genename[x] for x in data.index], sep='\n',
       file=open('topGO_gene_lists/deg32.universe.txt', 'w'))
print (*[full_genename[x] for x in sig_up_meth_data.index], sep='\n',
       file=open('topGO_gene_lists/deg32.up_meth.txt', 'w'))
print (*[full_genename[x] for x in sig_down_meth_data.index], sep='\n',
       file=open('topGO_gene_lists/deg32.down_meth.txt', 'w'))
print (*[full_genename[x] for x in sig_up_unmeth_data.index], sep='\n',
       file=open('topGO_gene_lists/deg32.up_unmeth.txt', 'w'))
print (*[full_genename[x] for x in sig_down_unmeth_data.index], sep='\n',
       file=open('topGO_gene_lists/deg32.down_unmeth.txt', 'w'))
